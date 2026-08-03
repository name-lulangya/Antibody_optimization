"""Path-safety checks and rollback-capable staged file replacement.

The helpers in this module are deliberately independent of any experiment
schema.  They protect immutable inputs from output-path collisions, constrain
outputs to a declared project root, and replace a set of staged files while
restoring the complete pre-transaction state if a replacement fails.

The multi-file transaction provides failure atomicity through backups and
rollback.  Like any sequence of filesystem operations, it does not make all
paths change at the same instant and does not guarantee crash durability.
"""

from __future__ import annotations

import os
import shutil
import tempfile
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import TypeAlias


Pathish: TypeAlias = str | os.PathLike[str]
ReplaceFunc: TypeAlias = Callable[[Pathish, Pathish], object]
StagedFilePairs: TypeAlias = Mapping[Pathish, Pathish] | Iterable[tuple[Pathish, Pathish]]


class PathSafetyError(ValueError):
    """Raised before mutation when a path plan is ambiguous or unsafe."""


class FileTransactionError(RuntimeError):
    """Raised when a staged replacement cannot be committed safely."""


@dataclass(frozen=True)
class ValidatedPaths:
    """Canonical paths returned by :func:`validate_file_paths`."""

    project_root: Path
    source_paths: tuple[Path, ...]
    target_paths: tuple[Path, ...]


def validate_file_paths(
    *,
    project_root: Pathish,
    source_paths: Iterable[Pathish] = (),
    target_paths: Iterable[Pathish],
) -> ValidatedPaths:
    """Validate read-only inputs and file outputs before a workflow mutates data.

    Relative paths are interpreted relative to ``project_root``.  Source paths
    may resolve outside the project, but every target must resolve inside it.
    A source and target, or two targets, may not be equal and neither may be an
    ancestor of the other.  Any existing target must be a regular file and may
    not itself be a symbolic link.

    The returned paths are absolute and resolved.  This also makes collisions
    through existing symbolic-link aliases visible to the comparison logic.
    """

    try:
        root = Path(project_root).expanduser().resolve(strict=True)
    except OSError as exc:
        raise PathSafetyError(f"Cannot resolve project root {project_root!r}: {exc}") from exc
    if not root.is_dir():
        raise PathSafetyError(f"Project root is not a directory: {root}")

    raw_sources = tuple(source_paths)
    raw_targets = tuple(target_paths)
    if not raw_targets:
        raise PathSafetyError("At least one target path is required")

    sources = tuple(_resolved_path(path, root) for path in raw_sources)
    target_lexical = tuple(_absolute_lexical_path(path, root) for path in raw_targets)
    targets = tuple(path.resolve(strict=False) for path in target_lexical)

    for lexical, target in zip(target_lexical, targets, strict=True):
        if lexical.is_symlink():
            raise PathSafetyError(f"Existing target must not be a symbolic link: {lexical}")
        if lexical.exists() and not lexical.is_file():
            raise PathSafetyError(f"Existing target is not a regular file: {lexical}")
        if not _is_within(target, root):
            raise PathSafetyError(
                f"Target resolves outside the project root: {lexical} -> {target}"
            )

    _reject_pairwise_collisions(targets, label="target")
    for source in sources:
        for target in targets:
            if _paths_collide(source, target):
                raise PathSafetyError(
                    "Source/target paths collide exactly or by ancestry: "
                    f"{source} <-> {target}"
                )

    return ValidatedPaths(root, sources, targets)


def replace_staged_files(
    staged_to_final: StagedFilePairs,
    *,
    project_root: Pathish,
    protected_source_paths: Iterable[Pathish] = (),
    replace_func: ReplaceFunc = os.replace,
) -> None:
    """Commit staged regular files to final paths with backup and rollback.

    Staged files may live outside ``project_root``.  Final paths may not.  The
    optional ``protected_source_paths`` are immutable workflow inputs (for
    example, the original DOCX) that must not collide exactly or by ancestry
    with any stage or final path.

    Each staged file is first copied to a unique candidate in the final file's
    parent directory.  This preserves the target directory's access-control
    inheritance on Windows and makes the final replacement same-filesystem.
    Windows candidates deliberately inherit the final parent ACL rather than
    copying the private staging ACL; POSIX candidates copy the staged file's
    permission mode.  Timestamps, extended attributes, and explicit ACLs are
    outside this content-artifact transaction contract.  Candidate preparation
    temporarily requires space for one additional copy of every staged file.
    Existing final files are then moved to unique same-directory backups, and
    the candidates are installed in input order.  If any replacement fails,
    installed candidates and all original targets are restored.  ``replace_func``
    must follow the ``os.replace`` contract: on a raised exception it must not
    have performed the replacement.  It is injectable so tests can force a
    failure at a precise operation.
    """

    pairs = _materialize_pairs(staged_to_final)
    if not pairs:
        raise PathSafetyError("At least one staged/final pair is required")

    raw_protected = tuple(protected_source_paths)
    raw_stages = tuple(stage for stage, _ in pairs)
    raw_finals = tuple(final for _, final in pairs)
    validated = validate_file_paths(
        project_root=project_root,
        source_paths=(*raw_protected, *raw_stages),
        target_paths=raw_finals,
    )
    protected_count = len(raw_protected)
    protected = validated.source_paths[:protected_count]
    stages = validated.source_paths[protected_count:]
    finals = validated.target_paths

    _reject_pairwise_collisions(stages, label="staged source")
    for raw_stage, stage in zip(raw_stages, stages, strict=True):
        for source in protected:
            if _paths_collide(stage, source):
                raise PathSafetyError(
                    "Staged/protected-source paths collide exactly or by ancestry: "
                    f"{stage} <-> {source}"
                )
        lexical_stage = _absolute_lexical_path(raw_stage, validated.project_root)
        if lexical_stage.is_symlink() or not stage.is_file():
            raise PathSafetyError(f"Staged source is not a regular non-symlink file: {stage}")

    for final in finals:
        if not final.parent.is_dir():
            raise PathSafetyError(f"Target parent directory does not exist: {final.parent}")

    canonical_pairs = list(zip(stages, finals, strict=True))
    install_candidates: list[Path] = []
    try:
        for stage, final in canonical_pairs:
            candidate = _reserve_install_path(final)
            install_candidates.append(candidate)
            shutil.copyfile(stage, candidate)
            if os.name != "nt":
                shutil.copymode(stage, candidate)
    except BaseException as exc:
        cleanup_errors = _remove_paths(install_candidates)
        detail = _format_errors(cleanup_errors)
        raise FileTransactionError(
            "Could not prepare same-directory installation candidates before "
            f"mutation{detail}"
        ) from exc

    install_pairs = list(zip(install_candidates, finals, strict=True))
    existed_before = {final: final.exists() for final in finals}
    backup_paths: dict[Path, Path] = {}
    placeholders: list[Path] = []

    try:
        for final in finals:
            if existed_before[final]:
                backup = _reserve_backup_path(final)
                backup_paths[final] = backup
                placeholders.append(backup)
    except BaseException as exc:
        cleanup_errors = _remove_paths([*placeholders, *install_candidates])
        detail = _format_errors(cleanup_errors)
        raise FileTransactionError(
            f"Could not reserve transaction backups before mutation{detail}"
        ) from exc

    backed_up: set[Path] = set()
    installed: list[tuple[Path, Path]] = []
    try:
        for _, final in canonical_pairs:
            if existed_before[final]:
                replace_func(final, backup_paths[final])
                backed_up.add(final)

        for candidate, final in install_pairs:
            replace_func(candidate, final)
            installed.append((candidate, final))
    except BaseException as exc:
        rollback_errors: list[str] = []

        for candidate, final in reversed(installed):
            try:
                replace_func(final, candidate)
            except BaseException as rollback_exc:
                rollback_errors.append(
                    "could not remove installed candidate "
                    f"{final} to {candidate}: {rollback_exc!r}"
                )

        for _, final in reversed(canonical_pairs):
            if final not in backed_up:
                continue
            backup = backup_paths[final]
            try:
                replace_func(backup, final)
            except BaseException as rollback_exc:
                rollback_errors.append(
                    f"could not restore original target {final} from {backup}: "
                    f"{rollback_exc!r}"
                )

        unused_placeholders = [
            path
            for final, path in backup_paths.items()
            if final not in backed_up
        ]
        rollback_errors.extend(
            _remove_paths([*unused_placeholders, *install_candidates])
        )
        if rollback_errors:
            raise FileTransactionError(
                "Staged replacement failed and rollback was incomplete; retained backup "
                f"files must not be discarded. Details: {'; '.join(rollback_errors)}"
            ) from exc
        raise FileTransactionError(
            "Staged replacement failed; all original targets and staged files were restored"
        ) from exc

    cleanup_errors = _remove_paths([*backup_paths.values(), *stages])
    if cleanup_errors:
        raise FileTransactionError(
            "Staged files were committed, but obsolete transaction files could not be "
            f"removed. Details: {'; '.join(cleanup_errors)}"
        )


def _materialize_pairs(staged_to_final: StagedFilePairs) -> list[tuple[Pathish, Pathish]]:
    if isinstance(staged_to_final, Mapping):
        return list(staged_to_final.items())
    return list(staged_to_final)


def _absolute_lexical_path(path: Pathish, root: Path) -> Path:
    candidate = Path(path).expanduser()
    if not candidate.is_absolute():
        candidate = root / candidate
    return Path(os.path.abspath(os.fspath(candidate)))


def _resolved_path(path: Pathish, root: Path) -> Path:
    return _absolute_lexical_path(path, root).resolve(strict=False)


def _is_within(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
    except ValueError:
        return False
    return True


def _paths_collide(first: Path, second: Path) -> bool:
    return first == second or first in second.parents or second in first.parents


def _reject_pairwise_collisions(paths: tuple[Path, ...], *, label: str) -> None:
    for index, first in enumerate(paths):
        for second in paths[index + 1 :]:
            if _paths_collide(first, second):
                raise PathSafetyError(
                    f"{label.capitalize()} paths collide exactly or by ancestry: "
                    f"{first} <-> {second}"
                )


def _reserve_backup_path(final: Path) -> Path:
    descriptor, name = tempfile.mkstemp(
        prefix=".file-transaction-",
        suffix=".backup",
        dir=final.parent,
    )
    os.close(descriptor)
    return Path(name)


def _reserve_install_path(final: Path) -> Path:
    descriptor, name = tempfile.mkstemp(
        prefix=".file-transaction-",
        suffix=".install",
        dir=final.parent,
    )
    os.close(descriptor)
    return Path(name)


def _remove_paths(paths: Iterable[Path]) -> list[str]:
    errors: list[str] = []
    for path in paths:
        try:
            path.unlink(missing_ok=True)
        except BaseException as exc:
            errors.append(f"could not remove {path}: {exc!r}")
    return errors


def _format_errors(errors: Iterable[str]) -> str:
    materialized = list(errors)
    return "" if not materialized else f"; cleanup errors: {'; '.join(materialized)}"
