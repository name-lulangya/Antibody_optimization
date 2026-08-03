#!/usr/bin/env python3
"""Export the already-open Nb252 ChimeraX session without changing it.

Run this file from ChimeraX 1.12's command line with ``runscript``.  The
script never opens or saves a session, changes a model transform, or changes a
display/color setting.  It exports each atomic model in its native model frame
and exports the two non-reference models in the experimental NK2R--Nb252 frame.
The output directory must not already exist.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import sys
import tempfile
import time
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Iterable


EXPORTER_VERSION = "1.0.0"
EXPECTED_MODEL_NAMES = (
    "NK2R-252.pdb",
    "NK2R-NKA.pdb",
    "fold_2r_252_nomg_model_0.cif",
)
REFERENCE_MODEL_NAME = "NK2R-252.pdb"
BUILTIN_ORANGE_RGB = (255, 165, 0)
COLOR_FIELDS = (
    "model_id",
    "model_name",
    "residue_atomspec",
    "chimerax_chain_id",
    "mmcif_chain_id",
    "auth_seq_id",
    "insertion_code",
    "residue_name",
    "ribbon_display",
    "ribbon_rgba",
    "displayed_atom_count",
    "atom_rgba_histogram_json",
    "surface_rgba_histogram_json",
    "selected_in_saved_session",
    "exact_builtin_orange_rgb_channels",
    "exact_builtin_orange_rgba_channels",
    "orange_annotation_status",
)


class SessionExportError(RuntimeError):
    """Raised before committing an export that cannot be tied to the CXS."""


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-cxs", type=Path, required=True)
    parser.add_argument("--expected-source-sha256", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--expected-chimerax-major-minor", default="1.12")
    parser.add_argument(
        "--generated-at",
        help="ISO-8601 timestamp; defaults to current local time",
    )
    return parser.parse_args(argv)


def main(chimerax_session: object, argv: list[str] | None = None) -> int:
    started = time.perf_counter()
    args = parse_args(argv)
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(
        timespec="seconds"
    )
    source = _regular_source(args.source_cxs)
    expected_hash = args.expected_source_sha256.strip().lower()
    if len(expected_hash) != 64 or any(c not in "0123456789abcdef" for c in expected_hash):
        raise SessionExportError("--expected-source-sha256 must be 64 lowercase hex digits")
    source_before = _file_identity(source)
    if source_before["sha256"] != expected_hash:
        raise SessionExportError(
            "Source CXS hash mismatch: "
            f"expected {expected_hash}, observed {source_before['sha256']}"
        )

    output_dir = args.output_dir.expanduser().absolute()
    if output_dir.exists() or output_dir.is_symlink():
        raise FileExistsError(f"Refusing existing export directory: {output_dir}")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    if output_dir.parent.is_symlink() or not output_dir.parent.is_dir():
        raise SessionExportError(f"Unsafe export parent directory: {output_dir.parent}")

    try:
        from chimerax.atomic import AtomicStructure
        from chimerax.core import __version__ as chimerax_version
        from chimerax.core.commands import StringArg, run
    except ImportError as exc:
        raise SessionExportError(
            "This script must be run inside ChimeraX 1.12 with `runscript`"
        ) from exc

    if not str(chimerax_version).startswith(
        f"{args.expected_chimerax_major_minor}."
    ) and str(chimerax_version) != args.expected_chimerax_major_minor:
        raise SessionExportError(
            f"Expected ChimeraX {args.expected_chimerax_major_minor}.x, found "
            f"{chimerax_version}"
        )

    session_path = getattr(chimerax_session, "session_file_path", None)
    if not session_path:
        raise SessionExportError(
            "The current ChimeraX session has no recorded session_file_path; "
            "open the supplied CXS and rerun explicitly"
        )
    session_path = Path(session_path).expanduser().resolve(strict=True)
    if not os.path.samefile(source, session_path):
        raise SessionExportError(
            f"Open session path {session_path} is not --source-cxs {source}"
        )

    all_session_models = list(chimerax_session.models.list())
    atomic_models = list(chimerax_session.models.list(type=AtomicStructure))
    models_by_name: dict[str, object] = {}
    duplicates: list[str] = []
    for model in atomic_models:
        name = str(model.name)
        if name in models_by_name:
            duplicates.append(name)
        models_by_name[name] = model
    if duplicates or set(models_by_name) != set(EXPECTED_MODEL_NAMES):
        raise SessionExportError(
            "Atomic model names must be exactly the expected three; "
            f"expected={list(EXPECTED_MODEL_NAMES)!r}, "
            f"observed={sorted(str(model.name) for model in atomic_models)!r}, "
            f"duplicates={sorted(set(duplicates))!r}"
        )
    for name, model in models_by_name.items():
        if int(model.num_coordsets) != 1:
            raise SessionExportError(
                f"Model {name!r} must have one coordinate set, found {model.num_coordsets}"
            )
        if len(model.atoms) == 0 or len(model.residues) == 0:
            raise SessionExportError(f"Atomic model {name!r} is empty")
        if not _finite_matrix(model.position.matrix) or not _finite_matrix(
            model.scene_position.matrix
        ):
            raise SessionExportError(f"Model {name!r} has a non-finite transform")

    stage = Path(
        tempfile.mkdtemp(prefix=f".{output_dir.name}-stage-", dir=output_dir.parent)
    )
    try:
        export_records: list[dict[str, object]] = []
        reference = models_by_name[REFERENCE_MODEL_NAME]
        for name in EXPECTED_MODEL_NAMES:
            model = models_by_name[name]
            native_path = stage / _native_filename(name)
            _save_mmcif(run, StringArg, chimerax_session, model, model, native_path)
            native_record = _export_record(
                path=native_path,
                output_dir=stage,
                source_model=model,
                frame_model=model,
                frame_kind="native_model_frame",
            )
            export_records.append(native_record)
            if name != REFERENCE_MODEL_NAME:
                reference_path = stage / _reference_filename(name)
                _save_mmcif(
                    run,
                    StringArg,
                    chimerax_session,
                    model,
                    reference,
                    reference_path,
                )
                export_records.append(
                    _export_record(
                        path=reference_path,
                        output_dir=stage,
                        source_model=model,
                        frame_model=reference,
                        frame_kind="experimental_reference_model_frame",
                    )
                )

        color_rows, color_notes, surface_overall_colors = collect_residue_color_rows(
            [models_by_name[name] for name in EXPECTED_MODEL_NAMES]
        )
        color_path = stage / "cxs_residue_colors.csv"
        _write_csv(color_path, COLOR_FIELDS, color_rows)
        source_after = _file_identity(source)
        if source_before != source_after:
            raise SessionExportError(
                "Source CXS identity changed during export; staged output was not committed"
            )

        manifest = {
            "schema_version": 1,
            "exporter_version": EXPORTER_VERSION,
            "status": "pass",
            "generated_at": generated_at,
            "execution_contract": {
                "entry_point": "ChimeraX GUI command: runscript",
                "opened_or_saved_session": False,
                "model_transforms_or_display_modified": False,
                "existing_output_directory_allowed": False,
            },
            "software": {"name": "UCSF ChimeraX", "version": str(chimerax_version)},
            "source_cxs": {
                "path": str(source),
                **source_before,
                "session_file_path": str(session_path),
                "unchanged_after_export": True,
            },
            "strict_model_identity": {
                "expected_names": list(EXPECTED_MODEL_NAMES),
                "observed_names": [str(models_by_name[name].name) for name in EXPECTED_MODEL_NAMES],
                "matched_exactly": True,
            },
            "session_model_inventory": {
                "model_count": len(all_session_models),
                "atomic_structure_count": len(atomic_models),
                "non_atomic_model_count": len(all_session_models) - len(atomic_models),
                "models": [
                    _session_model_record(
                        model,
                        is_atomic_structure=isinstance(model, AtomicStructure),
                    )
                    for model in all_session_models
                ],
                "semantics": (
                    "All models returned by session.models.list(), including surfaces, "
                    "groups, and other non-AtomicStructure children when present."
                ),
            },
            "reference_frame": {
                "model_name": REFERENCE_MODEL_NAME,
                "model_id": str(reference.id_string),
                "meaning": "experimental NK2R--Nb252 model frame",
            },
            "models": [_model_record(models_by_name[name]) for name in EXPECTED_MODEL_NAMES],
            "atom_site_count_semantics": {
                "atom_object": "one ChimeraX Atom object",
                "atom_site": (
                    "one serialized coordinate state; nonblank Atom.alt_locs are "
                    "enumerated and occupancy is read with get_alt_loc_occupancy"
                ),
                "no_altloc": (
                    "when alt_locs is empty and active alt_loc is blank, one blank "
                    "site is counted using Atom.occupancy"
                ),
                "occupancy_tolerance": 1e-6,
                "occupancy_classes": [
                    "negative", "zero", "partial", "unit", "above_unit", "nonfinite"
                ],
            },
            "exports": export_records,
            "color_inventory": {
                "path": color_path.name,
                "sha256": _sha256(color_path),
                "size_bytes": color_path.stat().st_size,
                "encoding": "utf-8-sig",
                "builtin_orange_rgb": list(BUILTIN_ORANGE_RGB),
                "semantics": (
                    "Exact RGB/RGBA channel matches are candidate color evidence only; "
                    "they do not prove the collaborator's orange selection."
                ),
                "notes": color_notes,
                "unmapped_or_overall_surface_colors": surface_overall_colors,
                "custom_colors": _custom_colors(chimerax_session),
            },
        }
        manifest_path = stage / "cxs_export_manifest.json"
        _write_json(manifest_path, manifest)
        recorded_argv = list(sys.argv[1:] if argv is None else argv)
        if not args.generated_at:
            recorded_argv.extend(["--generated-at", generated_at])
        _write_json(
            stage / "cxs_export_run_summary.json",
            {
                "schema_version": 1,
                "status": "pass",
                "generated_at": generated_at,
                "elapsed_seconds": time.perf_counter() - started,
                "software": {"name": "UCSF ChimeraX", "version": str(chimerax_version)},
                "script": str(Path(__file__).resolve()),
                "entry_point": "ChimeraX command: runscript",
                "argv": recorded_argv,
                "source_cxs": {"path": str(source), **source_before},
                "manifest": {
                    "path": manifest_path.name,
                    "sha256": _sha256(manifest_path),
                    "size_bytes": manifest_path.stat().st_size,
                },
                "output_counts": {
                    "session_models": len(all_session_models),
                    "atomic_models": len(EXPECTED_MODEL_NAMES),
                    "non_atomic_models": len(all_session_models) - len(atomic_models),
                    "mmcif_exports": len(export_records),
                    "residue_color_rows": len(color_rows),
                    "files_excluding_run_summary": len(list(stage.iterdir())),
                },
            },
        )
        os.replace(stage, output_dir)
    except BaseException:
        if stage.exists():
            shutil.rmtree(stage)
        raise
    return 0


def collect_residue_color_rows(
    models: Iterable[object],
) -> tuple[list[dict[str, object]], list[str], list[dict[str, object]]]:
    """Collect exact residue color channels without interpreting annotations."""

    rows: list[dict[str, object]] = []
    notes: list[str] = []
    overall_surface_colors: list[dict[str, object]] = []
    for model in models:
        surface_by_residue: dict[int, Counter[tuple[int, int, int, int]]] = defaultdict(Counter)
        for surface in model.surfaces():
            vertex_colors = getattr(surface, "vertex_colors", None)
            try:
                vertex_to_atom = surface.vertex_to_atom_map()
            except Exception as exc:  # ChimeraX extension can reject stale surfaces.
                notes.append(f"{model.name}: surface {surface.name!r} mapping failed: {exc}")
                continue
            if vertex_to_atom is None:
                rgba = _rgba(getattr(surface, "color", ()))
                overall_surface_colors.append(
                    {
                        "model_name": str(model.name),
                        "model_id": str(model.id_string),
                        "surface_name": str(surface.name),
                        "surface_id": str(getattr(surface, "id_string", "")),
                        "overall_rgba": [] if rgba is None else list(rgba),
                        "residue_mapping_status": "unavailable_no_vertex_to_atom_map",
                    }
                )
                notes.append(f"{model.name}: surface {surface.name!r} has no vertex-to-atom map")
                continue
            atoms = surface.atoms
            if vertex_colors is None:
                uniform = _rgba(getattr(surface, "color", ()))
                if uniform is None:
                    notes.append(f"{model.name}: surface {surface.name!r} has no readable color")
                    continue
                for atom_index in set(int(value) for value in vertex_to_atom):
                    surface_by_residue[id(atoms[atom_index].residue)][uniform] += 1
            else:
                if len(vertex_colors) != len(vertex_to_atom):
                    raise SessionExportError(
                        f"Surface color/map length mismatch in {model.name}/{surface.name}"
                    )
                for color, atom_index in zip(vertex_colors, vertex_to_atom, strict=True):
                    rgba = _rgba(color)
                    if rgba is not None:
                        surface_by_residue[id(atoms[int(atom_index)].residue)][rgba] += 1

        for residue in model.residues:
            atom_hist: Counter[tuple[int, int, int, int]] = Counter()
            displayed_atom_count = 0
            selected = bool(getattr(residue, "selected", False))
            for atom in residue.atoms:
                rgba = _rgba(atom.color)
                if rgba is not None:
                    atom_hist[rgba] += 1
                displayed_atom_count += int(bool(atom.display))
                selected = selected or bool(getattr(atom, "selected", False))
            ribbon = _rgba(residue.ribbon_color)
            surface_hist = surface_by_residue.get(id(residue), Counter())
            channels = list(atom_hist) + list(surface_hist)
            if ribbon is not None:
                channels.append(ribbon)
            exact_rgb = any(color[:3] == BUILTIN_ORANGE_RGB for color in channels)
            exact_rgba = any(color == (*BUILTIN_ORANGE_RGB, 255) for color in channels)
            rows.append(
                {
                    "model_id": str(model.id_string),
                    "model_name": str(model.name),
                    "residue_atomspec": str(residue.atomspec),
                    "chimerax_chain_id": str(residue.chain_id),
                    "mmcif_chain_id": str(getattr(residue, "mmcif_chain_id", "")),
                    "auth_seq_id": int(residue.number),
                    "insertion_code": _clean_code(residue.insertion_code),
                    "residue_name": str(residue.name),
                    "ribbon_display": bool(residue.ribbon_display),
                    "ribbon_rgba": _rgba_text(ribbon),
                    "displayed_atom_count": displayed_atom_count,
                    "atom_rgba_histogram_json": _histogram_json(atom_hist),
                    "surface_rgba_histogram_json": _histogram_json(surface_hist),
                    "selected_in_saved_session": selected,
                    "exact_builtin_orange_rgb_channels": exact_rgb,
                    "exact_builtin_orange_rgba_channels": exact_rgba,
                    "orange_annotation_status": "candidate_unconfirmed" if exact_rgb else "no_exact_match",
                }
            )
    return rows, notes, overall_surface_colors


def _save_mmcif(run, string_arg, session, model, frame_model, path: Path) -> None:
    command = " ".join(
        (
            "save",
            string_arg.unparse(str(path)),
            "format mmcif models",
            str(model.atomspec),
            "relModel",
            str(frame_model.atomspec),
            "selectedOnly false displayedOnly false allCoordsets false "
            "computedSheets false bestGuess false",
        )
    )
    run(session, command)
    if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
        raise SessionExportError(f"ChimeraX did not create a regular nonempty export: {path}")


def _export_record(*, path: Path, output_dir: Path, source_model, frame_model, frame_kind: str) -> dict[str, object]:
    return {
        "path": str(path.relative_to(output_dir)),
        "sha256": _sha256(path),
        "size_bytes": path.stat().st_size,
        "format": "mmCIF",
        "source_model_name": str(source_model.name),
        "source_model_id": str(source_model.id_string),
        "coordinate_frame_kind": frame_kind,
        "coordinate_frame_model_name": str(frame_model.name),
        "coordinate_frame_model_id": str(frame_model.id_string),
    }


def _model_record(model) -> dict[str, object]:
    parent = getattr(model, "parent", None)
    atom_site_counts = atom_site_classification_counts(model.atoms)
    return {
        "model_name": str(model.name),
        "model_id": str(model.id_string),
        "atomspec": str(model.atomspec),
        "atom_count": len(model.atoms),
        "atom_site_classification_counts": atom_site_counts,
        "residue_count": len(model.residues),
        "coordinate_set_count": int(model.num_coordsets),
        "chain_count": len(model.chains),
        "parent_model": (
            None
            if parent is None
            else {
                "model_id": str(getattr(parent, "id_string", "")),
                "model_name": str(getattr(parent, "name", "")),
                "atomspec": str(getattr(parent, "atomspec", "")),
            }
        ),
        "position_3x4": _matrix(model.position.matrix),
        "scene_position_3x4": _matrix(model.scene_position.matrix),
    }


def _session_model_record(model, *, is_atomic_structure: bool) -> dict[str, object]:
    """Return a conservative inventory row for any saved-session model."""

    parent = getattr(model, "parent", None)
    children = getattr(model, "child_models", ())
    try:
        child_count = len(children()) if callable(children) else len(children)
    except (TypeError, AttributeError):
        child_count = 0
    return {
        "model_id": str(getattr(model, "id_string", "")),
        "model_name": str(getattr(model, "name", "")),
        "atomspec": str(getattr(model, "atomspec", "")),
        "python_class": type(model).__name__,
        "is_atomic_structure": bool(is_atomic_structure),
        "parent_model_id": (
            "" if parent is None else str(getattr(parent, "id_string", ""))
        ),
        "parent_model_name": (
            "" if parent is None else str(getattr(parent, "name", ""))
        ),
        "child_model_count": child_count,
        "display": bool(getattr(model, "display", False)),
    }


def atom_site_classification_counts(atoms: Iterable[object]) -> dict[str, int]:
    """Count serialized atom sites, alternate locations, and occupancy classes."""

    counts = Counter()
    for atom in atoms:
        active_altloc = _clean_code(getattr(atom, "alt_loc", ""))
        raw_altlocs = getattr(atom, "alt_locs", ())
        altlocs = sorted(
            {
                cleaned
                for value in raw_altlocs
                if (cleaned := _clean_code(value))
            }
            | ({active_altloc} if active_altloc else set())
        )
        if altlocs:
            counts["atom_objects_with_nonblank_altlocs"] += 1
            for altloc in altlocs:
                occupancy = float(atom.get_alt_loc_occupancy(altloc))
                counts["atom_site_count"] += 1
                counts["nonblank_altloc_site_count"] += 1
                counts[f"occupancy_{_occupancy_class(occupancy)}_site_count"] += 1
        else:
            occupancy = float(getattr(atom, "occupancy"))
            counts["atom_site_count"] += 1
            counts["blank_altloc_site_count"] += 1
            counts[f"occupancy_{_occupancy_class(occupancy)}_site_count"] += 1
    fields = (
        "atom_site_count", "blank_altloc_site_count", "nonblank_altloc_site_count",
        "atom_objects_with_nonblank_altlocs", "occupancy_negative_site_count",
        "occupancy_zero_site_count", "occupancy_partial_site_count",
        "occupancy_unit_site_count", "occupancy_above_unit_site_count",
        "occupancy_nonfinite_site_count",
    )
    return {field: int(counts[field]) for field in fields}


def _occupancy_class(value: float) -> str:
    if not math.isfinite(value):
        return "nonfinite"
    if value < -1e-6:
        return "negative"
    if math.isclose(value, 0.0, abs_tol=1e-6):
        return "zero"
    if value < 1.0 - 1e-6:
        return "partial"
    if math.isclose(value, 1.0, abs_tol=1e-6):
        return "unit"
    return "above_unit"


def _custom_colors(session) -> list[dict[str, object]]:
    result: list[dict[str, object]] = []
    for name in session.user_colors.list():
        result.append({"name": str(name), "rgba": list(_rgba(session.user_colors[name].uint8x4()))})
    return result


def _native_filename(name: str) -> str:
    return f"{Path(name).stem}__native.cif"


def _reference_filename(name: str) -> str:
    return f"{Path(name).stem}__in_NK2R-252_frame.cif"


def _file_identity(path: Path) -> dict[str, object]:
    stat = path.stat()
    return {"sha256": _sha256(path), "size_bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def _regular_source(path: Path) -> Path:
    lexical = path.expanduser().absolute()
    if lexical.is_symlink() or not lexical.is_file():
        raise SessionExportError(f"Source CXS must be a regular non-symlink file: {lexical}")
    return lexical.resolve(strict=True)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _matrix(value) -> list[list[float]]:
    return [[float(item) for item in row] for row in value]


def _finite_matrix(value) -> bool:
    return all(math.isfinite(float(item)) for row in value for item in row)


def _rgba(value) -> tuple[int, int, int, int] | None:
    try:
        values = tuple(int(item) for item in value)
    except (TypeError, ValueError):
        return None
    return values if len(values) == 4 and all(0 <= item <= 255 for item in values) else None


def _rgba_text(value: tuple[int, int, int, int] | None) -> str:
    return "" if value is None else ",".join(str(item) for item in value)


def _histogram_json(histogram: Counter[tuple[int, int, int, int]]) -> str:
    serializable = {",".join(map(str, rgba)): count for rgba, count in sorted(histogram.items())}
    return json.dumps(serializable, sort_keys=True, separators=(",", ":"))


def _clean_code(value: object) -> str:
    text = str(value or "")
    return "" if text in {" ", "\x00", ".", "?"} else text


def _write_csv(path: Path, fields: Iterable[str], rows: Iterable[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            extrasaction="raise",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, value: object) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
        newline="\n",
    )


if __name__ == "__main__":
    try:
        active_session = session  # type: ignore[name-defined]  # provided by ChimeraX
    except NameError as exc:
        raise SystemExit(
            "Run this script inside ChimeraX 1.12 with the `runscript` command"
        ) from exc
    raise SystemExit(main(active_session, sys.argv[1:]))
