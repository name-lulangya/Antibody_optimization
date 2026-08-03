"""Tests for safe staged-to-final file transactions."""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import unittest
import uuid
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from antibody_optimization.file_transaction import (  # noqa: E402
    FileTransactionError,
    PathSafetyError,
    replace_staged_files,
    validate_file_paths,
)


class FilePathValidationTests(unittest.TestCase):
    def test_rejects_exact_source_target_collision(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            project = Path(temporary_dir) / "project"
            project.mkdir()
            source = project / "source.txt"
            source.write_text("source", encoding="utf-8")

            with self.assertRaisesRegex(PathSafetyError, "Source/target paths collide"):
                validate_file_paths(
                    project_root=project,
                    source_paths=[source],
                    target_paths=[source],
                )

    def test_rejects_source_target_ancestor_collision(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            project = Path(temporary_dir) / "project"
            source_directory = project / "source-tree"
            source_directory.mkdir(parents=True)
            nested_target = source_directory / "derived.txt"

            with self.assertRaisesRegex(PathSafetyError, "by ancestry"):
                validate_file_paths(
                    project_root=project,
                    source_paths=[source_directory],
                    target_paths=[nested_target],
                )

    def test_rejects_target_outside_project(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            base = Path(temporary_dir)
            project = base / "project"
            project.mkdir()

            with self.assertRaisesRegex(PathSafetyError, "outside the project root"):
                validate_file_paths(
                    project_root=project,
                    target_paths=[base / "outside.txt"],
                )

    def test_rejects_existing_target_directory(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            project = Path(temporary_dir) / "project"
            target_directory = project / "not-a-file"
            target_directory.mkdir(parents=True)

            with self.assertRaisesRegex(PathSafetyError, "not a regular file"):
                validate_file_paths(
                    project_root=project,
                    target_paths=[target_directory],
                )

    def test_rejects_existing_target_symlink(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            project = Path(temporary_dir) / "project"
            project.mkdir()
            real_file = project / "real.txt"
            real_file.write_text("real", encoding="utf-8")
            link = project / "linked-target.txt"
            try:
                link.symlink_to(real_file)
            except OSError as exc:
                self.skipTest(f"Symbolic links are unavailable: {exc}")

            with self.assertRaisesRegex(PathSafetyError, "symbolic link"):
                validate_file_paths(project_root=project, target_paths=[link])


class StagedFileTransactionTests(unittest.TestCase):
    def test_candidate_copy_failure_does_not_mutate_inputs_or_targets(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            base = Path(temporary_dir)
            project = base / "project"
            stage_directory = base / "private-staging"
            project.mkdir()
            stage_directory.mkdir()
            final = project / "result.txt"
            stage = stage_directory / "result.txt"
            final.write_text("old", encoding="utf-8")
            stage.write_text("new", encoding="utf-8")

            with mock.patch(
                "antibody_optimization.file_transaction.shutil.copyfile",
                side_effect=OSError("injected candidate-copy failure"),
            ):
                with self.assertRaisesRegex(
                    FileTransactionError, "installation candidates"
                ):
                    replace_staged_files([(stage, final)], project_root=project)

            self.assertEqual(final.read_text(encoding="utf-8"), "old")
            self.assertEqual(stage.read_text(encoding="utf-8"), "new")
            self.assertEqual(list(project.glob(".file-transaction-*")), [])

    def test_injected_mid_transaction_failure_restores_every_old_file(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            base = Path(temporary_dir)
            project = base / "project"
            stage_directory = base / "staging-outside-project"
            project.mkdir()
            stage_directory.mkdir()

            final_a = project / "a.txt"
            final_b = project / "b.txt"
            stage_a = stage_directory / "a.txt"
            stage_b = stage_directory / "b.txt"
            final_a.write_text("old-a", encoding="utf-8")
            final_b.write_text("old-b", encoding="utf-8")
            stage_a.write_text("new-a", encoding="utf-8")
            stage_b.write_text("new-b", encoding="utf-8")

            call_count = 0

            def fail_once_on_fourth_replace(source: str | os.PathLike[str], target: str | os.PathLike[str]) -> None:
                nonlocal call_count
                call_count += 1
                if call_count == 4:
                    raise OSError("injected replacement failure")
                os.replace(source, target)

            with self.assertRaisesRegex(FileTransactionError, "all original targets"):
                replace_staged_files(
                    [(stage_a, final_a), (stage_b, final_b)],
                    project_root=project,
                    replace_func=fail_once_on_fourth_replace,
                )

            self.assertEqual(final_a.read_text(encoding="utf-8"), "old-a")
            self.assertEqual(final_b.read_text(encoding="utf-8"), "old-b")
            self.assertEqual(stage_a.read_text(encoding="utf-8"), "new-a")
            self.assertEqual(stage_b.read_text(encoding="utf-8"), "new-b")
            self.assertEqual(list(project.glob(".file-transaction-*.backup")), [])
            self.assertEqual(list(project.glob(".file-transaction-*.install")), [])

    def test_successfully_replaces_existing_and_new_targets(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            base = Path(temporary_dir)
            project = base / "project"
            stage_directory = base / "staging-outside-project"
            project.mkdir()
            stage_directory.mkdir()

            final_existing = project / "existing.txt"
            final_new = project / "new.txt"
            stage_existing = stage_directory / "existing.txt"
            stage_new = stage_directory / "new.txt"
            protected_source = base / "immutable-source.docx"
            final_existing.write_text("old", encoding="utf-8")
            stage_existing.write_text("replacement", encoding="utf-8")
            stage_new.write_text("created", encoding="utf-8")
            protected_source.write_text("source", encoding="utf-8")

            replace_staged_files(
                {
                    stage_existing: final_existing,
                    stage_new: final_new,
                },
                project_root=project,
                protected_source_paths=[protected_source],
            )

            self.assertEqual(final_existing.read_text(encoding="utf-8"), "replacement")
            self.assertEqual(final_new.read_text(encoding="utf-8"), "created")
            self.assertFalse(stage_existing.exists())
            self.assertFalse(stage_new.exists())
            self.assertEqual(protected_source.read_text(encoding="utf-8"), "source")
            self.assertEqual(list(project.glob(".file-transaction-*.backup")), [])
            self.assertEqual(list(project.glob(".file-transaction-*.install")), [])

    def test_installs_from_candidate_in_final_parent(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            base = Path(temporary_dir)
            project = base / "project"
            stage_directory = base / "private-staging"
            project.mkdir()
            stage_directory.mkdir()
            stage = stage_directory / "result.txt"
            final = project / "result.txt"
            stage.write_text("result", encoding="utf-8")
            install_sources: list[Path] = []

            def observe_replace(
                source: str | os.PathLike[str], target: str | os.PathLike[str]
            ) -> None:
                if Path(target) == final:
                    install_sources.append(Path(source))
                os.replace(source, target)

            replace_staged_files(
                [(stage, final)],
                project_root=project,
                replace_func=observe_replace,
            )

            self.assertEqual(final.read_text(encoding="utf-8"), "result")
            self.assertFalse(stage.exists())
            self.assertEqual(len(install_sources), 1)
            self.assertEqual(install_sources[0].parent, final.parent)
            self.assertNotEqual(install_sources[0], stage)
            self.assertEqual(list(project.glob(".file-transaction-*.install")), [])

    @unittest.skipUnless(os.name == "nt", "Windows ACL regression test")
    def test_windows_final_file_inherits_target_parent_acl(self) -> None:
        target_directory = ROOT / f".acl-regression-{uuid.uuid4().hex}"
        target_directory.mkdir()
        final = target_directory / "result.txt"
        try:
            with tempfile.TemporaryDirectory() as private_stage_directory:
                stage = Path(private_stage_directory) / "result.txt"
                stage.write_text("result", encoding="utf-8")
                replace_staged_files([(stage, final)], project_root=ROOT)

            powershell = (
                "$parentAcl = Get-Acl -LiteralPath $env:ACL_REGRESSION_PARENT\n"
                "$fileAcl = Get-Acl -LiteralPath $env:ACL_REGRESSION_FILE\n"
                "[PSCustomObject]@{"
                "ParentProtected=$parentAcl.AreAccessRulesProtected;"
                "FileProtected=$fileAcl.AreAccessRulesProtected"
                "} | ConvertTo-Json -Compress"
            )
            result = subprocess.run(
                [
                    "powershell",
                    "-NoProfile",
                    "-NonInteractive",
                    "-Command",
                    powershell,
                ],
                capture_output=True,
                text=True,
                encoding="utf-8",
                env={
                    **os.environ,
                    "ACL_REGRESSION_PARENT": str(target_directory),
                    "ACL_REGRESSION_FILE": str(final),
                },
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            acl_status = json.loads(result.stdout)
            self.assertFalse(acl_status["ParentProtected"])
            self.assertFalse(acl_status["FileProtected"])
        finally:
            final.unlink(missing_ok=True)
            target_directory.rmdir()


if __name__ == "__main__":
    unittest.main()
