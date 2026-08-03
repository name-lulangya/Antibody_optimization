"""CLI safety and reproducibility tests for expression-data preparation."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest
from contextlib import contextmanager
from decimal import Decimal
from pathlib import Path
from types import ModuleType
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "nb序列及产量（1L）.docx"
SOURCE_SHA256 = "a6e4022f0978fbd70a0e04dc78f479140ab6f55caaa90b467fb77a62eb5db5d1"
PREPARE = ROOT / "scripts" / "data_preparation" / "prepare_nb_expression_dataset.py"
PLOT = ROOT / "scripts" / "data_preparation" / "plot_nb_expression_qc.py"
VERIFY = ROOT / "scripts" / "data_preparation" / "verify_nb_expression_outputs.py"
OUTPUT_NAMES = [
    "samples.csv",
    "yield_observations.csv",
    "assay_context.csv",
    "nb_expression_records.csv",
    "raw_transcription.csv",
    "nb_expression_sequences.fasta",
    "qc_plot_data.csv",
    "nb_expression_qc.svg",
    "validation_report.json",
    "manifest.json",
]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_script(path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(f"_test_{path.stem}", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load script for testing: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@contextmanager
def mocked_symlink_endpoint(link: Path, resolved_target: Path):
    """Model a leaf symlink without requiring OS symlink privileges."""

    original_resolve = Path.resolve
    original_is_file = Path.is_file
    original_is_symlink = Path.is_symlink

    def resolve(path: Path, strict: bool = False) -> Path:
        if path == link:
            return resolved_target
        return original_resolve(path, strict=strict)

    def is_symlink(path: Path) -> bool:
        if path == link:
            return True
        return original_is_symlink(path)

    def is_file(path: Path) -> bool:
        if path == link:
            return resolved_target.is_file()
        return original_is_file(path)

    with mock.patch.object(Path, "resolve", resolve), mock.patch.object(
        Path, "is_file", is_file
    ), mock.patch.object(Path, "is_symlink", is_symlink):
        yield


class ExpressionCliTests(unittest.TestCase):
    def run_command(self, *arguments: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, *arguments],
            cwd=ROOT,
            capture_output=True,
            text=True,
            encoding="utf-8",
            env={**os.environ, "PYTHONUTF8": "1"},
            check=False,
        )

    def prepare_arguments(self, output: Path, summary: Path) -> list[str]:
        return [
            str(PREPARE),
            "--source",
            str(SOURCE),
            "--output-dir",
            str(output),
            "--run-summary",
            str(summary),
            "--expected-source-sha256",
            SOURCE_SHA256,
            "--document-title-culture-volume-l",
            "1",
            "--generated-at",
            "2026-08-03T12:01:38+08:00",
        ]

    def test_prepare_rejects_run_summary_output_collision(self) -> None:
        source_before = digest(SOURCE)
        with tempfile.TemporaryDirectory(prefix=".cli-test-", dir=ROOT) as temporary_dir:
            output = Path(temporary_dir) / "artifacts"
            collision = output / "samples.csv"
            result = self.run_command(*self.prepare_arguments(output, collision))
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(collision.exists())
        self.assertEqual(digest(SOURCE), source_before)

    def test_plot_rejects_input_output_collision_without_modification(self) -> None:
        with tempfile.TemporaryDirectory(prefix=".cli-test-", dir=ROOT) as temporary_dir:
            plot_data = Path(temporary_dir) / "qc.csv"
            original = b"metric,category,count\nsource_records,LTT,23\n"
            plot_data.write_bytes(original)
            result = self.run_command(
                str(PLOT),
                "--plot-data",
                str(plot_data),
                "--output",
                str(plot_data),
                "--source-sha256",
                SOURCE_SHA256,
                "--overwrite",
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(plot_data.read_bytes(), original)

    def test_verifier_rejects_report_inside_artifact_directory(self) -> None:
        with tempfile.TemporaryDirectory(prefix=".cli-test-", dir=ROOT) as temporary_dir:
            artifact_dir = Path(temporary_dir) / "artifacts"
            artifact_dir.mkdir()
            collision = artifact_dir / "samples.csv"
            collision.write_text("sentinel", encoding="utf-8")
            result = self.run_command(
                str(VERIFY),
                "--source",
                str(SOURCE),
                "--artifact-dir",
                str(artifact_dir),
                "--report",
                str(collision),
                "--expected-source-sha256",
                SOURCE_SHA256,
                "--document-title-culture-volume-l",
                "1",
                "--overwrite",
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(collision.read_text(encoding="utf-8"), "sentinel")

    def test_existing_directory_target_is_rejected_before_changes(self) -> None:
        with tempfile.TemporaryDirectory(prefix=".cli-test-", dir=ROOT) as temporary_dir:
            base = Path(temporary_dir)
            output = base / "artifacts"
            output.mkdir()
            sentinel = output / "samples.csv"
            sentinel.write_text("old-samples", encoding="utf-8")
            (output / "nb_expression_records.csv").mkdir()
            summary = base / "summary.json"
            result = self.run_command(
                *self.prepare_arguments(output, summary), "--overwrite"
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(sentinel.read_text(encoding="utf-8"), "old-samples")
            self.assertFalse(summary.exists())

    def test_prepare_preserves_lexical_summary_path_for_symlink_rejection(self) -> None:
        module = load_script(PREPARE)
        with tempfile.TemporaryDirectory(prefix=".cli-test-", dir=ROOT) as temporary_dir:
            base = Path(temporary_dir)
            output = base / "artifacts"
            link = base / "summary-link.json"
            resolved_target = base / "real-summary.json"
            arguments = argparse.Namespace(
                source=SOURCE,
                output_dir=output,
                run_summary=link,
                expected_source_sha256=SOURCE_SHA256,
                document_title_culture_volume_l=Decimal("1"),
                generated_at="2026-08-03T12:01:38+08:00",
                overwrite=False,
            )
            with mocked_symlink_endpoint(link, resolved_target), mock.patch.object(
                module, "parse_args", return_value=arguments
            ), mock.patch.object(
                module,
                "parse_expression_docx",
                side_effect=AssertionError("parser must not run for a symlink target"),
            ):
                with self.assertRaisesRegex(ValueError, "symbolic link"):
                    module.main()
            self.assertFalse(output.exists())
            self.assertFalse(resolved_target.exists())

    def test_prepare_preserves_lexical_source_path_for_symlink_rejection(self) -> None:
        module = load_script(PREPARE)
        with tempfile.TemporaryDirectory(prefix=".cli-test-", dir=ROOT) as temporary_dir:
            base = Path(temporary_dir)
            source_link = base / "source-link.docx"
            arguments = argparse.Namespace(
                source=source_link,
                output_dir=base / "artifacts",
                run_summary=base / "summary.json",
                expected_source_sha256=SOURCE_SHA256,
                document_title_culture_volume_l=Decimal("1"),
                generated_at="2026-08-03T12:01:38+08:00",
                overwrite=False,
            )
            with mocked_symlink_endpoint(source_link, SOURCE), mock.patch.object(
                module, "parse_args", return_value=arguments
            ), mock.patch.object(
                module,
                "parse_expression_docx",
                side_effect=AssertionError("parser must not run for a symlink source"),
            ):
                with self.assertRaisesRegex(FileNotFoundError, "non-symlink file"):
                    module.main()
            self.assertFalse(arguments.output_dir.exists())
            self.assertFalse(arguments.run_summary.exists())

    def test_prepare_preserves_lexical_output_dir_for_symlink_rejection(self) -> None:
        module = load_script(PREPARE)
        with tempfile.TemporaryDirectory(prefix=".cli-test-", dir=ROOT) as temporary_dir:
            base = Path(temporary_dir)
            output_link = base / "artifacts-link"
            resolved_output = base / "real-artifacts"
            resolved_output.mkdir()
            arguments = argparse.Namespace(
                source=SOURCE,
                output_dir=output_link,
                run_summary=base / "summary.json",
                expected_source_sha256=SOURCE_SHA256,
                document_title_culture_volume_l=Decimal("1"),
                generated_at="2026-08-03T12:01:38+08:00",
                overwrite=False,
            )
            with mocked_symlink_endpoint(
                output_link, resolved_output
            ), mock.patch.object(
                module, "parse_args", return_value=arguments
            ), mock.patch.object(
                module,
                "parse_expression_docx",
                side_effect=AssertionError("parser must not run for a symlink output dir"),
            ):
                with self.assertRaisesRegex(ValueError, "Output directory"):
                    module.main()
            self.assertEqual(list(resolved_output.iterdir()), [])
            self.assertFalse(arguments.run_summary.exists())

    def test_plot_preserves_lexical_output_path_for_symlink_rejection(self) -> None:
        module = load_script(PLOT)
        with tempfile.TemporaryDirectory(prefix=".cli-test-", dir=ROOT) as temporary_dir:
            base = Path(temporary_dir)
            plot_data = base / "qc.csv"
            plot_data.write_text(
                "metric,category,count\nsource_records,LTT,23\n", encoding="utf-8"
            )
            link = base / "plot-link.svg"
            resolved_target = base / "real-plot.svg"
            arguments = argparse.Namespace(
                plot_data=plot_data,
                output=link,
                source_sha256=SOURCE_SHA256,
                overwrite=True,
            )
            with mocked_symlink_endpoint(link, resolved_target), mock.patch.object(
                module.argparse.ArgumentParser, "parse_args", return_value=arguments
            ), mock.patch.object(
                module,
                "render_qc_svg",
                side_effect=AssertionError("renderer must not run for a symlink target"),
            ):
                with self.assertRaisesRegex(ValueError, "symbolic link"):
                    module.main()
            self.assertFalse(resolved_target.exists())

    def test_verifier_preserves_lexical_report_path_for_symlink_rejection(self) -> None:
        module = load_script(VERIFY)
        with tempfile.TemporaryDirectory(prefix=".cli-test-", dir=ROOT) as temporary_dir:
            base = Path(temporary_dir)
            artifact_dir = base / "artifacts"
            artifact_dir.mkdir()
            link = base / "report-link.json"
            resolved_target = base / "real-report.json"
            arguments = argparse.Namespace(
                source=SOURCE,
                artifact_dir=artifact_dir,
                report=link,
                expected_source_sha256=SOURCE_SHA256,
                document_title_culture_volume_l=Decimal("1"),
                overwrite=True,
            )
            with mocked_symlink_endpoint(link, resolved_target), mock.patch.object(
                module, "arguments", return_value=arguments
            ), mock.patch.object(
                module,
                "source_records",
                side_effect=AssertionError("verifier must not parse for a symlink target"),
            ):
                with self.assertRaisesRegex(ValueError, "symbolic link"):
                    module.main()
            self.assertFalse(resolved_target.exists())

    def test_fixed_timestamp_runs_are_byte_reproducible_and_verifiable(self) -> None:
        with tempfile.TemporaryDirectory(prefix=".cli-test-", dir=ROOT) as temporary_dir:
            base = Path(temporary_dir)
            outputs = [base / "run1" / "artifacts", base / "run2" / "artifacts"]
            summaries = [base / "run1" / "summary.json", base / "run2" / "summary.json"]
            for output, summary in zip(outputs, summaries, strict=True):
                result = self.run_command(*self.prepare_arguments(output, summary))
                self.assertEqual(result.returncode, 0, result.stderr)
            for name in OUTPUT_NAMES:
                self.assertEqual(digest(outputs[0] / name), digest(outputs[1] / name), name)

            independent_report = base / "independent_validation.json"
            verify = self.run_command(
                str(VERIFY),
                "--source",
                str(SOURCE),
                "--artifact-dir",
                str(outputs[0]),
                "--report",
                str(independent_report),
                "--expected-source-sha256",
                SOURCE_SHA256,
                "--document-title-culture-volume-l",
                "1",
            )
            self.assertEqual(verify.returncode, 0, verify.stderr)
            self.assertTrue(independent_report.is_file())

    def test_implicit_timestamp_is_injected_into_replay_argv(self) -> None:
        with tempfile.TemporaryDirectory(prefix=".cli-test-", dir=ROOT) as temporary_dir:
            base = Path(temporary_dir)
            output = base / "artifacts"
            summary_path = base / "summary.json"
            arguments = self.prepare_arguments(output, summary_path)
            timestamp_index = arguments.index("--generated-at")
            del arguments[timestamp_index : timestamp_index + 2]
            arguments.append("--overwrite")

            first_run = self.run_command(*arguments)
            self.assertEqual(first_run.returncode, 0, first_run.stderr)
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            invocation = summary["invocation"]
            recorded_argv = invocation["argv"]
            self.assertEqual(recorded_argv.count("--generated-at"), 1)
            recorded_timestamp_index = recorded_argv.index("--generated-at")
            self.assertEqual(
                recorded_argv[recorded_timestamp_index + 1], summary["generated_at"]
            )
            self.assertIn(
                f"'--generated-at' '{summary['generated_at']}'",
                invocation["replay_command_powershell"],
            )

            hashes_before_replay = {
                name: digest(output / name) for name in OUTPUT_NAMES
            }
            replay = subprocess.run(
                [
                    invocation["python_executable"],
                    invocation["script_path"],
                    *recorded_argv,
                ],
                cwd=invocation["working_directory"],
                capture_output=True,
                text=True,
                encoding="utf-8",
                env={**os.environ, "PYTHONUTF8": "1"},
                check=False,
            )
            self.assertEqual(replay.returncode, 0, replay.stderr)
            self.assertEqual(
                {name: digest(output / name) for name in OUTPUT_NAMES},
                hashes_before_replay,
            )


if __name__ == "__main__":
    unittest.main()
