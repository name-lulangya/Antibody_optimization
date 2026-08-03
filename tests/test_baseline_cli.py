from __future__ import annotations

import csv
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import unittest
import uuid
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_ROOT / "scripts/input_baseline/finalize_input_baseline.py"
PLOT_SCRIPT = PROJECT_ROOT / "scripts/input_baseline/plot_input_baseline.py"
FIXED_TIME = "2026-08-03T20:00:00+08:00"
SOURCES = {
    "source_cxs": PROJECT_ROOT / "Nb252-optimization.cxs",
    "source_docx": PROJECT_ROOT / "nb序列及产量（1L）.docx",
    "expression_source_manifest": (
        PROJECT_ROOT / "docs/result_artifacts/nb_expression/manifest.json"
    ),
    "expression_records": (
        PROJECT_ROOT / "docs/result_artifacts/nb_expression/nb_expression_records.csv"
    ),
    "sequence_manifest": (
        PROJECT_ROOT
        / "docs/result_artifacts/input_baseline/sequence/sequence_numbering_manifest.json"
    ),
    "numbering_review": (
        PROJECT_ROOT
        / "docs/result_artifacts/input_baseline/sequence/sequence_numbering_review.csv"
    ),
    "numbering_positions": (
        PROJECT_ROOT
        / "docs/result_artifacts/input_baseline/sequence/sequence_numbering_positions.csv"
    ),
    "expression_audit_manifest": (
        PROJECT_ROOT
        / "docs/result_artifacts/input_baseline/expression/allowed_use_manifest.json"
    ),
    "sample_comparability": (
        PROJECT_ROOT
        / "docs/result_artifacts/input_baseline/expression/sample_comparability_review.csv"
    ),
}
EXPECTED_OUTPUTS = {
    "input_freeze_manifest.json",
    "nb252_baseline_plot_data.csv",
    "baseline_status_counts.csv",
    "stage1_gate.json",
    "input_baseline_qc.png",
    "input_baseline_qc.svg",
    "summary_manifest.json",
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class BaselineCliTests(unittest.TestCase):
    def command(self, output: Path, summary: Path) -> list[str]:
        arguments = [sys.executable, str(SCRIPT)]
        for key, path in SOURCES.items():
            arguments.extend(["--" + key.replace("_", "-"), str(path)])
        arguments.extend(
            [
                "--output-dir",
                str(output),
                "--run-summary",
                str(summary),
                "--generated-at",
                FIXED_TIME,
            ]
        )
        return arguments

    @staticmethod
    def execute(command: list[str]) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            command,
            cwd=PROJECT_ROOT,
            capture_output=True,
            text=True,
            encoding="utf-8",
            env={**os.environ, "PYTHONUTF8": "1"},
            check=False,
        )

    def test_real_finalizer_is_transactional_reproducible_and_blocked_honestly(self) -> None:
        source_hashes = {key: digest(path) for key, path in SOURCES.items()}
        with tempfile.TemporaryDirectory(
            prefix=".baseline-cli-test-", dir=PROJECT_ROOT
        ) as temporary_dir:
            root = Path(temporary_dir)
            output = root / "summary"
            run_summary = root / "run_summary.json"
            command = self.command(output, run_summary)

            first = self.execute(command)
            self.assertEqual(first.returncode, 0, first.stderr)
            self.assertEqual({path.name for path in output.iterdir()}, EXPECTED_OUTPUTS)
            self.assertTrue(run_summary.is_file())
            for name in (
                "nb252_baseline_plot_data.csv",
                "baseline_status_counts.csv",
            ):
                self.assertEqual((output / name).read_bytes()[:3], b"\xef\xbb\xbf")

            with (output / "nb252_baseline_plot_data.csv").open(
                "r", encoding="utf-8-sig", newline=""
            ) as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 128)
            self.assertEqual(rows[-2]["residue"] + rows[-1]["residue"], "GS")
            self.assertEqual(rows[-2]["imgt_region"], "UNNUMBERED")

            gate = json.loads((output / "stage1_gate.json").read_text(encoding="utf-8"))
            self.assertEqual(gate["local_baseline_build"], "blocked")
            self.assertEqual(gate["candidate_design_release"], "blocked")
            self.assertEqual(gate["pooled_expression_model_release"], "blocked")
            self.assertIn(
                "structure_export", gate["local_baseline_build_blockers"]
            )
            self.assertIn(
                "authoritative_nb252_sequence",
                gate["candidate_design_release_blockers"],
            )
            self.assertEqual(
                gate["pooled_expression_model_release_blockers"],
                ["cross_assay_pooling"],
            )

            freeze = json.loads(
                (output / "input_freeze_manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(freeze["status"], "pass")
            self.assertEqual(freeze["sequence_identity"]["record_count"], 47)
            artifact_hashes = {
                name: digest(output / name) for name in EXPECTED_OUTPUTS
            }

            rerender_png = root / "rerender.png"
            rerender_svg = root / "rerender.svg"
            rerender_summary = root / "rerender_summary.json"
            plot_command = [
                sys.executable,
                str(PLOT_SCRIPT),
                "--plot-data",
                str(output / "nb252_baseline_plot_data.csv"),
                "--status-counts",
                str(output / "baseline_status_counts.csv"),
                "--png",
                str(rerender_png),
                "--svg",
                str(rerender_svg),
                "--run-summary",
                str(rerender_summary),
                "--generated-at",
                FIXED_TIME,
            ]
            plotted = self.execute(plot_command)
            self.assertEqual(plotted.returncode, 0, plotted.stderr)
            self.assertEqual(digest(rerender_png), digest(output / "input_baseline_qc.png"))
            self.assertEqual(digest(rerender_svg), digest(output / "input_baseline_qc.svg"))
            plot_summary = json.loads(rerender_summary.read_text(encoding="utf-8"))
            self.assertEqual(plot_summary["outputs"]["png"]["sha256"], digest(rerender_png))
            self.assertNotEqual(self.execute(plot_command).returncode, 0)

            refused = self.execute(command)
            self.assertNotEqual(refused.returncode, 0)
            self.assertIn("Refusing to overwrite", refused.stderr)
            self.assertEqual(
                {name: digest(output / name) for name in EXPECTED_OUTPUTS},
                artifact_hashes,
            )

            overwritten = self.execute([*command, "--overwrite"])
            self.assertEqual(overwritten.returncode, 0, overwritten.stderr)
            self.assertEqual(
                {name: digest(output / name) for name in EXPECTED_OUTPUTS},
                artifact_hashes,
            )

        self.assertEqual(
            {key: digest(path) for key, path in SOURCES.items()}, source_hashes
        )

    def test_project_external_targets_are_rejected_before_directory_creation(self) -> None:
        outside = PROJECT_ROOT.parent / f".baseline-outside-{uuid.uuid4().hex}"
        command = self.command(outside / "artifacts", outside / "summary.json")
        completed = self.execute(command)
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("outside the project root", completed.stderr)
        self.assertFalse(outside.exists())

        plot_command = [
            sys.executable,
            str(PLOT_SCRIPT),
            "--plot-data",
            str(
                PROJECT_ROOT
                / "docs/result_artifacts/input_baseline/sequence/sequence_numbering_review.csv"
            ),
            "--status-counts",
            str(
                PROJECT_ROOT
                / "docs/result_artifacts/input_baseline/sequence/sequence_numbering_review.csv"
            ),
            "--png",
            str(outside / "plot.png"),
            "--svg",
            str(outside / "plot.svg"),
            "--run-summary",
            str(outside / "plot.json"),
            "--generated-at",
            FIXED_TIME,
        ]
        plotted = self.execute(plot_command)
        self.assertNotEqual(plotted.returncode, 0)
        self.assertIn("outside the project root", plotted.stderr)
        self.assertFalse(outside.exists())


if __name__ == "__main__":
    unittest.main()
