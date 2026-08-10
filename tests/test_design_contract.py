from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from antibody_optimization.design_contract import (  # noqa: E402
    DesignContractError,
    build_stage0_contract,
    resolve_and_validate_critical_bindings,
)


CRITICAL = PROJECT_ROOT / "docs/result_artifacts/input_baseline/reviews/nb252_critical_residue_sets.json"
GATE = PROJECT_ROOT / "docs/result_artifacts/input_baseline/summary/stage1_gate.json"
STRUCTURE = PROJECT_ROOT / "data/structures/cxs_exports/NK2R-252__native.cif"
SCRIPT = PROJECT_ROOT / "scripts/candidate_design/build_stage2_design_contract.py"


class DesignContractTests(unittest.TestCase):
    def test_real_stage0_contract_reproduces_critical_sets(self) -> None:
        contract, rows, preflight = build_stage0_contract(
            project_root=PROJECT_ROOT,
            critical_facts_path=CRITICAL,
            stage1_gate_path=GATE,
            experimental_structure_path=STRUCTURE,
            generated_at="2026-08-10T16:00:00+08:00",
        )
        self.assertEqual(len(rows), 128)
        self.assertEqual(contract["counts"]["experimental_missing_positions"], 13)
        self.assertEqual(contract["counts"]["experimental_interface_positions"], 24)
        self.assertEqual(contract["counts"]["hard_immutable_positions"], 6)
        self.assertEqual(contract["counts"]["first_round_affinity_allowed_positions"], 24)
        self.assertEqual(
            contract["hard_immutable"]["coordinate_supported_disulfide_indices_1based"],
            [22, 95],
        )
        self.assertEqual(
            contract["hard_immutable"]["terminal_SSGS_indices_1based"],
            [125, 126, 127, 128],
        )
        self.assertFalse(
            contract["affinity_structure_policy"]["bulk_completion_required_before_first_round"]
        )
        self.assertEqual(preflight["stage0_local_contract"], "pass")
        self.assertEqual(
            preflight["pyrosetta_affinity_scoring_release"],
            "blocked_pending_remote_gap_safe_import",
        )

    def test_stale_critical_binding_is_rejected(self) -> None:
        critical = json.loads(CRITICAL.read_text(encoding="utf-8"))
        critical["source_bindings"]["interface_manifest"]["sha256"] = "0" * 64
        with self.assertRaises(DesignContractError):
            resolve_and_validate_critical_bindings(
                project_root=PROJECT_ROOT, critical_facts=critical
            )

    def test_real_cli_outputs_are_parseable_bom_lf_and_non_overwriting(self) -> None:
        with tempfile.TemporaryDirectory(dir=PROJECT_ROOT) as temporary:
            root = Path(temporary)
            output = root / "artifacts"
            summary = root / "summary.json"
            command = [
                sys.executable,
                str(SCRIPT),
                "--critical-residue-facts",
                str(CRITICAL),
                "--stage1-gate",
                str(GATE),
                "--experimental-structure",
                str(STRUCTURE),
                "--output-dir",
                str(output),
                "--run-summary",
                str(summary),
                "--generated-at",
                "2026-08-10T16:00:00+08:00",
            ]
            completed = subprocess.run(command, capture_output=True, text=True, check=False)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            expected = {
                "stage2_design_contract.json",
                "mutable_position_inventory.csv",
                "stage2_preflight.json",
                "stage2_design_contract_qc.png",
                "stage2_design_contract_qc.svg",
                "stage2_stage0_manifest.json",
            }
            self.assertEqual({path.name for path in output.iterdir()}, expected)
            data = (output / "mutable_position_inventory.csv").read_bytes()
            self.assertTrue(data.startswith(b"\xef\xbb\xbf"))
            self.assertNotIn(b"\r\n", data)
            with (output / "mutable_position_inventory.csv").open(
                encoding="utf-8-sig", newline=""
            ) as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 128)
            manifest = json.loads(
                (output / "stage2_stage0_manifest.json").read_text(encoding="utf-8")
            )
            for key, record in manifest["outputs"].items():
                path = output / record["file"]
                self.assertEqual(
                    hashlib.sha256(path.read_bytes()).hexdigest(), record["sha256"]
                )
            second_output = root / "artifacts_second"
            second_summary = root / "summary_second.json"
            second_command = list(command)
            second_command[second_command.index("--output-dir") + 1] = str(second_output)
            second_command[second_command.index("--run-summary") + 1] = str(second_summary)
            second = subprocess.run(
                [str(item) for item in second_command],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(second.returncode, 0, second.stderr)
            for name in expected:
                self.assertEqual(
                    hashlib.sha256((output / name).read_bytes()).hexdigest(),
                    hashlib.sha256((second_output / name).read_bytes()).hexdigest(),
                    name,
                )
            repeated = subprocess.run(command, capture_output=True, text=True, check=False)
            self.assertNotEqual(repeated.returncode, 0)
            self.assertIn("Refusing to overwrite", repeated.stderr)


if __name__ == "__main__":
    unittest.main()
