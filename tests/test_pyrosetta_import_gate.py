from __future__ import annotations

import subprocess
import sys
import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from antibody_optimization.pyrosetta_import_gate import (  # noqa: E402
    ResidueRecord,
    build_import_gate,
    compare_pose_to_source,
    evaluate_breaks,
    load_released_stage_inputs,
)


STAGE0 = (
    PROJECT_ROOT
    / "docs/result_artifacts/candidate_design/stage0_contract_20260810"
)
STRUCTURE_BASELINE = (
    PROJECT_ROOT
    / "docs/result_artifacts/input_baseline/structure_released_20260810"
)
SCRIPT = (
    PROJECT_ROOT
    / "scripts/structure_preparation/validate_pyrosetta_wt_import.py"
)
SLURM = (
    PROJECT_ROOT
    / "scripts/structure_preparation/submit_pyrosetta_wt_import.slurm"
)


class PyRosettaImportGateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.inputs = load_released_stage_inputs(
            stage0_dir=STAGE0,
            structure_baseline_dir=STRUCTURE_BASELINE,
        )

    def test_real_released_inputs_define_only_meaningful_breaks(self) -> None:
        residues = self.inputs["source_residues"]
        self.assertEqual(len(residues), 396)
        self.assertEqual(
            [
                (
                    item.break_type,
                    item.left.index,
                    item.right.index,
                    item.left.chain_id,
                    item.left.auth_seq_id,
                    item.right.chain_id,
                    item.right.auth_seq_id,
                )
                for item in self.inputs["expected_breaks"]
            ],
            [
                ("missing_density", 8, 9, "C", 8, "C", 16),
                ("missing_density", 16, 17, "C", 23, "C", 30),
                ("chain_boundary", 115, 116, "C", 128, "R", 30),
                ("missing_density", 314, 315, "R", 228, "R", 240),
            ],
        )

    def test_mapping_and_safe_breaks_pass(self) -> None:
        source = self.inputs["source_residues"]
        pose = [
            ResidueRecord(
                record.index,
                record.chain_id,
                record.auth_seq_id,
                record.insertion_code,
                record.residue_name,
            )
            for record in source
        ]
        self.assertEqual(
            compare_pose_to_source(source_residues=source, pose_residues=pose), []
        )
        breaks, problems = evaluate_breaks(
            expected_breaks=self.inputs["expected_breaks"],
            fold_tree_cutpoints={8, 16, 115, 314},
            jump_cutpoints={8, 16, 115, 314},
            bonded_c_n_pairs=set(),
        )
        self.assertEqual(problems, [])
        self.assertTrue(all(row["status"] == "pass" for row in breaks))
        gate = build_import_gate(
            generated_at="2026-08-10T22:00:00+08:00",
            pyrosetta_version="PyRosetta 2026.03 commit 5e498f",
            score_function="ref2015",
            source_residues=source,
            break_rows=breaks,
            score_rows=[
                {
                    "score_term": "total_score",
                    "raw_value": -1.0,
                    "weight": 1.0,
                    "weighted_value": -1.0,
                }
            ],
            mapping_problems=[],
            break_problems=[],
            disulfide_bonded=True,
            stage0_run_id="stage0",
            structure_run_id="structure",
        )
        self.assertEqual(gate["pyrosetta_wt_import_release"], "pass")
        self.assertEqual(
            gate["pyrosetta_affinity_scoring_release"],
            "ready_for_scoring_protocol_calibration",
        )

    def test_wrong_bond_and_unexpected_cutpoint_block(self) -> None:
        rows, problems = evaluate_breaks(
            expected_breaks=self.inputs["expected_breaks"],
            fold_tree_cutpoints={8, 16, 100, 115, 314},
            jump_cutpoints={8, 16, 115, 314},
            bonded_c_n_pairs={(8, 9)},
        )
        self.assertTrue(any("unexpected FoldTree" in item for item in problems))
        self.assertTrue(any("unsafe missing_density" in item for item in problems))
        self.assertEqual(rows[0]["c_n_atoms_bonded"], True)
        self.assertEqual(rows[0]["status"], "blocked")

    def test_cli_help_does_not_require_local_pyrosetta(self) -> None:
        completed = subprocess.run(
            [sys.executable, str(SCRIPT), "--help"],
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("--stage0-dir", completed.stdout)

    def test_slurm_wrapper_matches_project_resource_policy(self) -> None:
        text = SLURM.read_text(encoding="utf-8")
        self.assertIn("#SBATCH --partition=batch", text)
        self.assertIn("#SBATCH --cpus-per-task=12", text)
        self.assertIn("#SBATCH --gres=gpu:1", text)
        self.assertNotIn("#SBATCH --mem", text)
        self.assertNotIn("#SBATCH --partition=gpu", text)
        self.assertIn("/data/software/env/luly25/multi_ligand", text)
        self.assertIn("validate_pyrosetta_wt_import.py", text)


if __name__ == "__main__":
    unittest.main()
