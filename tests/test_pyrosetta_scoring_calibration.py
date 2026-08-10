from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from antibody_optimization.pyrosetta_scoring_calibration import (  # noqa: E402
    CalibrationThresholds,
    build_calibration_gate,
    choose_representative_replicate,
    load_calibration_inputs,
    render_calibration_svg,
    select_protocol,
    summarize_protocol_rows,
)


STAGE0 = PROJECT_ROOT / "docs/result_artifacts/candidate_design/stage0_contract_20260810"
IMPORT_GATE = (
    PROJECT_ROOT
    / "docs/result_artifacts/structure_preparation/pyrosetta_wt_import_20260810"
)
SCRIPT = PROJECT_ROOT / "scripts/structure_preparation/calibrate_pyrosetta_scoring.py"
SLURM = (
    PROJECT_ROOT
    / "scripts/structure_preparation/submit_pyrosetta_scoring_calibration.slurm"
)


def _row(
    protocol: str,
    replicate: int,
    *,
    dg: float,
    interface_rep: float = 3.0,
    vhh_retention: float = 0.95,
    receptor_retention: float = 0.95,
    rmsd: float = 0.1,
    status: str = "pass",
) -> dict[str, object]:
    return {
        "protocol": protocol,
        "replicate": replicate,
        "seed": 1000 + replicate,
        "total_score": 100.0,
        "dG_separated": dg,
        "cross_interface_energy": -20.0,
        "interface_fa_atr": -30.0,
        "interface_fa_rep": interface_rep,
        "vhh_contact_count": 24,
        "receptor_epitope_count": 20,
        "vhh_contact_retention": vhh_retention,
        "receptor_epitope_retention": receptor_retention,
        "interface_ca_rmsd": rmsd,
        "minimum_interchain_distance": 2.5,
        "mapping_pass": status == "pass",
        "breaks_pass": status == "pass",
        "disulfide_pass": status == "pass",
        "finite_metrics": True,
        "status": status,
    }


class PyRosettaScoringCalibrationTests(unittest.TestCase):
    def test_real_released_inputs_are_bound_to_passed_import_gate(self) -> None:
        inputs = load_calibration_inputs(stage0_dir=STAGE0, import_gate_dir=IMPORT_GATE)
        self.assertEqual(len(inputs["vhh_interface_auth_positions"]), 24)
        self.assertEqual(inputs["vhh_interface_auth_positions"][0], 33)
        self.assertEqual(inputs["vhh_interface_auth_positions"][-1], 116)
        self.assertEqual(inputs["import_gate"]["status"], "pass")

    def test_simplest_passing_protocol_is_selected(self) -> None:
        rows = []
        for protocol in ("interface_repack", "interface_repack_constrained_min"):
            rows.extend(_row(protocol, i, dg=-10.0 + i * 0.2) for i in range(1, 5))
        summaries = summarize_protocol_rows(
            rows,
            raw_interface_fa_rep=10.0,
            thresholds=CalibrationThresholds(),
        )
        selected, blockers = select_protocol(summaries)
        self.assertEqual(selected, "interface_repack")
        self.assertEqual(blockers, [])
        self.assertEqual(
            choose_representative_replicate(rows, protocol=selected),
            2,
        )

    def test_failed_repack_can_select_constrained_minimization(self) -> None:
        rows = [
            _row("interface_repack", i, dg=-10.0, receptor_retention=0.5)
            for i in range(1, 4)
        ]
        rows.extend(
            _row("interface_repack_constrained_min", i, dg=-9.5 + i * 0.1)
            for i in range(1, 4)
        )
        summaries = summarize_protocol_rows(
            rows,
            raw_interface_fa_rep=10.0,
            thresholds=CalibrationThresholds(),
        )
        selected, blockers = select_protocol(summaries)
        self.assertEqual(selected, "interface_repack_constrained_min")
        self.assertEqual(blockers, [])

    def test_both_failed_protocols_keep_affinity_release_blocked(self) -> None:
        rows = []
        for protocol in ("interface_repack", "interface_repack_constrained_min"):
            rows.extend(
                _row(protocol, i, dg=float(i * 10), interface_rep=20.0)
                for i in range(1, 4)
            )
        summaries = summarize_protocol_rows(
            rows,
            raw_interface_fa_rep=10.0,
            thresholds=CalibrationThresholds(),
        )
        selected, blockers = select_protocol(summaries)
        self.assertIsNone(selected)
        self.assertTrue(blockers)
        gate = build_calibration_gate(
            generated_at="2026-08-10T23:00:00+08:00",
            pyrosetta_version="PyRosetta 2026.03",
            score_function="ref2015",
            thresholds=CalibrationThresholds(),
            raw_metrics={"interface_fa_rep": 10.0},
            protocol_summaries=summaries,
            selected_protocol=None,
            representative_replicate=None,
            stage0_run_id="stage0",
            import_gate_run_id="import",
        )
        self.assertEqual(gate["status"], "blocked")
        self.assertEqual(gate["pyrosetta_affinity_scoring_release"], "blocked")
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "gate.svg"
            render_calibration_svg(gate=gate, path=path)
            self.assertIn("BLOCKED", path.read_text(encoding="utf-8"))

    def test_cli_help_and_slurm_contract(self) -> None:
        completed = subprocess.run(
            [sys.executable, str(SCRIPT), "--help"],
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("--replicates", completed.stdout)
        text = SLURM.read_text(encoding="utf-8")
        self.assertIn("#SBATCH --partition=batch", text)
        self.assertIn("#SBATCH --cpus-per-task=12", text)
        self.assertIn("#SBATCH --gres=gpu:1", text)
        self.assertIn("#SBATCH --time=04:00:00", text)
        self.assertNotIn("#SBATCH --mem", text)
        self.assertIn("/data/software/env/luly25/multi_ligand", text)
        self.assertIn("calibrate_pyrosetta_scoring.py", text)


if __name__ == "__main__":
    unittest.main()
