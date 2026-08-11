from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from antibody_optimization.pyrosetta_scoring_calibration import (  # noqa: E402
    CalibrationThresholds,
    audit_source_incomplete_sidechains,
    build_calibration_gate,
    build_contact_change_rows,
    choose_representative_replicate,
    energy_edge_map,
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
STRUCTURE_BASELINE = (
    PROJECT_ROOT
    / "docs/result_artifacts/input_baseline/structure_released_20260810"
)
V1_CALIBRATION = (
    PROJECT_ROOT
    / "docs/result_artifacts/structure_preparation/pyrosetta_scoring_calibration_20260810"
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
        incomplete = audit_source_incomplete_sidechains(
            structure_baseline_dir=STRUCTURE_BASELINE,
            vhh_interface_auth_positions=inputs["vhh_interface_auth_positions"],
        )
        self.assertEqual(len(incomplete), 20)
        self.assertEqual(
            sum(int(row["missing_heavy_atom_count"]) for row in incomplete),
            90,
        )
        self.assertEqual(
            [
                (row["chain_id"], row["auth_seq_id"], row["residue_name"])
                for row in incomplete
                if row["vhh_experimental_interface"]
            ],
            [("C", 102, "TYR")],
        )

    def test_energy_edge_map_uses_pyrosetta_2026_zero_argument_api(self) -> None:
        class FakeEnergyEdge:
            def __init__(self) -> None:
                self.call_count = 0

            def fill_energy_map(self):
                self.call_count += 1
                return "emap"

        edge = FakeEnergyEdge()
        self.assertEqual(energy_edge_map(edge), "emap")
        self.assertEqual(edge.call_count, 1)

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

    def test_real_v1_metrics_block_repack_and_select_constrained_min(self) -> None:
        with (V1_CALIBRATION / "protocol_replicate_metrics.csv").open(
            "r", encoding="utf-8-sig", newline=""
        ) as handle:
            rows = list(csv.DictReader(handle))
        gate = json.loads(
            (V1_CALIBRATION / "pyrosetta_scoring_calibration_gate.json").read_text(
                encoding="utf-8"
            )
        )
        summaries = summarize_protocol_rows(
            rows,
            raw_interface_fa_rep=float(
                gate["raw_import_metrics"]["interface_fa_rep"]
            ),
            thresholds=CalibrationThresholds(),
        )
        by_protocol = {row["protocol"]: row for row in summaries}
        self.assertEqual(by_protocol["interface_repack"]["status"], "blocked")
        self.assertIn(
            "one_or_more_dg_separated_not_negative",
            by_protocol["interface_repack"]["blockers"],
        )
        self.assertIn(
            "one_or_more_cross_interface_energy_not_negative",
            by_protocol["interface_repack"]["blockers"],
        )
        self.assertEqual(
            by_protocol["interface_repack_constrained_min"]["status"], "pass"
        )
        self.assertEqual(
            select_protocol(summaries)[0],
            "interface_repack_constrained_min",
        )

    def test_contact_change_rows_preserve_exact_position_status(self) -> None:
        rows = build_contact_change_rows(
            molecule_side="Nb252_VHH",
            chain_id="C",
            reference_positions={33, 37},
            prepared_positions={37, 45},
            residue_names={33: "TYR", 37: "SER", 45: "ARG"},
            reference_minimum_distances={33: 3.5, 37: 3.2, 45: 5.0},
            prepared_minimum_distances={33: 4.2, 37: 3.3, 45: 3.8},
        )
        self.assertEqual(
            [(row["auth_seq_id"], row["contact_status"]) for row in rows],
            [(33, "lost"), (37, "retained"), (45, "gained")],
        )
        self.assertEqual(rows[0]["reference_minimum_distance_angstrom"], 3.5)
        self.assertEqual(rows[0]["prepared_minimum_distance_angstrom"], 4.2)

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
            contact_change_rows=[],
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
        self.assertIn("pyrosetta_scoring_calibration_v2_20260811", text)


if __name__ == "__main__":
    unittest.main()
