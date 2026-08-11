from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from antibody_optimization.affinity_candidates import (  # noqa: E402
    AffinityCandidateError,
    build_affinity_candidates,
    load_affinity_candidate_inputs,
    validate_affinity_candidates,
)


STAGE0 = PROJECT_ROOT / "docs/result_artifacts/candidate_design/stage0_contract_20260810"
CALIBRATION = (
    PROJECT_ROOT
    / "docs/result_artifacts/structure_preparation/pyrosetta_scoring_calibration_v2_20260811"
)
SCRIPT = PROJECT_ROOT / "scripts/candidate_design/build_affinity_single_mutants.py"


class AffinityCandidateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.inputs = load_affinity_candidate_inputs(
            project_root=PROJECT_ROOT,
            stage0_dir=STAGE0,
            calibration_dir=CALIBRATION,
        )
        cls.rows, cls.positions, cls.gate = build_affinity_candidates(cls.inputs)

    def test_real_inputs_generate_complete_reversible_space(self) -> None:
        self.assertEqual(len(self.rows), 456)
        self.assertEqual(len(self.positions), 24)
        self.assertEqual(Counter(row["sequence_index_1based"] for row in self.rows), Counter({position: 19 for position in self.inputs["interface_positions"]}))
        self.assertEqual(self.gate["candidate_manifest_release"], "pass")
        self.assertEqual(self.gate["selected_scoring_protocol"], "interface_repack_constrained_min")
        self.assertEqual({len(row["candidate_sequence"]) for row in self.rows}, {128})
        self.assertEqual({row["sequence_difference_count"] for row in self.rows}, {1})

    def test_numbering_and_sensitive_positions_are_explicit(self) -> None:
        sensitive = {
            int(row["sequence_index_1based"])
            for row in self.rows
            if row["prepared_contact_sensitive"]
        }
        self.assertEqual(sensitive, {46, 101, 103})
        row = next(
            row
            for row in self.rows
            if row["sequence_index_1based"] == 103 and row["mutant_residue"] == "N"
        )
        self.assertEqual(row["numbering_position_label"], "111A")
        self.assertEqual(row["experimental_auth_seq_id"], "103")
        self.assertEqual(row["prepared_wt_contact_status"], "lost")
        self.assertIn("reported_seq", row["mutation_reported_label"])
        self.assertIn("IMGT", row["mutation_numbering_label"])
        self.assertIn("chain C", row["mutation_source_auth_label"])

    def test_pilot_set_is_stratified_and_contains_sensitive_positions(self) -> None:
        pilot = [row for row in self.rows if row["pilot_selected"]]
        self.assertEqual(len(pilot), 12)
        self.assertEqual(
            {int(row["sequence_index_1based"]) for row in pilot} & {46, 101, 103},
            {46, 101, 103},
        )
        self.assertGreaterEqual(len({row["region"] for row in pilot}), 4)
        self.assertGreaterEqual(len({row["mutant_residue_class"] for row in pilot}), 4)
        self.assertTrue(all(row["pilot_selection_reason"] for row in pilot))

    def test_candidate_validation_rejects_multi_mutant_sequence(self) -> None:
        broken = [dict(row) for row in self.rows]
        sequence = list(str(broken[0]["candidate_sequence"]))
        extra = 0 if int(broken[0]["sequence_index_1based"]) != 1 else 1
        sequence[extra] = "A" if sequence[extra] != "A" else "V"
        broken[0]["candidate_sequence"] = "".join(sequence)
        with self.assertRaises(AffinityCandidateError):
            validate_affinity_candidates(
                broken,
                parent_sequence=str(self.inputs["parent_sequence"]),
                positions=self.inputs["interface_positions"],
            )

    def test_cli_outputs_are_deterministic_bom_lf_and_non_overwriting(self) -> None:
        with tempfile.TemporaryDirectory(dir=PROJECT_ROOT) as temporary:
            root = Path(temporary)
            output = root / "artifacts"
            summary = root / "summary.json"
            command = [
                sys.executable,
                str(SCRIPT),
                "--stage0-dir",
                str(STAGE0),
                "--calibration-dir",
                str(CALIBRATION),
                "--output-dir",
                str(output),
                "--run-summary",
                str(summary),
                "--generated-at",
                "2026-08-11T12:00:00+08:00",
            ]
            completed = subprocess.run(command, capture_output=True, text=True, check=False)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            expected = {
                "affinity_single_mutants.csv",
                "affinity_single_mutants.fasta",
                "affinity_position_summary.csv",
                "pilot_candidate_ids.txt",
                "affinity_candidate_gate.json",
                "affinity_candidate_space_qc.png",
                "affinity_candidate_space_qc.svg",
                "affinity_candidate_manifest.json",
            }
            self.assertEqual({path.name for path in output.iterdir()}, expected)
            data = (output / "affinity_single_mutants.csv").read_bytes()
            self.assertTrue(data.startswith(b"\xef\xbb\xbf"))
            self.assertNotIn(b"\r\n", data)
            with (output / "affinity_single_mutants.csv").open(
                "r", encoding="utf-8-sig", newline=""
            ) as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 456)
            self.assertEqual(sum(row["pilot_selected"] == "True" for row in rows), 12)
            manifest = json.loads(
                (output / "affinity_candidate_manifest.json").read_text(encoding="utf-8")
            )
            for record in manifest["outputs"].values():
                path = output / record["file"]
                self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), record["sha256"])

            output_two = root / "artifacts_two"
            summary_two = root / "summary_two.json"
            second = list(command)
            second[second.index("--output-dir") + 1] = str(output_two)
            second[second.index("--run-summary") + 1] = str(summary_two)
            completed_two = subprocess.run(second, capture_output=True, text=True, check=False)
            self.assertEqual(completed_two.returncode, 0, completed_two.stderr)
            for name in expected:
                self.assertEqual(
                    hashlib.sha256((output / name).read_bytes()).hexdigest(),
                    hashlib.sha256((output_two / name).read_bytes()).hexdigest(),
                    name,
                )
            repeated = subprocess.run(command, capture_output=True, text=True, check=False)
            self.assertNotEqual(repeated.returncode, 0)
            self.assertIn("Refusing to overwrite", repeated.stderr)


if __name__ == "__main__":
    unittest.main()
