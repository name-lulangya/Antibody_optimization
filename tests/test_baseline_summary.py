from __future__ import annotations

import sys
import tempfile
import unittest
import hashlib
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from antibody_optimization.baseline_summary import (  # noqa: E402
    NB252_UID,
    BaselineSummaryError,
    build_plot_rows,
    build_status_counts,
)
from antibody_optimization.baseline_plot import render_baseline_figure  # noqa: E402
from antibody_optimization.stage_gates import evaluate_stage1_gates  # noqa: E402


NB252_SEQUENCE = (
    "QVQLQESGGGLVQAGGSLRLSCAASGTIFFGYDMGWYRQAPGKEREFVASITTGSNTNYADSVKGRF"
    "TISRDNAKNTVYLQMNSLKPEDTAVYYCAVDTIDYIIEWNVYYYIFSYWGQGTQVTVSSGS"
)


class BaselineSummaryTests(unittest.TestCase):
    def test_unconfirmed_structure_stays_unavailable_and_blocks_release(self) -> None:
        expression = [{"sample_uid": NB252_UID, "sequence_raw": NB252_SEQUENCE}]
        review = [
            {
                "sample_uid": NB252_UID,
                "numbering_status": "pass",
                "query_start_0based_inclusive": "0",
                "query_end_0based_inclusive": "125",
                "provisional_numbered_span_sequence": NB252_SEQUENCE[:126],
                "unnumbered_n_sequence": "",
                "unnumbered_c_sequence": "GS",
            }
        ]
        positions = [
            {
                "sample_uid": NB252_UID,
                "sequence_index_1based": str(index),
                "residue": residue,
                "imgt_position": str(min(index, 128)),
                "imgt_insertion_code": "",
                "imgt_region": "FR1" if index <= 26 else "FR4",
            }
            for index, residue in enumerate(NB252_SEQUENCE[:126], start=1)
        ]
        rows = build_plot_rows(
            expression_records=expression,
            numbering_review=review,
            numbering_positions=positions,
        )
        self.assertEqual(len(rows), 128)
        self.assertEqual(rows[-2]["imgt_region"], "UNNUMBERED")
        self.assertEqual(rows[-1]["imgt_region"], "UNNUMBERED")
        self.assertTrue(
            all(row["experimental_coordinate_status"] == "not_available" for row in rows)
        )

    def test_status_counts_and_gate_do_not_require_all_sequences_to_number(self) -> None:
        numbering = [
            {
                "sample_uid": f"sample_{index}",
                "numbering_status": "failed" if index == 46 else "pass",
            }
            for index in range(47)
        ]
        comparability = [
            {
                "sample_uid": f"sample_{index}",
                "cross_assay_pooling_status": "blocked",
                "nb252_transfer_status": "blocked",
                "highest_allowed_use": "sequence_descriptive_only",
            }
            for index in range(47)
        ]
        counts = build_status_counts(
            numbering_review=numbering,
            sample_comparability=comparability,
        )
        count_map = {(row["metric"], row["category"]): row["count"] for row in counts}
        self.assertEqual(count_map[("numbering_status", "pass")], 46)
        self.assertEqual(count_map[("numbering_status", "failed")], 1)

        gate = evaluate_stage1_gates(
            numbering_review=numbering,
            sample_comparability=comparability,
            input_freeze_manifest={"status": "pass"},
            structure_manifest=None,
            interface_manifest=None,
        )
        self.assertEqual(gate["gates"]["sequence_numbering_inventory"]["status"], "pass")
        self.assertEqual(gate["gates"]["expression_audit"]["status"], "pass")
        self.assertEqual(gate["local_baseline_build"], "blocked")
        self.assertEqual(gate["candidate_design_release"], "blocked")
        self.assertEqual(gate["pooled_expression_model_release"], "blocked")

    def test_numbering_residue_mismatch_fails(self) -> None:
        with self.assertRaisesRegex(BaselineSummaryError, "residue mismatch"):
            build_plot_rows(
                expression_records=[
                    {"sample_uid": NB252_UID, "sequence_raw": NB252_SEQUENCE}
                ],
                numbering_review=[
                    {
                        "sample_uid": NB252_UID,
                        "numbering_status": "pass",
                        "query_start_0based_inclusive": "0",
                        "query_end_0based_inclusive": "0",
                        "provisional_numbered_span_sequence": "Q",
                        "unnumbered_n_sequence": "",
                        "unnumbered_c_sequence": NB252_SEQUENCE[1:],
                    }
                ],
                numbering_positions=[
                    {
                        "sample_uid": NB252_UID,
                        "sequence_index_1based": "1",
                        "residue": "A",
                    }
                ],
            )

    def test_candidate_release_cannot_bypass_local_structure_gates(self) -> None:
        numbering = [
            {"sample_uid": f"sample_{index}", "numbering_status": "pass"}
            for index in range(47)
        ]
        comparability = [
            {
                "sample_uid": f"sample_{index}",
                "cross_assay_pooling_status": "blocked",
                "nb252_transfer_status": "blocked",
            }
            for index in range(47)
        ]
        gate = evaluate_stage1_gates(
            numbering_review=numbering,
            sample_comparability=comparability,
            input_freeze_manifest={"status": "pass"},
            structure_manifest={
                "export_status": "blocked",
                "inventory_status": "pass",
                "chain_role_status": "pass",
                "residue_mapping_status": "pass",
                "authoritative_sequence_status": "pass",
            },
            interface_manifest={
                "interface_status": "pass",
                "orange_annotation_status": "pass",
            },
        )
        self.assertEqual(gate["local_baseline_build"], "blocked")
        self.assertEqual(gate["candidate_design_release"], "blocked")
        self.assertEqual(
            gate["candidate_design_release_blockers"], ["structure_export"]
        )

    def test_render_is_byte_stable_with_fixed_timestamp(self) -> None:
        expression = [{"sample_uid": NB252_UID, "sequence_raw": NB252_SEQUENCE}]
        review = [
            {
                "sample_uid": NB252_UID,
                "numbering_status": "pass",
                "query_start_0based_inclusive": "0",
                "query_end_0based_inclusive": "125",
                "provisional_numbered_span_sequence": NB252_SEQUENCE[:126],
                "unnumbered_n_sequence": "",
                "unnumbered_c_sequence": "GS",
            }
        ]
        positions = [
            {
                "sample_uid": NB252_UID,
                "sequence_index_1based": str(index),
                "residue_aa": residue,
                "numbering_position": str(index),
                "region": "FR1" if index <= 26 else "FR4",
            }
            for index, residue in enumerate(NB252_SEQUENCE[:126], start=1)
        ]
        rows = build_plot_rows(
            expression_records=expression,
            numbering_review=review,
            numbering_positions=positions,
        )
        counts = [
            {"metric": "numbering_status", "category": "pass", "count": 46},
            {"metric": "numbering_status", "category": "failed", "count": 1},
        ]
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            hashes = []
            for suffix in ("first", "second"):
                png = root / f"{suffix}.png"
                svg = root / f"{suffix}.svg"
                render_baseline_figure(
                    plot_rows=rows,
                    status_counts=counts,
                    png_path=png,
                    svg_path=svg,
                    generated_at="2026-08-03T20:00:00+08:00",
                )
                svg_bytes = svg.read_bytes()
                self.assertNotIn(b"\r", svg_bytes)
                self.assertTrue(
                    all(line == line.rstrip() for line in svg_bytes.splitlines())
                )
                hashes.append(
                    (
                        hashlib.sha256(png.read_bytes()).hexdigest(),
                        hashlib.sha256(svg_bytes).hexdigest(),
                    )
                )
        self.assertEqual(hashes[0], hashes[1])


if __name__ == "__main__":
    unittest.main()
