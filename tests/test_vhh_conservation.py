from __future__ import annotations

import csv
import json
from pathlib import Path

from antibody_optimization.vhh_conservation import (
    NumberedVhh,
    build_project_vhh_records,
    build_expression_constraints,
    classify_nb252_positions,
    cluster_and_weight,
    load_tnp_paper_sequences,
)
from antibody_optimization.vhh_conservation_plot import _region_runs


ROOT = Path(__file__).resolve().parents[1]
TNP_TABLE = (
    ROOT
    / "data/reference/tnp_natural_vhh/VHH_OAS_all_properties_FINAL.csv"
)
CRITICAL = ROOT / "docs/result_artifacts/input_baseline/reviews/nb252_critical_residue_sets.json"
REVIEW = ROOT / "docs/result_artifacts/input_baseline/sequence/sequence_numbering_review.csv"


def _csv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def test_real_tnp_paper_table_is_exact_unique_4059_sequence_set():
    records = load_tnp_paper_sequences(_csv(TNP_TABLE))
    assert len(records) == 4059
    assert len({record.seq_id for record in records}) == 4059
    assert len({record.sequence for record in records}) == 4059
    assert all(record.sequence == record.sequence.strip() for record in records)


def test_imgt_identity_components_receive_total_weight_one():
    records = [
        NumberedVhh("a", "", {str(i): "A" for i in range(1, 11)}),
        NumberedVhh("b", "", {**{str(i): "A" for i in range(1, 10)}, "10": "G"}),
        NumberedVhh("c", "", {str(i): "G" for i in range(1, 11)}),
    ]
    weighted, summary = cluster_and_weight(records, identity_threshold=0.90)
    assert summary["cluster_count"] == 2
    assert weighted[0].cluster_id == weighted[1].cluster_id
    assert weighted[2].cluster_id != weighted[0].cluster_id
    by_cluster = {}
    for record in weighted:
        by_cluster.setdefault(record.cluster_id, 0.0)
        by_cluster[record.cluster_id] += record.weight
    assert all(abs(total - 1.0) < 1e-12 for total in by_cluster.values())


def test_hard_conservation_and_all_hard_constraints_block_mutation():
    reference_rows = [
        {
            "sequence_index_1based": "1",
            "numbering_position_label": "1",
            "region": "FR1",
            "residue_aa": "Q",
            "is_gap": "false",
        }
    ]
    frequency_fields = {f"frequency_{aa}": 0.0 for aa in "ACDEFGHIKLMNPQRSTVWY"}
    global_rows = [
        {
            "imgt_position_label": "1",
            "dominant_aa": "Q",
            "dominant_frequency": 0.92,
            **frequency_fields,
            "frequency_Q": 0.92,
        }
    ]
    neighbor_rows = [
        {
            "imgt_position_label": "1",
            "dominant_aa": "Q",
            "dominant_frequency": 0.95,
            "coverage": 0.99,
            "effective_cluster_count": 60,
            "normalized_conservation": 0.9,
            **frequency_fields,
            "frequency_Q": 0.95,
        }
    ]
    parent = next(row for row in _csv(REVIEW) if row["sample_uid"] == "LTT__Nb252")[
        "sequence_raw"
    ]
    conservation = classify_nb252_positions(parent, reference_rows, global_rows, neighbor_rows)
    assert conservation[0]["conservation_class"] == "hard_conserved"
    assert conservation[126]["conservation_class"] == "insufficient_evidence"

    critical = json.loads(CRITICAL.read_text(encoding="utf-8"))
    contract, positions, candidates = build_expression_constraints(
        parent, critical, conservation
    )
    frozen = set(contract["hard_frozen_reported_indices_1based"])
    interface = set(
        critical["reproduced_experimental_interface"]["reported_sequence_indices_1based"]
    )
    assert {1, 22, 95, 125, 126, 127, 128} | interface <= frozen
    assert {row["reported_sequence_index_1based"] for row in candidates}.isdisjoint(frozen)
    assert all(row["mutant_residue"] != "C" for row in candidates)
    assert all(
        sum(left != right for left, right in zip(parent, row["sequence"], strict=True)) == 1
        for row in candidates
    )
    assert len(positions) == 128


def test_highly_conserved_nonconsensus_position_allows_only_consensus_reversion():
    reference_rows = [
        {
            "sequence_index_1based": "5",
            "numbering_position_label": "5",
            "region": "FR1",
            "residue_aa": "Q",
            "is_gap": "false",
        }
    ]
    frequency_fields = {f"frequency_{aa}": 0.0 for aa in "ACDEFGHIKLMNPQRSTVWY"}
    global_rows = [
        {
            "imgt_position_label": "5",
            "dominant_aa": "V",
            "dominant_frequency": 0.99,
            **frequency_fields,
            "frequency_Q": 0.01,
            "frequency_V": 0.99,
        }
    ]
    neighbor_rows = [
        {
            "imgt_position_label": "5",
            "dominant_aa": "V",
            "dominant_frequency": 0.99,
            "coverage": 0.99,
            "effective_cluster_count": 60,
            "normalized_conservation": 0.9,
            **frequency_fields,
            "frequency_Q": 0.01,
            "frequency_V": 0.99,
        }
    ]
    parent = next(row for row in _csv(REVIEW) if row["sample_uid"] == "LTT__Nb252")[
        "sequence_raw"
    ]
    conservation = classify_nb252_positions(parent, reference_rows, global_rows, neighbor_rows)
    assert conservation[4]["conservation_class"] == "conserved_nonconsensus"
    assert conservation[4]["classification_reason"] == (
        "high_conservation_parent_differs_from_consensus"
    )

    critical = json.loads(CRITICAL.read_text(encoding="utf-8"))
    contract, positions, candidates = build_expression_constraints(
        parent, critical, conservation
    )
    position5 = positions[4]
    assert position5["hard_frozen"] is False
    assert position5["allowed_substitution_rule"] == "natural_consensus_reversion_only"
    assert position5["allowed_mutant_residues"] == "V"
    position5_candidates = [
        row for row in candidates if row["reported_sequence_index_1based"] == 5
    ]
    assert [(row["wt_residue"], row["mutant_residue"]) for row in position5_candidates] == [
        ("Q", "V")
    ]
    assert contract["consensus_reversion_only"] == [
        {
            "reported_sequence_index_1based": 5,
            "wt_residue": "Q",
            "allowed_mutant_residue": "V",
        }
    ]


def test_project_logo_records_exclude_failed_and_non_heavy_sequences():
    review_rows = [
        {"sample_uid": "H1", "numbering_status": "pass", "chain_type": "H", "sequence_raw": "AA"},
        {"sample_uid": "L1", "numbering_status": "pass", "chain_type": "L", "sequence_raw": "CC"},
        {"sample_uid": "F1", "numbering_status": "failed", "chain_type": "F", "sequence_raw": "DD"},
    ]
    position_rows = [
        {"sample_uid": "H1", "is_gap": "false", "numbering_position_label": "1", "residue_aa": "A"},
        {"sample_uid": "L1", "is_gap": "false", "numbering_position_label": "1", "residue_aa": "C"},
    ]
    records, audit = build_project_vhh_records(review_rows, position_rows)
    assert [record.seq_id for record in records] == ["H1"]
    assert [row["logo_status"] for row in audit] == ["included", "excluded", "excluded"]
    assert [row["logo_reason"] for row in audit[1:]] == [
        "non_heavy_chain_assignment",
        "numbering_failed",
    ]


def test_region_runs_preserve_fr_cdr_boundaries():
    rows = [{"region": region} for region in ["FR1", "FR1", "CDR1", "FR2", "FR2"]]
    assert _region_runs(rows) == [(0, 1, "FR1"), (2, 2, "CDR1"), (3, 4, "FR2")]
