from __future__ import annotations

import csv
import json
from pathlib import Path

from antibody_optimization.vhh_conservation import (
    NumberedVhh,
    build_expression_constraints,
    classify_nb252_positions,
    cluster_and_weight,
    load_tnp_paper_sequences,
)


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
