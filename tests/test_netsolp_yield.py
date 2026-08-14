from __future__ import annotations

import csv
from pathlib import Path

import pytest

from antibody_optimization.netsolp_yield import (
    NetSolPYieldError,
    analyze_netsolp_associations,
    build_netsolp_validation_inputs,
    normalize_netsolp_scores,
    sequence_identity_clusters,
)


ROOT = Path(__file__).resolve().parents[1]


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _real_samples() -> list[dict[str, object]]:
    return build_netsolp_validation_inputs(
        _csv(ROOT / "docs/result_artifacts/nb_expression/nb_expression_records.csv"),
        _csv(ROOT / "docs/result_artifacts/input_baseline/sequence/sequence_numbering_review.csv"),
        _csv(ROOT / "docs/result_artifacts/input_baseline/sequence/sequence_numbering_positions.csv"),
    )["sample_rows"]


def test_real_47_netsolp_plan_preserves_failed_numbering_and_clusters() -> None:
    samples = _real_samples()
    assert len(samples) == 47
    assert len({row["sample_uid"] for row in samples}) == 47
    assert len({row["sequence_cluster_90"] for row in samples}) == 40
    failed = next(row for row in samples if row["sample_uid"] == "WCC__4-28")
    assert failed["numbering_status"] == "failed"
    assert failed["sequence_cluster_90"]


def test_normalize_netsolp_scores_requires_exact_ids_sequences_and_probabilities() -> None:
    samples = _real_samples()
    raw = [
        {
            "sid": row["sample_uid"],
            "fasta": row["sequence_raw"],
            "predicted_solubility": "0.4",
            "predicted_usability": "0.6",
        }
        for row in samples
    ]
    normalized = normalize_netsolp_scores(samples, raw)
    assert len(normalized) == 47
    assert normalized[0]["predicted_usability"] == 0.6
    raw[0]["predicted_usability"] = "1.1"
    with pytest.raises(NetSolPYieldError, match="probability"):
        normalize_netsolp_scores(samples, raw)


def test_netsolp_association_keeps_u_primary_s_secondary_and_llj_ordinal() -> None:
    samples = _real_samples()
    raw = []
    for index, row in enumerate(samples):
        if row["observation_semantics"] == "individual_approximate":
            usability = min(float(row["numeric_yield_value"]) / 30.0, 0.99)
        else:
            usability = int(row["llj_ordinal_level"]) / 4.0
        raw.append(
            {
                "sid": row["sample_uid"],
                "fasta": row["sequence_raw"],
                "predicted_solubility": 0.25 + index / 200.0,
                "predicted_usability": usability,
            }
        )
    result = analyze_netsolp_associations(samples, normalize_netsolp_scores(samples, raw))
    assert result["primary"]["feature"] == "predicted_usability"
    assert next(row for row in result["metric_rows"] if row["feature"] == "predicted_solubility")
    assert result["primary"]["llj_ordinal_n"] == 16
    assert result["primary"]["sequence_cluster_count_90"] <= 31
    assert result["evidence_level"] in {
        "weak_ranking_evidence",
        "compatibility_filter_only",
        "no_supported_use",
    }


def test_sequence_cluster_definition_is_deterministic_and_single_linkage() -> None:
    clusters = sequence_identity_clusters(
        {"a": "AAAAAAAAAA", "b": "AAAAAAAATA", "c": "AAAAAAATTA", "z": "CCCCCCCCCC"},
        threshold=0.9,
    )
    assert clusters["a"] == clusters["b"] == clusters["c"]
    assert clusters["z"] != clusters["a"]
