from __future__ import annotations

import csv
import json
from pathlib import Path

from antibody_optimization.nanobert_yield import (
    build_validation_inputs,
    classify_primary_evidence,
    sequence_features,
)


ROOT = Path(__file__).resolve().parents[1]


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def test_real_47_sequence_validation_plan_preserves_semantics() -> None:
    result = build_validation_inputs(
        _csv(ROOT / "docs/result_artifacts/nb_expression/nb_expression_records.csv"),
        _csv(ROOT / "docs/result_artifacts/input_baseline/sequence/sequence_numbering_review.csv"),
        _csv(ROOT / "docs/result_artifacts/input_baseline/sequence/sequence_numbering_positions.csv"),
    )
    samples = result["sample_rows"]
    assert len(samples) == 47
    assert sum(row["observation_semantics"] == "individual_approximate" for row in samples) == 31
    assert sum(row["provider_code"] == "LLJ" and row["numeric_yield_value"] == "" for row in samples) == 16
    failed = next(row for row in samples if row["sample_uid"] == "WCC__4-28")
    assert failed["numbering_status"] == "failed"
    assert failed["numbered_start_0based"] == ""


def test_physicochemical_features_are_finite_and_predeclared() -> None:
    features = sequence_features("QVQLVESGGGLVQAGGSLRLSCAAS")
    assert set(features) == {
        "sequence_length_aa", "molecular_weight_da", "theoretical_pi", "charge_at_ph7_4",
        "gravy", "aromaticity", "instability_index", "hydrophobic_fraction", "positive_fraction",
        "negative_fraction", "nxs_t_motif_count", "deamidation_motif_count", "oxidation_residue_count",
    }
    assert features["sequence_length_aa"] == 25.0


def test_evidence_gate_requires_direction_uncertainty_and_loocv() -> None:
    passing = {
        "stratified_spearman_rho": 0.4,
        "bootstrap_95ci_low": 0.05,
        "stratified_permutation_p": 0.02,
        "ltt_spearman_rho": 0.3,
        "wcc_spearman_rho": 0.4,
        "length_adjusted_partial_spearman": 0.35,
        "loocv_increment_over_provider": 0.12,
    }
    assert classify_primary_evidence(passing)[0] == "weak_ranking_evidence"
    passing["bootstrap_95ci_low"] = -0.1
    assert classify_primary_evidence(passing)[0] == "compatibility_filter_only"
    passing["stratified_spearman_rho"] = -0.2
    assert classify_primary_evidence(passing)[0] == "no_supported_use"


def test_allowed_use_manifest_releases_only_exploratory_pooling() -> None:
    manifest = json.loads(
        (ROOT / "docs/result_artifacts/input_baseline/expression/allowed_use_manifest.json").read_text(encoding="utf-8")
    )
    assert manifest["gates"]["cross_assay_pooling_gate"] == "pass"
    assert manifest["gates"]["nb252_transfer_gate"] == "blocked"
    assert "pooled_continuous_yield_model" in manifest["blocked_uses"]
