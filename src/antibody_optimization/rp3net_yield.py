"""RP3Net score normalization and BL21 reported-yield validation."""

from __future__ import annotations

from typing import Mapping, Sequence

from .netsolp_yield import build_netsolp_validation_inputs, yield_metric_row
from .nanobert_yield import (
    classify_primary_evidence,
    sequence_features,
    stratified_bootstrap_ci,
    stratified_permutation_p,
)
from .yield_classification import nested_yield_classification


PRIMARY_FEATURE = "rp3net_expression_probability"
EXPECTED_SAMPLE_COUNT = 47


class RP3NetYieldError(ValueError):
    """Raised when RP3Net inputs or outputs violate the frozen contract."""


def build_rp3net_validation_inputs(
    expression_rows: Sequence[Mapping[str, object]],
    numbering_rows: Sequence[Mapping[str, object]],
    position_rows: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """Build the same 47-sample phenotype and 90%-identity cluster plan used by NetSolP."""

    return build_netsolp_validation_inputs(expression_rows, numbering_rows, position_rows)


def normalize_rp3net_scores(
    samples: Sequence[Mapping[str, object]], raw_rows: Sequence[Mapping[str, object]]
) -> list[dict[str, object]]:
    """Validate official ``id,score`` output and preserve validation-plan order."""

    if len(samples) != EXPECTED_SAMPLE_COUNT or len(raw_rows) != EXPECTED_SAMPLE_COUNT:
        raise RP3NetYieldError("Expected 47 plan rows and 47 RP3Net output rows")
    if not raw_rows or not {"id", "score"}.issubset(raw_rows[0]):
        raise RP3NetYieldError("RP3Net output must contain id and score columns")
    by_id = {str(row["id"]): row for row in raw_rows}
    expected = {str(row["sample_uid"]) for row in samples}
    if len(by_id) != EXPECTED_SAMPLE_COUNT or set(by_id) != expected:
        raise RP3NetYieldError("RP3Net output IDs do not match the validation plan")
    normalized = []
    for sample in samples:
        uid = str(sample["sample_uid"])
        score = float(by_id[uid]["score"])
        if not 0 <= score <= 1:
            raise RP3NetYieldError(f"RP3Net score is outside [0, 1] for {uid}")
        normalized.append(
            {
                "sample_uid": uid,
                "sequence_raw": str(sample["sequence_raw"]),
                PRIMARY_FEATURE: score,
                "scoring_status": "pass",
            }
        )
    return normalized


def analyze_rp3net_associations(
    sample_rows: Sequence[Mapping[str, object]],
    score_rows: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """Evaluate RP3Net continuously and with leakage-controlled binary validation."""

    if len(sample_rows) != EXPECTED_SAMPLE_COUNT or len(score_rows) != EXPECTED_SAMPLE_COUNT:
        raise RP3NetYieldError("RP3Net validation requires exactly 47 samples and scores")
    scores = {str(row["sample_uid"]): row for row in score_rows}
    expected = {str(row["sample_uid"]) for row in sample_rows}
    if len(scores) != EXPECTED_SAMPLE_COUNT or set(scores) != expected:
        raise RP3NetYieldError("RP3Net score identities do not match the plan")
    combined = []
    for sample in sample_rows:
        uid = str(sample["sample_uid"])
        score = scores[uid]
        if score.get("scoring_status") != "pass" or score.get("sequence_raw") != sample.get("sequence_raw"):
            raise RP3NetYieldError(f"RP3Net score provenance mismatch for {uid}")
        row = dict(sample)
        row.update(sequence_features(str(sample["sequence_raw"])))
        row[PRIMARY_FEATURE] = float(score[PRIMARY_FEATURE])
        combined.append(row)

    numeric = [row for row in combined if row["observation_semantics"] == "individual_approximate"]
    llj = [row for row in combined if row["provider_code"] == "LLJ"]
    primary = yield_metric_row(numeric, llj, PRIMARY_FEATURE)
    low, high = stratified_bootstrap_ci(numeric, PRIMARY_FEATURE)
    primary["bootstrap_95ci_low"] = low
    primary["bootstrap_95ci_high"] = high
    primary["stratified_permutation_p"] = stratified_permutation_p(
        numeric, PRIMARY_FEATURE, float(primary["stratified_spearman_rho"])
    )
    continuous_level, continuous_reasons = classify_primary_evidence(primary)
    loo = nested_yield_classification(numeric, PRIMARY_FEATURE, outer_scheme="leave_one_out")
    cluster = nested_yield_classification(numeric, PRIMARY_FEATURE, outer_scheme="leave_one_cluster_out")
    evidence_level, reasons = _combined_evidence(continuous_level, loo["summary"], cluster["summary"])
    prediction_rows = loo["prediction_rows"] + cluster["prediction_rows"]
    classification_rows = [
        {"outer_scheme": "leave_one_out", **loo["summary"]},
        {"outer_scheme": "leave_one_cluster_out", **cluster["summary"]},
    ]
    primary["continuous_evidence_level"] = continuous_level
    return {
        "sample_rows": combined,
        "metric_rows": [primary],
        "classification_rows": classification_rows,
        "classification_prediction_rows": prediction_rows,
        "primary": primary,
        "evidence_level": evidence_level,
        "decision_reasons": continuous_reasons + reasons,
    }


def _combined_evidence(
    continuous_level: str,
    loo: Mapping[str, object],
    cluster: Mapping[str, object],
) -> tuple[str, list[str]]:
    loo_pass = (
        float(loo["roc_auc"]) >= 0.65
        and float(loo["pr_auc_average_precision"]) >= float(loo["prevalence"]) + 0.10
        and float(loo["mcc"]) >= 0.25
        and float(loo["balanced_accuracy"]) >= 0.60
        and min(float(loo["sensitivity"]), float(loo["specificity"])) >= 0.50
    )
    cluster_support = float(cluster["mcc"]) > 0 and float(cluster["balanced_accuracy"]) > 0.50
    threshold_stable = float(loo["score_threshold_q3"]) - float(loo["score_threshold_q1"]) <= 0.25
    if continuous_level == "weak_ranking_evidence" and loo_pass and cluster_support and threshold_stable:
        return "weak_ranking_evidence", ["continuous_and_nested_classification_gates_passed"]
    if continuous_level != "no_supported_use" and loo_pass:
        return "compatibility_filter_only", ["classification_supported_but_full_cross_validation_gate_not_met"]
    return "no_supported_use", ["independent_continuous_and_classification_support_not_established"]
