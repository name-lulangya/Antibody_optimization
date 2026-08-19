"""PLM_Sol score normalization and BL21 reported-yield validation."""

from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np


PRIMARY_FEATURE = "plm_sol_solubility_score"
EXPECTED_SAMPLE_COUNT = 47


class PLMSolYieldError(ValueError):
    """Raised when PLM_Sol evidence violates the fixed validation contract."""


def normalize_plm_sol_scores(
    samples: Sequence[Mapping[str, object]], raw_rows: Sequence[Mapping[str, object]]
) -> list[dict[str, object]]:
    """Map official classifier rows to samples by the exact, unique input sequence.

    PLM_Sol returns the bio-embeddings hash as ``protein_ID``.  The complete
    sequence is therefore the reversible identity key for this fixed panel;
    duplicate input sequences are rejected rather than assigned arbitrarily.
    """

    if len(samples) != EXPECTED_SAMPLE_COUNT or len(raw_rows) != EXPECTED_SAMPLE_COUNT:
        raise PLMSolYieldError("Expected 47 plan rows and 47 PLM_Sol output rows")
    required = {"protein_ID", "sequence", "predict_result"}
    if not raw_rows or not required.issubset(raw_rows[0]):
        raise PLMSolYieldError("PLM_Sol output lacks required columns")
    planned = {str(row["sequence_raw"]): row for row in samples}
    if len(planned) != EXPECTED_SAMPLE_COUNT:
        raise PLMSolYieldError("PLM_Sol validation requires 47 unique input sequences")
    observed = {str(row["sequence"]): row for row in raw_rows}
    if len(observed) != EXPECTED_SAMPLE_COUNT or set(observed) != set(planned):
        raise PLMSolYieldError("PLM_Sol output sequences do not match the validation plan")
    normalized = []
    for sample in samples:
        sequence = str(sample["sequence_raw"])
        raw = observed[sequence]
        score = float(raw["predict_result"])
        if not np.isfinite(score) or not 0 <= score <= 1:
            raise PLMSolYieldError("PLM_Sol scores must be finite values in [0, 1]")
        normalized.append(
            {
                "sample_uid": str(sample["sample_uid"]),
                "sequence_raw": sequence,
                "embedding_key": str(raw["protein_ID"]),
                "sequence_length_aa": len(sequence),
                PRIMARY_FEATURE: score,
                "scoring_status": "pass",
            }
        )
    return normalized


def analyze_plm_sol_associations(
    sample_rows: Sequence[Mapping[str, object]],
    score_rows: Sequence[Mapping[str, object]],
    netsolp_rows: Sequence[Mapping[str, object]],
    rp3net_rows: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """Evaluate PLM_Sol continuously, discretely, and against fixed comparators."""

    from .nanobert_yield import (
        _spearman,
        classify_primary_evidence,
        sequence_features,
        stratified_bootstrap_ci,
        stratified_permutation_p,
    )
    from .netsolp_yield import yield_metric_row
    from .yield_classification import (
        fixed_yield_apparent_classification,
        fixed_yield_nested_classification,
        nested_yield_classification,
    )

    combined = _join_scores(sample_rows, score_rows)
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

    classifications = []
    predictions = []
    for scheme in ("leave_one_out", "leave_one_cluster_out"):
        result = nested_yield_classification(numeric, PRIMARY_FEATURE, outer_scheme=scheme)
        classifications.append({"outer_scheme": scheme, **result["summary"]})
        predictions.extend(result["prediction_rows"])
    loo = classifications[0]
    cluster = classifications[1]
    evidence_level, reasons = _combined_evidence(continuous_level, loo, cluster)

    comparison_rows, independent_increment = _compare_predictors(
        combined, netsolp_rows, rp3net_rows, _spearman
    )
    if evidence_level == "weak_ranking_evidence" and independent_increment <= 0:
        evidence_level = "compatibility_filter_only"
        reasons = ["full_single_predictor_gate_passed_but_no_cluster_cv_increment_over_netsolp_s"]

    fixed_metrics = []
    fixed_predictions = []
    for scheme in ("leave_one_out", "leave_one_cluster_out"):
        result = fixed_yield_nested_classification(
            numeric, PRIMARY_FEATURE, outer_scheme=scheme, yield_threshold=5.0
        )
        fixed_metrics.append({"outer_scheme": scheme, **result["summary"]})
        fixed_predictions.extend(result["prediction_rows"])
    apparent = fixed_yield_apparent_classification(
        numeric, PRIMARY_FEATURE, yield_threshold=5.0
    )
    fixed_metrics.append({"outer_scheme": "apparent_full_sample", **apparent["summary"]})
    fixed_predictions.extend(apparent["prediction_rows"])
    primary["continuous_evidence_level"] = continuous_level
    return {
        "sample_rows": combined,
        "metric_rows": [primary],
        "classification_rows": classifications,
        "classification_prediction_rows": predictions,
        "comparison_rows": comparison_rows,
        "fixed5_metric_rows": fixed_metrics,
        "fixed5_prediction_rows": fixed_predictions,
        "primary": primary,
        "evidence_level": evidence_level,
        "decision_reasons": continuous_reasons + reasons,
        "independent_cluster_cv_increment": independent_increment,
    }


def _join_scores(
    sample_rows: Sequence[Mapping[str, object]], score_rows: Sequence[Mapping[str, object]]
) -> list[dict[str, object]]:
    from .nanobert_yield import sequence_features

    if len(sample_rows) != EXPECTED_SAMPLE_COUNT or len(score_rows) != EXPECTED_SAMPLE_COUNT:
        raise PLMSolYieldError("PLM_Sol validation requires exactly 47 samples and scores")
    scores = {str(row["sample_uid"]): row for row in score_rows}
    if len(scores) != EXPECTED_SAMPLE_COUNT or set(scores) != {
        str(row["sample_uid"]) for row in sample_rows
    }:
        raise PLMSolYieldError("PLM_Sol score identities do not match the plan")
    combined = []
    for sample in sample_rows:
        uid = str(sample["sample_uid"])
        score = scores[uid]
        if score.get("scoring_status") != "pass" or score.get("sequence_raw") != sample.get("sequence_raw"):
            raise PLMSolYieldError(f"PLM_Sol score provenance mismatch for {uid}")
        row = dict(sample)
        row.update(sequence_features(str(sample["sequence_raw"])))
        row[PRIMARY_FEATURE] = float(score[PRIMARY_FEATURE])
        row["embedding_key"] = str(score["embedding_key"])
        combined.append(row)
    return combined


def _compare_predictors(combined, netsolp_rows, rp3net_rows, spearman):
    net = {str(row["sample_uid"]): row for row in netsolp_rows}
    rp3 = {str(row["sample_uid"]): row for row in rp3net_rows}
    identifiers = {str(row["sample_uid"]) for row in combined}
    if len(net) != 47 or len(rp3) != 47 or set(net) != identifiers or set(rp3) != identifiers:
        raise PLMSolYieldError("Comparator inputs must contain the same 47 samples")
    joined = []
    for row in combined:
        uid = str(row["sample_uid"])
        if net[uid]["sequence_raw"] != row["sequence_raw"] or rp3[uid]["sequence_raw"] != row["sequence_raw"]:
            raise PLMSolYieldError(f"Comparator sequence mismatch for {uid}")
        item = dict(row)
        item["predicted_solubility"] = float(net[uid]["predicted_solubility"])
        item["predicted_usability"] = float(net[uid]["predicted_usability"])
        item["rp3net_expression_probability"] = float(rp3[uid]["rp3net_expression_probability"])
        joined.append(item)
    numeric = [row for row in joined if row["observation_semantics"] == "individual_approximate"]
    plm_all = np.asarray([float(row[PRIMARY_FEATURE]) for row in joined])
    plm_numeric = np.asarray([float(row[PRIMARY_FEATURE]) for row in numeric])
    rows = []
    for feature, display in (
        ("predicted_solubility", "NetSolP S"),
        ("predicted_usability", "NetSolP U"),
        ("rp3net_expression_probability", "RP3Net"),
    ):
        rows.append(
            {
                "comparison": display,
                "all_47_spearman": spearman(plm_all, np.asarray([float(row[feature]) for row in joined]))[0],
                "numeric_31_spearman": spearman(plm_numeric, np.asarray([float(row[feature]) for row in numeric]))[0],
                "partial_spearman_provider_length_netsolp_s": "",
                "cluster_cv_spearman": "",
                "increment_over_netsolp_s": "",
            }
        )
    outcomes = np.asarray([float(row["numeric_yield_value"]) for row in numeric])
    providers = [str(row["provider_code"]) for row in numeric]
    groups = [str(row["sequence_cluster_90"]) for row in numeric]
    net_s = np.asarray([float(row["predicted_solubility"]) for row in numeric])
    conditional = _partial_spearman_with_netsolp_s(
        plm_numeric, outcomes,
        np.asarray([float(row["sequence_length_aa"]) for row in numeric]),
        net_s, providers,
    )
    base = _group_cv_multifeature([net_s], outcomes, providers, groups, spearman)
    combined_cv = _group_cv_multifeature([net_s, plm_numeric], outcomes, providers, groups, spearman)
    rows.append(
        {
            "comparison": "NetSolP S + PLM_Sol incremental model",
            "all_47_spearman": "",
            "numeric_31_spearman": "",
            "partial_spearman_provider_length_netsolp_s": conditional,
            "cluster_cv_spearman": combined_cv,
            "increment_over_netsolp_s": combined_cv - base,
        }
    )
    rows.append(
        {
            "comparison": "NetSolP S baseline model",
            "all_47_spearman": "",
            "numeric_31_spearman": "",
            "partial_spearman_provider_length_netsolp_s": "",
            "cluster_cv_spearman": base,
            "increment_over_netsolp_s": 0.0,
        }
    )
    return rows, combined_cv - base


def _partial_spearman_with_netsolp_s(values, outcomes, lengths, netsolp_s, providers):
    from scipy.stats import rankdata

    provider_array = np.asarray(providers)
    columns = [np.ones(len(values)), rankdata(lengths), rankdata(netsolp_s)]
    for provider in sorted(set(providers))[1:]:
        columns.append((provider_array == provider).astype(float))
    design = np.column_stack(columns)
    ranked_values = rankdata(values)
    ranked_outcomes = rankdata(outcomes)
    value_residual = ranked_values - design @ np.linalg.lstsq(design, ranked_values, rcond=None)[0]
    outcome_residual = ranked_outcomes - design @ np.linalg.lstsq(design, ranked_outcomes, rcond=None)[0]
    correlation = float(np.corrcoef(value_residual, outcome_residual)[0, 1])
    return correlation if np.isfinite(correlation) else 0.0


def _group_cv_multifeature(features, outcomes, providers, groups, spearman):
    provider_array = np.asarray(providers)
    group_array = np.asarray(groups)
    provider_levels = sorted(set(providers))
    predictions = np.empty(len(outcomes), dtype=float)
    for group in sorted(set(groups)):
        held = group_array == group
        train = ~held
        columns = [np.ones(int(train.sum()))]
        held_columns = [np.ones(int(held.sum()))]
        for provider in provider_levels[1:]:
            columns.append((provider_array[train] == provider).astype(float))
            held_columns.append((provider_array[held] == provider).astype(float))
        for feature in features:
            mean = float(feature[train].mean())
            sd = float(feature[train].std()) or 1.0
            columns.append((feature[train] - mean) / sd)
            held_columns.append((feature[held] - mean) / sd)
        coefficients = np.linalg.lstsq(
            np.column_stack(columns), np.log1p(outcomes[train]), rcond=None
        )[0]
        predictions[held] = np.column_stack(held_columns) @ coefficients
    return float(spearman(predictions, outcomes)[0])


def _combined_evidence(continuous_level, loo, cluster):
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
