"""TNP developability metrics versus collaborator-reported BL21 yield.

The module normalizes the six official TNP outputs, preserves the existing
47-sample phenotype semantics, and performs low-capacity association tests.
It does not treat TNP metrics as measured expression or train a production
yield model.
"""

from __future__ import annotations

import math
from typing import Mapping, Sequence

import numpy as np

from .nanobert_yield import (
    _kendall,
    _optional_float,
    _spearman,
    partial_spearman_provider_length,
    sequence_features,
    stratified_spearman,
)
from .netsolp_yield import (
    build_netsolp_validation_inputs,
)


TNP_FEATURES = (
    "tnp_total_cdr_length",
    "tnp_cdr3_length",
    "tnp_cdr3_compactness",
    "tnp_psh",
    "tnp_ppc",
    "tnp_pnc",
)
PRIMARY_FEATURE = "tnp_psh"
BOOTSTRAP_SEED = 252041
RESAMPLING_REPLICATES = 5000
TNP_NOT_APPLICABLE = {
    "WCC__4-1": "tnp_anarci_rejected_incomplete_heavy_domain",
    "WCC__4-28": "tnp_anarci_rejected_unrecognized_antibody_domain",
    "WCC__4-11": "tnp_anarci_rejected_project_anarcii_light_chain",
    "WCC__4-42": "nanobodybuilder2_rejected_too_many_missing_residues",
}
ELIGIBLE_COUNT = 43
ELIGIBLE_NUMERIC_COUNT = 27


class TNPYieldError(ValueError):
    """Raised when TNP plan, score, or phenotype evidence violates the contract."""


def verify_immune_builder_refine_patch(source_text: str) -> None:
    """Require the two intended OpenMM CPU calls and reject the buggy set literal."""

    broken = "platform, {'Threads', str(n_threads)})"
    corrected = "platform, {'Threads': str(n_threads)})"
    if broken in source_text or source_text.count(corrected) != 2:
        raise TNPYieldError("ImmuneBuilder OpenMM Threads mapping patch is absent or unexpected")


def build_tnp_validation_inputs(
    expression_rows: Sequence[Mapping[str, object]],
    numbering_rows: Sequence[Mapping[str, object]],
    position_rows: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """Return the frozen phenotype rows with the existing 90% sequence clusters."""

    result = build_netsolp_validation_inputs(expression_rows, numbering_rows, position_rows)
    samples = [dict(row) for row in result["sample_rows"]]
    for row in samples:
        uid = str(row["sample_uid"])
        reason = TNP_NOT_APPLICABLE.get(uid, "")
        row["tnp_applicability"] = "not_applicable" if reason else "eligible"
        row["tnp_inapplicability_reason"] = reason
        row["tnp_scoring_requested"] = not bool(reason)
    if sum(row["tnp_applicability"] == "eligible" for row in samples) != ELIGIBLE_COUNT:
        raise TNPYieldError("Unexpected TNP-eligible sample count")
    return {"sample_rows": samples}


def normalize_tnp_result(
    sample: Mapping[str, object],
    result: Mapping[str, object],
    *,
    modelled_sequence: str,
    elapsed_seconds: float,
) -> dict[str, object]:
    """Normalize one official TNP JSON record and document any terminal trimming."""

    uid = str(sample["sample_uid"])
    sequence = str(sample["sequence_raw"])
    if len(result) != 1:
        raise TNPYieldError(f"TNP result for {uid} must contain exactly one record")
    record = next(iter(result.values()))
    if not isinstance(record, Mapping) or str(record.get("name")) != uid:
        raise TNPYieldError(f"TNP result identity mismatch for {uid}")
    start = sequence.find(modelled_sequence)
    if not modelled_sequence or start < 0 or sequence.find(modelled_sequence, start + 1) >= 0:
        raise TNPYieldError(f"TNP modelled sequence is not a unique subsequence for {uid}")
    flags = record.get("Flags")
    if not isinstance(flags, Mapping):
        raise TNPYieldError(f"TNP flags missing for {uid}")
    values = {
        "tnp_total_cdr_length": _finite(record.get("Total CDR Length")),
        "tnp_cdr3_length": _finite(record.get("CDR3 Length")),
        "tnp_cdr3_compactness": _finite(record.get("CDR3 Compactness")),
        "tnp_psh": _finite(record.get("PSH")),
        "tnp_ppc": _finite(record.get("PPC")),
        "tnp_pnc": _finite(record.get("PNC")),
    }
    return {
        "sample_uid": uid,
        "sequence_raw": sequence,
        "input_length_aa": len(sequence),
        "modelled_sequence": modelled_sequence,
        "modelled_length_aa": len(modelled_sequence),
        "trimmed_n_terminal": sequence[:start],
        "trimmed_c_terminal": sequence[start + len(modelled_sequence) :],
        **values,
        "tnp_flag_total_cdr_length": _flag(flags, "L"),
        "tnp_flag_cdr3_length": _flag(flags, "L3"),
        "tnp_flag_cdr3_compactness": _flag(flags, "C"),
        "tnp_flag_psh": _flag(flags, "PSH"),
        "tnp_flag_ppc": _flag(flags, "PPC"),
        "tnp_flag_pnc": _flag(flags, "PNC"),
        "elapsed_seconds": round(float(elapsed_seconds), 6),
        "scoring_status": "pass",
        "failure_reason": "",
    }


def failed_tnp_result(sample: Mapping[str, object], reason: str, elapsed_seconds: float) -> dict[str, object]:
    """Return a schema-compatible explicit failure row for one real sample."""

    return {
        "sample_uid": sample["sample_uid"],
        "sequence_raw": sample["sequence_raw"],
        "input_length_aa": len(str(sample["sequence_raw"])),
        "modelled_sequence": "",
        "modelled_length_aa": "",
        "trimmed_n_terminal": "",
        "trimmed_c_terminal": "",
        **{feature: "" for feature in TNP_FEATURES},
        **{f"tnp_flag_{name}": "" for name in ("total_cdr_length", "cdr3_length", "cdr3_compactness", "psh", "ppc", "pnc")},
        "elapsed_seconds": round(float(elapsed_seconds), 6),
        "scoring_status": "failed",
        "failure_reason": reason,
    }


def not_applicable_tnp_result(sample: Mapping[str, object]) -> dict[str, object]:
    """Preserve one explicit unscored row for a documented out-of-domain sequence."""

    row = failed_tnp_result(sample, str(sample["tnp_inapplicability_reason"]), 0.0)
    row["scoring_status"] = "not_applicable"
    return row


def analyze_tnp_associations(
    sample_rows: Sequence[Mapping[str, object]],
    score_rows: Sequence[Mapping[str, object]],
    netsolp_rows: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """Analyze TNP metrics and compare PSH increment with the existing NetSolP U signal."""

    if len(sample_rows) != 47 or len(score_rows) != 47:
        raise TNPYieldError("TNP analysis requires 47 planned samples and 47 status rows")
    scores = {str(row["sample_uid"]): row for row in score_rows}
    netsolp = {str(row["sample_uid"]): row for row in netsolp_rows}
    expected = {str(row["sample_uid"]) for row in sample_rows}
    if set(scores) != expected or set(netsolp) != expected:
        raise TNPYieldError("TNP/NetSolP identities do not match the validation plan")
    eligible = [row for row in sample_rows if str(row["tnp_applicability"]) == "eligible"]
    not_applicable = [row for row in sample_rows if str(row["tnp_applicability"]) == "not_applicable"]
    if len(eligible) != ELIGIBLE_COUNT or {str(row["sample_uid"]) for row in not_applicable} != set(TNP_NOT_APPLICABLE):
        raise TNPYieldError("TNP V2 applicability identities do not match the fixed contract")
    for sample in sample_rows:
        expected_status = "not_applicable" if sample["tnp_applicability"] == "not_applicable" else "pass"
        observed_status = str(scores[str(sample["sample_uid"])]["scoring_status"])
        if observed_status != expected_status:
            raise TNPYieldError(f'TNP V2 status mismatch for {sample["sample_uid"]}: {observed_status}')
    passed = [row for row in eligible if str(scores[str(row["sample_uid"])]["scoring_status"]) == "pass"]
    numeric_passed = [row for row in passed if row["observation_semantics"] == "individual_approximate"]
    providers = {str(row["provider_code"]) for row in numeric_passed}
    if len(passed) != ELIGIBLE_COUNT or len(numeric_passed) != ELIGIBLE_NUMERIC_COUNT or not {"LTT", "WCC"}.issubset(providers):
        raise TNPYieldError("TNP V2 coverage gate requires 43/43 eligible and 27/27 numeric samples")

    combined: list[dict[str, object]] = []
    for sample in sample_rows:
        uid = str(sample["sample_uid"])
        row = dict(sample)
        row.update(sequence_features(str(sample["sequence_raw"])))
        score = scores[uid]
        row["tnp_scoring_status"] = score["scoring_status"]
        row["tnp_failure_reason"] = score.get("failure_reason", "")
        for field in TNP_FEATURES:
            row[field] = _optional_float(score.get(field, ""))
        for field in (
            "tnp_flag_total_cdr_length", "tnp_flag_cdr3_length", "tnp_flag_cdr3_compactness",
            "tnp_flag_psh", "tnp_flag_ppc", "tnp_flag_pnc",
        ):
            row[field] = score.get(field, "")
        row["predicted_usability"] = _optional_float(netsolp[uid].get("predicted_usability"))
        combined.append(row)

    numeric = [row for row in combined if row[PRIMARY_FEATURE] is not None and row["observation_semantics"] == "individual_approximate"]
    llj = [row for row in combined if row[PRIMARY_FEATURE] is not None and row["provider_code"] == "LLJ"]
    metrics = [_metric_row(numeric, llj, feature) for feature in TNP_FEATURES]
    _add_bh_fdr(metrics)
    primary = next(row for row in metrics if row["feature"] == PRIMARY_FEATURE)
    low, high = _stratified_bootstrap_ci(numeric)
    primary["bootstrap_95ci_low"] = low
    primary["bootstrap_95ci_high"] = high
    primary["stratified_permutation_p"] = _stratified_permutation_p(numeric, float(primary["stratified_spearman_rho"]))
    cv_rows = _comparison_cv(numeric)
    primary["loocv_increment_over_provider"] = next(row["increment_over_provider"] for row in cv_rows if row["model"] == "TNP PSH")
    primary["cluster_cv_increment_over_provider"] = next(row["cluster_increment_over_provider"] for row in cv_rows if row["model"] == "TNP PSH")
    primary_only_fields = (
        "bootstrap_95ci_low",
        "bootstrap_95ci_high",
        "stratified_permutation_p",
        "loocv_increment_over_provider",
        "cluster_cv_increment_over_provider",
    )
    for row in metrics:
        for field in primary_only_fields:
            row.setdefault(field, "")
    level, reasons = _classify_psh(primary)
    return {
        "sample_rows": combined,
        "metric_rows": metrics,
        "cv_rows": cv_rows,
        "primary": primary,
        "evidence_level": level,
        "decision_reasons": reasons,
        "coverage": {
            "planned": 47,
            "eligible": len(eligible),
            "not_applicable": len(not_applicable),
            "eligible_passed": len(passed),
            "numeric_eligible_passed": len(numeric),
            "llj_eligible_passed": len(llj),
        },
    }


def _metric_row(numeric: Sequence[Mapping[str, object]], llj: Sequence[Mapping[str, object]], feature: str) -> dict[str, object]:
    usable = [row for row in numeric if row[feature] is not None]
    llj_usable = [row for row in llj if row[feature] is not None]
    values = np.asarray([float(row[feature]) for row in usable])
    outcomes = np.asarray([float(row["numeric_yield_value"]) for row in usable])
    providers = [str(row["provider_code"]) for row in usable]
    constant = len(set(values.tolist())) < 2
    pooled, pooled_p = (0.0, 1.0) if constant else _spearman(values, outcomes)
    source = {
        p: 0.0 if len(set(values[np.asarray(providers) == p].tolist())) < 2 else _spearman(values[np.asarray(providers) == p], outcomes[np.asarray(providers) == p])[0]
        for p in sorted(set(providers))
    }
    llj_values = np.asarray([float(row[feature]) for row in llj_usable])
    return {
        "feature": feature,
        "numeric_n": len(usable),
        "pooled_spearman_rho": pooled,
        "pooled_spearman_p": pooled_p,
        "stratified_spearman_rho": 0.0 if constant else stratified_spearman(values, outcomes, providers),
        "length_adjusted_partial_spearman": 0.0 if constant else partial_spearman_provider_length(
            values, outcomes, np.asarray([float(row["sequence_length_aa"]) for row in usable]), providers
        ),
        "ltt_spearman_rho": source.get("LTT"),
        "wcc_spearman_rho": source.get("WCC"),
        "llj_ordinal_n": len(llj_usable),
        "llj_kendall_tau_b": 0.0 if len(set(llj_values.tolist())) < 2 else _kendall(
            llj_values,
            np.asarray([int(row["llj_ordinal_level"]) for row in llj_usable]),
        )[0],
    }


def _comparison_cv(rows: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    outcomes = np.asarray([float(row["numeric_yield_value"]) for row in rows])
    providers = [str(row["provider_code"]) for row in rows]
    groups = [str(row["sequence_cluster_90"]) for row in rows]
    features = {
        "Provider only": [],
        "NetSolP U": ["predicted_usability"],
        "TNP PSH": [PRIMARY_FEATURE],
        "NetSolP U + TNP PSH": ["predicted_usability", PRIMARY_FEATURE],
    }
    ordinary_base = _multi_feature_cv(rows, outcomes, providers, groups, [], grouped=False)
    cluster_base = _multi_feature_cv(rows, outcomes, providers, groups, [], grouped=True)
    result = []
    for name, fields in features.items():
        ordinary = _multi_feature_cv(rows, outcomes, providers, groups, fields, grouped=False)
        clustered = _multi_feature_cv(rows, outcomes, providers, groups, fields, grouped=True)
        result.append({
            "model": name,
            "loocv_spearman": ordinary,
            "increment_over_provider": ordinary - ordinary_base,
            "cluster_cv_spearman": clustered,
            "cluster_increment_over_provider": clustered - cluster_base,
        })
    return result


def _multi_feature_cv(rows, outcomes, providers, groups, fields, *, grouped: bool) -> float:
    provider_array = np.asarray(providers)
    group_array = np.asarray(groups)
    folds = sorted(set(groups)) if grouped else list(range(len(rows)))
    predictions = np.empty(len(rows), dtype=float)
    for fold in folds:
        held = group_array == fold if grouped else np.arange(len(rows)) == fold
        train = ~held
        train_columns = [np.ones(int(train.sum()))]
        held_columns = [np.ones(int(held.sum()))]
        for provider in sorted(set(providers))[1:]:
            train_columns.append((provider_array[train] == provider).astype(float))
            held_columns.append((provider_array[held] == provider).astype(float))
        for field in fields:
            values = np.asarray([float(row[field]) for row in rows])
            mean, sd = float(values[train].mean()), float(values[train].std()) or 1.0
            train_columns.append((values[train] - mean) / sd)
            held_columns.append((values[held] - mean) / sd)
        coefficients = np.linalg.lstsq(np.column_stack(train_columns), np.log1p(outcomes[train]), rcond=None)[0]
        predictions[held] = np.column_stack(held_columns) @ coefficients
    return _spearman(predictions, outcomes)[0]


def _stratified_bootstrap_ci(rows):
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    grouped = {p: [row for row in rows if row["provider_code"] == p] for p in ("LTT", "WCC")}
    estimates = []
    for _ in range(RESAMPLING_REPLICATES):
        sample = []
        for members in grouped.values():
            sample.extend(members[i] for i in rng.integers(0, len(members), len(members)))
        estimates.append(stratified_spearman(
            np.asarray([float(row[PRIMARY_FEATURE]) for row in sample]),
            np.asarray([float(row["numeric_yield_value"]) for row in sample]),
            [str(row["provider_code"]) for row in sample],
        ))
    return tuple(float(value) for value in np.quantile(estimates, [0.025, 0.975]))


def _stratified_permutation_p(rows, observed):
    rng = np.random.default_rng(BOOTSTRAP_SEED + 1)
    values = np.asarray([float(row[PRIMARY_FEATURE]) for row in rows])
    outcomes = np.asarray([float(row["numeric_yield_value"]) for row in rows])
    providers = [str(row["provider_code"]) for row in rows]
    masks = [np.asarray(providers) == provider for provider in sorted(set(providers))]
    exceed = 0
    for _ in range(RESAMPLING_REPLICATES):
        permuted = outcomes.copy()
        for mask in masks:
            permuted[mask] = rng.permutation(permuted[mask])
        exceed += abs(stratified_spearman(values, permuted, providers)) >= abs(observed)
    return (exceed + 1) / (RESAMPLING_REPLICATES + 1)


def _classify_psh(metric):
    strict = (
        float(metric["stratified_spearman_rho"]) <= -0.30
        and float(metric["bootstrap_95ci_high"]) < 0
        and float(metric["stratified_permutation_p"]) <= 0.05
        and float(metric["ltt_spearman_rho"]) < 0
        and float(metric["wcc_spearman_rho"]) < 0
        and float(metric["length_adjusted_partial_spearman"]) < 0
        and float(metric["cluster_cv_increment_over_provider"]) > 0
    )
    if strict:
        return "weak_ranking_evidence", ["all_predeclared_psh_direction_uncertainty_cross_provider_and_cluster_cv_criteria_passed"]
    negative = sum(float(metric[field]) < 0 for field in (
        "stratified_spearman_rho", "length_adjusted_partial_spearman", "ltt_spearman_rho", "wcc_spearman_rho"
    ))
    if negative >= 3 and float(metric["cluster_cv_increment_over_provider"]) > 0:
        return "compatibility_filter_only", ["psh_has_partial_expected_direction_but_full_weak_ranking_gate_not_met"]
    return "no_supported_yield_use", ["psh_direction_cross_provider_uncertainty_or_cluster_cv_not_supported"]


def _add_bh_fdr(rows):
    order = sorted(range(len(rows)), key=lambda index: float(rows[index]["pooled_spearman_p"]))
    adjusted = [1.0] * len(rows)
    running = 1.0
    for reverse_rank, index in enumerate(reversed(order), 1):
        rank = len(rows) - reverse_rank + 1
        running = min(running, float(rows[index]["pooled_spearman_p"]) * len(rows) / rank)
        adjusted[index] = min(1.0, running)
    for row, value in zip(rows, adjusted, strict=True):
        row["pooled_spearman_bh_fdr"] = value


def _finite(value: object) -> float:
    number = _optional_float(value)
    if number is None or not math.isfinite(number):
        raise TNPYieldError("TNP metric must be finite")
    return number


def _flag(flags: Mapping[str, object], key: str) -> str:
    value = str(flags.get(key, "")).lower()
    if value not in {"green", "amber", "red"}:
        raise TNPYieldError(f"Invalid TNP flag for {key}")
    return value
