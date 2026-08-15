"""NanoMelt predicted-apparent-Tm association with reported BL21 yield."""

from __future__ import annotations

import math
from typing import Mapping, Sequence

import numpy as np

from .nanobert_yield import (
    _kendall,
    _spearman,
    loocv_spearman,
    partial_spearman_provider_length,
    stratified_spearman,
)
from .netsolp_yield import (
    build_netsolp_validation_inputs,
    group_cv_spearman,
)


PRIMARY_FEATURE = "nanomelt_predicted_apparent_tm_c"
RESAMPLING_REPLICATES = 10_000
BOOTSTRAP_SEED = 252041
NB252_UID = "LTT__Nb252"


class NanoMeltYieldError(ValueError):
    """Raised when NanoMelt inputs or scores violate the validation contract."""


def build_nanomelt_validation_inputs(
    expression_rows: Sequence[Mapping[str, object]],
    numbering_rows: Sequence[Mapping[str, object]],
    position_rows: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """Return the fixed 47 samples with the established 90% identity clusters."""

    return build_netsolp_validation_inputs(expression_rows, numbering_rows, position_rows)


def normalize_nanomelt_scores(
    samples: Sequence[Mapping[str, object]],
    raw_rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    """Validate official NanoMelt CSV and map its scored domain to each raw sequence."""

    if len(samples) != 47 or len(raw_rows) != 47:
        raise NanoMeltYieldError("Expected 47 plan rows and 47 NanoMelt output rows")
    required = {"ID", "Aligned Sequence", "Sequence", "NanoMelt Tm (C)"}
    if not raw_rows or not required.issubset(raw_rows[0]):
        raise NanoMeltYieldError("NanoMelt output lacks required columns")
    by_id = {str(row["ID"]): row for row in raw_rows}
    expected = {str(row["sample_uid"]) for row in samples}
    if len(by_id) != 47 or set(by_id) != expected:
        raise NanoMeltYieldError("NanoMelt output IDs do not match the 47-sample plan")

    normalized = []
    for sample in samples:
        uid = str(sample["sample_uid"])
        raw_sequence = str(sample["sequence_raw"])
        raw = by_id[uid]
        aligned = str(raw["Aligned Sequence"]).strip().upper()
        scored = str(raw["Sequence"]).strip().upper()
        ungapped = "".join(residue for residue in aligned if residue.isalpha())
        if ungapped != scored:
            raise NanoMeltYieldError(f"Aligned/scored sequence mismatch for {uid}")
        starts = [index for index in range(len(raw_sequence) - len(scored) + 1) if raw_sequence.startswith(scored, index)]
        if len(starts) != 1:
            raise NanoMeltYieldError(f"NanoMelt scored domain is not a unique raw-sequence segment for {uid}")
        start = starts[0]
        tm = _finite_float(raw["NanoMelt Tm (C)"], f"NanoMelt Tm for {uid}")
        normalized.append(
            {
                "sample_uid": uid,
                "sequence_raw": raw_sequence,
                "aligned_sequence": aligned,
                "scored_ungapped_sequence": scored,
                "scored_length_aa": len(scored),
                "trimmed_n_terminal": raw_sequence[:start],
                "trimmed_c_terminal": raw_sequence[start + len(scored) :],
                PRIMARY_FEATURE: tm,
                "scoring_status": "pass",
            }
        )
    return normalized


def analyze_nanomelt_associations(
    sample_rows: Sequence[Mapping[str, object]],
    score_rows: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """Evaluate predicted apparent Tm without treating it as measured Tm or yield."""

    if len(sample_rows) != 47 or len(score_rows) != 47:
        raise NanoMeltYieldError("NanoMelt association analysis requires 47 samples and scores")
    scores = {str(row["sample_uid"]): row for row in score_rows}
    if len(scores) != 47 or set(scores) != {str(row["sample_uid"]) for row in sample_rows}:
        raise NanoMeltYieldError("NanoMelt score identities do not match the plan")

    combined = []
    for sample in sample_rows:
        uid = str(sample["sample_uid"])
        score = scores[uid]
        if score.get("scoring_status") != "pass" or score.get("sequence_raw") != sample["sequence_raw"]:
            raise NanoMeltYieldError(f"NanoMelt score identity/status mismatch for {uid}")
        row = dict(sample)
        row.update(score)
        combined.append(row)

    numeric = [row for row in combined if row["observation_semantics"] == "individual_approximate"]
    llj = [row for row in combined if row["provider_code"] == "LLJ"]
    if len(numeric) != 31 or len(llj) != 16:
        raise NanoMeltYieldError("Expected 31 numeric and 16 LLJ ordinal records")
    primary = _primary_metric(numeric, llj)
    low, high = _stratified_bootstrap_ci(numeric)
    primary.update(
        {
            "bootstrap_95ci_low": low,
            "bootstrap_95ci_high": high,
            "stratified_permutation_p": _stratified_permutation_p(
                numeric, float(primary["stratified_spearman_rho"])
            ),
        }
    )
    leave_one_out = _leave_one_out_rhos(numeric)
    primary["leave_one_out_rho_min"] = min(leave_one_out.values())
    primary["leave_one_out_rho_max"] = max(leave_one_out.values())
    primary["without_nb252_stratified_spearman_rho"] = leave_one_out[NB252_UID]
    level, reasons = _classify_evidence(primary)
    primary["evidence_level"] = level
    cv_rows = [
        {
            "model": "provider_only",
            "loocv_spearman": primary["provider_only_loocv_spearman"],
            "leave_cluster_out_spearman": primary["cluster_provider_only_cv_spearman"],
            "increment_over_provider": 0.0,
            "cluster_increment_over_provider": 0.0,
        },
        {
            "model": "provider_plus_nanomelt_tm",
            "loocv_spearman": primary["loocv_log1p_yield_spearman"],
            "leave_cluster_out_spearman": primary["cluster_cv_log1p_yield_spearman"],
            "increment_over_provider": primary["loocv_increment_over_provider"],
            "cluster_increment_over_provider": primary["cluster_cv_increment_over_provider"],
        },
    ]
    return {
        "sample_rows": combined,
        "metric_rows": [primary],
        "cv_rows": cv_rows,
        "leave_one_out_rows": [
            {"held_out_sample_uid": uid, "stratified_spearman_rho": rho}
            for uid, rho in leave_one_out.items()
        ],
        "primary": primary,
        "evidence_level": level,
        "decision_reasons": reasons,
    }


def _primary_metric(numeric, llj) -> dict[str, object]:
    values = np.asarray([float(row[PRIMARY_FEATURE]) for row in numeric])
    outcomes = np.asarray([float(row["numeric_yield_value"]) for row in numeric])
    providers = [str(row["provider_code"]) for row in numeric]
    groups = [str(row["sequence_cluster_90"]) for row in numeric]
    pooled_rho, pooled_p = _spearman(values, outcomes)
    source = {
        provider: _spearman(values[np.asarray(providers) == provider], outcomes[np.asarray(providers) == provider])[0]
        for provider in sorted(set(providers))
    }
    loocv = loocv_spearman(values, outcomes, providers, include_feature=True)
    baseline = loocv_spearman(values, outcomes, providers, include_feature=False)
    cluster_cv = group_cv_spearman(values, outcomes, providers, groups, include_feature=True)
    cluster_baseline = group_cv_spearman(values, outcomes, providers, groups, include_feature=False)
    llj_values = np.asarray([float(row[PRIMARY_FEATURE]) for row in llj])
    llj_levels = np.asarray([int(row["llj_ordinal_level"]) for row in llj])
    return {
        "feature": PRIMARY_FEATURE,
        "numeric_n": len(numeric),
        "sequence_cluster_count_90": len(set(groups)),
        "pooled_spearman_rho": pooled_rho,
        "pooled_spearman_p": pooled_p,
        "stratified_spearman_rho": stratified_spearman(values, outcomes, providers),
        "scored_length_adjusted_partial_spearman": partial_spearman_provider_length(
            values, outcomes, np.asarray([float(row["scored_length_aa"]) for row in numeric]), providers
        ),
        "raw_length_adjusted_partial_spearman": partial_spearman_provider_length(
            values, outcomes, np.asarray([float(row["sequence_length_aa"]) for row in numeric]), providers
        ),
        "ltt_spearman_rho": source["LTT"],
        "wcc_spearman_rho": source["WCC"],
        "loocv_log1p_yield_spearman": loocv,
        "provider_only_loocv_spearman": baseline,
        "loocv_increment_over_provider": loocv - baseline,
        "cluster_cv_log1p_yield_spearman": cluster_cv,
        "cluster_provider_only_cv_spearman": cluster_baseline,
        "cluster_cv_increment_over_provider": cluster_cv - cluster_baseline,
        "llj_ordinal_n": len(llj),
        "llj_kendall_tau_b": _kendall(llj_values, llj_levels)[0],
    }


def _classify_evidence(metric: Mapping[str, object]) -> tuple[str, list[str]]:
    strict = (
        float(metric["stratified_spearman_rho"]) >= 0.30
        and float(metric["bootstrap_95ci_low"]) > 0
        and float(metric["stratified_permutation_p"]) <= 0.05
        and float(metric["ltt_spearman_rho"]) > 0
        and float(metric["wcc_spearman_rho"]) > 0
        and float(metric["scored_length_adjusted_partial_spearman"]) > 0
        and float(metric["cluster_cv_increment_over_provider"]) > 0
        and float(metric["without_nb252_stratified_spearman_rho"]) > 0
    )
    if strict:
        return "weak_ranking_evidence", ["all_predeclared_direction_uncertainty_cluster_cv_and_nb252_influence_criteria_passed"]
    positive = sum(
        float(metric[field]) > 0
        for field in (
            "stratified_spearman_rho",
            "ltt_spearman_rho",
            "wcc_spearman_rho",
            "scored_length_adjusted_partial_spearman",
            "cluster_cv_increment_over_provider",
            "without_nb252_stratified_spearman_rho",
        )
    )
    if float(metric["stratified_spearman_rho"]) > 0 and positive >= 5:
        return "compatibility_filter_only", ["positive_direction_but_full_weak_ranking_gate_not_met"]
    return "no_supported_use", ["predicted_tm_direction_uncertainty_transfer_or_nb252_influence_not_supported"]


def _leave_one_out_rhos(rows) -> dict[str, float]:
    result = {}
    for held in rows:
        kept = [row for row in rows if row["sample_uid"] != held["sample_uid"]]
        result[str(held["sample_uid"])] = stratified_spearman(
            np.asarray([float(row[PRIMARY_FEATURE]) for row in kept]),
            np.asarray([float(row["numeric_yield_value"]) for row in kept]),
            [str(row["provider_code"]) for row in kept],
        )
    return result


def _stratified_bootstrap_ci(rows) -> tuple[float, float]:
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    grouped = {provider: [row for row in rows if row["provider_code"] == provider] for provider in ("LTT", "WCC")}
    estimates = []
    for _ in range(RESAMPLING_REPLICATES):
        sample = []
        for members in grouped.values():
            sample.extend(members[index] for index in rng.integers(0, len(members), len(members)))
        estimates.append(
            stratified_spearman(
                np.asarray([float(row[PRIMARY_FEATURE]) for row in sample]),
                np.asarray([float(row["numeric_yield_value"]) for row in sample]),
                [str(row["provider_code"]) for row in sample],
            )
        )
    return tuple(float(value) for value in np.quantile(estimates, [0.025, 0.975]))


def _stratified_permutation_p(rows, observed: float) -> float:
    rng = np.random.default_rng(BOOTSTRAP_SEED + 1)
    values = np.asarray([float(row[PRIMARY_FEATURE]) for row in rows])
    outcomes = np.asarray([float(row["numeric_yield_value"]) for row in rows])
    providers = np.asarray([str(row["provider_code"]) for row in rows])
    masks = [providers == provider for provider in sorted(set(providers))]
    exceed = 0
    for _ in range(RESAMPLING_REPLICATES):
        permuted = outcomes.copy()
        for mask in masks:
            permuted[mask] = rng.permutation(permuted[mask])
        exceed += abs(stratified_spearman(values, permuted, providers.tolist())) >= abs(observed)
    return (exceed + 1) / (RESAMPLING_REPLICATES + 1)


def _finite_float(value: object, label: str) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise NanoMeltYieldError(f"{label} must be finite")
    return number
