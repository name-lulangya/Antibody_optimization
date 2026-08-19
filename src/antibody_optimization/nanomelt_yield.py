"""NanoMelt predicted-apparent-Tm association with reported BL21 yield."""

from __future__ import annotations

import json
import math
from pathlib import Path
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
from .yield_classification import nested_yield_classification


PRIMARY_FEATURE = "nanomelt_predicted_apparent_tm_c"
RESAMPLING_REPLICATES = 10_000
BOOTSTRAP_SEED = 252041
NB252_UID = "LTT__Nb252"


class NanoMeltYieldError(ValueError):
    """Raised when NanoMelt inputs or scores violate the validation contract."""


def verify_required_openmm_platforms(
    observed: Sequence[str], required: Sequence[str]
) -> list[str]:
    """Require baseline OpenMM platforms while allowing node-specific extras."""

    observed_names = [str(name) for name in observed]
    missing = sorted(set(required) - set(observed_names))
    if missing:
        raise NanoMeltYieldError(f"OpenMM runtime lacks required platforms: {missing}")
    return observed_names


def verify_anarci_runtime(
    module: object,
    environment_root: Path,
    *,
    expected_conda_version: str,
) -> dict[str, object]:
    """Verify the imported ANARCI code and HMM data without trusting stale dist metadata."""

    module_file = Path(str(getattr(module, "__file__", ""))).resolve()
    environment_root = environment_root.resolve()
    if not module_file.is_relative_to(environment_root):
        raise NanoMeltYieldError(f"ANARCI imported outside the active environment: {module_file}")
    required_api = ("anarci", "run_anarci", "validate_sequence", "scheme_short_to_long")
    missing_api = [name for name in required_api if not hasattr(module, name)]
    if missing_api:
        raise NanoMeltYieldError(f"ANARCI runtime lacks required API: {missing_api}")
    hmm_dir = module_file.parent / "dat" / "HMMs"
    required_hmms = [hmm_dir / f"ALL.hmm{suffix}" for suffix in ("", ".h3f", ".h3i", ".h3m", ".h3p")]
    missing_hmms = [str(path) for path in required_hmms if not path.is_file()]
    if missing_hmms:
        raise NanoMeltYieldError(f"ANARCI runtime lacks required HMM data: {missing_hmms}")
    conda_records = [
        json.loads(path.read_text(encoding="utf-8"))
        for path in (environment_root / "conda-meta").glob("anarci-*.json")
    ]
    conda_versions = {str(record.get("version")) for record in conda_records if record.get("name") == "anarci"}
    if conda_versions != {expected_conda_version}:
        raise NanoMeltYieldError(
            f"ANARCI conda package mismatch: {sorted(conda_versions)} != {[expected_conda_version]}"
        )
    return {
        "module_file": str(module_file),
        "conda_package_version": expected_conda_version,
        "required_api": list(required_api),
        "hmm_database": str(required_hmms[0]),
        "hmm_pressed_index_count": 4,
    }


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
    *,
    expected_pass_count: int | None = None,
    expected_plan_count: int = 47,
) -> list[dict[str, object]]:
    """Validate official NanoMelt CSV and map its scored domain to each raw sequence."""

    if expected_plan_count < 1 or len(samples) != expected_plan_count:
        raise NanoMeltYieldError(f"Expected {expected_plan_count} plan rows")
    if expected_pass_count is not None and len(raw_rows) != expected_pass_count:
        raise NanoMeltYieldError(
            f"Expected {expected_pass_count} NanoMelt output rows, found {len(raw_rows)}"
        )
    required = {"ID", "Aligned Sequence", "Sequence", "NanoMelt Tm (C)"}
    if not raw_rows or not required.issubset(raw_rows[0]):
        raise NanoMeltYieldError("NanoMelt output lacks required columns")
    by_id = {str(row["ID"]): row for row in raw_rows}
    expected = {str(row["sample_uid"]) for row in samples}
    if len(by_id) != len(raw_rows) or not set(by_id).issubset(expected):
        raise NanoMeltYieldError("NanoMelt output contains duplicate or unknown IDs")

    normalized = []
    for sample in samples:
        uid = str(sample["sample_uid"])
        raw_sequence = str(sample["sequence_raw"])
        if uid not in by_id:
            normalized.append(
                {
                    "sample_uid": uid,
                    "sequence_raw": raw_sequence,
                    "aligned_sequence": "",
                    "scored_ungapped_sequence": "",
                    "scored_length_aa": "",
                    "trimmed_n_terminal": "",
                    "trimmed_c_terminal": "",
                    PRIMARY_FEATURE: "",
                    "scoring_status": "nanomelt_not_scored",
                    "not_scored_reason": "not_returned_by_nanomelt_anarci_alignment",
                }
            )
            continue
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
                "not_scored_reason": "",
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
        if score.get("scoring_status") not in {"pass", "nanomelt_not_scored"} or score.get("sequence_raw") != sample["sequence_raw"]:
            raise NanoMeltYieldError(f"NanoMelt score identity/status mismatch for {uid}")
        row = dict(sample)
        row.update(score)
        combined.append(row)

    passed = [row for row in combined if row["scoring_status"] == "pass"]
    not_scored = [row for row in combined if row["scoring_status"] == "nanomelt_not_scored"]
    if len(passed) != 43 or len(not_scored) != 4:
        raise NanoMeltYieldError("Expected 43 scored and 4 NanoMelt-not-scored records")
    numeric = [row for row in passed if row["observation_semantics"] == "individual_approximate"]
    llj = [row for row in passed if row["provider_code"] == "LLJ"]
    if NB252_UID not in {str(row["sample_uid"]) for row in numeric} or {row["provider_code"] for row in numeric} != {"LTT", "WCC"}:
        raise NanoMeltYieldError("Scored numeric subset must retain Nb252 and both numeric providers")
    primary = _primary_metric(numeric, llj)
    primary.update(
        {
            "planned_n": 47,
            "scored_n": len(passed),
            "not_scored_count": len(not_scored),
            "not_scored_sample_uids": "|".join(str(row["sample_uid"]) for row in not_scored),
        }
    )
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
    continuous_level, continuous_reasons = _classify_evidence(primary)
    classifications = []
    classification_predictions = []
    for scheme in ("leave_one_out", "leave_one_cluster_out"):
        classification = nested_yield_classification(
            numeric, PRIMARY_FEATURE, outer_scheme=scheme
        )
        classifications.append({"outer_scheme": scheme, **classification["summary"]})
        classification_predictions.extend(classification["prediction_rows"])
    level, classification_reasons = _combined_evidence(
        continuous_level, classifications[0], classifications[1]
    )
    reasons = continuous_reasons + classification_reasons
    reasons.append("association_scope_is_nanomelt_scored_standard_vhh_domains_only")
    primary["continuous_evidence_level"] = continuous_level
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
        "classification_rows": classifications,
        "classification_prediction_rows": classification_predictions,
        "primary": primary,
        "evidence_level": level,
        "decision_reasons": reasons,
        "not_scored_uids": [str(row["sample_uid"]) for row in not_scored],
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


def _combined_evidence(
    continuous_level: str,
    leave_one_out: Mapping[str, object],
    leave_cluster_out: Mapping[str, object],
) -> tuple[str, list[str]]:
    classification_support = (
        float(leave_one_out["roc_auc"]) >= 0.65
        and float(leave_one_out["pr_auc_average_precision"])
        >= float(leave_one_out["prevalence"]) + 0.10
        and float(leave_one_out["mcc"]) >= 0.25
        and float(leave_one_out["balanced_accuracy"]) >= 0.60
        and min(
            float(leave_one_out["sensitivity"]),
            float(leave_one_out["specificity"]),
        )
        >= 0.50
    )
    cluster_support = (
        float(leave_cluster_out["mcc"]) > 0
        and float(leave_cluster_out["balanced_accuracy"]) > 0.50
    )
    if continuous_level == "weak_ranking_evidence" and classification_support and cluster_support:
        return "weak_ranking_evidence", ["continuous_and_nested_classification_gates_passed"]
    if continuous_level != "no_supported_use" and classification_support:
        return "compatibility_filter_only", ["classification_supported_but_full_gate_not_met"]
    return "no_supported_use", ["independent_continuous_and_classification_support_not_established"]


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
