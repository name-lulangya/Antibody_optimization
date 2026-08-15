"""NetSolP–BL21 yield validation using the existing fixed phenotype semantics."""

from __future__ import annotations

from collections import defaultdict
from typing import Mapping, Sequence

import numpy as np

from .nanobert_yield import (
    _kendall,
    _optional_float,
    _spearman,
    build_validation_inputs,
    classify_primary_evidence,
    loocv_spearman,
    partial_spearman_provider_length,
    sequence_features,
    stratified_bootstrap_ci,
    stratified_permutation_p,
    stratified_spearman,
)


PRIMARY_FEATURE = "predicted_usability"
SECONDARY_FEATURE = "predicted_solubility"
SEQUENCE_CLUSTER_IDENTITY = 0.90


class NetSolPYieldError(ValueError):
    """Raised when NetSolP plan or score evidence violates the fixed contract."""


def build_netsolp_validation_inputs(
    expression_rows: Sequence[Mapping[str, object]],
    numbering_rows: Sequence[Mapping[str, object]],
    position_rows: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """Reuse the frozen 47-sample phenotype semantics and add sequence clusters."""

    result = build_validation_inputs(expression_rows, numbering_rows, position_rows)
    sample_rows = [dict(row) for row in result["sample_rows"]]
    clusters = sequence_identity_clusters(
        {str(row["sample_uid"]): str(row["sequence_raw"]) for row in sample_rows},
        threshold=SEQUENCE_CLUSTER_IDENTITY,
    )
    for row in sample_rows:
        row["sequence_cluster_90"] = clusters[str(row["sample_uid"])]
    return {"sample_rows": sample_rows}


def analyze_netsolp_associations(
    sample_rows: Sequence[Mapping[str, object]],
    score_rows: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """Evaluate Distilled U/S against reported yield without training a new model."""

    if len(sample_rows) != 47 or len(score_rows) != 47:
        raise NetSolPYieldError("NetSolP association analysis requires 47 samples and 47 scores")
    scores = {str(row["sample_uid"]): row for row in score_rows}
    expected = {str(row["sample_uid"]) for row in sample_rows}
    if len(scores) != 47 or set(scores) != expected:
        raise NetSolPYieldError("NetSolP score identities do not match the validation plan")

    combined: list[dict[str, object]] = []
    for sample in sample_rows:
        uid = str(sample["sample_uid"])
        score = scores[uid]
        if str(score.get("scoring_status")) != "pass":
            raise NetSolPYieldError(f"NetSolP scoring failed for {uid}")
        if str(score.get("sequence_raw")) != str(sample["sequence_raw"]):
            raise NetSolPYieldError(f"NetSolP sequence identity mismatch for {uid}")
        row = dict(sample)
        row.update(sequence_features(str(sample["sequence_raw"])))
        row[PRIMARY_FEATURE] = _probability(score.get(PRIMARY_FEATURE), PRIMARY_FEATURE)
        row[SECONDARY_FEATURE] = _probability(score.get(SECONDARY_FEATURE), SECONDARY_FEATURE)
        combined.append(row)

    features = [
        PRIMARY_FEATURE,
        SECONDARY_FEATURE,
        "sequence_length_aa",
        "molecular_weight_da",
        "theoretical_pi",
        "charge_at_ph7_4",
        "gravy",
        "aromaticity",
        "instability_index",
        "hydrophobic_fraction",
        "positive_fraction",
        "negative_fraction",
        "nxs_t_motif_count",
        "deamidation_motif_count",
        "oxidation_residue_count",
    ]
    numeric = [row for row in combined if row["observation_semantics"] == "individual_approximate"]
    llj = [row for row in combined if str(row["provider_code"]) == "LLJ"]
    metric_rows = [_metric_row(numeric, llj, feature) for feature in features]
    primary = next(row for row in metric_rows if row["feature"] == PRIMARY_FEATURE)
    low, high = stratified_bootstrap_ci(numeric, PRIMARY_FEATURE)
    primary["bootstrap_95ci_low"] = low
    primary["bootstrap_95ci_high"] = high
    primary["stratified_permutation_p"] = stratified_permutation_p(
        numeric, PRIMARY_FEATURE, float(primary["stratified_spearman_rho"])
    )
    level, reasons = classify_primary_evidence(primary)
    cluster_increment = float(primary["cluster_cv_increment_over_provider"])
    if level == "weak_ranking_evidence" and cluster_increment <= 0:
        level = "compatibility_filter_only"
        reasons = ["ordinary_gate_passed_but_sequence_cluster_cv_increment_not_positive"]
    elif level == "compatibility_filter_only" and cluster_increment <= 0:
        level = "no_supported_use"
        reasons = ["positive_partial_signals_but_sequence_cluster_cv_increment_not_positive"]
    primary["evidence_level"] = level
    return {
        "sample_rows": combined,
        "metric_rows": metric_rows,
        "primary": primary,
        "evidence_level": level,
        "decision_reasons": reasons,
    }


def normalize_netsolp_scores(
    samples: Sequence[Mapping[str, object]],
    raw_rows: Sequence[Mapping[str, object]],
    *,
    expected_count: int = 47,
) -> list[dict[str, object]]:
    """Validate official Distilled SU output and preserve plan order."""

    if expected_count < 1 or len(samples) != expected_count or len(raw_rows) != expected_count:
        raise NetSolPYieldError(
            f"Expected {expected_count} plan and NetSolP raw-output rows"
        )
    required = {"sid", "fasta", PRIMARY_FEATURE, SECONDARY_FEATURE}
    if not raw_rows or not required.issubset(raw_rows[0]):
        raise NetSolPYieldError("NetSolP raw output lacks required Distilled SU columns")
    by_id = {str(row["sid"]): row for row in raw_rows}
    expected = {str(row["sample_uid"]) for row in samples}
    if len(by_id) != expected_count or set(by_id) != expected:
        raise NetSolPYieldError("NetSolP raw output IDs do not match the 47-sample plan")
    normalized = []
    for sample in samples:
        uid = str(sample["sample_uid"])
        raw = by_id[uid]
        sequence = str(sample["sequence_raw"])
        if str(raw["fasta"]) != sequence:
            raise NetSolPYieldError(f"NetSolP raw sequence mismatch for {uid}")
        normalized.append(
            {
                "sample_uid": uid,
                "sequence_raw": sequence,
                PRIMARY_FEATURE: _probability(raw[PRIMARY_FEATURE], PRIMARY_FEATURE),
                SECONDARY_FEATURE: _probability(raw[SECONDARY_FEATURE], SECONDARY_FEATURE),
                "scoring_status": "pass",
            }
        )
    return normalized


def sequence_identity_clusters(sequences: Mapping[str, str], *, threshold: float) -> dict[str, str]:
    """Single-linkage clusters using 1-normalized-Levenshtein identity."""

    if not 0 < threshold <= 1:
        raise NetSolPYieldError("Sequence-cluster identity threshold must be in (0, 1]")
    identifiers = sorted(sequences)
    parent = {identifier: identifier for identifier in identifiers}

    def find(identifier: str) -> str:
        while parent[identifier] != identifier:
            parent[identifier] = parent[parent[identifier]]
            identifier = parent[identifier]
        return identifier

    def union(first: str, second: str) -> None:
        root_first, root_second = find(first), find(second)
        if root_first != root_second:
            parent[max(root_first, root_second)] = min(root_first, root_second)

    for index, first in enumerate(identifiers):
        for second in identifiers[index + 1 :]:
            identity = 1.0 - _levenshtein(sequences[first], sequences[second]) / max(
                len(sequences[first]), len(sequences[second])
            )
            if identity >= threshold:
                union(first, second)
    groups: dict[str, list[str]] = defaultdict(list)
    for identifier in identifiers:
        groups[find(identifier)].append(identifier)
    ordered = sorted(groups.values(), key=lambda members: members[0])
    return {identifier: f"cluster_{index:03d}" for index, members in enumerate(ordered, 1) for identifier in members}


def group_cv_spearman(
    values: np.ndarray,
    outcomes: np.ndarray,
    providers: Sequence[str],
    groups: Sequence[str],
    *,
    include_feature: bool,
) -> float:
    """Leave-one-sequence-cluster-out provider-intercept regression."""

    provider_array = np.asarray(providers)
    group_array = np.asarray(groups)
    provider_levels = sorted(set(providers))
    predictions = np.empty(len(values), dtype=float)
    for group in sorted(set(groups)):
        held = group_array == group
        train = ~held
        if int(train.sum()) < 3:
            raise NetSolPYieldError("Insufficient training rows for sequence-cluster CV")
        columns = [np.ones(int(train.sum()))]
        held_columns = [np.ones(int(held.sum()))]
        for provider in provider_levels[1:]:
            columns.append((provider_array[train] == provider).astype(float))
            held_columns.append((provider_array[held] == provider).astype(float))
        if include_feature:
            mean = float(values[train].mean())
            sd = float(values[train].std()) or 1.0
            columns.append((values[train] - mean) / sd)
            held_columns.append((values[held] - mean) / sd)
        coefficients = np.linalg.lstsq(np.column_stack(columns), np.log1p(outcomes[train]), rcond=None)[0]
        predictions[held] = np.column_stack(held_columns) @ coefficients
    return _spearman(predictions, outcomes)[0]


def _metric_row(
    numeric: Sequence[Mapping[str, object]],
    llj: Sequence[Mapping[str, object]],
    feature: str,
) -> dict[str, object]:
    values = np.asarray([float(row[feature]) for row in numeric])
    outcomes = np.asarray([float(row["numeric_yield_value"]) for row in numeric])
    providers = [str(row["provider_code"]) for row in numeric]
    groups = [str(row["sequence_cluster_90"]) for row in numeric]
    pooled_rho, pooled_p = _spearman(values, outcomes)
    tau, tau_p = _kendall(values, outcomes)
    source_stats = {
        provider: _spearman(values[np.asarray(providers) == provider], outcomes[np.asarray(providers) == provider])[0]
        for provider in sorted(set(providers))
    }
    loocv = loocv_spearman(values, outcomes, providers, include_feature=True)
    baseline = loocv_spearman(values, outcomes, providers, include_feature=False)
    cluster_cv = group_cv_spearman(values, outcomes, providers, groups, include_feature=True)
    cluster_baseline = group_cv_spearman(values, outcomes, providers, groups, include_feature=False)
    llj_values = np.asarray([float(row[feature]) for row in llj])
    llj_levels = np.asarray([int(row["llj_ordinal_level"]) for row in llj])
    return {
        "feature": feature,
        "numeric_n": len(numeric),
        "sequence_cluster_count_90": len(set(groups)),
        "pooled_spearman_rho": pooled_rho,
        "pooled_spearman_p": pooled_p,
        "pooled_kendall_tau_b": tau,
        "pooled_kendall_p": tau_p,
        "stratified_spearman_rho": stratified_spearman(values, outcomes, providers),
        "length_adjusted_partial_spearman": partial_spearman_provider_length(
            values,
            outcomes,
            np.asarray([float(row["sequence_length_aa"]) for row in numeric]),
            providers,
        ),
        "ltt_spearman_rho": source_stats.get("LTT"),
        "wcc_spearman_rho": source_stats.get("WCC"),
        "loocv_log1p_yield_spearman": loocv,
        "provider_only_loocv_spearman": baseline,
        "loocv_increment_over_provider": loocv - baseline,
        "cluster_cv_log1p_yield_spearman": cluster_cv,
        "cluster_provider_only_cv_spearman": cluster_baseline,
        "cluster_cv_increment_over_provider": cluster_cv - cluster_baseline,
        "llj_ordinal_n": len(llj),
        "llj_kendall_tau_b": _kendall(llj_values, llj_levels)[0],
    }


def _probability(value: object, field: str) -> float:
    number = _optional_float(value)
    if number is None or not 0 <= number <= 1:
        raise NetSolPYieldError(f"{field} must be a finite probability in [0, 1]")
    return number


def _levenshtein(first: str, second: str) -> int:
    previous = list(range(len(second) + 1))
    for first_index, first_residue in enumerate(first, 1):
        current = [first_index]
        for second_index, second_residue in enumerate(second, 1):
            current.append(
                min(
                    current[-1] + 1,
                    previous[second_index] + 1,
                    previous[second_index - 1] + (first_residue != second_residue),
                )
            )
        previous = current
    return previous[-1]
