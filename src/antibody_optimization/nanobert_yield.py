"""Contracts and statistics for nanoBERT–reported-yield validation."""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from typing import Mapping, Sequence

import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from scipy.stats import kendalltau, rankdata, spearmanr


AA20 = frozenset("ACDEFGHIKLMNPQRSTVWY")
HYDROPHOBIC = frozenset("AVILMFWY")
POSITIVE = frozenset("KRH")
NEGATIVE = frozenset("DE")
PRIMARY_FEATURE = "nanobert_mean_pll_raw"
BOOTSTRAP_SEED = 252031
BOOTSTRAP_REPLICATES = 5000
PERMUTATION_REPLICATES = 5000


class NanoBertYieldError(ValueError):
    """Raised when source or model-score evidence violates the fixed contract."""


def build_validation_inputs(
    expression_rows: Sequence[Mapping[str, object]],
    numbering_rows: Sequence[Mapping[str, object]],
    position_rows: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """Join frozen expression semantics to provisional IMGT regions."""

    if len(expression_rows) != 47 or len(numbering_rows) != 47:
        raise NanoBertYieldError("Expected exactly 47 expression and numbering rows")
    numbering = {str(row["sample_uid"]): row for row in numbering_rows}
    regions: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in position_rows:
        if str(row["is_gap"]).lower() == "true":
            continue
        regions[str(row["sample_uid"])].append(
            {
                "sample_uid": row["sample_uid"],
                "sequence_index_1based": int(row["sequence_index_1based"]),
                "region": row["region"],
                "numbering_position_label": row["numbering_position_label"],
                "residue_aa": row["residue_aa"],
            }
        )
    samples = []
    for row in expression_rows:
        uid = str(row["sample_uid"])
        sequence = str(row["sequence_raw"])
        if not sequence or set(sequence) - AA20:
            raise NanoBertYieldError(f"Invalid amino-acid sequence for {uid}")
        nr = numbering.get(uid)
        if nr is None or str(nr["sequence_raw"]) != sequence:
            raise NanoBertYieldError(f"Sequence/numbering identity mismatch for {uid}")
        semantics = str(row["observation_semantics"])
        numeric = semantics == "individual_approximate"
        if numeric:
            phenotype_value = float(row["point_estimate_mg"])
            ordinal_level = ""
        elif semantics == "group_approximate":
            anchor = float(row["group_anchor_mg"])
            phenotype_value = ""
            ordinal_level = 1 if anchor == 2 else 2 if anchor == 10 else None
        elif semantics == "group_lower_bound" and float(row["lower_bound_mg"]) == 20:
            phenotype_value = ""
            ordinal_level = 3
        else:
            raise NanoBertYieldError(f"Unsupported phenotype semantics for {uid}")
        if ordinal_level is None:
            raise NanoBertYieldError(f"Unsupported LLJ ordinal anchor for {uid}")
        samples.append(
            {
                "sample_uid": uid,
                "provider_code": row["provider_code"],
                "source_sample_id": row["source_sample_id"],
                "sequence_raw": sequence,
                "sequence_length_aa": len(sequence),
                "numbering_status": nr["numbering_status"],
                "numbered_start_0based": nr["query_start_0based_inclusive"],
                "numbered_end_0based": nr["query_end_0based_inclusive"],
                "observation_semantics": semantics,
                "reported_text": row["reported_text"],
                "numeric_yield_value": phenotype_value,
                "llj_ordinal_level": ordinal_level,
            }
        )
    counts = Counter(str(row["observation_semantics"]) for row in samples)
    if counts != Counter({"individual_approximate": 31, "group_lower_bound": 9, "group_approximate": 7}):
        raise NanoBertYieldError(f"Unexpected phenotype counts: {dict(counts)}")
    if Counter(str(row["numbering_status"]) for row in samples) != Counter({"pass": 46, "failed": 1}):
        raise NanoBertYieldError("Unexpected numbering status counts")
    flat_regions = [entry for uid in sorted(regions) for entry in regions[uid]]
    return {"sample_rows": samples, "region_rows": flat_regions}


def sequence_features(sequence: str) -> dict[str, float]:
    """Return a small predeclared physicochemical baseline."""

    analysis = ProteinAnalysis(sequence)
    length = len(sequence)
    return {
        "sequence_length_aa": float(length),
        "molecular_weight_da": float(analysis.molecular_weight()),
        "theoretical_pi": float(analysis.isoelectric_point()),
        "charge_at_ph7_4": float(analysis.charge_at_pH(7.4)),
        "gravy": float(analysis.gravy()),
        "aromaticity": float(analysis.aromaticity()),
        "instability_index": float(analysis.instability_index()),
        "hydrophobic_fraction": sum(aa in HYDROPHOBIC for aa in sequence) / length,
        "positive_fraction": sum(aa in POSITIVE for aa in sequence) / length,
        "negative_fraction": sum(aa in NEGATIVE for aa in sequence) / length,
        "nxs_t_motif_count": float(sum(sequence[i] == "N" and sequence[i + 1] != "P" and sequence[i + 2] in "ST" for i in range(length - 2))),
        "deamidation_motif_count": float(sum(sequence[i] == "N" and sequence[i + 1] in "GST" for i in range(length - 1))),
        "oxidation_residue_count": float(sum(aa in "MW" for aa in sequence)),
    }


def analyze_associations(
    sample_rows: Sequence[Mapping[str, object]],
    score_rows: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """Evaluate nanoBERT and baseline features without fitting a high-capacity model."""

    if len(sample_rows) != 47 or len(score_rows) != 47:
        raise NanoBertYieldError("Association analysis requires 47 samples and 47 scores")
    scores = {str(row["sample_uid"]): row for row in score_rows}
    if len(scores) != 47 or set(scores) != {str(row["sample_uid"]) for row in sample_rows}:
        raise NanoBertYieldError("nanoBERT score identities do not match the validation plan")
    combined = []
    for sample in sample_rows:
        score = scores[str(sample["sample_uid"])]
        if str(score["scoring_status"]) != "pass":
            raise NanoBertYieldError(f'nanoBERT scoring failed for {sample["sample_uid"]}')
        row = dict(sample)
        row.update(sequence_features(str(sample["sequence_raw"])))
        for field in ("nanobert_mean_pll_raw", "nanobert_sum_pll_raw", "nanobert_mean_pll_numbered", "nanobert_mean_pll_fr", "nanobert_mean_pll_cdr"):
            row[field] = _optional_float(score.get(field, ""))
        combined.append(row)
    feature_names = [
        PRIMARY_FEATURE, "nanobert_mean_pll_numbered", "nanobert_mean_pll_fr", "nanobert_mean_pll_cdr",
        "sequence_length_aa", "molecular_weight_da", "theoretical_pi", "charge_at_ph7_4", "gravy",
        "aromaticity", "instability_index", "hydrophobic_fraction", "positive_fraction", "negative_fraction",
        "nxs_t_motif_count", "deamidation_motif_count", "oxidation_residue_count",
    ]
    numeric = [row for row in combined if row["observation_semantics"] == "individual_approximate"]
    llj = [row for row in combined if str(row["provider_code"]) == "LLJ"]
    metrics = []
    for feature in feature_names:
        usable = [row for row in numeric if row[feature] is not None]
        values = np.asarray([float(row[feature]) for row in usable])
        yields = np.asarray([float(row["numeric_yield_value"]) for row in usable])
        providers = [str(row["provider_code"]) for row in usable]
        pooled_rho, pooled_p = _spearman(values, yields)
        tau, tau_p = _kendall(values, yields)
        stratified = stratified_spearman(values, yields, providers)
        length_adjusted = partial_spearman_provider_length(
            values,
            yields,
            np.asarray([float(row["sequence_length_aa"]) for row in usable]),
            providers,
        )
        source_stats = {provider: _spearman(values[np.asarray(providers) == provider], yields[np.asarray(providers) == provider])[0] for provider in sorted(set(providers))}
        loocv = loocv_spearman(values, yields, providers, include_feature=True)
        baseline = loocv_spearman(values, yields, providers, include_feature=False)
        metric = {
            "feature": feature, "numeric_n": len(usable), "pooled_spearman_rho": pooled_rho,
            "pooled_spearman_p": pooled_p, "pooled_kendall_tau_b": tau, "pooled_kendall_p": tau_p,
            "stratified_spearman_rho": stratified, "length_adjusted_partial_spearman": length_adjusted,
            "ltt_spearman_rho": source_stats.get("LTT"),
            "wcc_spearman_rho": source_stats.get("WCC"), "loocv_log1p_yield_spearman": loocv,
            "provider_only_loocv_spearman": baseline, "loocv_increment_over_provider": loocv - baseline,
        }
        llj_usable = [row for row in llj if row[feature] is not None]
        metric["llj_ordinal_n"] = len(llj_usable)
        metric["llj_kendall_tau_b"] = _kendall(
            np.asarray([float(row[feature]) for row in llj_usable]),
            np.asarray([int(row["llj_ordinal_level"]) for row in llj_usable]),
        )[0]
        metrics.append(metric)
    primary = next(row for row in metrics if row["feature"] == PRIMARY_FEATURE)
    ci_low, ci_high = stratified_bootstrap_ci(numeric, PRIMARY_FEATURE)
    permutation_p = stratified_permutation_p(numeric, PRIMARY_FEATURE, float(primary["stratified_spearman_rho"]))
    primary.update({"bootstrap_95ci_low": ci_low, "bootstrap_95ci_high": ci_high, "stratified_permutation_p": permutation_p})
    evidence_level, reasons = classify_primary_evidence(primary)
    return {"sample_rows": combined, "metric_rows": metrics, "primary": primary, "evidence_level": evidence_level, "decision_reasons": reasons}


def stratified_spearman(values: np.ndarray, outcomes: np.ndarray, providers: Sequence[str]) -> float:
    """Correlate within-provider ranks after removing provider mean ranks."""

    xres, yres = [], []
    provider_array = np.asarray(providers)
    for provider in sorted(set(providers)):
        mask = provider_array == provider
        xr = rankdata(values[mask], method="average")
        yr = rankdata(outcomes[mask], method="average")
        xres.extend(xr - xr.mean())
        yres.extend(yr - yr.mean())
    correlation = np.corrcoef(xres, yres)[0, 1]
    return float(correlation) if math.isfinite(correlation) else 0.0


def loocv_spearman(values: np.ndarray, outcomes: np.ndarray, providers: Sequence[str], *, include_feature: bool) -> float:
    """LOOCV a provider-intercept plus optional single-feature log-yield model."""

    provider_levels = sorted(set(providers))
    predictions = []
    for held_out in range(len(values)):
        train = np.arange(len(values)) != held_out
        columns = [np.ones(int(train.sum()))]
        held_columns = [1.0]
        for provider in provider_levels[1:]:
            columns.append(np.asarray(providers)[train] == provider)
            held_columns.append(float(providers[held_out] == provider))
        if include_feature:
            mean, sd = float(values[train].mean()), float(values[train].std()) or 1.0
            columns.append((values[train] - mean) / sd)
            held_columns.append((float(values[held_out]) - mean) / sd)
        design = np.column_stack(columns).astype(float)
        coefficients = np.linalg.lstsq(design, np.log1p(outcomes[train]), rcond=None)[0]
        predictions.append(float(np.dot(held_columns, coefficients)))
    return _spearman(np.asarray(predictions), outcomes)[0]


def partial_spearman_provider_length(
    values: np.ndarray,
    outcomes: np.ndarray,
    lengths: np.ndarray,
    providers: Sequence[str],
) -> float:
    """Partial rank correlation after provider and sequence-length adjustment."""

    provider_levels = sorted(set(providers))
    design_columns = [np.ones(len(values)), rankdata(lengths, method="average")]
    for provider in provider_levels[1:]:
        design_columns.append((np.asarray(providers) == provider).astype(float))
    design = np.column_stack(design_columns)
    ranked_values = rankdata(values, method="average")
    ranked_outcomes = rankdata(outcomes, method="average")
    value_residual = ranked_values - design @ np.linalg.lstsq(design, ranked_values, rcond=None)[0]
    outcome_residual = ranked_outcomes - design @ np.linalg.lstsq(design, ranked_outcomes, rcond=None)[0]
    correlation = np.corrcoef(value_residual, outcome_residual)[0, 1]
    return float(correlation) if math.isfinite(correlation) else 0.0


def stratified_bootstrap_ci(rows: Sequence[Mapping[str, object]], feature: str) -> tuple[float, float]:
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    grouped = {provider: [row for row in rows if row["provider_code"] == provider] for provider in ("LTT", "WCC")}
    estimates = []
    for _ in range(BOOTSTRAP_REPLICATES):
        sample = []
        for provider, members in grouped.items():
            sample.extend(members[index] for index in rng.integers(0, len(members), len(members)))
        estimates.append(stratified_spearman(np.asarray([float(row[feature]) for row in sample]), np.asarray([float(row["numeric_yield_value"]) for row in sample]), [str(row["provider_code"]) for row in sample]))
    return tuple(float(value) for value in np.quantile(estimates, [0.025, 0.975]))


def stratified_permutation_p(rows: Sequence[Mapping[str, object]], feature: str, observed: float) -> float:
    rng = np.random.default_rng(BOOTSTRAP_SEED + 1)
    values = np.asarray([float(row[feature]) for row in rows])
    outcomes = np.asarray([float(row["numeric_yield_value"]) for row in rows])
    providers = [str(row["provider_code"]) for row in rows]
    masks = [np.asarray(providers) == provider for provider in sorted(set(providers))]
    exceed = 0
    for _ in range(PERMUTATION_REPLICATES):
        permuted = outcomes.copy()
        for mask in masks:
            permuted[mask] = rng.permutation(permuted[mask])
        exceed += abs(stratified_spearman(values, permuted, providers)) >= abs(observed)
    return (exceed + 1) / (PERMUTATION_REPLICATES + 1)


def classify_primary_evidence(metric: Mapping[str, object]) -> tuple[str, list[str]]:
    """Apply the predeclared conservative evidence gate."""

    strong = (
        float(metric["stratified_spearman_rho"]) >= 0.30
        and float(metric["bootstrap_95ci_low"]) > 0
        and float(metric["stratified_permutation_p"]) <= 0.05
        and float(metric["ltt_spearman_rho"]) > 0
        and float(metric["wcc_spearman_rho"]) > 0
        and float(metric["length_adjusted_partial_spearman"]) > 0
        and float(metric["loocv_increment_over_provider"]) >= 0.10
    )
    positive_signals = sum(
        float(metric[field]) > 0
        for field in (
            "stratified_spearman_rho", "length_adjusted_partial_spearman", "ltt_spearman_rho",
            "wcc_spearman_rho", "loocv_increment_over_provider",
        )
    )
    if strong:
        return "weak_ranking_evidence", ["all_predeclared_direction_uncertainty_and_loocv_criteria_passed"]
    if float(metric["stratified_spearman_rho"]) > 0 and positive_signals >= 4:
        return "compatibility_filter_only", ["positive_direction_but_full_weak_ranking_gate_not_met"]
    return "no_supported_use", ["primary_direction_or_cross_source_stability_not_supported"]


def _optional_float(value: object) -> float | None:
    if value in (None, ""):
        return None
    number = float(value)
    if not math.isfinite(number):
        raise NanoBertYieldError("Non-finite nanoBERT score")
    return number


def _spearman(first: np.ndarray, second: np.ndarray) -> tuple[float, float]:
    result = spearmanr(first, second)
    return _finite_or_zero(result.statistic), _finite_or_one(result.pvalue)


def _kendall(first: np.ndarray, second: np.ndarray) -> tuple[float, float]:
    result = kendalltau(first, second, variant="b")
    return _finite_or_zero(result.statistic), _finite_or_one(result.pvalue)


def _finite_or_zero(value: object) -> float:
    number = float(value)
    return number if math.isfinite(number) else 0.0


def _finite_or_one(value: object) -> float:
    number = float(value)
    return number if math.isfinite(number) else 1.0
