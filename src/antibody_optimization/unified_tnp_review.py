"""Magnitude-aware TNP review for the released Nb252 single-mutant pools.

The module combines the 49 property-discovery Pareto-front-1 records with the
46 non-cysteine candidates from the 50-candidate Flex ddG production review.
It preserves raw predictions and prior evidence, adds conservative operational
magnitude labels, and compares TNP developability metrics with the WT control.
It does not predict yield or select final candidates.
"""

from __future__ import annotations

from collections import Counter
from typing import Mapping, Sequence


WT_SCORE_ID = "LTT__Nb252__WT"
PROPERTY_COUNT = 49
AFFINITY_COUNT = 46
CANDIDATE_COUNT = PROPERTY_COUNT + AFFINITY_COUNT
SCORE_COUNT = CANDIDATE_COUNT + 1
BLOCKED_PRODUCTION_CYS_IDS = {
    "Nb252_aff_seq037_Y37C",
    "Nb252_aff_seq045_R45C",
    "Nb252_aff_seq098_D98C",
    "Nb252_aff_seq105_E105C",
}
MAGNITUDE_THRESHOLDS = {
    "netsolp_delta_usability_vs_wt": 0.01,
    "netsolp_delta_solubility_vs_wt": 0.02,
    "nanomelt_delta_predicted_apparent_tm_c_vs_wt": 1.0,
}
TNP_METRICS = (
    "tnp_total_cdr_length",
    "tnp_cdr3_length",
    "tnp_cdr3_compactness",
    "tnp_psh",
    "tnp_ppc",
    "tnp_pnc",
)
TNP_FLAGS = (
    "tnp_flag_total_cdr_length",
    "tnp_flag_cdr3_length",
    "tnp_flag_cdr3_compactness",
    "tnp_flag_psh",
    "tnp_flag_ppc",
    "tnp_flag_pnc",
)


class UnifiedTNPReviewError(ValueError):
    """Raised when the unified TNP plan or evidence is inconsistent."""


def magnitude_label(value: object, threshold: float) -> str:
    """Classify a WT-relative value with inclusive neutral boundaries."""

    number = float(value)
    if number > threshold:
        return "favorable"
    if number < -threshold:
        return "adverse"
    return "neutral"


def build_unified_tnp_samples(
    property_rows: Sequence[Mapping[str, object]],
    property_samples: Sequence[Mapping[str, object]],
    production_rows: Sequence[Mapping[str, object]],
    unified_candidate_rows: Sequence[Mapping[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Return WT plus 95 eligible candidates and four blocked audit records."""

    evidence = _unique(property_rows, "candidate_id", 1962, "property evidence")
    plan_samples = _unique(property_samples, "score_id", 1963, "property samples")
    unified = _unique(unified_candidate_rows, "candidate_id", 2318, "unified candidates")
    production = _unique(production_rows, "candidate_id", 50, "Flex ddG production")

    property_ids = {
        identifier
        for identifier, row in evidence.items()
        if str(row["design_track"]) == "stability_developability_discovery"
        and int(row["property_pareto_layer"]) == 1
    }
    if len(property_ids) != PROPERTY_COUNT:
        raise UnifiedTNPReviewError("Expected 49 property-discovery Pareto-front-1 candidates")

    blocked = set(production) - set(evidence)
    if blocked != BLOCKED_PRODUCTION_CYS_IDS:
        raise UnifiedTNPReviewError("Flex ddG candidates missing from property evidence are not the four fixed Cys records")
    for identifier in blocked:
        row = unified[identifier]
        if str(row["design_status"]) != "blocked_new_unpaired_cys":
            raise UnifiedTNPReviewError(f"Blocked Flex ddG identity is not an unpaired-Cys record: {identifier}")
    affinity_ids = set(production) - blocked
    if len(affinity_ids) != AFFINITY_COUNT or property_ids & affinity_ids:
        raise UnifiedTNPReviewError("Unified TNP source pools have unexpected size or overlap")

    wt = plan_samples.get(WT_SCORE_ID)
    if wt is None:
        raise UnifiedTNPReviewError("WT control is missing from the property scoring plan")
    parent = str(wt["sequence_raw"])
    samples: list[dict[str, object]] = [{
        "sample_uid": WT_SCORE_ID,
        "candidate_id": "WT",
        "sequence_raw": parent,
        "candidate_source": "wild_type_control",
        "design_track": "wild_type_control",
        "sequence_index_1based": "",
        "wt_residue": "",
        "mutant_residue": "",
        "region": "",
        "property_pareto_layer": "",
        "netsolp_delta_usability_vs_wt": 0.0,
        "netsolp_usability_magnitude": "neutral",
        "netsolp_delta_solubility_vs_wt": 0.0,
        "netsolp_solubility_magnitude": "neutral",
        "nanomelt_delta_predicted_apparent_tm_c_vs_wt": 0.0,
        "nanomelt_tm_magnitude": "neutral",
        "experimental_complex_context_delta_log_probability": 0.0,
        "property_magnitude_class": "wild_type_control",
        "material_favorable_count": 0,
        "material_adverse_count": 0,
        "affinity_core_module": False,
        "chemical_risk_count": 0,
        "new_liability_flags": "",
        "tnp_scoring_requested": True,
    }]

    for source_name, identifiers in (
        ("property_pareto_front_1", property_ids),
        ("affinity_flex_ddg_20_sample_pool", affinity_ids),
    ):
        for identifier in sorted(identifiers, key=lambda item: _candidate_sort_key(evidence[item])):
            source = evidence[identifier]
            sequence = str(source["sequence"])
            index = int(source["sequence_index_1based"])
            wt_residue = str(source["wt_residue"])
            mutant_residue = str(source["mutant_residue"])
            if len(sequence) != 128 or sequence[-4:] != "SSGS" or sequence[21] != "C" or sequence[94] != "C":
                raise UnifiedTNPReviewError(f"Invalid authoritative candidate sequence: {identifier}")
            reconstructed = sequence[: index - 1] + wt_residue + sequence[index:]
            if reconstructed != parent or sequence[index - 1] != mutant_residue:
                raise UnifiedTNPReviewError(f"Candidate is not a single mutation of WT: {identifier}")
            labels = {
                name: magnitude_label(source[name], threshold)
                for name, threshold in MAGNITUDE_THRESHOLDS.items()
            }
            favorable = sum(value == "favorable" for value in labels.values())
            adverse = sum(value == "adverse" for value in labels.values())
            magnitude_class = (
                "tradeoff_material_adverse" if adverse
                else "multi_signal_favorable" if favorable >= 2
                else "single_signal_favorable" if favorable == 1
                else "no_material_change"
            )
            samples.append({
                "sample_uid": identifier,
                "candidate_id": identifier,
                "sequence_raw": sequence,
                "candidate_source": source_name,
                "design_track": str(source["design_track"]),
                "sequence_index_1based": index,
                "wt_residue": wt_residue,
                "mutant_residue": mutant_residue,
                "region": str(source["region"]),
                "property_pareto_layer": int(source["property_pareto_layer"]),
                "netsolp_delta_usability_vs_wt": float(source["netsolp_delta_usability_vs_wt"]),
                "netsolp_usability_magnitude": labels["netsolp_delta_usability_vs_wt"],
                "netsolp_delta_solubility_vs_wt": float(source["netsolp_delta_solubility_vs_wt"]),
                "netsolp_solubility_magnitude": labels["netsolp_delta_solubility_vs_wt"],
                "nanomelt_delta_predicted_apparent_tm_c_vs_wt": float(source["nanomelt_delta_predicted_apparent_tm_c_vs_wt"]),
                "nanomelt_tm_magnitude": labels["nanomelt_delta_predicted_apparent_tm_c_vs_wt"],
                "experimental_complex_context_delta_log_probability": float(source["experimental_complex_context_delta_log_probability"]),
                "property_magnitude_class": magnitude_class,
                "material_favorable_count": favorable,
                "material_adverse_count": adverse,
                "affinity_core_module": _truth(source.get("affinity_core_module")),
                "chemical_risk_count": int(source["chemical_risk_count"]),
                "new_liability_flags": _clean(source.get("new_liability_flags")),
                "tnp_scoring_requested": True,
            })

    if len(samples) != SCORE_COUNT or len({str(row["sample_uid"]) for row in samples}) != SCORE_COUNT:
        raise UnifiedTNPReviewError("Unified TNP plan does not contain 96 unique records")

    audit = [{
        "candidate_id": identifier,
        "sequence_index_1based": int(unified[identifier]["sequence_index_1based"]),
        "wt_residue": str(unified[identifier]["wt_residue"]),
        "mutant_residue": str(unified[identifier]["mutant_residue"]),
        "design_status": str(unified[identifier]["design_status"]),
        "design_status_reason": str(unified[identifier]["design_status_reason"]),
        "tnp_scoring_requested": False,
        "audit_reason": "blocked_new_unpaired_cys_not_a_design_candidate",
    } for identifier in sorted(blocked)]
    return samples, audit


def analyze_unified_tnp_scores(
    sample_rows: Sequence[Mapping[str, object]],
    score_rows: Sequence[Mapping[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, object]]:
    """Join 96 TNP scores, calculate WT deltas, and retain risk semantics."""

    samples = _unique(sample_rows, "sample_uid", SCORE_COUNT, "TNP samples")
    scores = _unique(score_rows, "sample_uid", SCORE_COUNT, "TNP scores")
    if set(samples) != set(scores) or WT_SCORE_ID not in samples:
        raise UnifiedTNPReviewError("TNP score identities do not match the plan")
    for identifier, score in scores.items():
        if str(score["scoring_status"]) != "pass":
            raise UnifiedTNPReviewError(f"TNP scoring did not pass for {identifier}")
        if str(score["sequence_raw"]) != str(samples[identifier]["sequence_raw"]):
            raise UnifiedTNPReviewError(f"TNP sequence mismatch for {identifier}")
        if str(score["trimmed_n_terminal"]):
            raise UnifiedTNPReviewError(f"Unexpected TNP N-terminal trim for {identifier}")
        if str(score["trimmed_c_terminal"]) != "GS" or int(score["modelled_length_aa"]) != 126:
            raise UnifiedTNPReviewError(f"Unexpected TNP VHH modeling domain for {identifier}")

    wt_score = scores[WT_SCORE_ID]
    flag_severity = {"green": 0, "amber": 1, "red": 2}
    output: list[dict[str, object]] = []
    for identifier, sample in samples.items():
        if identifier == WT_SCORE_ID:
            continue
        score = scores[identifier]
        row = dict(sample)
        for metric in TNP_METRICS:
            row[metric] = float(score[metric])
            row[f"{metric}_delta_vs_wt"] = float(score[metric]) - float(wt_score[metric])
        regressions = 0
        improvements = 0
        new_red = 0
        for flag in TNP_FLAGS:
            current = str(score[flag]).lower()
            wt_value = str(wt_score[flag]).lower()
            if current not in flag_severity or wt_value not in flag_severity:
                raise UnifiedTNPReviewError(f"Invalid TNP flag for {identifier}: {flag}")
            row[flag] = current
            row[f"{flag}_vs_wt"] = _flag_transition(wt_value, current, flag_severity)
            regressions += flag_severity[current] > flag_severity[wt_value]
            improvements += flag_severity[current] < flag_severity[wt_value]
            new_red += current == "red" and wt_value != "red"
        row["tnp_flag_regression_count"] = regressions
        row["tnp_flag_improvement_count"] = improvements
        row["tnp_new_red_flag_count"] = new_red
        row["tnp_developability_review"] = (
            "new_red_flag" if new_red else
            "flag_regression" if regressions else
            "flag_improvement" if improvements else
            "no_flag_change"
        )
        row["tnp_yield_prediction_performed"] = False
        row["candidate_selection_performed"] = False
        output.append(row)

    output.sort(key=lambda row: (str(row["candidate_source"]), int(row["sequence_index_1based"]), str(row["mutant_residue"])))
    summary_counter = Counter(
        (str(row["candidate_source"]), str(row["property_magnitude_class"]), str(row["tnp_developability_review"]))
        for row in output
    )
    summary = [{
        "candidate_source": source,
        "property_magnitude_class": magnitude,
        "tnp_developability_review": review,
        "candidate_count": count,
    } for (source, magnitude, review), count in sorted(summary_counter.items())]
    gate = {
        "schema_version": 1,
        "gate_name": "nb252_unified_tnp_candidate_review",
        "status": "pass",
        "planned_count": SCORE_COUNT,
        "candidate_count": CANDIDATE_COUNT,
        "wt_control_count": 1,
        "pass_count": len(scores),
        "source_counts": dict(sorted(Counter(str(row["candidate_source"]) for row in output).items())),
        "flag_regression_count": sum(int(row["tnp_flag_regression_count"]) > 0 for row in output),
        "new_red_flag_candidate_count": sum(int(row["tnp_new_red_flag_count"]) > 0 for row in output),
        "yield_prediction_performed": False,
        "candidate_selection_performed": False,
        "release": "ready_for_magnitude_aware_multitool_shortlisting",
        "interpretation": "TNP developability risk relative to WT; not measured expression, yield, affinity, or a standalone selection score.",
    }
    return output, summary, gate


def _candidate_sort_key(row: Mapping[str, object]) -> tuple[int, str]:
    return int(row["sequence_index_1based"]), str(row["mutant_residue"])


def _unique(
    rows: Sequence[Mapping[str, object]], key: str, expected: int, label: str
) -> dict[str, Mapping[str, object]]:
    result = {str(row[key]): row for row in rows}
    if len(rows) != expected or len(result) != expected:
        raise UnifiedTNPReviewError(f"{label} must contain {expected} unique rows")
    return result


def _truth(value: object) -> bool:
    return value is True or str(value).strip().lower() == "true"


def _clean(value: object) -> str:
    text = str(value).strip()
    return "" if text.lower() == "nan" else text


def _flag_transition(wt_value: str, current: str, severity: Mapping[str, int]) -> str:
    if severity[current] > severity[wt_value]:
        return "worse"
    if severity[current] < severity[wt_value]:
        return "better"
    return "same"
