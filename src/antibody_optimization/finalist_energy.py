"""Paired Rosetta energy decomposition for the Nb252 finalist pool.

The functions in this module reuse completed, protocol-matched mutant and WT
records.  They do not run PyRosetta.  Absolute scores are never compared
between protocol families; every reported value is mutant minus its recorded
paired WT.  The separated-state value is a Rosetta ranking proxy, not a
measured folding free energy, melting temperature, or monomer stability.
"""

from __future__ import annotations

import math
import statistics
from collections import Counter
from typing import Mapping, Sequence


class FinalistEnergyError(ValueError):
    """Raised when finalist identities or paired Rosetta records disagree."""


EXPECTED_PANEL = 30
EXPECTED_RESERVES = 6
EXPECTED_REPLICATES = 3


def build_finalist_energy_review(
    panel_rows: Sequence[Mapping[str, object]],
    reserve_rows: Sequence[Mapping[str, object]],
    *,
    affinity_paired: Sequence[Mapping[str, object]],
    affinity_wt: Sequence[Mapping[str, object]],
    property_paired: Sequence[Mapping[str, object]],
    property_wt: Sequence[Mapping[str, object]],
    double_paired: Sequence[Mapping[str, object]],
    double_wt: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """Build replicate, summary, and decision-template rows for 36 finalists.

    Affinity singles deliberately use their completed three-replicate full-scan
    records for this decomposition.  Their 20-sample Flex ddG evidence remains
    the primary affinity evidence in the preliminary panel and is not replaced.
    """

    panel = _unique(panel_rows, EXPECTED_PANEL, "preliminary panel")
    reserves = _unique(reserve_rows, EXPECTED_RESERVES, "reserves")
    overlap = set(panel) & set(reserves)
    if overlap:
        raise FinalistEnergyError(f"Panel/reserve overlap: {sorted(overlap)}")
    candidates = {**panel, **reserves}

    sources = {
        "affinity_single_full_scan_3rep": _paired_source(
            affinity_paired, affinity_wt, "affinity single"
        ),
        "property_single_local_3rep": _paired_source(
            property_paired, property_wt, "property single"
        ),
        "double_position_pair_local_3rep": _paired_source(
            double_paired, double_wt, "double mutant"
        ),
    }
    replicate_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    decision_rows: list[dict[str, object]] = []

    for identifier, candidate in candidates.items():
        family = _protocol_family(candidate)
        paired_by_candidate, wt_by_id = sources[family]
        selected = paired_by_candidate.get(identifier, [])
        if len(selected) != EXPECTED_REPLICATES:
            raise FinalistEnergyError(
                f"{identifier} has {len(selected)} {family} records; expected 3"
            )
        selected = sorted(selected, key=lambda row: int(row["replicate"]))
        expected_replicates = list(range(1, EXPECTED_REPLICATES + 1))
        if [int(row["replicate"]) for row in selected] != expected_replicates:
            raise FinalistEnergyError(f"Invalid replicate identities for {identifier}")

        candidate_replicates = []
        for paired in selected:
            wt_id = str(paired["wt_control_id"])
            if wt_id not in wt_by_id:
                raise FinalistEnergyError(f"Missing paired WT {wt_id} for {identifier}")
            wt = wt_by_id[wt_id]
            replicate = int(paired["replicate"])
            seed = int(paired["seed"])
            if replicate != int(wt["replicate"]) or seed != int(wt["seed"]):
                raise FinalistEnergyError(f"WT replicate/seed mismatch for {identifier}")
            if str(paired.get("status")) != "pass" or str(wt.get("status")) != "pass":
                raise FinalistEnergyError(f"Non-pass paired record for {identifier}")

            mutant_total = _finite(paired, "mutant_total_score")
            mutant_dg = _finite(paired, "mutant_dG_separated")
            wt_total = _finite(wt, "total_score")
            wt_dg = _finite(wt, "dG_separated")
            delta_complex = mutant_total - wt_total
            mutant_separated = mutant_total - mutant_dg
            wt_separated = wt_total - wt_dg
            delta_separated = mutant_separated - wt_separated
            delta_dg = _finite(paired, "delta_dG_separated")
            identity_error = delta_dg - (delta_complex - delta_separated)
            if abs(identity_error) > 1e-8:
                raise FinalistEnergyError(
                    f"Energy identity failed for {identifier} replicate {replicate}"
                )
            row = {
                "candidate_id": identifier,
                "mutation_set": str(candidate["mutation_set"]),
                "candidate_kind": str(candidate["candidate_kind"]),
                "panel_category": str(candidate["panel_category"]),
                "current_pool_status": (
                    "preliminary_panel" if identifier in panel else "reserve"
                ),
                "energy_protocol_family": family,
                "replicate": replicate,
                "seed": seed,
                "wt_control_id": wt_id,
                "mutant_complex_total_score": mutant_total,
                "paired_wt_complex_total_score": wt_total,
                "mutant_separated_score": mutant_separated,
                "paired_wt_separated_score": wt_separated,
                "delta_complex_total_score": delta_complex,
                "delta_separated_state_score": delta_separated,
                "delta_dG_separated": delta_dg,
                "delta_cross_interface_energy": _finite(
                    paired, "delta_cross_interface_energy"
                ),
                "energy_identity_error": identity_error,
                "separated_state_interpretation": (
                    "paired_separated_system_proxy_not_measured_monomer_stability"
                ),
            }
            replicate_rows.append(row)
            candidate_replicates.append(row)

        summary = _summarize(candidate, candidate_replicates)
        summary_rows.append(summary)
        decision_rows.append(
            {
                **summary,
                "sequence": str(candidate["sequence"]),
                "expert_risk_level": str(candidate.get("expert_risk_level", "")),
                "risk_flags": str(candidate.get("risk_flags", "")),
                "pyrosetta_contact_change_status": str(
                    candidate.get("pyrosetta_contact_change_status", "")
                ),
                "property_material_favorable_count": str(
                    candidate.get("property_material_favorable_count", "")
                ),
                "property_material_adverse_count": str(
                    candidate.get("property_material_adverse_count", "")
                ),
                "antifold_delta_log_probability": str(
                    candidate.get("antifold_delta_log_probability", "")
                ),
                "tnp_flag_regression_count": str(
                    candidate.get("tnp_flag_regression_count", "")
                ),
                "review_decision": "pending",
                "review_rationale": "",
                "final_candidate_selection_performed": False,
            }
        )

    replicate_rows.sort(key=lambda row: (str(row["candidate_id"]), int(row["replicate"])))
    summary_rows.sort(key=lambda row: str(row["candidate_id"]))
    decision_rows.sort(
        key=lambda row: (
            0 if row["current_pool_status"] == "preliminary_panel" else 1,
            str(row["candidate_id"]),
        )
    )
    counts = Counter(str(row["energy_origin_class"]) for row in summary_rows)
    return {
        "replicate_rows": replicate_rows,
        "summary_rows": summary_rows,
        "decision_rows": decision_rows,
        "facts": {
            "candidate_count": len(candidates),
            "preliminary_panel_count": len(panel),
            "reserve_count": len(reserves),
            "replicate_row_count": len(replicate_rows),
            "energy_origin_class_counts": dict(sorted(counts.items())),
            "strong_separated_destabilization_caution_count": sum(
                row["energy_origin_review_flag"]
                == "strong_separated_destabilization_caution"
                for row in summary_rows
            ),
        },
    }


def _summarize(
    candidate: Mapping[str, object], rows: Sequence[Mapping[str, object]]
) -> dict[str, object]:
    complex_values = [float(row["delta_complex_total_score"]) for row in rows]
    separated_values = [float(row["delta_separated_state_score"]) for row in rows]
    dg_values = [float(row["delta_dG_separated"]) for row in rows]
    cross_values = [float(row["delta_cross_interface_energy"]) for row in rows]
    complex_negative = sum(value < 0 for value in complex_values)
    separated_negative = sum(value < 0 for value in separated_values)
    separated_positive = sum(value > 0 for value in separated_values)
    dg_negative = sum(value < 0 for value in dg_values)
    cross_negative = sum(value < 0 for value in cross_values)

    if complex_negative == 3 and separated_negative == 3:
        origin = "complex_and_separated_state_stabilization"
    elif complex_negative == 3 and separated_positive < 3:
        origin = "complex_stabilization_without_consistent_separated_destabilization"
    elif dg_negative == 3 and separated_positive == 3 and complex_negative == 0:
        origin = "apparent_binding_gain_driven_by_separated_destabilization"
    elif separated_positive == 3:
        origin = "consistent_separated_destabilization_caution"
    else:
        origin = "mixed_or_noisy_energy_origin"
    review_flag = (
        "strong_separated_destabilization_caution"
        if dg_negative == 3 and separated_positive == 3 and complex_negative == 0
        else "review_with_existing_multitool_evidence"
    )
    return {
        "candidate_id": str(candidate["candidate_id"]),
        "mutation_set": str(candidate["mutation_set"]),
        "candidate_kind": str(candidate["candidate_kind"]),
        "panel_category": str(candidate["panel_category"]),
        "current_pool_status": str(candidate["selection_status"]),
        "energy_protocol_family": str(rows[0]["energy_protocol_family"]),
        "replicate_count": len(rows),
        "delta_complex_total_score_median": statistics.median(complex_values),
        "delta_complex_total_score_mad": _mad(complex_values),
        "delta_separated_state_score_median": statistics.median(separated_values),
        "delta_separated_state_score_mad": _mad(separated_values),
        "delta_dG_separated_median": statistics.median(dg_values),
        "delta_dG_separated_mad": _mad(dg_values),
        "delta_cross_interface_energy_median": statistics.median(cross_values),
        "delta_cross_interface_energy_mad": _mad(cross_values),
        "complex_favorable_replicate_count": complex_negative,
        "separated_state_favorable_replicate_count": separated_negative,
        "separated_state_adverse_replicate_count": separated_positive,
        "binding_favorable_replicate_count": dg_negative,
        "cross_interface_favorable_replicate_count": cross_negative,
        "energy_origin_class": origin,
        "energy_origin_review_flag": review_flag,
        "absolute_rosetta_scores_compared_across_protocols": False,
        "separated_state_is_measured_monomer_stability": False,
    }


def _protocol_family(candidate: Mapping[str, object]) -> str:
    kind = str(candidate["candidate_kind"])
    category = str(candidate["panel_category"])
    if kind == "double_mutant":
        return "double_position_pair_local_3rep"
    if kind != "single_mutant":
        raise FinalistEnergyError(f"Unexpected candidate kind: {kind}")
    if category == "affinity_focused_single":
        return "affinity_single_full_scan_3rep"
    if category == "property_focused_single":
        return "property_single_local_3rep"
    raise FinalistEnergyError(f"Unexpected single category: {category}")


def _paired_source(
    paired_rows: Sequence[Mapping[str, object]],
    wt_rows: Sequence[Mapping[str, object]],
    label: str,
) -> tuple[dict[str, list[Mapping[str, object]]], dict[str, Mapping[str, object]]]:
    by_candidate: dict[str, list[Mapping[str, object]]] = {}
    for row in paired_rows:
        by_candidate.setdefault(str(row["candidate_id"]), []).append(row)
    wt: dict[str, Mapping[str, object]] = {}
    for row in wt_rows:
        identifier = str(row["wt_control_id"])
        if identifier in wt and dict(wt[identifier]) != dict(row):
            raise FinalistEnergyError(f"Conflicting {label} WT: {identifier}")
        wt[identifier] = row
    return by_candidate, wt


def _unique(
    rows: Sequence[Mapping[str, object]], expected: int, label: str
) -> dict[str, Mapping[str, object]]:
    result: dict[str, Mapping[str, object]] = {}
    for row in rows:
        identifier = str(row["candidate_id"])
        if identifier in result:
            raise FinalistEnergyError(f"Duplicate {label} candidate: {identifier}")
        result[identifier] = row
    if len(result) != expected:
        raise FinalistEnergyError(f"Expected {expected} {label} rows, found {len(result)}")
    return result


def _finite(row: Mapping[str, object], field: str) -> float:
    try:
        value = float(row[field])
    except (KeyError, TypeError, ValueError) as exc:
        raise FinalistEnergyError(f"Missing/non-numeric {field}") from exc
    if not math.isfinite(value):
        raise FinalistEnergyError(f"Non-finite {field}")
    return value


def _mad(values: Sequence[float]) -> float:
    center = statistics.median(values)
    return statistics.median(abs(value - center) for value in values)
