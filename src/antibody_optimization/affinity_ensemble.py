"""Select interpretable Nb252 affinity modules from the 20-sample ensemble."""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from typing import Mapping, Sequence


CORE_SUPPORT_MINIMUM = 18
EXPECTED_CANDIDATES = 50
EXPECTED_CORE_COUNT = 8


class AffinityEnsembleError(ValueError):
    """Raised when ensemble evidence or selection provenance is invalid."""


def select_affinity_core_modules(
    ensemble_rows: Sequence[Mapping[str, object]],
    post_scan_rows: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """Apply a repeat-support gate, retain risks, and assign global Pareto fronts.

    The only selection gate requires both affinity metrics to be negative in at
    least 18 of 20 samples and to have negative medians.  Contact retention,
    fa_rep, prepared-structure sensitivity, and chemical/site diversity remain
    explicit evidence and are never collapsed into a weighted score.
    """

    if len(ensemble_rows) != EXPECTED_CANDIDATES:
        raise AffinityEnsembleError("Expected exactly 50 ensemble candidates")
    if len({str(row["candidate_id"]) for row in ensemble_rows}) != EXPECTED_CANDIDATES:
        raise AffinityEnsembleError("Ensemble candidate IDs must be unique")
    post_by_id = {str(row["candidate_id"]): row for row in post_scan_rows}
    if len(post_by_id) != 456:
        raise AffinityEnsembleError("Post-scan evidence must contain 456 candidates")

    evidence: list[dict[str, object]] = []
    for source in ensemble_rows:
        candidate_id = str(source["candidate_id"])
        if candidate_id not in post_by_id:
            raise AffinityEnsembleError(f"Missing post-scan evidence for {candidate_id}")
        post = post_by_id[candidate_id]
        if source["candidate_selection_performed"] not in (False, "False"):
            raise AffinityEnsembleError("Production input must be unfiltered")
        numeric = {
            field: float(source[field])
            for field in (
                "delta_dG_separated_median",
                "delta_dG_separated_mad",
                "delta_cross_interface_energy_median",
                "delta_cross_interface_energy_mad",
                "delta_interface_fa_rep_median",
                "minimum_vhh_contact_retention",
                "minimum_receptor_epitope_retention",
            )
        }
        if not all(math.isfinite(value) for value in numeric.values()):
            raise AffinityEnsembleError(f"Non-finite evidence for {candidate_id}")
        dg_support = int(source["negative_delta_dG_count"])
        cross_support = int(source["negative_delta_cross_interface_count"])
        if not 0 <= dg_support <= 20 or not 0 <= cross_support <= 20:
            raise AffinityEnsembleError(f"Invalid support count for {candidate_id}")
        selected = (
            dg_support >= CORE_SUPPORT_MINIMUM
            and cross_support >= CORE_SUPPORT_MINIMUM
            and numeric["delta_dG_separated_median"] < 0
            and numeric["delta_cross_interface_energy_median"] < 0
        )
        risks = _tokens(str(source["risk_flags"]))
        if str(post["prepared_contact_sensitive"]).lower() == "true":
            risks.add("prepared_contact_sensitive")
        if numeric["delta_interface_fa_rep_median"] > 0:
            risks.add("ensemble_fa_rep_increase")
        if numeric["minimum_vhh_contact_retention"] < 1:
            risks.add("ensemble_vhh_contact_change")
        if numeric["minimum_receptor_epitope_retention"] < 1:
            risks.add("ensemble_receptor_epitope_contact_change")
        evidence.append(
            {
                "candidate_id": candidate_id,
                "source_tier": source["tier"],
                "sequence_index_1based": int(source["sequence_index_1based"]),
                "wt_residue": source["wt_residue"],
                "mutant_residue": source["mutant_residue"],
                "region": post["region"],
                "mutation_reported_label": source["mutation_reported_label"],
                "mutation_numbering_label": source["mutation_numbering_label"],
                "mutation_source_auth_label": source["mutation_source_auth_label"],
                "negative_delta_dG_count": dg_support,
                "negative_delta_cross_interface_count": cross_support,
                **numeric,
                "prepared_contact_sensitive": str(post["prepared_contact_sensitive"]).lower() == "true",
                "risk_flags": ";".join(sorted(risks)),
                "core_support_gate": "pass" if selected else "not_selected",
                "core_module_selected": selected,
                "selection_reason": (
                    "both_metrics_negative_in_at_least_18_of_20_and_negative_medians"
                    if selected
                    else "repeat_support_or_median_direction_below_core_gate"
                ),
                "pareto_front": 0,
                "candidate_selection_performed": True,
                "combination_mutations_generated": False,
            }
        )
    _assign_pareto_fronts(evidence)
    cores = [row.copy() for row in evidence if row["core_module_selected"]]
    if len(cores) != EXPECTED_CORE_COUNT:
        raise AffinityEnsembleError(f"Expected 8 core modules, observed {len(cores)}")
    position_groups = _position_groups(cores)
    return {
        "evidence_rows": evidence,
        "core_rows": cores,
        "position_rows": position_groups,
        "counts": {
            "candidate_count": len(evidence),
            "core_module_count": len(cores),
            "core_position_count": len(position_groups),
            "source_tier_counts": dict(Counter(str(row["source_tier"]) for row in cores)),
        },
    }


def _assign_pareto_fronts(rows: Sequence[dict[str, object]]) -> None:
    remaining = list(rows)
    front = 1
    while remaining:
        current = [
            row
            for row in remaining
            if not any(_dominates(other, row) for other in remaining if other is not row)
        ]
        if not current:
            raise AffinityEnsembleError("Could not resolve Pareto fronts")
        for row in current:
            row["pareto_front"] = front
        selected_ids = {id(row) for row in current}
        remaining = [row for row in remaining if id(row) not in selected_ids]
        front += 1


def _dominates(first: Mapping[str, object], second: Mapping[str, object]) -> bool:
    first_values = _objectives(first)
    second_values = _objectives(second)
    return all(a <= b for a, b in zip(first_values, second_values)) and any(
        a < b for a, b in zip(first_values, second_values)
    )


def _objectives(row: Mapping[str, object]) -> tuple[float, ...]:
    return (
        -float(row["negative_delta_dG_count"]),
        -float(row["negative_delta_cross_interface_count"]),
        float(row["delta_dG_separated_median"]),
        float(row["delta_cross_interface_energy_median"]),
        float(row["delta_dG_separated_mad"]),
        float(row["delta_cross_interface_energy_mad"]),
        float(row["delta_interface_fa_rep_median"]),
        -float(row["minimum_vhh_contact_retention"]),
        -float(row["minimum_receptor_epitope_retention"]),
    )


def _position_groups(rows: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    grouped: dict[int, list[Mapping[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[int(row["sequence_index_1based"])].append(row)
    output = []
    for position, members in sorted(grouped.items()):
        output.append(
            {
                "sequence_index_1based": position,
                "wt_residue": members[0]["wt_residue"],
                "region": members[0]["region"],
                "core_module_count": len(members),
                "core_candidate_ids": ";".join(str(row["candidate_id"]) for row in members),
                "mutant_residues": ";".join(str(row["mutant_residue"]) for row in members),
                "same_position_modules_mutually_exclusive": len(members) > 1,
            }
        )
    return output


def _tokens(value: str) -> set[str]:
    return {token for token in value.split(";") if token}
