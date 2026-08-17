"""Build the traceable Nb252 30-sequence preliminary experimental panel.

The module keeps the 14 already released single mutants as individually
interpretable hypotheses and selects 16 balanced double mutants without a
weighted composite score.  Double-mutant Pareto comparison is restricted to
the shared V2.1 protocol; single- and double-mutant Rosetta magnitudes are
never compared directly.  The output is a computational preselection, not a
claim of measured affinity, expression, stability, or yield.
"""

from __future__ import annotations

import math
import re
from collections import Counter
from typing import Mapping, Sequence

from .unified_tnp_review import MAGNITUDE_THRESHOLDS, magnitude_label


class PreliminaryPanelError(ValueError):
    """Raised when released evidence violates the preliminary-panel contract."""


EXPECTED_SINGLES = 14
EXPECTED_DOUBLES = 86
PRIMARY_DOUBLE_CLASSES = {
    "balanced_supported": "balanced_combination",
    "affinity_supported_property_nonadverse": "affinity_supported_double",
    "property_supported_affinity_nonadverse": "property_supported_double",
}
MAIN_QUOTAS = {
    "affinity_focused_single": 8,
    "property_focused_single": 6,
    "balanced_combination": 16,
}
RESERVE_QUOTAS = {
    "balanced_combination": 2,
    "affinity_supported_double": 2,
    "property_supported_double": 2,
}
MUTATION_RE = re.compile(r"^([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWY])$")


def build_preliminary_panel(
    single_rows: Sequence[Mapping[str, object]],
    double_rows: Sequence[Mapping[str, object]],
    affinity_rows: Sequence[Mapping[str, object]],
    property_rows: Sequence[Mapping[str, object]],
    single_gate: Mapping[str, object],
    double_gate: Mapping[str, object],
    stage2_contract: Mapping[str, object],
) -> dict[str, object]:
    """Return the 100-row audit, 30-row panel, six reserves, and facts.

    All 14 released singles are retained because they provide the component
    controls needed to interpret selected combinations.  Sixteen of the 28
    ``balanced_supported`` doubles are selected by four-objective Pareto
    peeling followed by a deterministic mutation-diversity tie-break.  The
    property objectives use the already released non-micro favorable count,
    so tiny U/S/Tm differences do not drive selection.
    """

    _validate_gates(single_gate, double_gate)
    parent, immutable = _parent_contract(stage2_contract)
    singles = _unique(single_rows, "candidate_id", EXPECTED_SINGLES, "singles")
    doubles = _unique(double_rows, "candidate_id", EXPECTED_DOUBLES, "doubles")
    affinity = _unique_any(affinity_rows, "candidate_id", "affinity evidence")
    properties = _unique_any(property_rows, "candidate_id", "property evidence")

    audit: list[dict[str, object]] = []
    for identifier, source in singles.items():
        if str(source.get("shortlist_decision")) != "retain_active":
            raise PreliminaryPanelError(f"Single is not active: {identifier}")
        mutation = str(source["mutation"])
        _validate_candidate_sequence(parent, str(source["sequence"]), mutation, immutable)
        track = str(source["design_track"])
        if track == "affinity":
            if identifier not in affinity:
                raise PreliminaryPanelError(f"Missing affinity evidence for {identifier}")
            evidence = affinity[identifier]
            panel_category = "affinity_focused_single"
            protocol = "flex_ddg_20_sample_single_mutant"
            dg = _float(evidence, "delta_dG_separated_median")
            cross = _float(evidence, "delta_cross_interface_energy_median")
            dg_count = _int(evidence, "negative_delta_dG_count")
            cross_count = _int(evidence, "negative_delta_cross_interface_count")
            replicate_count = 20
            affinity_class = "affinity_core" if str(evidence["core_support_gate"]) == "pass" else "affinity_alternative"
        elif track == "property":
            if identifier not in properties:
                raise PreliminaryPanelError(f"Missing property affinity evidence for {identifier}")
            evidence = properties[identifier]
            panel_category = "property_focused_single"
            protocol = "property_local_repack_3_replicate_single_mutant"
            dg = _float(evidence, "delta_dG_separated_median")
            cross = _float(evidence, "delta_cross_interface_energy_median")
            dg_count = _int(evidence, "delta_dg_negative_replicate_count")
            cross_count = _int(evidence, "delta_cross_negative_replicate_count")
            replicate_count = 3
            affinity_class = str(evidence["affinity_direction_class"])
        else:
            raise PreliminaryPanelError(f"Unexpected single design track: {track}")

        du = _float(source, "netsolp_delta_usability_vs_wt")
        ds = _float(source, "netsolp_delta_solubility_vs_wt")
        dt = _float(source, "nanomelt_delta_predicted_tm_c_vs_wt")
        labels = _property_labels(du, ds, dt)
        flags = _tokens(source.get("hard_risk_flags")) | _tokens(source.get("structural_review_flags"))
        audit.append(
            {
                "candidate_id": identifier,
                "candidate_kind": "single_mutant",
                "mutation_count": 1,
                "mutation_set": mutation,
                "panel_category": panel_category,
                "source_evidence_class": affinity_class,
                "sequence": str(source["sequence"]),
                "primary_pool_eligible": True,
                "primary_pool_reason": "released_active_single_component_hypothesis",
                "hard_constraint_status": "pass",
                "hard_constraint_blockers": "",
                "expert_risk_level": str(source.get("expert_risk_level", "")),
                "risk_flags": ";".join(sorted(flags)),
                "position_pair": str(source["sequence_index_1based"]),
                "component_mutations": mutation,
                "pyrosetta_protocol_family": protocol,
                "pyrosetta_replicate_count": replicate_count,
                "pyrosetta_negative_dg_count": dg_count,
                "pyrosetta_negative_cross_count": cross_count,
                "pyrosetta_delta_dg_median": dg,
                "pyrosetta_delta_cross_median": cross,
                "pyrosetta_contact_change_status": "reviewed_in_single_mutant_source",
                "paired_contact_changed_replicate_count": "",
                "experimental_reference_sensitivity_status": "",
                "antifold_delta_log_probability": _float(source, "antifold_complex_delta_log_probability"),
                "antifold_evidence_scope": "experimental_complex_single_position",
                "netsolp_delta_usability_vs_wt": du,
                "netsolp_usability_magnitude": labels[0],
                "netsolp_delta_solubility_vs_wt": ds,
                "netsolp_solubility_magnitude": labels[1],
                "nanomelt_delta_predicted_tm_c_vs_wt": dt,
                "nanomelt_tm_magnitude": labels[2],
                "property_material_favorable_count": labels.count("favorable"),
                "property_material_adverse_count": labels.count("adverse"),
                "tnp_psh_delta_vs_wt": _float(source, "tnp_psh_delta_vs_wt"),
                "tnp_flag_regression_count": _int(source, "tnp_flag_regression_count"),
                "new_liability_flags": str(source.get("hard_risk_flags", "")),
                "paired_wt_vhh_lost_auth_positions_union": "",
                "paired_wt_vhh_gained_auth_positions_union": "",
                "paired_wt_receptor_lost_auth_positions_union": "",
                "paired_wt_receptor_gained_auth_positions_union": "",
                "pareto_front_within_double_class": "",
                "selection_status": "preliminary_panel",
                "selection_reason": "retain_released_single_for_component_interpretability",
                "selection_order_within_category": 0,
                "final_candidate_selection_performed": False,
            }
        )

    for identifier, source in doubles.items():
        mutation = str(source["mutation_set"])
        _validate_candidate_sequence(parent, str(source["sequence"]), mutation, immutable)
        evidence_class = str(source["joint_evidence_class"])
        category = PRIMARY_DOUBLE_CLASSES.get(evidence_class, "tradeoff_or_unclear_double")
        eligible = (
            evidence_class in PRIMARY_DOUBLE_CLASSES
            and str(source["pyrosetta_structural_safety_status"]) == "pass"
            and str(source["hard_constraint_status"]) == "pass"
        )
        audit.append(
            {
                "candidate_id": identifier,
                "candidate_kind": "double_mutant",
                "mutation_count": 2,
                "mutation_set": mutation,
                "panel_category": category,
                "source_evidence_class": evidence_class,
                "sequence": str(source["sequence"]),
                "primary_pool_eligible": eligible,
                "primary_pool_reason": "v2_1_supported_double" if eligible else "v2_1_tradeoff_or_unclear",
                "hard_constraint_status": str(source["hard_constraint_status"]),
                "hard_constraint_blockers": str(source.get("pyrosetta_structural_safety_blockers", "")),
                "expert_risk_level": "contact_review" if str(source["paired_contact_change_status"]) == "changed" else "no_contact_change",
                "risk_flags": str(source.get("new_liability_flags", "")),
                "position_pair": str(source["position_pair"]),
                "component_mutations": mutation,
                "pyrosetta_protocol_family": "paired_local_repack_3_replicate_double_mutant",
                "pyrosetta_replicate_count": _int(source, "pyrosetta_replicate_count"),
                "pyrosetta_negative_dg_count": "",
                "pyrosetta_negative_cross_count": "",
                "pyrosetta_delta_dg_median": _float(source, "pyrosetta_delta_dG_separated_median"),
                "pyrosetta_delta_cross_median": _float(source, "pyrosetta_delta_cross_interface_energy_median"),
                "pyrosetta_contact_change_status": str(source["paired_contact_change_status"]),
                "paired_contact_changed_replicate_count": _int(source, "paired_contact_changed_replicate_count"),
                "experimental_reference_sensitivity_status": str(source["experimental_reference_sensitivity_status"]),
                "antifold_delta_log_probability": _float(source, "antifold_additive_fixed_backbone_delta_log_probability"),
                "antifold_evidence_scope": "additive_fixed_backbone_single_position_identity",
                "netsolp_delta_usability_vs_wt": _float(source, "netsolp_delta_usability_vs_wt"),
                "netsolp_usability_magnitude": str(source["netsolp_usability_magnitude"]),
                "netsolp_delta_solubility_vs_wt": _float(source, "netsolp_delta_solubility_vs_wt"),
                "netsolp_solubility_magnitude": str(source["netsolp_solubility_magnitude"]),
                "nanomelt_delta_predicted_tm_c_vs_wt": _float(source, "nanomelt_delta_predicted_apparent_tm_c_vs_wt"),
                "nanomelt_tm_magnitude": str(source["nanomelt_tm_magnitude"]),
                "property_material_favorable_count": _int(source, "property_material_favorable_count"),
                "property_material_adverse_count": _int(source, "property_material_adverse_count"),
                "tnp_psh_delta_vs_wt": _float(source, "tnp_psh_delta_vs_wt"),
                "tnp_flag_regression_count": _int(source, "tnp_flag_regression_count"),
                "new_liability_flags": str(source.get("new_liability_flags", "")),
                "paired_wt_vhh_lost_auth_positions_union": str(source.get("paired_wt_vhh_lost_auth_positions_union", "")),
                "paired_wt_vhh_gained_auth_positions_union": str(source.get("paired_wt_vhh_gained_auth_positions_union", "")),
                "paired_wt_receptor_lost_auth_positions_union": str(source.get("paired_wt_receptor_lost_auth_positions_union", "")),
                "paired_wt_receptor_gained_auth_positions_union": str(source.get("paired_wt_receptor_gained_auth_positions_union", "")),
                "pareto_front_within_double_class": "",
                "selection_status": "not_selected_current_round",
                "selection_reason": "not_in_preliminary_panel_or_reserve",
                "selection_order_within_category": "",
                "final_candidate_selection_performed": False,
            }
        )

    by_category: dict[str, list[dict[str, object]]] = {}
    for row in audit:
        if row["candidate_kind"] == "double_mutant" and row["primary_pool_eligible"]:
            by_category.setdefault(str(row["panel_category"]), []).append(row)
    expected_double_counts = {
        "balanced_combination": 28,
        "affinity_supported_double": 12,
        "property_supported_double": 2,
    }
    if {key: len(value) for key, value in by_category.items()} != expected_double_counts:
        raise PreliminaryPanelError("Unexpected V2.1 supported-double class counts")

    for rows in by_category.values():
        _assign_pareto_fronts(rows)

    selected_balanced = _select_diverse_balanced(
        by_category["balanced_combination"], MAIN_QUOTAS["balanced_combination"]
    )
    selected_ids = {str(row["candidate_id"]) for row in selected_balanced}
    reserve_balanced = _diversity_order(
        [row for row in by_category["balanced_combination"] if str(row["candidate_id"]) not in selected_ids]
    )[: RESERVE_QUOTAS["balanced_combination"]]
    reserve_affinity = _diversity_order(by_category["affinity_supported_double"])[
        : RESERVE_QUOTAS["affinity_supported_double"]
    ]
    reserve_property = _diversity_order(by_category["property_supported_double"])[
        : RESERVE_QUOTAS["property_supported_double"]
    ]

    for category in ("affinity_focused_single", "property_focused_single"):
        rows = sorted(
            (row for row in audit if row["panel_category"] == category),
            key=lambda row: str(row["candidate_id"]),
        )
        if len(rows) != MAIN_QUOTAS[category]:
            raise PreliminaryPanelError(f"Unexpected single count for {category}")
        for order, row in enumerate(rows, start=1):
            row["selection_order_within_category"] = order

    for order, row in enumerate(selected_balanced, start=1):
        row["selection_status"] = "preliminary_panel"
        row["selection_reason"] = "balanced_double_selected_by_pareto_then_mutation_diversity"
        row["selection_order_within_category"] = order
    reserves = [*reserve_balanced, *reserve_affinity, *reserve_property]
    for category in RESERVE_QUOTAS:
        category_rows = [row for row in reserves if row["panel_category"] == category]
        for order, row in enumerate(category_rows, start=1):
            row["selection_status"] = "reserve"
            row["selection_reason"] = "class_stratified_reserve_after_preliminary_panel"
            row["selection_order_within_category"] = order

    panel = sorted(
        (row for row in audit if row["selection_status"] == "preliminary_panel"),
        key=lambda row: (
            list(MAIN_QUOTAS).index(str(row["panel_category"])),
            int(row["selection_order_within_category"]),
        ),
    )
    reserves = sorted(
        (row for row in audit if row["selection_status"] == "reserve"),
        key=lambda row: (
            list(RESERVE_QUOTAS).index(str(row["panel_category"])),
            int(row["selection_order_within_category"]),
        ),
    )
    panel_counts = Counter(str(row["panel_category"]) for row in panel)
    reserve_counts = Counter(str(row["panel_category"]) for row in reserves)
    if panel_counts != MAIN_QUOTAS or reserve_counts != RESERVE_QUOTAS:
        raise PreliminaryPanelError("Panel or reserve quotas were not satisfied")
    if len({str(row["sequence"]) for row in panel}) != 30:
        raise PreliminaryPanelError("Preliminary panel sequences are not unique")

    component_counts = Counter(
        token
        for row in selected_balanced
        for token in str(row["mutation_set"]).split(";")
    )
    facts = {
        "reviewed_candidate_count": len(audit),
        "active_single_count": EXPECTED_SINGLES,
        "supported_double_count": sum(expected_double_counts.values()),
        "tradeoff_or_unclear_double_count": EXPECTED_DOUBLES - sum(expected_double_counts.values()),
        "primary_pool_count": sum(bool(row["primary_pool_eligible"]) for row in audit),
        "preliminary_panel_count": len(panel),
        "reserve_count": len(reserves),
        "preliminary_panel_category_counts": dict(panel_counts),
        "reserve_category_counts": dict(reserve_counts),
        "selected_single_count": sum(row["candidate_kind"] == "single_mutant" for row in panel),
        "selected_double_count": sum(row["candidate_kind"] == "double_mutant" for row in panel),
        "selected_double_contact_status_counts": dict(
            Counter(str(row["pyrosetta_contact_change_status"]) for row in selected_balanced)
        ),
        "selected_double_component_counts": dict(sorted(component_counts.items())),
        "selected_double_position_pair_counts": dict(
            sorted(Counter(str(row["position_pair"]) for row in selected_balanced).items())
        ),
        "preliminary_panel_selection_performed": True,
        "final_candidate_selection_performed": False,
    }
    return {"audit_rows": audit, "panel_rows": panel, "reserve_rows": reserves, "facts": facts}


def _validate_gates(single_gate: Mapping[str, object], double_gate: Mapping[str, object]) -> None:
    if single_gate.get("status") != "pass" or single_gate.get("release") != "ready_for_small_combination_contract":
        raise PreliminaryPanelError("Single-mutant shortlist is not released")
    if double_gate.get("status") != "pass" or double_gate.get("release") != "ready_for_scientific_shortlist_definition":
        raise PreliminaryPanelError("Double-mutant V2.1 evidence is not released")
    if int(double_gate.get("candidate_count", -1)) != EXPECTED_DOUBLES:
        raise PreliminaryPanelError("Double-mutant gate candidate count mismatch")


def _parent_contract(stage2: Mapping[str, object]) -> tuple[str, set[int]]:
    if stage2.get("status") != "pass":
        raise PreliminaryPanelError("Stage-2 contract is not released")
    parent_block = stage2.get("authoritative_parent")
    hard = stage2.get("hard_immutable")
    if not isinstance(parent_block, Mapping) or not isinstance(hard, Mapping):
        raise PreliminaryPanelError("Stage-2 parent or immutable contract is missing")
    parent = str(parent_block.get("sequence", ""))
    immutable = {int(value) for value in hard.get("reported_sequence_indices_1based", [])}
    if len(parent) != 128 or parent[-4:] != "SSGS" or immutable != {22, 95, 125, 126, 127, 128}:
        raise PreliminaryPanelError("Authoritative parent or hard immutable set mismatch")
    return parent, immutable


def _validate_candidate_sequence(parent: str, sequence: str, mutation_set: str, immutable: set[int]) -> None:
    if len(sequence) != len(parent) or any(residue not in "ACDEFGHIKLMNPQRSTVWY" for residue in sequence):
        raise PreliminaryPanelError(f"Invalid candidate sequence for {mutation_set}")
    tokens = mutation_set.split(";")
    parsed: list[tuple[str, int, str]] = []
    for token in tokens:
        match = MUTATION_RE.fullmatch(token)
        if match is None:
            raise PreliminaryPanelError(f"Invalid mutation label: {token}")
        wt, position_text, mutant = match.groups()
        position = int(position_text)
        if not 1 <= position <= len(parent) or parent[position - 1] != wt or sequence[position - 1] != mutant:
            raise PreliminaryPanelError(f"Mutation identity mismatch: {token}")
        if position in immutable:
            raise PreliminaryPanelError(f"Mutation touches immutable position: {token}")
        parsed.append((wt, position, mutant))
    if len({position for _, position, _ in parsed}) != len(parsed):
        raise PreliminaryPanelError(f"Repeated mutation position: {mutation_set}")
    expected_differences = {position for _, position, _ in parsed}
    observed_differences = {
        index for index, (wt, mutant) in enumerate(zip(parent, sequence, strict=True), start=1) if wt != mutant
    }
    if observed_differences != expected_differences:
        raise PreliminaryPanelError(f"Sequence differences do not match label: {mutation_set}")
    if sequence.count("C") > parent.count("C"):
        raise PreliminaryPanelError(f"Candidate introduces an unpaired cysteine: {mutation_set}")


def _assign_pareto_fronts(rows: Sequence[dict[str, object]]) -> None:
    remaining = list(rows)
    front = 1
    while remaining:
        current = [row for row in remaining if not any(_dominates(other, row) for other in remaining if other is not row)]
        if not current:
            raise PreliminaryPanelError("Could not resolve double-mutant Pareto fronts")
        for row in current:
            row["pareto_front_within_double_class"] = front
            remaining.remove(row)
        front += 1


def _dominates(left: Mapping[str, object], right: Mapping[str, object]) -> bool:
    left_values = _pareto_values(left)
    right_values = _pareto_values(right)
    return all(a <= b for a, b in zip(left_values, right_values, strict=True)) and any(
        a < b for a, b in zip(left_values, right_values, strict=True)
    )


def _pareto_values(row: Mapping[str, object]) -> tuple[float, float, float, float]:
    values = (
        float(row["pyrosetta_delta_dg_median"]),
        float(row["pyrosetta_delta_cross_median"]),
        -float(row["property_material_favorable_count"]),
        -float(row["antifold_delta_log_probability"]),
    )
    if not all(math.isfinite(value) for value in values):
        raise PreliminaryPanelError("Non-finite Pareto objective")
    return values


def _diversity_order(rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    remaining = list(rows)
    ordered: list[dict[str, object]] = []
    component_use: Counter[str] = Counter()
    pair_use: Counter[str] = Counter()
    while remaining:
        def key(row: Mapping[str, object]) -> tuple[object, ...]:
            components = str(row["mutation_set"]).split(";")
            contact = str(row["pyrosetta_contact_change_status"])
            contact_rank = 0 if contact == "unchanged" else 1
            projected = [component_use[item] + 1 for item in components]
            return (
                int(row["pareto_front_within_double_class"]),
                max(projected),
                sum(projected),
                pair_use[str(row["position_pair"])],
                contact_rank,
                -int(row["property_material_favorable_count"]),
                max(float(row["pyrosetta_delta_dg_median"]), float(row["pyrosetta_delta_cross_median"])),
                str(row["candidate_id"]),
            )

        chosen = min(remaining, key=key)
        ordered.append(chosen)
        remaining.remove(chosen)
        for component in str(chosen["mutation_set"]).split(";"):
            component_use[component] += 1
        pair_use[str(chosen["position_pair"])] += 1
    return ordered


def _select_diverse_balanced(
    rows: Sequence[dict[str, object]], count: int
) -> list[dict[str, object]]:
    """Select balanced doubles with explicit, feasibility-checked diversity caps.

    No component mutation may occupy more than five of the sixteen slots and
    no reported position pair may occupy more than two. These caps prevent a
    single R45 or 30/45 hypothesis from consuming the panel while still
    allowing repeated substitutions where the evidence is strong.
    """

    remaining = list(rows)
    selected: list[dict[str, object]] = []
    component_use: Counter[str] = Counter()
    pair_use: Counter[str] = Counter()
    while len(selected) < count:
        eligible = []
        for row in remaining:
            components = str(row["mutation_set"]).split(";")
            if any(component_use[item] >= 5 for item in components):
                continue
            if pair_use[str(row["position_pair"])] >= 2:
                continue
            eligible.append(row)
        if not eligible:
            raise PreliminaryPanelError("Balanced-double diversity caps cannot fill the panel")

        def key(row: Mapping[str, object]) -> tuple[object, ...]:
            components = str(row["mutation_set"]).split(";")
            projected = [component_use[item] + 1 for item in components]
            contact_rank = 0 if str(row["pyrosetta_contact_change_status"]) == "unchanged" else 1
            return (
                int(row["pareto_front_within_double_class"]),
                max(projected),
                sum(projected),
                pair_use[str(row["position_pair"])],
                contact_rank,
                -int(row["property_material_favorable_count"]),
                max(float(row["pyrosetta_delta_dg_median"]), float(row["pyrosetta_delta_cross_median"])),
                str(row["candidate_id"]),
            )

        chosen = min(eligible, key=key)
        selected.append(chosen); remaining.remove(chosen)
        for component in str(chosen["mutation_set"]).split(";"):
            component_use[component] += 1
        pair_use[str(chosen["position_pair"])] += 1
    return selected


def _property_labels(du: float, ds: float, dt: float) -> tuple[str, str, str]:
    return (
        magnitude_label(du, MAGNITUDE_THRESHOLDS["netsolp_delta_usability_vs_wt"]),
        magnitude_label(ds, MAGNITUDE_THRESHOLDS["netsolp_delta_solubility_vs_wt"]),
        magnitude_label(dt, MAGNITUDE_THRESHOLDS["nanomelt_delta_predicted_apparent_tm_c_vs_wt"]),
    )


def _unique(rows: Sequence[Mapping[str, object]], key: str, expected: int, label: str) -> dict[str, Mapping[str, object]]:
    result = _unique_any(rows, key, label)
    if len(result) != expected:
        raise PreliminaryPanelError(f"Expected {expected} {label}, found {len(result)}")
    return result


def _unique_any(rows: Sequence[Mapping[str, object]], key: str, label: str) -> dict[str, Mapping[str, object]]:
    result: dict[str, Mapping[str, object]] = {}
    for row in rows:
        identifier = str(row.get(key, ""))
        if not identifier or identifier in result:
            raise PreliminaryPanelError(f"Missing or duplicate {label} identifier: {identifier}")
        result[identifier] = row
    return result


def _tokens(value: object) -> set[str]:
    return {token for token in str(value or "").split(";") if token}


def _float(row: Mapping[str, object], field: str) -> float:
    value = float(row[field])
    if not math.isfinite(value):
        raise PreliminaryPanelError(f"Non-finite {field}")
    return value


def _int(row: Mapping[str, object], field: str) -> int:
    return int(row[field])
