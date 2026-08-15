"""Unified Nb252 single-mutant space and AntiFold landscape joins.

The module enumerates every non-WT amino acid at each non-immutable reported
sequence position.  It records hard constraints and lightweight sequence
liability deltas without selecting candidates.  AntiFold joins reuse existing
per-position WT-backbone logits; they do not rerun inverse folding or interpret
compatibility as affinity, stability, expression, or measured yield.
"""

from __future__ import annotations

import re
from collections import Counter
from typing import Mapping, Sequence

from antibody_optimization.antifold_validation import AA_COLUMNS


class UnifiedSingleMutantError(ValueError):
    """Raised when a unified-design input or result is inconsistent."""


HYDROPHOBIC = frozenset("AVILMFWY")
RESIDUE_CLASS = {
    **{aa: "hydrophobic" for aa in "AVILMFWY"},
    **{aa: "polar" for aa in "STNQ"},
    **{aa: "positive" for aa in "KRH"},
    **{aa: "negative" for aa in "DE"},
    **{aa: "special" for aa in "CGP"},
}


def build_unified_space(
    stage2: Mapping[str, object],
    position_rows: Sequence[Mapping[str, object]],
    critical: Mapping[str, object],
    core_rows: Sequence[Mapping[str, object]],
    affinity_candidate_rows: Sequence[Mapping[str, object]],
    affinity_summary_rows: Sequence[Mapping[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Build position and complete single-substitution tables.

    The input position table must contain one row for every reported parent
    residue.  Missing experimental coordinates remain enumerable but deferred;
    immutable positions have no candidates; new cysteine substitutions remain
    visible but blocked.  No FR/CDR or AntiFold score threshold is applied.
    """

    parent = stage2.get("authoritative_parent", {})
    if not isinstance(parent, Mapping):
        raise UnifiedSingleMutantError("Missing authoritative parent")
    sequence = str(parent.get("sequence", ""))
    if len(sequence) != 128 or parent.get("sequence_sha256") != critical.get("authoritative_parent", {}).get("sequence_sha256"):
        raise UnifiedSingleMutantError("Authoritative parent identity mismatch")
    if stage2.get("status") != "pass" or critical.get("status") != "pass":
        raise UnifiedSingleMutantError("Released upstream contracts must pass")

    hard = set(map(int, stage2["hard_immutable"]["reported_sequence_indices_1based"]))
    missing = set(map(int, critical["experimental_missing_coordinates"]["reported_sequence_indices_1based"]))
    interface = set(map(int, critical["reproduced_experimental_interface"]["reported_sequence_indices_1based"]))
    by_position = {int(row["sequence_index_1based"]): row for row in position_rows}
    if set(by_position) != set(range(1, 129)):
        raise UnifiedSingleMutantError("Position inventory must contain exactly reported positions 1-128")

    core_by_key: dict[tuple[int, str], Mapping[str, object]] = {}
    for row in core_rows:
        if not _bool(row.get("core_module_selected")):
            continue
        key = (int(row["sequence_index_1based"]), str(row["mutant_residue"]))
        if key in core_by_key:
            raise UnifiedSingleMutantError(f"Duplicate affinity core substitution: {key}")
        core_by_key[key] = row
    if len(core_by_key) != 8:
        raise UnifiedSingleMutantError("Expected eight released affinity-core substitutions")

    affinity_by_key: dict[tuple[int, str], Mapping[str, object]] = {}
    for row in affinity_candidate_rows:
        key = (int(row["sequence_index_1based"]), str(row["mutant_residue"]))
        if key in affinity_by_key:
            raise UnifiedSingleMutantError(f"Duplicate existing affinity substitution: {key}")
        affinity_by_key[key] = row
    if len(affinity_by_key) != 456:
        raise UnifiedSingleMutantError("Expected the complete existing 456-member interface scan")
    summary_by_id = {str(row["candidate_id"]): row for row in affinity_summary_rows}
    summary_ids = set(summary_by_id)
    if len(summary_ids) != 456 or summary_ids != {str(row["candidate_id"]) for row in affinity_candidate_rows}:
        raise UnifiedSingleMutantError("Existing PyRosetta summary does not exactly cover the 456 interface mutants")

    positions: list[dict[str, object]] = []
    candidates: list[dict[str, object]] = []
    for index in range(1, 129):
        source = by_position[index]
        wt = sequence[index - 1]
        if str(source["residue_aa"]) != wt:
            raise UnifiedSingleMutantError(f"WT mismatch at reported position {index}")
        immutable = index in hard
        missing_coordinates = index in missing
        interface_position = index in interface
        positions.append({
            "sample_uid": parent["sample_uid"],
            "sequence_index_1based": index,
            "wt_residue": wt,
            "numbering_scheme": str(source["numbering_scheme"]),
            "numbering_position_label": str(source["numbering_position_label"]),
            "region": str(source["region"]),
            "experimental_coordinate_status": str(source["experimental_coordinate_status"]),
            "experimental_missing_coordinates": missing_coordinates,
            "experimental_interface": interface_position,
            "interface_mutation_semantics": "cautious_not_forbidden" if interface_position else "not_interface_annotated",
            "hard_immutable": immutable,
            "hard_immutable_reasons": str(source.get("hard_immutable_reasons", "")),
            "enumerated_substitution_count": 0 if immutable else 19,
            "design_scope": "hard_immutable" if immutable else "complete_non_wt_single_substitution_space",
        })
        if immutable:
            continue
        for mutant in AA_COLUMNS:
            if mutant == wt:
                continue
            key = (index, mutant)
            core = core_by_key.get(key)
            affinity = affinity_by_key.get(key)
            affinity_summary = summary_by_id.get(str(affinity["candidate_id"])) if affinity else None
            if interface_position and affinity is None:
                raise UnifiedSingleMutantError(f"Interface substitution lacks existing PyRosetta identity: {key}")
            if not interface_position and affinity is not None:
                raise UnifiedSingleMutantError(f"Non-interface substitution unexpectedly appears in affinity scan: {key}")
            status = (
                "blocked_new_unpaired_cys" if mutant == "C"
                else "deferred_missing_experimental_coordinates" if missing_coordinates
                else "eligible_current_round"
            )
            mutated = sequence[: index - 1] + mutant + sequence[index:]
            liabilities = _liability_deltas(sequence, mutated)
            candidates.append({
                "candidate_id": str(affinity["candidate_id"]) if affinity else f"Nb252_uni_seq{index:03d}_{wt}{index}{mutant}",
                "parent_sample_uid": parent["sample_uid"],
                "sequence_index_1based": index,
                "wt_residue": wt,
                "mutant_residue": mutant,
                "mutation_reported_label": f"Nb252 reported_seq {wt}{index}{mutant}",
                "numbering_scheme": str(source["numbering_scheme"]),
                "numbering_position_label": str(source["numbering_position_label"]),
                "region": str(source["region"]),
                "experimental_missing_coordinates": missing_coordinates,
                "experimental_interface": interface_position,
                "interface_mutation_semantics": "cautious_not_forbidden" if interface_position else "not_interface_annotated",
                "design_status": status,
                "design_status_reason": {
                    "blocked_new_unpaired_cys": "introduces_extra_unpaired_cysteine",
                    "deferred_missing_experimental_coordinates": "experimental_structure_not_evaluable_current_round",
                    "eligible_current_round": "passes_unified_single_mutant_hard_constraints",
                }[status],
                "design_track": (
                    "affinity_existing_interface_scan" if interface_position
                    else "audit_only_missing_coordinates" if missing_coordinates
                    else "stability_developability_discovery"
                ),
                "pyrosetta_evidence_status": (
                    "reuse_existing_three_replicate_full_scan" if interface_position
                    else "not_evaluable_missing_experimental_coordinates" if missing_coordinates
                    else "not_run_shortlist_before_affinity_noninferiority_check"
                ),
                "pyrosetta_rescoring_required_now": False,
                "existing_affinity_scan_candidate": affinity is not None,
                "existing_affinity_scan_candidate_id": str(affinity["candidate_id"]) if affinity else "",
                "pyrosetta_full_scan_status": str(affinity_summary["status"]) if affinity_summary else "not_run",
                "pyrosetta_replicate_count": int(affinity_summary["replicate_count"]) if affinity_summary else "",
                "pyrosetta_delta_dG_separated_median": float(affinity_summary["delta_dG_separated_median"]) if affinity_summary else "",
                "pyrosetta_delta_cross_interface_energy_median": float(affinity_summary["delta_cross_interface_energy_median"]) if affinity_summary else "",
                "pyrosetta_delta_interface_fa_rep_median": float(affinity_summary["delta_interface_fa_rep_median"]) if affinity_summary else "",
                "pyrosetta_minimum_paired_wt_vhh_contact_retention": float(affinity_summary["minimum_candidate_vs_paired_wt_vhh_contact_retention"]) if affinity_summary else "",
                "pyrosetta_minimum_paired_wt_receptor_epitope_retention": float(affinity_summary["minimum_candidate_vs_paired_wt_receptor_epitope_retention"]) if affinity_summary else "",
                "affinity_core_module": core is not None,
                "affinity_core_source_tier": str(core.get("source_tier", "")) if core else "",
                "affinity_core_risk_flags": str(core.get("risk_flags", "")) if core else "",
                "wt_residue_class": RESIDUE_CLASS[wt],
                "mutant_residue_class": RESIDUE_CLASS[mutant],
                "formal_charge_delta": _formal_charge(mutant) - _formal_charge(wt),
                "hydrophobic_fraction_delta": round((int(mutant in HYDROPHOBIC) - int(wt in HYDROPHOBIC)) / len(sequence), 8),
                **liabilities,
                "sequence": mutated,
                "candidate_selection_performed": False,
            })
    _validate_space(sequence, positions, candidates, hard, missing)
    return positions, candidates


def evaluate_antifold_landscape(
    candidates: Sequence[Mapping[str, object]],
    indexed_views: Mapping[str, Mapping[str, Mapping[str, object]]],
) -> tuple[list[dict[str, object]], dict[str, object]]:
    """Join every candidate to three existing AntiFold per-position tables."""

    required = ("experimental_vhh_only", "experimental_complex_context", "af3_vhh_only")
    if tuple(indexed_views) != required:
        raise UnifiedSingleMutantError(f"Expected AntiFold views {required}")
    output: list[dict[str, object]] = []
    for candidate in candidates:
        label = str(candidate["numbering_position_label"])
        wt = str(candidate["wt_residue"])
        mutant = str(candidate["mutant_residue"])
        row = dict(candidate)
        evaluable = 0
        directions: list[str] = []
        for view in required:
            score = indexed_views[view].get(label)
            prefix = view
            if score is None:
                row[f"{prefix}_evaluation_status"] = "not_evaluable"
                row[f"{prefix}_wt_log_probability"] = ""
                row[f"{prefix}_mutant_log_probability"] = ""
                row[f"{prefix}_delta_log_probability"] = ""
                row[f"{prefix}_perplexity"] = ""
                row[f"{prefix}_direction"] = "not_evaluable"
                directions.append("not_evaluable")
                continue
            if str(score["pdb_res"]) != wt:
                raise UnifiedSingleMutantError(f"AntiFold WT mismatch for {candidate['candidate_id']} in {view}")
            wt_logp = float(score[wt]); mutant_logp = float(score[mutant]); delta = mutant_logp - wt_logp
            direction = "positive" if delta > 0 else "negative" if delta < 0 else "zero"
            row[f"{prefix}_evaluation_status"] = "pass"
            row[f"{prefix}_wt_log_probability"] = wt_logp
            row[f"{prefix}_mutant_log_probability"] = mutant_logp
            row[f"{prefix}_delta_log_probability"] = delta
            row[f"{prefix}_perplexity"] = float(score["perplexity"])
            row[f"{prefix}_direction"] = direction
            evaluable += 1; directions.append(direction)
        row["antifold_evaluable_view_count"] = evaluable
        row["antifold_evaluation_scope"] = "three_views" if evaluable == 3 else "af3_only" if evaluable == 1 else "partial"
        row["all_view_directions_concordant"] = evaluable == 3 and len(set(directions)) == 1
        row["antifold_candidate_filtering_applied"] = False
        output.append(row)

    counts = Counter(str(row["antifold_evaluation_scope"]) for row in output)
    eligible = [row for row in output if row["design_status"] == "eligible_current_round"]
    passed = (
        len(output) == 2318 and counts == {"three_views": 2071, "af3_only": 247}
        and len(eligible) == 1962 and all(row["antifold_evaluation_scope"] == "three_views" for row in eligible)
    )
    gate = {
        "schema_version": 1,
        "gate_name": "nb252_unified_single_mutant_antifold_landscape",
        "status": "pass" if passed else "blocked",
        "candidate_count": len(output),
        "evaluation_scope_counts": dict(sorted(counts.items())),
        "eligible_current_round_count": len(eligible),
        "eligible_current_round_all_three_views": all(row["antifold_evaluation_scope"] == "three_views" for row in eligible),
        "candidate_selection_performed": False,
        "scientific_score_threshold_applied": False,
        "release": "ready_for_unified_property_scoring" if passed else "blocked",
        "interpretation": "AntiFold structure-conditioned compatibility only; not affinity, stability, expression, yield, or experimental validation.",
    }
    return output, gate


def _validate_space(sequence: str, positions: Sequence[Mapping[str, object]], candidates: Sequence[Mapping[str, object]], hard: set[int], missing: set[int]) -> None:
    if len(positions) != 128 or len(candidates) != 2318:
        raise UnifiedSingleMutantError("Unexpected unified-space size")
    counts = Counter(str(row["design_status"]) for row in candidates)
    if counts != {"eligible_current_round": 1962, "deferred_missing_experimental_coordinates": 234, "blocked_new_unpaired_cys": 122}:
        raise UnifiedSingleMutantError(f"Unexpected design-status counts: {counts}")
    for row in candidates:
        index = int(row["sequence_index_1based"])
        candidate_sequence = str(row["sequence"])
        differences = [i for i, (left, right) in enumerate(zip(sequence, candidate_sequence, strict=True), 1) if left != right]
        if differences != [index] or index in hard:
            raise UnifiedSingleMutantError(f"Invalid single mutant: {row['candidate_id']}")
        if index in missing and row["mutant_residue"] != "C" and row["design_status"] != "deferred_missing_experimental_coordinates":
            raise UnifiedSingleMutantError("Missing-coordinate candidate was not deferred")


def _liability_deltas(wt: str, mutant: str) -> dict[str, object]:
    patterns = {
        "n_linked_glycosylation_motif": r"N[^P][ST]",
        "deamidation_motif": r"N[GST]",
        "isomerization_motif": r"D[GST]",
    }
    values: dict[str, object] = {}
    flags: list[str] = []
    for name, pattern in patterns.items():
        before = len(re.findall(f"(?=({pattern}))", wt))
        after = len(re.findall(f"(?=({pattern}))", mutant))
        delta = after - before
        values[f"{name}_delta"] = delta
        if delta > 0:
            flags.append(f"new_{name}")
    oxidation_delta = sum(mutant.count(aa) - wt.count(aa) for aa in "MW")
    values["oxidation_susceptible_residue_delta"] = oxidation_delta
    if oxidation_delta > 0:
        flags.append("more_M_or_W")
    values["new_liability_flags"] = "|".join(flags)
    return values


def _formal_charge(aa: str) -> int:
    return -1 if aa in "DE" else 1 if aa in "KR" else 0


def _bool(value: object) -> bool:
    return value is True or str(value).lower() == "true"
