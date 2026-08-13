"""Build the position-level WT discovery contract for stability/expression design."""

from __future__ import annotations

from collections import Counter
from typing import Mapping, Sequence


EXPECTED_POSITIONS = 128


class StabilityExpressionContractError(ValueError):
    """Raised when the design-position contract is internally inconsistent."""


def build_stability_expression_contract(
    position_rows: Sequence[Mapping[str, object]],
    stage2_contract: Mapping[str, object],
) -> dict[str, object]:
    """Classify all 128 positions without generating mutations.

    WT module discovery is limited to non-interface framework positions.
    Framework positions missing from the experimental model remain cautiously
    designable only with full-VHH predicted-structure evidence.  CDRs,
    experimental interface positions, disulfides, and terminal SSGS are frozen.
    """

    if stage2_contract.get("status") != "pass":
        raise StabilityExpressionContractError("Stage-2 contract is not released")
    if len(position_rows) != EXPECTED_POSITIONS:
        raise StabilityExpressionContractError("Position inventory must contain 128 rows")
    indices = [int(row["sequence_index_1based"]) for row in position_rows]
    if indices != list(range(1, EXPECTED_POSITIONS + 1)):
        raise StabilityExpressionContractError("Position inventory is not contiguous")

    output = []
    for source in position_rows:
        region = str(source["region"])
        hard = _as_bool(source["hard_immutable"])
        interface = _as_bool(source["experimental_interface"])
        coordinate = str(source["experimental_coordinate_status"])
        if hard:
            status, reason, evidence = "frozen_hard", str(source["hard_immutable_reasons"]), "experimental_or_construct_constraint"
        elif interface:
            status, reason, evidence = "frozen_interface", "preserve_experimental_epitope_and_binding_pose", "experimental_complex"
        elif not region.startswith("FR"):
            status, reason, evidence = "frozen_cdr_or_flank", "WT_stability_expression_discovery_restricted_to_framework", "numbering_region"
        elif coordinate == "observed":
            status, reason, evidence = "allowed_observed_framework", "noninterface_framework_with_experimental_coordinates", "experimental_VHH_coordinates"
        elif coordinate == "missing_coordinates":
            status, reason, evidence = "allowed_cautious_predicted_framework", "noninterface_framework_missing_experimental_coordinates", "AF3_full_VHH_required"
        else:
            status, reason, evidence = "frozen_unresolved", "unsupported_coordinate_state", "none"
        designable = status.startswith("allowed_")
        output.append(
            {
                "sequence_index_1based": int(source["sequence_index_1based"]),
                "sequence_index_0based": int(source["sequence_index_0based"]),
                "residue_aa": source["residue_aa"],
                "numbering_scheme": source["numbering_scheme"],
                "numbering_position_label": source["numbering_position_label"],
                "region": region or "terminal_flank",
                "experimental_coordinate_status": coordinate,
                "experimental_interface": interface,
                "hard_immutable": hard,
                "wt_discovery_status": status,
                "wt_discovery_designable": designable,
                "status_reason": reason,
                "required_structure_evidence": evidence,
                "antifold_sampling_allowed": designable,
                "nanobert_scoring_required_for_generated_sequences": designable,
                "conditional_affinity_background_rule": (
                    "freeze_installed_affinity_mutations_and_recompute_this_position_set"
                    if designable
                    else "remain_frozen"
                ),
            }
        )
    counts = Counter(str(row["wt_discovery_status"]) for row in output)
    expected = {
        "allowed_observed_framework": 72,
        "allowed_cautious_predicted_framework": 9,
        "frozen_interface": 24,
        "frozen_cdr_or_flank": 17,
        "frozen_hard": 6,
    }
    if dict(counts) != expected:
        raise StabilityExpressionContractError(f"Unexpected position counts: {dict(counts)}")
    return {
        "position_rows": output,
        "counts": {
            **expected,
            "designable_positions": sum(bool(row["wt_discovery_designable"]) for row in output),
            "frozen_positions": sum(not bool(row["wt_discovery_designable"]) for row in output),
        },
        "module_rules": {
            "primary_mutations_per_module_minimum": 1,
            "primary_mutations_per_module_maximum": 3,
            "exploratory_mutations_per_module_maximum": 4,
            "generated_mutations": False,
            "antifold_primary_generator": True,
            "abmpnn_parallel_generation": False,
            "nanobert_evidence_level": "pending_yield_association_validation",
        },
    }


def _as_bool(value: object) -> bool:
    return value is True or str(value).lower() == "true"
