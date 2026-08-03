"""Evaluate stage-1 engineering and scientific release gates.

This module consumes explicit sequence, expression, structure, interface, and
input-freeze evidence. It never converts a blocked or absent upstream result
into a pass and reports release-specific blocker lists.
"""

from __future__ import annotations

from typing import Mapping, Sequence


VALID_GATE_STATUSES = {"pass", "blocked", "pending", "not_applicable"}


def evaluate_stage1_gates(
    *,
    numbering_review: Sequence[Mapping[str, str]],
    sample_comparability: Sequence[Mapping[str, str]],
    input_freeze_manifest: Mapping[str, object] | None,
    structure_manifest: Mapping[str, object] | None,
    interface_manifest: Mapping[str, object] | None,
) -> dict[str, object]:
    """Evaluate engineering completion separately from scientific release.

    A failed ANARCII classification does not fail the numbering inventory gate
    when the failure is retained explicitly.  In contrast, missing structural
    exports, unconfirmed roles, or an unconfirmed construct block candidate
    design.  Pooled-expression release is read only from explicit audit facts.
    """

    numbered_ids = [row.get("sample_uid", "") for row in numbering_review]
    reviewed_ids = [row.get("sample_uid", "") for row in sample_comparability]
    unique_numbered = set(numbered_ids)
    unique_reviewed = set(reviewed_ids)
    row_identity_ok = (
        len(numbering_review) == 47
        and len(unique_numbered) == 47
        and len(sample_comparability) == 47
        and len(unique_reviewed) == 47
        and all(numbered_ids)
        and all(reviewed_ids)
        and unique_numbered == unique_reviewed
    )
    input_ok = (
        row_identity_ok
        and _manifest_status(input_freeze_manifest, "status") == "pass"
    )
    numbering_closed = row_identity_ok and all(
        _first(row, "numbering_status", "status") in {"pass", "failed"}
        for row in numbering_review
    )
    expression_closed = row_identity_ok and all(
        row.get("cross_assay_pooling_status") in {"blocked", "pass"}
        and row.get("nb252_transfer_status") in {"blocked", "pass"}
        for row in sample_comparability
    )

    structure_export = _manifest_status(structure_manifest, "export_status", "status")
    structure_inventory = _manifest_status(
        structure_manifest, "inventory_status", "structure_inventory_status"
    )
    structure_identity = _manifest_status(
        structure_manifest, "chain_role_status", "structure_identity_status"
    )
    residue_mapping = _manifest_status(
        structure_manifest, "residue_mapping_status", "mapping_status"
    )
    authoritative_sequence = _manifest_status(
        structure_manifest, "authoritative_sequence_status"
    )
    interface_status = _manifest_status(
        interface_manifest, "interface_status", "status"
    )
    orange_status = _manifest_status(
        interface_manifest, "orange_annotation_status"
    )

    gates = {
        "input_integrity": _gate(input_ok, "Frozen 47-sequence identity is consistent"),
        "sequence_numbering_inventory": _gate(
            numbering_closed,
            "All 47 ANARCII outcomes are retained, including explicit failures",
        ),
        "expression_audit": _gate(
            expression_closed,
            "All sample-level use decisions are closed or explicitly blocked",
        ),
        "structure_export": _gate_from_status(
            structure_export,
            "Verified ChimeraX export is required",
            structure_manifest,
        ),
        "structure_inventory": _gate_from_status(
            structure_inventory,
            "Gemmi read-back inventory is required",
            structure_manifest,
        ),
        "structure_identity": _gate_from_status(
            structure_identity,
            "VHH and NK2R chain roles require explicit review",
            structure_manifest,
        ),
        "residue_mapping": _gate_from_status(
            residue_mapping,
            "Unique reversible sequence/IMGT/structure mapping is required",
            structure_manifest,
        ),
        "authoritative_nb252_sequence": _gate_from_status(
            authoritative_sequence,
            "The reported Nb252 sequence remains provisional until construct confirmation",
            structure_manifest,
        ),
        "interface_safety": _combined_gate(
            interface_status,
            orange_status,
            "Both the temporary <4 A interface and orange annotation require confirmation",
            interface_manifest,
        ),
        "cross_assay_pooling": _gate(
            expression_closed
            and all(
                row.get("cross_assay_pooling_status") == "pass"
                for row in sample_comparability
            ),
            "All 47 samples require an explicit pass for cross-assay pooling",
        ),
        "nb252_expression_transfer": _gate(
            expression_closed
            and all(
                row.get("nb252_transfer_status") == "pass"
                for row in sample_comparability
            ),
            "All 47 samples require an explicit pass before Nb252 expression transfer",
        ),
    }

    local_required = (
        "input_integrity",
        "sequence_numbering_inventory",
        "expression_audit",
        "structure_export",
        "structure_inventory",
        "structure_identity",
        "residue_mapping",
        "interface_safety",
    )
    local_pass = all(gates[name]["status"] == "pass" for name in local_required)
    candidate_required = (
        "authoritative_nb252_sequence",
        "structure_identity",
        "residue_mapping",
        "interface_safety",
    )
    candidate_pass = local_pass and all(
        gates[name]["status"] == "pass" for name in candidate_required
    )
    pooled_required = (
        "input_integrity",
        "expression_audit",
        "cross_assay_pooling",
    )
    pooled_pass = all(gates[name]["status"] == "pass" for name in pooled_required)

    blockers = [
        name for name, value in gates.items() if value["status"] != "pass"
    ]
    local_blockers = [
        name for name in local_required if gates[name]["status"] != "pass"
    ]
    candidate_dependencies = tuple(dict.fromkeys((*local_required, *candidate_required)))
    candidate_blockers = [
        name for name in candidate_dependencies if gates[name]["status"] != "pass"
    ]
    pooled_blockers = [
        name for name in pooled_required if gates[name]["status"] != "pass"
    ]
    return {
        "schema_version": 1,
        "local_baseline_build": "pass" if local_pass else "blocked",
        "candidate_design_release": "pass" if candidate_pass else "blocked",
        "pooled_expression_model_release": "pass" if pooled_pass else "blocked",
        "gates": gates,
        "blocking_gates": blockers,
        "local_baseline_build_blockers": local_blockers,
        "candidate_design_release_blockers": candidate_blockers,
        "pooled_expression_model_release_blockers": pooled_blockers,
        "interpretation": (
            "A blocked scientific release is an evidence state, not an implementation failure."
        ),
    }


def structure_evidence_is_verified(
    manifest: Mapping[str, object] | None,
) -> bool:
    """Return true only for a passed export, inventory, identity, and mapping."""

    return all(
        _manifest_status(manifest, *keys) == "pass"
        for keys in (
            ("export_status", "status"),
            ("inventory_status", "structure_inventory_status"),
            ("chain_role_status", "structure_identity_status"),
            ("residue_mapping_status", "mapping_status"),
        )
    )


def interface_evidence_statuses(
    manifest: Mapping[str, object] | None,
) -> tuple[bool, bool]:
    """Return separate verified flags for distance contacts and orange labels."""

    return (
        _manifest_status(manifest, "interface_status", "status") == "pass",
        _manifest_status(manifest, "orange_annotation_status") == "pass",
    )


def _first(row: Mapping[str, object], *keys: str) -> str:
    for key in keys:
        value = row.get(key)
        if value is not None and str(value) != "":
            return str(value)
    return ""


def _manifest_status(manifest: Mapping[str, object] | None, *keys: str) -> str:
    if manifest is None:
        return "blocked"
    for key in keys:
        value = manifest.get(key)
        if isinstance(value, Mapping):
            value = value.get("status")
        if isinstance(value, str) and value in VALID_GATE_STATUSES:
            return value
    gates = manifest.get("gates")
    if isinstance(gates, Mapping):
        for key in keys:
            value = gates.get(key)
            if isinstance(value, Mapping):
                value = value.get("status")
            if isinstance(value, str) and value in VALID_GATE_STATUSES:
                return value
    return "blocked"


def _gate(passed: bool, reason: str) -> dict[str, str]:
    return {"status": "pass" if passed else "blocked", "reason": reason}


def _gate_from_status(
    status: str,
    reason: str,
    manifest: Mapping[str, object] | None = None,
) -> dict[str, str]:
    details = _manifest_blockers(manifest)
    if manifest is None:
        details = ["upstream manifest absent"]
    suffix = f" Observed={status}."
    if details:
        suffix += " Upstream blockers: " + "; ".join(details)
    return {
        "status": "pass" if status == "pass" else "blocked",
        "reason": reason + suffix,
    }


def _combined_gate(
    first: str,
    second: str,
    reason: str,
    manifest: Mapping[str, object] | None = None,
) -> dict[str, str]:
    details = _manifest_blockers(manifest)
    if manifest is None:
        details = ["upstream manifest absent"]
    suffix = f" Observed interface={first}, orange={second}."
    if details:
        suffix += " Upstream blockers: " + "; ".join(details)
    return {
        "status": "pass" if first == second == "pass" else "blocked",
        "reason": reason + suffix,
    }


def _manifest_blockers(manifest: Mapping[str, object] | None) -> list[str]:
    if manifest is None:
        return []
    value = manifest.get("blockers")
    if not isinstance(value, list):
        return []
    return [str(item) for item in value if str(item)]
