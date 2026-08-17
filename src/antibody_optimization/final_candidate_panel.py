"""Validate an explicitly reviewed final Nb252 30-sequence panel.

This module performs identity and decision-contract validation only.  It does
not rank candidates, infer missing review decisions, or run prediction tools.
An explicit, independently reviewed decision contract can be joined to the
energy-review template before the final 30-sequence panel is frozen.
"""

from __future__ import annotations

import re
from collections import Counter
from typing import Mapping, Sequence


class FinalCandidatePanelError(ValueError):
    """Raised when the explicit final-panel review is incomplete or invalid."""


MUTATION_RE = re.compile(r"^([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWY])$")
ALLOWED_DECISIONS = {"select", "reserve", "exclude"}


def apply_explicit_finalist_decisions(
    template_rows: Sequence[Mapping[str, object]],
    decision_contract: Mapping[str, object],
) -> list[dict[str, object]]:
    """Join a complete explicit decision contract to the 36-row review template.

    ``decision_contract`` must contain a ``decisions`` list with one unique row
    per template candidate.  Each row supplies ``candidate_id``,
    ``review_decision`` and a non-empty candidate-specific ``review_rationale``.
    The function does not infer, rank, or fill decisions and requires exactly
    30 selections.  All other template evidence is retained unchanged.
    """

    if len(template_rows) != 36:
        raise FinalCandidatePanelError(
            f"Expected 36 energy-review template rows, found {len(template_rows)}"
        )
    raw_decisions = decision_contract.get("decisions")
    if not isinstance(raw_decisions, list) or len(raw_decisions) != 36:
        count = len(raw_decisions) if isinstance(raw_decisions, list) else 0
        raise FinalCandidatePanelError(
            f"Expected 36 explicit decision rows, found {count}"
        )
    decisions: dict[str, Mapping[str, object]] = {}
    for row in raw_decisions:
        if not isinstance(row, Mapping):
            raise FinalCandidatePanelError("Explicit decision rows must be mappings")
        identifier = str(row.get("candidate_id", "")).strip()
        if not identifier or identifier in decisions:
            raise FinalCandidatePanelError(
                f"Missing or duplicate explicit candidate_id: {identifier}"
            )
        decision = str(row.get("review_decision", "")).strip().lower()
        rationale = str(row.get("review_rationale", "")).strip()
        if decision not in ALLOWED_DECISIONS:
            raise FinalCandidatePanelError(
                f"Invalid explicit decision for {identifier}: {decision}"
            )
        if not rationale:
            raise FinalCandidatePanelError(
                f"Missing explicit rationale for {identifier}"
            )
        decisions[identifier] = row

    template_ids = {str(row["candidate_id"]) for row in template_rows}
    if set(decisions) != template_ids:
        raise FinalCandidatePanelError(
            "Explicit decision candidate identities do not match the review template"
        )
    if sum(
        str(row["review_decision"]).strip().lower() == "select"
        for row in decisions.values()
    ) != 30:
        raise FinalCandidatePanelError("Explicit decision contract must select exactly 30")

    reviewer = str(decision_contract.get("reviewer", "")).strip()
    basis = str(decision_contract.get("decision_basis", "")).strip()
    if not reviewer or not basis:
        raise FinalCandidatePanelError(
            "Explicit decision contract requires reviewer and decision_basis"
        )
    return [
        {
            **dict(row),
            "review_decision": str(decisions[str(row["candidate_id"])]["review_decision"])
            .strip()
            .lower(),
            "review_rationale": str(
                decisions[str(row["candidate_id"])]["review_rationale"]
            ).strip(),
            "reviewer": reviewer,
            "decision_basis": basis,
        }
        for row in template_rows
    ]


def finalize_candidate_panel(
    rows: Sequence[Mapping[str, object]],
    parent_sequence: str,
    expected_candidates: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """Return validated final, reserve, and complete decision-audit rows."""

    if len(rows) != 36 or len(expected_candidates) != 36:
        raise FinalCandidatePanelError(f"Expected 36 reviewed rows, found {len(rows)}")
    expected = {str(row["candidate_id"]): row for row in expected_candidates}
    if len(expected) != 36:
        raise FinalCandidatePanelError("Expected-candidate identities are not unique")
    by_id: dict[str, Mapping[str, object]] = {}
    for row in rows:
        identifier = str(row["candidate_id"])
        if identifier in by_id:
            raise FinalCandidatePanelError(f"Duplicate candidate: {identifier}")
        if identifier not in expected:
            raise FinalCandidatePanelError(f"Unknown reviewed candidate: {identifier}")
        for field in ("mutation_set", "sequence", "candidate_kind", "panel_category"):
            if str(row[field]) != str(expected[identifier][field]):
                raise FinalCandidatePanelError(
                    f"Reviewed identity mismatch for {identifier}: {field}"
                )
        decision = str(row.get("review_decision", "")).strip().lower()
        rationale = str(row.get("review_rationale", "")).strip()
        if decision not in ALLOWED_DECISIONS:
            raise FinalCandidatePanelError(f"Unreviewed decision for {identifier}")
        if not rationale:
            raise FinalCandidatePanelError(f"Missing review rationale for {identifier}")
        by_id[identifier] = row
    missing = set(expected) - set(by_id)
    if missing:
        raise FinalCandidatePanelError(f"Missing reviewed candidates: {sorted(missing)}")

    selected = [row for row in rows if str(row["review_decision"]).strip().lower() == "select"]
    reserves = [row for row in rows if str(row["review_decision"]).strip().lower() == "reserve"]
    if len(selected) != 30:
        raise FinalCandidatePanelError(f"Expected 30 selected candidates, found {len(selected)}")
    sequences: set[str] = set()
    final_rows = []
    for order, row in enumerate(sorted(selected, key=_selection_key), start=1):
        sequence = str(row["sequence"])
        mutation_set = str(row["mutation_set"])
        _validate_sequence(parent_sequence, sequence, mutation_set)
        if sequence in sequences:
            raise FinalCandidatePanelError("Selected sequences are not unique")
        sequences.add(sequence)
        final_rows.append(
            {
                **dict(row),
                "final_panel_order": order,
                "final_candidate_selection_performed": True,
                "experimental_validation_status": "recommended_for_testing_not_validated",
            }
        )
    reserve_rows = [
        {**dict(row), "final_candidate_selection_performed": True}
        for row in sorted(reserves, key=_selection_key)
    ]
    audit_rows = [
        {**dict(row), "final_candidate_selection_performed": True}
        for row in sorted(rows, key=lambda value: str(value["candidate_id"]))
    ]
    return {
        "final_rows": final_rows,
        "reserve_rows": reserve_rows,
        "audit_rows": audit_rows,
        "facts": {
            "reviewed_candidate_count": len(rows),
            "final_candidate_count": len(final_rows),
            "final_unique_sequence_count": len(sequences),
            "reserve_count": len(reserve_rows),
            "excluded_count": sum(
                str(row["review_decision"]).strip().lower() == "exclude" for row in rows
            ),
            "final_category_counts": dict(
                sorted(Counter(str(row["panel_category"]) for row in final_rows).items())
            ),
        },
    }


def _validate_sequence(parent: str, sequence: str, mutation_set: str) -> None:
    if len(parent) != 128 or len(sequence) != 128:
        raise FinalCandidatePanelError("Parent and candidates must be 128 aa")
    expected = list(parent)
    tokens = mutation_set.split(";")
    if len(tokens) not in {1, 2}:
        raise FinalCandidatePanelError(f"Invalid mutation set: {mutation_set}")
    positions = set()
    for token in tokens:
        match = MUTATION_RE.fullmatch(token)
        if not match:
            raise FinalCandidatePanelError(f"Invalid mutation token: {token}")
        wt, position_text, mutant = match.groups()
        position = int(position_text)
        if position in positions or position < 1 or position > len(parent):
            raise FinalCandidatePanelError(f"Invalid mutation position: {token}")
        positions.add(position)
        if parent[position - 1] != wt:
            raise FinalCandidatePanelError(f"WT mismatch: {token}")
        expected[position - 1] = mutant
    if sequence != "".join(expected):
        raise FinalCandidatePanelError(f"Sequence does not match {mutation_set}")
    if sequence[21] != "C" or sequence[94] != "C" or sequence[124:] != "SSGS":
        raise FinalCandidatePanelError(f"Immutable residues changed in {mutation_set}")
    if sequence.count("C") != parent.count("C"):
        raise FinalCandidatePanelError(f"Cysteine count changed in {mutation_set}")


def _selection_key(row: Mapping[str, object]) -> tuple[str, str]:
    return str(row["panel_category"]), str(row["candidate_id"])
