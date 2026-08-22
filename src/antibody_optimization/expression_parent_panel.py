"""Build the approved 19-single-mutant parent set for double-mutant design.

The module only subsets an already released 30-member computational trial
panel.  It does not rescore candidates, alter magnitude bands, or enumerate
double mutants.  Selection decisions are supplied explicitly by the caller so
that mentor/user judgement remains visible and auditable.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping, Sequence
from math import comb


class ParentPanelError(ValueError):
    """Raised when the source panel or explicit selection contract is invalid."""


def mutation_code(row: Mapping[str, object]) -> str:
    """Return the compact mutation code from a trial-panel record."""
    label = str(row.get("mutation_reported_label", ""))
    prefix = "Nb252 reported_seq "
    if not label.startswith(prefix):
        raise ParentPanelError(f"Unexpected mutation label: {label!r}")
    code = label[len(prefix) :]
    if not code:
        raise ParentPanelError("Empty compact mutation code")
    return code


def build_parent19_panel(
    trial_rows: Sequence[Mapping[str, object]],
    selected_mutations: Sequence[str],
    selected_reasons: Mapping[str, str],
    *,
    focal_positions: Sequence[int] = (30, 1, 27),
) -> dict[str, object]:
    """Subset a 30-row trial panel under the approved 3/3/3 + 10x1 rule.

    Args:
        trial_rows: Existing 30-member single-mutant trial panel.  Rows are
            treated as immutable evidence and copied into returned records.
        selected_mutations: Ordered compact mutation codes for the new panel.
        selected_reasons: One decision-relevant rationale per selected code.
        focal_positions: Positions that must contribute exactly three retained
            substitutions each.  Every other source position contributes one.

    Returns:
        A dictionary containing the ordered 19-row panel, a 30-row audit, and
        exact combinatorial counts for the next double-enumeration stage.

    Out of scope:
        Predictor reruns, mutation generation outside the supplied trial panel,
        double-mutant enumeration, and final experimental release.
    """
    rows = [dict(row) for row in trial_rows]
    if len(rows) != 30:
        raise ParentPanelError(f"Expected 30 source rows, observed {len(rows)}")

    by_code: dict[str, dict[str, object]] = {}
    candidate_ids: set[str] = set()
    sequences: set[str] = set()
    source_position_counts: Counter[int] = Counter()
    for row in rows:
        code = mutation_code(row)
        candidate_id = str(row.get("candidate_id", ""))
        sequence = str(row.get("sequence", ""))
        try:
            position = int(row["reported_sequence_index_1based"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ParentPanelError(f"Invalid reported position for {code}") from exc
        if code in by_code or not candidate_id or candidate_id in candidate_ids:
            raise ParentPanelError(f"Duplicate or empty candidate identity for {code}")
        if len(sequence) != 128 or sequence in sequences:
            raise ParentPanelError(f"Invalid or duplicate 128-aa sequence for {code}")
        if str(row.get("trial_selection_status")) != "trial_final30":
            raise ParentPanelError(f"Source row is not in the frozen trial30 panel: {code}")
        by_code[code] = row
        candidate_ids.add(candidate_id)
        sequences.add(sequence)
        source_position_counts[position] += 1

    selected = list(selected_mutations)
    if len(selected) != 19 or len(set(selected)) != 19:
        raise ParentPanelError("The parent panel must contain 19 unique mutation codes")
    missing = [code for code in selected if code not in by_code]
    if missing:
        raise ParentPanelError(f"Selected mutations are absent from trial30: {missing}")
    if set(selected_reasons) != set(selected):
        raise ParentPanelError("Selected reasons must cover exactly the 19 selected mutations")
    if any(not str(selected_reasons[code]).strip() for code in selected):
        raise ParentPanelError("Every selected mutation requires a non-empty rationale")

    focal = tuple(int(position) for position in focal_positions)
    if len(focal) != 3 or len(set(focal)) != 3:
        raise ParentPanelError("Exactly three unique focal positions are required")
    selected_position_counts = Counter(
        int(by_code[code]["reported_sequence_index_1based"]) for code in selected
    )
    if any(selected_position_counts[position] != 3 for position in focal):
        raise ParentPanelError("Each focal position must retain exactly three substitutions")
    nonfocal_positions = set(source_position_counts) - set(focal)
    if len(nonfocal_positions) != 10:
        raise ParentPanelError(
            f"Expected 10 nonfocal source positions, observed {len(nonfocal_positions)}"
        )
    if any(selected_position_counts[position] != 1 for position in nonfocal_positions):
        raise ParentPanelError("Each nonfocal source position must retain exactly one substitution")
    if set(selected_position_counts) != set(source_position_counts):
        raise ParentPanelError("The selected panel must preserve all 13 source positions")

    panel_rows: list[dict[str, object]] = []
    selected_order = {code: index for index, code in enumerate(selected, start=1)}
    for code in selected:
        row = dict(by_code[code])
        row.update(
            {
                "parent19_selection_status": "selected_parent19",
                "parent19_selection_order": selected_order[code],
                "parent19_selection_reason": selected_reasons[code],
                "approved_as_double_mutant_parent": True,
            }
        )
        panel_rows.append(row)

    audit_rows: list[dict[str, object]] = []
    for row in sorted(rows, key=lambda value: int(value["trial_selection_order"])):
        code = mutation_code(row)
        position = int(row["reported_sequence_index_1based"])
        selected_status = code in selected_order
        if selected_status:
            reason = selected_reasons[code]
        elif position in focal:
            reason = "not_retained_within_three_substitution_focal_position_quota"
        else:
            reason = "not_retained_within_one_substitution_nonfocal_position_quota"
        audit = dict(row)
        audit.update(
            {
                "parent19_selection_status": (
                    "selected_parent19" if selected_status else "not_selected_parent19"
                ),
                "parent19_selection_order": selected_order.get(code, ""),
                "parent19_selection_reason": reason,
                "approved_as_double_mutant_parent": selected_status,
            }
        )
        audit_rows.append(audit)

    all_pairs = comb(len(selected), 2)
    invalid_same_position_pairs = sum(comb(count, 2) for count in selected_position_counts.values())
    facts = {
        "source_trial_single_mutant_count": len(rows),
        "source_position_count": len(source_position_counts),
        "selected_parent_single_mutant_count": len(panel_rows),
        "selected_position_count": len(selected_position_counts),
        "focal_position_retained_counts": {
            str(position): selected_position_counts[position] for position in focal
        },
        "nonfocal_position_count": len(nonfocal_positions),
        "nonfocal_retained_per_position": 1,
        "theoretical_all_pair_count": all_pairs,
        "invalid_same_position_pair_count": invalid_same_position_pairs,
        "planned_valid_double_mutant_count": all_pairs - invalid_same_position_pairs,
        "double_mutant_enumeration_performed": False,
        "existing_trial30_artifacts_modified": False,
    }
    if facts["planned_valid_double_mutant_count"] != 162:
        raise ParentPanelError("Expected 162 valid future double-mutant combinations")
    return {"panel_rows": panel_rows, "audit_rows": audit_rows, "facts": facts}
