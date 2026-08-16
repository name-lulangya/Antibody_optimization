"""Audit paired-WT contact changes in Nb252 double-mutant PyRosetta results.

This module only compares contact sets already written by the released scoring
run. It does not read or modify structures, recalculate geometric contacts, or
apply candidate-selection rules.
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from typing import Mapping, Sequence


class DoubleMutantContactError(ValueError):
    """Raised when paired contact records violate the released run contract."""


def audit_paired_contacts(
    paired_rows: Sequence[Mapping[str, object]],
    wt_controls: Sequence[Mapping[str, object]],
    *,
    expected_candidate_count: int = 86,
    expected_replicates: int = 3,
    expected_wt_control_count: int = 135,
) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, object]]:
    """Recompute per-repeat lost/gained sets and candidate-level summaries.

    Inputs are the unfiltered candidate-repeat and position-pair WT tables from
    ``score_double_mutants_pyrosetta.py``. Retention values in the paired table
    are independently recomputed and required to agree within floating-point
    tolerance. Source-auth residue identifiers remain integer strings.
    """

    if len(paired_rows) != expected_candidate_count * expected_replicates:
        raise DoubleMutantContactError("Unexpected paired-row coverage")
    if len(wt_controls) != expected_wt_control_count:
        raise DoubleMutantContactError("Unexpected WT-control coverage")
    wt_map = _unique(wt_controls, "wt_control_id", "WT controls")
    candidate_counts = Counter(str(row["candidate_id"]) for row in paired_rows)
    if len(candidate_counts) != expected_candidate_count or set(
        candidate_counts.values()
    ) != {expected_replicates}:
        raise DoubleMutantContactError("Candidate replicate coverage is incomplete")

    replicate_rows: list[dict[str, object]] = []
    seen_candidate_replicates: set[tuple[str, int]] = set()
    referenced_wt_ids: set[str] = set()
    for source in paired_rows:
        candidate_id = str(source["candidate_id"])
        replicate = int(source["replicate"])
        key = (candidate_id, replicate)
        if key in seen_candidate_replicates:
            raise DoubleMutantContactError(f"Duplicate candidate replicate: {key}")
        seen_candidate_replicates.add(key)
        wt_id = str(source["wt_control_id"])
        if wt_id not in wt_map:
            raise DoubleMutantContactError(f"Unknown WT control: {wt_id}")
        referenced_wt_ids.add(wt_id)
        wt = wt_map[wt_id]
        if int(wt["replicate"]) != replicate or int(wt["seed"]) != int(
            source["seed"]
        ):
            raise DoubleMutantContactError(
                f"WT replicate/seed mismatch for {candidate_id} replicate {replicate}"
            )
        if str(source["status"]) != "pass" or str(wt["status"]) != "pass":
            raise DoubleMutantContactError("Contact audit requires pass records")

        wt_vhh = _position_set(wt["vhh_contact_auth_positions"])
        wt_receptor = _position_set(wt["receptor_contact_auth_positions"])
        mutant_vhh = _position_set(source["mutant_vhh_contact_auth_positions"])
        mutant_receptor = _position_set(
            source["mutant_receptor_contact_auth_positions"]
        )
        if len(wt_vhh) != int(source["paired_wt_vhh_contact_count"]) or len(
            wt_receptor
        ) != int(source["paired_wt_receptor_epitope_count"]):
            raise DoubleMutantContactError("Paired WT contact count is inconsistent")
        vhh_retention = _retention(wt_vhh, mutant_vhh)
        receptor_retention = _retention(wt_receptor, mutant_receptor)
        if not math.isclose(
            vhh_retention,
            float(source["candidate_vs_paired_wt_vhh_contact_retention"]),
            rel_tol=0,
            abs_tol=1e-12,
        ) or not math.isclose(
            receptor_retention,
            float(source["candidate_vs_paired_wt_receptor_epitope_retention"]),
            rel_tol=0,
            abs_tol=1e-12,
        ):
            raise DoubleMutantContactError("Stored paired retention is inconsistent")

        vhh_lost = wt_vhh - mutant_vhh
        vhh_gained = mutant_vhh - wt_vhh
        receptor_lost = wt_receptor - mutant_receptor
        receptor_gained = mutant_receptor - wt_receptor
        replicate_rows.append(
            {
                "candidate_id": candidate_id,
                "mutation_reported_label": source["mutation_reported_label"],
                "position_pair": source["position_pair"],
                "replicate": replicate,
                "seed": int(source["seed"]),
                "wt_control_id": wt_id,
                "paired_wt_vhh_contact_auth_positions": _position_text(wt_vhh),
                "mutant_vhh_contact_auth_positions": _position_text(mutant_vhh),
                "vhh_lost_auth_positions": _position_text(vhh_lost),
                "vhh_gained_auth_positions": _position_text(vhh_gained),
                "paired_wt_vhh_contact_retention": vhh_retention,
                "paired_wt_receptor_contact_auth_positions": _position_text(
                    wt_receptor
                ),
                "mutant_receptor_contact_auth_positions": _position_text(
                    mutant_receptor
                ),
                "receptor_lost_auth_positions": _position_text(receptor_lost),
                "receptor_gained_auth_positions": _position_text(receptor_gained),
                "paired_wt_receptor_epitope_retention": receptor_retention,
                "paired_contact_change": bool(
                    vhh_lost or vhh_gained or receptor_lost or receptor_gained
                ),
            }
        )
    if referenced_wt_ids != set(wt_map):
        raise DoubleMutantContactError("WT controls are not referenced exactly in scope")

    candidate_groups: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in replicate_rows:
        candidate_groups[str(row["candidate_id"])].append(row)
    summaries = [
        _summarize_candidate(candidate_id, selected, expected_replicates)
        for candidate_id, selected in candidate_groups.items()
    ]
    summaries.sort(key=lambda row: str(row["candidate_id"]))
    replicate_rows.sort(
        key=lambda row: (str(row["candidate_id"]), int(row["replicate"]))
    )
    changed = [row for row in summaries if row["paired_contact_change_status"] == "changed"]
    facts = {
        "paired_row_count": len(replicate_rows),
        "wt_control_count": len(wt_controls),
        "candidate_count": len(summaries),
        "paired_contact_changed_candidate_count": len(changed),
        "paired_contact_unchanged_candidate_count": len(summaries) - len(changed),
        "vhh_changed_candidate_count": sum(
            bool(row["paired_wt_vhh_contact_change"]) for row in summaries
        ),
        "receptor_changed_candidate_count": sum(
            bool(row["paired_wt_receptor_contact_change"]) for row in summaries
        ),
    }
    return replicate_rows, summaries, facts


def _summarize_candidate(
    candidate_id: str,
    rows: Sequence[Mapping[str, object]],
    expected_replicates: int,
) -> dict[str, object]:
    if len(rows) != expected_replicates:
        raise DoubleMutantContactError(
            f"Unexpected replicate count for {candidate_id}"
        )
    if sorted(int(row["replicate"]) for row in rows) != list(
        range(1, expected_replicates + 1)
    ):
        raise DoubleMutantContactError(
            f"Unexpected replicate identities for {candidate_id}"
        )
    if len({str(row["mutation_reported_label"]) for row in rows}) != 1 or len(
        {str(row["position_pair"]) for row in rows}
    ) != 1:
        raise DoubleMutantContactError(
            f"Mutation or position-pair identity varies for {candidate_id}"
        )
    first = rows[0]
    vhh_lost = _union_field(rows, "vhh_lost_auth_positions")
    vhh_gained = _union_field(rows, "vhh_gained_auth_positions")
    receptor_lost = _union_field(rows, "receptor_lost_auth_positions")
    receptor_gained = _union_field(rows, "receptor_gained_auth_positions")
    vhh_changed = bool(vhh_lost or vhh_gained)
    receptor_changed = bool(receptor_lost or receptor_gained)
    patterns = {
        (
            str(row["vhh_lost_auth_positions"]),
            str(row["vhh_gained_auth_positions"]),
            str(row["receptor_lost_auth_positions"]),
            str(row["receptor_gained_auth_positions"]),
        )
        for row in rows
    }
    return {
        "candidate_id": candidate_id,
        "mutation_reported_label": first["mutation_reported_label"],
        "position_pair": first["position_pair"],
        "replicate_count": len(rows),
        "minimum_recomputed_paired_wt_vhh_contact_retention": min(
            float(row["paired_wt_vhh_contact_retention"]) for row in rows
        ),
        "minimum_recomputed_paired_wt_receptor_epitope_retention": min(
            float(row["paired_wt_receptor_epitope_retention"]) for row in rows
        ),
        "paired_contact_changed_replicate_count": sum(
            bool(row["paired_contact_change"]) for row in rows
        ),
        "paired_wt_vhh_contact_change": vhh_changed,
        "paired_wt_vhh_lost_auth_positions_union": _position_text(vhh_lost),
        "paired_wt_vhh_gained_auth_positions_union": _position_text(vhh_gained),
        "paired_wt_receptor_contact_change": receptor_changed,
        "paired_wt_receptor_lost_auth_positions_union": _position_text(
            receptor_lost
        ),
        "paired_wt_receptor_gained_auth_positions_union": _position_text(
            receptor_gained
        ),
        "paired_contact_change_pattern_count": len(patterns),
        "paired_contact_change_replicate_concordance": (
            "all_replicates_same" if len(patterns) == 1 else "replicate_variable"
        ),
        "paired_contact_change_status": (
            "changed" if vhh_changed or receptor_changed else "unchanged"
        ),
    }


def _union_field(
    rows: Sequence[Mapping[str, object]], field: str
) -> set[int]:
    result: set[int] = set()
    for row in rows:
        result.update(_position_set(row[field]))
    return result


def _position_set(value: object) -> set[int]:
    return {int(item) for item in str(value).split(";") if item}


def _position_text(values: set[int]) -> str:
    return ";".join(str(value) for value in sorted(values))


def _retention(reference: set[int], candidate: set[int]) -> float:
    if not reference:
        raise DoubleMutantContactError("Paired WT contact set must not be empty")
    return len(reference & candidate) / len(reference)


def _unique(
    rows: Sequence[Mapping[str, object]], key: str, label: str
) -> dict[str, Mapping[str, object]]:
    result = {str(row[key]): row for row in rows}
    if len(result) != len(rows):
        raise DoubleMutantContactError(f"Duplicate identities in {label}")
    return result
