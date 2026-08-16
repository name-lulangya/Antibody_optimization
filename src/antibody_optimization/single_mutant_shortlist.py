"""Magnitude- and risk-aware narrowing of the Nb252 V2 single-mutant pool.

The V2 review remains the immutable evidence layer.  This module adds a
decision layer that removes property-only test candidates when strong negative
or decision-changing structural evidence outweighs their exploratory property
signal.  It does not rescore candidates, alter affinity candidates, or generate
multi-mutants.
"""

from __future__ import annotations

from collections import Counter
from typing import Mapping, Sequence


class SingleMutantShortlistError(ValueError):
    """Raised when the released V2 review violates the shortlist contract."""


EXPECTED_ROWS = 80
EXPECTED_PROPERTY_ACTIVE_BEFORE = 22
EXPECTED_AFFINITY_ACTIVE_BEFORE = 8
ACTIVE_V2_STATUSES = {
    "combination_ready",
    "single_mutant_test_only",
    "targeted_alternative_review",
}


def build_single_mutant_shortlist(
    rows: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """Return an 80-row audit table and the narrowed active single-mutant pool.

    Property candidates already released as ``combination_ready`` are retained.
    A property ``single_mutant_test_only`` candidate is deprioritized when it
    has a receptor-contact change, failed the AF3 local non-adverse gate while a
    same-position alternative passed, directionally adverse paired affinity,
    or a pre-existing strong-negative AntiFold/exposed-hydrophobe flag.  The
    rules consume existing V2 evidence only; no thresholds are recomputed.
    """

    if len(rows) != EXPECTED_ROWS:
        raise SingleMutantShortlistError(f"Expected {EXPECTED_ROWS} V2 rows")
    if len({str(row["candidate_id"]) for row in rows}) != EXPECTED_ROWS:
        raise SingleMutantShortlistError("V2 candidate identifiers are not unique")

    active_before = [row for row in rows if str(row["v2_qualification_status"]) in ACTIVE_V2_STATUSES]
    before_by_track = Counter(str(row["design_track"]) for row in active_before)
    if before_by_track != {
        "property": EXPECTED_PROPERTY_ACTIVE_BEFORE,
        "affinity": EXPECTED_AFFINITY_ACTIVE_BEFORE,
    }:
        raise SingleMutantShortlistError(f"Unexpected V2 active pool: {dict(before_by_track)}")

    review_rows: list[dict[str, object]] = []
    for source in rows:
        row = dict(source)
        track = str(row["design_track"])
        status = str(row["v2_qualification_status"])
        triggers: list[str] = []
        primary = ""

        if track == "property" and status == "single_mutant_test_only":
            flags = _tokens(row.get("structural_review_flags"))
            affinity_class = str(row.get("property_affinity_direction_class", ""))
            v2_reason = str(row.get("v2_qualification_reason", ""))
            if "paired_receptor_contact_change" in flags:
                triggers.append("paired_receptor_contact_change")
                primary = primary or "paired_receptor_contact_change"
            if v2_reason == "af3_complete_vhh_local_nonadverse_gate_not_met":
                triggers.append("af3_local_nonadverse_gate_not_met")
                primary = primary or "af3_local_nonadverse_gate_not_met"
            if affinity_class == "directionally_adverse":
                triggers.append("paired_affinity_directionally_adverse")
                primary = primary or "paired_affinity_directionally_adverse"
            if {
                "strong_negative_antifold_complex_signal",
                "exposed_hydrophobic_substitution",
            }.issubset(flags):
                triggers.append("strong_negative_antifold_and_exposed_hydrophobe")
                primary = primary or "strong_negative_antifold_and_exposed_hydrophobe"
            elif "strong_negative_antifold_complex_signal" in flags:
                triggers.append("strong_negative_antifold_complex_signal")
                primary = primary or "strong_negative_antifold_complex_signal"

        if triggers:
            decision = "deprioritize"
            role = "none"
            reason = primary
        elif status == "combination_ready":
            decision = "retain_active"
            role = "combination_module_review"
            reason = "v2_combination_ready_retained"
        elif status == "single_mutant_test_only":
            decision = "retain_active"
            role = "single_mutant_test"
            reason = "no_decision_changing_integrated_risk_trigger"
        elif status == "targeted_alternative_review":
            decision = "retain_active"
            role = "affinity_alternative_review"
            reason = "v2_targeted_affinity_alternative_retained"
        elif status == "do_not_advance":
            decision = "do_not_advance"
            role = "none"
            reason = str(row["v2_qualification_reason"])
        else:
            decision = "retain_audit_only"
            role = "none"
            reason = "v2_not_prioritized_retained_for_audit"

        row["shortlist_decision"] = decision
        row["shortlist_role"] = role
        row["shortlist_reason"] = reason
        row["shortlist_trigger_count"] = len(triggers)
        row["shortlist_triggers"] = ";".join(triggers)
        review_rows.append(row)

    shortlist_rows = [row for row in review_rows if row["shortlist_decision"] == "retain_active"]
    after_by_track = Counter(str(row["design_track"]) for row in shortlist_rows)
    property_deprioritized = [
        row for row in review_rows
        if row["design_track"] == "property" and row["shortlist_decision"] == "deprioritize"
    ]
    facts = {
        "candidate_count": len(review_rows),
        "active_before_count": len(active_before),
        "active_after_count": len(shortlist_rows),
        "active_before_by_track": dict(before_by_track),
        "active_after_by_track": dict(after_by_track),
        "property_deprioritized_count": len(property_deprioritized),
        "decision_counts": dict(Counter(str(row["shortlist_decision"]) for row in review_rows)),
        "role_counts": dict(Counter(str(row["shortlist_role"]) for row in shortlist_rows)),
        "primary_deprioritization_reason_counts": dict(
            Counter(str(row["shortlist_reason"]) for row in property_deprioritized)
        ),
        "retained_mutations": [str(row["mutation"]) for row in shortlist_rows],
        "deprioritized_property_mutations": [str(row["mutation"]) for row in property_deprioritized],
        "combination_generated": False,
    }
    if len(shortlist_rows) != 14 or after_by_track != {"affinity": 8, "property": 6}:
        raise SingleMutantShortlistError(f"Unexpected narrowed pool: {facts}")
    return {"review_rows": review_rows, "shortlist_rows": shortlist_rows, "facts": facts}


def _tokens(value: object) -> set[str]:
    return {token for token in str(value or "").split(";") if token}
