"""Build the ordered, non-ranking V3 parent-single expert-review pool.

This module intentionally uses only the Python standard library so both the
project environment and ChimeraX's bundled Python can import it.  It preserves
the immutable 30-row V3 shortlist and appends explicitly requested candidates
from a complete upstream audit.  Review-pool order is display/provenance order,
not a scientific rank or a parent-selection result.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence


V3_BASE_SHORTLIST_COUNT = 30
V3_SUPPLEMENTAL_T99F_ID = "Nb252_expr_seq099_T99F"
V3_REVIEW_POOL_COUNT = 31


class V3ExpertReviewPoolError(ValueError):
    """Raised when shortlist and supplemental-candidate identities disagree."""


def build_v3_review_candidate_pool(
    shortlist_rows: Sequence[Mapping[str, object]],
    supplemental_source_rows: Sequence[Mapping[str, object]],
    *,
    supplemental_candidate_ids: Sequence[str] = (V3_SUPPLEMENTAL_T99F_ID,),
) -> list[dict[str, object]]:
    """Return the immutable shortlist plus requested supplemental candidates."""

    shortlist = [dict(row) for row in shortlist_rows]
    if len(shortlist) != V3_BASE_SHORTLIST_COUNT:
        raise V3ExpertReviewPoolError(
            f"Expected {V3_BASE_SHORTLIST_COUNT} immutable V3 shortlist rows"
        )
    shortlist_ids = [str(row.get("candidate_id", "")).strip() for row in shortlist]
    if len(set(shortlist_ids)) != V3_BASE_SHORTLIST_COUNT or any(
        not value for value in shortlist_ids
    ):
        raise V3ExpertReviewPoolError(
            "V3 shortlist candidate identities are not unique"
        )
    try:
        shortlist_orders = [int(row["selection_order_v3"]) for row in shortlist]
    except (KeyError, TypeError, ValueError) as error:
        raise V3ExpertReviewPoolError(
            "V3 shortlist lacks valid selection_order_v3 values"
        ) from error
    if sorted(shortlist_orders) != list(range(1, V3_BASE_SHORTLIST_COUNT + 1)):
        raise V3ExpertReviewPoolError(
            "V3 shortlist selection_order_v3 must contain integers 1 through 30"
        )

    supplemental_lookup: dict[str, dict[str, object]] = {}
    for row in supplemental_source_rows:
        identifier = str(row.get("candidate_id", "")).strip()
        if not identifier:
            continue
        if identifier in supplemental_lookup:
            raise V3ExpertReviewPoolError(
                f"Supplemental source contains duplicate candidate: {identifier}"
            )
        supplemental_lookup[identifier] = dict(row)

    requested_ids = [str(value).strip() for value in supplemental_candidate_ids]
    if not requested_ids or any(not value for value in requested_ids):
        raise V3ExpertReviewPoolError(
            "At least one supplemental candidate is required"
        )
    if len(set(requested_ids)) != len(requested_ids):
        raise V3ExpertReviewPoolError(
            "Supplemental candidate identities are not unique"
        )
    overlap = sorted(set(requested_ids) & set(shortlist_ids))
    if overlap:
        raise V3ExpertReviewPoolError(
            f"Supplemental candidates already occur in the V3 shortlist: {overlap}"
        )
    missing = [
        identifier
        for identifier in requested_ids
        if identifier not in supplemental_lookup
    ]
    if missing:
        raise V3ExpertReviewPoolError(
            f"Supplemental source lacks requested candidates: {missing}"
        )

    pool: list[dict[str, object]] = []
    for row in sorted(shortlist, key=lambda value: int(value["selection_order_v3"])):
        row["review_pool_order"] = int(row["selection_order_v3"])
        row["review_pool_role"] = "immutable_v3_upstream_shortlist"
        pool.append(row)
    for offset, identifier in enumerate(requested_ids, start=1):
        row = dict(supplemental_lookup[identifier])
        row["review_pool_order"] = V3_BASE_SHORTLIST_COUNT + offset
        row["review_pool_role"] = "user_added_stable_word_exploratory"
        pool.append(row)

    if len({str(row["candidate_id"]) for row in pool}) != len(pool):
        raise V3ExpertReviewPoolError(
            "Combined expert-review pool identities are not unique"
        )
    return pool


__all__ = [
    "V3_BASE_SHORTLIST_COUNT",
    "V3_REVIEW_POOL_COUNT",
    "V3_SUPPLEMENTAL_T99F_ID",
    "V3ExpertReviewPoolError",
    "build_v3_review_candidate_pool",
]
