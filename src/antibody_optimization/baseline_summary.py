"""Join validated stage-1 artifacts into compact tabular summaries.

This module joins already validated sequence, expression, and optional
structure/interface artifacts.  It never infers a chain role, construct
boundary, or interface annotation.  Missing or unconfirmed structural inputs
remain explicit blocked gates instead of being replaced with placeholders.

Inputs are paths to Git-tracked CSV/JSON artifacts.  Returned values are plain
rows and dictionaries suitable for deterministic CSV/JSON serialization.  The
plot consumes only the compact rows returned by :func:`build_plot_rows` and
:func:`build_status_counts`.
"""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path
from typing import Iterable, Mapping, Sequence


NB252_UID = "LTT__Nb252"
PLOT_FIELDS = [
    "sample_uid",
    "sequence_index_1based",
    "residue",
    "imgt_position",
    "imgt_insertion_code",
    "imgt_position_label",
    "imgt_region",
    "numbering_status",
    "experimental_coordinate_status",
    "af3_coordinate_status",
    "collaborator_orange_annotation",
    "temporary_interface_lt4A",
]

COUNT_FIELDS = ["metric", "category", "count"]


class BaselineSummaryError(ValueError):
    """Raised when upstream artifacts violate the stage-1 data contract."""


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    """Read a human-facing UTF-8 CSV and require a non-empty unique header."""

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise BaselineSummaryError(f"Invalid CSV header: {path}")
        return list(reader)


def read_json_object(path: Path) -> dict[str, object]:
    """Read a JSON object; arrays and scalar top-level values are rejected."""

    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise BaselineSummaryError(f"JSON top level must be an object: {path}")
    return value


def build_plot_rows(
    *,
    expression_records: Sequence[Mapping[str, str]],
    numbering_review: Sequence[Mapping[str, str]],
    numbering_positions: Sequence[Mapping[str, str]],
    structure_mapping: Sequence[Mapping[str, str]] = (),
    interface_rows: Sequence[Mapping[str, str]] = (),
    structure_evidence_verified: bool = False,
    orange_annotation_verified: bool = False,
    interface_evidence_verified: bool = False,
) -> list[dict[str, object]]:
    """Build the exact 128-row Nb252 table consumed by the baseline figure.

    The reported ``sequence_raw`` is retained as the provisional reference.
    Structure and interface columns are populated only from supplied artifacts;
    absent evidence is represented by ``not_available`` or an empty flag.
    Several accepted aliases make this join tolerant of auth/label-oriented
    structure tables without changing the scientific semantics.
    """

    records = [row for row in expression_records if row.get("sample_uid") == NB252_UID]
    if len(records) != 1:
        raise BaselineSummaryError(f"Expected one {NB252_UID} expression row")
    sequence = records[0].get("sequence_raw", "")
    if len(sequence) != 128:
        raise BaselineSummaryError(f"Expected 128-aa Nb252 sequence, found {len(sequence)}")

    review_rows = [row for row in numbering_review if row.get("sample_uid") == NB252_UID]
    if len(review_rows) != 1:
        raise BaselineSummaryError(f"Expected one {NB252_UID} numbering review row")
    numbering_status = _first(review_rows[0], "numbering_status", "status") or "unknown"
    review = review_rows[0]
    if numbering_status != "pass":
        raise BaselineSummaryError("Nb252 must have a passed provisional numbering result")
    start_0based = _required_nonnegative_int(
        review, "query_start_0based_inclusive"
    )
    end_0based = _required_nonnegative_int(review, "query_end_0based_inclusive")
    if end_0based < start_0based or end_0based >= len(sequence):
        raise BaselineSummaryError("Invalid Nb252 provisional numbering query bounds")
    expected_span = sequence[start_0based : end_0based + 1]
    if _first(review, "provisional_numbered_span_sequence") != expected_span:
        raise BaselineSummaryError("Nb252 numbered span does not reconstruct sequence_raw")
    if _first(review, "unnumbered_n_sequence") != sequence[:start_0based]:
        raise BaselineSummaryError("Nb252 unnumbered N terminus does not reconstruct")
    if _first(review, "unnumbered_c_sequence") != sequence[end_0based + 1 :]:
        raise BaselineSummaryError("Nb252 unnumbered C terminus does not reconstruct")

    position_by_index: dict[int, Mapping[str, str]] = {}
    for row in numbering_positions:
        if row.get("sample_uid") != NB252_UID:
            continue
        index = _optional_index(
            row, "sequence_index_1based", "raw_sequence_index_1based"
        )
        if index is None:
            if _normal_bool(_first(row, "is_gap")) != "true":
                raise BaselineSummaryError(
                    "A numbering row without a sequence index is not marked as a gap"
                )
            continue
        if index in position_by_index:
            raise BaselineSummaryError(f"Duplicate Nb252 numbering index: {index}")
        position_by_index[index] = row
    expected_numbered_indices = set(range(start_0based + 1, end_0based + 2))
    if set(position_by_index) != expected_numbered_indices:
        missing = sorted(expected_numbered_indices - set(position_by_index))
        extra = sorted(set(position_by_index) - expected_numbered_indices)
        raise BaselineSummaryError(
            f"Nb252 numbering-index coverage mismatch: missing={missing}, extra={extra}"
        )

    structure_by_index: dict[int, dict[str, str]] = {}
    for row in structure_mapping if structure_evidence_verified else ():
        if not row.get("sample_uid"):
            raise BaselineSummaryError("Structure mapping requires explicit sample_uid")
        if row.get("sample_uid") != NB252_UID:
            continue
        index = _optional_index(row, "sequence_index_1based", "raw_sequence_index_1based")
        if index is None:
            continue
        if index > len(sequence):
            raise BaselineSummaryError(f"Structure mapping index exceeds Nb252: {index}")
        mapped_residue = _first(
            row, "sequence_residue", "raw_sequence_residue", "residue_aa"
        )
        if mapped_residue and mapped_residue != sequence[index - 1]:
            raise BaselineSummaryError(
                f"Structure mapping WT mismatch at {index}: {mapped_residue}"
            )
        role = _first(row, "source_model_role", "model_role", "evidence_class").lower()
        status = _first(row, "coordinate_status", "mapping_status") or "unmapped"
        target = structure_by_index.setdefault(index, {})
        if "experimental" in role or "nk2r_nb252" in role:
            if "experimental" in target:
                raise BaselineSummaryError(
                    f"Duplicate experimental structure mapping at {index}"
                )
            target["experimental"] = status
        elif "af3" in role or "prediction" in role:
            if "af3" in target:
                raise BaselineSummaryError(f"Duplicate AF3 structure mapping at {index}")
            target["af3"] = status
        else:
            raise BaselineSummaryError(f"Unrecognized structure model role: {role!r}")

    interface_by_index: dict[int, dict[str, str]] = {}
    for row in interface_rows:
        if not row.get("sample_uid"):
            raise BaselineSummaryError("Interface summary requires explicit sample_uid")
        if row.get("sample_uid") != NB252_UID:
            continue
        index = _optional_index(row, "sequence_index_1based", "raw_sequence_index_1based")
        if index is None:
            continue
        if index > len(sequence):
            raise BaselineSummaryError(f"Interface mapping index exceeds Nb252: {index}")
        if index in interface_by_index:
            raise BaselineSummaryError(f"Duplicate interface summary at {index}")
        mapped_residue = _first(row, "residue_aa", "sequence_residue")
        if mapped_residue and mapped_residue != sequence[index - 1]:
            raise BaselineSummaryError(
                f"Interface mapping WT mismatch at {index}: {mapped_residue}"
            )
        values = interface_by_index.setdefault(index, {})
        if orange_annotation_verified:
            values["orange"] = _normal_bool(
                _first(
                    row,
                    "confirmed_orange",
                    "collaborator_orange_annotation",
                    "session_orange",
                )
            )
        if interface_evidence_verified:
            values["interface"] = _normal_bool(
                _first(row, "temporary_interface_lt4A", "interface_lt_4A")
            )

    result: list[dict[str, object]] = []
    for index, residue in enumerate(sequence, start=1):
        numbered = position_by_index.get(index, {})
        numbered_residue = _first(numbered, "residue", "residue_aa", "amino_acid")
        if numbered_residue and numbered_residue != residue:
            raise BaselineSummaryError(
                f"Nb252 numbering residue mismatch at {index}: {numbered_residue} != {residue}"
            )
        imgt_position = _first(
            numbered, "imgt_position", "numbering_position", "position"
        )
        insertion = _first(numbered, "imgt_insertion_code", "insertion_code")
        region = _first(numbered, "imgt_region", "region") or "UNNUMBERED"
        position_label = _first(numbered, "numbering_position_label")
        structure = structure_by_index.get(index, {})
        interface = interface_by_index.get(index, {})
        result.append(
            {
                "sample_uid": NB252_UID,
                "sequence_index_1based": index,
                "residue": residue,
                "imgt_position": imgt_position,
                "imgt_insertion_code": insertion,
                "imgt_position_label": (
                    position_label
                    or (f"{imgt_position}{insertion}" if imgt_position else "")
                ),
                "imgt_region": region,
                "numbering_status": numbering_status,
                "experimental_coordinate_status": structure.get(
                    "experimental", "not_available"
                ),
                "af3_coordinate_status": structure.get("af3", "not_available"),
                "collaborator_orange_annotation": interface.get("orange", ""),
                "temporary_interface_lt4A": interface.get("interface", ""),
            }
        )
    return result


def build_status_counts(
    *,
    numbering_review: Sequence[Mapping[str, str]],
    sample_comparability: Sequence[Mapping[str, str]],
) -> list[dict[str, object]]:
    """Count numbering outcomes and maximum currently allowed data uses."""

    if len(numbering_review) != 47 or len(sample_comparability) != 47:
        raise BaselineSummaryError("Status summaries require exactly 47 rows per table")
    numbering = Counter(
        (_first(row, "numbering_status", "status") or "unknown")
        for row in numbering_review
    )
    allowed = Counter(
        (row.get("highest_allowed_use") or "unknown") for row in sample_comparability
    )
    rows = [
        {"metric": "numbering_status", "category": key, "count": numbering[key]}
        for key in sorted(numbering)
    ]
    rows.extend(
        {"metric": "highest_allowed_use", "category": key, "count": allowed[key]}
        for key in sorted(allowed)
    )
    return rows


def _first(row: Mapping[str, object], *keys: str) -> str:
    for key in keys:
        value = row.get(key)
        if value is not None and str(value) != "":
            return str(value)
    return ""


def _parse_index(row: Mapping[str, object], *keys: str) -> int:
    value = _first(row, *keys)
    try:
        index = int(value)
    except ValueError as exc:
        raise BaselineSummaryError(f"Invalid sequence index: {value!r}") from exc
    if index < 1:
        raise BaselineSummaryError(f"Sequence index must be positive: {index}")
    return index


def _optional_index(row: Mapping[str, object], *keys: str) -> int | None:
    value = _first(row, *keys)
    if not value:
        return None
    return _parse_index(row, *keys)


def _required_nonnegative_int(row: Mapping[str, object], key: str) -> int:
    value = _first(row, key)
    try:
        parsed = int(value)
    except ValueError as exc:
        raise BaselineSummaryError(f"Invalid {key}: {value!r}") from exc
    if parsed < 0:
        raise BaselineSummaryError(f"{key} must be nonnegative")
    return parsed


def _normal_bool(value: str) -> str:
    if not value:
        return ""
    lowered = value.strip().lower()
    if lowered in {"true", "1", "yes"}:
        return "true"
    if lowered in {"false", "0", "no"}:
        return "false"
    raise BaselineSummaryError(f"Invalid boolean value: {value!r}")
