"""Validated provisional IMGT numbering for the expression-sequence baseline.

This module keeps three concepts separate: the collaborator-transcribed input
sequence, the span that ANARCII can number, and the unresolved biological VHH
construct boundary.  A successful ANARCII result therefore produces a
``provisional_numbered_domain`` span; it never populates or replaces the source
``vhh_region_sequence`` field.

The pure validation and transformation functions accept ordinary Python data
and do not read files, invoke ANARCII, or mutate their inputs.  File loading,
tool execution, and artifact writing are thin explicit wrappers around those
functions.
"""

from __future__ import annotations

import hashlib
import math
from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass


ANARCII_VERSION = "2.0.8"
ANARCII_PARAMETERS: dict[str, object] = {
    "seq_type": "antibody",
    "mode": "accuracy",
    "scheme": "imgt",
    "cpu": True,
    "ncpu": 1,
    "batch_size": 8,
    "scfv": False,
}
EXPECTED_RECORD_COUNT = 47
STANDARD_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")

REQUIRED_INPUT_FIELDS = (
    "sample_uid",
    "provider_code",
    "source_sample_id",
    "sequence_raw",
    "sequence_length_aa",
    "sequence_sha256",
    "sequence_scope",
    "vhh_region_sequence",
    "sequence_review_flags",
)

class SequenceNumberingError(ValueError):
    """Raised when input provenance or an ANARCII result is inconsistent."""


@dataclass(frozen=True)
class InputSequence:
    """One immutable sequence record copied from the validated input table."""

    sample_uid: str
    provider_code: str
    source_sample_id: str
    sequence_raw: str
    sequence_length_aa: int
    sequence_sha256: str
    source_sequence_scope: str
    source_vhh_region_sequence: str
    source_sequence_review_flags: str


@dataclass(frozen=True)
class NumberingPosition:
    """One ANARCII IMGT output position and its source-sequence mapping."""

    position_order: int
    numbering_position: int
    insertion_code: str
    residue_aa: str
    sequence_index_0based: int | None
    region: str

    @property
    def label(self) -> str:
        return f"{self.numbering_position}{self.insertion_code}"

    @property
    def is_gap(self) -> bool:
        return self.residue_aa == "-"


@dataclass(frozen=True)
class NumberingAudit:
    """Validated sequence-level interpretation of one raw ANARCII result."""

    input_record: InputSequence
    numbering_status: str
    sequence_scope_status: str
    chain_type: str
    score: float
    scheme: str
    query_start_0based_inclusive: int | None
    query_end_0based_inclusive: int | None
    provisional_numbered_span_sequence: str
    provisional_numbered_span_sha256: str
    unnumbered_n_sequence: str
    unnumbered_c_sequence: str
    first_numbered_imgt_position: str
    last_numbered_imgt_position: str
    numbering_review_flags: tuple[str, ...]
    error: str
    positions: tuple[NumberingPosition, ...]


def sha256_bytes(value: bytes) -> str:
    """Return the lowercase SHA-256 hex digest of ``value``."""

    return hashlib.sha256(value).hexdigest()


def validate_sequence_rows(
    fieldnames: Sequence[str],
    rows: Sequence[Mapping[str, object]],
    *,
    expected_count: int = EXPECTED_RECORD_COUNT,
) -> tuple[InputSequence, ...]:
    """Validate input identities, literal lengths, alphabets, and sequence hashes.

    Args:
        fieldnames: Observed CSV header.
        rows: Parsed row mappings, which are never modified.
        expected_count: Exact number of required sequence rows.

    Returns:
        An immutable tuple of normalized ``InputSequence`` records in source
        row order.  Sequence strings themselves are not normalized.

    Scope:
        This pure function checks transcription integrity only.  It deliberately
        does not decide whether a sequence is a VHH, trim terminal residues, or
        populate the unresolved ``vhh_region_sequence`` field.
    """

    missing_fields = [field for field in REQUIRED_INPUT_FIELDS if field not in fieldnames]
    if missing_fields:
        raise SequenceNumberingError(
            f"Input records table lacks required fields: {missing_fields}"
        )
    if len(rows) != expected_count:
        raise SequenceNumberingError(
            f"Input records table must contain {expected_count} rows, found {len(rows)}"
        )

    records: list[InputSequence] = []
    seen_uids: set[str] = set()
    for row_number, row in enumerate(rows, start=2):
        sample_uid = _row_text(row, "sample_uid")
        if not sample_uid:
            raise SequenceNumberingError(f"Blank sample_uid at CSV row {row_number}")
        if sample_uid in seen_uids:
            raise SequenceNumberingError(f"Duplicate sample_uid: {sample_uid}")
        seen_uids.add(sample_uid)

        sequence = _row_text(row, "sequence_raw")
        if not sequence or sequence != sequence.strip():
            raise SequenceNumberingError(
                f"Sequence must be non-empty with no outer whitespace: {sample_uid}"
            )
        unsupported = sorted(set(sequence) - STANDARD_AMINO_ACIDS)
        if unsupported:
            raise SequenceNumberingError(
                f"Sequence contains unsupported residues for {sample_uid}: {unsupported}"
            )

        length_text = _row_text(row, "sequence_length_aa")
        try:
            recorded_length = int(length_text)
        except ValueError as exc:
            raise SequenceNumberingError(
                f"Invalid sequence_length_aa for {sample_uid}: {length_text!r}"
            ) from exc
        if recorded_length != len(sequence):
            raise SequenceNumberingError(
                f"Sequence length mismatch for {sample_uid}: "
                f"recorded {recorded_length}, observed {len(sequence)}"
            )

        recorded_sha256 = _row_text(row, "sequence_sha256")
        observed_sha256 = sha256_bytes(sequence.encode("ascii"))
        if recorded_sha256 != observed_sha256:
            raise SequenceNumberingError(
                f"Sequence SHA-256 mismatch for {sample_uid}"
            )

        records.append(
            InputSequence(
                sample_uid=sample_uid,
                provider_code=_row_text(row, "provider_code"),
                source_sample_id=_row_text(row, "source_sample_id"),
                sequence_raw=sequence,
                sequence_length_aa=recorded_length,
                sequence_sha256=recorded_sha256,
                source_sequence_scope=_row_text(row, "sequence_scope"),
                source_vhh_region_sequence=_row_text(row, "vhh_region_sequence"),
                source_sequence_review_flags=_row_text(row, "sequence_review_flags"),
            )
        )
    return tuple(records)


def imgt_region(numbering_position: int) -> str:
    """Map an IMGT V-domain numeric position to its provisional region.

    The fixed ranges are FR1 1-26, CDR1 27-38, FR2 39-55, CDR2 56-65,
    FR3 66-104, CDR3 105-117, and FR4 118-128.  An insertion inherits the
    region of its numeric base position.  This pure mapping is a reporting aid;
    it does not establish a biological VHH construct boundary.
    """

    if not 1 <= numbering_position <= 128:
        raise SequenceNumberingError(
            f"IMGT V-domain position must be within 1..128: {numbering_position}"
        )
    if numbering_position <= 26:
        return "FR1"
    if numbering_position <= 38:
        return "CDR1"
    if numbering_position <= 55:
        return "FR2"
    if numbering_position <= 65:
        return "CDR2"
    if numbering_position <= 104:
        return "FR3"
    if numbering_position <= 117:
        return "CDR3"
    return "FR4"


def audit_numbering_result(
    record: InputSequence, raw_result: Mapping[str, object]
) -> NumberingAudit:
    """Validate and map one raw ANARCII result without changing scientific scope.

    Args:
        record: Authoritative literal input sequence and provenance.
        raw_result: One ANARCII result mapping in IMGT scheme.

    Returns:
        A ``NumberingAudit``.  Successful results map every non-gap ANARCII
        residue to its zero- and one-based input-sequence index.  Failed results
        retain the tool error and contain no invented positions or span.

    Scope:
        ``query_start`` and ``query_end`` are validated as zero-based inclusive
        indices.  The reconstructed span must exactly equal the input slice.
        The function does not reinterpret chain type ``L`` or failure marker
        ``F`` as experimentally established molecular identity.
    """

    scheme = _required_text(raw_result.get("scheme"), record.sample_uid, "scheme")
    if scheme != "imgt":
        raise SequenceNumberingError(
            f"ANARCII scheme must be imgt for {record.sample_uid}, found {scheme!r}"
        )
    chain_type = _required_text(
        raw_result.get("chain_type"), record.sample_uid, "chain_type"
    )
    score = _required_finite_float(raw_result.get("score"), record.sample_uid, "score")
    raw_numbering = raw_result.get("numbering")
    raw_error = raw_result.get("error")

    if raw_numbering is None:
        error = _required_text(raw_error, record.sample_uid, "error")
        if raw_result.get("query_start") is not None or raw_result.get("query_end") is not None:
            raise SequenceNumberingError(
                f"Failed ANARCII result has query bounds for {record.sample_uid}"
            )
        return NumberingAudit(
            input_record=record,
            numbering_status="failed",
            sequence_scope_status="unresolved",
            chain_type=chain_type,
            score=score,
            scheme=scheme,
            query_start_0based_inclusive=None,
            query_end_0based_inclusive=None,
            provisional_numbered_span_sequence="",
            provisional_numbered_span_sha256="",
            unnumbered_n_sequence="",
            unnumbered_c_sequence="",
            first_numbered_imgt_position="",
            last_numbered_imgt_position="",
            numbering_review_flags=("anarcii_numbering_failed",),
            error=error,
            positions=(),
        )

    if raw_error not in (None, ""):
        raise SequenceNumberingError(
            f"Successful ANARCII numbering also has an error for {record.sample_uid}: "
            f"{raw_error!r}"
        )
    if not isinstance(raw_numbering, Sequence) or isinstance(raw_numbering, (str, bytes)):
        raise SequenceNumberingError(
            f"ANARCII numbering must be a sequence for {record.sample_uid}"
        )

    query_start = _required_index(
        raw_result.get("query_start"), record.sample_uid, "query_start"
    )
    query_end = _required_index(raw_result.get("query_end"), record.sample_uid, "query_end")
    if not 0 <= query_start <= query_end < record.sequence_length_aa:
        raise SequenceNumberingError(
            f"ANARCII query bounds are outside {record.sample_uid}: "
            f"{query_start}..{query_end} for length {record.sequence_length_aa}"
        )

    positions: list[NumberingPosition] = []
    seen_labels: set[tuple[int, str]] = set()
    non_gap_residues: list[str] = []
    next_sequence_index = query_start
    for order, raw_position in enumerate(raw_numbering, start=1):
        if not isinstance(raw_position, Sequence) or isinstance(raw_position, (str, bytes)):
            raise SequenceNumberingError(
                f"Malformed ANARCII position {order} for {record.sample_uid}"
            )
        if len(raw_position) != 2:
            raise SequenceNumberingError(
                f"Malformed ANARCII position {order} for {record.sample_uid}"
            )
        raw_label, raw_residue = raw_position
        if not isinstance(raw_label, Sequence) or isinstance(raw_label, (str, bytes)):
            raise SequenceNumberingError(
                f"Malformed ANARCII label {order} for {record.sample_uid}"
            )
        if len(raw_label) != 2:
            raise SequenceNumberingError(
                f"Malformed ANARCII label {order} for {record.sample_uid}"
            )
        numbering_position = _required_index(
            raw_label[0], record.sample_uid, f"numbering[{order}].position"
        )
        insertion_raw = raw_label[1]
        if not isinstance(insertion_raw, str):
            raise SequenceNumberingError(
                f"Non-text insertion code at position {order} for {record.sample_uid}"
            )
        insertion_code = insertion_raw.strip()
        if insertion_code and (
            len(insertion_code) != 1 or not insertion_code.isalpha() or not insertion_code.isupper()
        ):
            raise SequenceNumberingError(
                f"Invalid insertion code {insertion_raw!r} for {record.sample_uid}"
            )
        label_key = (numbering_position, insertion_code)
        if label_key in seen_labels:
            raise SequenceNumberingError(
                f"Duplicate IMGT label {numbering_position}{insertion_code} "
                f"for {record.sample_uid}"
            )
        seen_labels.add(label_key)

        if not isinstance(raw_residue, str) or len(raw_residue) != 1:
            raise SequenceNumberingError(
                f"Invalid residue at ANARCII position {order} for {record.sample_uid}"
            )
        if raw_residue != "-" and raw_residue not in STANDARD_AMINO_ACIDS:
            raise SequenceNumberingError(
                f"Unsupported ANARCII residue {raw_residue!r} for {record.sample_uid}"
            )
        sequence_index: int | None = None
        if raw_residue != "-":
            sequence_index = next_sequence_index
            next_sequence_index += 1
            non_gap_residues.append(raw_residue)

        positions.append(
            NumberingPosition(
                position_order=order,
                numbering_position=numbering_position,
                insertion_code=insertion_code,
                residue_aa=raw_residue,
                sequence_index_0based=sequence_index,
                region=imgt_region(numbering_position),
            )
        )

    numbered_span = "".join(non_gap_residues)
    expected_span = record.sequence_raw[query_start : query_end + 1]
    if next_sequence_index != query_end + 1 or numbered_span != expected_span:
        raise SequenceNumberingError(
            f"ANARCII numbering does not reconstruct inclusive query span for "
            f"{record.sample_uid}"
        )
    non_gap_positions = [position for position in positions if not position.is_gap]
    if not non_gap_positions:
        raise SequenceNumberingError(
            f"Successful ANARCII result has no residues for {record.sample_uid}"
        )

    first_label = non_gap_positions[0].label
    last_label = non_gap_positions[-1].label
    review_flags: list[str] = []
    if chain_type != "H":
        review_flags.append(f"anarcii_chain_type_{chain_type}")
    if query_start:
        review_flags.append("unnumbered_n_terminal_residues")
    if query_end != record.sequence_length_aa - 1:
        review_flags.append("unnumbered_c_terminal_residues")
    if first_label != "1":
        review_flags.append("first_numbered_imgt_position_not_1")
    if last_label != "128":
        review_flags.append("last_numbered_imgt_position_not_128")

    return NumberingAudit(
        input_record=record,
        numbering_status="pass",
        sequence_scope_status="provisional_numbered_domain",
        chain_type=chain_type,
        score=score,
        scheme=scheme,
        query_start_0based_inclusive=query_start,
        query_end_0based_inclusive=query_end,
        provisional_numbered_span_sequence=numbered_span,
        provisional_numbered_span_sha256=sha256_bytes(numbered_span.encode("ascii")),
        unnumbered_n_sequence=record.sequence_raw[:query_start],
        unnumbered_c_sequence=record.sequence_raw[query_end + 1 :],
        first_numbered_imgt_position=first_label,
        last_numbered_imgt_position=last_label,
        numbering_review_flags=tuple(review_flags),
        error="",
        positions=tuple(positions),
    )


def build_numbering_audits(
    records: Sequence[InputSequence], raw_results: Mapping[str, object]
) -> tuple[NumberingAudit, ...]:
    """Validate one and only one raw ANARCII result for every input record."""

    expected_keys = [record.sample_uid for record in records]
    if set(raw_results) != set(expected_keys) or len(raw_results) != len(expected_keys):
        missing = sorted(set(expected_keys) - set(raw_results))
        unexpected = sorted(set(raw_results) - set(expected_keys))
        raise SequenceNumberingError(
            "ANARCII output identities differ from input identities; "
            f"missing={missing}, unexpected={unexpected}"
        )
    audits: list[NumberingAudit] = []
    for record in records:
        raw_result = raw_results[record.sample_uid]
        if not isinstance(raw_result, Mapping):
            raise SequenceNumberingError(
                f"ANARCII result is not a mapping for {record.sample_uid}"
            )
        audits.append(audit_numbering_result(record, raw_result))
    return tuple(audits)


def validate_expected_baseline_outcome(
    audits: Sequence[NumberingAudit],
) -> dict[str, object]:
    """Apply the fixed acceptance gate for the validated 47-sequence baseline.

    Inputs are already provenance-checked ``NumberingAudit`` objects.  The
    function returns concrete categorical/boundary facts after requiring 46
    passes, the single observed ``WCC__4-28`` failure, successful-chain counts
    H=45/L=1, and the Nb252 zero-based inclusive bounds 0..125 with IMGT span
    1..128.  It intentionally does not assert an ANARCII score, infer that the
    L result is a biological light chain, or upgrade provisional spans to
    authoritative VHH boundaries.
    """

    if len(audits) != EXPECTED_RECORD_COUNT:
        raise SequenceNumberingError(
            f"Baseline acceptance requires {EXPECTED_RECORD_COUNT} audits"
        )
    status_counts = Counter(audit.numbering_status for audit in audits)
    if status_counts != Counter({"pass": 46, "failed": 1}):
        raise SequenceNumberingError(
            f"Unexpected baseline numbering counts: {dict(status_counts)}"
        )
    failed = [audit for audit in audits if audit.numbering_status == "failed"]
    if [audit.input_record.sample_uid for audit in failed] != ["WCC__4-28"]:
        raise SequenceNumberingError(
            "Expected WCC__4-28 to be the only failed baseline numbering"
        )
    success_chain_counts = Counter(
        audit.chain_type for audit in audits if audit.numbering_status == "pass"
    )
    if success_chain_counts != Counter({"H": 45, "L": 1}):
        raise SequenceNumberingError(
            f"Unexpected successful chain-type counts: {dict(success_chain_counts)}"
        )

    matches = [
        audit for audit in audits if audit.input_record.sample_uid == "LTT__Nb252"
    ]
    if len(matches) != 1:
        raise SequenceNumberingError("Baseline must contain exactly one LTT__Nb252")
    nb252 = matches[0]
    nb252_non_gap_count = sum(not position.is_gap for position in nb252.positions)
    observed_nb252 = (
        nb252.numbering_status,
        nb252.chain_type,
        nb252.scheme,
        nb252.query_start_0based_inclusive,
        nb252.query_end_0based_inclusive,
        nb252_non_gap_count,
        nb252.first_numbered_imgt_position,
        nb252.last_numbered_imgt_position,
        nb252.unnumbered_n_sequence,
        nb252.unnumbered_c_sequence,
    )
    expected_nb252 = (
        "pass",
        "H",
        "imgt",
        0,
        125,
        126,
        "1",
        "128",
        "",
        "GS",
    )
    if observed_nb252 != expected_nb252:
        raise SequenceNumberingError(
            f"Unexpected Nb252 numbering boundary facts: {observed_nb252!r}"
        )

    return {
        "status": "pass",
        "expected_record_count": EXPECTED_RECORD_COUNT,
        "numbering_pass_count": 46,
        "numbering_failed_count": 1,
        "failed_sample_uid": "WCC__4-28",
        "successful_chain_type_counts": {"H": 45, "L": 1},
        "nb252": {
            "sample_uid": "LTT__Nb252",
            "numbering_status": "pass",
            "chain_type": "H",
            "scheme": "imgt",
            "query_start_0based_inclusive": 0,
            "query_end_0based_inclusive": 125,
            "numbered_non_gap_count": 126,
            "first_numbered_imgt_position": "1",
            "last_numbered_imgt_position": "128",
            "unnumbered_n_sequence": "",
            "unnumbered_c_sequence": "GS",
            "provisional_numbered_span_sha256": nb252.provisional_numbered_span_sha256,
        },
        "score_assertion_policy": "scores recorded but no exact score asserted",
    }


def _row_text(row: Mapping[str, object], field: str) -> str:
    value = row.get(field, "")
    if value is None:
        return ""
    return str(value)


def _required_text(value: object, sample_uid: str, field: str) -> str:
    if not isinstance(value, str) or not value:
        raise SequenceNumberingError(
            f"ANARCII {field} must be non-empty text for {sample_uid}"
        )
    return value


def _required_index(value: object, sample_uid: str, field: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise SequenceNumberingError(
            f"ANARCII {field} must be an integer for {sample_uid}: {value!r}"
        )
    return value


def _required_finite_float(value: object, sample_uid: str, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise SequenceNumberingError(
            f"ANARCII {field} must be numeric for {sample_uid}: {value!r}"
        )
    result = float(value)
    if not math.isfinite(result):
        raise SequenceNumberingError(
            f"ANARCII {field} must be finite for {sample_uid}: {value!r}"
        )
    return result
