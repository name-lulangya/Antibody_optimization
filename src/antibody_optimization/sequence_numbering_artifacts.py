"""Artifact I/O and compact summaries for provisional sequence numbering.

The core :mod:`antibody_optimization.sequence_numbering` module owns sequence
and ANARCII semantic validation.  This module owns the stable tabular/JSON
handoff contract, upstream manifest verification, and deterministic UTF-8
serialization.  It does not invoke ANARCII or infer sequence boundaries.
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

from .sequence_numbering import (
    EXPECTED_RECORD_COUNT,
    InputSequence,
    NumberingAudit,
    SequenceNumberingError,
    sha256_bytes,
    validate_sequence_rows,
)


INPUT_TABLE_NAME = "nb_expression_records.csv"

SEQUENCE_REVIEW_FIELDS = (
    "sample_uid",
    "provider_code",
    "source_sample_id",
    "sequence_raw",
    "sequence_length_aa",
    "sequence_sha256",
    "source_sequence_scope",
    "source_sequence_review_flags",
    "numbering_status",
    "sequence_scope_status",
    "provisional_numbered_span_sequence",
    "provisional_numbered_span_sha256",
    "chain_type",
    "score",
    "scheme",
    "query_start_0based_inclusive",
    "query_end_0based_inclusive",
    "numbered_span_length_aa",
    "numbered_non_gap_count",
    "numbering_position_row_count",
    "unnumbered_n_length_aa",
    "unnumbered_n_sequence",
    "unnumbered_c_length_aa",
    "unnumbered_c_sequence",
    "first_numbered_imgt_position",
    "last_numbered_imgt_position",
    "numbering_review_flags",
    "error",
)

POSITION_FIELDS = (
    "sample_uid",
    "provider_code",
    "source_sample_id",
    "sequence_sha256",
    "provisional_numbered_span_sha256",
    "numbering_status",
    "sequence_scope_status",
    "chain_type",
    "scheme",
    "query_start_0based_inclusive",
    "query_end_0based_inclusive",
    "position_order",
    "numbering_position",
    "insertion_code",
    "numbering_position_label",
    "region",
    "residue_aa",
    "is_gap",
    "sequence_index_0based",
    "sequence_index_1based",
)


def sha256_file(path: Path) -> str:
    """Return the lowercase SHA-256 digest of a concrete file."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_manifest_document(
    manifest: Mapping[str, object],
    *,
    records_size_bytes: int,
    records_sha256: str,
    csv_fieldnames: Sequence[str],
    expected_count: int = EXPECTED_RECORD_COUNT,
) -> dict[str, object]:
    """Validate the upstream manifest against the exact records-table bytes.

    Args:
        manifest: Parsed upstream JSON mapping.  It is treated as read-only.
        records_size_bytes: Byte size of the records CSV that will be numbered.
        records_sha256: SHA-256 of those exact CSV bytes.
        csv_fieldnames: CSV header in its observed order.
        expected_count: Required manifest record count.

    Returns:
        A compact mapping of validated upstream provenance facts.

    Scope:
        This pure function validates the existing manifest/table contract.  It
        does not inspect individual sequences, infer VHH boundaries, or read
        either file from disk.
    """

    record_count = manifest.get("record_count")
    if record_count != expected_count:
        raise SequenceNumberingError(
            f"Manifest record_count must be {expected_count}, found {record_count!r}"
        )

    outputs = _require_mapping(manifest.get("outputs"), "manifest.outputs")
    records_metadata = _require_mapping(
        outputs.get(INPUT_TABLE_NAME), f"manifest.outputs[{INPUT_TABLE_NAME!r}]"
    )
    if records_metadata.get("sha256") != records_sha256:
        raise SequenceNumberingError(
            "Input records SHA-256 does not match the upstream manifest"
        )
    if records_metadata.get("size_bytes") != records_size_bytes:
        raise SequenceNumberingError(
            "Input records byte size does not match the upstream manifest"
        )

    tables = _require_mapping(manifest.get("tables"), "manifest.tables")
    recorded_fields = tables.get(INPUT_TABLE_NAME)
    if not isinstance(recorded_fields, list) or recorded_fields != list(csv_fieldnames):
        raise SequenceNumberingError(
            "Input records header does not match the schema recorded in the manifest"
        )

    integrity = _require_mapping(
        manifest.get("sequence_integrity"), "manifest.sequence_integrity"
    )
    if integrity.get("per_record_sha256") is not True:
        raise SequenceNumberingError(
            "Upstream manifest does not attest per-record sequence SHA-256 values"
        )

    return {
        "dataset": manifest.get("dataset"),
        "parser_version": manifest.get("parser_version"),
        "record_count": record_count,
        "records_size_bytes": records_size_bytes,
        "records_sha256": records_sha256,
        "per_record_sha256": True,
    }


def load_validated_sequence_input(
    records_path: Path,
    manifest_path: Path,
    *,
    expected_count: int = EXPECTED_RECORD_COUNT,
) -> tuple[tuple[InputSequence, ...], dict[str, object]]:
    """Read and validate the upstream manifest plus records CSV."""

    records_bytes = records_path.read_bytes()
    records_sha256 = sha256_bytes(records_bytes)
    with records_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise SequenceNumberingError("Input records CSV has no header")
        fieldnames = tuple(reader.fieldnames)
        rows = list(reader)

    try:
        manifest_value = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
    except json.JSONDecodeError as exc:
        raise SequenceNumberingError(
            f"Input manifest is not valid JSON: {manifest_path}"
        ) from exc
    if not isinstance(manifest_value, dict):
        raise SequenceNumberingError("Input manifest JSON must contain an object")

    provenance = validate_manifest_document(
        manifest_value,
        records_size_bytes=len(records_bytes),
        records_sha256=records_sha256,
        csv_fieldnames=fieldnames,
        expected_count=expected_count,
    )
    records = validate_sequence_rows(fieldnames, rows, expected_count=expected_count)
    return records, provenance


def sequence_review_rows(audits: Sequence[NumberingAudit]) -> list[dict[str, object]]:
    """Build deterministic sequence-level CSV rows from validated audits."""

    rows: list[dict[str, object]] = []
    for audit in audits:
        record = audit.input_record
        rows.append(
            {
                "sample_uid": record.sample_uid,
                "provider_code": record.provider_code,
                "source_sample_id": record.source_sample_id,
                "sequence_raw": record.sequence_raw,
                "sequence_length_aa": record.sequence_length_aa,
                "sequence_sha256": record.sequence_sha256,
                "source_sequence_scope": record.source_sequence_scope,
                "source_sequence_review_flags": record.source_sequence_review_flags,
                "numbering_status": audit.numbering_status,
                "sequence_scope_status": audit.sequence_scope_status,
                "provisional_numbered_span_sequence": audit.provisional_numbered_span_sequence,
                "provisional_numbered_span_sha256": audit.provisional_numbered_span_sha256,
                "chain_type": audit.chain_type,
                "score": repr(audit.score),
                "scheme": audit.scheme,
                "query_start_0based_inclusive": _csv_optional_int(
                    audit.query_start_0based_inclusive
                ),
                "query_end_0based_inclusive": _csv_optional_int(
                    audit.query_end_0based_inclusive
                ),
                "numbered_span_length_aa": (
                    len(audit.provisional_numbered_span_sequence)
                    if audit.numbering_status == "pass"
                    else ""
                ),
                "numbered_non_gap_count": (
                    sum(not position.is_gap for position in audit.positions)
                    if audit.numbering_status == "pass"
                    else ""
                ),
                "numbering_position_row_count": (
                    len(audit.positions) if audit.numbering_status == "pass" else ""
                ),
                "unnumbered_n_length_aa": (
                    len(audit.unnumbered_n_sequence)
                    if audit.numbering_status == "pass"
                    else ""
                ),
                "unnumbered_n_sequence": audit.unnumbered_n_sequence,
                "unnumbered_c_length_aa": (
                    len(audit.unnumbered_c_sequence)
                    if audit.numbering_status == "pass"
                    else ""
                ),
                "unnumbered_c_sequence": audit.unnumbered_c_sequence,
                "first_numbered_imgt_position": audit.first_numbered_imgt_position,
                "last_numbered_imgt_position": audit.last_numbered_imgt_position,
                "numbering_review_flags": ";".join(audit.numbering_review_flags),
                "error": audit.error,
            }
        )
    return rows


def numbering_position_rows(
    audits: Sequence[NumberingAudit],
) -> list[dict[str, object]]:
    """Build gap-preserving IMGT position rows with authoritative input indices."""

    rows: list[dict[str, object]] = []
    for audit in audits:
        record = audit.input_record
        for position in audit.positions:
            rows.append(
                {
                    "sample_uid": record.sample_uid,
                    "provider_code": record.provider_code,
                    "source_sample_id": record.source_sample_id,
                    "sequence_sha256": record.sequence_sha256,
                    "provisional_numbered_span_sha256": audit.provisional_numbered_span_sha256,
                    "numbering_status": audit.numbering_status,
                    "sequence_scope_status": audit.sequence_scope_status,
                    "chain_type": audit.chain_type,
                    "scheme": audit.scheme,
                    "query_start_0based_inclusive": audit.query_start_0based_inclusive,
                    "query_end_0based_inclusive": audit.query_end_0based_inclusive,
                    "position_order": position.position_order,
                    "numbering_position": position.numbering_position,
                    "insertion_code": position.insertion_code,
                    "numbering_position_label": position.label,
                    "region": position.region,
                    "residue_aa": position.residue_aa,
                    "is_gap": str(position.is_gap).lower(),
                    "sequence_index_0based": _csv_optional_int(
                        position.sequence_index_0based
                    ),
                    "sequence_index_1based": (
                        position.sequence_index_0based + 1
                        if position.sequence_index_0based is not None
                        else ""
                    ),
                }
            )
    return rows


def sample_summaries(audits: Sequence[NumberingAudit]) -> list[dict[str, object]]:
    """Build compact JSON sample objects for downstream audit gates."""

    return [
        {
            "sample_uid": audit.input_record.sample_uid,
            "sequence_sha256": audit.input_record.sequence_sha256,
            "numbering_status": audit.numbering_status,
            "sequence_scope_status": audit.sequence_scope_status,
            "provisional_numbered_span_sha256": (
                audit.provisional_numbered_span_sha256 or None
            ),
            "chain_type": audit.chain_type,
            "scheme": audit.scheme,
            "query_start_0based_inclusive": audit.query_start_0based_inclusive,
            "query_end_0based_inclusive": audit.query_end_0based_inclusive,
            "numbering_review_flags": list(audit.numbering_review_flags),
            "error": audit.error or None,
        }
        for audit in audits
    ]


def result_statistics(audits: Sequence[NumberingAudit]) -> dict[str, object]:
    """Aggregate categorical and boundary facts from validated audit objects."""

    status_counts = Counter(audit.numbering_status for audit in audits)
    scope_counts = Counter(audit.sequence_scope_status for audit in audits)
    all_chain_counts = Counter(audit.chain_type for audit in audits)
    success_chain_counts = Counter(
        audit.chain_type for audit in audits if audit.numbering_status == "pass"
    )
    provider_status: dict[str, Counter[str]] = {}
    for audit in audits:
        provider_status.setdefault(audit.input_record.provider_code, Counter()).update(
            [audit.numbering_status]
        )
    successful = [audit for audit in audits if audit.numbering_status == "pass"]
    return {
        "record_count": len(audits),
        "numbering_status_counts": dict(sorted(status_counts.items())),
        "sequence_scope_status_counts": dict(sorted(scope_counts.items())),
        "chain_type_counts_all_results": dict(sorted(all_chain_counts.items())),
        "chain_type_counts_successful_results": dict(
            sorted(success_chain_counts.items())
        ),
        "provider_numbering_status_counts": {
            provider: dict(sorted(counts.items()))
            for provider, counts in sorted(provider_status.items())
        },
        "query_start_0based_inclusive_counts": _counter_of_optional_ints(
            audit.query_start_0based_inclusive for audit in successful
        ),
        "query_end_0based_inclusive_counts": _counter_of_optional_ints(
            audit.query_end_0based_inclusive for audit in successful
        ),
        "unnumbered_n_length_counts": _counter_of_optional_ints(
            len(audit.unnumbered_n_sequence) for audit in successful
        ),
        "unnumbered_c_length_counts": _counter_of_optional_ints(
            len(audit.unnumbered_c_sequence) for audit in successful
        ),
        "first_numbered_imgt_position_counts": dict(
            sorted(Counter(audit.first_numbered_imgt_position for audit in successful).items())
        ),
        "last_numbered_imgt_position_counts": dict(
            sorted(Counter(audit.last_numbered_imgt_position for audit in successful).items())
        ),
        "position_row_count": sum(len(audit.positions) for audit in audits),
        "numbered_non_gap_position_count": sum(
            sum(not position.is_gap for position in audit.positions) for audit in audits
        ),
        "failed_samples": [
            {
                "sample_uid": audit.input_record.sample_uid,
                "chain_type": audit.chain_type,
                "error": audit.error,
            }
            for audit in audits
            if audit.numbering_status == "failed"
        ],
        "multi_domain_output_status": "not_observed_with_scfv_false",
    }


def write_csv_utf8_bom(
    path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, object]]
) -> None:
    """Write a stable RFC-4180-style CSV with UTF-8 BOM and LF endings."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: object) -> None:
    """Write sorted, indented UTF-8 JSON with a final LF."""

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def _require_mapping(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise SequenceNumberingError(f"{label} must be a JSON object")
    return value


def _csv_optional_int(value: int | None) -> int | str:
    return "" if value is None else value


def _counter_of_optional_ints(values: Iterable[int | None]) -> dict[str, int]:
    counter = Counter("null" if value is None else str(value) for value in values)
    return dict(sorted(counter.items(), key=lambda item: item[0]))
