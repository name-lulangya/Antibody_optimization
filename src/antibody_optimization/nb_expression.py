"""Parse and validate collaborator nanobody-yield records from DOCX XML.

The source stores data as body paragraphs rather than Word tables. Raw paragraph
text is retained separately from whitespace-trimmed parser text. Amino-acid
sequences are accepted only when raw and parser text are identical, so sequence
whitespace is never normalized silently.
"""

from __future__ import annotations

import hashlib
import re
from collections import Counter
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path
from typing import Sequence
from xml.etree import ElementTree
from zipfile import BadZipFile, ZipFile


PARSER_VERSION = "1.1.0"
WORDPROCESSINGML_NS = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
STANDARD_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
EXPECTED_SECTION_COUNTS = {"LTT": 23, "WCC": 8, "LLJ": 16}
EXPECTED_SEMANTICS_COUNTS = {
    "individual_approximate": 31,
    "group_lower_bound": 9,
    "group_approximate": 7,
}
EXPECTED_NONEMPTY_PARAGRAPHS = 100

LTT_HEADER_RE = re.compile(
    r"^(?P<clone_id>[A-Za-z0-9-]+)\s+~\s*(?P<yield>[0-9]+(?:\.[0-9]+)?)\s*mg$"
)
WCC_HEADER_RE = re.compile(
    r"^>(?P<clone_id>[A-Za-z0-9-]+)\s+1L\s+TB纯化得到约\s*"
    r"(?P<yield>[0-9]+(?:\.[0-9]+)?)\s*mg$"
)
LLJ_BIN_RE = re.compile(
    r"^(?P<relation>[>~])\s*(?P<yield>[0-9]+(?:\.[0-9]+)?)\s*mg$"
)
CLONE_ID_RE = re.compile(r"^[A-Za-z0-9-]+$")


class ParseError(ValueError):
    """Raised when source content violates the verified document contract."""


@dataclass(frozen=True)
class DocxParagraph:
    """One nonempty body paragraph with both literal and parser text."""

    nonempty_index: int
    raw_text: str
    parse_text: str


@dataclass(frozen=True)
class ExpressionRecord:
    """One source sequence plus its scientifically explicit yield assignment."""

    sample_uid: str
    provider_code: str
    source_sample_id: str
    source_header_raw: str
    sequence_raw: str
    sequence_length_aa: int
    sequence_sha256: str
    sequence_scope: str
    sequence_review_flags: str
    observation_id: str
    assay_id: str
    reported_text: str
    observation_semantics: str
    value_relation: str
    point_estimate_mg: Decimal | None
    group_anchor_mg: Decimal | None
    lower_bound_mg: Decimal | None
    lower_bound_inclusive: bool | None
    assignment_level: str
    group_id: str
    individual_numeric_available: bool
    censoring_type: str
    censoring_reason: str
    culture_volume_l: Decimal | None
    volume_evidence: str
    medium: str
    yield_stage: str
    source_yield_paragraph_index: int
    source_clone_paragraph_index: int
    source_sequence_paragraph_index: int
    source_document: str
    source_sha256: str

    def samples_row(self) -> dict[str, str | int]:
        return {
            "sample_uid": self.sample_uid,
            "provider_code": self.provider_code,
            "source_sample_id": self.source_sample_id,
            "source_section": self.provider_code,
            "source_header_raw": self.source_header_raw,
            "sequence_raw": self.sequence_raw,
            "sequence_length_aa": self.sequence_length_aa,
            "sequence_sha256": self.sequence_sha256,
            "sequence_scope": self.sequence_scope,
            "vhh_region_sequence": "",
            "sequence_review_flags": self.sequence_review_flags,
            "source_document": self.source_document,
            "source_sha256": self.source_sha256,
            "source_clone_paragraph_index": self.source_clone_paragraph_index,
            "source_sequence_paragraph_index": self.source_sequence_paragraph_index,
        }

    def yield_row(self) -> dict[str, str | int | bool]:
        return {
            "observation_id": self.observation_id,
            "sample_uid": self.sample_uid,
            "assay_id": self.assay_id,
            "phenotype_name": "reported_yield",
            "reported_text": self.reported_text,
            "observation_semantics": self.observation_semantics,
            "value_relation": self.value_relation,
            "point_estimate_mg": optional_decimal(self.point_estimate_mg),
            "group_anchor_mg": optional_decimal(self.group_anchor_mg),
            "lower_bound_mg": optional_decimal(self.lower_bound_mg),
            "lower_bound_inclusive": (
                "" if self.lower_bound_inclusive is None else self.lower_bound_inclusive
            ),
            "upper_bound_mg": "",
            "assignment_level": self.assignment_level,
            "group_id": self.group_id,
            "individual_numeric_available": self.individual_numeric_available,
            "censoring_type": self.censoring_type,
            "censoring_reason": self.censoring_reason,
            "replicate_count": "",
            "uncertainty_value": "",
            "uncertainty_type": "",
            "source_yield_paragraph_index": self.source_yield_paragraph_index,
        }

    def wide_row(self) -> dict[str, str | int | bool]:
        row = self.samples_row()
        row.update(
            {key: value for key, value in self.yield_row().items() if key != "sample_uid"}
        )
        row.update(
            {
                "culture_volume_l": optional_decimal(self.culture_volume_l),
                "volume_evidence": self.volume_evidence,
                "medium": self.medium,
                "yield_stage": self.yield_stage,
            }
        )
        return row


def sha256_file(path: Path) -> str:
    """Return a lowercase SHA-256 hex digest."""

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def extract_nonempty_docx_paragraphs(source_path: Path) -> list[DocxParagraph]:
    """Extract literal nonempty body paragraphs and reject unsupported content.

    The verified collaborator format has no tables, text boxes, tabs, line
    breaks, tracked changes, or hidden runs. These constructs are rejected rather
    than silently skipped because any of them could hide or split a sequence.
    """

    try:
        with ZipFile(source_path) as archive:
            document_xml = archive.read("word/document.xml")
    except (BadZipFile, KeyError, OSError) as exc:
        raise ParseError(f"Cannot read DOCX body XML from {source_path}: {exc}") from exc
    try:
        root = ElementTree.fromstring(document_xml)
    except ElementTree.ParseError as exc:
        raise ParseError(f"Invalid word/document.xml in {source_path}: {exc}") from exc

    namespace = {"w": WORDPROCESSINGML_NS}
    unsupported = {
        "tables": root.findall(".//w:body//w:tbl", namespace),
        "text_boxes": root.findall(".//w:body//w:txbxContent", namespace),
        "tabs": root.findall(".//w:body//w:tab", namespace),
        "line_breaks": root.findall(".//w:body//w:br", namespace),
        "tracked_insertions": root.findall(".//w:body//w:ins", namespace),
        "tracked_deletions": root.findall(".//w:body//w:del", namespace),
        "hidden_text_markers": root.findall(".//w:body//w:vanish", namespace),
    }
    present = {name: len(nodes) for name, nodes in unsupported.items() if nodes}
    if present:
        raise ParseError(f"Unsupported DOCX structures are present: {present}")

    body = root.find(".//w:body", namespace)
    if body is None:
        raise ParseError("DOCX has no Word body")
    paragraphs: list[DocxParagraph] = []
    for paragraph in body.findall("w:p", namespace):
        raw_text = "".join(
            node.text or "" for node in paragraph.findall(".//w:t", namespace)
        )
        parse_text = raw_text.strip()
        if parse_text:
            paragraphs.append(
                DocxParagraph(
                    nonempty_index=len(paragraphs) + 1,
                    raw_text=raw_text,
                    parse_text=parse_text,
                )
            )
    return paragraphs


def parse_expression_docx(
    source_path: Path,
    *,
    document_title_culture_volume_l: Decimal | None,
) -> tuple[list[ExpressionRecord], list[DocxParagraph]]:
    """Parse all provider sections using explicit document-title metadata.

    ``document_title_culture_volume_l`` records the 1 L information conveyed by
    the original title and must be supplied independently of the mutable file
    path. Passing ``None`` leaves LTT/LLJ volume unknown. WCC volume remains 1 L
    because each WCC entry states it directly.
    """

    source_path = Path(source_path)
    paragraphs = extract_nonempty_docx_paragraphs(source_path)
    source_sha256 = sha256_file(source_path)
    records: list[ExpressionRecord] = []
    section: str | None = None
    llj_bin: tuple[str, Decimal, str, int] | None = None
    index = 0

    while index < len(paragraphs):
        paragraph = paragraphs[index]
        text = paragraph.parse_text
        if text.startswith("Data from "):
            section = text.removeprefix("Data from ")
            if section not in EXPECTED_SECTION_COUNTS:
                raise ParseError(f"Unknown section {section!r} at paragraph {index + 1}")
            llj_bin = None
            index += 1
            continue
        if section is None:
            raise ParseError(f"Record before section heading at paragraph {index + 1}")

        if section == "LLJ":
            bin_match = LLJ_BIN_RE.fullmatch(text)
            if bin_match:
                relation = "gt" if bin_match.group("relation") == ">" else "approx"
                llj_bin = (
                    relation,
                    Decimal(bin_match.group("yield")),
                    paragraph.raw_text,
                    paragraph.nonempty_index,
                )
                index += 1
                continue
            if llj_bin is None:
                raise ParseError(f"LLJ clone without yield bin at paragraph {index + 1}")
            if not CLONE_ID_RE.fullmatch(text):
                raise ParseError(f"Invalid LLJ clone ID {text!r} at paragraph {index + 1}")
            sequence, sequence_index = _next_sequence(paragraphs, index + 1)
            relation, value, reported_text, yield_index = llj_bin
            is_lower_bound = relation == "gt"
            records.append(
                _make_record(
                    provider_code=section,
                    source_sample_id=text,
                    source_header_raw=paragraph.raw_text,
                    sequence=sequence,
                    reported_text=reported_text,
                    observation_semantics=(
                        "group_lower_bound" if is_lower_bound else "group_approximate"
                    ),
                    value_relation=relation,
                    point_estimate_mg=None,
                    group_anchor_mg=None if is_lower_bound else value,
                    lower_bound_mg=value if is_lower_bound else None,
                    lower_bound_inclusive=False if is_lower_bound else None,
                    assignment_level="group",
                    group_id=(
                        f"LLJ_GT{format_decimal(value)}"
                        if is_lower_bound
                        else f"LLJ_APPROX{format_decimal(value)}"
                    ),
                    individual_numeric_available=False,
                    censoring_type="right_censored" if is_lower_bound else "",
                    censoring_reason="unknown" if is_lower_bound else "",
                    culture_volume_l=document_title_culture_volume_l,
                    volume_evidence=(
                        "document_title"
                        if document_title_culture_volume_l is not None
                        else "unavailable"
                    ),
                    medium="",
                    yield_stage="",
                    source_yield_paragraph_index=yield_index,
                    source_clone_paragraph_index=paragraph.nonempty_index,
                    source_sequence_paragraph_index=sequence_index,
                    source_path=source_path,
                    source_sha256=source_sha256,
                )
            )
            index += 2
            continue

        header_match = (
            LTT_HEADER_RE.fullmatch(text)
            if section == "LTT"
            else WCC_HEADER_RE.fullmatch(text)
        )
        if header_match is None:
            raise ParseError(f"Invalid {section} header {text!r} at paragraph {index + 1}")
        sequence, sequence_index = _next_sequence(paragraphs, index + 1)
        records.append(
            _make_record(
                provider_code=section,
                source_sample_id=header_match.group("clone_id"),
                source_header_raw=paragraph.raw_text,
                sequence=sequence,
                reported_text=paragraph.raw_text,
                observation_semantics="individual_approximate",
                value_relation="approx",
                point_estimate_mg=Decimal(header_match.group("yield")),
                group_anchor_mg=None,
                lower_bound_mg=None,
                lower_bound_inclusive=None,
                assignment_level="individual",
                group_id="",
                individual_numeric_available=True,
                censoring_type="",
                censoring_reason="",
                culture_volume_l=(
                    Decimal("1")
                    if section == "WCC"
                    else document_title_culture_volume_l
                ),
                volume_evidence=(
                    "entry_text"
                    if section == "WCC"
                    else (
                        "document_title"
                        if document_title_culture_volume_l is not None
                        else "unavailable"
                    )
                ),
                medium="TB" if section == "WCC" else "",
                yield_stage="post_purification" if section == "WCC" else "",
                source_yield_paragraph_index=paragraph.nonempty_index,
                source_clone_paragraph_index=paragraph.nonempty_index,
                source_sequence_paragraph_index=sequence_index,
                source_path=source_path,
                source_sha256=source_sha256,
            )
        )
        index += 2

    return records, paragraphs


def validate_records(
    records: Sequence[ExpressionRecord],
    paragraphs: Sequence[DocxParagraph],
    *,
    expected_source_sha256: str | None = None,
) -> dict[str, object]:
    """Fail on any count, identity, alphabet, hash, or raw-text mismatch."""

    if not records:
        raise ParseError("No expression records were parsed")
    source_hashes = {record.source_sha256 for record in records}
    if len(source_hashes) != 1:
        raise ParseError(f"Multiple source hashes in records: {sorted(source_hashes)}")
    source_sha256 = next(iter(source_hashes))
    if expected_source_sha256 and source_sha256.lower() != expected_source_sha256.lower():
        raise ParseError(
            f"Source SHA-256 mismatch: expected {expected_source_sha256}, got {source_sha256}"
        )
    section_counts = Counter(record.provider_code for record in records)
    if dict(section_counts) != EXPECTED_SECTION_COUNTS:
        raise ParseError(
            f"Section counts mismatch: expected {EXPECTED_SECTION_COUNTS}, got {dict(section_counts)}"
        )
    if len(paragraphs) != EXPECTED_NONEMPTY_PARAGRAPHS:
        raise ParseError(
            f"Nonempty paragraph count mismatch: expected {EXPECTED_NONEMPTY_PARAGRAPHS}, "
            f"got {len(paragraphs)}"
        )

    uids = [record.sample_uid for record in records]
    duplicate_uids = sorted(uid for uid, count in Counter(uids).items() if count > 1)
    if duplicate_uids:
        raise ParseError(f"Duplicate sample UIDs: {duplicate_uids}")
    sequence_hashes = [record.sequence_sha256 for record in records]
    duplicate_sequences = sorted(
        value for value, count in Counter(sequence_hashes).items() if count > 1
    )
    if duplicate_sequences:
        raise ParseError(f"Duplicate sequences: {duplicate_sequences}")

    for record in records:
        _validate_sequence(record.sequence_raw, record.source_sequence_paragraph_index)
        calculated_hash = hashlib.sha256(record.sequence_raw.encode("ascii")).hexdigest()
        if calculated_hash != record.sequence_sha256:
            raise ParseError(f"Sequence hash mismatch for {record.sample_uid}")
        if len(record.sequence_raw) != record.sequence_length_aa:
            raise ParseError(f"Sequence length mismatch for {record.sample_uid}")
        sequence_paragraph = paragraphs[record.source_sequence_paragraph_index - 1]
        if sequence_paragraph.raw_text != record.sequence_raw:
            raise ParseError(f"Raw source sequence mismatch for {record.sample_uid}")
        yield_paragraph = paragraphs[record.source_yield_paragraph_index - 1]
        if yield_paragraph.raw_text != record.reported_text:
            raise ParseError(f"Raw source yield-text mismatch for {record.sample_uid}")
        clone_paragraph = paragraphs[record.source_clone_paragraph_index - 1]
        if clone_paragraph.raw_text != record.source_header_raw:
            raise ParseError(f"Raw source header mismatch for {record.sample_uid}")
        if record.provider_code == "LLJ":
            if clone_paragraph.parse_text != record.source_sample_id:
                raise ParseError(f"Source clone-ID mismatch for {record.sample_uid}")
        elif record.source_sample_id not in clone_paragraph.parse_text:
            raise ParseError(f"Source header lost clone ID for {record.sample_uid}")
        if record.provider_code == "LLJ" and record.point_estimate_mg is not None:
            raise ParseError(f"LLJ point estimate must be blank for {record.sample_uid}")

    semantics_counts = Counter(record.observation_semantics for record in records)
    if dict(semantics_counts) != EXPECTED_SEMANTICS_COUNTS:
        raise ParseError(
            f"Semantics counts mismatch: expected {EXPECTED_SEMANTICS_COUNTS}, "
            f"got {dict(semantics_counts)}"
        )
    nb252 = [record for record in records if record.sample_uid == "LTT__Nb252"]
    if len(nb252) != 1 or nb252[0].point_estimate_mg != Decimal("0.5"):
        raise ParseError("Expected one LTT__Nb252 record with individual ~0.5 mg")

    flag_counts = Counter(
        flag
        for record in records
        for flag in record.sequence_review_flags.split(";")
        if flag
    )
    expected_flag_counts = {
        "short_literal_sequence_lt115": 6,
        "single_cysteine_literal_sequence": 2,
        "wgqgt_motif_absent": 9,
    }
    if dict(flag_counts) != expected_flag_counts:
        raise ParseError(
            f"Sequence flag counts mismatch: expected {expected_flag_counts}, "
            f"got {dict(flag_counts)}"
        )

    whitespace_paragraphs = [
        paragraph.nonempty_index
        for paragraph in paragraphs
        if paragraph.raw_text != paragraph.parse_text
    ]
    return {
        "status": "pass",
        "parser_version": PARSER_VERSION,
        "source_file": records[0].source_document,
        "source_file_sha256": source_sha256,
        "source_nonempty_paragraph_count": len(paragraphs),
        "paragraphs_with_outer_whitespace": whitespace_paragraphs,
        "record_count": len(records),
        "record_counts_by_source": dict(section_counts),
        "yield_semantics_counts": dict(semantics_counts),
        "sequence_review_flag_counts": dict(flag_counts),
        "unique_sample_uid_count": len(set(uids)),
        "unique_sequence_count": len(set(sequence_hashes)),
        "invalid_sequence_count": 0,
        "sequence_round_trip_mismatch_count": 0,
        "source_text_round_trip_mismatch_count": 0,
        "checks": [
            "unsupported_docx_structure_absence",
            "source_sha256",
            "nonempty_paragraph_count",
            "section_record_counts",
            "sample_uid_uniqueness",
            "sequence_uniqueness",
            "standard_amino_acid_alphabet",
            "sequence_raw_equals_parser_text",
            "sequence_length",
            "per_sequence_sha256",
            "raw_source_paragraph_round_trip",
            "yield_semantics_counts",
            "llj_numeric_field_separation",
            "target_clone_nb252_presence",
            "sequence_review_flag_counts",
        ],
    }


def format_decimal(value: Decimal) -> str:
    text = format(value, "f")
    return text.rstrip("0").rstrip(".") if "." in text else text


def optional_decimal(value: Decimal | None) -> str:
    return "" if value is None else format_decimal(value)


def _make_record(
    *,
    provider_code: str,
    source_sample_id: str,
    source_header_raw: str,
    sequence: str,
    reported_text: str,
    observation_semantics: str,
    value_relation: str,
    point_estimate_mg: Decimal | None,
    group_anchor_mg: Decimal | None,
    lower_bound_mg: Decimal | None,
    lower_bound_inclusive: bool | None,
    assignment_level: str,
    group_id: str,
    individual_numeric_available: bool,
    censoring_type: str,
    censoring_reason: str,
    culture_volume_l: Decimal | None,
    volume_evidence: str,
    medium: str,
    yield_stage: str,
    source_yield_paragraph_index: int,
    source_clone_paragraph_index: int,
    source_sequence_paragraph_index: int,
    source_path: Path,
    source_sha256: str,
) -> ExpressionRecord:
    sample_uid = f"{provider_code}__{source_sample_id}"
    return ExpressionRecord(
        sample_uid=sample_uid,
        provider_code=provider_code,
        source_sample_id=source_sample_id,
        source_header_raw=source_header_raw,
        sequence_raw=sequence,
        sequence_length_aa=len(sequence),
        sequence_sha256=hashlib.sha256(sequence.encode("ascii")).hexdigest(),
        sequence_scope="unknown",
        sequence_review_flags=_sequence_review_flags(sequence),
        observation_id=f"OBS__{sample_uid}",
        assay_id=f"ASSAY__{provider_code}",
        reported_text=reported_text,
        observation_semantics=observation_semantics,
        value_relation=value_relation,
        point_estimate_mg=point_estimate_mg,
        group_anchor_mg=group_anchor_mg,
        lower_bound_mg=lower_bound_mg,
        lower_bound_inclusive=lower_bound_inclusive,
        assignment_level=assignment_level,
        group_id=group_id,
        individual_numeric_available=individual_numeric_available,
        censoring_type=censoring_type,
        censoring_reason=censoring_reason,
        culture_volume_l=culture_volume_l,
        volume_evidence=volume_evidence,
        medium=medium,
        yield_stage=yield_stage,
        source_yield_paragraph_index=source_yield_paragraph_index,
        source_clone_paragraph_index=source_clone_paragraph_index,
        source_sequence_paragraph_index=source_sequence_paragraph_index,
        source_document=source_path.name,
        source_sha256=source_sha256,
    )


def _next_sequence(
    paragraphs: Sequence[DocxParagraph], next_index: int
) -> tuple[str, int]:
    if next_index >= len(paragraphs):
        raise ParseError("Record header is missing its following sequence paragraph")
    paragraph = paragraphs[next_index]
    if paragraph.raw_text != paragraph.parse_text:
        raise ParseError(
            f"Sequence paragraph {paragraph.nonempty_index} has outer whitespace; "
            "refusing to normalize it"
        )
    _validate_sequence(paragraph.raw_text, paragraph.nonempty_index)
    return paragraph.raw_text, paragraph.nonempty_index


def _validate_sequence(sequence: str, paragraph_index: int) -> None:
    invalid = sorted(set(sequence) - STANDARD_AMINO_ACIDS)
    if not sequence or invalid:
        raise ParseError(
            f"Invalid sequence at nonempty paragraph {paragraph_index}; characters={invalid}"
        )


def _sequence_review_flags(sequence: str) -> str:
    flags: list[str] = []
    if len(sequence) < 115:
        flags.append("short_literal_sequence_lt115")
    if sequence.count("C") == 1:
        flags.append("single_cysteine_literal_sequence")
    if "WGQGT" not in sequence:
        flags.append("wgqgt_motif_absent")
    return ";".join(flags)
