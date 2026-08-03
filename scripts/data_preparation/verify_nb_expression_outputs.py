#!/usr/bin/env python3
"""Independently verify every exported sequence and yield field against DOCX XML.

This verifier intentionally does not import the production parser or artifact
writers. Its separate state machine guards against shared sequence-copy errors.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
import tempfile
from collections import Counter
from decimal import Decimal
from pathlib import Path
from xml.etree import ElementTree
from zipfile import ZipFile


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from antibody_optimization.file_transaction import validate_file_paths  # noqa: E402


AA_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")
LTT_RE = re.compile(
    r"^(?P<id>[A-Za-z0-9-]+)\s+~\s*(?P<value>[0-9]+(?:\.[0-9]+)?)\s*mg$"
)
WCC_RE = re.compile(
    r"^>(?P<id>[A-Za-z0-9-]+)\s+1L\s+TB纯化得到约\s*"
    r"(?P<value>[0-9]+(?:\.[0-9]+)?)\s*mg$"
)
LLJ_RE = re.compile(r"^(?P<relation>[>~])\s*(?P<value>[0-9]+(?:\.[0-9]+)?)\s*mg$")
EXPECTED_COUNTS = {"LTT": 23, "WCC": 8, "LLJ": 16}


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--artifact-dir", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--expected-source-sha256", required=True)
    parser.add_argument(
        "--document-title-culture-volume-l",
        required=True,
        type=Decimal,
        help="Explicit metadata read from the original document title",
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def source_records(
    path: Path, title_volume_l: Decimal
) -> tuple[list[dict[str, str]], list[dict[str, str]], str]:
    with ZipFile(path) as archive:
        document_xml = archive.read("word/document.xml")
    namespace = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
    root = ElementTree.fromstring(document_xml)
    unsupported = {
        "tables": root.findall(".//w:body//w:tbl", namespace),
        "text_boxes": root.findall(".//w:body//w:txbxContent", namespace),
        "tabs": root.findall(".//w:body//w:tab", namespace),
        "line_breaks": root.findall(".//w:body//w:br", namespace),
        "tracked_insertions": root.findall(".//w:body//w:ins", namespace),
        "tracked_deletions": root.findall(".//w:body//w:del", namespace),
        "hidden_text": root.findall(".//w:body//w:vanish", namespace),
    }
    present = {name: len(nodes) for name, nodes in unsupported.items() if nodes}
    if present:
        raise AssertionError(f"Unsupported DOCX structures: {present}")
    body = root.find(".//w:body", namespace)
    if body is None:
        raise AssertionError("DOCX has no Word body")
    paragraphs: list[dict[str, str]] = []
    for paragraph in body.findall("w:p", namespace):
        raw_text = "".join(
            node.text or "" for node in paragraph.findall(".//w:t", namespace)
        )
        parse_text = raw_text.strip()
        if parse_text:
            paragraphs.append({"raw_text": raw_text, "parse_text": parse_text})

    records: list[dict[str, str]] = []
    section = ""
    llj_bin: tuple[str, str, str] | None = None
    index = 0
    while index < len(paragraphs):
        paragraph = paragraphs[index]
        text = paragraph["parse_text"]
        if text.startswith("Data from "):
            section = text[10:]
            llj_bin = None
            index += 1
            continue
        if section == "LLJ":
            bin_match = LLJ_RE.fullmatch(text)
            if bin_match:
                llj_bin = (
                    bin_match.group("relation"),
                    bin_match.group("value"),
                    paragraph["raw_text"],
                )
                index += 1
                continue
            if llj_bin is None:
                raise AssertionError(f"LLJ clone lacks bin at paragraph {index + 1}")
            clone_id = text
            source_header_raw = paragraph["raw_text"]
            sequence_paragraph = paragraphs[index + 1]
            relation_symbol, value, reported_text = llj_bin
            semantics = "group_lower_bound" if relation_symbol == ">" else "group_approximate"
            point_estimate = ""
            group_anchor = value if relation_symbol == "~" else ""
            lower_bound = value if relation_symbol == ">" else ""
            relation = "gt" if relation_symbol == ">" else "approx"
            culture_volume_l = format(title_volume_l, "f")
            volume_evidence = "document_title"
            medium = ""
            yield_stage = ""
            index += 2
        elif section in {"LTT", "WCC"}:
            match = LTT_RE.fullmatch(text) if section == "LTT" else WCC_RE.fullmatch(text)
            if match is None:
                raise AssertionError(f"Invalid {section} header at paragraph {index + 1}: {text}")
            clone_id = match.group("id")
            source_header_raw = paragraph["raw_text"]
            sequence_paragraph = paragraphs[index + 1]
            reported_text = paragraph["raw_text"]
            semantics = "individual_approximate"
            point_estimate = match.group("value")
            group_anchor = ""
            lower_bound = ""
            relation = "approx"
            culture_volume_l = "1" if section == "WCC" else format(title_volume_l, "f")
            volume_evidence = "entry_text" if section == "WCC" else "document_title"
            medium = "TB" if section == "WCC" else ""
            yield_stage = "post_purification" if section == "WCC" else ""
            index += 2
        else:
            raise AssertionError(f"Unexpected paragraph {index + 1}: {text}")
        if sequence_paragraph["raw_text"] != sequence_paragraph["parse_text"]:
            raise AssertionError(f"Sequence whitespace for {section}__{clone_id}")
        sequence = sequence_paragraph["raw_text"]
        if not AA_RE.fullmatch(sequence):
            raise AssertionError(f"Invalid sequence for {section}__{clone_id}")
        records.append(
            {
                "sample_uid": f"{section}__{clone_id}",
                "provider_code": section,
                "source_sample_id": clone_id,
                "source_header_raw": source_header_raw,
                "sequence_raw": sequence,
                "sequence_length_aa": str(len(sequence)),
                "sequence_sha256": sha256_bytes(sequence.encode("ascii")),
                "reported_text": reported_text,
                "observation_semantics": semantics,
                "value_relation": relation,
                "point_estimate_mg": point_estimate,
                "group_anchor_mg": group_anchor,
                "lower_bound_mg": lower_bound,
                "culture_volume_l": culture_volume_l.rstrip("0").rstrip(".")
                if "." in culture_volume_l
                else culture_volume_l,
                "volume_evidence": volume_evidence,
                "medium": medium,
                "yield_stage": yield_stage,
            }
        )
    return records, paragraphs, sha256_bytes(document_xml)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def read_fasta(path: Path) -> list[tuple[str, str]]:
    entries: list[tuple[str, str]] = []
    header = ""
    sequence = ""
    for line in path.read_text(encoding="ascii").splitlines():
        if line.startswith(">"):
            if header:
                entries.append((header, sequence))
            header, sequence = line[1:], ""
        else:
            sequence += line
    if header:
        entries.append((header, sequence))
    return entries


def main() -> int:
    args = arguments()
    source_lexical = args.source.expanduser().absolute()
    artifact_dir_lexical = args.artifact_dir.expanduser().absolute()
    report_lexical = args.report.expanduser().absolute()
    validated_paths = validate_file_paths(
        project_root=PROJECT_ROOT,
        source_paths=[source_lexical, artifact_dir_lexical],
        target_paths=[report_lexical],
    )
    source, artifact_dir = validated_paths.source_paths
    report_path = validated_paths.target_paths[0]
    source_hash = sha256_file(source)
    if source_hash.lower() != args.expected_source_sha256.lower():
        raise AssertionError(
            f"Source SHA-256 mismatch: {source_hash} != {args.expected_source_sha256}"
        )
    if report_path.exists() and not args.overwrite:
        raise FileExistsError(f"Refusing to overwrite {report_path}")

    expected, paragraphs, document_xml_hash = source_records(
        source, args.document_title_culture_volume_l
    )
    counts = dict(Counter(record["provider_code"] for record in expected))
    if counts != EXPECTED_COUNTS or len(expected) != 47 or len(paragraphs) != 100:
        raise AssertionError(f"Unexpected source structure: {counts}, {len(paragraphs)} paragraphs")

    samples = read_csv(artifact_dir / "samples.csv")
    yields = read_csv(artifact_dir / "yield_observations.csv")
    wide = read_csv(artifact_dir / "nb_expression_records.csv")
    raw = read_csv(artifact_dir / "raw_transcription.csv")
    fasta = read_fasta(artifact_dir / "nb_expression_sequences.fasta")
    expected_raw = [
        {
            "source_nonempty_paragraph_index": str(index),
            "raw_text": paragraph["raw_text"],
            "parse_text": paragraph["parse_text"],
            "leading_or_trailing_whitespace": str(
                paragraph["raw_text"] != paragraph["parse_text"]
            ),
        }
        for index, paragraph in enumerate(paragraphs, start=1)
    ]
    if raw != expected_raw:
        raise AssertionError("raw_transcription.csv differs from literal source paragraphs")
    if not (len(samples) == len(yields) == len(wide) == len(fasta) == len(expected)):
        raise AssertionError("Output record counts differ")

    sample_fields = [
        "sample_uid",
        "provider_code",
        "source_sample_id",
        "source_header_raw",
        "sequence_raw",
        "sequence_length_aa",
        "sequence_sha256",
    ]
    yield_fields = [
        "sample_uid",
        "reported_text",
        "observation_semantics",
        "value_relation",
        "point_estimate_mg",
        "group_anchor_mg",
        "lower_bound_mg",
    ]
    context_fields = ["culture_volume_l", "volume_evidence", "medium", "yield_stage"]
    for position, (source_row, sample, observation, combined, fasta_record) in enumerate(
        zip(expected, samples, yields, wide, fasta, strict=True), start=1
    ):
        for field in sample_fields:
            if sample[field] != source_row[field] or combined[field] != source_row[field]:
                raise AssertionError(f"Row {position} {field} mismatch")
        for field in yield_fields:
            if observation[field] != source_row[field] or combined[field] != source_row[field]:
                raise AssertionError(f"Row {position} {field} mismatch")
        for field in context_fields:
            if combined[field] != source_row[field]:
                raise AssertionError(f"Row {position} {field} mismatch")
        expected_header = f"{source_row['sample_uid']}|sha256={source_row['sequence_sha256']}"
        if fasta_record != (expected_header, source_row["sequence_raw"]):
            raise AssertionError(f"Row {position} FASTA mismatch")

    report = {
        "status": "pass",
        "verification_implementation": "independent_stdlib_parser_v2",
        "production_parser_imported": False,
        "source_file": source.name,
        "source_sha256": source_hash,
        "word_document_xml_sha256": document_xml_hash,
        "source_nonempty_paragraph_count": len(paragraphs),
        "paragraphs_with_outer_whitespace": [
            index
            for index, paragraph in enumerate(paragraphs, start=1)
            if paragraph["raw_text"] != paragraph["parse_text"]
        ],
        "record_count": len(expected),
        "record_counts_by_source": counts,
        "unique_sample_uid_count": len({row["sample_uid"] for row in expected}),
        "unique_sequence_count": len({row["sequence_sha256"] for row in expected}),
        "sequence_mismatch_count": 0,
        "sequence_hash_mismatch_count": 0,
        "sequence_length_mismatch_count": 0,
        "yield_field_mismatch_count": 0,
        "context_field_mismatch_count": 0,
        "raw_transcription_mismatch_count": 0,
        "fasta_mismatch_count": 0,
        "artifacts_checked": {
            name: {
                "sha256": sha256_file(artifact_dir / name),
                "size_bytes": (artifact_dir / name).stat().st_size,
            }
            for name in (
                "samples.csv",
                "yield_observations.csv",
                "nb_expression_records.csv",
                "raw_transcription.csv",
                "nb_expression_sequences.fasta",
            )
        },
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", newline="\n", dir=report_path.parent, delete=False
    ) as handle:
        handle.write(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n")
        staged_report = Path(handle.name)
    try:
        os.replace(staged_report, report_path)
    except Exception:
        staged_report.unlink(missing_ok=True)
        raise
    print(json.dumps(report, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
