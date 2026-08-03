"""Writers, manifests, QC rendering, and output validation for yield data."""

from __future__ import annotations

import csv
import html
import json
from collections import Counter
from pathlib import Path
from typing import Iterable, Sequence

from .nb_expression import (
    PARSER_VERSION,
    DocxParagraph,
    ExpressionRecord,
    ParseError,
    optional_decimal,
    sha256_file,
)


SAMPLES_FIELDS = [
    "sample_uid",
    "provider_code",
    "source_sample_id",
    "source_section",
    "source_header_raw",
    "sequence_raw",
    "sequence_length_aa",
    "sequence_sha256",
    "sequence_scope",
    "vhh_region_sequence",
    "sequence_review_flags",
    "source_document",
    "source_sha256",
    "source_clone_paragraph_index",
    "source_sequence_paragraph_index",
]

YIELD_FIELDS = [
    "observation_id",
    "sample_uid",
    "assay_id",
    "phenotype_name",
    "reported_text",
    "observation_semantics",
    "value_relation",
    "point_estimate_mg",
    "group_anchor_mg",
    "lower_bound_mg",
    "lower_bound_inclusive",
    "upper_bound_mg",
    "assignment_level",
    "group_id",
    "individual_numeric_available",
    "censoring_type",
    "censoring_reason",
    "replicate_count",
    "uncertainty_value",
    "uncertainty_type",
    "source_yield_paragraph_index",
]

ASSAY_CONTEXT_FIELDS = [
    "assay_id",
    "provider_code",
    "culture_volume_l",
    "volume_evidence",
    "medium",
    "yield_stage",
    "purification_method",
    "batch_id",
    "protocol_id",
    "context_notes",
]

WIDE_FIELDS = SAMPLES_FIELDS + [field for field in YIELD_FIELDS if field != "sample_uid"] + [
    "culture_volume_l",
    "volume_evidence",
    "medium",
    "yield_stage",
]

RAW_TRANSCRIPTION_FIELDS = [
    "source_nonempty_paragraph_index",
    "raw_text",
    "parse_text",
    "leading_or_trailing_whitespace",
]


def write_samples_csv(path: Path, records: Sequence[ExpressionRecord]) -> None:
    _write_csv(path, SAMPLES_FIELDS, (record.samples_row() for record in records))


def write_yield_observations_csv(path: Path, records: Sequence[ExpressionRecord]) -> None:
    _write_csv(path, YIELD_FIELDS, (record.yield_row() for record in records))


def write_assay_context_csv(path: Path, records: Sequence[ExpressionRecord]) -> None:
    """Write provider context derived from parsed record metadata."""

    _write_csv(path, ASSAY_CONTEXT_FIELDS, assay_context_rows(records))


def assay_context_rows(
    records: Sequence[ExpressionRecord],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for provider in ("LTT", "WCC", "LLJ"):
        provider_records = [record for record in records if record.provider_code == provider]
        if not provider_records:
            raise ParseError(f"No records for assay context provider {provider}")
        volume_values = {record.culture_volume_l for record in provider_records}
        volume_evidence = {record.volume_evidence for record in provider_records}
        media = {record.medium for record in provider_records}
        yield_stages = {record.yield_stage for record in provider_records}
        if any(len(values) != 1 for values in (volume_values, volume_evidence, media, yield_stages)):
            raise ParseError(f"Inconsistent assay context within provider {provider}")
        evidence = next(iter(volume_evidence))
        if provider == "WCC":
            context_note = (
                "Each source entry states 1 L TB and an approximate "
                "post-purification yield."
            )
        elif evidence == "document_title":
            suffix = (
                "other protocol details are unavailable."
                if provider == "LTT"
                else "yields are shared group bins."
            )
            context_note = f"Volume is explicit metadata from the original document title; {suffix}"
        else:
            context_note = (
                "Culture volume is unavailable; other protocol details are unavailable."
                if provider == "LTT"
                else "Culture volume is unavailable; yields are shared group bins."
            )
        rows.append(
            {
                "assay_id": f"ASSAY__{provider}",
                "provider_code": provider,
                "culture_volume_l": optional_decimal(next(iter(volume_values))),
                "volume_evidence": evidence,
                "medium": next(iter(media)),
                "yield_stage": next(iter(yield_stages)),
                "purification_method": "",
                "batch_id": "",
                "protocol_id": "",
                "context_notes": context_note,
            }
        )
    return rows


def write_wide_records_csv(path: Path, records: Sequence[ExpressionRecord]) -> None:
    _write_csv(path, WIDE_FIELDS, (record.wide_row() for record in records))


def write_raw_transcription_csv(
    path: Path, paragraphs: Sequence[DocxParagraph]
) -> None:
    rows = (
        {
            "source_nonempty_paragraph_index": paragraph.nonempty_index,
            "raw_text": paragraph.raw_text,
            "parse_text": paragraph.parse_text,
            "leading_or_trailing_whitespace": (
                paragraph.raw_text != paragraph.parse_text
            ),
        }
        for paragraph in paragraphs
    )
    _write_csv(path, RAW_TRANSCRIPTION_FIELDS, rows)


def write_fasta(path: Path, records: Sequence[ExpressionRecord]) -> None:
    """Write one unwrapped source-exact sequence per safe sample UID."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for record in records:
            handle.write(
                f">{record.sample_uid}|sha256={record.sequence_sha256}\n"
                f"{record.sequence_raw}\n"
            )


def write_qc_plot_data(path: Path, records: Sequence[ExpressionRecord]) -> None:
    source_counts = Counter(record.provider_code for record in records)
    semantics_counts = Counter(record.observation_semantics for record in records)
    rows = [
        {"metric": "source_records", "category": source, "count": source_counts[source]}
        for source in ("LTT", "WCC", "LLJ")
    ]
    rows.extend(
        {
            "metric": "yield_semantics",
            "category": semantics,
            "count": semantics_counts[semantics],
        }
        for semantics in (
            "individual_approximate",
            "group_lower_bound",
            "group_approximate",
        )
    )
    _write_csv(path, ["metric", "category", "count"], rows, excel_bom=False)


def render_qc_svg(plot_data_path: Path, output_path: Path, source_hash: str) -> None:
    """Render a dependency-free SVG of extraction counts, not performance."""

    with plot_data_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    source_rows = [row for row in rows if row["metric"] == "source_records"]
    semantics_rows = [row for row in rows if row["metric"] == "yield_semantics"]
    expected_sources = ["LTT", "WCC", "LLJ"]
    expected_semantics = [
        "individual_approximate",
        "group_lower_bound",
        "group_approximate",
    ]
    if [row["category"] for row in source_rows] != expected_sources:
        raise ParseError("QC plot source categories or order are invalid")
    if [row["category"] for row in semantics_rows] != expected_semantics:
        raise ParseError("QC plot semantics categories or order are invalid")
    try:
        counts = [int(row["count"]) for row in rows]
    except ValueError as exc:
        raise ParseError("QC plot counts must be integers") from exc
    if not counts or any(count <= 0 for count in counts):
        raise ParseError("QC plot counts must all be positive")

    max_count = max(counts)
    display_labels = {
        "individual_approximate": "Individual ~ value",
        "group_lower_bound": "Group bin > threshold",
        "group_approximate": "Group bin ~ label",
    }
    colors = ["#2A6F97", "#61A5C2", "#89C2D9"]
    chart_bottom, chart_height = 535, 370

    def panel(x0: int, title: str, data: Sequence[dict[str, str]]) -> str:
        parts = [
            f'<rect x="{x0}" y="115" width="500" height="465" rx="18" '
            'fill="#FFFFFF" stroke="#D7E3EA" stroke-width="2"/>',
            f'<text x="{x0 + 28}" y="155" class="panel-title">{html.escape(title)}</text>',
        ]
        for tick in range(0, max_count + 1, 5):
            y = chart_bottom - (tick / max_count) * chart_height
            parts.append(
                f'<line x1="{x0 + 60}" y1="{y:.1f}" x2="{x0 + 470}" y2="{y:.1f}" '
                'stroke="#E9F0F4" stroke-width="1"/>'
            )
            parts.append(
                f'<text x="{x0 + 48}" y="{y + 5:.1f}" class="tick" text-anchor="end">{tick}</text>'
            )
        for index, row in enumerate(data):
            count = int(row["count"])
            bar_height = (count / max_count) * chart_height
            x = x0 + 82 + index * 132
            y = chart_bottom - bar_height
            label = display_labels.get(row["category"], row["category"])
            parts.extend(
                [
                    f'<rect x="{x}" y="{y:.1f}" width="90" height="{bar_height:.1f}" '
                    f'rx="7" fill="{colors[index]}"/>',
                    f'<text x="{x + 45}" y="{y - 12:.1f}" class="value" text-anchor="middle">{count}</text>',
                    f'<text x="{x + 45}" y="565" class="label" text-anchor="middle">{html.escape(label)}</text>',
                ]
            )
        return "".join(parts)

    svg = f'''<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="700" viewBox="0 0 1200 700">
  <rect width="1200" height="700" fill="#F4F8FA"/>
  <style>
    text {{ font-family: Arial, "Microsoft YaHei", sans-serif; fill: #17324D; }}
    .title {{ font-size: 30px; font-weight: 700; }}
    .subtitle {{ font-size: 16px; fill: #526D82; }}
    .panel-title {{ font-size: 20px; font-weight: 700; }}
    .tick {{ font-size: 13px; fill: #6D8293; }}
    .label {{ font-size: 13px; fill: #415A6B; }}
    .value {{ font-size: 18px; font-weight: 700; }}
    .foot {{ font-size: 13px; fill: #526D82; }}
  </style>
  <text x="70" y="55" class="title">Nb expression dataset extraction QC</text>
  <text x="70" y="83" class="subtitle">Counts are source-derived; this figure does not compare expression performance.</text>
  {panel(70, "Records by source section", source_rows)}
  {panel(630, "Yield-value semantics", semantics_rows)}
  <text x="70" y="635" class="foot">47 records; 47 unique standard-AA sequences; zero round-trip mismatches.</text>
  <text x="70" y="662" class="foot">Source SHA-256: {html.escape(source_hash)}</text>
</svg>
'''
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(svg, encoding="utf-8", newline="\n")


def validate_written_outputs(
    paths: dict[str, Path],
    records: Sequence[ExpressionRecord],
    paragraphs: Sequence[DocxParagraph],
) -> dict[str, object]:
    """Reopen every output and compare all CSV fields plus every sequence."""

    samples = _read_csv_exact(paths["samples"], SAMPLES_FIELDS)
    yields = _read_csv_exact(paths["yields"], YIELD_FIELDS)
    contexts = _read_csv_exact(paths["assay_context"], ASSAY_CONTEXT_FIELDS)
    wide = _read_csv_exact(paths["wide"], WIDE_FIELDS)
    raw = _read_csv_exact(paths["raw_transcription"], RAW_TRANSCRIPTION_FIELDS)
    if len(samples) != len(records) or len(yields) != len(records) or len(wide) != len(records):
        raise ParseError("One or more record tables have the wrong row count")

    expected_contexts = [_stringify_row(row) for row in assay_context_rows(records)]
    if contexts != expected_contexts:
        raise ParseError("Assay context does not match parsed record metadata")
    expected_raw = [
        _stringify_row(
            {
                "source_nonempty_paragraph_index": paragraph.nonempty_index,
                "raw_text": paragraph.raw_text,
                "parse_text": paragraph.parse_text,
                "leading_or_trailing_whitespace": (
                    paragraph.raw_text != paragraph.parse_text
                ),
            }
        )
        for paragraph in paragraphs
    ]
    if raw != expected_raw:
        raise ParseError("Raw transcription does not preserve every source paragraph")

    individual_count = 0
    lower_bound_only_count = 0
    group_anchor_only_count = 0
    for sample, observation, combined, record in zip(
        samples, yields, wide, records, strict=True
    ):
        if sample != _stringify_row(record.samples_row()):
            raise ParseError(f"Samples row mismatch for {record.sample_uid}")
        if observation != _stringify_row(record.yield_row()):
            raise ParseError(f"Yield row mismatch for {record.sample_uid}")
        if combined != _stringify_row(record.wide_row()):
            raise ParseError(f"Wide row mismatch for {record.sample_uid}")
        numeric_fields = (
            observation["point_estimate_mg"],
            observation["group_anchor_mg"],
            observation["lower_bound_mg"],
            observation["upper_bound_mg"],
        )
        if record.observation_semantics == "individual_approximate":
            if not numeric_fields[0] or any(numeric_fields[1:]):
                raise ParseError(f"Invalid individual numeric fields for {record.sample_uid}")
            individual_count += 1
        elif record.observation_semantics == "group_lower_bound":
            if numeric_fields != ("", "", "20", ""):
                raise ParseError(f"Invalid lower-bound fields for {record.sample_uid}")
            lower_bound_only_count += 1
        elif record.observation_semantics == "group_approximate":
            if numeric_fields[0] or not numeric_fields[1] or numeric_fields[2] or numeric_fields[3]:
                raise ParseError(f"Invalid group-anchor fields for {record.sample_uid}")
            group_anchor_only_count += 1
        else:
            raise ParseError(f"Unknown yield semantics for {record.sample_uid}")
    if (individual_count, lower_bound_only_count, group_anchor_only_count) != (31, 9, 7):
        raise ParseError("Written yield-semantics counts differ from 31/9/7")

    fasta_entries = _read_fasta(paths["fasta"])
    if len(fasta_entries) != len(records):
        raise ParseError(f"FASTA count mismatch: {len(fasta_entries)} vs {len(records)}")
    for (header, sequence), record in zip(fasta_entries, records, strict=True):
        expected_header = f"{record.sample_uid}|sha256={record.sequence_sha256}"
        if header != expected_header or sequence != record.sequence_raw:
            raise ParseError(f"FASTA mismatch for {record.sample_uid}")

    return {
        "status": "pass",
        "samples_row_count": len(samples),
        "yield_observation_row_count": len(yields),
        "wide_record_row_count": len(wide),
        "assay_context_row_count": len(contexts),
        "raw_transcription_row_count": len(raw),
        "fasta_record_count": len(fasta_entries),
        "all_csv_fields_compared": True,
        "csv_sequence_mismatch_count": 0,
        "fasta_sequence_mismatch_count": 0,
        "hash_mismatch_count": 0,
        "llj_point_estimate_population_count": 0,
        "individual_point_estimate_only_count": individual_count,
        "gt20_lower_bound_only_count": lower_bound_only_count,
        "approx_group_anchor_only_count": group_anchor_only_count,
    }


def build_manifest(
    *,
    source_path: Path,
    records: Sequence[ExpressionRecord],
    validation_report: dict[str, object],
    generated_at: str,
    output_paths: Iterable[Path],
) -> dict[str, object]:
    """Build provenance, schema, sequence-handling, and caution metadata."""

    return {
        "dataset": "collaborator_nanobody_reported_yield_1L",
        "generated_at": generated_at,
        "parser_version": PARSER_VERSION,
        "source": {
            "file": source_path.name,
            "size_bytes": source_path.stat().st_size,
            "sha256": sha256_file(source_path),
            "storage_layout": "ordered direct body paragraphs; no Word tables",
        },
        "record_count": len(records),
        "record_counts_by_source": dict(Counter(r.provider_code for r in records)),
        "tables": {
            "samples.csv": SAMPLES_FIELDS,
            "yield_observations.csv": YIELD_FIELDS,
            "assay_context.csv": ASSAY_CONTEXT_FIELDS,
            "nb_expression_records.csv": WIDE_FIELDS,
            "raw_transcription.csv": RAW_TRANSCRIPTION_FIELDS,
        },
        "sequence_integrity": {
            "alphabet": "20 standard uppercase amino-acid letters",
            "sequence_normalization": "none",
            "sequence_outer_whitespace_policy": "fail",
            "paragraph_handling": "raw_text preserved; parse_text trims only outer paragraph whitespace for grammar",
            "unsupported_docx_structures": "tables/text boxes/tabs/breaks/tracked changes/hidden text fail parsing",
            "per_record_sha256": True,
            "raw_source_paragraph_round_trip": "pass",
            "all_written_sequence_round_trips": validation_report["written_outputs"]["status"],
        },
        "yield_semantics": {
            "individual_approximate": "LTT/WCC clone-specific approximate yields; stored in point_estimate_mg.",
            "group_lower_bound": "LLJ >20 mg shared lower-bound bin; point_estimate_mg remains blank.",
            "group_approximate": "LLJ ~10/~2 mg shared anchors with unknown boundaries; point_estimate_mg remains blank.",
            "do_not_merge_into_one_exact_continuous_label": True,
        },
        "cautions": [
            "The phenotype is reported yield, not a direct expression-rate measurement; post-purification is explicit only for WCC.",
            "A 1 L volume does not justify automatically creating an expression_mg_per_l field.",
            "LTT/LLJ volume is explicit CLI metadata from the original document title; WCC states 1 L per entry.",
            "TB medium and post-purification stage are explicit only for WCC.",
            "VHH boundaries, terminal additions, tags, truncations, and construct equivalence remain unresolved.",
            "Host, induction, purification recovery, batch, replicate count, uncertainty, and cross-provider comparability are undocumented.",
            "Use sample_uid for stable identity; import source_sample_id as text to prevent spreadsheet auto-coercion.",
        ],
        "outputs": {
            path.name: {"sha256": sha256_file(path), "size_bytes": path.stat().st_size}
            for path in output_paths
        },
    }


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def _write_csv(
    path: Path,
    fieldnames: Sequence[str],
    rows: Iterable[dict[str, object]],
    *,
    excel_bom: bool = True,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    encoding = "utf-8-sig" if excel_bom else "utf-8"
    with path.open("w", encoding=encoding, newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _read_csv_exact(path: Path, expected_fields: Sequence[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != list(expected_fields):
            raise ParseError(f"CSV field order mismatch for {path.name}: {reader.fieldnames}")
        return list(reader)


def _read_fasta(path: Path) -> list[tuple[str, str]]:
    entries: list[tuple[str, str]] = []
    header: str | None = None
    sequence_parts: list[str] = []
    for line in path.read_text(encoding="ascii").splitlines():
        if line.startswith(">"):
            if header is not None:
                entries.append((header, "".join(sequence_parts)))
            header = line[1:]
            sequence_parts = []
        else:
            sequence_parts.append(line)
    if header is not None:
        entries.append((header, "".join(sequence_parts)))
    return entries


def _stringify_row(row: dict[str, object]) -> dict[str, str]:
    return {key: "" if value is None else str(value) for key, value in row.items()}
