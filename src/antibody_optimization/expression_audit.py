"""Build a conservative comparability audit for collaborator-reported VHH yields.

The audit preserves the serialized source semantics and assigns allowed-use gates.  It
does not normalize yields, infer missing protocol metadata, train a model, or decide
that a provisional antibody-numbering span is a confirmed mature VHH construct.
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Iterable, Mapping, Sequence


AUDIT_VERSION = "1.0.0"
EXPECTED_SOURCE_COUNTS = {"LTT": 23, "WCC": 8, "LLJ": 16}
EXPECTED_SEMANTICS_COUNTS = {
    "individual_approximate": 31,
    "group_lower_bound": 9,
    "group_approximate": 7,
}
PROVIDER_ORDER = ("LTT", "WCC", "LLJ")

ASSAY_CONTEXT_FIELDS = (
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
)
ASSAY_METADATA_FIELDS = (
    "assay_id",
    "provider_code",
    "field_name",
    "field_value",
    "unit",
    "evidence_status",
    "evidence_source",
    "source_locator",
    "review_note",
    "reviewed_at",
)
METADATA_FIELD_ORDER = (
    "culture_volume_l",
    "medium",
    "host_species",
    "host_strain",
    "vector_backbone",
    "promoter",
    "signal_peptide",
    "construct_boundary_definition",
    "terminal_tags_or_linkers",
    "expression_compartment",
    "induction_agent",
    "induction_concentration",
    "induction_temperature_c",
    "induction_duration_h",
    "harvest_definition",
    "yield_stage",
    "purification_method",
    "purification_recovery",
    "quantification_method",
    "batch_id",
    "protocol_id",
    "replicate_design",
    "replicate_count",
    "uncertainty_definition",
    "system_equivalence_claim",
    "cross_provider_protocol_equivalence",
)
EVIDENCE_STATUSES = {
    "source_explicit",
    "user_provided",
    "collaborator_confirmed",
    "derived_exact",
    "unknown_not_reported",
    "conflicting",
    "not_applicable",
}
GATE_STATUSES = {"pass", "conditional", "blocked", "pending", "not_applicable"}

SAMPLE_REVIEW_FIELDS = (
    "sample_uid",
    "observation_id",
    "assay_id",
    "provider_code",
    "source_sample_id",
    "source_sha256",
    "sequence_sha256",
    "sequence_scope",
    "sequence_scope_status",
    "numbering_status",
    "provisional_numbered_span_sha256",
    "sequence_review_flags",
    "construct_comparability_status",
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
    "individual_numeric_available",
    "censoring_type",
    "protocol_completeness_status",
    "within_assay_numeric_use_status",
    "within_assay_ordinal_use_status",
    "cross_assay_pooling_status",
    "nb252_transfer_status",
    "highest_allowed_use",
    "blocking_fields",
    "decision_reasons",
    "review_version",
)
VIEW_FIELDS = (
    "sample_uid",
    "provider_code",
    "source_sample_id",
    "assay_id",
    "observation_id",
    "phenotype_name",
    "reported_text",
    "observation_semantics",
    "value_relation",
    "point_estimate_mg",
    "group_anchor_mg",
    "lower_bound_mg",
    "censoring_type",
    "culture_volume_l",
    "culture_volume_evidence_status",
    "medium",
    "medium_evidence_status",
    "yield_stage",
    "yield_stage_evidence_status",
    "system_equivalence_claim",
    "system_equivalence_claim_evidence_status",
    "cross_provider_protocol_equivalence_status",
    "sequence_scope_status",
    "numbering_status",
    "provisional_numbered_span_sha256",
    "sequence_review_flags",
    "within_assay_numeric_use_status",
    "within_assay_ordinal_use_status",
    "cross_assay_pooling_status",
    "nb252_transfer_status",
    "highest_allowed_use",
    "blocking_fields",
)


class ExpressionAuditError(ValueError):
    """Raised when an input or derived audit violates the frozen data contract."""


@dataclass(frozen=True)
class AuditInputs:
    records: tuple[dict[str, str], ...]
    assay_contexts: tuple[dict[str, str], ...]
    source_manifest: dict[str, object]
    sequence_samples: Mapping[str, Mapping[str, str]]
    input_sha256: Mapping[str, str]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_and_validate_inputs(
    records_path: Path,
    assay_context_path: Path,
    manifest_path: Path,
    sequence_audit_summary_path: Path | None = None,
) -> AuditInputs:
    """Read frozen artifacts, validate hashes/schema/semantics, and join no data yet."""

    manifest = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
    if manifest.get("dataset") != "collaborator_nanobody_reported_yield_1L":
        raise ExpressionAuditError("Unexpected source manifest dataset")
    tables = manifest.get("tables")
    if not isinstance(tables, dict):
        raise ExpressionAuditError("Source manifest is missing table schemas")
    record_fields = tables.get("nb_expression_records.csv")
    if not isinstance(record_fields, list):
        raise ExpressionAuditError("Source manifest lacks the wide-record schema")
    records = _read_csv_exact(records_path, record_fields)
    contexts = _read_csv_exact(assay_context_path, ASSAY_CONTEXT_FIELDS)

    _validate_manifest_file_hash(manifest, "nb_expression_records.csv", records_path)
    _validate_manifest_file_hash(manifest, "assay_context.csv", assay_context_path)
    _validate_records(records, manifest)
    _validate_contexts(contexts, records)

    sequence_samples = _read_sequence_summary(sequence_audit_summary_path, records)
    hashes = {
        "nb_expression_records.csv": sha256_file(records_path),
        "assay_context.csv": sha256_file(assay_context_path),
        "manifest.json": sha256_file(manifest_path),
    }
    if sequence_audit_summary_path is not None:
        hashes["sequence_audit_summary"] = sha256_file(sequence_audit_summary_path)
    return AuditInputs(
        tuple(records), tuple(contexts), manifest, sequence_samples, hashes
    )


def build_assay_metadata_rows(
    contexts: Sequence[Mapping[str, str]],
    *,
    generated_at: str,
) -> list[dict[str, str]]:
    """Expand three source-level contexts into evidence-bearing metadata rows."""

    by_provider = {row["provider_code"]: row for row in contexts}
    rows: list[dict[str, str]] = []
    for provider in PROVIDER_ORDER:
        context = by_provider[provider]
        for field_name in METADATA_FIELD_ORDER:
            value = ""
            unit = ""
            status = "unknown_not_reported"
            source = ""
            locator = ""
            note = "Not reported in the current source artifacts."
            if field_name == "culture_volume_l":
                value, unit, status = context[field_name], "L", "source_explicit"
                source, locator = "source_document", context["volume_evidence"]
                note = "Volume evidence is retained at its original document scope."
            elif field_name in {"medium", "yield_stage"} and context[field_name]:
                value, status = context[field_name], "source_explicit"
                source, locator = "source_document", "entry_text"
                note = "Explicit only for WCC entries; not propagated to other sources."
            elif field_name in {"purification_method", "batch_id", "protocol_id"}:
                value = context[field_name]
                if value:
                    status, source, locator = (
                        "source_explicit",
                        "source_document",
                        "assay_context.csv",
                    )
            elif field_name == "system_equivalence_claim":
                value, status = "same_system", "user_provided"
                source, locator = "user_statement", "project conversation"
                note = (
                    "High-level claim only; it does not establish field-level protocol "
                    "or construct equivalence."
                )
            elif field_name == "cross_provider_protocol_equivalence":
                note = (
                    "Unknown until host, construct, induction, yield stage, purification, "
                    "quantification, batch, and replicate metadata are reconciled."
                )
            rows.append(
                {
                    "assay_id": context["assay_id"],
                    "provider_code": provider,
                    "field_name": field_name,
                    "field_value": value,
                    "unit": unit,
                    "evidence_status": status,
                    "evidence_source": source,
                    "source_locator": locator,
                    "review_note": note,
                    "reviewed_at": generated_at,
                }
            )
    return rows


def build_sample_review_rows(inputs: AuditInputs) -> list[dict[str, str]]:
    """Assign conservative per-sample use gates without changing source values."""

    rows: list[dict[str, str]] = []
    for record in inputs.records:
        sequence = inputs.sequence_samples.get(record["sample_uid"], {})
        scope_status = sequence.get("sequence_scope_status") or "unknown_not_reported"
        if scope_status == "unknown":
            scope_status = "unknown_not_reported"
        numbering_status = sequence.get("numbering_status") or "pending"
        span_hash = sequence.get("provisional_numbered_span_sha256") or ""
        common_blockers = {
            "batch_id",
            "construct_boundary_definition",
            "cross_provider_protocol_equivalence",
            "host_species",
            "protocol_id",
            "replicate_count",
            "uncertainty_definition",
        }
        if numbering_status != "pass":
            common_blockers.add("numbering_status")
        if scope_status not in {"confirmed", "confirmed_mature_vhh"}:
            common_blockers.add("sequence_scope")
        if record["provider_code"] in {"LTT", "LLJ"}:
            common_blockers.update({"medium", "yield_stage"})

        if record["provider_code"] in {"LTT", "WCC"}:
            numeric_status, ordinal_status = "conditional", "not_applicable"
            highest = "within_assay_numeric_exploratory"
            reasons = (
                "individual_approximate_numeric_but_protocol_incomplete;"
                "cross_source_pooling_not_authorized;"
                "nb252_transfer_not_authorized"
            )
        else:
            numeric_status, ordinal_status = "blocked", "conditional"
            highest = "within_assay_ordinal_exploratory"
            reasons = (
                "group_level_label_not_continuous_target;"
                "ordinal_or_censored_use_only;"
                "cross_source_pooling_not_authorized;"
                "nb252_transfer_not_authorized"
            )
        row = {field: record.get(field, "") for field in SAMPLE_REVIEW_FIELDS}
        row.update(
            {
                "sequence_scope_status": scope_status,
                "numbering_status": numbering_status,
                "provisional_numbered_span_sha256": span_hash,
                "construct_comparability_status": "pending",
                "protocol_completeness_status": "blocked",
                "within_assay_numeric_use_status": numeric_status,
                "within_assay_ordinal_use_status": ordinal_status,
                "cross_assay_pooling_status": "blocked",
                "nb252_transfer_status": "blocked",
                "highest_allowed_use": highest,
                "blocking_fields": ";".join(sorted(common_blockers)),
                "decision_reasons": reasons,
                "review_version": AUDIT_VERSION,
            }
        )
        rows.append(row)
    return rows


def build_comparability_view(
    sample_rows: Sequence[Mapping[str, str]],
    metadata_rows: Sequence[Mapping[str, str]],
) -> list[dict[str, str]]:
    """Materialize selected source evidence beside each sample decision."""

    metadata = {
        (row["assay_id"], row["field_name"]): row for row in metadata_rows
    }
    rows: list[dict[str, str]] = []
    for sample in sample_rows:
        assay_id = sample["assay_id"]
        volume = metadata[(assay_id, "culture_volume_l")]
        medium = metadata[(assay_id, "medium")]
        stage = metadata[(assay_id, "yield_stage")]
        system = metadata[(assay_id, "system_equivalence_claim")]
        equivalence = metadata[(assay_id, "cross_provider_protocol_equivalence")]
        combined = dict(sample)
        combined.update(
            {
                "culture_volume_l": volume["field_value"],
                "culture_volume_evidence_status": volume["evidence_status"],
                "medium": medium["field_value"],
                "medium_evidence_status": medium["evidence_status"],
                "yield_stage": stage["field_value"],
                "yield_stage_evidence_status": stage["evidence_status"],
                "system_equivalence_claim": system["field_value"],
                "system_equivalence_claim_evidence_status": system["evidence_status"],
                "cross_provider_protocol_equivalence_status": equivalence[
                    "evidence_status"
                ],
            }
        )
        rows.append({field: combined.get(field, "") for field in VIEW_FIELDS})
    return rows


def build_allowed_use_manifest(
    *,
    inputs: AuditInputs,
    metadata_rows: Sequence[Mapping[str, str]],
    sample_rows: Sequence[Mapping[str, str]],
    view_rows: Sequence[Mapping[str, str]],
    generated_at: str,
    output_hashes: Mapping[str, str],
) -> dict[str, object]:
    """Describe evidence, allowed uses, blocked uses, and the expression sub-gate."""

    return {
        "audit": "nb_expression_comparability",
        "audit_version": AUDIT_VERSION,
        "generated_at": generated_at,
        "inputs": dict(inputs.input_sha256),
        "sequence_audit_summary": {
            "status": "provided" if "sequence_audit_summary" in inputs.input_sha256 else "not_provided",
            "sample_coverage_count": len(inputs.sequence_samples),
            "numbering_status_counts": dict(
                Counter(row["numbering_status"] for row in sample_rows)
            ),
            "sequence_scope_status_counts": dict(
                Counter(row["sequence_scope_status"] for row in sample_rows)
            ),
            "provisional_numbering_is_not_construct_confirmation": True,
        },
        "counts": {
            "samples": len(sample_rows),
            "assay_metadata_rows": len(metadata_rows),
            "view_rows": len(view_rows),
            "by_provider": dict(Counter(row["provider_code"] for row in sample_rows)),
            "by_observation_semantics": dict(
                Counter(row["observation_semantics"] for row in sample_rows)
            ),
            "by_highest_allowed_use": dict(
                Counter(row["highest_allowed_use"] for row in sample_rows)
            ),
        },
        "evidence_policy": {
            "same_system_claim": "user_provided_summary_only",
            "cross_provider_protocol_equivalence": "unknown_not_reported",
            "reported_yield_is_expression_rate": False,
            "one_liter_implies_mg_per_l_expression": False,
        },
        "allowed_uses": {
            "LTT_WCC": "conditional within-assay numeric exploration only",
            "LLJ": "conditional within-assay ordinal/censored exploration only",
            "sequence": "descriptive numbering and leakage audit only",
        },
        "blocked_uses": [
            "cross_assay_pooling",
            "pooled_continuous_yield_model",
            "nb252_candidate_transfer_or_ranking",
            "conversion_to_expression_mg_per_l",
            "conversion_of_LLJ_group_bins_to_exact_individual_values",
        ],
        "gates": {
            "expression_audit_gate": "pass",
            "cross_assay_pooling_gate": "blocked",
            "nb252_transfer_gate": "blocked",
            "stage_1_baseline_gate": (
                "pending_structure_baseline"
                if "sequence_audit_summary" in inputs.input_sha256
                else "pending_sequence_and_structure_baselines"
            ),
        },
        "schemas": {
            "assay_metadata_review.csv": list(ASSAY_METADATA_FIELDS),
            "sample_comparability_review.csv": list(SAMPLE_REVIEW_FIELDS),
            "expression_comparability_view.csv": list(VIEW_FIELDS),
        },
        "outputs": dict(output_hashes),
    }


def write_csv(
    path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, object]]
) -> None:
    """Write an Excel-readable UTF-8 BOM CSV with a fixed schema."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def validate_written_audit(paths: Mapping[str, Path]) -> dict[str, object]:
    metadata = _read_csv_exact(paths["assay_metadata"], ASSAY_METADATA_FIELDS)
    samples = _read_csv_exact(paths["sample_review"], SAMPLE_REVIEW_FIELDS)
    view = _read_csv_exact(paths["view"], VIEW_FIELDS)
    for name in ("assay_metadata", "sample_review", "view"):
        if not paths[name].read_bytes().startswith(b"\xef\xbb\xbf"):
            raise ExpressionAuditError(f"CSV lacks UTF-8 BOM: {paths[name].name}")
    if len(metadata) != len(PROVIDER_ORDER) * len(METADATA_FIELD_ORDER):
        raise ExpressionAuditError("Assay metadata row count is not complete")
    if len(samples) != 47 or len(view) != 47:
        raise ExpressionAuditError("Sample audit and materialized view must each have 47 rows")
    if Counter(row["provider_code"] for row in samples) != Counter(EXPECTED_SOURCE_COUNTS):
        raise ExpressionAuditError("Written provider counts differ from 23/8/16")
    if Counter(row["observation_semantics"] for row in samples) != Counter(
        EXPECTED_SEMANTICS_COUNTS
    ):
        raise ExpressionAuditError("Written semantics counts differ from 31/9/7")
    if any(row["cross_assay_pooling_status"] != "blocked" for row in samples):
        raise ExpressionAuditError("Cross-assay pooling must remain blocked")
    if any(row["nb252_transfer_status"] != "blocked" for row in samples):
        raise ExpressionAuditError("Nb252 transfer must remain blocked")
    if any(row["evidence_status"] not in EVIDENCE_STATUSES for row in metadata):
        raise ExpressionAuditError("Unknown evidence status in assay metadata")
    for row in samples:
        for field in (
            "construct_comparability_status",
            "protocol_completeness_status",
            "within_assay_numeric_use_status",
            "within_assay_ordinal_use_status",
            "cross_assay_pooling_status",
            "nb252_transfer_status",
        ):
            if row[field] not in GATE_STATUSES:
                raise ExpressionAuditError(f"Unknown gate status in {field}")
    if [row["sample_uid"] for row in samples] != [row["sample_uid"] for row in view]:
        raise ExpressionAuditError("Materialized view changed sample order or identity")
    return {
        "status": "pass",
        "assay_metadata_row_count": len(metadata),
        "sample_review_row_count": len(samples),
        "view_row_count": len(view),
        "provider_counts": dict(Counter(row["provider_code"] for row in samples)),
        "semantics_counts": dict(
            Counter(row["observation_semantics"] for row in samples)
        ),
        "cross_assay_pooling_blocked_count": 47,
        "nb252_transfer_blocked_count": 47,
    }


def _read_csv_exact(path: Path, fields: Sequence[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != list(fields):
            raise ExpressionAuditError(
                f"CSV field order mismatch for {path.name}: {reader.fieldnames}"
            )
        return list(reader)


def _validate_manifest_file_hash(
    manifest: Mapping[str, object], name: str, path: Path
) -> None:
    outputs = manifest.get("outputs")
    entry = outputs.get(name) if isinstance(outputs, dict) else None
    expected = entry.get("sha256") if isinstance(entry, dict) else None
    observed = sha256_file(path)
    if expected != observed:
        raise ExpressionAuditError(
            f"Input hash mismatch for {name}: expected {expected}, observed {observed}"
        )


def _validate_records(
    records: Sequence[Mapping[str, str]], manifest: Mapping[str, object]
) -> None:
    if len(records) != 47:
        raise ExpressionAuditError(f"Expected 47 records, found {len(records)}")
    if len({row["sample_uid"] for row in records}) != 47:
        raise ExpressionAuditError("sample_uid values are not unique")
    if len({row["observation_id"] for row in records}) != 47:
        raise ExpressionAuditError("observation_id values are not unique")
    providers = Counter(row["provider_code"] for row in records)
    semantics = Counter(row["observation_semantics"] for row in records)
    if providers != Counter(EXPECTED_SOURCE_COUNTS):
        raise ExpressionAuditError(f"Provider counts differ from 23/8/16: {providers}")
    if semantics != Counter(EXPECTED_SEMANTICS_COUNTS):
        raise ExpressionAuditError(f"Semantics counts differ from 31/9/7: {semantics}")
    if manifest.get("record_count") != 47 or manifest.get("record_counts_by_source") != EXPECTED_SOURCE_COUNTS:
        raise ExpressionAuditError("Source manifest record counts are inconsistent")

    source = manifest.get("source")
    source_hash = source.get("sha256") if isinstance(source, dict) else None
    for row in records:
        if row["phenotype_name"] != "reported_yield":
            raise ExpressionAuditError("Unexpected phenotype name")
        if row["source_sha256"] != source_hash:
            raise ExpressionAuditError(f"Source hash mismatch in {row['sample_uid']}")
        observed_sequence_hash = hashlib.sha256(row["sequence_raw"].encode("ascii")).hexdigest()
        if observed_sequence_hash != row["sequence_sha256"]:
            raise ExpressionAuditError(f"Sequence hash mismatch in {row['sample_uid']}")
        if row["sequence_scope"] != "unknown" or row["vhh_region_sequence"]:
            raise ExpressionAuditError("Source sequences must remain untrimmed with unknown scope")
        if row["culture_volume_l"] != "1":
            raise ExpressionAuditError("All records must retain source-reported 1 L metadata")
        if any(row[field] for field in ("replicate_count", "uncertainty_value", "uncertainty_type")):
            raise ExpressionAuditError("Replicate and uncertainty fields must remain unknown")
        _validate_record_semantics(row)


def _validate_record_semantics(row: Mapping[str, str]) -> None:
    semantics = row["observation_semantics"]
    provider = row["provider_code"]
    if semantics == "individual_approximate":
        if provider not in {"LTT", "WCC"} or row["assignment_level"] != "individual":
            raise ExpressionAuditError("Individual numeric semantics/provider mismatch")
        _decimal(row["point_estimate_mg"], row["sample_uid"])
        if any(row[field] for field in ("group_anchor_mg", "lower_bound_mg", "upper_bound_mg")):
            raise ExpressionAuditError("Individual yield populated a group/bound field")
        if row["individual_numeric_available"] != "True":
            raise ExpressionAuditError("Individual numeric availability mismatch")
    elif semantics == "group_lower_bound":
        if provider != "LLJ" or row["assignment_level"] != "group":
            raise ExpressionAuditError("Lower-bound semantics/provider mismatch")
        if (row["point_estimate_mg"], row["group_anchor_mg"], row["lower_bound_mg"], row["upper_bound_mg"]) != ("", "", "20", ""):
            raise ExpressionAuditError("LLJ lower bound was changed or made exact")
        if row["value_relation"] != "gt" or row["censoring_type"] != "right_censored":
            raise ExpressionAuditError("LLJ lower-bound relation/censoring mismatch")
        if row["lower_bound_inclusive"] != "False" or row["individual_numeric_available"] != "False":
            raise ExpressionAuditError("LLJ lower-bound inclusivity/availability mismatch")
    elif semantics == "group_approximate":
        if provider != "LLJ" or row["assignment_level"] != "group":
            raise ExpressionAuditError("Group-anchor semantics/provider mismatch")
        if row["point_estimate_mg"] or row["lower_bound_mg"] or row["upper_bound_mg"]:
            raise ExpressionAuditError("LLJ group anchor was made into an individual value")
        if row["group_anchor_mg"] not in {"2", "10"}:
            raise ExpressionAuditError("Unexpected LLJ group anchor")
        if row["value_relation"] != "approx" or row["individual_numeric_available"] != "False":
            raise ExpressionAuditError("LLJ group-anchor relation/availability mismatch")
    else:
        raise ExpressionAuditError(f"Unknown observation semantics: {semantics}")


def _validate_contexts(
    contexts: Sequence[Mapping[str, str]], records: Sequence[Mapping[str, str]]
) -> None:
    if len(contexts) != 3 or {row["provider_code"] for row in contexts} != set(PROVIDER_ORDER):
        raise ExpressionAuditError("Assay context must contain exactly LTT, WCC, and LLJ")
    by_provider = {row["provider_code"]: row for row in contexts}
    for provider in PROVIDER_ORDER:
        context = by_provider[provider]
        if context["assay_id"] != f"ASSAY__{provider}" or context["culture_volume_l"] != "1":
            raise ExpressionAuditError(f"Invalid assay identity/volume for {provider}")
        if provider == "WCC":
            expected = ("entry_text", "TB", "post_purification")
        else:
            expected = ("document_title", "", "")
        if (context["volume_evidence"], context["medium"], context["yield_stage"]) != expected:
            raise ExpressionAuditError(f"Assay context changed for {provider}")
    if any(row["assay_id"] != f"ASSAY__{row['provider_code']}" for row in records):
        raise ExpressionAuditError("Record assay/provider mapping is inconsistent")
    for record in records:
        context = by_provider[record["provider_code"]]
        for field in ("culture_volume_l", "volume_evidence", "medium", "yield_stage"):
            if record[field] != context[field]:
                raise ExpressionAuditError(
                    f"Record/context {field} mismatch for {record['sample_uid']}"
                )


def _read_sequence_summary(
    path: Path | None, records: Sequence[Mapping[str, str]]
) -> dict[str, dict[str, str]]:
    if path is None:
        return {}
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    entries = payload.get("samples") if isinstance(payload, dict) else None
    if not isinstance(entries, list):
        raise ExpressionAuditError("Sequence audit summary must contain a samples array")
    known = {row["sample_uid"]: row for row in records}
    result: dict[str, dict[str, str]] = {}
    for entry in entries:
        if not isinstance(entry, dict) or not isinstance(entry.get("sample_uid"), str):
            raise ExpressionAuditError("Sequence audit sample lacks sample_uid")
        uid = entry["sample_uid"]
        if uid not in known or uid in result:
            raise ExpressionAuditError(f"Unknown or duplicate sequence-audit sample: {uid}")
        sequence_hash = entry.get("sequence_sha256")
        if sequence_hash is not None and sequence_hash != known[uid]["sequence_sha256"]:
            raise ExpressionAuditError(f"Stale sequence-audit input for {uid}")
        item: dict[str, str] = {}
        for field in (
            "numbering_status",
            "sequence_scope_status",
            "provisional_numbered_span_sha256",
        ):
            value = entry.get(field, "")
            if value is not None and not isinstance(value, str):
                raise ExpressionAuditError(f"Non-string {field} for {uid}")
            item[field] = value or ""
        span_hash = item["provisional_numbered_span_sha256"]
        if item["numbering_status"] and item["numbering_status"] not in {"pass", "failed"}:
            raise ExpressionAuditError(f"Unknown numbering status for {uid}")
        if span_hash and (len(span_hash) != 64 or any(c not in "0123456789abcdef" for c in span_hash.lower())):
            raise ExpressionAuditError(f"Invalid provisional span SHA-256 for {uid}")
        result[uid] = item
    return result


def _decimal(value: str, uid: str) -> Decimal:
    try:
        return Decimal(value)
    except InvalidOperation as exc:
        raise ExpressionAuditError(f"Invalid numeric value for {uid}: {value!r}") from exc
