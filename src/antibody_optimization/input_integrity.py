"""Validate and freeze the stage-1 source identities.

The functions in this module bind the collaborator CXS, source DOCX, the
validated 47-row expression table, and the two stage-1 review manifests into
one machine-readable record.  They validate literal sequences and hashes but
do not infer VHH construct boundaries, chain roles, assay equivalence, or
structure provenance that is absent from the supplied inputs.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Mapping, Sequence


EXPECTED_CXS_SHA256 = (
    "1bc636c28f66ae60edc658d2e1c4aad0b07f4141ca5411c78662aa19da793c4d"
)
EXPECTED_RECORD_COUNT = 47
EXPECTED_PYTHON_MAJOR_MINOR = "3.11"
EXPECTED_ANARCII_VERSION = "2.0.8"
EXPECTED_GEMMI_VERSION = "0.7.5"
STANDARD_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")


class InputIntegrityError(ValueError):
    """Raised when a frozen source or its recorded provenance disagrees."""


def sha256_file(path: Path) -> str:
    """Return a lowercase SHA-256 digest without modifying ``path``."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_identity(path: Path) -> dict[str, object]:
    """Return content and filesystem identity fields used for change checks."""

    if path.is_symlink() or not path.is_file():
        raise InputIntegrityError(f"Expected a regular non-symlink file: {path}")
    stat = path.stat()
    return {
        "path": str(path),
        "sha256": sha256_file(path),
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def assert_same_identity(
    path: Path, expected: Mapping[str, object], *, label: str
) -> None:
    """Fail if content, size, or mtime changed since ``expected`` was recorded."""

    observed = file_identity(path)
    for field in ("sha256", "size_bytes", "mtime_ns"):
        if observed[field] != expected.get(field):
            raise InputIntegrityError(
                f"{label} changed during the run: {field} "
                f"{expected.get(field)!r} -> {observed[field]!r}"
            )


def build_input_freeze_manifest(
    *,
    source_cxs: Path,
    source_docx: Path,
    expression_manifest_path: Path,
    expression_manifest: Mapping[str, object],
    expression_records_path: Path,
    expression_records: Sequence[Mapping[str, str]],
    sequence_manifest_path: Path,
    sequence_manifest: Mapping[str, object],
    numbering_review_path: Path,
    numbering_positions_path: Path,
    expression_audit_manifest_path: Path,
    expression_audit_manifest: Mapping[str, object],
    sample_comparability_path: Path,
    generated_at: str,
    python_version: str,
    gemmi_version: str,
    expected_cxs_sha256: str = EXPECTED_CXS_SHA256,
) -> dict[str, object]:
    """Build a validated, reversible identity record for stage-1 inputs.

    ``expression_records`` must preserve every literal ``sequence_raw``.  The
    returned per-sample list is sorted by ``sample_uid`` only for stable
    serialization; source order remains available in the frozen CSV itself.
    """

    cxs_identity = file_identity(source_cxs)
    docx_identity = file_identity(source_docx)
    expression_manifest_identity = file_identity(expression_manifest_path)
    records_identity = file_identity(expression_records_path)
    sequence_manifest_identity = file_identity(sequence_manifest_path)
    numbering_review_identity = file_identity(numbering_review_path)
    numbering_positions_identity = file_identity(numbering_positions_path)
    audit_manifest_identity = file_identity(expression_audit_manifest_path)
    sample_comparability_identity = file_identity(sample_comparability_path)

    normalized_expected_cxs = _normalized_digest(expected_cxs_sha256)
    if normalized_expected_cxs != EXPECTED_CXS_SHA256:
        raise InputIntegrityError(
            "The requested CXS contract is not the project-frozen Nb252 session hash"
        )
    if cxs_identity["sha256"] != normalized_expected_cxs:
        raise InputIntegrityError(
            "Nb252 CXS SHA-256 mismatch: "
            f"expected {normalized_expected_cxs}, observed {cxs_identity['sha256']}"
        )

    source_record = _mapping(expression_manifest.get("source"), "expression source")
    _require_equal(
        docx_identity["sha256"],
        source_record.get("sha256"),
        "DOCX SHA-256 versus expression manifest",
    )
    _require_equal(
        docx_identity["size_bytes"],
        source_record.get("size_bytes"),
        "DOCX size versus expression manifest",
    )
    if Path(str(source_record.get("file", ""))).name != source_docx.name:
        raise InputIntegrityError("DOCX filename disagrees with the expression manifest")

    outputs = _mapping(expression_manifest.get("outputs"), "expression outputs")
    records_output = _mapping(
        outputs.get(expression_records_path.name), "expression records output"
    )
    _require_equal(
        records_identity["sha256"],
        records_output.get("sha256"),
        "Expression-record CSV SHA-256 versus manifest",
    )
    _require_equal(
        records_identity["size_bytes"],
        records_output.get("size_bytes"),
        "Expression-record CSV size versus manifest",
    )

    per_sequence = _validate_literal_sequences(expression_records, docx_identity)
    _validate_sequence_manifest(
        sequence_manifest,
        records_identity,
        per_sequence,
        numbering_review_path=numbering_review_path,
        numbering_review_identity=numbering_review_identity,
        numbering_positions_path=numbering_positions_path,
        numbering_positions_identity=numbering_positions_identity,
        project_root=source_cxs.parent,
    )
    _validate_expression_audit_manifest(
        expression_audit_manifest,
        records_identity,
        len(per_sequence),
        sample_comparability_path=sample_comparability_path,
        sample_comparability_identity=sample_comparability_identity,
        project_root=source_cxs.parent,
    )

    if ".".join(python_version.split(".")[:2]) != EXPECTED_PYTHON_MAJOR_MINOR:
        raise InputIntegrityError(
            f"Expected Python {EXPECTED_PYTHON_MAJOR_MINOR}.x, observed {python_version}"
        )
    if gemmi_version != EXPECTED_GEMMI_VERSION:
        raise InputIntegrityError(
            f"Expected Gemmi {EXPECTED_GEMMI_VERSION}, observed {gemmi_version}"
        )

    return {
        "schema_version": 1,
        "status": "pass",
        "generated_at": generated_at,
        "sources": {
            "collaborator_cxs": cxs_identity,
            "reported_yield_docx": docx_identity,
            "expression_manifest": expression_manifest_identity,
            "expression_records": records_identity,
            "sequence_numbering_manifest": sequence_manifest_identity,
            "sequence_numbering_review": numbering_review_identity,
            "sequence_numbering_positions": numbering_positions_identity,
            "expression_audit_manifest": audit_manifest_identity,
            "sample_comparability_review": sample_comparability_identity,
        },
        "sequence_identity": {
            "record_count": len(per_sequence),
            "unique_sample_uid_count": len(per_sequence),
            "literal_sequence_policy": "no normalization, trimming, or substitution",
            "per_sample": per_sequence,
        },
        "software_contract": {
            "python": python_version,
            "python_required_major_minor": EXPECTED_PYTHON_MAJOR_MINOR,
            "anarcii": EXPECTED_ANARCII_VERSION,
            "gemmi": gemmi_version,
            "chimerax_export_required": "1.12",
            "source_session_generator": {
                "chimerax": "1.9",
                "platform": "macOS",
                "evidence_status": "user_provided",
            },
        },
        "scientific_scope": {
            "reported_nb252_is_provisional": True,
            "construct_boundary_confirmed": False,
            "chain_roles_inferred": False,
            "assay_equivalence_inferred": False,
        },
    }


def _validate_literal_sequences(
    rows: Sequence[Mapping[str, str]], docx_identity: Mapping[str, object]
) -> list[dict[str, object]]:
    if len(rows) != EXPECTED_RECORD_COUNT:
        raise InputIntegrityError(
            f"Expected {EXPECTED_RECORD_COUNT} expression records, found {len(rows)}"
        )
    seen: set[str] = set()
    result: list[dict[str, object]] = []
    for row in rows:
        sample_uid = row.get("sample_uid", "")
        sequence = row.get("sequence_raw", "")
        if not sample_uid or sample_uid in seen:
            raise InputIntegrityError(f"Missing or duplicate sample_uid: {sample_uid!r}")
        seen.add(sample_uid)
        if not sequence or set(sequence) - STANDARD_AMINO_ACIDS:
            raise InputIntegrityError(f"Invalid literal sequence for {sample_uid}")
        try:
            recorded_length = int(row.get("sequence_length_aa", ""))
        except ValueError as exc:
            raise InputIntegrityError(
                f"Invalid recorded sequence length for {sample_uid}"
            ) from exc
        observed_digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
        if recorded_length != len(sequence):
            raise InputIntegrityError(f"Sequence length mismatch for {sample_uid}")
        if row.get("sequence_sha256") != observed_digest:
            raise InputIntegrityError(f"Sequence SHA-256 mismatch for {sample_uid}")
        if row.get("source_sha256") != docx_identity["sha256"]:
            raise InputIntegrityError(f"Source SHA-256 mismatch for {sample_uid}")
        result.append(
            {
                "sample_uid": sample_uid,
                "sequence_length_aa": len(sequence),
                "sequence_sha256": observed_digest,
            }
        )
    return sorted(result, key=lambda row: str(row["sample_uid"]))


def _validate_sequence_manifest(
    manifest: Mapping[str, object],
    records_identity: Mapping[str, object],
    per_sequence: Sequence[Mapping[str, object]],
    *,
    numbering_review_path: Path,
    numbering_review_identity: Mapping[str, object],
    numbering_positions_path: Path,
    numbering_positions_identity: Mapping[str, object],
    project_root: Path,
) -> None:
    if manifest.get("status") != "pass":
        raise InputIntegrityError("Sequence-numbering manifest is not pass")
    tool = _mapping(manifest.get("tool"), "ANARCII tool metadata")
    if tool.get("version") != EXPECTED_ANARCII_VERSION:
        raise InputIntegrityError("ANARCII version disagrees with the frozen contract")
    inputs = _mapping(manifest.get("input"), "sequence-numbering input")
    records = _mapping(inputs.get("records"), "sequence-numbering records input")
    _require_equal(
        records.get("sha256"),
        records_identity["sha256"],
        "Sequence-numbering input versus expression-record CSV",
    )
    samples = manifest.get("samples")
    if not isinstance(samples, list) or len(samples) != EXPECTED_RECORD_COUNT:
        raise InputIntegrityError("Sequence-numbering manifest must contain 47 samples")
    expected_hashes = {
        str(row["sample_uid"]): row["sequence_sha256"] for row in per_sequence
    }
    observed_hashes: dict[str, object] = {}
    for row in samples:
        if not isinstance(row, Mapping):
            raise InputIntegrityError("Invalid sequence-numbering sample record")
        uid = str(row.get("sample_uid", ""))
        if not uid or uid in observed_hashes:
            raise InputIntegrityError("Duplicate sequence-numbering sample identity")
        observed_hashes[uid] = row.get("sequence_sha256")
    if observed_hashes != expected_hashes:
        raise InputIntegrityError("Per-sequence hashes disagree with numbering manifest")
    validate_manifest_output_file(
        numbering_review_path,
        numbering_review_identity,
        manifest,
        output_key="sequence_review",
        label="sequence-numbering review",
        project_root=project_root,
    )
    validate_manifest_output_file(
        numbering_positions_path,
        numbering_positions_identity,
        manifest,
        output_key="positions",
        label="sequence-numbering positions",
        project_root=project_root,
    )


def _validate_expression_audit_manifest(
    manifest: Mapping[str, object],
    records_identity: Mapping[str, object],
    record_count: int,
    *,
    sample_comparability_path: Path,
    sample_comparability_identity: Mapping[str, object],
    project_root: Path,
) -> None:
    gates = _mapping(manifest.get("gates"), "expression-audit gates")
    if gates.get("expression_audit_gate") != "pass":
        raise InputIntegrityError("Expression audit gate is not pass")
    counts = _mapping(manifest.get("counts"), "expression-audit counts")
    _require_equal(counts.get("samples"), record_count, "Expression-audit count")
    inputs = _mapping(manifest.get("inputs"), "expression-audit inputs")
    _require_equal(
        inputs.get("nb_expression_records.csv"),
        records_identity["sha256"],
        "Expression-audit input versus expression-record CSV",
    )
    validate_manifest_output_file(
        sample_comparability_path,
        sample_comparability_identity,
        manifest,
        output_key=sample_comparability_path.name,
        label="sample comparability review",
        project_root=project_root,
    )


def validate_manifest_output_file(
    path: Path,
    identity: Mapping[str, object],
    manifest: Mapping[str, object],
    *,
    output_key: str,
    label: str,
    project_root: Path | None = None,
) -> None:
    """Bind a concrete derived file to one manifest output record.

    The output record must contain a matching SHA-256 and, when present, a
    matching size.  A recorded path/file field must name the same file.  This
    prevents a passed manifest from being paired with an edited replacement.
    """

    outputs = _mapping(manifest.get("outputs"), f"{label} outputs")
    raw_record = outputs.get(output_key)
    if isinstance(raw_record, str):
        record: Mapping[str, object] = {"sha256": raw_record}
    else:
        record = _mapping(raw_record, f"{label} output record")
    _require_equal(identity.get("sha256"), record.get("sha256"), f"{label} SHA-256")
    if "size_bytes" in record:
        _require_equal(identity.get("size_bytes"), record.get("size_bytes"), f"{label} size")
    recorded_path = record.get("path", record.get("file"))
    if recorded_path is not None:
        recorded = Path(str(recorded_path))
        if recorded.name != path.name:
            raise InputIntegrityError(
                f"{label} filename mismatch: {recorded.name!r} != {path.name!r}"
            )
        if project_root is not None:
            recorded_absolute = (
                recorded if recorded.is_absolute() else project_root / recorded
            ).resolve(strict=False)
            if recorded_absolute != path.resolve(strict=True):
                raise InputIntegrityError(
                    f"{label} path mismatch: {recorded_absolute} != {path.resolve(strict=True)}"
                )


def _mapping(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise InputIntegrityError(f"Missing or invalid {label}")
    return value


def _normalized_digest(value: str) -> str:
    normalized = value.strip().lower()
    if len(normalized) != 64 or any(c not in "0123456789abcdef" for c in normalized):
        raise InputIntegrityError(f"Invalid SHA-256 value: {value!r}")
    return normalized


def _require_equal(observed: object, expected: object, label: str) -> None:
    if observed != expected:
        raise InputIntegrityError(
            f"{label} mismatch: observed {observed!r}, expected {expected!r}"
        )
