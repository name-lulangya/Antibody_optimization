"""Reusable validation helpers for the Nb252 structure-baseline workflow.

The functions here validate ChimeraX export contracts and preserve the
distinction between source-declared mmCIF polymer evidence and metadata that
Gemmi adds only to an in-memory analysis copy.  They do not assign chain roles
or modify a source/exported structure.
"""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

import gemmi

from .polymer_mapping import read_source_polymer_evidence
from .structure_inventory import (
    ChainSelector,
    StructureBaselineError,
    StructureResidue,
    extract_confirmed_chain_residues,
    file_sha256,
)


class StructureBaselineSupportError(ValueError):
    """Raised when an exported baseline input violates its declared contract."""


def annotate_source_polymer_inventory(
    rows: Sequence[dict[str, object]],
    *,
    native_path: Path,
    reported_sequence: str,
) -> None:
    """Add source-mmCIF polymer evidence without relabeling heuristic metadata.

    ``rows`` is updated in place because it is the inventory being assembled by
    the caller.  Missing source categories remain explicitly marked ``absent``;
    analysis-copy entity metadata is never promoted to source evidence.
    """

    for row in rows:
        row.update(
            {
                "entity_full_sequence_provenance": (
                    "gemmi_analysis_copy_after_setup_entities"
                ),
                "source_entity_id": "",
                "source_entity_description": "",
                "source_entity_type": "",
                "source_polymer_type": "",
                "source_polymer_sequence": "",
                "source_polymer_sequence_sha256": "",
                "source_polymer_sequence_provenance": "absent",
                "source_struct_asym_provenance": "absent",
                "source_poly_seq_scheme_provenance": "absent",
                "source_poly_seq_scheme_row_count": 0,
                "exact_source_polymer_match_to_reported": False,
            }
        )
        if str(row.get("entity_type", "")).lower() != "polymer":
            continue
        label_asym_id = str(row.get("label_asym_id", ""))
        if not label_asym_id:
            continue
        evidence = read_source_polymer_evidence(
            native_path,
            label_asym_id=label_asym_id,
        )
        if evidence is None:
            continue
        row.update(
            {
                "source_entity_id": evidence.entity_id,
                "source_entity_description": evidence.entity_description,
                "source_entity_type": evidence.entity_type,
                "source_polymer_type": evidence.polymer_type,
                "source_polymer_sequence": evidence.polymer_sequence,
                "source_polymer_sequence_sha256": evidence.polymer_sequence_sha256,
                "source_polymer_sequence_provenance": (
                    evidence.polymer_sequence_source
                ),
                "source_struct_asym_provenance": evidence.struct_asym_source,
                "source_poly_seq_scheme_provenance": evidence.scheme_source,
                "source_poly_seq_scheme_row_count": len(evidence.scheme_rows),
                "exact_source_polymer_match_to_reported": (
                    evidence.polymer_sequence == reported_sequence
                ),
            }
        )


def mapping_residues_and_label_source(
    *,
    raw_structure: gemmi.Structure,
    analysis_structure: gemmi.Structure,
    selector: ChainSelector,
) -> tuple[list[StructureResidue], str]:
    """Select observed residues and label their label-sequence-ID provenance."""

    raw_selector_exists = any(
        chain.name == selector.auth_asym_id
        and residue.subchain == selector.label_asym_id
        and gemmi.find_tabulated_residue(residue.name).is_amino_acid()
        for chain in raw_structure[0]
        for residue in chain
    )
    if raw_selector_exists:
        residues = extract_confirmed_chain_residues(raw_structure, selector)
        source_kind = "source_mmcif_atom_site"
    else:
        residues = extract_confirmed_chain_residues(analysis_structure, selector)
        source_kind = "gemmi_heuristic"
    label_ids = [residue.key.label_seq_id for residue in residues]
    present = [value is not None for value in label_ids]
    if any(present) and not all(present):
        raise StructureBaselineError(
            f"Selected chain mixes present and absent label_seq_id values: {selector}"
        )
    if not any(present):
        source_kind = "absent"
    return residues, source_kind


def validate_export_records(
    manifest: Mapping[str, object],
    directory: Path,
    *,
    expected_models: Sequence[str],
    reference_model: str,
) -> dict[tuple[str, str], Path]:
    """Validate the exact native/reference-frame export set and file hashes."""

    records = manifest.get("exports")
    if not isinstance(records, list):
        raise StructureBaselineSupportError("CXS manifest exports is not a list")
    result: dict[tuple[str, str], Path] = {}
    for record in records:
        if not isinstance(record, dict):
            raise StructureBaselineSupportError(
                "CXS export record is not an object"
            )
        key = (
            str(record.get("source_model_name", "")),
            str(record.get("coordinate_frame_kind", "")),
        )
        if key in result:
            raise StructureBaselineSupportError(
                f"duplicate CXS export record {key!r}"
            )
        path = (directory / str(record.get("path", ""))).resolve(strict=True)
        if (
            directory.resolve(strict=True) not in path.parents
            or path.is_symlink()
            or not path.is_file()
        ):
            raise StructureBaselineSupportError(f"unsafe CXS export path: {path}")
        if file_sha256(path) != record.get("sha256"):
            raise StructureBaselineSupportError(
                f"CXS export hash mismatch: {path}"
            )
        result[key] = path
    expected = {(name, "native_model_frame") for name in expected_models}
    expected |= {
        (name, "experimental_reference_model_frame")
        for name in expected_models
        if name != reference_model
    }
    if set(result) != expected:
        raise StructureBaselineSupportError(
            "CXS export record set mismatch: "
            f"missing={sorted(expected-set(result))!r}, "
            f"extra={sorted(set(result)-expected)!r}"
        )
    return result


def export_model_site_counts(
    manifest: Mapping[str, object], *, expected_models: Sequence[str]
) -> dict[str, dict[str, int]]:
    """Read the exact per-model logical atom-site count contract."""

    models = manifest.get("models")
    if not isinstance(models, list):
        raise StructureBaselineSupportError("CXS manifest models is not a list")
    result: dict[str, dict[str, int]] = {}
    for record in models:
        if not isinstance(record, dict):
            raise StructureBaselineSupportError(
                "CXS model record is not an object"
            )
        name = str(record.get("model_name", ""))
        raw_counts = record.get("atom_site_classification_counts")
        if name in result or not isinstance(raw_counts, dict):
            raise StructureBaselineSupportError(
                f"missing or duplicate atom-site counts for model {name!r}"
            )
        try:
            result[name] = {
                str(key): int(value) for key, value in raw_counts.items()
            }
        except (TypeError, ValueError) as exc:
            raise StructureBaselineSupportError(
                f"invalid atom-site count value for model {name!r}"
            ) from exc
    if set(result) != set(expected_models):
        raise StructureBaselineSupportError(
            "CXS model atom-site count set is incomplete"
        )
    return result


def export_model_structure_counts(
    manifest: Mapping[str, object], *, expected_models: Sequence[str]
) -> dict[str, dict[str, int]]:
    """Read the model/chain/residue/logical-atom cross-tool count contract."""

    models = manifest.get("models")
    if not isinstance(models, list):
        raise StructureBaselineSupportError("CXS manifest models is not a list")
    result: dict[str, dict[str, int]] = {}
    for record in models:
        if not isinstance(record, dict):
            raise StructureBaselineSupportError(
                "CXS model record is not an object"
            )
        name = str(record.get("model_name", ""))
        site_counts = record.get("atom_site_classification_counts")
        if name in result or not isinstance(site_counts, dict):
            raise StructureBaselineSupportError(
                f"missing or duplicate structure counts for model {name!r}"
            )
        try:
            result[name] = {
                "model_count": int(record["coordinate_set_count"]),
                "chain_count": int(record["chain_count"]),
                "residue_count": int(record["residue_count"]),
                "atom_object_count": int(record["atom_count"]),
                "atom_site_count": int(site_counts["atom_site_count"]),
            }
        except (KeyError, TypeError, ValueError) as exc:
            raise StructureBaselineSupportError(
                f"invalid model/chain/residue/atom count for model {name!r}"
            ) from exc
    if set(result) != set(expected_models):
        raise StructureBaselineSupportError(
            "CXS model structure-count set is incomplete"
        )
    return result


def validated_color_inventory(
    export_manifest: Mapping[str, object], export_dir: Path
) -> Path:
    """Resolve and hash-check the color inventory declared by an export."""

    record = export_manifest.get("color_inventory")
    if not isinstance(record, dict):
        raise StructureBaselineSupportError("CXS manifest lacks color_inventory")
    lexical = export_dir / str(record.get("path", ""))
    if lexical.is_symlink() or not lexical.is_file():
        raise StructureBaselineSupportError(
            f"CXS color inventory is not a regular file: {lexical}"
        )
    path = lexical.resolve(strict=True)
    if file_sha256(path) != record.get("sha256"):
        raise StructureBaselineSupportError("CXS color inventory hash mismatch")
    return path
