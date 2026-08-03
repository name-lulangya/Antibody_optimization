"""Read and inventory traceable structure exports with Gemmi.

The helpers in this module preserve both author and label identifiers and never
assign molecular roles from chain order.  They are intended for ChimeraX mmCIF
exports that have already been recorded in an export manifest.  The module does
not mutate, renumber, complete, or write structures.
"""

from __future__ import annotations

import hashlib
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import gemmi
import numpy as np


STRUCTURE_INVENTORY_VERSION = "1.0.0"


class StructureBaselineError(ValueError):
    """Raised when a structure violates the baseline identity contract."""


@dataclass(frozen=True, order=True)
class ChainSelector:
    """An explicit author/label chain selection within one named source model."""

    model_name: str
    auth_asym_id: str
    label_asym_id: str


@dataclass(frozen=True, order=True)
class ResidueKey:
    """A reversible structure-residue identifier."""

    model_name: str
    auth_asym_id: str
    label_asym_id: str
    auth_seq_id: int
    insertion_code: str
    label_seq_id: int | None
    residue_name: str


@dataclass(frozen=True)
class StructureResidue:
    """One observed amino-acid residue from a confirmed chain."""

    key: ResidueKey
    one_letter_code: str
    observed_index_1based: int
    ca_coordinate: tuple[float, float, float] | None
    ca_altloc: str
    ca_status: str


def file_sha256(path: Path) -> str:
    """Return the lowercase SHA-256 digest of a concrete file."""

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_single_model_structure(path: Path) -> gemmi.Structure:
    """Read a structure and require one nonempty coordinate model.

    Entity/subchain information is read as written.  No ``setup_entities`` or
    label-sequence assignment is performed because those operations can invent
    identifiers that were absent from the exported file.
    """

    path = Path(path)
    if not path.is_file() or path.is_symlink():
        raise StructureBaselineError(
            f"Structure must be a regular non-symlink file: {path}"
        )
    try:
        structure = gemmi.read_structure(str(path))
    except (RuntimeError, OSError) as exc:
        raise StructureBaselineError(f"Gemmi could not read {path}: {exc}") from exc
    if len(structure) != 1:
        raise StructureBaselineError(
            f"Expected exactly one coordinate model in {path}, found {len(structure)}"
        )
    if structure[0].count_atom_sites() == 0:
        raise StructureBaselineError(f"Structure contains no atom sites: {path}")
    for chain in structure[0]:
        for residue in chain:
            for atom in residue:
                coordinates = (atom.pos.x, atom.pos.y, atom.pos.z)
                if not all(math.isfinite(value) for value in coordinates):
                    raise StructureBaselineError(
                        "Non-finite coordinate at "
                        f"{chain.name}/{residue.seqid}/{atom.name} in {path}"
                    )
    return structure


def prepare_heuristic_analysis_copy(
    exported_structure: gemmi.Structure,
) -> tuple[gemmi.Structure, dict[str, object]]:
    """Clone an export and let Gemmi establish analysis-only entity metadata.

    ``setup_entities`` and ``assign_label_seq_id`` are heuristic conveniences,
    not source annotations.  They are deliberately applied only to an in-memory
    clone; the ChimeraX export and the initially read structure remain unchanged.
    Callers must retain the returned metadata in their manifest.
    """

    analysis = exported_structure.clone()
    before = topology_signature(analysis)
    physical_before = physical_atom_site_signature(analysis)
    analysis.setup_entities()
    analysis.assign_label_seq_id()
    after = topology_signature(analysis)
    physical_after = physical_atom_site_signature(analysis)
    if physical_before != physical_after:
        raise StructureBaselineError(
            "Gemmi analysis setup changed physical atom-site identity or coordinates"
        )
    return analysis, {
        "analysis_copy_only": True,
        "source_export_modified": False,
        "gemmi_methods": ["Structure.setup_entities", "Structure.assign_label_seq_id"],
        "classification": "heuristic_entity_subchain_and_label_sequence_setup",
        "atom_site_count_before": len(before),
        "atom_site_count_after": len(after),
        "atom_site_count_preserved": len(before) == len(after),
        "physical_atom_site_signature_preserved": True,
        "identifier_signature_changed": before != after,
    }


def topology_signature(structure: gemmi.Structure) -> tuple[tuple[object, ...], ...]:
    """Return a coordinate-independent atom-site signature.

    The signature retains model, author/label chain, author/label residue,
    insertion, component, atom, element, alternate-location, and occupancy
    identity.  It is used to prove that native-frame and reference-frame
    exports differ only in coordinates.
    """

    rows: list[tuple[object, ...]] = []
    for model in structure:
        for chain in model:
            for residue in chain:
                label_seq = _optional_label_seq(residue)
                insertion_code = _clean_icode(residue.seqid.icode)
                for atom in residue:
                    rows.append(
                        (
                            str(model.num),
                            chain.name,
                            residue.subchain,
                            residue.entity_id,
                            residue.seqid.num,
                            insertion_code,
                            label_seq,
                            residue.name,
                            residue.het_flag,
                            atom.name,
                            _clean_altloc(atom.altloc),
                            atom.element.name,
                            round(float(atom.occ), 6),
                        )
                    )
    return tuple(rows)


def physical_atom_site_signature(
    structure: gemmi.Structure,
) -> tuple[tuple[object, ...], ...]:
    """Return physical atom identity while excluding entity/label metadata.

    This signature is used only to prove that analysis-only Gemmi entity setup
    did not alter residue order, atom identity, altloc, occupancy, or Cartesian
    coordinates.  Fields that the setup operation is expected to establish
    (entity ID, subchain/label asym, and label sequence ID) are omitted.
    """

    rows: list[tuple[object, ...]] = []
    for model_index, model in enumerate(structure):
        for chain_index, chain in enumerate(model):
            for residue_index, residue in enumerate(chain):
                for atom_index, atom in enumerate(residue):
                    rows.append(
                        (
                            model_index,
                            str(model.num),
                            chain_index,
                            chain.name,
                            residue_index,
                            residue.seqid.num,
                            _clean_icode(residue.seqid.icode),
                            residue.name,
                            residue.het_flag,
                            atom_index,
                            atom.name,
                            atom.element.name,
                            _clean_altloc(atom.altloc),
                            float(atom.occ),
                            float(atom.pos.x),
                            float(atom.pos.y),
                            float(atom.pos.z),
                        )
                    )
    return tuple(rows)


def require_matching_topology(
    native: gemmi.Structure,
    reference_frame: gemmi.Structure,
    *,
    native_label: str,
    reference_label: str,
) -> None:
    """Fail unless two exports have identical atom-site topology."""

    native_signature = topology_signature(native)
    reference_signature = topology_signature(reference_frame)
    if native_signature != reference_signature:
        mismatch_index = next(
            (
                index
                for index, (left, right) in enumerate(
                    zip(native_signature, reference_signature, strict=False)
                )
                if left != right
            ),
            min(len(native_signature), len(reference_signature)),
        )
        raise StructureBaselineError(
            "Native/reference-frame topology mismatch at atom-site index "
            f"{mismatch_index}: {native_label} ({len(native_signature)} sites) vs "
            f"{reference_label} ({len(reference_signature)} sites)"
        )


def atom_site_classification_counts(
    structure: gemmi.Structure,
) -> dict[str, int]:
    """Count Gemmi atom sites using the exporter occupancy categories."""

    counts: dict[str, int] = defaultdict(int)
    atom_objects_with_altlocs: set[tuple[int, int, int, str]] = set()
    for model_index, model in enumerate(structure):
        for chain_index, chain in enumerate(model):
            for residue_index, residue in enumerate(chain):
                for atom in residue:
                    altloc = _clean_altloc(atom.altloc)
                    counts["atom_site_count"] += 1
                    counts[
                        "nonblank_altloc_site_count" if altloc else "blank_altloc_site_count"
                    ] += 1
                    if altloc:
                        atom_objects_with_altlocs.add(
                            (model_index, chain_index, residue_index, atom.name)
                        )
                    occupancy_class = _occupancy_class(float(atom.occ))
                    counts[f"occupancy_{occupancy_class}_site_count"] += 1
    counts["atom_objects_with_nonblank_altlocs"] = len(atom_objects_with_altlocs)
    fields = (
        "atom_site_count", "blank_altloc_site_count", "nonblank_altloc_site_count",
        "atom_objects_with_nonblank_altlocs", "occupancy_negative_site_count",
        "occupancy_zero_site_count", "occupancy_partial_site_count",
        "occupancy_unit_site_count", "occupancy_above_unit_site_count",
        "occupancy_nonfinite_site_count",
    )
    return {field: int(counts[field]) for field in fields}


def structure_count_summary(structure: gemmi.Structure) -> dict[str, int]:
    """Return counts that can be compared with one ChimeraX source model.

    Gemmi stores each alternate-location coordinate as an atom site.  The
    logical atom-object count therefore groups sites by model, author/label
    chain, author/label residue identity, atom name, and element while omitting
    only the altloc identifier.  This mirrors ChimeraX's one ``Atom`` object
    with zero or more alternate-location coordinate states.
    """

    logical_atoms: set[tuple[object, ...]] = set()
    residue_count = 0
    chain_count = 0
    atom_site_count = 0
    for model_index, model in enumerate(structure):
        chain_count += len(model)
        for chain in model:
            for residue in chain:
                residue_count += 1
                label_seq_id = _optional_label_seq(residue)
                for atom in residue:
                    atom_site_count += 1
                    logical_atoms.add(
                        (
                            model_index,
                            chain.name,
                            residue.subchain,
                            residue.entity_id,
                            residue.seqid.num,
                            _clean_icode(residue.seqid.icode),
                            label_seq_id,
                            residue.name,
                            residue.het_flag,
                            atom.name,
                            atom.element.name,
                        )
                    )
    return {
        "model_count": len(structure),
        "chain_count": chain_count,
        "residue_count": residue_count,
        "atom_object_count": len(logical_atoms),
        "atom_site_count": atom_site_count,
    }


def rigid_coordinate_relationship(
    native: gemmi.Structure,
    reference_frame: gemmi.Structure,
    *,
    maximum_residual_tolerance_angstrom: float = 1e-3,
) -> dict[str, object]:
    """Validate and report the finite rigid transform between two exports."""

    require_matching_topology(
        native,
        reference_frame,
        native_label="native structure",
        reference_label="reference-frame structure",
    )
    native_coordinates = _all_atom_coordinates(native)
    reference_coordinates = _all_atom_coordinates(reference_frame)
    if len(native_coordinates) < 3:
        raise StructureBaselineError("Rigid-frame validation needs at least three atom sites")
    native_center = native_coordinates.mean(axis=0)
    reference_center = reference_coordinates.mean(axis=0)
    native_centered = native_coordinates - native_center
    reference_centered = reference_coordinates - reference_center
    if np.linalg.matrix_rank(native_centered) < 2:
        raise StructureBaselineError("Atom coordinates are collinear; rigid transform is not unique")
    covariance = native_centered.T @ reference_centered
    left, _, right_transpose = np.linalg.svd(covariance)
    rotation = right_transpose.T @ left.T
    if np.linalg.det(rotation) < 0:
        right_transpose[-1, :] *= -1
        rotation = right_transpose.T @ left.T
    translation = reference_center - rotation @ native_center
    transformed = (rotation @ native_coordinates.T).T + translation
    residuals = np.linalg.norm(transformed - reference_coordinates, axis=1)
    rmsd = math.sqrt(float(np.mean(residuals * residuals)))
    maximum = float(np.max(residuals))
    if not math.isfinite(rmsd) or not math.isfinite(maximum):
        raise StructureBaselineError("Rigid-frame validation produced non-finite residuals")
    if maximum > maximum_residual_tolerance_angstrom:
        raise StructureBaselineError(
            "Native/reference exports are not related by one rigid transform: "
            f"maximum residual {maximum:.6g} A exceeds "
            f"{maximum_residual_tolerance_angstrom:.6g} A"
        )
    return {
        "status": "pass",
        "algorithm": "all_atom_site_kabsch_no_outlier_rejection",
        "atom_site_count": len(native_coordinates),
        "rmsd_angstrom": rmsd,
        "maximum_residual_angstrom": maximum,
        "maximum_residual_tolerance_angstrom": maximum_residual_tolerance_angstrom,
        "rotation_3x3": [[float(value) for value in row] for row in rotation],
        "translation_3": [float(value) for value in translation],
    }


def chain_inventory_rows(
    structure: gemmi.Structure,
    *,
    source_model_name: str,
    source_file: Path,
    authoritative_sequence: str | None = None,
) -> list[dict[str, object]]:
    """Build subchain-level inventory rows without assigning molecule roles."""

    source_file = Path(source_file)
    model = structure[0]
    grouped: dict[tuple[str, str, str], list[gemmi.Residue]] = defaultdict(list)
    for chain in model:
        for residue in chain:
            grouped[(chain.name, residue.subchain, residue.entity_id)].append(residue)

    rows: list[dict[str, object]] = []
    for (auth_asym_id, label_asym_id, entity_id), residues in sorted(grouped.items()):
        entity_types = {_enum_name(residue.entity_type) for residue in residues}
        if len(entity_types) != 1:
            raise StructureBaselineError(
                "A single author/label/entity group contains multiple entity types: "
                f"{source_model_name} {auth_asym_id}/{label_asym_id}/{entity_id}"
            )
        amino_acids = [residue for residue in residues if _one_letter_code(residue.name)]
        observed_sequence = "".join(
            _one_letter_code(residue.name) for residue in amino_acids
        )
        entity = structure.get_entity(entity_id) if entity_id else None
        entity_full_sequence = (
            "".join(_one_letter_code(name) or "X" for name in entity.full_sequence)
            if entity is not None
            else ""
        )
        polymer_type = _enum_name(entity.polymer_type) if entity is not None else "Unknown"
        atom_count = sum(len(residue) for residue in residues)
        heavy_atom_count = sum(
            1 for residue in residues for atom in residue if _is_heavy_atom(atom)
        )
        first = residues[0]
        last = residues[-1]
        rows.append(
            {
                "source_model_name": source_model_name,
                "source_file": str(source_file),
                "source_file_sha256": file_sha256(source_file),
                "gemmi_model_num": str(model.num),
                "auth_asym_id": auth_asym_id,
                "label_asym_id": label_asym_id,
                "entity_id": entity_id,
                "entity_type": next(iter(entity_types)),
                "polymer_type": polymer_type,
                "entity_full_sequence": entity_full_sequence,
                "observed_sequence": observed_sequence,
                "observed_sequence_sha256": (
                    hashlib.sha256(observed_sequence.encode("ascii")).hexdigest()
                    if observed_sequence
                    else ""
                ),
                "residue_count": len(residues),
                "observed_amino_acid_count": len(amino_acids),
                "atom_count": atom_count,
                "heavy_atom_count": heavy_atom_count,
                "first_auth_seq_id": first.seqid.num,
                "first_insertion_code": _clean_icode(first.seqid.icode),
                "last_auth_seq_id": last.seqid.num,
                "last_insertion_code": _clean_icode(last.seqid.icode),
                "exact_observed_match_to_authoritative": (
                    bool(authoritative_sequence)
                    and observed_sequence == authoritative_sequence
                ),
                "exact_entity_match_to_authoritative": (
                    bool(authoritative_sequence)
                    and entity_full_sequence == authoritative_sequence
                ),
                "confirmed_role": "",
                "confirmation_status": "pending_user_review",
                "confirmed_by": "",
                "confirmed_at": "",
                "confirmation_note": "",
            }
        )
    if not rows:
        raise StructureBaselineError(f"No chain/subchain rows found in {source_file}")
    return rows


def residue_inventory_rows(
    structure: gemmi.Structure,
    *,
    source_model_name: str,
    source_file: Path | None = None,
) -> list[dict[str, object]]:
    """Return one provenance-preserving row per observed residue."""

    rows: list[dict[str, object]] = []
    source_path_text = "" if source_file is None else str(source_file)
    source_digest = "" if source_file is None else file_sha256(source_file)
    model = structure[0]
    for chain in model:
        for residue in chain:
            altlocs = sorted(
                {_clean_altloc(atom.altloc) for atom in residue if _clean_altloc(atom.altloc)}
            )
            occupancies = [float(atom.occ) for atom in residue]
            rows.append(
                {
                    "source_model_name": source_model_name,
                    "source_file": source_path_text,
                    "source_file_sha256": source_digest,
                    "gemmi_model_num": str(model.num),
                    "auth_asym_id": chain.name,
                    "label_asym_id": residue.subchain,
                    "entity_id": residue.entity_id,
                    "entity_type": _enum_name(residue.entity_type),
                    "auth_seq_id": residue.seqid.num,
                    "insertion_code": _clean_icode(residue.seqid.icode),
                    "label_seq_id": _optional_label_seq(residue),
                    "residue_name": residue.name,
                    "one_letter_code": _one_letter_code(residue.name),
                    "het_flag": residue.het_flag,
                    "is_water": residue.is_water(),
                    "atom_count": len(residue),
                    "heavy_atom_count": sum(
                        1 for atom in residue if _is_heavy_atom(atom)
                    ),
                    "altlocs": ";".join(altlocs),
                    "minimum_occupancy": min(occupancies) if occupancies else "",
                    "maximum_occupancy": max(occupancies) if occupancies else "",
                    "coordinate_status": "observed",
                }
            )
    return rows


def extract_confirmed_chain_residues(
    structure: gemmi.Structure,
    selector: ChainSelector,
) -> list[StructureResidue]:
    """Extract amino-acid residues from one explicit author/label chain.

    Residues whose entity type is not explicitly ``Polymer`` are rejected.  A
    C-alpha coordinate is returned only when exactly one positive-occupancy CA
    site exists; alternative CA sites are marked ambiguous instead of selected
    silently.
    """

    observed: list[StructureResidue] = []
    seen_keys: set[ResidueKey] = set()
    for chain in structure[0]:
        if chain.name != selector.auth_asym_id:
            continue
        for residue in chain:
            if residue.subchain != selector.label_asym_id:
                continue
            one_letter = _one_letter_code(residue.name)
            if not one_letter:
                continue
            if residue.entity_type != gemmi.EntityType.Polymer:
                raise StructureBaselineError(
                    "Confirmed VHH/NK2R chain contains an amino acid whose entity "
                    f"type is not Polymer: {selector} {residue.seqid} {residue.name}"
                )
            key = ResidueKey(
                model_name=selector.model_name,
                auth_asym_id=chain.name,
                label_asym_id=residue.subchain,
                auth_seq_id=residue.seqid.num,
                insertion_code=_clean_icode(residue.seqid.icode),
                label_seq_id=_optional_label_seq(residue),
                residue_name=residue.name,
            )
            if key in seen_keys:
                raise StructureBaselineError(f"Duplicate residue key: {key}")
            seen_keys.add(key)
            ca_sites = [
                atom
                for atom in residue
                if atom.name.strip() == "CA"
                and atom.element.name.upper() == "C"
                and float(atom.occ) > 0
            ]
            if len(ca_sites) == 1:
                ca = ca_sites[0]
                ca_coordinate = (float(ca.pos.x), float(ca.pos.y), float(ca.pos.z))
                ca_altloc = _clean_altloc(ca.altloc)
                ca_status = "observed_unique"
            elif not ca_sites:
                ca_coordinate = None
                ca_altloc = ""
                ca_status = "missing"
            else:
                ca_coordinate = None
                ca_altloc = ";".join(
                    sorted({_clean_altloc(atom.altloc) for atom in ca_sites})
                )
                ca_status = "ambiguous_altloc"
            observed.append(
                StructureResidue(
                    key=key,
                    one_letter_code=one_letter,
                    observed_index_1based=len(observed) + 1,
                    ca_coordinate=ca_coordinate,
                    ca_altloc=ca_altloc,
                    ca_status=ca_status,
                )
            )
    if not observed:
        raise StructureBaselineError(
            "Explicit chain selector matched no polymer amino-acid residues: "
            f"{selector}"
        )
    return observed


def observed_sequence(residues: Sequence[StructureResidue]) -> str:
    """Return the literal coordinate-observed amino-acid sequence."""

    return "".join(residue.one_letter_code for residue in residues)


def _one_letter_code(residue_name: str) -> str:
    info = gemmi.find_tabulated_residue(residue_name)
    if not info.found() or not info.is_amino_acid():
        return ""
    code = info.one_letter_code
    return code if len(code) == 1 and code.isalpha() else ""


def _is_heavy_atom(atom: gemmi.Atom) -> bool:
    return atom.element.name.upper() not in {"H", "D"}


def _optional_label_seq(residue: gemmi.Residue) -> int | None:
    value = residue.label_seq
    return int(value) if value is not None else None


def _clean_icode(value: str) -> str:
    return "" if not value or value in {" ", "\x00", "?", "."} else value


def _clean_altloc(value: str) -> str:
    return "" if not value or value in {" ", "\x00", ".", "?"} else value


def _enum_name(value: object) -> str:
    return str(value).rsplit(".", 1)[-1]


def _occupancy_class(value: float) -> str:
    if not math.isfinite(value):
        return "nonfinite"
    if value < -1e-6:
        return "negative"
    if math.isclose(value, 0.0, abs_tol=1e-6):
        return "zero"
    if value < 1.0 - 1e-6:
        return "partial"
    if math.isclose(value, 1.0, abs_tol=1e-6):
        return "unit"
    return "above_unit"


def _all_atom_coordinates(structure: gemmi.Structure) -> np.ndarray:
    coordinates = np.asarray(
        [
            (float(atom.pos.x), float(atom.pos.y), float(atom.pos.z))
            for model in structure
            for chain in model
            for residue in chain
            for atom in residue
        ],
        dtype=float,
    )
    if coordinates.ndim != 2 or coordinates.shape[1:] != (3,):
        raise StructureBaselineError("Atom coordinates must form an N x 3 array")
    if not np.isfinite(coordinates).all():
        raise StructureBaselineError("Atom coordinates contain non-finite values")
    return coordinates
