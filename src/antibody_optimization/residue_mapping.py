"""Conservative sequence/structure mapping and explicit framework alignment.

Sequence mapping is exact-first and mismatch-free.  A coordinate-observed
sequence may omit residues.  When its exact order-preserving embedding is not
unique, a fixed BLOSUM62 global alignment may resolve the index placement only
if every optimal alignment gives the same complete mapping and all mapped
wild-type residues agree.  Other ambiguity or disagreement blocks the workflow.
"""

from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass
from typing import Iterable, Mapping, Protocol, Sequence

import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices

from .structure_inventory import ChainSelector, StructureResidue


RESIDUE_MAPPING_VERSION = "1.0.0"
STANDARD_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
GLOBAL_ALIGNMENT_PARAMETERS = {
    "mode": "global",
    "substitution_matrix": "BLOSUM62",
    "open_gap_score": -10.0,
    "extend_gap_score": -0.5,
    "terminal_gap_score": 0.0,
}


class ResidueMappingError(ValueError):
    """Raised when a residue mapping is missing, inconsistent, or ambiguous."""


@dataclass(frozen=True)
class ExactSequenceMapping:
    """Unique authoritative indices for coordinate-observed residues."""

    authoritative_index_1based_by_observed_index: tuple[int, ...]
    method: str
    optimal_global_alignment_count: int = 0


class ResidueIndexMapping(Protocol):
    """Minimal mapping contract shared by direct and source-aware routes."""

    authoritative_index_1based_by_observed_index: tuple[int, ...]

    @property
    def method(self) -> str: ...


@dataclass(frozen=True)
class RigidAlignment:
    """Kabsch transform that maps mobile coordinates into the reference frame."""

    rotation: tuple[tuple[float, float, float], ...]
    translation: tuple[float, float, float]
    rmsd_angstrom: float
    fitted_atom_count: int
    fitted_authoritative_indices_1based: tuple[int, ...]
    algorithm: str = "kabsch_framework_ca_no_outlier_rejection"

    def homogeneous_matrix(self) -> tuple[tuple[float, ...], ...]:
        return (
            (*self.rotation[0], self.translation[0]),
            (*self.rotation[1], self.translation[1]),
            (*self.rotation[2], self.translation[2]),
            (0.0, 0.0, 0.0, 1.0),
        )

    def as_dict(self) -> dict[str, object]:
        return {
            "algorithm": self.algorithm,
            "atom_selection": "explicit mapped framework C-alpha atoms",
            "outlier_rejection": "none",
            "fitted_atom_count": self.fitted_atom_count,
            "fitted_authoritative_indices_1based": list(
                self.fitted_authoritative_indices_1based
            ),
            "rmsd_angstrom": self.rmsd_angstrom,
            "rotation_3x3": [list(row) for row in self.rotation],
            "translation_3": list(self.translation),
            "homogeneous_transform_4x4": [
                list(row) for row in self.homogeneous_matrix()
            ],
        }


def unique_exact_order_mapping(
    authoritative_sequence: str,
    observed_sequence: str,
) -> ExactSequenceMapping:
    """Map ``observed_sequence`` uniquely as a mismatch-free subsequence.

    This supports unresolved structural residues while refusing residue
    substitutions and ambiguous repeated embeddings.  Both sequences must use
    the 20 standard uppercase amino-acid letters.
    """

    _validate_sequence(authoritative_sequence, label="authoritative sequence")
    _validate_sequence(observed_sequence, label="observed structure sequence")
    if len(observed_sequence) > len(authoritative_sequence):
        raise ResidueMappingError(
            "Observed structure sequence is longer than the authoritative sequence"
        )

    reference_length = len(authoritative_sequence)
    observed_length = len(observed_sequence)
    counts = [[0] * (observed_length + 1) for _ in range(reference_length + 1)]
    for reference_index in range(reference_length + 1):
        counts[reference_index][observed_length] = 1
    for reference_index in range(reference_length - 1, -1, -1):
        for observed_index in range(observed_length - 1, -1, -1):
            ways = counts[reference_index + 1][observed_index]
            if (
                authoritative_sequence[reference_index]
                == observed_sequence[observed_index]
            ):
                ways += counts[reference_index + 1][observed_index + 1]
            counts[reference_index][observed_index] = min(2, ways)

    embedding_count = counts[0][0]
    if embedding_count == 0:
        raise ResidueMappingError(
            "Observed residues are not a mismatch-free subsequence of the "
            "authoritative sequence; an explicit reviewed mapping is required"
        )
    if embedding_count > 1:
        raise ResidueMappingError(
            "Observed residues have multiple order-preserving embeddings in the "
            "authoritative sequence; an explicit reviewed mapping is required"
        )

    mapping: list[int] = []
    reference_index = 0
    observed_index = 0
    while observed_index < observed_length:
        can_match = (
            reference_index < reference_length
            and authoritative_sequence[reference_index]
            == observed_sequence[observed_index]
            and counts[reference_index + 1][observed_index + 1] > 0
        )
        can_skip = (
            reference_index < reference_length
            and counts[reference_index + 1][observed_index] > 0
        )
        if can_match and can_skip:
            raise ResidueMappingError(
                "Internal error: unique mapping reconstruction encountered a branch"
            )
        if can_match:
            mapping.append(reference_index + 1)
            reference_index += 1
            observed_index += 1
        elif can_skip:
            reference_index += 1
        else:
            raise ResidueMappingError(
                "Internal error: unique mapping could not be reconstructed"
            )

    if authoritative_sequence == observed_sequence:
        method = "exact_full_sequence"
    else:
        first = authoritative_sequence.find(observed_sequence)
        second = authoritative_sequence.find(observed_sequence, first + 1)
        method = (
            "unique_contiguous_subsequence"
            if first >= 0 and second < 0
            else "unique_order_preserving_subsequence"
        )
    return ExactSequenceMapping(tuple(mapping), method)


def exact_first_sequence_mapping(
    authoritative_sequence: str,
    observed_sequence: str,
) -> ExactSequenceMapping:
    """Use exact order mapping first, then a fixed BLOSUM62 global fallback.

    The fallback accepts only a single *mapping result* across every optimal
    alignment.  Every observed residue must be aligned to the identical
    authoritative residue; insertions, substitutions, and equally optimal
    alternatives with different observed-to-reference indices are rejected.
    Terminal gaps are unpenalized because unbuilt structure termini are not
    sequence evidence.
    """

    exact_failure: ResidueMappingError | None = None
    try:
        return unique_exact_order_mapping(authoritative_sequence, observed_sequence)
    except ResidueMappingError as exact_error:
        exact_failure = exact_error
        _validate_sequence(authoritative_sequence, label="authoritative sequence")
        _validate_sequence(observed_sequence, label="observed structure sequence")

    aligner = PairwiseAligner(mode="global")
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = GLOBAL_ALIGNMENT_PARAMETERS["open_gap_score"]
    aligner.extend_gap_score = GLOBAL_ALIGNMENT_PARAMETERS["extend_gap_score"]
    aligner.end_gap_score = GLOBAL_ALIGNMENT_PARAMETERS["terminal_gap_score"]
    alignments = aligner.align(authoritative_sequence, observed_sequence)
    accepted_mapping: tuple[int, ...] | None = None
    alignment_count = 0
    for alignment in alignments:
        alignment_count += 1
        mapping = _observed_mapping_from_alignment(
            alignment.coordinates,
            authoritative_sequence=authoritative_sequence,
            observed_sequence=observed_sequence,
        )
        if mapping is None:
            raise ResidueMappingError(
                "An optimal BLOSUM62 global alignment contains an observed insertion "
                "or a wild-type mismatch; exact source identity is not established"
            ) from exact_failure
        if accepted_mapping is None:
            accepted_mapping = mapping
        elif mapping != accepted_mapping:
            raise ResidueMappingError(
                "Optimal BLOSUM62 global alignments imply different observed-to-"
                "authoritative residue mappings"
            ) from exact_failure
    if accepted_mapping is None or alignment_count == 0:
        raise ResidueMappingError("BLOSUM62 global alignment produced no optimum") from exact_failure
    return ExactSequenceMapping(
        accepted_mapping,
        "blosum62_global_unique_optimal_mapping",
        optimal_global_alignment_count=alignment_count,
    )


def validate_numbering_rows(
    numbering_rows: Sequence[Mapping[str, object]],
    *,
    sample_uid: str,
    authoritative_sequence: str,
    authoritative_sequence_sha256: str,
    required_scheme: str,
) -> dict[int, dict[str, str]]:
    """Validate provisional numbering rows and index them by source position.

    Gap rows without an authoritative sequence index are retained by the source
    artifact but are intentionally not join keys for structure residues.
    """

    calculated_hash = hashlib.sha256(authoritative_sequence.encode("ascii")).hexdigest()
    if calculated_hash != authoritative_sequence_sha256.lower():
        raise ResidueMappingError(
            "Authoritative sequence SHA-256 does not match the supplied sequence"
        )
    indexed: dict[int, dict[str, str]] = {}
    matched_rows = 0
    for raw_row in numbering_rows:
        row = {str(key): "" if value is None else str(value) for key, value in raw_row.items()}
        if row.get("sample_uid") != sample_uid:
            continue
        matched_rows += 1
        if row.get("sequence_sha256", "").lower() != authoritative_sequence_sha256.lower():
            raise ResidueMappingError(
                f"Numbering row has the wrong sequence hash for {sample_uid}"
            )
        if row.get("scheme", "").lower() != required_scheme.lower():
            raise ResidueMappingError(
                f"Numbering scheme mismatch: expected {required_scheme!r}, "
                f"got {row.get('scheme')!r}"
            )
        index_text = row.get("sequence_index_1based", "").strip()
        if not index_text:
            if row.get("is_gap", "").strip().lower() not in {"true", "1", "yes"}:
                raise ResidueMappingError(
                    "A numbering row without sequence_index_1based is not marked as a gap"
                )
            continue
        try:
            index = int(index_text)
        except ValueError as exc:
            raise ResidueMappingError(
                f"Invalid numbering sequence index: {index_text!r}"
            ) from exc
        if not 1 <= index <= len(authoritative_sequence):
            raise ResidueMappingError(f"Numbering sequence index is out of range: {index}")
        if index in indexed:
            raise ResidueMappingError(f"Duplicate numbering row for sequence index {index}")
        if row.get("residue_aa") != authoritative_sequence[index - 1]:
            raise ResidueMappingError(
                "Numbering wild-type residue mismatch at authoritative index "
                f"{index}: {row.get('residue_aa')!r} vs "
                f"{authoritative_sequence[index - 1]!r}"
            )
        if row.get("numbering_status") not in {"success", "provisional", "pass"}:
            raise ResidueMappingError(
                f"Numbering row {index} is not successful/provisional: "
                f"{row.get('numbering_status')!r}"
            )
        indexed[index] = row
    if matched_rows == 0:
        raise ResidueMappingError(f"No numbering rows found for {sample_uid}")
    if not indexed:
        raise ResidueMappingError(f"No residue-bearing numbering rows found for {sample_uid}")
    return indexed


def build_mapping_rows(
    *,
    sample_uid: str,
    authoritative_sequence: str,
    authoritative_sequence_sha256: str,
    selector: ChainSelector,
    structure_residues: Sequence[StructureResidue],
    numbering_by_index: Mapping[int, Mapping[str, str]],
    numbering_scheme: str,
    sequence_mapping: ResidueIndexMapping | None = None,
) -> tuple[list[dict[str, object]], ResidueIndexMapping]:
    """Build one reversible mapping row per authoritative sequence residue."""

    observed = "".join(residue.one_letter_code for residue in structure_residues)
    if sequence_mapping is None:
        sequence_mapping = exact_first_sequence_mapping(
            authoritative_sequence, observed
        )
    if len(sequence_mapping.authoritative_index_1based_by_observed_index) != len(
        structure_residues
    ):
        raise ResidueMappingError(
            "Supplied residue mapping length does not match observed structure residues"
        )
    mapped_indices = sequence_mapping.authoritative_index_1based_by_observed_index
    if (
        any(index < 1 or index > len(authoritative_sequence) for index in mapped_indices)
        or len(set(mapped_indices)) != len(mapped_indices)
        or any(left >= right for left, right in zip(mapped_indices, mapped_indices[1:]))
    ):
        raise ResidueMappingError(
            "Supplied residue mapping indices are out of range, duplicate, or nonmonotonic"
        )
    for index, residue in zip(mapped_indices, structure_residues, strict=True):
        if authoritative_sequence[index - 1] != residue.one_letter_code:
            raise ResidueMappingError(
                "Supplied residue mapping has a wild-type mismatch at authoritative "
                f"index {index}"
            )
    observed_by_authoritative = {
        authoritative_index: residue
        for authoritative_index, residue in zip(
            sequence_mapping.authoritative_index_1based_by_observed_index,
            structure_residues,
            strict=True,
        )
    }
    rows: list[dict[str, object]] = []
    for index, residue_aa in enumerate(authoritative_sequence, start=1):
        structure_residue = observed_by_authoritative.get(index)
        numbering = numbering_by_index.get(index, {})
        key = structure_residue.key if structure_residue else None
        terminal_flank = not numbering or numbering.get(
            "numbering_status", ""
        ) == "outside_numbered_domain"
        coordinate_status = (
            "terminal_flank"
            if terminal_flank
            else "observed"
            if structure_residue
            else "missing_coordinates"
        )
        rows.append(
            {
                "sample_uid": sample_uid,
                "authoritative_sequence_sha256": authoritative_sequence_sha256,
                "sequence_index_1based": index,
                "sequence_index_0based": index - 1,
                "residue_aa": residue_aa,
                "numbering_scheme": numbering_scheme,
                "numbering_position": numbering.get("numbering_position", ""),
                "numbering_insertion_code": numbering.get("insertion_code", ""),
                "numbering_position_label": numbering.get(
                    "numbering_position_label", ""
                ),
                "region": numbering.get("region", ""),
                "numbering_status": numbering.get(
                    "numbering_status", "outside_numbered_domain"
                ),
                "source_model_name": selector.model_name,
                "auth_asym_id": selector.auth_asym_id,
                "label_asym_id": selector.label_asym_id,
                "auth_seq_id": key.auth_seq_id if key else "",
                "insertion_code": key.insertion_code if key else "",
                "label_seq_id": (
                    "" if key is None or key.label_seq_id is None else key.label_seq_id
                ),
                "structure_residue_name": key.residue_name if key else "",
                "structure_residue_aa": (
                    structure_residue.one_letter_code if structure_residue else ""
                ),
                "coordinate_status": coordinate_status,
                "coordinate_evaluable": bool(structure_residue),
                "ca_status": structure_residue.ca_status if structure_residue else "",
                "ca_altloc": structure_residue.ca_altloc if structure_residue else "",
                "mapping_method": sequence_mapping.method,
                "optimal_global_alignment_count": (
                    getattr(sequence_mapping, "optimal_global_alignment_count", "")
                ),
                "mapping_status": "provisional_identity_consistent",
            }
        )
    return rows, sequence_mapping


def ca_coordinates_by_authoritative_index(
    mapping: ResidueIndexMapping,
    structure_residues: Sequence[StructureResidue],
) -> dict[int, tuple[float, float, float]]:
    """Return uniquely observed C-alpha coordinates keyed by source sequence index."""

    coordinates: dict[int, tuple[float, float, float]] = {}
    for index, residue in zip(
        mapping.authoritative_index_1based_by_observed_index,
        structure_residues,
        strict=True,
    ):
        if residue.ca_coordinate is not None:
            coordinates[index] = residue.ca_coordinate
    return coordinates


def fit_explicit_framework_ca(
    reference_coordinates: Mapping[int, Sequence[float]],
    mobile_coordinates: Mapping[int, Sequence[float]],
    *,
    framework_authoritative_indices_1based: Iterable[int],
) -> RigidAlignment:
    """Fit mobile onto reference using exactly the supplied framework C-alpha set."""

    requested = tuple(sorted(set(framework_authoritative_indices_1based)))
    fitted = tuple(
        index
        for index in requested
        if index in reference_coordinates and index in mobile_coordinates
    )
    if len(fitted) < 3:
        raise ResidueMappingError(
            "Framework C-alpha alignment requires at least three mapped coordinates; "
            f"found {len(fitted)}"
        )
    reference = np.asarray([reference_coordinates[index] for index in fitted], dtype=float)
    mobile = np.asarray([mobile_coordinates[index] for index in fitted], dtype=float)
    if reference.shape != mobile.shape or reference.shape[1:] != (3,):
        raise ResidueMappingError("Framework coordinate arrays must both have shape (N, 3)")
    if not np.isfinite(reference).all() or not np.isfinite(mobile).all():
        raise ResidueMappingError("Framework coordinates contain non-finite values")

    reference_center = reference.mean(axis=0)
    mobile_center = mobile.mean(axis=0)
    reference_centered = reference - reference_center
    mobile_centered = mobile - mobile_center
    if np.linalg.matrix_rank(reference_centered) < 2 or np.linalg.matrix_rank(
        mobile_centered
    ) < 2:
        raise ResidueMappingError(
            "Framework coordinates are collinear and do not define a unique 3D fit"
        )
    covariance = mobile_centered.T @ reference_centered
    left, _, right_transpose = np.linalg.svd(covariance)
    rotation = right_transpose.T @ left.T
    if np.linalg.det(rotation) < 0:
        right_transpose[-1, :] *= -1
        rotation = right_transpose.T @ left.T
    translation = reference_center - rotation @ mobile_center
    aligned = (rotation @ mobile.T).T + translation
    rmsd = math.sqrt(float(np.mean(np.sum((aligned - reference) ** 2, axis=1))))
    return RigidAlignment(
        rotation=tuple(tuple(float(value) for value in row) for row in rotation),
        translation=tuple(float(value) for value in translation),
        rmsd_angstrom=rmsd,
        fitted_atom_count=len(fitted),
        fitted_authoritative_indices_1based=fitted,
    )


def apply_rigid_alignment(
    coordinates: Sequence[Sequence[float]], alignment: RigidAlignment
) -> np.ndarray:
    """Apply a recorded mobile-to-reference transform to an ``N x 3`` array."""

    points = np.asarray(coordinates, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ResidueMappingError("Coordinates must have shape (N, 3)")
    rotation = np.asarray(alignment.rotation, dtype=float)
    translation = np.asarray(alignment.translation, dtype=float)
    return (rotation @ points.T).T + translation


def aligned_ca_displacement_statistics(
    reference_coordinates: Mapping[int, Sequence[float]],
    mobile_coordinates: Mapping[int, Sequence[float]],
    alignment: RigidAlignment,
    *,
    region_by_authoritative_index: Mapping[int, str],
    region_order: Sequence[str] = (
        "FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"
    ),
) -> dict[str, object]:
    """Summarize fitted C-alpha displacement by explicit antibody regions."""

    fitted_mobile: dict[int, np.ndarray] = {}
    for index, coordinate in mobile_coordinates.items():
        point = apply_rigid_alignment([coordinate], alignment)[0]
        fitted_mobile[index] = point
    distances_by_region: dict[str, list[float]] = {
        region: [] for region in region_order
    }
    for index, region in region_by_authoritative_index.items():
        if region not in distances_by_region:
            continue
        if index not in reference_coordinates or index not in fitted_mobile:
            continue
        reference = np.asarray(reference_coordinates[index], dtype=float)
        distance = float(np.linalg.norm(fitted_mobile[index] - reference))
        distances_by_region[region].append(distance)
    per_region = {
        region: _displacement_summary(distances_by_region[region])
        for region in region_order
    }
    framework = [
        distance
        for region, distances in distances_by_region.items()
        if region.startswith("FR")
        for distance in distances
    ]
    cdr = [
        distance
        for region, distances in distances_by_region.items()
        if region.startswith("CDR")
        for distance in distances
    ]
    return {
        "metric": "post-fit paired C-alpha Euclidean displacement in angstrom",
        "per_region": per_region,
        "aggregates": {
            "FR": _displacement_summary(framework),
            "CDR": _displacement_summary(cdr),
        },
    }


def _validate_sequence(sequence: str, *, label: str) -> None:
    if not sequence:
        raise ResidueMappingError(f"{label.capitalize()} is empty")
    invalid = sorted(set(sequence) - STANDARD_AMINO_ACIDS)
    if invalid:
        raise ResidueMappingError(
            f"{label.capitalize()} contains unsupported residues: {invalid}"
        )


def _displacement_summary(distances: Sequence[float]) -> dict[str, object]:
    if not distances:
        return {
            "status": "not_evaluable",
            "count": 0,
            "rmsd_angstrom": None,
            "mean_displacement_angstrom": None,
            "maximum_displacement_angstrom": None,
        }
    values = np.asarray(distances, dtype=float)
    return {
        "status": "evaluable",
        "count": len(distances),
        "rmsd_angstrom": math.sqrt(float(np.mean(values * values))),
        "mean_displacement_angstrom": float(np.mean(values)),
        "maximum_displacement_angstrom": float(np.max(values)),
    }


def _observed_mapping_from_alignment(
    coordinates: Sequence[Sequence[int]],
    *,
    authoritative_sequence: str,
    observed_sequence: str,
) -> tuple[int, ...] | None:
    """Extract a complete identity-consistent mapping from alignment coordinates."""

    coordinate_array = np.asarray(coordinates, dtype=int)
    if coordinate_array.shape[0] != 2 or coordinate_array.shape[1] < 2:
        raise ResidueMappingError("Unexpected Biopython alignment coordinate shape")
    mapping: list[int | None] = [None] * len(observed_sequence)
    for segment in range(coordinate_array.shape[1] - 1):
        reference_start, reference_end = coordinate_array[0, segment : segment + 2]
        observed_start, observed_end = coordinate_array[1, segment : segment + 2]
        reference_span = int(reference_end - reference_start)
        observed_span = int(observed_end - observed_start)
        if reference_span and observed_span:
            if reference_span != observed_span:
                raise ResidueMappingError("Unexpected unequal diagonal alignment segment")
            for offset in range(reference_span):
                reference_index = int(reference_start + offset)
                observed_index = int(observed_start + offset)
                if observed_sequence[observed_index] != authoritative_sequence[reference_index]:
                    return None
                mapping[observed_index] = reference_index + 1
        elif observed_span:
            return None
    if any(index is None for index in mapping):
        return None
    return tuple(int(index) for index in mapping)
