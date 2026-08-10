"""Source-aware mmCIF polymer-to-authoritative residue mapping.

The source mmCIF polymer sequence and its label/auth numbering scheme are
stronger evidence than the coordinate-observed sequence.  This module keeps
those evidence layers separate and permits a sequence-alignment fallback only
when the corresponding source evidence is absent.  A present but inconsistent
source annotation is always an error.

No Gemmi entity or label-sequence inference is performed here. Raw category
parsing is isolated in ``source_mmcif_evidence`` and re-exported by this module.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from .residue_mapping import ResidueMappingError, exact_first_sequence_mapping
from .source_mmcif_evidence import (
    SOURCE_POLYMER_SEQUENCE_SOURCES,
    PolySeqSchemeRow,
    PolymerMappingError,
    SourcePolymerEvidence,
    read_source_polymer_evidence,
)
from .structure_inventory import StructureResidue


POLYMER_MAPPING_VERSION = "1.1.0"
STANDARD_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
LABEL_SEQ_ID_SOURCES = frozenset(
    {"source_mmcif_atom_site", "gemmi_heuristic", "absent"}
)
@dataclass(frozen=True)
class ObservedResidueEvidence:
    """One coordinate-observed residue with explicit identifier provenance."""

    residue_aa: str
    auth_asym_id: str = ""
    auth_seq_id: str = ""
    insertion_code: str = ""
    label_seq_id: int | None = None
    label_seq_id_source: str = "absent"


@dataclass(frozen=True)
class PolymerMappingResult:
    """Deterministic composed observed-to-authoritative mapping."""

    authoritative_index_1based_by_observed_index: tuple[int, ...]
    polymer_index_1based_by_observed_index: tuple[int, ...] | None
    polymer_sequence_source: str
    polymer_sequence_sha256: str
    entity_id: str
    entity_description: str
    label_seq_id_source: str
    polymer_to_authoritative_method: str
    observed_to_polymer_method: str
    observed_to_authoritative_method: str
    mapping_status: str
    fallback_reason: str
    source_path: str
    source_block_name: str

    @property
    def method(self) -> str:
        """Compatibility-friendly name for the final composed method."""

        return self.observed_to_authoritative_method

    def as_dict(self) -> dict[str, object]:
        """Return compact JSON-serializable method and provenance metadata."""

        return {
            "authoritative_index_1based_by_observed_index": list(
                self.authoritative_index_1based_by_observed_index
            ),
            "polymer_index_1based_by_observed_index": (
                None
                if self.polymer_index_1based_by_observed_index is None
                else list(self.polymer_index_1based_by_observed_index)
            ),
            "polymer_sequence_source": self.polymer_sequence_source,
            "polymer_sequence_sha256": self.polymer_sequence_sha256,
            "entity_id": self.entity_id,
            "entity_description": self.entity_description,
            "label_seq_id_source": self.label_seq_id_source,
            "polymer_to_authoritative_method": self.polymer_to_authoritative_method,
            "observed_to_polymer_method": self.observed_to_polymer_method,
            "observed_to_authoritative_method": self.observed_to_authoritative_method,
            "mapping_status": self.mapping_status,
            "fallback_reason": self.fallback_reason,
            "source_path": self.source_path,
            "source_block_name": self.source_block_name,
        }


def observations_from_structure_residues(
    residues: Sequence[StructureResidue],
    *,
    label_seq_id_source: str,
) -> tuple[ObservedResidueEvidence, ...]:
    """Copy observed structure residues while explicitly classifying label IDs."""

    if label_seq_id_source not in LABEL_SEQ_ID_SOURCES:
        raise PolymerMappingError(
            f"Unsupported label_seq_id_source: {label_seq_id_source!r}"
        )
    observations = tuple(
        ObservedResidueEvidence(
            residue_aa=residue.one_letter_code,
            auth_asym_id=residue.key.auth_asym_id,
            auth_seq_id=str(residue.key.auth_seq_id),
            insertion_code=residue.key.insertion_code,
            label_seq_id=residue.key.label_seq_id,
            label_seq_id_source=label_seq_id_source,
        )
        for residue in residues
    )
    _validate_observations(observations)
    return observations


def compose_polymer_mapping(
    *,
    authoritative_sequence: str,
    observed_residues: Sequence[ObservedResidueEvidence],
    source_evidence: SourcePolymerEvidence | None,
) -> PolymerMappingResult:
    """Map observed residues through source polymer evidence when available.

    Priority is fixed:

    1. source polymer sequence to authoritative sequence;
    2. source ``label_seq_id`` or source auth/insertion scheme to polymer index;
    3. observed-sequence alignment to the source polymer only if both identifier
       routes are absent;
    4. when the source polymer sequence is absent, exact integer source auth
       numbering may map directly to the authoritative sequence, but only if
       every observed WT residue agrees at that index;
    5. otherwise, observed-sequence alignment directly to the authoritative
       sequence.

    Lower-priority evidence never overrides a present conflicting source value.
    """

    _validate_standard_sequence(authoritative_sequence, "authoritative sequence")
    observations = tuple(observed_residues)
    _validate_observations(observations)
    _validate_source_evidence(source_evidence)
    observed_sequence = "".join(row.residue_aa for row in observations)
    label_source = observations[0].label_seq_id_source

    if source_evidence is None or not source_evidence.polymer_sequence:
        if source_evidence is not None and source_evidence.scheme_rows:
            raise PolymerMappingError(
                "Source scheme rows are present but no source polymer sequence can be established"
            )
        reason = (
            "all_source_polymer_categories_absent"
            if source_evidence is None
            else "source_polymer_sequence_absent"
        )
        auth_direct = _mapping_from_direct_auth_indices_if_consistent(
            observations, authoritative_sequence
        )
        if auth_direct is None:
            direct = _exact_mapping(
                authoritative_sequence,
                observed_sequence,
                context="observed sequence to authoritative fallback",
            )
            authoritative_indices = (
                direct.authoritative_index_1based_by_observed_index
            )
            direct_method = f"observed_sequence_to_authoritative_fallback:{direct.method}"
            mapping_status = "observed_sequence_only_fallback"
        else:
            authoritative_indices = auth_direct
            direct_method = "source_atom_site.auth_seq_id_direct_exact_wt"
            mapping_status = "source_auth_numbering_direct_exact_wt"
        return PolymerMappingResult(
            authoritative_index_1based_by_observed_index=authoritative_indices,
            polymer_index_1based_by_observed_index=None,
            polymer_sequence_source=(
                "absent" if source_evidence is None else source_evidence.polymer_sequence_source
            ),
            polymer_sequence_sha256="",
            entity_id="" if source_evidence is None else source_evidence.entity_id,
            entity_description=(
                "" if source_evidence is None else source_evidence.entity_description
            ),
            label_seq_id_source=label_source,
            polymer_to_authoritative_method="not_applicable_no_source_polymer_sequence",
            observed_to_polymer_method="not_applicable_no_source_polymer_sequence",
            observed_to_authoritative_method=direct_method,
            mapping_status=mapping_status,
            fallback_reason=reason,
            source_path="" if source_evidence is None else source_evidence.source_path,
            source_block_name=(
                "" if source_evidence is None else source_evidence.source_block_name
            ),
        )

    polymer_sequence = source_evidence.polymer_sequence
    _validate_standard_sequence(polymer_sequence, "source mmCIF polymer sequence")
    source_label_mapping: tuple[int, ...] | None = None
    if label_source == "source_mmcif_atom_site":
        source_label_mapping = _mapping_from_label_seq_ids(
            observations, polymer_sequence, evidence_label="source label_seq_id"
        )

    scheme_mapping = _mapping_from_auth_scheme_if_available(
        observations, source_evidence, polymer_sequence
    )
    if source_label_mapping is not None:
        if scheme_mapping is not None and source_label_mapping != scheme_mapping:
            raise PolymerMappingError(
                "Source label_seq_id and source auth/insertion scheme imply different "
                "observed-to-polymer mappings"
            )
        observed_to_polymer = source_label_mapping
        observed_method = "source_mmcif_atom_site.label_seq_id"
        fallback_reason = ""
        mapping_status = "source_mmcif_evidence_consistent"
    elif scheme_mapping is not None:
        observed_to_polymer = scheme_mapping
        observed_method = (
            "source_mmcif_pdbx_poly_seq_scheme.auth_seq_num+pdb_ins_code"
        )
        fallback_reason = ""
        mapping_status = "source_mmcif_evidence_consistent"
    else:
        fallback = _exact_mapping(
            polymer_sequence,
            observed_sequence,
            context="observed sequence to source polymer fallback",
        )
        observed_to_polymer = fallback.authoritative_index_1based_by_observed_index
        observed_method = f"observed_sequence_to_polymer_fallback:{fallback.method}"
        fallback_reason = "source_label_seq_id_and_auth_scheme_absent"
        mapping_status = "source_polymer_with_observed_sequence_fallback"

    if label_source == "gemmi_heuristic":
        heuristic_mapping = _mapping_from_label_seq_ids(
            observations, polymer_sequence, evidence_label="Gemmi heuristic label_seq_id"
        )
        if heuristic_mapping != observed_to_polymer:
            raise PolymerMappingError(
                "Gemmi heuristic label_seq_id conflicts with the source-aware mapping"
            )

    authoritative_indices, polymer_to_authoritative_method = (
        _authoritative_indices_for_observed_polymer_positions(
            authoritative_sequence,
            polymer_sequence,
            observed_to_polymer,
        )
    )
    final_method = (
        f"composed:{observed_method}+polymer_to_authoritative:"
        f"{polymer_to_authoritative_method}"
    )
    return PolymerMappingResult(
        authoritative_index_1based_by_observed_index=authoritative_indices,
        polymer_index_1based_by_observed_index=observed_to_polymer,
        polymer_sequence_source=source_evidence.polymer_sequence_source,
        polymer_sequence_sha256=source_evidence.polymer_sequence_sha256,
        entity_id=source_evidence.entity_id,
        entity_description=source_evidence.entity_description,
        label_seq_id_source=label_source,
        polymer_to_authoritative_method=polymer_to_authoritative_method,
        observed_to_polymer_method=observed_method,
        observed_to_authoritative_method=final_method,
        mapping_status=mapping_status,
        fallback_reason=fallback_reason,
        source_path=source_evidence.source_path,
        source_block_name=source_evidence.source_block_name,
    )


def _validate_observations(rows: Sequence[ObservedResidueEvidence]) -> None:
    if not rows:
        raise PolymerMappingError("Observed residue evidence is empty")
    label_sources = {row.label_seq_id_source for row in rows}
    unknown = label_sources - LABEL_SEQ_ID_SOURCES
    if unknown:
        raise PolymerMappingError(f"Unsupported label_seq_id source(s): {sorted(unknown)!r}")
    if len(label_sources) != 1:
        raise PolymerMappingError(
            f"Observed residues mix label_seq_id provenance: {sorted(label_sources)!r}"
        )
    for index, row in enumerate(rows, start=1):
        _validate_standard_sequence(row.residue_aa, f"observed residue {index}")
        if len(row.residue_aa) != 1:
            raise PolymerMappingError(
                f"Observed residue {index} must contain exactly one amino-acid letter"
            )
        if row.label_seq_id_source == "absent" and row.label_seq_id is not None:
            raise PolymerMappingError(
                f"Observed residue {index} has label_seq_id but marks its source absent"
            )
        if row.label_seq_id_source != "absent" and row.label_seq_id is None:
            raise PolymerMappingError(
                f"Observed residue {index} lacks label_seq_id despite source "
                f"{row.label_seq_id_source!r}"
            )


def _validate_source_evidence(evidence: SourcePolymerEvidence | None) -> None:
    if evidence is None:
        return
    if evidence.polymer_sequence_source not in SOURCE_POLYMER_SEQUENCE_SOURCES:
        raise PolymerMappingError(
            "Polymer sequence is not classified as a raw source mmCIF value: "
            f"{evidence.polymer_sequence_source!r}"
        )
    if evidence.polymer_sequence and evidence.polymer_sequence_source == "absent":
        raise PolymerMappingError(
            "A polymer sequence is present but its source is classified as absent"
        )
    if not evidence.polymer_sequence and evidence.polymer_sequence_source != "absent":
        raise PolymerMappingError(
            "Polymer sequence source is declared but the sequence is empty"
        )
    if evidence.struct_asym_source not in {"source_mmcif_struct_asym", "absent"}:
        raise PolymerMappingError(
            f"Unsupported struct_asym provenance: {evidence.struct_asym_source!r}"
        )
    if evidence.scheme_source not in {
        "source_mmcif_pdbx_poly_seq_scheme",
        "absent",
    }:
        raise PolymerMappingError(
            f"Unsupported poly-seq scheme provenance: {evidence.scheme_source!r}"
        )


def _mapping_from_label_seq_ids(
    observations: Sequence[ObservedResidueEvidence],
    polymer_sequence: str,
    *,
    evidence_label: str,
) -> tuple[int, ...]:
    indices: list[int] = []
    for row in observations:
        value = row.label_seq_id
        if isinstance(value, bool) or not isinstance(value, int):
            raise PolymerMappingError(f"{evidence_label} is not an integer: {value!r}")
        indices.append(value)
    return _validate_polymer_indices(
        tuple(indices), observations, polymer_sequence, evidence_label=evidence_label
    )


def _mapping_from_direct_auth_indices_if_consistent(
    observations: Sequence[ObservedResidueEvidence],
    authoritative_sequence: str,
) -> tuple[int, ...] | None:
    """Use source auth numbers as sequence indices only under exact WT agreement."""

    indices: list[int] = []
    for row in observations:
        if row.insertion_code or not row.auth_seq_id.isdecimal():
            return None
        index = int(row.auth_seq_id)
        if not 1 <= index <= len(authoritative_sequence):
            return None
        if authoritative_sequence[index - 1] != row.residue_aa:
            return None
        indices.append(index)
    if len(set(indices)) != len(indices):
        return None
    if any(left >= right for left, right in zip(indices, indices[1:])):
        return None
    return tuple(indices)


def _authoritative_indices_for_observed_polymer_positions(
    authoritative_sequence: str,
    polymer_sequence: str,
    observed_to_polymer: tuple[int, ...],
) -> tuple[tuple[int, ...], str]:
    """Map observed polymer indices into one exact authoritative segment.

    A longer source polymer is accepted only when it contains the complete
    authoritative sequence exactly once and every coordinate-observed residue
    lies inside that segment. This preserves source flanks without treating
    them as part of the reviewed VHH.
    """

    if len(polymer_sequence) > len(authoritative_sequence):
        starts = tuple(
            index
            for index in range(len(polymer_sequence) - len(authoritative_sequence) + 1)
            if polymer_sequence.startswith(authoritative_sequence, index)
        )
        if len(starts) == 1:
            start_0based = starts[0]
            start_1based = start_0based + 1
            end_1based = start_0based + len(authoritative_sequence)
            if any(
                index < start_1based or index > end_1based
                for index in observed_to_polymer
            ):
                raise PolymerMappingError(
                    "Observed residues extend outside the unique exact authoritative "
                    "segment in the source polymer sequence"
                )
            return (
                tuple(index - start_0based for index in observed_to_polymer),
                f"unique_exact_authoritative_segment_in_source_polymer:start={start_1based}",
            )

    polymer_to_authoritative = _exact_mapping(
        authoritative_sequence,
        polymer_sequence,
        context="source polymer sequence to authoritative sequence",
    )
    return (
        tuple(
            polymer_to_authoritative.authoritative_index_1based_by_observed_index[
                index - 1
            ]
            for index in observed_to_polymer
        ),
        polymer_to_authoritative.method,
    )


def _mapping_from_auth_scheme_if_available(
    observations: Sequence[ObservedResidueEvidence],
    evidence: SourcePolymerEvidence,
    polymer_sequence: str,
) -> tuple[int, ...] | None:
    usable_rows = [
        row for row in evidence.scheme_rows if row.auth_asym_id and row.auth_seq_id
    ]
    if not usable_rows:
        return None
    by_key: dict[tuple[str, str, str], int] = {}
    for row in usable_rows:
        key = (row.auth_asym_id, row.auth_seq_id, row.insertion_code)
        if key in by_key:
            raise PolymerMappingError(
                "Duplicate source auth/insertion key in _pdbx_poly_seq_scheme: "
                f"{key!r}"
            )
        by_key[key] = row.sequence_index_1based
    indices: list[int] = []
    for row in observations:
        key = (
            _clean_identifier(row.auth_asym_id),
            _clean_identifier(str(row.auth_seq_id)),
            _clean_identifier(row.insertion_code),
        )
        if not key[0] or not key[1]:
            raise PolymerMappingError(
                "Source auth scheme is available, but an observed residue lacks an "
                "auth_asym_id/auth_seq_id key"
            )
        if key not in by_key:
            raise PolymerMappingError(
                f"Observed auth/insertion key {key!r} is absent from the source scheme"
            )
        indices.append(by_key[key])
    return _validate_polymer_indices(
        tuple(indices),
        observations,
        polymer_sequence,
        evidence_label="source auth/insertion scheme",
    )


def _validate_polymer_indices(
    indices: tuple[int, ...],
    observations: Sequence[ObservedResidueEvidence],
    polymer_sequence: str,
    *,
    evidence_label: str,
) -> tuple[int, ...]:
    if len(indices) != len(observations):
        raise PolymerMappingError(f"{evidence_label} mapping length mismatch")
    if any(index < 1 or index > len(polymer_sequence) for index in indices):
        raise PolymerMappingError(
            f"{evidence_label} contains an out-of-range polymer sequence index: {indices!r}"
        )
    if len(set(indices)) != len(indices):
        raise PolymerMappingError(f"{evidence_label} contains duplicate sequence indices")
    if any(left >= right for left, right in zip(indices, indices[1:])):
        raise PolymerMappingError(
            f"{evidence_label} is not strictly increasing in observed residue order"
        )
    for observation, index in zip(observations, indices, strict=True):
        expected = polymer_sequence[index - 1]
        if observation.residue_aa != expected:
            raise PolymerMappingError(
                f"{evidence_label} wild-type mismatch at polymer index {index}: "
                f"observed {observation.residue_aa!r}, source {expected!r}"
            )
    return indices


def _exact_mapping(reference: str, query: str, *, context: str):
    try:
        return exact_first_sequence_mapping(reference, query)
    except ResidueMappingError as exc:
        raise PolymerMappingError(f"{context} failed: {exc}") from exc


def _validate_standard_sequence(sequence: str, label: str) -> None:
    if not sequence:
        raise PolymerMappingError(f"{label} is empty")
    invalid = sorted(set(sequence) - STANDARD_AMINO_ACIDS)
    if invalid:
        raise PolymerMappingError(
            f"{label} contains unsupported residue code(s): {invalid!r}"
        )


def _clean_identifier(value: object) -> str:
    text = "" if value is None else str(value).strip()
    return "" if text in {".", "?"} else text


__all__ = [
    "ObservedResidueEvidence",
    "POLYMER_MAPPING_VERSION",
    "PolySeqSchemeRow",
    "PolymerMappingError",
    "PolymerMappingResult",
    "SourcePolymerEvidence",
    "compose_polymer_mapping",
    "observations_from_structure_residues",
    "read_source_polymer_evidence",
]
