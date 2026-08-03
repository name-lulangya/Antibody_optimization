"""Read source-declared polymer evidence directly from raw mmCIF categories.

This module never calls Gemmi entity or label-sequence inference.  Values
classified as source evidence come only from categories present in the input
file; partial or contradictory source annotations raise instead of falling
back to inferred metadata.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import gemmi


SOURCE_POLYMER_SEQUENCE_SOURCES = frozenset(
    {
        "source_mmcif_entity_poly_seq",
        "source_mmcif_entity_poly.pdbx_seq_one_letter_code_can",
        "source_mmcif_pdbx_poly_seq_scheme",
        "absent",
    }
)


class PolymerMappingError(ValueError):
    """Raised when source polymer evidence is malformed or inconsistent."""


@dataclass(frozen=True)
class PolySeqSchemeRow:
    """One source ``_pdbx_poly_seq_scheme`` row.

    ``auth_seq_id`` deliberately remains text.  Author numbering can contain
    conventions that are not sequence indices; the insertion code is retained
    as a separate part of the reversible key.
    """

    label_asym_id: str
    entity_id: str
    sequence_index_1based: int
    monomer_id: str
    residue_aa: str
    auth_asym_id: str
    auth_seq_id: str
    insertion_code: str
    pdb_seq_id: str
    pdb_monomer_id: str
    auth_monomer_id: str

    def as_dict(self) -> dict[str, object]:
        """Return a JSON-serializable provenance record."""

        return {
            "label_asym_id": self.label_asym_id,
            "entity_id": self.entity_id,
            "sequence_index_1based": self.sequence_index_1based,
            "monomer_id": self.monomer_id,
            "residue_aa": self.residue_aa,
            "auth_asym_id": self.auth_asym_id,
            "auth_seq_id": self.auth_seq_id,
            "insertion_code": self.insertion_code,
            "pdb_seq_id": self.pdb_seq_id,
            "pdb_monomer_id": self.pdb_monomer_id,
            "auth_monomer_id": self.auth_monomer_id,
        }


@dataclass(frozen=True)
class SourcePolymerEvidence:
    """Source-declared polymer identity for one label asym ID."""

    source_path: str
    source_block_name: str
    label_asym_id: str
    entity_id: str
    entity_description: str
    entity_type: str
    polymer_type: str
    polymer_sequence: str
    polymer_sequence_source: str
    struct_asym_source: str
    scheme_source: str
    scheme_rows: tuple[PolySeqSchemeRow, ...]

    @property
    def polymer_sequence_sha256(self) -> str:
        """Return the sequence digest, or an empty string when it is absent."""

        if not self.polymer_sequence:
            return ""
        return hashlib.sha256(self.polymer_sequence.encode("ascii")).hexdigest()

    def as_dict(self) -> dict[str, object]:
        """Return all source evidence in a JSON-serializable form."""

        return {
            "source_path": self.source_path,
            "source_block_name": self.source_block_name,
            "label_asym_id": self.label_asym_id,
            "entity_id": self.entity_id,
            "entity_description": self.entity_description,
            "entity_type": self.entity_type,
            "polymer_type": self.polymer_type,
            "polymer_sequence": self.polymer_sequence,
            "polymer_sequence_sha256": self.polymer_sequence_sha256,
            "polymer_sequence_source": self.polymer_sequence_source,
            "struct_asym_source": self.struct_asym_source,
            "scheme_source": self.scheme_source,
            "scheme_rows": [row.as_dict() for row in self.scheme_rows],
        }


def read_source_polymer_evidence(
    path: Path,
    *,
    label_asym_id: str,
) -> SourcePolymerEvidence | None:
    """Read source-declared polymer evidence for one mmCIF label asym ID.

    ``None`` is returned only when the file contains none of the relevant
    source categories.  A partially present or contradictory annotation raises
    :class:`PolymerMappingError` so it cannot be silently replaced with
    Gemmi-generated metadata.
    """

    path = Path(path)
    if not path.is_file() or path.is_symlink():
        raise PolymerMappingError(f"mmCIF must be a regular non-symlink file: {path}")
    label_asym_id = _required_text(label_asym_id, "requested label_asym_id")
    try:
        document = gemmi.cif.read_file(str(path))
    except (RuntimeError, OSError) as exc:
        raise PolymerMappingError(f"Gemmi could not read mmCIF {path}: {exc}") from exc
    if len(document) != 1:
        raise PolymerMappingError(
            f"Expected exactly one mmCIF data block in {path}, found {len(document)}"
        )
    block = document.sole_block()

    entity_rows = _category_rows(block, "_entity.")
    entity_poly_rows = _category_rows(block, "_entity_poly.")
    entity_poly_seq_rows = _category_rows(block, "_entity_poly_seq.")
    struct_asym_rows = _category_rows(block, "_struct_asym.")
    scheme_rows_raw = _category_rows(block, "_pdbx_poly_seq_scheme.")
    relevant_categories_present = any(
        (entity_rows, entity_poly_rows, entity_poly_seq_rows, struct_asym_rows, scheme_rows_raw)
    )
    if not relevant_categories_present:
        return None

    descriptions, entity_types = _entity_metadata(entity_rows)
    struct_by_label = _struct_asym_entities(struct_asym_rows)
    matching_scheme_raw = [
        row
        for row in scheme_rows_raw
        if _cif_text(row.get("asym_id", "")) == label_asym_id
    ]
    scheme_entity_ids = {
        _required_text(row.get("entity_id", ""), "_pdbx_poly_seq_scheme.entity_id")
        for row in matching_scheme_raw
    }
    if len(scheme_entity_ids) > 1:
        raise PolymerMappingError(
            f"Scheme rows for label asym {label_asym_id!r} reference multiple entities: "
            f"{sorted(scheme_entity_ids)!r}"
        )

    struct_entity_id = struct_by_label.get(label_asym_id, "")
    scheme_entity_id = next(iter(scheme_entity_ids), "")
    if struct_entity_id and scheme_entity_id and struct_entity_id != scheme_entity_id:
        raise PolymerMappingError(
            "_struct_asym and _pdbx_poly_seq_scheme disagree for label asym "
            f"{label_asym_id!r}: {struct_entity_id!r} vs {scheme_entity_id!r}"
        )
    entity_id = struct_entity_id or scheme_entity_id
    if not entity_id:
        raise PolymerMappingError(
            f"Relevant mmCIF source categories are present, but label asym "
            f"{label_asym_id!r} has no source entity association"
        )

    entity_type = entity_types.get(entity_id, "")
    if entity_type and entity_type.lower() != "polymer":
        raise PolymerMappingError(
            f"Selected entity {entity_id!r} is source-classified as {entity_type!r}, not polymer"
        )

    sequence_from_poly_seq = _entity_poly_seq_sequence(
        entity_poly_seq_rows, entity_id=entity_id
    )
    sequence_from_entity_poly, polymer_type = _entity_poly_sequence(
        entity_poly_rows, entity_id=entity_id
    )
    if (
        sequence_from_poly_seq
        and sequence_from_entity_poly
        and sequence_from_poly_seq != sequence_from_entity_poly
    ):
        raise PolymerMappingError(
            f"_entity_poly_seq and _entity_poly canonical sequence disagree for entity {entity_id!r}"
        )

    parsed_scheme_rows = _parse_scheme_rows(
        matching_scheme_raw,
        label_asym_id=label_asym_id,
        entity_id=entity_id,
    )
    scheme_sequence = _sequence_from_scheme(parsed_scheme_rows)
    polymer_sequence = sequence_from_poly_seq or sequence_from_entity_poly or scheme_sequence
    if sequence_from_poly_seq:
        sequence_source = "source_mmcif_entity_poly_seq"
    elif sequence_from_entity_poly:
        sequence_source = "source_mmcif_entity_poly.pdbx_seq_one_letter_code_can"
    elif scheme_sequence:
        sequence_source = "source_mmcif_pdbx_poly_seq_scheme"
    else:
        sequence_source = "absent"

    if polymer_sequence and scheme_sequence and polymer_sequence != scheme_sequence:
        raise PolymerMappingError(
            f"Polymer sequence and _pdbx_poly_seq_scheme disagree for label asym {label_asym_id!r}"
        )
    if parsed_scheme_rows and polymer_sequence:
        expected_indices = tuple(range(1, len(polymer_sequence) + 1))
        observed_indices = tuple(row.sequence_index_1based for row in parsed_scheme_rows)
        if observed_indices != expected_indices:
            raise PolymerMappingError(
                "_pdbx_poly_seq_scheme does not cover the complete source polymer "
                f"sequence for {label_asym_id!r}: expected {expected_indices!r}, "
                f"observed {observed_indices!r}"
            )

    return SourcePolymerEvidence(
        source_path=str(path),
        source_block_name=str(block.name),
        label_asym_id=label_asym_id,
        entity_id=entity_id,
        entity_description=descriptions.get(entity_id, ""),
        entity_type=entity_type,
        polymer_type=polymer_type,
        polymer_sequence=polymer_sequence,
        polymer_sequence_source=sequence_source,
        struct_asym_source=("source_mmcif_struct_asym" if struct_entity_id else "absent"),
        scheme_source=(
            "source_mmcif_pdbx_poly_seq_scheme" if parsed_scheme_rows else "absent"
        ),
        scheme_rows=parsed_scheme_rows,
    )


def _category_rows(block: gemmi.cif.Block, prefix: str) -> list[dict[str, str]]:
    table = block.find_mmcif_category(prefix)
    if not table:
        return []
    column_names = [str(tag)[len(prefix) :] for tag in table.tags]
    return [
        {
            name: "" if gemmi.cif.as_string(value) is None else gemmi.cif.as_string(value)
            for name, value in zip(column_names, row, strict=True)
        }
        for row in table
    ]


def _entity_metadata(
    rows: Sequence[Mapping[str, str]],
) -> tuple[dict[str, str], dict[str, str]]:
    descriptions: dict[str, str] = {}
    entity_types: dict[str, str] = {}
    for row in rows:
        entity_id = _required_text(row.get("id", ""), "_entity.id")
        if entity_id in descriptions:
            raise PolymerMappingError(f"Duplicate _entity row for entity {entity_id!r}")
        descriptions[entity_id] = _cif_text(row.get("pdbx_description", ""))
        entity_types[entity_id] = _cif_text(row.get("type", ""))
    return descriptions, entity_types


def _struct_asym_entities(rows: Sequence[Mapping[str, str]]) -> dict[str, str]:
    result: dict[str, str] = {}
    for row in rows:
        label_asym_id = _required_text(row.get("id", ""), "_struct_asym.id")
        entity_id = _required_text(row.get("entity_id", ""), "_struct_asym.entity_id")
        if label_asym_id in result:
            raise PolymerMappingError(
                f"Duplicate _struct_asym row for label asym {label_asym_id!r}"
            )
        result[label_asym_id] = entity_id
    return result


def _entity_poly_seq_sequence(
    rows: Sequence[Mapping[str, str]],
    *,
    entity_id: str,
) -> str:
    selected = [
        row for row in rows if _cif_text(row.get("entity_id", "")) == entity_id
    ]
    if not selected:
        return ""
    numbered: dict[int, tuple[str, str]] = {}
    for row in selected:
        index = _positive_int(row.get("num", ""), "_entity_poly_seq.num")
        monomer_id = _required_text(row.get("mon_id", ""), "_entity_poly_seq.mon_id")
        residue_aa = _monomer_to_aa(monomer_id)
        prior = numbered.get(index)
        if prior is not None and prior != (monomer_id, residue_aa):
            raise PolymerMappingError(
                f"Microheterogeneous _entity_poly_seq position {index} for entity {entity_id!r} "
                "does not define one deterministic sequence"
            )
        if prior is not None:
            raise PolymerMappingError(
                f"Duplicate _entity_poly_seq position {index} for entity {entity_id!r}"
            )
        numbered[index] = (monomer_id, residue_aa)
    expected = list(range(1, len(numbered) + 1))
    if sorted(numbered) != expected:
        raise PolymerMappingError(
            f"_entity_poly_seq numbering for entity {entity_id!r} is not contiguous from 1"
        )
    return "".join(numbered[index][1] for index in expected)


def _entity_poly_sequence(
    rows: Sequence[Mapping[str, str]],
    *,
    entity_id: str,
) -> tuple[str, str]:
    selected = [
        row for row in rows if _cif_text(row.get("entity_id", "")) == entity_id
    ]
    if len(selected) > 1:
        raise PolymerMappingError(f"Duplicate _entity_poly row for entity {entity_id!r}")
    if not selected:
        return "", ""
    row = selected[0]
    raw_sequence = _cif_text(row.get("pdbx_seq_one_letter_code_can", ""))
    sequence = "".join(raw_sequence.split()).upper()
    if sequence and not all(character.isalpha() for character in sequence):
        raise PolymerMappingError(
            f"Invalid canonical one-letter polymer sequence for entity {entity_id!r}"
        )
    return sequence, _cif_text(row.get("type", ""))


def _parse_scheme_rows(
    rows: Sequence[Mapping[str, str]],
    *,
    label_asym_id: str,
    entity_id: str,
) -> tuple[PolySeqSchemeRow, ...]:
    parsed: dict[int, PolySeqSchemeRow] = {}
    for row in rows:
        row_entity = _required_text(
            row.get("entity_id", ""), "_pdbx_poly_seq_scheme.entity_id"
        )
        if row_entity != entity_id:
            raise PolymerMappingError(
                f"Scheme entity changed within label asym {label_asym_id!r}"
            )
        index = _positive_int(
            row.get("seq_id", ""), "_pdbx_poly_seq_scheme.seq_id"
        )
        if index in parsed:
            raise PolymerMappingError(
                f"Duplicate _pdbx_poly_seq_scheme seq_id {index} for {label_asym_id!r}"
            )
        monomer_id = _required_text(
            row.get("mon_id", ""), "_pdbx_poly_seq_scheme.mon_id"
        )
        parsed[index] = PolySeqSchemeRow(
            label_asym_id=label_asym_id,
            entity_id=entity_id,
            sequence_index_1based=index,
            monomer_id=monomer_id,
            residue_aa=_monomer_to_aa(monomer_id),
            auth_asym_id=_cif_text(row.get("pdb_strand_id", "")),
            auth_seq_id=_cif_text(row.get("auth_seq_num", "")),
            insertion_code=_cif_text(row.get("pdb_ins_code", "")),
            pdb_seq_id=_cif_text(row.get("pdb_seq_num", "")),
            pdb_monomer_id=_cif_text(row.get("pdb_mon_id", "")),
            auth_monomer_id=_cif_text(row.get("auth_mon_id", "")),
        )
    return tuple(parsed[index] for index in sorted(parsed))


def _sequence_from_scheme(rows: Sequence[PolySeqSchemeRow]) -> str:
    if not rows:
        return ""
    expected = tuple(range(1, len(rows) + 1))
    observed = tuple(row.sequence_index_1based for row in rows)
    if observed != expected:
        raise PolymerMappingError(
            "_pdbx_poly_seq_scheme sequence indices are not contiguous from 1"
        )
    return "".join(row.residue_aa for row in rows)


def _monomer_to_aa(monomer_id: str) -> str:
    info = gemmi.find_tabulated_residue(monomer_id)
    code = str(info.one_letter_code).strip().upper()
    return code if len(code) == 1 and code.isalpha() else "X"


def _positive_int(value: object, label: str) -> int:
    text = _required_text(value, label)
    try:
        number = int(text)
    except ValueError as exc:
        raise PolymerMappingError(f"{label} is not an integer: {text!r}") from exc
    if number < 1:
        raise PolymerMappingError(f"{label} must be positive: {number}")
    return number


def _required_text(value: object, label: str) -> str:
    text = _cif_text(value)
    if not text:
        raise PolymerMappingError(f"Missing required mmCIF value {label}")
    return text


def _cif_text(value: object) -> str:
    text = "" if value is None else str(value).strip()
    return "" if text in {".", "?"} else text


__all__ = [
    "PolymerMappingError",
    "PolySeqSchemeRow",
    "SOURCE_POLYMER_SEQUENCE_SOURCES",
    "SourcePolymerEvidence",
    "read_source_polymer_evidence",
]
