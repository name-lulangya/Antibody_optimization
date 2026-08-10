from __future__ import annotations

import sys
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC = PROJECT_ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from antibody_optimization.polymer_mapping import (
    ObservedResidueEvidence,
    PolymerMappingError,
    compose_polymer_mapping,
    read_source_polymer_evidence,
)


ONE_TO_THREE = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "N": "ASN",
}


def write_source_cif(
    path: Path,
    sequence: str,
    *,
    description: str | None = "Nb252 VHH",
    include_scheme: bool = True,
    auth_ids: list[str] | None = None,
    insertion_codes: list[str] | None = None,
) -> None:
    auth_ids = auth_ids or [str(index) for index in range(1, len(sequence) + 1)]
    insertion_codes = insertion_codes or ["" for _ in sequence]
    lines = [
        "data_source",
        "loop_",
        "_entity.id",
        "_entity.type",
        "_entity.pdbx_description",
        f"1 polymer {'?' if description is None else repr(description)}",
        "#",
        "loop_",
        "_entity_poly.entity_id",
        "_entity_poly.type",
        "_entity_poly.pdbx_seq_one_letter_code_can",
        f"1 'polypeptide(L)' {sequence}",
        "#",
        "loop_",
        "_entity_poly_seq.entity_id",
        "_entity_poly_seq.num",
        "_entity_poly_seq.mon_id",
        "_entity_poly_seq.hetero",
    ]
    lines.extend(
        f"1 {index} {ONE_TO_THREE[aa]} n"
        for index, aa in enumerate(sequence, start=1)
    )
    lines.extend(
        [
            "#",
            "loop_",
            "_struct_asym.id",
            "_struct_asym.entity_id",
            "V 1",
            "#",
        ]
    )
    if include_scheme:
        lines.extend(
            [
                "loop_",
                "_pdbx_poly_seq_scheme.asym_id",
                "_pdbx_poly_seq_scheme.entity_id",
                "_pdbx_poly_seq_scheme.seq_id",
                "_pdbx_poly_seq_scheme.mon_id",
                "_pdbx_poly_seq_scheme.pdb_seq_num",
                "_pdbx_poly_seq_scheme.auth_seq_num",
                "_pdbx_poly_seq_scheme.pdb_mon_id",
                "_pdbx_poly_seq_scheme.auth_mon_id",
                "_pdbx_poly_seq_scheme.pdb_strand_id",
                "_pdbx_poly_seq_scheme.pdb_ins_code",
            ]
        )
        for index, (aa, auth_id, insertion) in enumerate(
            zip(sequence, auth_ids, insertion_codes, strict=True), start=1
        ):
            monomer = ONE_TO_THREE[aa]
            lines.append(
                f"V 1 {index} {monomer} {auth_id} {auth_id} {monomer} "
                f"{monomer} H {insertion or '?'}"
            )
        lines.append("#")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def observations(
    sequence: str,
    *,
    labels: list[int] | None = None,
    label_source: str = "absent",
    auth_ids: list[str] | None = None,
    insertion_codes: list[str] | None = None,
) -> tuple[ObservedResidueEvidence, ...]:
    labels_or_none = labels or [None for _ in sequence]
    auth_ids = auth_ids or [str(index) for index in range(1, len(sequence) + 1)]
    insertion_codes = insertion_codes or ["" for _ in sequence]
    return tuple(
        ObservedResidueEvidence(
            residue_aa=aa,
            auth_asym_id="H",
            auth_seq_id=auth_id,
            insertion_code=insertion,
            label_seq_id=label,
            label_seq_id_source=label_source,
        )
        for aa, auth_id, insertion, label in zip(
            sequence, auth_ids, insertion_codes, labels_or_none, strict=True
        )
    )


class SourceAwarePolymerMappingTests(unittest.TestCase):
    def read_evidence(self, path: Path):
        evidence = read_source_polymer_evidence(path, label_asym_id="V")
        self.assertIsNotNone(evidence)
        return evidence

    def test_source_label_seq_disambiguates_repeated_motif(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "source.cif"
            write_source_cif(path, "AAAC", auth_ids=["10", "11", "12", "13"])
            evidence = self.read_evidence(path)
            result = compose_polymer_mapping(
                authoritative_sequence="AAAC",
                observed_residues=observations(
                    "AC",
                    labels=[1, 4],
                    label_source="source_mmcif_atom_site",
                    auth_ids=["10", "13"],
                ),
                source_evidence=evidence,
            )

        self.assertEqual(result.polymer_index_1based_by_observed_index, (1, 4))
        self.assertEqual(result.authoritative_index_1based_by_observed_index, (1, 4))
        self.assertEqual(
            result.observed_to_polymer_method,
            "source_mmcif_atom_site.label_seq_id",
        )
        self.assertEqual(result.fallback_reason, "")

    def test_full_entity_conflict_blocks_even_when_observed_subset_matches(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "source.cif"
            write_source_cif(path, "ACNEFG")
            evidence = self.read_evidence(path)
            with self.assertRaisesRegex(
                PolymerMappingError, "source polymer sequence to authoritative"
            ):
                compose_polymer_mapping(
                    authoritative_sequence="ACDEFG",
                    observed_residues=observations(
                        "ACEFG",
                        labels=[1, 2, 4, 5, 6],
                        label_source="source_mmcif_atom_site",
                        auth_ids=["1", "2", "4", "5", "6"],
                    ),
                    source_evidence=evidence,
                )

    def test_unique_authoritative_segment_allows_extra_source_polymer_flank(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "source.cif"
            write_source_cif(path, "ACDEFG")
            evidence = self.read_evidence(path)
            result = compose_polymer_mapping(
                authoritative_sequence="ACD",
                observed_residues=observations(
                    "ACD",
                    labels=[1, 2, 3],
                    label_source="source_mmcif_atom_site",
                ),
                source_evidence=evidence,
            )

        self.assertEqual(result.authoritative_index_1based_by_observed_index, (1, 2, 3))
        self.assertEqual(
            result.polymer_to_authoritative_method,
            "unique_exact_authoritative_segment_in_source_polymer:start=1",
        )

    def test_internal_missing_coordinates_follow_source_label_positions(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "source.cif"
            write_source_cif(path, "ACDEFG")
            evidence = self.read_evidence(path)
            result = compose_polymer_mapping(
                authoritative_sequence="ACDEFG",
                observed_residues=observations(
                    "ACFG",
                    labels=[1, 2, 5, 6],
                    label_source="source_mmcif_atom_site",
                    auth_ids=["1", "2", "5", "6"],
                ),
                source_evidence=evidence,
            )

        self.assertEqual(result.polymer_index_1based_by_observed_index, (1, 2, 5, 6))
        self.assertEqual(
            result.authoritative_index_1based_by_observed_index, (1, 2, 5, 6)
        )
        self.assertEqual(result.mapping_status, "source_mmcif_evidence_consistent")

    def test_label_seq_anomalies_block_without_alignment_fallback(self) -> None:
        cases = (
            (
                "duplicate",
                "AA",
                [1, 1],
                "duplicate",
            ),
            (
                "out_of_range",
                "AC",
                [1, 4],
                "out-of-range",
            ),
            (
                "nonmonotonic",
                "CA",
                [2, 1],
                "not strictly increasing",
            ),
            (
                "wild_type_mismatch",
                "C",
                [1],
                "wild-type mismatch",
            ),
        )
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "source.cif"
            write_source_cif(path, "AAC", include_scheme=False)
            evidence = self.read_evidence(path)
            for name, observed_sequence, labels, message in cases:
                with self.subTest(name=name):
                    with self.assertRaisesRegex(PolymerMappingError, message):
                        compose_polymer_mapping(
                            authoritative_sequence="AAC",
                            observed_residues=observations(
                                observed_sequence,
                                labels=labels,
                                label_source="source_mmcif_atom_site",
                            ),
                            source_evidence=evidence,
                        )

    def test_auth_scheme_preserves_insertion_and_is_not_a_sequence_index(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "source.cif"
            write_source_cif(
                path,
                "ACD",
                description="Nb252 source entity",
                auth_ids=["100", "100", "101"],
                insertion_codes=["", "A", ""],
            )
            evidence = self.read_evidence(path)
            result = compose_polymer_mapping(
                authoritative_sequence="ACD",
                observed_residues=observations(
                    "ACD",
                    auth_ids=["100", "100", "101"],
                    insertion_codes=["", "A", ""],
                ),
                source_evidence=evidence,
            )

        self.assertEqual(result.polymer_index_1based_by_observed_index, (1, 2, 3))
        self.assertIn("auth_seq_num+pdb_ins_code", result.observed_to_polymer_method)
        self.assertEqual(evidence.entity_description, "Nb252 source entity")
        self.assertEqual(
            [(row.auth_seq_id, row.insertion_code) for row in evidence.scheme_rows],
            [("100", ""), ("100", "A"), ("101", "")],
        )
        self.assertEqual(evidence.as_dict()["entity_description"], "Nb252 source entity")

    def test_source_label_and_scheme_conflict_blocks(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "source.cif"
            write_source_cif(path, "AAC", auth_ids=["10", "11", "12"])
            evidence = self.read_evidence(path)
            with self.assertRaisesRegex(PolymerMappingError, "imply different"):
                compose_polymer_mapping(
                    authoritative_sequence="AAC",
                    observed_residues=observations(
                        "AC",
                        labels=[2, 3],
                        label_source="source_mmcif_atom_site",
                        auth_ids=["10", "12"],
                    ),
                    source_evidence=evidence,
                )

    def test_source_polymer_uses_observed_fallback_only_when_ids_absent(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "source.cif"
            write_source_cif(path, "AAAC", include_scheme=False)
            evidence = self.read_evidence(path)
            result = compose_polymer_mapping(
                authoritative_sequence="AAAC",
                observed_residues=observations("AC"),
                source_evidence=evidence,
            )

        self.assertEqual(result.polymer_index_1based_by_observed_index, (3, 4))
        self.assertTrue(
            result.observed_to_polymer_method.startswith(
                "observed_sequence_to_polymer_fallback:"
            )
        )
        self.assertEqual(
            result.fallback_reason, "source_label_seq_id_and_auth_scheme_absent"
        )

    def test_direct_observed_fallback_requires_all_source_categories_absent(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "empty.cif"
            path.write_text("data_empty\n_entry.id test\n", encoding="utf-8")
            evidence = read_source_polymer_evidence(path, label_asym_id="V")
            self.assertIsNone(evidence)
            result = compose_polymer_mapping(
                authoritative_sequence="AAAC",
                observed_residues=observations("AC"),
                source_evidence=evidence,
            )

        self.assertEqual(result.authoritative_index_1based_by_observed_index, (3, 4))
        self.assertIsNone(result.polymer_index_1based_by_observed_index)
        self.assertEqual(result.mapping_status, "observed_sequence_only_fallback")
        self.assertEqual(result.fallback_reason, "all_source_polymer_categories_absent")

    def test_exact_source_auth_numbers_resolve_repeated_sequence_without_polymer(self) -> None:
        result = compose_polymer_mapping(
            authoritative_sequence="AAAA",
            observed_residues=observations("AA", auth_ids=["1", "4"]),
            source_evidence=None,
        )

        self.assertEqual(result.authoritative_index_1based_by_observed_index, (1, 4))
        self.assertEqual(
            result.observed_to_authoritative_method,
            "source_atom_site.auth_seq_id_direct_exact_wt",
        )
        self.assertEqual(result.mapping_status, "source_auth_numbering_direct_exact_wt")

    def test_missing_entity_description_is_preserved_as_empty(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "source.cif"
            write_source_cif(path, "ACD", description=None)
            evidence = self.read_evidence(path)
        self.assertEqual(evidence.entity_description, "")

    def test_heuristic_label_id_is_not_promoted_to_source_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "source.cif"
            write_source_cif(path, "AAAC", include_scheme=False)
            evidence = self.read_evidence(path)
            with self.assertRaisesRegex(PolymerMappingError, "Gemmi heuristic"):
                compose_polymer_mapping(
                    authoritative_sequence="AAAC",
                    observed_residues=observations(
                        "AC",
                        labels=[1, 4],
                        label_source="gemmi_heuristic",
                    ),
                    source_evidence=evidence,
                )

    def test_heuristic_entity_sequence_is_rejected_as_source(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "source.cif"
            write_source_cif(path, "ACD", include_scheme=False)
            evidence = replace(
                self.read_evidence(path),
                polymer_sequence_source="gemmi_heuristic_entity",
            )
            with self.assertRaisesRegex(PolymerMappingError, "raw source mmCIF"):
                compose_polymer_mapping(
                    authoritative_sequence="ACD",
                    observed_residues=observations("ACD"),
                    source_evidence=evidence,
                )


if __name__ == "__main__":
    unittest.main()
