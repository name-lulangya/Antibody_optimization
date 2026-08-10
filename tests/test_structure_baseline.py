from __future__ import annotations

import importlib.util
import hashlib
import json
import runpy
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import gemmi
import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC = PROJECT_ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from antibody_optimization.baseline_review import (
    BaselineReviewError,
    build_review_template,
    chain_identity_sha256,
    load_authoritative_construct_confirmation,
    load_confirmed_review,
)
from antibody_optimization.interface_contacts import (
    AtomSite,
    ContactPair,
    InterfaceContactError,
    neighbor_search_heavy_atom_contacts,
    strict_heavy_atom_contacts,
)
from antibody_optimization.residue_mapping import (
    ResidueMappingError,
    aligned_ca_displacement_statistics,
    apply_rigid_alignment,
    build_mapping_rows,
    exact_first_sequence_mapping,
    fit_explicit_framework_ca,
    unique_exact_order_mapping,
    validate_numbering_rows,
)
from antibody_optimization.structure_inventory import (
    ChainSelector,
    ResidueKey,
    atom_site_classification_counts,
    chain_inventory_rows,
    extract_confirmed_chain_residues,
    prepare_heuristic_analysis_copy,
    read_single_model_structure,
    residue_inventory_rows,
    rigid_coordinate_relationship,
    structure_count_summary,
    require_matching_topology,
)


def make_structure(
    chains: list[tuple[str, str, list[tuple[str, int, str, tuple[float, float, float]]]]]
) -> gemmi.Structure:
    structure = gemmi.Structure()
    structure.name = "test_only_fixture"
    model = gemmi.Model("1")
    for entity_number, (auth_chain, label_chain, residues) in enumerate(chains, start=1):
        chain = gemmi.Chain(auth_chain)
        sequence: list[str] = []
        for label_index, (name, auth_number, insertion, xyz) in enumerate(residues, start=1):
            residue = gemmi.Residue()
            residue.name = name
            residue.seqid = gemmi.SeqId(auth_number, insertion or " ")
            residue.subchain = label_chain
            residue.entity_id = str(entity_number)
            residue.entity_type = gemmi.EntityType.Polymer
            residue.label_seq = label_index
            atom = gemmi.Atom()
            atom.name = "CA"
            atom.element = gemmi.Element("C")
            atom.pos = gemmi.Position(*xyz)
            atom.occ = 1.0
            residue.add_atom(atom)
            chain.add_residue(residue)
            sequence.append(name)
        model.add_chain(chain)
        entity = gemmi.Entity(str(entity_number))
        entity.entity_type = gemmi.EntityType.Polymer
        entity.polymer_type = gemmi.PolymerType.PeptideL
        entity.subchains = [label_chain]
        entity.full_sequence = sequence
        structure.entities.append(entity)
    structure.add_model(model)
    return structure


def write_cif(path: Path, structure: gemmi.Structure) -> None:
    structure.make_mmcif_document().write_file(str(path))


class StructureInventoryTests(unittest.TestCase):
    def test_gemmi_roundtrip_preserves_auth_label_and_insertion(self) -> None:
        structure = make_structure(
            [("A", "L1", [("ALA", 10, "A", (0.0, 0.0, 0.0)), ("CYS", 11, "", (1.0, 0.0, 0.0))])]
        )
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "fixture.cif"
            write_cif(path, structure)
            read = read_single_model_structure(path)
            residues = extract_confirmed_chain_residues(
                read, ChainSelector("fixture", "A", "L1")
            )
            self.assertEqual(residues[0].key.auth_seq_id, 10)
            self.assertEqual(residues[0].key.insertion_code, "A")
            self.assertEqual(residues[0].key.label_seq_id, 1)
            inventory = chain_inventory_rows(
                read, source_model_name="fixture", source_file=path,
                authoritative_sequence="AC",
            )
            self.assertEqual(inventory[0]["auth_asym_id"], "A")
            self.assertEqual(inventory[0]["label_asym_id"], "L1")
            self.assertEqual(inventory[0]["confirmation_status"], "pending_user_review")

    def test_analysis_setup_is_clone_only_and_topology_check_detects_change(self) -> None:
        structure = make_structure([("A", "L1", [("ALA", 1, "", (0, 0, 0))])])
        original_subchain = structure[0][0][0].subchain
        analysis, metadata = prepare_heuristic_analysis_copy(structure)
        self.assertIsNot(analysis, structure)
        self.assertEqual(structure[0][0][0].subchain, original_subchain)
        self.assertTrue(metadata["source_export_modified"] is False)
        moved = structure.clone()
        moved[0][0][0][0].pos.x = 100
        require_matching_topology(structure, moved, native_label="a", reference_label="b")
        moved[0][0][0][0].name = "CB"
        with self.assertRaises(Exception):
            require_matching_topology(structure, moved, native_label="a", reference_label="b")

    def test_rigid_coordinate_relationship_reports_and_rejects_residual(self) -> None:
        native = make_structure([
            ("A", "L1", [
                ("ALA", 1, "", (0, 0, 0)), ("CYS", 2, "", (1, 0, 0)),
                ("ASP", 3, "", (0, 1, 0)), ("GLU", 4, "", (0, 0, 1)),
            ])
        ])
        reference = native.clone()
        for chain in reference[0]:
            for residue in chain:
                for atom in residue:
                    x, y, z = atom.pos.x, atom.pos.y, atom.pos.z
                    atom.pos = gemmi.Position(-y + 4, x - 2, z + 3)
        result = rigid_coordinate_relationship(native, reference)
        self.assertEqual(result["atom_site_count"], 4)
        self.assertLess(result["maximum_residual_angstrom"], 1e-12)
        reference[0][0][0][0].pos.x += 0.1
        with self.assertRaises(Exception):
            rigid_coordinate_relationship(native, reference)

    def test_residue_inventory_retains_water_hetero_altloc_and_occupancy(self) -> None:
        structure = make_structure([("A", "LP", [("ALA", 1, "", (0, 0, 0))])])
        chain = structure[0][0]
        water = gemmi.Residue()
        water.name = "HOH"
        water.seqid = gemmi.SeqId(2, " ")
        water.subchain = "LW"
        water.entity_type = gemmi.EntityType.Water
        oxygen = gemmi.Atom()
        oxygen.name = "O"
        oxygen.element = gemmi.Element("O")
        oxygen.altloc = "A"
        oxygen.occ = 0.5
        oxygen.pos = gemmi.Position(1, 1, 1)
        water.add_atom(oxygen)
        chain.add_residue(water)
        ligand = gemmi.Residue()
        ligand.name = "LIG"
        ligand.seqid = gemmi.SeqId(3, " ")
        ligand.subchain = "LX"
        ligand.entity_type = gemmi.EntityType.NonPolymer
        carbon = gemmi.Atom()
        carbon.name = "C1"
        carbon.element = gemmi.Element("C")
        carbon.altloc = "B"
        carbon.occ = 0.75
        carbon.pos = gemmi.Position(2, 2, 2)
        ligand.add_atom(carbon)
        chain.add_residue(ligand)
        rows = residue_inventory_rows(structure, source_model_name="fixture")
        by_name = {row["residue_name"]: row for row in rows}
        self.assertTrue(by_name["HOH"]["is_water"])
        self.assertEqual(by_name["HOH"]["altlocs"], "A")
        self.assertEqual(by_name["HOH"]["minimum_occupancy"], 0.5)
        self.assertEqual(by_name["LIG"]["entity_type"], "NonPolymer")
        self.assertEqual(by_name["LIG"]["altlocs"], "B")
        counts = atom_site_classification_counts(structure)
        self.assertEqual(counts["atom_site_count"], 3)
        self.assertEqual(counts["nonblank_altloc_site_count"], 2)
        self.assertEqual(counts["occupancy_partial_site_count"], 2)

    def test_structure_counts_group_altloc_sites_into_one_atom_object(self) -> None:
        structure = make_structure(
            [("A", "L1", [("ALA", 1, "", (0.0, 0.0, 0.0))])]
        )
        first = structure[0][0][0][0]
        first.altloc = "A"
        first.occ = 0.5
        alternate = gemmi.Atom()
        alternate.name = "CA"
        alternate.element = gemmi.Element("C")
        alternate.altloc = "B"
        alternate.occ = 0.5
        alternate.pos = gemmi.Position(0.1, 0.0, 0.0)
        structure[0][0][0].add_atom(alternate)
        counts = structure_count_summary(structure)
        self.assertEqual(
            counts,
            {
                "model_count": 1,
                "chain_count": 1,
                "residue_count": 1,
                "atom_object_count": 1,
                "atom_site_count": 2,
            },
        )


class ResidueMappingTests(unittest.TestCase):
    def test_exact_full_and_unique_missing_residue_mapping(self) -> None:
        full = unique_exact_order_mapping("ACDE", "ACDE")
        self.assertEqual(full.authoritative_index_1based_by_observed_index, (1, 2, 3, 4))
        missing = unique_exact_order_mapping("ACDE", "ACE")
        self.assertEqual(missing.authoritative_index_1based_by_observed_index, (1, 2, 4))

    def test_fixed_blosum62_branch_resolves_only_unique_optimal_mapping(self) -> None:
        mapping = exact_first_sequence_mapping("AAAC", "AC")
        self.assertEqual(mapping.method, "blosum62_global_unique_optimal_mapping")
        self.assertEqual(mapping.authoritative_index_1based_by_observed_index, (3, 4))
        self.assertEqual(mapping.optimal_global_alignment_count, 1)
        with self.assertRaises(ResidueMappingError):
            exact_first_sequence_mapping("AAAA", "AA")
        with self.assertRaises(ResidueMappingError):
            exact_first_sequence_mapping("AC", "AD")

    def test_numbering_contract_accepts_pass_and_validates_wild_type(self) -> None:
        rows = [
            {
                "sample_uid": "x", "sequence_sha256": "",
                "scheme": "imgt", "sequence_index_1based": "1",
                "residue_aa": "A", "numbering_status": "pass", "is_gap": "false",
            }
        ]
        import hashlib
        rows[0]["sequence_sha256"] = hashlib.sha256(b"A").hexdigest()
        indexed = validate_numbering_rows(
            rows, sample_uid="x", authoritative_sequence="A",
            authoritative_sequence_sha256=rows[0]["sequence_sha256"], required_scheme="imgt",
        )
        self.assertEqual(indexed[1]["numbering_status"], "pass")
        rows[0]["residue_aa"] = "C"
        with self.assertRaises(ResidueMappingError):
            validate_numbering_rows(
                rows, sample_uid="x", authoritative_sequence="A",
                authoritative_sequence_sha256=rows[0]["sequence_sha256"], required_scheme="imgt",
            )

    def test_explicit_framework_kabsch(self) -> None:
        mobile = {1: (0, 0, 0), 2: (1, 0, 0), 3: (0, 1, 0), 4: (0, 0, 2)}
        rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], dtype=float)
        translation = np.array([5, -2, 3], dtype=float)
        reference = {index: tuple(rotation @ np.asarray(point) + translation) for index, point in mobile.items()}
        fit = fit_explicit_framework_ca(reference, mobile, framework_authoritative_indices_1based=[1, 2, 3, 4])
        transformed = apply_rigid_alignment(list(mobile.values()), fit)
        np.testing.assert_allclose(transformed, list(reference.values()), atol=1e-12)
        self.assertLess(fit.rmsd_angstrom, 1e-12)
        statistics = aligned_ca_displacement_statistics(
            reference,
            mobile,
            fit,
            region_by_authoritative_index={1: "FR1", 2: "FR1", 3: "CDR1", 4: "FR2"},
        )
        self.assertEqual(statistics["per_region"]["FR1"]["count"], 2)
        self.assertEqual(statistics["aggregates"]["CDR"]["count"], 1)
        self.assertEqual(statistics["per_region"]["CDR3"]["status"], "not_evaluable")
        self.assertIsNone(statistics["per_region"]["CDR3"]["rmsd_angstrom"])

    def test_mapping_coordinate_statuses_distinguish_missing_and_terminal(self) -> None:
        structure = make_structure(
            [("V", "LV", [("ALA", 1, "", (0, 0, 0)), ("GLY", 4, "", (3, 0, 0))])]
        )
        residues = extract_confirmed_chain_residues(
            structure, ChainSelector("model", "V", "LV")
        )
        numbering = {
            index: {
                "numbering_status": "pass",
                "numbering_position": str(index),
                "numbering_position_label": str(index),
                "region": "FR1",
            }
            for index in (1, 2, 3)
        }
        import hashlib
        sequence = "ACDG"
        rows, _ = build_mapping_rows(
            sample_uid="x",
            authoritative_sequence=sequence,
            authoritative_sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
            selector=ChainSelector("model", "V", "LV"),
            structure_residues=residues,
            numbering_by_index=numbering,
            numbering_scheme="imgt",
        )
        self.assertEqual(
            [row["coordinate_status"] for row in rows],
            ["observed", "missing_coordinates", "missing_coordinates", "terminal_flank"],
        )
        self.assertEqual(
            [row["coordinate_evaluable"] for row in rows],
            [True, False, False, True],
        )


class InterfaceTests(unittest.TestCase):
    @staticmethod
    def site(
        x: float,
        *,
        element: str = "C",
        occupancy: float = 1.0,
        altloc: str = "",
        chain: str = "V",
        is_polymer: bool = True,
    ) -> AtomSite:
        return AtomSite(
            residue=ResidueKey("model", chain, chain, 1, "", 1, "ALA"),
            atom_name="CA", element=element, altloc=altloc, occupancy=occupancy,
            x=x, y=0.0, z=0.0, is_polymer=is_polymer,
        )

    def test_strict_boundary_and_filters(self) -> None:
        vhh = [self.site(0)]
        receptor = [
            self.site(3.999999, chain="R"), self.site(4.0, chain="R"),
            self.site(4.000001, chain="R"), self.site(1.0, element="H", chain="R"),
            self.site(1.0, occupancy=0.0, chain="R"), self.site(1.0, altloc="B", chain="R"),
            self.site(1.0, element="D", chain="R"),
            self.site(1.0, chain="R", is_polymer=False),
        ]
        vhh[0] = self.site(0, altloc="A")
        contacts = strict_heavy_atom_contacts(vhh, receptor, cutoff_angstrom=4.0)
        self.assertEqual(len(contacts), 1)
        self.assertAlmostEqual(contacts[0].distance_angstrom, 3.999999)

    def test_neighbor_search_is_real_nonperiodic_candidate_path(self) -> None:
        structure = make_structure(
            [
                ("V", "LV", [("ALA", 1, "", (0, 0, 0))]),
                ("R", "LR", [("GLY", 2, "", (3.9, 0, 0)), ("SER", 3, "", (4.0, 0, 0))]),
            ]
        )
        contacts, metadata = neighbor_search_heavy_atom_contacts(
            structure,
            vhh_selector=ChainSelector("model", "V", "LV"),
            receptor_selectors=[ChainSelector("model", "R", "LR")],
            cutoff_angstrom=4.0,
        )
        self.assertEqual(len(contacts), 1)
        self.assertEqual(metadata["candidate_search"], "gemmi.NeighborSearch")
        self.assertFalse(metadata["periodic_or_symmetry_images"])

    def test_neighbor_search_excludes_crystal_mates(self) -> None:
        structure = make_structure(
            [
                ("V", "LV", [("ALA", 1, "", (0.1, 0, 0))]),
                ("R", "LR", [("GLY", 2, "", (9.9, 0, 0))]),
            ]
        )
        structure.cell = gemmi.UnitCell(10, 10, 10, 90, 90, 90)
        structure.spacegroup_hm = "P 1"
        contacts, metadata = neighbor_search_heavy_atom_contacts(
            structure,
            vhh_selector=ChainSelector("model", "V", "LV"),
            receptor_selectors=[ChainSelector("model", "R", "LR")],
            cutoff_angstrom=4.0,
        )
        self.assertEqual(contacts, [])
        self.assertFalse(metadata["periodic_or_symmetry_images"])

    def test_neighbor_search_excludes_water_and_ligand_residues(self) -> None:
        structure = make_structure(
            [
                ("V", "LV", [("ALA", 1, "", (0.0, 0, 0))]),
                ("R", "LR", [("GLY", 2, "", (5.0, 0, 0))]),
            ]
        )
        receptor_chain = structure[0][1]
        for residue_name, entity_type, x in (
            ("HOH", gemmi.EntityType.Water, 1.0),
            ("LIG", gemmi.EntityType.NonPolymer, 2.0),
        ):
            residue = gemmi.Residue()
            residue.name = residue_name
            residue.seqid = gemmi.SeqId(10 + len(receptor_chain), " ")
            residue.subchain = "LR"
            residue.entity_type = entity_type
            atom = gemmi.Atom()
            atom.name = "O" if residue_name == "HOH" else "C1"
            atom.element = gemmi.Element("O" if residue_name == "HOH" else "C")
            atom.pos = gemmi.Position(x, 0, 0)
            atom.occ = 1.0
            residue.add_atom(atom)
            receptor_chain.add_residue(residue)
        contacts, _ = neighbor_search_heavy_atom_contacts(
            structure,
            vhh_selector=ChainSelector("model", "V", "LV"),
            receptor_selectors=[ChainSelector("model", "R", "LR")],
            cutoff_angstrom=4.0,
        )
        self.assertEqual(contacts, [])


class EntryPointTests(unittest.TestCase):
    @staticmethod
    def load_script_module(name: str, filename: str):
        path = PROJECT_ROOT / "scripts" / "input_baseline" / filename
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        spec.loader.exec_module(module)
        return module

    def test_exporter_imports_without_chimerax_and_has_best_guess_false(self) -> None:
        path = PROJECT_ROOT / "scripts" / "input_baseline" / "export_cxs_session_chimerax.py"
        module = self.load_script_module("test_cxs_exporter", path.name)
        structure_module = self.load_script_module(
            "test_structure_model_identity", "build_structure_baseline.py"
        )
        self.assertEqual(
            module.EXPECTED_MODEL_NAMES,
            (
                "NK2R-252.pdb",
                "NK2R-NKA",
                "fold_2r_252_nomg_model_0.cif",
            ),
        )
        self.assertEqual(module.EXPECTED_MODEL_NAMES, structure_module.EXPECTED_MODELS)
        self.assertEqual(module._rgba((255, 165, 0, 255)), (255, 165, 0, 255))
        self.assertIn("bestGuess false", path.read_text(encoding="utf-8"))

    def test_exporter_runscript_sandbox_invokes_main(self) -> None:
        path = PROJECT_ROOT / "scripts" / "input_baseline" / "export_cxs_session_chimerax.py"
        session_sentinel = object()
        invocation: dict[str, object] = {}

        class MainInvoked(Exception):
            pass

        def trace_main_call(frame, event, arg):
            if (
                event == "call"
                and frame.f_code.co_name == "main"
                and Path(frame.f_code.co_filename) == path
            ):
                invocation["session"] = frame.f_locals["chimerax_session"]
                invocation["argv"] = frame.f_locals["argv"]
                raise MainInvoked
            return trace_main_call

        original_argv = sys.argv
        original_trace = sys.gettrace()
        try:
            sys.argv = [str(path)]
            sys.settrace(trace_main_call)
            with self.assertRaises(MainInvoked):
                runpy.run_path(
                    str(path),
                    init_globals={"session": session_sentinel},
                    run_name="chimerax_runscript_sandbox",
                )
        finally:
            sys.settrace(original_trace)
            sys.argv = original_argv

        self.assertIs(invocation["session"], session_sentinel)
        self.assertEqual(invocation["argv"], [])

    def test_exporter_normal_import_does_not_invoke_main(self) -> None:
        path = PROJECT_ROOT / "scripts" / "input_baseline" / "export_cxs_session_chimerax.py"
        invocations: list[object] = []

        def trace_main_call(frame, event, arg):
            if (
                event == "call"
                and frame.f_code.co_name == "main"
                and Path(frame.f_code.co_filename) == path
            ):
                invocations.append(frame.f_locals.get("chimerax_session"))
            return trace_main_call

        original_trace = sys.gettrace()
        try:
            sys.settrace(trace_main_call)
            self.load_script_module("test_cxs_exporter_normal_import", path.name)
        finally:
            sys.settrace(original_trace)

        self.assertEqual(invocations, [])

    def test_interface_contract_has_atom_pairs_and_complete_128_row_comparison(self) -> None:
        module = self.load_script_module(
            "test_interface_entry", "calculate_temporary_interface.py"
        )
        self.assertEqual(
            module.OUTPUT_NAMES["atom_contacts"],
            "temporary_interface_atom_contacts.csv",
        )
        self.assertEqual(
            module.OUTPUT_NAMES["orange_vs_4A_comparison"], "orange_vs_4A.csv"
        )
        self.assertIn("temporary_protected_union", module.COMPARISON_FIELDS)
        self.assertEqual(
            module._comparison_class(evaluable=True, orange=False, distance=True),
            "distance_only",
        )
        self.assertEqual(
            module._comparison_class(evaluable=False, orange=True, distance=True),
            "not_evaluable",
        )
        self.assertEqual(
            module._comparison_evaluability(
                coordinate_status="terminal_flank", coordinate_evaluable=True
            ),
            (True, ""),
        )
        self.assertEqual(
            module._comparison_evaluability(
                coordinate_status="terminal_flank", coordinate_evaluable=False
            ),
            (False, "terminal_flank_missing_coordinates"),
        )
        complete = [
            {"sample_uid": "LTT__Nb252", "sequence_index_1based": str(index)}
            for index in range(1, 129)
        ]
        module._validate_complete_nb252_mapping(complete)
        with self.assertRaises(module.InterfaceBuildBlocked):
            module._validate_complete_nb252_mapping(complete[:-1])

        vhh = InterfaceTests.site(0.0, chain="V")
        partner = InterfaceTests.site(3.0, chain="R")
        row = ContactPair(vhh, partner, 3.0).as_row()
        self.assertEqual(set(row), set(module.ATOM_CONTACT_FIELDS))

    def test_export_atom_site_counter_enumerates_altloc_occupancies(self) -> None:
        module = self.load_script_module(
            "test_cxs_count_entry", "export_cxs_session_chimerax.py"
        )

        class FakeAtom:
            def __init__(self, occupancy, alt_locs=(), alt_loc="", alt_occ=None):
                self.occupancy = occupancy
                self.alt_locs = alt_locs
                self.alt_loc = alt_loc
                self.alt_occ = alt_occ or {}

            def get_alt_loc_occupancy(self, altloc):
                return self.alt_occ[altloc]

        counts = module.atom_site_classification_counts([
            FakeAtom(1.0),
            FakeAtom(0.4, ("A", "B"), "A", {"A": 0.4, "B": 0.6}),
        ])
        self.assertEqual(counts["atom_site_count"], 3)
        self.assertEqual(counts["blank_altloc_site_count"], 1)
        self.assertEqual(counts["nonblank_altloc_site_count"], 2)
        self.assertEqual(counts["occupancy_partial_site_count"], 2)

    def test_exporter_inventories_non_atomic_session_models(self) -> None:
        module = self.load_script_module(
            "test_cxs_session_inventory", "export_cxs_session_chimerax.py"
        )

        class Parent:
            id_string = "1"
            name = "parent"

        class Surface:
            id_string = "1.1"
            name = "surface"
            atomspec = "#1.1"
            parent = Parent()
            display = True

            @staticmethod
            def child_models():
                return []

        record = module._session_model_record(
            Surface(), is_atomic_structure=False
        )
        self.assertFalse(record["is_atomic_structure"])
        self.assertEqual(record["parent_model_id"], "1")
        self.assertEqual(record["python_class"], "Surface")

    def test_single_review_template_and_confirmed_schema(self) -> None:
        module = self.load_script_module(
            "test_structure_build_entry", "build_structure_baseline.py"
        )
        inventory = [{
            "source_model_name": "NK2R-252.pdb", "auth_asym_id": "V",
            "label_asym_id": "LV", "entity_id": "1", "entity_type": "Polymer",
            "observed_sequence_sha256": "abc", "residue_count": 1, "atom_count": 1,
        }]
        color_rows = [{
            "model_name": "NK2R-252.pdb", "chimerax_chain_id": "V",
            "mmcif_chain_id": "LV", "auth_seq_id": "1", "insertion_code": "",
            "residue_name": "ALA", "ribbon_rgba": "255,165,0,255",
            "atom_rgba_histogram_json": "{}", "surface_rgba_histogram_json": "{}",
        }]
        binding = {
            "cxs_export_manifest_sha256": "a", "cxs_residue_colors_sha256": "b",
            "chain_identity_sha256": chain_identity_sha256(inventory),
        }
        review = build_review_template(
            inventory_rows=inventory, color_rows=color_rows, source_binding=binding
        )
        self.assertEqual(review["status"], "pending_user_review")
        review["status"] = "confirmed"
        chain = review["chain_reviews"][0]
        chain.update({
            "confirmed_role": "Nb252_VHH", "confirmation_status": "confirmed",
            "confirmed_by": "reviewer", "confirmed_at": "2026-08-03T00:00:00+08:00",
            "confirmation_note": "test-only confirmation",
        })
        orange = review["orange_annotation_review"]
        orange.update({
            "status": "confirmed", "confirmed_rgb": [255, 165, 0],
            "confirmed_rgba": [255, 165, 0, 255], "confirmed_channels": ["ribbon"],
            "confirmed_by": "reviewer", "confirmed_at": "2026-08-03T00:00:00+08:00",
            "confirmation_note": "test-only confirmation", "evidence": "test-only fixture",
        })
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "baseline_review.json"
            path.write_text(json.dumps(review), encoding="utf-8")
            roles, confirmed_orange = load_confirmed_review(
                path, inventory_rows=inventory, color_rows=color_rows,
                expected_binding=binding,
            )
        self.assertEqual(next(iter(roles.values()))["confirmed_role"], "Nb252_VHH")
        self.assertEqual(confirmed_orange["confirmed_rgba"], [255, 165, 0, 255])

    def test_review_template_preserves_transparent_exact_orange_classes(self) -> None:
        inventory = [{
            "source_model_name": "NK2R-252.pdb", "auth_asym_id": "V",
            "label_asym_id": "LV", "entity_id": "1", "entity_type": "Polymer",
            "observed_sequence_sha256": "abc", "residue_count": 1, "atom_count": 1,
        }]
        color_rows = [{
            "model_name": "NK2R-252.pdb", "chimerax_chain_id": "V",
            "mmcif_chain_id": "LV", "auth_seq_id": "1", "insertion_code": "",
            "residue_name": "ALA", "ribbon_rgba": "255,165,0,96",
            "atom_rgba_histogram_json": json.dumps({
                "255,165,0,64": 2,
                "255,164,0,255": 1,
            }),
            "surface_rgba_histogram_json": json.dumps({"255,165,0,96": 3}),
        }]
        review = build_review_template(
            inventory_rows=inventory,
            color_rows=color_rows,
            source_binding={"test_only": "binding"},
        )
        orange = review["orange_annotation_review"]
        self.assertEqual(orange["candidate_exact_builtin_orange_residue_count"], 1)
        classes = {
            tuple(candidate["exact_rgba"]): candidate["candidate_channels"]
            for candidate in orange["candidates"]
        }
        self.assertEqual(classes[(255, 165, 0, 64)], ["atom"])
        self.assertEqual(
            classes[(255, 165, 0, 96)],
            ["ribbon", "surface_vertex_or_uniform_patch"],
        )
        self.assertNotIn((255, 164, 0, 255), classes)

    def test_authoritative_construct_requires_hash_boundaries_gs_and_evidence(self) -> None:
        reported = "A" * 126 + "GS"
        reported_hash = hashlib.sha256(reported.encode("ascii")).hexdigest()
        review = build_review_template(
            inventory_rows=[],
            color_rows=[],
            source_binding={"test_only": "binding"},
            reported_sequence=reported,
            reported_sequence_sha256=reported_hash,
        )
        construct = review["authoritative_construct_review"]
        authoritative = reported[:126]
        construct.update(
            {
                "status": "confirmed",
                "authoritative_sequence": authoritative,
                "authoritative_sequence_sha256": hashlib.sha256(
                    authoritative.encode("ascii")
                ).hexdigest(),
                "reported_start_1based_inclusive": 1,
                "reported_end_1based_inclusive": 126,
                "construct_scope": "mature_vhh",
                "terminal_gs_interpretation": "excluded_expression_flank",
                "confirmed_by": "test reviewer",
                "confirmed_at": "2026-08-03T00:00:00+08:00",
                "confirmation_note": "test-only exact boundary confirmation",
                "evidence": "test-only collaborator record",
            }
        )
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "baseline_review.json"
            path.write_text(json.dumps(review), encoding="utf-8")
            confirmed = load_authoritative_construct_confirmation(
                path,
                reported_sequence=reported,
                reported_sequence_sha256=reported_hash,
            )
            self.assertEqual(confirmed["reported_end_1based_inclusive"], 126)
            construct["authoritative_sequence_sha256"] = "0" * 64
            path.write_text(json.dumps(review), encoding="utf-8")
            with self.assertRaises(BaselineReviewError):
                load_authoritative_construct_confirmation(
                    path,
                    reported_sequence=reported,
                    reported_sequence_sha256=reported_hash,
                )

        build_script = (
            PROJECT_ROOT
            / "scripts"
            / "input_baseline"
            / "build_structure_baseline.py"
        ).read_text(encoding="utf-8")
        self.assertNotIn("--authoritative-sequence-confirmed-by", build_script)

    def test_export_structure_count_contract_includes_all_cross_tool_levels(self) -> None:
        module = self.load_script_module(
            "test_structure_count_contract", "build_structure_baseline.py"
        )
        models = []
        for name in module.EXPECTED_MODELS:
            models.append(
                {
                    "model_name": name,
                    "coordinate_set_count": 1,
                    "chain_count": 2,
                    "residue_count": 3,
                    "atom_count": 4,
                    "atom_site_classification_counts": {
                        "atom_site_count": 5
                    },
                }
            )
        counts = module._export_model_structure_counts({"models": models})
        self.assertEqual(
            counts[module.EXPECTED_MODELS[0]],
            {
                "model_count": 1,
                "chain_count": 2,
                "residue_count": 3,
                "atom_object_count": 4,
                "atom_site_count": 5,
            },
        )
        del models[0]["chain_count"]
        with self.assertRaises(module.StructureBuildBlocked):
            module._export_model_structure_counts({"models": models})

    def test_structure_csv_writers_use_bom_and_lf(self) -> None:
        exporter = self.load_script_module(
            "test_cxs_csv_writer", "export_cxs_session_chimerax.py"
        )
        structure = self.load_script_module(
            "test_structure_csv_writer", "build_structure_baseline.py"
        )
        interface = self.load_script_module(
            "test_interface_csv_writer", "calculate_temporary_interface.py"
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            targets = [root / "export.csv", root / "structure.csv", root / "interface.csv"]
            exporter._write_csv(targets[0], ["x"], [{"x": "1"}])
            structure._write_csv(targets[1], [{"x": "1"}])
            interface._write_csv(targets[2], ["x"], [{"x": "1"}])
            for target in targets:
                content = target.read_bytes()
                self.assertTrue(content.startswith(b"\xef\xbb\xbf"))
                self.assertNotIn(b"\r\n", content)
                self.assertIn(b"x\n1\n", content)

    def test_structure_json_writers_use_lf(self) -> None:
        exporter = self.load_script_module(
            "test_cxs_json_writer", "export_cxs_session_chimerax.py"
        )
        structure = self.load_script_module(
            "test_structure_json_writer", "build_structure_baseline.py"
        )
        interface = self.load_script_module(
            "test_interface_json_writer", "calculate_temporary_interface.py"
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            targets = [root / "export.json", root / "structure.json", root / "interface.json"]
            for module, target in zip(
                (exporter, structure, interface), targets, strict=True
            ):
                module._write_json(target, {"x": 1})
                content = target.read_bytes()
                self.assertNotIn(b"\r", content)
                self.assertTrue(content.endswith(b"\n"))
                self.assertEqual(json.loads(content), {"x": 1})

    def test_missing_export_writes_only_blocked_structure_manifest(self) -> None:
        with tempfile.TemporaryDirectory(dir=PROJECT_ROOT) as temporary:
            root = Path(temporary)
            output = root / "structure"
            summary = root / "summary.json"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(PROJECT_ROOT / "scripts" / "input_baseline" / "build_structure_baseline.py"),
                    "--cxs-export-dir", str(root / "missing_export"),
                    "--expression-records", str(PROJECT_ROOT / "docs" / "result_artifacts" / "nb_expression" / "nb_expression_records.csv"),
                    "--numbering-positions", str(PROJECT_ROOT / "docs" / "result_artifacts" / "input_baseline" / "sequence" / "sequence_numbering_positions.csv"),
                    "--output-dir", str(output), "--run-summary", str(summary),
                ],
                capture_output=True, text=True, check=False,
            )
            self.assertEqual(completed.returncode, 2, completed.stderr)
            self.assertEqual({path.name for path in output.iterdir()}, {"structure_baseline_manifest.json"})
            manifest = json.loads((output / "structure_baseline_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["status"], "blocked")
            self.assertFalse((output / "model_chain_inventory.csv").exists())

    def test_missing_baseline_writes_only_blocked_interface_manifest(self) -> None:
        with tempfile.TemporaryDirectory(dir=PROJECT_ROOT) as temporary:
            root = Path(temporary)
            output = root / "interface"
            summary = root / "summary.json"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(PROJECT_ROOT / "scripts" / "input_baseline" / "calculate_temporary_interface.py"),
                    "--structure-baseline-manifest", str(root / "missing.json"),
                    "--confirmed-review", str(root / "missing_review.json"),
                    "--output-dir", str(output), "--run-summary", str(summary),
                ],
                capture_output=True, text=True, check=False,
            )
            self.assertEqual(completed.returncode, 2, completed.stderr)
            self.assertEqual({path.name for path in output.iterdir()}, {"interface_manifest.json"})
            manifest = json.loads((output / "interface_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["interface_status"], "blocked")


if __name__ == "__main__":
    unittest.main()
