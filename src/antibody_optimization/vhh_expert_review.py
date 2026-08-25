"""Structure-aware expert review helpers for the active V3 Nb252 shortlist.

The module derives descriptive residue environments from the released
experimental NK2R--Nb252 complex and the independently predicted AF3 VHH.
It then joins those facts to curated, mutation-specific expert judgements.
The judgements are hypotheses for triage: they do not select parents, model
mutant coordinates, predict affinity, or replace experimental measurements.
"""

from __future__ import annotations

import copy
import math
from collections import defaultdict
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
from Bio.PDB import MMCIFParser, PPBuilder, ShrakeRupley


MAX_ASA = {
    "A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0,
    "Q": 225.0, "E": 223.0, "G": 104.0, "H": 224.0, "I": 197.0,
    "L": 201.0, "K": 236.0, "M": 224.0, "F": 240.0, "P": 159.0,
    "S": 155.0, "T": 172.0, "W": 285.0, "Y": 263.0, "V": 174.0,
}
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
HYDROPHOBIC = frozenset("AVILMFWY")
AROMATIC = frozenset("FWY")
POSITIVE = frozenset("KR")
NEGATIVE = frozenset("DE")
POLAR = frozenset("STNQH")
SIDECHAIN_VOLUME = {
    "A": 67.0, "R": 148.0, "N": 96.0, "D": 91.0, "C": 86.0,
    "Q": 114.0, "E": 109.0, "G": 48.0, "H": 118.0, "I": 124.0,
    "L": 124.0, "K": 135.0, "M": 124.0, "F": 135.0, "P": 90.0,
    "S": 73.0, "T": 93.0, "W": 163.0, "Y": 141.0, "V": 105.0,
}
KYTE_DOOLITTLE = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O", "OXT"})


class VHHExpertReviewError(ValueError):
    """Raised when shortlist, mapping, or structure identities disagree."""


def derive_position_contexts(
    mapping_rows: Sequence[Mapping[str, object]],
    positions: Sequence[int],
    experimental_cif: Path,
    af3_cif: Path,
    *,
    alignment_summary: Mapping[str, object] | None = None,
    interface_positions: Sequence[int] = (),
) -> dict[int, dict[str, object]]:
    """Return source-aware structural context for requested reported positions.

    Experimental coordinates are primary when observed.  AF3 is always
    calculated as a sensitivity view and becomes primary only when the
    experimental residue has no coordinates.  Solvent accessibility uses
    Shrake--Rupley with 200 points; relative accessibility is residue SASA
    divided by the Tien et al. 2013 theoretical maximum-accessibility scale.
    Project-descriptive exposure bins are RSA <=0.10 buried, 0.10--<0.25
    partially buried, and >=0.25 exposed.  They are not predictor thresholds.
    The phi/psi class is a descriptive heuristic, not a DSSP assignment.
    AF3 ``B_iso_or_equiv`` is retained without assuming that it is pLDDT.
    """

    requested = sorted(set(int(value) for value in positions))
    if not requested:
        raise VHHExpertReviewError("At least one reported position is required")
    mapping = _mapping_by_position(mapping_rows)
    missing = [position for position in requested if position not in mapping]
    if missing:
        raise VHHExpertReviewError(f"Structure mapping lacks positions: {missing}")

    experimental = _load_model(experimental_cif)
    af3 = _load_model(af3_cif)
    if "C" not in experimental or "R" not in experimental:
        raise VHHExpertReviewError("Experimental model must contain VHH C and NK2R R")
    if "A" not in af3:
        raise VHHExpertReviewError("AF3 model must contain VHH chain A")

    exp_context = _chain_context(
        experimental,
        "C",
        requested,
        mapping,
        source_name="experimental_complex",
        receptor_chain_id="R",
        interface_positions=interface_positions,
    )
    af3_context = _chain_context(
        af3,
        "A",
        requested,
        mapping,
        source_name="af3_vhh_only",
        receptor_chain_id=None,
        interface_positions=(),
    )
    transform = _alignment_transform(alignment_summary)
    experimental_missing = sorted(
        position
        for position, model_rows in mapping.items()
        if (
            (row := model_rows.get("NK2R-252.pdb")) is not None
            and str(row["coordinate_status"]) != "observed"
        )
    )
    output: dict[int, dict[str, object]] = {}
    for position in requested:
        exp = exp_context.get(position)
        predicted = af3_context.get(position)
        if predicted is None:
            raise VHHExpertReviewError(f"AF3 lacks reported position {position}")
        observed = exp is not None
        primary = exp if observed else predicted
        assert primary is not None
        gap_distance = min(
            (abs(position - missing_position) for missing_position in experimental_missing),
            default=None,
        )
        if not observed:
            gap_proximity = "within_experimental_missing_coordinate_segment"
        elif gap_distance == 1:
            gap_proximity = "adjacent_to_experimental_missing_coordinate_segment"
        elif gap_distance is not None and gap_distance <= 2:
            gap_proximity = "within_two_positions_of_experimental_missing_coordinate_segment"
        else:
            gap_proximity = "not_near_experimental_missing_coordinate_segment"
        ca_displacement = ""
        if exp is not None and transform is not None:
            rotation, translation = transform
            predicted_ca = rotation @ np.asarray(predicted["ca_xyz"], dtype=float) + translation
            ca_displacement = round(
                float(np.linalg.norm(predicted_ca - np.asarray(exp["ca_xyz"], dtype=float))),
                5,
            )
        interface_distance = "" if exp is None else exp["minimum_interface_residue_distance_a"]
        if exp is None:
            near_interface = "not_evaluable_experimental_coordinates_missing"
        elif interface_distance != "" and float(interface_distance) <= 4.0:
            near_interface = "within_4A_of_hard_interface_residue"
        elif interface_distance != "" and float(interface_distance) <= 4.5:
            near_interface = "borderline_4_to_4p5A_from_hard_interface_residue"
        else:
            near_interface = "outside_4p5A_hard_interface_shell"
        output[position] = {
            "reported_sequence_index_1based": position,
            "experimental_coordinate_status": "observed" if observed else "missing_coordinates",
            "primary_structure_source": "experimental_complex" if observed else "af3_only_due_missing_experimental_coordinates",
            "structure_evidence_limit": (
                "experimental_coordinate_observation_with_af3_sensitivity_view"
                if observed
                else "predicted_context_only_not_experimentally_observed"
            ),
            "primary_chain_id": "C" if observed else "A",
            "primary_auth_seq_id": primary["auth_seq_id"],
            "primary_relative_sasa": primary["relative_sasa"],
            "primary_exposure_class": primary["exposure_class"],
            "primary_phi_degrees": primary["phi_degrees"],
            "primary_psi_degrees": primary["psi_degrees"],
            "primary_backbone_class_heuristic": primary["backbone_class_heuristic"],
            "primary_intra_vhh_neighbor_count_4p5a": primary["intra_vhh_neighbor_count_4p5a"],
            "primary_intra_vhh_neighbors_4p5a": primary["intra_vhh_neighbors_4p5a"],
            "primary_nearest_disulfide_sulfur_distance_a": primary["nearest_disulfide_sulfur_distance_a"],
            "primary_nearest_disulfide_cys": primary["nearest_disulfide_cys"],
            "experimental_relative_sasa": "" if exp is None else exp["relative_sasa"],
            "experimental_exposure_class": "not_evaluable" if exp is None else exp["exposure_class"],
            "experimental_minimum_receptor_distance_a": "" if exp is None else exp["minimum_receptor_distance_a"],
            "experimental_nearest_receptor_residue": "" if exp is None else exp["nearest_receptor_residue"],
            "experimental_sasa_buried_by_receptor_a2": "" if exp is None else exp["sasa_buried_by_receptor_a2"],
            "experimental_minimum_interface_residue_distance_a": interface_distance,
            "experimental_nearest_interface_residue": "" if exp is None else exp["nearest_interface_residue"],
            "near_interface_shell_status": near_interface,
            "experimental_missing_coordinate_proximity": gap_proximity,
            "experimental_missing_coordinate_distance_positions": "" if gap_distance is None else gap_distance,
            "af3_relative_sasa": predicted["relative_sasa"],
            "af3_exposure_class": predicted["exposure_class"],
            "af3_backbone_class_heuristic": predicted["backbone_class_heuristic"],
            "af3_mean_b_iso_or_equiv_uninterpreted": predicted["mean_atom_b_factor"],
            "experimental_af3_exposure_agreement": (
                "not_applicable" if exp is None else (
                    "same_class" if exp["exposure_class"] == predicted["exposure_class"] else "different_class"
                )
            ),
            "experimental_af3_backbone_class_agreement": (
                "not_applicable" if exp is None else (
                    "same_class"
                    if exp["backbone_class_heuristic"] == predicted["backbone_class_heuristic"]
                    else "different_class"
                )
            ),
            "experimental_af3_ca_displacement_after_fr_fit_a": ca_displacement,
            "position_overview_view_id": f"reported_pos_{position:03d}",
        }
    return output


def mutation_chemistry(wt: str, mutant: str) -> dict[str, object]:
    """Return transparent amino-acid-class changes for one substitution."""

    if wt not in MAX_ASA or mutant not in MAX_ASA or wt == mutant:
        raise VHHExpertReviewError(f"Invalid substitution {wt}->{mutant}")
    return {
        "sidechain_volume_delta_a3": round(SIDECHAIN_VOLUME[mutant] - SIDECHAIN_VOLUME[wt], 3),
        "kyte_doolittle_delta": round(KYTE_DOOLITTLE[mutant] - KYTE_DOOLITTLE[wt], 3),
        "formal_charge_class_wt": _charge_class(wt),
        "formal_charge_class_mutant": _charge_class(mutant),
        "charge_class_change": f"{_charge_class(wt)}_to_{_charge_class(mutant)}",
        "hydrophobic_class_change": f"{_bool_class(wt in HYDROPHOBIC)}_to_{_bool_class(mutant in HYDROPHOBIC)}",
        "aromatic_class_change": f"{_bool_class(wt in AROMATIC)}_to_{_bool_class(mutant in AROMATIC)}",
        "glycine_introduced": mutant == "G" and wt != "G",
        "proline_introduced": mutant == "P" and wt != "P",
        "methionine_introduced": mutant == "M" and wt != "M",
        "asparagine_introduced": mutant == "N" and wt != "N",
    }


def build_expert_review_rows(
    candidate_rows: Sequence[Mapping[str, object]],
    contexts: Mapping[int, Mapping[str, object]],
    assessments: Mapping[str, Mapping[str, object]],
    *,
    visual_view_ids: Mapping[str, str] | None = None,
) -> list[dict[str, object]]:
    """Join unchanged V3 evidence, structural facts, and curated judgements."""

    if len(candidate_rows) != 30 or len({str(row["candidate_id"]) for row in candidate_rows}) != 30:
        raise VHHExpertReviewError("Expected exactly 30 unique V3 shortlist candidates")
    identifiers = {str(row["candidate_id"]) for row in candidate_rows}
    if identifiers != set(assessments):
        missing = sorted(identifiers - set(assessments))
        extra = sorted(set(assessments) - identifiers)
        raise VHHExpertReviewError(f"Assessment identity mismatch: missing={missing}, extra={extra}")

    position_counts: dict[int, int] = defaultdict(int)
    for row in candidate_rows:
        position_counts[int(row["reported_sequence_index_1based"])] += 1
    output: list[dict[str, object]] = []
    for source in sorted(candidate_rows, key=lambda row: int(row["selection_order_v3"])):
        identifier = str(source["candidate_id"])
        position = int(source["reported_sequence_index_1based"])
        wt = str(source["wt_residue"])
        mutant = str(source["mutant_residue"])
        sequence = str(source["sequence"])
        if len(sequence) != 128 or sequence[position - 1] != mutant:
            raise VHHExpertReviewError(f"Candidate sequence mismatch: {identifier}")
        context = contexts[position]
        assessment = assessments[identifier]
        required = {
            "structural_facts_cn",
            "chimerax_single_rotamer_observation_cn",
            "vhh_expert_inference_cn",
            "expert_structural_assessment",
            "expert_solubility_expectation",
            "expert_thermal_stability_expectation",
            "expert_confidence",
            "expert_primary_concern",
            "expert_rationale_cn",
            "expert_uncertainty_cn",
            "expert_rule_flags",
        }
        if not required.issubset(assessment) or any(not str(assessment[key]).strip() for key in required):
            raise VHHExpertReviewError(f"Incomplete expert assessment: {identifier}")
        expected_mutation = f"{wt}{position}{mutant}"
        if assessment.get("mutation", expected_mutation) != expected_mutation:
            raise VHHExpertReviewError(f"Assessment mutation mismatch: {identifier}")
        assessment_payload = {
            key: value for key, value in assessment.items() if key != "mutation"
        }
        row = {
            "upstream_shortlist_row_order_not_expert_rank": int(source["selection_order_v3"]),
            "candidate_id": identifier,
            "mutation_reported_label": source["mutation_reported_label"],
            "reported_sequence_index_1based": position,
            "imgt_position_label": source["imgt_position_label"],
            "region": source["region"],
            "wt_residue": wt,
            "mutant_residue": mutant,
            "same_position_candidate_count": position_counts[position],
            "conservation_class": source["conservation_class"],
            "selection_tier_v3_upstream": source["selection_tier_v3"],
            "netsolp_delta_u": float(source["netsolp_delta_usability_vs_current_wt"]),
            "netsolp_u_band_v3": source["netsolp_u_magnitude_band_v3"],
            "netsolp_delta_s": float(source["netsolp_delta_solubility_vs_current_wt"]),
            "netsolp_s_band_v3": source["netsolp_s_magnitude_band_v3"],
            "nanomelt_delta_tm_c": float(source["nanomelt_delta_predicted_apparent_tm_c_vs_current_wt"]),
            "nanomelt_tm_band_v3": source["nanomelt_tm_magnitude_band_v3"],
            "antifold_selection_source": source["antifold_selection_source"],
            "antifold_delta_logp": float(source["antifold_selection_delta_log_probability"]),
            "antifold_mutant_rank_worst_first": int(source["antifold_mutant_rank_worst_first"]),
            "antifold_veto_status": source["antifold_veto_status"],
            "upstream_hard_sequence_risk_flags": source["hard_sequence_risk_flags_v3"],
            "upstream_soft_sequence_risk_flags": source["soft_sequence_risk_flags_v3"],
            "stable_word_effect": source["stable_word_effect"],
            **context,
            **mutation_chemistry(wt, mutant),
            **assessment_payload,
            "visual_review_candidate_view_id": (
                "" if visual_view_ids is None else visual_view_ids.get(identifier, "")
            ),
            "manual_visual_review_status": (
                "not_performed" if visual_view_ids is None
                else "reviewed_in_chimerax_1_12_single_rotamer_view"
            ),
            "parent_single_selection_status": "not_performed",
            "sequence": sequence,
        }
        if (
            row["primary_structure_source"] == "af3_only_due_missing_experimental_coordinates"
            and row["expert_confidence"] != "low"
        ):
            raise VHHExpertReviewError(f"AF3-only assessment must have low confidence: {identifier}")
        if isinstance(row["expert_rule_flags"], (list, tuple)):
            row["expert_rule_flags"] = ";".join(str(value) for value in row["expert_rule_flags"])
        if visual_view_ids is not None and not row["visual_review_candidate_view_id"]:
            raise VHHExpertReviewError(f"Missing visual-review view identity: {identifier}")
        output.append(row)
    return output


def summarize_review(rows: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    """Return compact counts without interpreting them as candidate ranks."""

    categories = (
        "expert_structural_assessment",
        "expert_solubility_expectation",
        "expert_thermal_stability_expectation",
        "expert_confidence",
        "primary_structure_source",
    )
    summary: list[dict[str, object]] = []
    for field in categories:
        counts: dict[str, int] = defaultdict(int)
        for row in rows:
            counts[str(row[field])] += 1
        summary.extend(
            {"field": field, "category": category, "count": count}
            for category, count in sorted(counts.items())
        )
    return summary


def _load_model(path: Path):
    if not path.is_file():
        raise VHHExpertReviewError(f"Structure file is missing: {path}")
    structure = MMCIFParser(QUIET=True, auth_chains=True, auth_residues=True).get_structure(path.stem, str(path))
    models = list(structure.get_models())
    if len(models) != 1:
        raise VHHExpertReviewError(f"Expected one model in {path}, found {len(models)}")
    return models[0]


def _mapping_by_position(rows: Sequence[Mapping[str, object]]) -> dict[int, dict[str, Mapping[str, object]]]:
    output: dict[int, dict[str, Mapping[str, object]]] = defaultdict(dict)
    for row in rows:
        position = int(row["sequence_index_1based"])
        model = str(row["source_model_name"])
        output[position][model] = row
    return dict(output)


def _chain_context(
    model,
    chain_id: str,
    positions: Sequence[int],
    mapping: Mapping[int, Mapping[str, Mapping[str, object]]],
    *,
    source_name: str,
    receptor_chain_id: str | None,
    interface_positions: Sequence[int],
) -> dict[int, dict[str, object]]:
    chain = model[chain_id]
    model_name = "NK2R-252.pdb" if source_name == "experimental_complex" else "fold_2r_252_nomg_model_0.cif"
    residues_by_auth = {
        (int(residue.id[1]), str(residue.id[2]).strip()): residue
        for residue in chain.get_residues()
        if residue.id[0] == " " and residue.resname in THREE_TO_ONE
    }
    phi_psi = _phi_psi(chain)

    complex_model = copy.deepcopy(model)
    surface = ShrakeRupley(n_points=200)
    surface.compute(complex_model, level="R")
    complex_chain = complex_model[chain_id]
    complex_sasa_by_auth = {
        (int(residue.id[1]), str(residue.id[2]).strip()): float(residue.sasa)
        for residue in complex_chain.get_residues()
        if residue.id[0] == " " and residue.resname in THREE_TO_ONE
    }
    isolated_chain = copy.deepcopy(chain)
    surface.compute(isolated_chain, level="R")
    isolated_by_auth = {
        (int(residue.id[1]), str(residue.id[2]).strip()): residue
        for residue in isolated_chain.get_residues()
        if residue.id[0] == " " and residue.resname in THREE_TO_ONE
    }

    output: dict[int, dict[str, object]] = {}
    all_residues = [residue for residue in chain.get_residues() if residue.id[0] == " " and residue.resname in THREE_TO_ONE]
    receptor = model[receptor_chain_id] if receptor_chain_id else None
    interface_residues: list[tuple[int, object]] = []
    for interface_position in sorted(set(int(value) for value in interface_positions)):
        interface_row = mapping.get(interface_position, {}).get(model_name)
        if interface_row is None or str(interface_row["coordinate_status"]) != "observed":
            continue
        interface_key = (
            int(interface_row["auth_seq_id"]),
            str(interface_row.get("insertion_code", "")).strip(),
        )
        interface_residue = residues_by_auth.get(interface_key)
        if interface_residue is not None:
            interface_residues.append((interface_position, interface_residue))
    for position in positions:
        row = mapping[position].get(model_name)
        if row is None:
            raise VHHExpertReviewError(f"Missing {model_name} mapping at reported position {position}")
        if str(row["coordinate_status"]) != "observed":
            continue
        auth = int(row["auth_seq_id"])
        insertion = str(row.get("insertion_code", "")).strip()
        residue = residues_by_auth.get((auth, insertion))
        if residue is None:
            raise VHHExpertReviewError(f"Mapped residue absent in {model_name} {chain_id}:{auth}{insertion}")
        wt = str(row["residue_aa"])
        if THREE_TO_ONE[residue.resname] != wt:
            raise VHHExpertReviewError(f"WT mismatch at reported position {position} in {model_name}")
        isolated_residue = isolated_by_auth[(auth, insertion)]
        sasa = float(isolated_residue.sasa)
        relative = sasa / MAX_ASA[wt]
        target_atoms = _heavy_atoms(residue, sidechain=True) or _heavy_atoms(residue)
        neighbors = []
        for other in all_residues:
            if other is residue:
                continue
            distance = _minimum_distance(target_atoms, _heavy_atoms(other))
            if distance < 4.5:
                neighbors.append((distance, other))
        neighbors.sort(key=lambda item: (item[0], item[1].id[1]))
        cysteine_distances = []
        for cys_position in (22, 95):
            cys_row = mapping[cys_position].get(model_name)
            if cys_row is None or str(cys_row["coordinate_status"]) != "observed":
                continue
            cys_key = (int(cys_row["auth_seq_id"]), str(cys_row.get("insertion_code", "")).strip())
            cys = residues_by_auth.get(cys_key)
            if cys is not None and "SG" in cys:
                cysteine_distances.append((_minimum_distance(target_atoms, [cys["SG"]]), cys_position))
        closest_cys = min(cysteine_distances, default=(math.inf, None))
        phi, psi = phi_psi.get((auth, insertion), (None, None))
        receptor_distance = ""
        receptor_name = ""
        if receptor is not None:
            receptor_distance_value, receptor_residue = _nearest_residue(target_atoms, receptor.get_residues())
            receptor_distance = round(receptor_distance_value, 4)
            receptor_name = f"{THREE_TO_ONE.get(receptor_residue.resname, receptor_residue.resname)}{receptor_residue.id[1]}"
        interface_distance = ""
        interface_name = ""
        eligible_interface = [
            (interface_position, interface_residue)
            for interface_position, interface_residue in interface_residues
            if interface_position != position
        ]
        if eligible_interface:
            interface_distance_value, interface_position, interface_residue = min(
                (
                    _minimum_distance(target_atoms, _heavy_atoms(interface_residue)),
                    interface_position,
                    interface_residue,
                )
                for interface_position, interface_residue in eligible_interface
            )
            interface_distance = round(interface_distance_value, 4)
            interface_name = f"{THREE_TO_ONE[interface_residue.resname]}{interface_position}"
        atoms = list(residue.get_atoms())
        mean_b = sum(float(atom.bfactor) for atom in atoms) / len(atoms)
        output[position] = {
            "auth_seq_id": f"{auth}{insertion}",
            "relative_sasa": round(relative, 5),
            "exposure_class": _exposure_class(relative),
            "phi_degrees": "" if phi is None else round(math.degrees(phi), 3),
            "psi_degrees": "" if psi is None else round(math.degrees(psi), 3),
            "backbone_class_heuristic": _backbone_class(phi, psi),
            "intra_vhh_neighbor_count_4p5a": len(neighbors),
            "intra_vhh_neighbors_4p5a": ";".join(
                f"{THREE_TO_ONE[other.resname]}{other.id[1]}:{distance:.2f}"
                for distance, other in neighbors[:12]
            ),
            "nearest_disulfide_sulfur_distance_a": "" if closest_cys[1] is None else round(closest_cys[0], 4),
            "nearest_disulfide_cys": "" if closest_cys[1] is None else f"C{closest_cys[1]}",
            "minimum_receptor_distance_a": receptor_distance,
            "nearest_receptor_residue": receptor_name,
            "sasa_buried_by_receptor_a2": (
                "" if receptor is None else round(sasa - complex_sasa_by_auth[(auth, insertion)], 4)
            ),
            "minimum_interface_residue_distance_a": interface_distance,
            "nearest_interface_residue": interface_name,
            "mean_atom_b_factor": round(mean_b, 4),
            "ca_xyz": tuple(float(value) for value in residue["CA"].coord),
        }
    return output


def _phi_psi(chain) -> dict[tuple[int, str], tuple[float | None, float | None]]:
    output: dict[tuple[int, str], tuple[float | None, float | None]] = {}
    for peptide in PPBuilder().build_peptides(chain, aa_only=True):
        for residue, angles in zip(peptide, peptide.get_phi_psi_list()):
            output[(int(residue.id[1]), str(residue.id[2]).strip())] = angles
    return output


def _backbone_class(phi: float | None, psi: float | None) -> str:
    if phi is None or psi is None:
        return "terminal_or_not_assignable"
    phi_d = math.degrees(phi)
    psi_d = math.degrees(psi)
    if -180 <= phi_d <= -60 and (90 <= psi_d <= 180 or -180 <= psi_d <= -120):
        return "beta_like_phi_psi"
    if -100 <= phi_d <= -30 and -80 <= psi_d <= -5:
        return "helix_like_phi_psi"
    return "loop_or_turn_like_phi_psi"


def _exposure_class(relative_sasa: float) -> str:
    if relative_sasa <= 0.10:
        return "buried"
    if relative_sasa < 0.25:
        return "partially_buried"
    return "exposed"


def _alignment_transform(
    summary: Mapping[str, object] | None,
) -> tuple[np.ndarray, np.ndarray] | None:
    if summary is None:
        return None
    reference = summary.get("reference")
    mobile = summary.get("mobile")
    if not isinstance(reference, Mapping) or reference.get("model_name") != "NK2R-252.pdb":
        raise VHHExpertReviewError("Alignment reference must be NK2R-252.pdb")
    if not isinstance(mobile, Mapping) or mobile.get("model_name") != "fold_2r_252_nomg_model_0.cif":
        raise VHHExpertReviewError("Alignment mobile model must be the AF3 VHH")
    rotation = np.asarray(summary.get("rotation_3x3"), dtype=float)
    translation = np.asarray(summary.get("translation_3"), dtype=float)
    if rotation.shape != (3, 3) or translation.shape != (3,):
        raise VHHExpertReviewError("Alignment transform must contain a 3x3 rotation and 3-vector")
    if not np.isfinite(rotation).all() or not np.isfinite(translation).all():
        raise VHHExpertReviewError("Alignment transform contains non-finite values")
    return rotation, translation


def _nearest_residue(atoms: Sequence[object], residues) -> tuple[float, object]:
    best_distance = math.inf
    best_residue = None
    for residue in residues:
        if residue.id[0] != " " or residue.resname not in THREE_TO_ONE:
            continue
        distance = _minimum_distance(atoms, _heavy_atoms(residue))
        if distance < best_distance:
            best_distance = distance
            best_residue = residue
    if best_residue is None:
        raise VHHExpertReviewError("No polymer residue available for distance calculation")
    return best_distance, best_residue


def _minimum_distance(first_atoms: Sequence[object], second_atoms: Sequence[object]) -> float:
    return min(float(first - second) for first in first_atoms for second in second_atoms)


def _heavy_atoms(residue, *, sidechain: bool = False) -> list[object]:
    return [
        atom for atom in residue.get_atoms()
        if atom.element not in {"H", "D"} and (not sidechain or atom.name not in BACKBONE_ATOMS)
    ]


def _charge_class(residue: str) -> str:
    if residue in POSITIVE:
        return "positive"
    if residue in NEGATIVE:
        return "negative"
    if residue in POLAR:
        return "polar_neutral"
    return "neutral"


def _bool_class(value: bool) -> str:
    return "yes" if value else "no"
