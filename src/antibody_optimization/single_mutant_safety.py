"""Structure-aware combination qualification for Nb252 single-mutant modules.

This module reviews the completed 50-candidate affinity ensemble and the
30-candidate property pool.  It does not predict expression or measured
affinity, build mutant coordinates, or generate combinations.  Operational
flags are conservative triage rules for this project and are not validated
universal developability thresholds.
"""

from __future__ import annotations

import math
import re
from collections import Counter
from pathlib import Path
from typing import Mapping, Sequence

from Bio.PDB import PDBParser, ShrakeRupley


AFFINITY_COUNT = 50
PROPERTY_COUNT = 30
HYDROPHOBIC = frozenset("AVILMFWY")
MAX_ASA = {
    "A": 129,
    "R": 274,
    "N": 195,
    "D": 193,
    "C": 167,
    "Q": 225,
    "E": 223,
    "G": 104,
    "H": 224,
    "I": 197,
    "L": 201,
    "K": 236,
    "M": 224,
    "F": 240,
    "P": 159,
    "S": 155,
    "T": 172,
    "W": 285,
    "Y": 263,
    "V": 174,
}
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


class SingleMutantSafetyError(ValueError):
    """Raised when the fixed review inputs or residue identities disagree."""


def review_single_mutant_modules(
    affinity_rows: Sequence[Mapping[str, object]],
    property_rows: Sequence[Mapping[str, object]],
    candidate_rows: Sequence[Mapping[str, object]],
    property_evidence_rows: Sequence[Mapping[str, object]],
    antifold_core_rows: Sequence[Mapping[str, object]],
    tnp_rows: Sequence[Mapping[str, object]],
    prepared_wt_pdb: Path,
    missing_positions: set[int],
) -> dict[str, object]:
    """Return candidate-level safety evidence and conservative qualification.

    Affinity modules retain the published 18/20 core gate.  Non-core affinity
    rows with favorable medians and at least 13/20 support for each metric are
    exposed only as ``targeted_alternative_review``; this descriptive 65%
    review floor does not promote them to combination-ready status.  Property
    modules require no directionally adverse paired PyRosetta result and no
    major structural flag.  Positions directly adjacent to unresolved
    experimental coordinates are held for targeted structural review.
    """

    _validate_unique(affinity_rows, "candidate_id", AFFINITY_COUNT, "affinity")
    _validate_unique(property_rows, "candidate_id", PROPERTY_COUNT, "property")
    if {str(row["candidate_id"]) for row in affinity_rows} & {
        str(row["candidate_id"]) for row in property_rows
    }:
        raise SingleMutantSafetyError("Affinity and property pools overlap")
    candidates = {str(row["candidate_id"]): row for row in candidate_rows}
    properties = {str(row["candidate_id"]): row for row in property_evidence_rows}
    antifold = {str(row["candidate_id"]): row for row in antifold_core_rows}
    tnp = {str(row["candidate_id"]): row for row in tnp_rows}
    requested = [*(str(row["candidate_id"]) for row in affinity_rows), *(str(row["candidate_id"]) for row in property_rows)]
    missing = [identifier for identifier in requested if identifier not in candidates]
    if missing:
        raise SingleMutantSafetyError(f"Missing unified candidate identities: {missing}")

    structure = _structure_features(prepared_wt_pdb, requested, candidates)
    output: list[dict[str, object]] = []
    for track, source_rows in (("affinity", affinity_rows), ("property", property_rows)):
        for source in source_rows:
            identifier = str(source["candidate_id"])
            candidate = candidates[identifier]
            position = int(source["sequence_index_1based"])
            wt = str(source["wt_residue"])
            mutant = str(source["mutant_residue"])
            sequence = str(candidate["sequence"])
            if sequence[position - 1] != mutant:
                raise SingleMutantSafetyError(f"Mutant sequence mismatch for {identifier}")
            wt_sequence = sequence[: position - 1] + wt + sequence[position:]
            feature = structure[position]
            evidence = properties.get(identifier, candidate)
            anti = _optional_float(
                evidence.get("experimental_complex_context_delta_log_probability")
                or antifold.get(identifier, {}).get("experimental_complex_context_delta_log_probability")
            )
            before_hydrophobic, after_hydrophobic, local_context = _local_hydrophobicity(
                wt_sequence, sequence, position
            )
            new_nglycan = _motif_count(sequence, r"N[^P][ST]") - _motif_count(wt_sequence, r"N[^P][ST]")
            flags: list[str] = []
            hard_flags: list[str] = []
            if str(candidate["design_status"]) == "blocked_new_unpaired_cys":
                hard_flags.append("new_unpaired_cysteine")
            if new_nglycan > 0:
                hard_flags.append("new_n_linked_glycosylation_motif")
            gap_distance = min(abs(position - item) for item in missing_positions)
            if gap_distance <= 1:
                flags.append("adjacent_to_experimental_coordinate_gap")
            if mutant == "P" and gap_distance <= 1:
                flags.append("proline_at_coordinate_gap_boundary")
            if mutant == "P" and wt != "P":
                flags.append("new_proline_backbone_constraint")
            if mutant == "G" and wt != "G" and str(source["region"]).startswith("CDR"):
                flags.append("cdr_glycine_flexibility_change")
            if anti is not None and anti <= -3.0:
                flags.append("strong_negative_antifold_complex_signal")
            exposed_hydrophobic = (
                mutant in HYDROPHOBIC
                and wt not in HYDROPHOBIC
                and feature["vhh_alone_relative_sasa"] >= 0.25
            )
            if exposed_hydrophobic:
                flags.append("exposed_hydrophobic_substitution")
            if after_hydrophobic >= 6 and after_hydrophobic > before_hydrophobic:
                flags.append("dense_local_hydrophobic_window")
            liability_tokens = _tokens(candidate.get("new_liability_flags"), "|")
            if "more_M_or_W" in liability_tokens:
                flags.append("increased_oxidation_susceptibility")
            source_risks = _tokens(source.get("risk_flags"), ";")
            if any("contact_change" in token or token == "vhh_contact_reorganization" for token in source_risks):
                flags.append("ensemble_contact_change")
            if "prepared_contact_sensitive" in source_risks or _bool(source.get("prepared_contact_sensitive")):
                flags.append("prepared_structure_sensitive")
            if "ensemble_fa_rep_increase" in source_risks or "fa_rep_increase" in source_risks:
                flags.append("ensemble_fa_rep_increase")
            if position == 45 and mutant in HYDROPHOBIC | {"C"}:
                flags.append("vhh_hallmark_hydrophobization")
            paired_contact_status = str(source.get("paired_contact_status", ""))
            if paired_contact_status and paired_contact_status != "preserved_all":
                flags.append("paired_receptor_contact_change")

            qualification, reason = _qualification(
                track, source, hard_flags, flags
            )
            all_flags = [*hard_flags, *flags]
            output.append(
                {
                    "candidate_id": identifier,
                    "design_track": track,
                    "mutation": f"{wt}{position}{mutant}",
                    "sequence_index_1based": position,
                    "wt_residue": wt,
                    "mutant_residue": mutant,
                    "region": source["region"],
                    "sequence": sequence,
                    "qualification_status": qualification,
                    "qualification_reason": reason,
                    "expert_risk_level": _risk_level(hard_flags, flags),
                    "hard_risk_flags": ";".join(hard_flags),
                    "structural_review_flags": ";".join(flags),
                    "risk_flag_count": len(all_flags),
                    "experimental_gap_distance_residues": gap_distance,
                    **feature,
                    "local_sequence_wt": local_context[0],
                    "local_sequence_mutant": local_context[1],
                    "local_hydrophobic_count_wt": before_hydrophobic,
                    "local_hydrophobic_count_mutant": after_hydrophobic,
                    "new_n_linked_glycosylation_motif_count": new_nglycan,
                    "cysteine_count_delta": sequence.count("C") - wt_sequence.count("C"),
                    "formal_charge_delta": int(float(candidate["formal_charge_delta"])),
                    "antifold_complex_delta_log_probability": "" if anti is None else anti,
                    "netsolp_delta_usability_vs_wt": evidence.get("netsolp_delta_usability_vs_wt", ""),
                    "netsolp_delta_solubility_vs_wt": evidence.get("netsolp_delta_solubility_vs_wt", ""),
                    "nanomelt_delta_predicted_tm_c_vs_wt": evidence.get("nanomelt_delta_predicted_apparent_tm_c_vs_wt", ""),
                    "tnp_psh_delta_vs_wt": tnp.get(identifier, {}).get("tnp_psh_delta_vs_wt", ""),
                    "tnp_flag_regression_count": tnp.get(identifier, {}).get("tnp_flag_regression_count", ""),
                    "affinity_core_support_gate": source.get("core_support_gate", "not_applicable"),
                    "affinity_negative_dg_count": source.get("negative_delta_dG_count", ""),
                    "affinity_negative_cross_count": source.get("negative_delta_cross_interface_count", ""),
                    "property_affinity_direction_class": source.get("affinity_direction_class", ""),
                    "combination_generated": False,
                }
            )

    output.sort(key=lambda row: (str(row["design_track"]), int(row["sequence_index_1based"]), str(row["mutant_residue"])))
    counts = Counter(str(row["qualification_status"]) for row in output)
    track_counts = {
        track: dict(Counter(str(row["qualification_status"]) for row in output if row["design_track"] == track))
        for track in ("affinity", "property")
    }
    combination_ready = [str(row["mutation"]) for row in output if row["qualification_status"] == "combination_ready"]
    return {
        "review_rows": output,
        "facts": {
            "candidate_count": len(output),
            "track_counts": {"affinity": AFFINITY_COUNT, "property": PROPERTY_COUNT},
            "qualification_counts": dict(counts),
            "qualification_counts_by_track": track_counts,
            "combination_ready_mutations": combination_ready,
            "combination_generated": False,
        },
    }


def _qualification(
    track: str,
    source: Mapping[str, object],
    hard_flags: Sequence[str],
    flags: Sequence[str],
) -> tuple[str, str]:
    if hard_flags:
        return "blocked", "hard_sequence_or_disulfide_risk"
    if "adjacent_to_experimental_coordinate_gap" in flags:
        return "blocked_pending_structure", "requires_targeted_complete_loop_structure_review"
    major = {
        "strong_negative_antifold_complex_signal",
        "exposed_hydrophobic_substitution",
        "dense_local_hydrophobic_window",
        "vhh_hallmark_hydrophobization",
        "new_proline_backbone_constraint",
        "cdr_glycine_flexibility_change",
        "ensemble_contact_change",
        "paired_receptor_contact_change",
        "prepared_structure_sensitive",
    }
    if track == "property":
        if str(source["affinity_direction_class"]) == "directionally_adverse":
            return "not_prioritized", "paired_pyrosetta_directionally_adverse"
        if major.intersection(flags):
            return "single_mutant_test_only", "property_signal_with_structural_or_compatibility_risk"
        return "combination_ready", "property_signal_with_affinity_nonadverse_and_no_major_structural_flag"

    if _bool(source["core_module_selected"]):
        if major.intersection(flags) or "increased_oxidation_susceptibility" in flags:
            return "single_mutant_test_only", "affinity_core_with_developability_or_structure_risk"
        return "combination_ready", "affinity_core_without_major_structural_flag"
    support = min(int(source["negative_delta_dG_count"]), int(source["negative_delta_cross_interface_count"]))
    favorable_medians = (
        float(source["delta_dG_separated_median"]) < 0
        and float(source["delta_cross_interface_energy_median"]) < 0
    )
    if support >= 13 and favorable_medians and not major.intersection(flags):
        return "targeted_alternative_review", "noncore_favorable_medians_with_at_least_13_of_20_support"
    return "not_prioritized", "does_not_meet_combination_or_targeted_alternative_review_contract"


def _risk_level(hard_flags: Sequence[str], flags: Sequence[str]) -> str:
    if hard_flags:
        return "blocked"
    if (
        "proline_at_coordinate_gap_boundary" in flags
        or "new_proline_backbone_constraint" in flags
        or "vhh_hallmark_hydrophobization" in flags
    ):
        return "high"
    major = sum(
        flag in {
            "strong_negative_antifold_complex_signal",
            "exposed_hydrophobic_substitution",
            "dense_local_hydrophobic_window",
            "ensemble_contact_change",
            "paired_receptor_contact_change",
            "prepared_structure_sensitive",
            "cdr_glycine_flexibility_change",
        }
        for flag in flags
    )
    if major >= 2:
        return "high"
    if major == 1 or "adjacent_to_experimental_coordinate_gap" in flags:
        return "medium_high"
    if flags:
        return "medium"
    return "low"


def _structure_features(
    pdb_path: Path,
    requested: Sequence[str],
    candidates: Mapping[str, Mapping[str, object]],
) -> dict[int, dict[str, object]]:
    if not pdb_path.is_file():
        raise SingleMutantSafetyError(f"Prepared WT PDB is missing: {pdb_path}")
    structure = PDBParser(QUIET=True).get_structure("prepared_wt", str(pdb_path))
    model = next(structure.get_models())
    if "C" not in model or "R" not in model:
        raise SingleMutantSafetyError("Prepared WT must contain VHH chain C and receptor chain R")
    vhh = model["C"]
    receptor_atoms = [atom for atom in model["R"].get_atoms() if atom.element != "H"]
    positions = sorted({int(candidates[identifier]["sequence_index_1based"]) for identifier in requested})
    residues = {}
    for position in positions:
        residue = next((item for item in vhh.get_residues() if item.id[1] == position), None)
        if residue is None:
            raise SingleMutantSafetyError(f"Prepared WT lacks chain C residue {position}")
        expected = str(next(candidates[identifier]["wt_residue"] for identifier in requested if int(candidates[identifier]["sequence_index_1based"]) == position))
        observed = THREE_TO_ONE.get(residue.resname)
        if observed != expected:
            raise SingleMutantSafetyError(f"Prepared WT residue mismatch at {position}: {observed} != {expected}")
        residues[position] = residue

    surface = ShrakeRupley(n_points=200)
    surface.compute(model, level="R")
    complex_sasa = {position: float(residue.sasa) for position, residue in residues.items()}
    surface.compute(vhh, level="R")
    output: dict[int, dict[str, object]] = {}
    all_vhh = list(vhh.get_residues())
    for position, residue in residues.items():
        sidechain = _heavy_atoms(residue, sidechain=True) or _heavy_atoms(residue)
        receptor_distance, receptor_residue = _nearest_residue(sidechain, model["R"].get_residues())
        neighbors = []
        for other in all_vhh:
            if other is residue:
                continue
            distance, _ = _nearest_residue(sidechain, [other])
            if distance < 5.0:
                neighbors.append((distance, other))
        neighbors.sort(key=lambda item: item[0])
        wt = THREE_TO_ONE[residue.resname]
        output[position] = {
            "complex_sasa_a2": round(complex_sasa[position], 6),
            "vhh_alone_sasa_a2": round(float(residue.sasa), 6),
            "vhh_alone_relative_sasa": round(float(residue.sasa) / MAX_ASA[wt], 6),
            "sasa_buried_by_receptor_a2": round(float(residue.sasa) - complex_sasa[position], 6),
            "minimum_receptor_distance_a": round(receptor_distance, 6),
            "nearest_receptor_residue": f"{THREE_TO_ONE.get(receptor_residue.resname, receptor_residue.resname)}{receptor_residue.id[1]}",
            "vhh_neighbors_within_5a": ";".join(
                f"{THREE_TO_ONE.get(other.resname, other.resname)}{other.id[1]}:{distance:.2f}"
                for distance, other in neighbors[:12]
            ),
        }
    return output


def _nearest_residue(atoms: Sequence[object], residues: Sequence[object]) -> tuple[float, object]:
    best_distance = math.inf
    best_residue = None
    for residue in residues:
        for first in atoms:
            for second in _heavy_atoms(residue):
                distance = float(first - second)
                if distance < best_distance:
                    best_distance = distance
                    best_residue = residue
    if best_residue is None:
        raise SingleMutantSafetyError("Could not calculate a residue distance")
    return best_distance, best_residue


def _heavy_atoms(residue: object, sidechain: bool = False) -> list[object]:
    backbone = {"N", "CA", "C", "O", "OXT"}
    return [
        atom
        for atom in residue.get_atoms()
        if atom.element != "H" and (not sidechain or atom.name not in backbone)
    ]


def _local_hydrophobicity(wt: str, mutant: str, position: int) -> tuple[int, int, tuple[str, str]]:
    start = max(0, position - 4)
    stop = min(len(wt), position + 3)
    wt_window = wt[start:stop]
    mutant_window = mutant[start:stop]
    return (
        sum(residue in HYDROPHOBIC for residue in wt_window),
        sum(residue in HYDROPHOBIC for residue in mutant_window),
        (wt_window, mutant_window),
    )


def _motif_count(sequence: str, pattern: str) -> int:
    return len(re.findall(f"(?=({pattern}))", sequence))


def _optional_float(value: object) -> float | None:
    if value in (None, "", "nan", "NaN"):
        return None
    result = float(value)
    return result if math.isfinite(result) else None


def _tokens(value: object, separator: str) -> set[str]:
    return {token for token in str(value or "").split(separator) if token and token.lower() != "nan"}


def _bool(value: object) -> bool:
    return value is True or str(value).lower() == "true"


def _validate_unique(
    rows: Sequence[Mapping[str, object]], key: str, expected: int, label: str
) -> None:
    if len(rows) != expected or len({str(row[key]) for row in rows}) != expected:
        raise SingleMutantSafetyError(f"Expected {expected} unique {label} rows")
