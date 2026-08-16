"""Nonredundant structural review for selected Nb252 single mutants.

The existing unified safety table already carries paired-complex PyRosetta or
Flex ddG evidence.  New computation is therefore restricted to nine non-Pro
A23/F30 mutants adjacent to the experimental coordinate gap and uses only the
complete AF3 VHH parent.  Intrinsic hard risks are excluded before runtime
work.  This module does not predict affinity/expression or make combinations.
"""

from __future__ import annotations

import math
import statistics
from collections import Counter
from pathlib import Path
from typing import Mapping, Sequence

from Bio.PDB import MMCIFParser, ShrakeRupley
from Bio.PDB.vectors import calc_dihedral


class TargetedStructureReviewError(ValueError):
    """Raised when targeted-review inputs violate the released contract."""


REVIEW_POOL_COUNT = 30
TARGET_COUNT = 9
REPLICATES = 3
MUTATION_NEIGHBORHOOD_ANGSTROM = 8.0
GAP_BOUNDARY_POSITIONS = {23, 30}
DIRECT_EXCLUSION_REASONS = {
    "F30P": "noncompensable_proline_backbone_risk_without_compelling_benefit",
    "Q5P": "new_proline_constraint_without_compelling_benefit",
    "R45P": "new_proline_constraint_without_compelling_benefit",
}
NONRESOLVABLE_COMBINATION_FLAGS = {
    "strong_negative_antifold_complex_signal",
    "exposed_hydrophobic_substitution",
    "dense_local_hydrophobic_window",
    "vhh_hallmark_hydrophobization",
    "increased_oxidation_susceptibility",
    "cdr_glycine_flexibility_change",
    "new_proline_backbone_constraint",
}

PLAN_FIELDS = [
    "candidate_id", "design_track", "mutation", "sequence_index_1based",
    "wt_residue", "mutant_residue", "region", "sequence",
    "prior_qualification_status", "prior_qualification_reason",
    "prior_expert_risk_level", "hard_risk_flags", "structural_review_flags",
    "review_group", "review_route", "experimental_auth_asym_id",
    "experimental_auth_seq_id", "experimental_insertion_code",
    "experimental_coordinate_status", "af3_auth_asym_id", "af3_auth_seq_id",
    "af3_insertion_code", "af3_coordinate_status", "requires_af3_vhh_branch",
    "af3_parent_relative_sasa", "af3_parent_neighbors_within_5a",
    "af3_parent_phi_degrees", "af3_parent_psi_degrees",
    "prior_complex_evidence_reused", "combination_generated",
]

HARD_EXCLUSION_FIELDS = [
    "candidate_id", "mutation", "sequence_index_1based", "design_track",
    "prior_qualification_status", "hard_exclusion_reason", "next_step_status",
]


def build_targeted_plan(
    safety_rows: Sequence[Mapping[str, object]],
    mapping_rows: Sequence[Mapping[str, object]],
    *,
    af3_cif: Path,
) -> dict[str, object]:
    """Build the 30-row evidence pool, hard exclusions, and nine AF3 jobs."""

    if len(safety_rows) != 80:
        raise TargetedStructureReviewError("Expected the exact 80-row safety review")
    review_pool = [
        row for row in safety_rows
        if str(row["qualification_status"]) in {
            "blocked_pending_structure", "single_mutant_test_only",
            "targeted_alternative_review",
        }
    ]
    if len(review_pool) != REVIEW_POOL_COUNT:
        raise TargetedStructureReviewError(
            f"Expected {REVIEW_POOL_COUNT} review-pool candidates, found {len(review_pool)}"
        )
    if len({str(row["candidate_id"]) for row in review_pool}) != REVIEW_POOL_COUNT:
        raise TargetedStructureReviewError("Review-pool candidate identifiers are not unique")

    hard_exclusions = _hard_exclusions(safety_rows)
    expected_exclusions = {"Y37C", "R45C", "D98C", "E105C", "Q5P", "R45P", "F30P"}
    observed_exclusions = {str(row["mutation"]) for row in hard_exclusions}
    if observed_exclusions != expected_exclusions:
        raise TargetedStructureReviewError(
            f"Unexpected hard-exclusion set: {sorted(observed_exclusions)}"
        )

    evidence_rows: list[dict[str, object]] = []
    for source in review_pool:
        row = dict(source)
        mutation = str(row["mutation"])
        position = int(row["sequence_index_1based"])
        row["prior_complex_evidence_reused"] = True
        row["next_step_status"] = (
            "do_not_advance" if mutation == "F30P"
            else "af3_gap_review" if position in GAP_BOUNDARY_POSITIONS
            else "retain_existing_evidence"
        )
        row["next_step_reason"] = DIRECT_EXCLUSION_REASONS.get(
            mutation,
            "complete_af3_local_gap_review" if row["next_step_status"] == "af3_gap_review"
            else "no_decision_changing_new_rosetta_evidence_expected",
        )
        evidence_rows.append(row)

    selected = [
        row for row in review_pool
        if str(row["qualification_status"]) == "blocked_pending_structure"
        and str(row["mutation"]) != "F30P"
    ]
    if len(selected) != TARGET_COUNT:
        raise TargetedStructureReviewError(
            f"Expected {TARGET_COUNT} non-Pro gap-boundary candidates, found {len(selected)}"
        )

    mapping = _mapping_index(mapping_rows)
    context = structure_parent_context(af3_cif, chain_id="A")
    plan_rows: list[dict[str, object]] = []
    for row in selected:
        position = int(row["sequence_index_1based"])
        exp = mapping[("NK2R-252.pdb", position)]
        af3 = mapping[("fold_2r_252_nomg_model_0.cif", position)]
        if str(exp["coordinate_status"]) != "observed":
            raise TargetedStructureReviewError(
                f"Target mutation itself lacks experimental coordinates: {row['candidate_id']}"
            )
        if str(af3["coordinate_status"]) != "observed":
            raise TargetedStructureReviewError(
                f"Target mutation lacks AF3 coordinates: {row['candidate_id']}"
            )
        plan_rows.append({
            "candidate_id": row["candidate_id"],
            "design_track": row["design_track"],
            "mutation": row["mutation"],
            "sequence_index_1based": position,
            "wt_residue": row["wt_residue"],
            "mutant_residue": row["mutant_residue"],
            "region": row["region"],
            "sequence": row["sequence"],
            "prior_qualification_status": row["qualification_status"],
            "prior_qualification_reason": row["qualification_reason"],
            "prior_expert_risk_level": row["expert_risk_level"],
            "hard_risk_flags": row["hard_risk_flags"],
            "structural_review_flags": row["structural_review_flags"],
            "review_group": "gap_boundary_nonproline",
            "review_route": "af3_complete_vhh_local_repack_only",
            "experimental_auth_asym_id": exp["auth_asym_id"],
            "experimental_auth_seq_id": exp["auth_seq_id"],
            "experimental_insertion_code": exp["insertion_code"],
            "experimental_coordinate_status": exp["coordinate_status"],
            "af3_auth_asym_id": af3["auth_asym_id"],
            "af3_auth_seq_id": af3["auth_seq_id"],
            "af3_insertion_code": af3["insertion_code"],
            "af3_coordinate_status": af3["coordinate_status"],
            "requires_af3_vhh_branch": True,
            **context[position],
            "prior_complex_evidence_reused": True,
            "combination_generated": False,
        })
    plan_rows.sort(key=lambda item: (int(item["sequence_index_1based"]), str(item["mutant_residue"])))
    mutations = {str(row["mutation"]) for row in plan_rows}
    expected_mutations = {
        "A23Q", "A23R", "A23S", "F30A", "F30K", "F30Q", "F30R", "F30S", "F30T",
    }
    if mutations != expected_mutations:
        raise TargetedStructureReviewError(f"Unexpected runtime mutation set: {sorted(mutations)}")
    return {
        "plan_rows": plan_rows,
        "evidence_rows": evidence_rows,
        "hard_exclusion_rows": hard_exclusions,
        "facts": {
            "candidate_count": TARGET_COUNT,
            "review_pool_count": REVIEW_POOL_COUNT,
            "hard_exclusion_count": len(hard_exclusions),
            "review_group_counts": {"gap_boundary_nonproline": TARGET_COUNT},
            "position_count": 2,
            "af3_branch_candidate_count": TARGET_COUNT,
            "prior_complex_evidence_reused_count": REVIEW_POOL_COUNT,
            "combination_generated": False,
        },
    }


def build_runtime_gate(
    plan_rows: Sequence[Mapping[str, object]],
    replicate_rows: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """Validate complete three-replicate AF3 local-review output."""

    expected_ids = {str(row["candidate_id"]) for row in plan_rows}
    grouped: dict[str, list[Mapping[str, object]]] = {}
    for row in replicate_rows:
        grouped.setdefault(str(row["candidate_id"]), []).append(row)
    failures: list[str] = []
    if set(grouped) != expected_ids:
        failures.append("candidate_identity_mismatch")
    for plan in plan_rows:
        identifier = str(plan["candidate_id"])
        rows = grouped.get(identifier, [])
        if len(rows) != REPLICATES:
            failures.append(f"{identifier}:replicate_count")
            continue
        if {int(row["replicate"]) for row in rows} != {1, 2, 3}:
            failures.append(f"{identifier}:replicate_identity")
        for row in rows:
            if str(row.get("status")) != "pass" or not _bool(row.get("af3_branch_pass")):
                failures.append(f"{identifier}:af3_runtime")
    return {
        "status": "pass" if not failures else "blocked",
        "release": "ready_for_targeted_qualification_v2" if not failures else "blocked_runtime_incomplete",
        "candidate_count": len(expected_ids),
        "replicate_count": len(replicate_rows),
        "expected_replicate_count": len(expected_ids) * REPLICATES,
        "failure_reasons": sorted(set(failures)),
        "candidate_filtering_applied_during_scoring": False,
        "combination_generated": False,
    }


def qualify_v2(
    safety_rows: Sequence[Mapping[str, object]],
    plan_rows: Sequence[Mapping[str, object]],
    replicate_rows: Sequence[Mapping[str, object]],
    runtime_gate: Mapping[str, object],
) -> dict[str, object]:
    """Apply hard exclusions and AF3 gap review without erasing soft risks."""

    if runtime_gate.get("status") != "pass":
        raise TargetedStructureReviewError("Runtime gate does not release V2 qualification")
    plan_by_id = {str(row["candidate_id"]): row for row in plan_rows}
    hard_by_mutation = {
        "Y37C": "new_unpaired_cysteine",
        "R45C": "new_unpaired_cysteine",
        "D98C": "new_unpaired_cysteine",
        "E105C": "new_unpaired_cysteine",
        **DIRECT_EXCLUSION_REASONS,
    }
    reps: dict[str, list[Mapping[str, object]]] = {}
    for row in replicate_rows:
        reps.setdefault(str(row["candidate_id"]), []).append(row)

    output: list[dict[str, object]] = []
    for source in safety_rows:
        row = dict(source)
        identifier = str(row["candidate_id"])
        mutation = str(row["mutation"])
        if mutation in hard_by_mutation:
            row["v2_qualification_status"] = "do_not_advance"
            row["v2_qualification_reason"] = hard_by_mutation[mutation]
            row["targeted_runtime_reviewed"] = False
            output.append(row)
            continue
        if identifier not in plan_by_id:
            row["v2_qualification_status"] = row["qualification_status"]
            row["v2_qualification_reason"] = "existing_complex_evidence_retained_without_redundant_rerun"
            row["targeted_runtime_reviewed"] = False
            output.append(row)
            continue
        plan = plan_by_id[identifier]
        summary = summarize_runtime_rows(reps[identifier])
        row["review_group"] = plan["review_group"]
        row["review_route"] = plan["review_route"]
        residual = NONRESOLVABLE_COMBINATION_FLAGS.intersection(
            _tokens(row.get("structural_review_flags"))
        ) - {"new_proline_backbone_constraint"}
        af3_nonadverse = (
            summary["af3_pass_count"] == REPLICATES
            and summary["median_af3_vhh_delta_local_fa_rep"] <= 0.0
            and summary["median_af3_vhh_delta_total_score"] <= 0.0
        )
        if residual:
            status, reason = "single_mutant_test_only", "nonresolvable_sequence_or_developability_risk_persists"
        elif af3_nonadverse:
            status, reason = "combination_ready", "af3_complete_vhh_local_review_nonadverse"
        else:
            status, reason = "single_mutant_test_only", "af3_complete_vhh_local_nonadverse_gate_not_met"
        row.update(summary)
        row["v2_qualification_status"] = status
        row["v2_qualification_reason"] = reason
        row["targeted_runtime_reviewed"] = True
        output.append(row)

    counts = Counter(str(row["v2_qualification_status"]) for row in output)
    return {
        "review_rows": output,
        "facts": {
            "candidate_count": len(output),
            "qualification_counts": dict(counts),
            "combination_ready_mutations": [
                str(row["mutation"]) for row in output
                if row["v2_qualification_status"] == "combination_ready"
            ],
            "do_not_advance_mutations": [
                str(row["mutation"]) for row in output
                if row["v2_qualification_status"] == "do_not_advance"
            ],
            "combination_generated": False,
        },
    }


def summarize_runtime_rows(rows: Sequence[Mapping[str, object]]) -> dict[str, object]:
    """Summarize three paired WT/mutant AF3 VHH local replicates."""

    if len(rows) != REPLICATES or {int(row["replicate"]) for row in rows} != {1, 2, 3}:
        raise TargetedStructureReviewError("Targeted runtime summary requires replicates 1,2,3")
    return {
        "af3_pass_count": sum(_bool(row.get("af3_branch_pass")) for row in rows),
        "median_af3_vhh_delta_total_score": statistics.median(
            float(row["af3_vhh_delta_total_score"]) for row in rows
        ),
        "median_af3_vhh_delta_local_fa_rep": statistics.median(
            float(row["af3_vhh_delta_local_fa_rep"]) for row in rows
        ),
    }


def local_pose_energy_metrics(pose, scorefxn, local_indices: set[int]) -> dict[str, float]:
    """Measure paired total and local weighted fa_rep scores for one Pose."""

    import pyrosetta

    total = float(scorefxn(pose))
    fa_rep = pyrosetta.rosetta.core.scoring.ScoreType.fa_rep
    weight = float(scorefxn.weights()[fa_rep])
    local_rep = sum(
        float(pose.energies().residue_total_energies(index)[fa_rep]) * weight
        for index in local_indices
    )
    if not math.isfinite(total) or not math.isfinite(local_rep):
        raise TargetedStructureReviewError("Non-finite local pose energy metric")
    return {"total_score": total, "local_fa_rep": local_rep}


def structure_parent_context(cif_path: Path, *, chain_id: str) -> dict[int, dict[str, object]]:
    """Describe parent AF3 VHH coordinates without modeling mutant side chains."""

    structure = MMCIFParser(QUIET=True).get_structure("af3", str(cif_path))
    model = next(structure.get_models())
    if chain_id not in model:
        raise TargetedStructureReviewError(f"AF3 chain {chain_id} is absent")
    chain = model[chain_id]
    residues = [residue for residue in chain.get_residues() if residue.id[0] == " "]
    surface = ShrakeRupley(n_points=200)
    surface.compute(chain, level="R")
    output: dict[int, dict[str, object]] = {}
    for index, residue in enumerate(residues):
        position = int(residue.id[1])
        atoms = [atom for atom in residue.get_atoms() if atom.element != "H"]
        neighbors: list[tuple[float, object]] = []
        for other in residues:
            if other is residue:
                continue
            distance = min(
                float(first - second)
                for first in atoms
                for second in other.get_atoms()
                if second.element != "H"
            )
            if distance < 5.0:
                neighbors.append((distance, other))
        neighbors.sort(key=lambda item: item[0])
        phi = _phi(residues, index)
        psi = _psi(residues, index)
        output[position] = {
            "af3_parent_relative_sasa": round(float(residue.sasa) / _max_asa(residue.resname), 6),
            "af3_parent_neighbors_within_5a": ";".join(
                f"{other.resname}{other.id[1]}:{distance:.2f}" for distance, other in neighbors[:12]
            ),
            "af3_parent_phi_degrees": "" if phi is None else round(math.degrees(phi), 6),
            "af3_parent_psi_degrees": "" if psi is None else round(math.degrees(psi), 6),
        }
    return output


def _hard_exclusions(rows: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for row in rows:
        mutation = str(row["mutation"])
        reason = DIRECT_EXCLUSION_REASONS.get(mutation)
        if "new_unpaired_cysteine" in _tokens(row.get("hard_risk_flags")):
            reason = "new_unpaired_cysteine"
        if reason:
            output.append({
                "candidate_id": row["candidate_id"],
                "mutation": mutation,
                "sequence_index_1based": row["sequence_index_1based"],
                "design_track": row["design_track"],
                "prior_qualification_status": row["qualification_status"],
                "hard_exclusion_reason": reason,
                "next_step_status": "do_not_advance",
            })
    return sorted(output, key=lambda row: (int(row["sequence_index_1based"]), str(row["mutation"])))


def _mapping_index(rows: Sequence[Mapping[str, object]]) -> dict[tuple[str, int], Mapping[str, object]]:
    index = {(str(row["source_model_name"]), int(row["sequence_index_1based"])): row for row in rows}
    if len(index) != 256:
        raise TargetedStructureReviewError("Expected 128 mappings for each of two VHH structures")
    return index


def _phi(residues: Sequence[object], index: int) -> float | None:
    if index == 0:
        return None
    try:
        return float(calc_dihedral(residues[index - 1]["C"].get_vector(), residues[index]["N"].get_vector(), residues[index]["CA"].get_vector(), residues[index]["C"].get_vector()))
    except KeyError:
        return None


def _psi(residues: Sequence[object], index: int) -> float | None:
    if index + 1 >= len(residues):
        return None
    try:
        return float(calc_dihedral(residues[index]["N"].get_vector(), residues[index]["CA"].get_vector(), residues[index]["C"].get_vector(), residues[index + 1]["N"].get_vector()))
    except KeyError:
        return None


def _max_asa(name: str) -> float:
    values = {
        "ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167,
        "GLN": 225, "GLU": 223, "GLY": 104, "HIS": 224, "ILE": 197,
        "LEU": 201, "LYS": 236, "MET": 224, "PHE": 240, "PRO": 159,
        "SER": 155, "THR": 172, "TRP": 285, "TYR": 263, "VAL": 174,
    }
    if name not in values:
        raise TargetedStructureReviewError(f"Unsupported residue for relative SASA: {name}")
    return float(values[name])


def _tokens(value: object) -> set[str]:
    return {token for token in str(value or "").split(";") if token}


def _bool(value: object) -> bool:
    return value is True or str(value).lower() == "true"
