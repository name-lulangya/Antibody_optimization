"""Shared PyRosetta runtime for calibrated Nb252 local-interface scoring.

The functions here implement the exact pose preparation and measurement
semantics selected by the WT calibration workflow.  They operate only inside
the pinned PyRosetta environment and deliberately exclude missing-region
completion, global relaxation, docking, membrane absolute stability, and any
interpretation of Rosetta units as measured affinity.
"""

from __future__ import annotations

import math
import platform
from typing import Mapping

from .pyrosetta_import_gate import compare_pose_to_source, evaluate_breaks
from .pyrosetta_scoring_calibration import energy_edge_map


INIT_OPTIONS = (
    "-missing_density_to_jump true "
    "-detect_disulf true "
    "-ignore_unrecognized_res false "
    "-constant_seed -jran 8102026"
)
SCORE_FUNCTION = "ref2015"

ONE_TO_THREE = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "K": "LYS",
    "L": "LEU",
    "M": "MET",
    "N": "ASN",
    "P": "PRO",
    "Q": "GLN",
    "R": "ARG",
    "S": "SER",
    "T": "THR",
    "V": "VAL",
    "W": "TRP",
    "Y": "TYR",
}


def initialize_pyrosetta(*, expected_version: str = "2026.03"):
    """Initialize and validate the pinned PyRosetta runtime."""

    try:
        import pyrosetta
    except ImportError as exc:
        raise RuntimeError(
            "PyRosetta is required; use /data/software/env/luly25/multi_ligand"
        ) from exc
    pyrosetta.init(INIT_OPTIONS)
    version = pyrosetta.version()
    if expected_version not in version:
        raise RuntimeError(f"PyRosetta version does not contain {expected_version!r}")
    if platform.python_version_tuple()[:2] != ("3", "10"):
        raise RuntimeError("The pinned PyRosetta environment must use Python 3.10")
    return pyrosetta


def prepare_interface_pose(
    starting_pose,
    scorefxn,
    *,
    local_indices: set[int],
    protocol: str,
    seed: int,
    coordinate_constraint_sd: float,
):
    """Apply one calibrated local-interface preparation protocol."""

    import pyrosetta

    pose = starting_pose.clone()
    pyrosetta.rosetta.numeric.random.rg().set_seed(seed)
    repack_indices(pose, scorefxn, local_indices)
    if protocol == "interface_repack":
        return pose
    if protocol != "interface_repack_constrained_min":
        raise ValueError(f"Unknown calibrated protocol: {protocol}")

    from pyrosetta.rosetta.core.id import AtomID
    from pyrosetta.rosetta.core.kinematics import MoveMap
    from pyrosetta.rosetta.core.scoring import ScoreType
    from pyrosetta.rosetta.core.scoring.constraints import CoordinateConstraint
    from pyrosetta.rosetta.core.scoring.func import HarmonicFunc
    from pyrosetta.rosetta.protocols.minimization_packing import MinMover

    anchor_index = next(
        index
        for index in range(1, pose.total_residue() + 1)
        if index not in local_indices
    )
    anchor_atom = AtomID(pose.residue(anchor_index).atom_index("CA"), anchor_index)
    harmonic = HarmonicFunc(0.0, coordinate_constraint_sd)
    for index in sorted(local_indices):
        residue = pose.residue(index)
        for atom_name in ("N", "CA", "C", "O"):
            atom_id = AtomID(residue.atom_index(atom_name), index)
            pose.add_constraint(
                CoordinateConstraint(atom_id, anchor_atom, pose.xyz(atom_id), harmonic)
            )
    move_map = MoveMap()
    move_map.set_bb(False)
    move_map.set_chi(False)
    move_map.set_jump(False)
    for index in local_indices:
        move_map.set_bb(index, True)
        move_map.set_chi(index, True)
    constrained_scorefxn = scorefxn.clone()
    constrained_scorefxn.set_weight(ScoreType.coordinate_constraint, 1.0)
    minimizer = MinMover(
        move_map,
        constrained_scorefxn,
        "lbfgs_armijo_nonmonotone",
        0.01,
        True,
    )
    minimizer.max_iter(200)
    minimizer.apply(pose)
    pose.remove_constraints()
    scorefxn(pose)
    return pose


def mutate_pose_residue(
    pose,
    *,
    chain_id: str,
    auth_seq_id: int,
    insertion_code: str,
    wt_residue: str,
    mutant_residue: str,
) -> tuple[int, str]:
    """Apply exactly one declared source-auth substitution to a Pose."""

    from pyrosetta.rosetta.protocols.simple_moves import MutateResidue

    if mutant_residue not in ONE_TO_THREE or wt_residue not in ONE_TO_THREE:
        raise ValueError("Mutation residues must be standard one-letter amino acids")
    pdb_info = pose.pdb_info()
    matches = [
        index
        for index in range(1, pose.total_residue() + 1)
        if str(pdb_info.chain(index)).strip() == chain_id
        and int(pdb_info.number(index)) == auth_seq_id
        and normalize_icode(pdb_info.icode(index)) == normalize_icode(insertion_code)
    ]
    if len(matches) != 1:
        raise RuntimeError(
            f"Mutation target must map uniquely: {chain_id} {auth_seq_id}{insertion_code!s}"
        )
    index = matches[0]
    observed = str(pose.residue(index).name1()).strip().upper()
    if observed != wt_residue:
        raise RuntimeError(
            f"Mutation WT mismatch at pose {index}: expected {wt_residue}, observed {observed}"
        )
    MutateResidue(index, ONE_TO_THREE[mutant_residue]).apply(pose)
    if str(pose.residue(index).name1()).strip().upper() != mutant_residue:
        raise RuntimeError("PyRosetta mutation did not install the requested residue")
    return index, ONE_TO_THREE[mutant_residue]


def repack_indices(pose, scorefxn, indices: set[int]) -> None:
    """Repack only the supplied pose indices."""

    import pyrosetta

    task = pyrosetta.standard_packer_task(pose)
    task.restrict_to_repacking()
    task.or_include_current(True)
    for index in range(1, pose.total_residue() + 1):
        if index not in indices:
            task.nonconst_residue_task(index).prevent_repacking()
    mover = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(
        scorefxn, task
    )
    mover.apply(pose)


def measure_interface_pose(
    pose,
    scorefxn,
    *,
    structure_inputs,
    local_indices: set[int],
    reference_ca: Mapping[int, tuple[float, float, float]],
    reference_contacts: Mapping[str, set[int]],
    protocol: str,
    replicate: int,
    seed: int,
    contact_cutoff: float,
    allowed_mutations: Mapping[tuple[str, int, str], str] | None = None,
    include_contact_sets: bool = False,
) -> dict[str, object]:
    """Measure one prepared WT or mutant using calibrated semantics."""

    mapping_pass, breaks_pass, disulfide_pass = pose_safety(
        pose,
        structure_inputs,
        allowed_mutations=allowed_mutations,
    )
    total = float(scorefxn(pose))
    contacts = contact_sets(pose, "C", "R", contact_cutoff)
    cross = cross_interface_energy(pose, scorefxn, "C", "R")
    separated = score_separated_with_repack(
        pose,
        scorefxn,
        moving_chain="C",
        repack_pose_indices=local_indices,
        seed=seed + 1000000,
    )
    vhh_retention = set_retention(
        reference_contacts["chain_a_auth_positions"],
        contacts["chain_a_auth_positions"],
    )
    receptor_retention = set_retention(
        reference_contacts["chain_b_auth_positions"],
        contacts["chain_b_auth_positions"],
    )
    rmsd = ca_rmsd(reference_ca, ca_coordinates(pose, local_indices))
    values = [
        total,
        total - separated,
        cross["total"],
        cross["fa_atr"],
        cross["fa_rep"],
        vhh_retention,
        receptor_retention,
        rmsd,
        contacts["minimum_distance"],
    ]
    finite = all(math.isfinite(value) for value in values)
    status = "pass" if mapping_pass and breaks_pass and disulfide_pass and finite else "blocked"
    result = {
        "protocol": protocol,
        "replicate": replicate,
        "seed": seed,
        "total_score": total,
        "dG_separated": total - separated,
        "cross_interface_energy": cross["total"],
        "interface_fa_atr": cross["fa_atr"],
        "interface_fa_rep": cross["fa_rep"],
        "vhh_contact_count": len(contacts["chain_a_auth_positions"]),
        "receptor_epitope_count": len(contacts["chain_b_auth_positions"]),
        "vhh_contact_retention": vhh_retention,
        "receptor_epitope_retention": receptor_retention,
        "interface_ca_rmsd": rmsd,
        "minimum_interchain_distance": contacts["minimum_distance"],
        "mapping_pass": mapping_pass,
        "breaks_pass": breaks_pass,
        "disulfide_pass": disulfide_pass,
        "finite_metrics": finite,
        "status": status,
    }
    if include_contact_sets:
        result.update(
            {
                "vhh_contact_auth_positions": sorted(
                    contacts["chain_a_auth_positions"]
                ),
                "receptor_contact_auth_positions": sorted(
                    contacts["chain_b_auth_positions"]
                ),
            }
        )
    return result


def pose_safety(
    pose,
    structure_inputs,
    *,
    allowed_mutations: Mapping[tuple[str, int, str], str] | None = None,
) -> tuple[bool, bool, bool]:
    """Validate source mapping, expected breaks, and the VHH disulfide."""

    from .pyrosetta_import_gate import ResidueRecord

    expected_mutations = {
        (chain, int(auth), normalize_icode(icode)): ONE_TO_THREE[residue.upper()]
        for (chain, auth, icode), residue in (allowed_mutations or {}).items()
    }
    source_residues = []
    for record in structure_inputs["source_residues"]:
        replacement = expected_mutations.get(
            (record.chain_id, record.auth_seq_id, normalize_icode(record.insertion_code))
        )
        source_residues.append(
            ResidueRecord(
                index=record.index,
                chain_id=record.chain_id,
                auth_seq_id=record.auth_seq_id,
                insertion_code=record.insertion_code,
                residue_name=replacement or record.residue_name,
            )
        )
    if len(expected_mutations) != sum(
        (record.chain_id, record.auth_seq_id, normalize_icode(record.insertion_code))
        in expected_mutations
        for record in structure_inputs["source_residues"]
    ):
        return False, False, False
    pdb_info = pose.pdb_info()
    pose_records = [
        ResidueRecord(
            index=index,
            chain_id=str(pdb_info.chain(index)).strip(),
            auth_seq_id=int(pdb_info.number(index)),
            insertion_code=normalize_icode(pdb_info.icode(index)),
            residue_name=str(pose.residue(index).name3()).strip().upper(),
        )
        for index in range(1, pose.total_residue() + 1)
    ]
    mapping_pass = not compare_pose_to_source(
        source_residues=source_residues,
        pose_residues=pose_records,
    )
    fold_tree = pose.fold_tree()
    fold_tree_cutpoints = {int(value) for value in fold_tree.cutpoints()}
    jump_cutpoints = {
        int(fold_tree.cutpoint_by_jump(number))
        for number in range(1, int(fold_tree.num_jump()) + 1)
    }
    bonded = bonded_break_pairs(pose, structure_inputs["expected_breaks"])
    _, break_problems = evaluate_breaks(
        expected_breaks=structure_inputs["expected_breaks"],
        fold_tree_cutpoints=fold_tree_cutpoints,
        jump_cutpoints=jump_cutpoints,
        bonded_c_n_pairs=bonded,
    )
    disulfide_pass = disulfide_is_bonded(pose, "C", {22, 95})
    return mapping_pass, not break_problems, disulfide_pass


def assert_pose_safety(
    pose,
    structure_inputs,
    *,
    allowed_mutations: Mapping[tuple[str, int, str], str] | None = None,
) -> None:
    flags = pose_safety(pose, structure_inputs, allowed_mutations=allowed_mutations)
    if flags != (True, True, True):
        raise RuntimeError(f"Pose failed mapping/break/disulfide safety: {flags}")


def bonded_break_pairs(pose, expected_breaks) -> set[tuple[int, int]]:
    from pyrosetta.rosetta.core.id import AtomID

    bonded: set[tuple[int, int]] = set()
    conformation = pose.conformation()
    for item in expected_breaks:
        left = pose.residue(item.left.index)
        right = pose.residue(item.right.index)
        c_atom = AtomID(left.atom_index("C"), item.left.index)
        n_atom = AtomID(right.atom_index("N"), item.right.index)
        if conformation.atoms_are_bonded(c_atom, n_atom):
            bonded.add((item.left.index, item.right.index))
    return bonded


def disulfide_is_bonded(pose, chain_id: str, auth_positions: set[int]) -> bool:
    pdb_info = pose.pdb_info()
    indices = [
        index
        for index in range(1, pose.total_residue() + 1)
        if str(pdb_info.chain(index)).strip() == chain_id
        and int(pdb_info.number(index)) in auth_positions
    ]
    return len(indices) == 2 and bool(pose.residue(indices[0]).is_bonded(indices[1]))


def interface_neighborhood(
    pose, chain_a: str, chain_b: str, cutoff: float
) -> set[int]:
    pairs = contact_pose_indices(pose, chain_a, chain_b, cutoff)
    return {index for pair in pairs for index in pair}


def contact_sets(pose, chain_a: str, chain_b: str, cutoff: float) -> dict[str, object]:
    pairs = contact_pose_indices(pose, chain_a, chain_b, cutoff)
    pdb_info = pose.pdb_info()
    return {
        "chain_a_auth_positions": {int(pdb_info.number(left)) for left, _ in pairs},
        "chain_b_auth_positions": {int(pdb_info.number(right)) for _, right in pairs},
        "minimum_distance": minimum_interchain_distance(pose, chain_a, chain_b),
    }


def contact_pose_indices(
    pose, chain_a: str, chain_b: str, cutoff: float
) -> set[tuple[int, int]]:
    pdb_info = pose.pdb_info()
    left_indices = [
        index
        for index in range(1, pose.total_residue() + 1)
        if str(pdb_info.chain(index)).strip() == chain_a
    ]
    right_indices = [
        index
        for index in range(1, pose.total_residue() + 1)
        if str(pdb_info.chain(index)).strip() == chain_b
    ]
    cutoff_squared = cutoff * cutoff
    return {
        (left, right)
        for left in left_indices
        for right in right_indices
        if residues_within(pose.residue(left), pose.residue(right), cutoff_squared)
    }


def residues_within(left, right, cutoff_squared: float) -> bool:
    for left_atom in range(1, left.nheavyatoms() + 1):
        left_xyz = left.xyz(left_atom)
        for right_atom in range(1, right.nheavyatoms() + 1):
            if (left_xyz - right.xyz(right_atom)).length_squared() < cutoff_squared:
                return True
    return False


def minimum_interchain_distance(pose, chain_a: str, chain_b: str) -> float:
    pdb_info = pose.pdb_info()
    left_indices = [
        index
        for index in range(1, pose.total_residue() + 1)
        if str(pdb_info.chain(index)).strip() == chain_a
    ]
    right_indices = [
        index
        for index in range(1, pose.total_residue() + 1)
        if str(pdb_info.chain(index)).strip() == chain_b
    ]
    minimum = math.inf
    for left in left_indices:
        for right in right_indices:
            left_residue = pose.residue(left)
            right_residue = pose.residue(right)
            for left_atom in range(1, left_residue.nheavyatoms() + 1):
                for right_atom in range(1, right_residue.nheavyatoms() + 1):
                    distance = (
                        left_residue.xyz(left_atom) - right_residue.xyz(right_atom)
                    ).norm()
                    minimum = min(minimum, distance)
    return minimum


def auth_minimum_partner_distances(
    pose, chain_id: str, partner_chain_id: str
) -> dict[int, float]:
    pdb_info = pose.pdb_info()
    indices = [
        index
        for index in range(1, pose.total_residue() + 1)
        if str(pdb_info.chain(index)).strip() == chain_id
    ]
    partner_indices = [
        index
        for index in range(1, pose.total_residue() + 1)
        if str(pdb_info.chain(index)).strip() == partner_chain_id
    ]
    result: dict[int, float] = {}
    for index in indices:
        residue = pose.residue(index)
        minimum = math.inf
        for partner_index in partner_indices:
            partner = pose.residue(partner_index)
            for atom in range(1, residue.nheavyatoms() + 1):
                for partner_atom in range(1, partner.nheavyatoms() + 1):
                    minimum = min(
                        minimum,
                        (residue.xyz(atom) - partner.xyz(partner_atom)).norm(),
                    )
        auth_seq_id = int(pdb_info.number(index))
        if not math.isfinite(minimum) or auth_seq_id in result:
            raise RuntimeError(f"Invalid partner distance for {chain_id} auth {auth_seq_id}")
        result[auth_seq_id] = minimum
    return result


def cross_interface_energy(pose, scorefxn, chain_a: str, chain_b: str) -> dict[str, float]:
    import pyrosetta

    scorefxn(pose)
    pdb_info = pose.pdb_info()
    left_indices = [
        index
        for index in range(1, pose.total_residue() + 1)
        if str(pdb_info.chain(index)).strip() == chain_a
    ]
    right_indices = [
        index
        for index in range(1, pose.total_residue() + 1)
        if str(pdb_info.chain(index)).strip() == chain_b
    ]
    graph = pose.energies().energy_graph()
    weights = scorefxn.weights()
    fa_atr = pyrosetta.rosetta.core.scoring.ScoreType.fa_atr
    fa_rep = pyrosetta.rosetta.core.scoring.ScoreType.fa_rep
    result = {"total": 0.0, "fa_atr": 0.0, "fa_rep": 0.0}
    for left in left_indices:
        for right in right_indices:
            edge = graph.find_energy_edge(left, right)
            if edge is None:
                continue
            result["total"] += float(edge.dot(weights))
            energy_map = energy_edge_map(edge)
            result["fa_atr"] += float(energy_map[fa_atr] * weights[fa_atr])
            result["fa_rep"] += float(energy_map[fa_rep] * weights[fa_rep])
    return result


def score_separated_without_repack(pose, scorefxn, moving_chain: str) -> float:
    separated = pose.clone()
    translate_chain(separated, moving_chain, 1000.0)
    return float(scorefxn(separated))


def score_separated_with_repack(
    pose,
    scorefxn,
    *,
    moving_chain: str,
    repack_pose_indices: set[int],
    seed: int,
) -> float:
    import pyrosetta

    separated = pose.clone()
    translate_chain(separated, moving_chain, 1000.0)
    pyrosetta.rosetta.numeric.random.rg().set_seed(seed)
    repack_indices(separated, scorefxn, repack_pose_indices)
    return float(scorefxn(separated))


def translate_chain(pose, chain_id: str, distance: float) -> None:
    from pyrosetta.rosetta.core.id import AtomID
    from pyrosetta.rosetta.numeric import xyzVector_double_t

    pdb_info = pose.pdb_info()
    shift = xyzVector_double_t(distance, 0.0, 0.0)
    for index in range(1, pose.total_residue() + 1):
        if str(pdb_info.chain(index)).strip() != chain_id:
            continue
        residue = pose.residue(index)
        for atom in range(1, residue.natoms() + 1):
            atom_id = AtomID(atom, index)
            pose.set_xyz(atom_id, pose.xyz(atom_id) + shift)


def ca_coordinates(pose, indices: set[int]) -> dict[int, tuple[float, float, float]]:
    result = {}
    for index in indices:
        xyz = pose.residue(index).xyz("CA")
        result[index] = (float(xyz.x), float(xyz.y), float(xyz.z))
    return result


def ca_rmsd(
    reference: Mapping[int, tuple[float, float, float]],
    observed: Mapping[int, tuple[float, float, float]],
) -> float:
    if reference.keys() != observed.keys() or not reference:
        return math.nan
    squared = sum(
        sum((a - b) ** 2 for a, b in zip(ref, observed[index], strict=True))
        for index, ref in reference.items()
    )
    return math.sqrt(squared / len(reference))


def set_retention(reference: set[int], observed: set[int]) -> float:
    return len(reference & observed) / len(reference) if reference else math.nan


def auth_residue_names(pose) -> dict[str, dict[int, str]]:
    pdb_info = pose.pdb_info()
    result: dict[str, dict[int, str]] = {"C": {}, "R": {}}
    for index in range(1, pose.total_residue() + 1):
        chain_id = str(pdb_info.chain(index)).strip()
        if chain_id not in result:
            continue
        auth_seq_id = int(pdb_info.number(index))
        residue_name = str(pose.residue(index).name3()).strip().upper()
        existing = result[chain_id].get(auth_seq_id)
        if existing is not None and existing != residue_name:
            raise RuntimeError(f"Conflicting residue identity at {chain_id} auth {auth_seq_id}")
        result[chain_id][auth_seq_id] = residue_name
    return result


def per_residue_rows(
    pose,
    scorefxn,
    *,
    structure_state: str,
    protocol: str,
    replicate: int,
    local_indices: set[int],
) -> list[dict[str, object]]:
    """Return the calibrated compact per-residue energy table."""

    import pyrosetta

    scorefxn(pose)
    weights = scorefxn.weights()
    score_types = {
        "fa_atr": pyrosetta.rosetta.core.scoring.ScoreType.fa_atr,
        "fa_rep": pyrosetta.rosetta.core.scoring.ScoreType.fa_rep,
        "fa_sol": pyrosetta.rosetta.core.scoring.ScoreType.fa_sol,
        "fa_dun": pyrosetta.rosetta.core.scoring.ScoreType.fa_dun,
    }
    pdb_info = pose.pdb_info()
    rows = []
    for index in range(1, pose.total_residue() + 1):
        energies = pose.energies().residue_total_energies(index)
        chain_id = str(pdb_info.chain(index)).strip()
        if index in local_indices:
            region = "local_interface_neighborhood"
        elif chain_id == "C":
            region = "Nb252_other"
        elif chain_id == "R":
            region = "NK2R_other"
        else:
            region = "other"
        rows.append(
            {
                "structure_state": structure_state,
                "protocol": protocol,
                "replicate": replicate,
                "pose_index": index,
                "chain_id": chain_id,
                "auth_seq_id": int(pdb_info.number(index)),
                "residue_name": str(pose.residue(index).name3()).strip().upper(),
                "region": region,
                **{
                    label: float(energies[score_type] * weights[score_type])
                    for label, score_type in score_types.items()
                },
                "residue_total_score": float(pose.energies().residue_total_energy(index)),
            }
        )
    return rows


def normalize_icode(value: object) -> str:
    text = str(value or "").strip()
    return "" if text in {"?", "."} else text
