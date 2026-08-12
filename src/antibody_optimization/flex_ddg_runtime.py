"""PyRosetta runtime for one Nb252 Flex ddG timing-pilot task.

The implementation follows the published Flex ddG ordering: constrain and
minimize the WT complex, sample a WT local-backbone ensemble with backrub,
then branch the same sampled backbone into independently packed/minimized WT
and mutant poses.  Existing project runtime functions provide mutation,
packing, interface energies, contacts, RMSD, and structural safety checks.
"""

from __future__ import annotations

import time
from typing import Mapping

from . import pyrosetta_runtime as shared


ATOM_PAIR_CONSTRAINT_SD = 0.5
ATOM_PAIR_CONSTRAINT_CUTOFF = 9.0
BRANCH_MINIMIZER_MAX_ITER = 200


def validate_backrub_api() -> dict[str, bool]:
    """Fail before a long task unless the pinned build exposes required APIs."""

    import pyrosetta

    backrub = pyrosetta.rosetta.protocols.backrub.BackrubMover()
    required = {
        "BackrubMover.set_pivot_residues": hasattr(backrub, "set_pivot_residues"),
        "BackrubMover.set_min_atoms": hasattr(backrub, "set_min_atoms"),
        "BackrubMover.set_max_atoms": hasattr(backrub, "set_max_atoms"),
        "BackrubMover.set_require_mm_bend": hasattr(
            backrub, "set_require_mm_bend"
        ),
        "BackrubMover.set_preserve_detailed_balance": hasattr(
            backrub, "set_preserve_detailed_balance"
        ),
        "BackrubMover.add_mainchain_segments": hasattr(
            backrub, "add_mainchain_segments"
        ),
        "BackrubMover.apply": hasattr(backrub, "apply"),
        "MonteCarlo.boltzmann": hasattr(
            pyrosetta.rosetta.protocols.moves.MonteCarlo, "boltzmann"
        ),
    }
    missing = [name for name, available in required.items() if not available]
    if missing:
        raise RuntimeError(f"Pinned PyRosetta lacks required backrub APIs: {missing}")
    return required


def run_paired_sample(
    *,
    starting_pose,
    scorefxn,
    structure_inputs,
    candidate: Mapping[str, object],
    seed: int,
    backrub_trials: int,
    backrub_temperature: float,
    backrub_neighborhood_angstrom: float,
    local_interface_indices: set[int],
    reference_ca: Mapping[int, tuple[float, float, float]],
    reference_contacts: Mapping[str, set[int]],
    contact_cutoff: float,
) -> dict[str, object]:
    """Generate one WT backbone and score paired WT/mutant branches."""

    if backrub_trials <= 0 or backrub_temperature <= 0:
        raise ValueError("Backrub trials and temperature must be positive")
    chain_id = str(candidate["experimental_auth_asym_id"])
    auth_seq_id = int(candidate["experimental_auth_seq_id"])
    insertion_code = str(candidate.get("experimental_insertion_code", ""))
    wt_residue = str(candidate["wt_residue"])
    mutant_residue = str(candidate["mutant_residue"])
    mutation_index = _pose_index(
        starting_pose,
        chain_id=chain_id,
        auth_seq_id=auth_seq_id,
        insertion_code=insertion_code,
    )
    backrub_indices = mutation_neighborhood(
        starting_pose, mutation_index, backrub_neighborhood_angstrom
    )

    started = time.perf_counter()
    phase = time.perf_counter()
    minimized = constrained_minimize(
        starting_pose,
        scorefxn,
        movable_indices=backrub_indices,
        max_iter=BRANCH_MINIMIZER_MAX_ITER,
    )
    initial_minimization_seconds = time.perf_counter() - phase
    shared.assert_pose_safety(minimized, structure_inputs)

    phase = time.perf_counter()
    backbone = sample_backrub(
        minimized,
        scorefxn,
        pivot_indices=backrub_indices,
        seed=seed,
        trials=backrub_trials,
        temperature=backrub_temperature,
    )
    backrub_seconds = time.perf_counter() - phase
    shared.assert_pose_safety(backbone, structure_inputs)

    phase = time.perf_counter()
    wt_pose = optimize_branch(
        backbone,
        scorefxn,
        pack_indices=backrub_indices,
        minimize_indices=backrub_indices,
        seed=seed + 100_000,
    )
    wt_branch_seconds = time.perf_counter() - phase
    shared.assert_pose_safety(wt_pose, structure_inputs)

    phase = time.perf_counter()
    mutant_pose = backbone.clone()
    shared.mutate_pose_residue(
        mutant_pose,
        chain_id=chain_id,
        auth_seq_id=auth_seq_id,
        insertion_code=insertion_code,
        wt_residue=wt_residue,
        mutant_residue=mutant_residue,
    )
    allowed_mutations = {(chain_id, auth_seq_id, insertion_code): mutant_residue}
    shared.assert_pose_safety(
        mutant_pose, structure_inputs, allowed_mutations=allowed_mutations
    )
    mutant_pose = optimize_branch(
        mutant_pose,
        scorefxn,
        pack_indices=backrub_indices,
        minimize_indices=backrub_indices,
        seed=seed + 100_000,
    )
    mutant_branch_seconds = time.perf_counter() - phase
    shared.assert_pose_safety(
        mutant_pose, structure_inputs, allowed_mutations=allowed_mutations
    )

    phase = time.perf_counter()
    wt_metrics = shared.measure_interface_pose(
        wt_pose,
        scorefxn,
        structure_inputs=structure_inputs,
        local_indices=local_interface_indices,
        reference_ca=reference_ca,
        reference_contacts=reference_contacts,
        protocol="flex_ddg_paired_backbone",
        replicate=int(candidate["sample_index"]),
        seed=seed,
        contact_cutoff=contact_cutoff,
        include_contact_sets=True,
    )
    mutant_metrics = shared.measure_interface_pose(
        mutant_pose,
        scorefxn,
        structure_inputs=structure_inputs,
        local_indices=local_interface_indices,
        reference_ca=reference_ca,
        reference_contacts=reference_contacts,
        protocol="flex_ddg_paired_backbone",
        replicate=int(candidate["sample_index"]),
        seed=seed,
        contact_cutoff=contact_cutoff,
        allowed_mutations=allowed_mutations,
        include_contact_sets=True,
    )
    measurement_seconds = time.perf_counter() - phase
    return {
        "backbone_pose": backbone,
        "wt_pose": wt_pose,
        "mutant_pose": mutant_pose,
        "wt_metrics": wt_metrics,
        "mutant_metrics": mutant_metrics,
        "backrub_neighborhood_pose_indices": sorted(backrub_indices),
        "initial_minimization_seconds": initial_minimization_seconds,
        "backrub_seconds": backrub_seconds,
        "wt_branch_seconds": wt_branch_seconds,
        "mutant_branch_seconds": mutant_branch_seconds,
        "measurement_seconds": measurement_seconds,
        "total_elapsed_seconds": time.perf_counter() - started,
    }


def mutation_neighborhood(pose, mutation_index: int, cutoff: float) -> set[int]:
    """Return the 8-A mutation neighborhood plus sequence-adjacent residues."""

    target = pose.residue(mutation_index)
    target_atom = "CA" if target.name1() == "G" else "CB"
    target_xyz = target.xyz(target_atom)
    cutoff_squared = cutoff * cutoff
    selected: set[int] = set()
    for index in range(1, pose.total_residue() + 1):
        residue = pose.residue(index)
        if not residue.is_protein():
            continue
        atom_name = "CA" if residue.name1() == "G" else "CB"
        if residue.has(atom_name) and residue.xyz(atom_name).distance_squared(target_xyz) <= cutoff_squared:
            selected.add(index)
    pdb_info = pose.pdb_info()
    expanded = set(selected)
    for index in selected:
        for adjacent in (index - 1, index + 1):
            if (
                1 <= adjacent <= pose.total_residue()
                and str(pdb_info.chain(adjacent)).strip()
                == str(pdb_info.chain(index)).strip()
            ):
                expanded.add(adjacent)
    if mutation_index not in expanded or len(expanded) < 3:
        raise RuntimeError("Backrub mutation neighborhood is unexpectedly small")
    return expanded


def sample_backrub(
    pose,
    scorefxn,
    *,
    pivot_indices: set[int],
    seed: int,
    trials: int,
    temperature: float,
):
    """Run one Monte Carlo backrub trajectory and return its final accepted pose."""

    import pyrosetta

    pyrosetta.rosetta.numeric.random.rg().set_seed(seed)
    sampled = pose.clone()
    pivot_vector = pyrosetta.rosetta.utility.vector1_unsigned_long()
    for index in sorted(pivot_indices):
        pivot_vector.append(index)
    mover = pyrosetta.rosetta.protocols.backrub.BackrubMover()
    mover.set_pivot_residues(pivot_vector)
    mover.set_preserve_detailed_balance(True)
    # Rosetta counts main-chain atoms; 3..34 atoms approximates the published
    # short 3..12-residue segment range used by the Flex ddG protocol.
    mover.set_min_atoms(3)
    mover.set_max_atoms(34)
    mover.set_require_mm_bend(False)
    segment_count = mover.add_mainchain_segments(sampled)
    if segment_count <= 0:
        raise RuntimeError("No valid backrub segments were generated")
    monte_carlo = pyrosetta.rosetta.protocols.moves.MonteCarlo(
        sampled, scorefxn, temperature
    )
    progress_interval = max(1, trials // 20)
    for trial in range(1, trials + 1):
        mover.apply(sampled)
        monte_carlo.boltzmann(sampled)
        if trial % progress_interval == 0 or trial == trials:
            print(f"  backrub {trial}/{trials}", flush=True)
    scorefxn(sampled)
    return sampled


def optimize_branch(
    pose,
    scorefxn,
    *,
    pack_indices: set[int],
    minimize_indices: set[int],
    seed: int,
):
    """Independently repack and constrained-minimize one sequence branch."""

    import pyrosetta

    optimized = pose.clone()
    pyrosetta.rosetta.numeric.random.rg().set_seed(seed)
    shared.repack_indices(optimized, scorefxn, pack_indices)
    return constrained_minimize(
        optimized,
        scorefxn,
        movable_indices=minimize_indices,
        max_iter=BRANCH_MINIMIZER_MAX_ITER,
    )


def constrained_minimize(pose, scorefxn, *, movable_indices: set[int], max_iter: int):
    """Minimize backbone/chi under current-conformation C-alpha pair restraints."""

    from pyrosetta.rosetta.core.id import AtomID
    from pyrosetta.rosetta.core.kinematics import MoveMap
    from pyrosetta.rosetta.core.scoring import ScoreType
    from pyrosetta.rosetta.core.scoring.constraints import AtomPairConstraint
    from pyrosetta.rosetta.core.scoring.func import HarmonicFunc
    from pyrosetta.rosetta.protocols.minimization_packing import MinMover

    minimized = pose.clone()
    ca_atoms = []
    for index in range(1, minimized.total_residue() + 1):
        residue = minimized.residue(index)
        if residue.is_protein() and residue.has("CA"):
            ca_atoms.append(
                (index, AtomID(residue.atom_index("CA"), index), residue.xyz("CA"))
            )
    cutoff_squared = ATOM_PAIR_CONSTRAINT_CUTOFF**2
    for left_number, (left_index, left_atom, left_xyz) in enumerate(ca_atoms):
        for right_index, right_atom, right_xyz in ca_atoms[left_number + 1 :]:
            if right_xyz.distance_squared(left_xyz) > cutoff_squared:
                continue
            distance = left_xyz.distance(right_xyz)
            minimized.add_constraint(
                AtomPairConstraint(
                    left_atom,
                    right_atom,
                    HarmonicFunc(distance, ATOM_PAIR_CONSTRAINT_SD),
                )
            )
    move_map = MoveMap()
    move_map.set_bb(False)
    move_map.set_chi(False)
    move_map.set_jump(False)
    for index in movable_indices:
        move_map.set_bb(index, True)
        move_map.set_chi(index, True)
    constrained = scorefxn.clone()
    constrained.set_weight(ScoreType.atom_pair_constraint, 1.0)
    minimizer = MinMover(
        move_map,
        constrained,
        "lbfgs_armijo_nonmonotone",
        0.000001,
        True,
    )
    minimizer.max_iter(max_iter)
    minimizer.apply(minimized)
    minimized.remove_constraints()
    scorefxn(minimized)
    return minimized


def _pose_index(
    pose, *, chain_id: str, auth_seq_id: int, insertion_code: str
) -> int:
    pdb_info = pose.pdb_info()
    matches = [
        index
        for index in range(1, pose.total_residue() + 1)
        if str(pdb_info.chain(index)).strip() == chain_id
        and int(pdb_info.number(index)) == auth_seq_id
        and shared.normalize_icode(pdb_info.icode(index))
        == shared.normalize_icode(insertion_code)
    ]
    if len(matches) != 1:
        raise RuntimeError(
            f"Mutation target must map uniquely: {chain_id} {auth_seq_id}{insertion_code}"
        )
    return matches[0]
