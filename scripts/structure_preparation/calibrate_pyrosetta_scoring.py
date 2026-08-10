#!/usr/bin/env python3
"""Calibrate one reproducible local-interface PyRosetta scoring protocol."""

from __future__ import annotations

import argparse
import csv
import json
import math
import platform
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from antibody_optimization.file_transaction import (  # noqa: E402
    replace_staged_files,
    validate_file_paths,
)
from antibody_optimization.pyrosetta_import_gate import (  # noqa: E402
    compare_pose_to_source,
    evaluate_breaks,
    load_released_stage_inputs,
)
from antibody_optimization.pyrosetta_scoring_calibration import (  # noqa: E402
    PER_RESIDUE_FIELDS,
    PROTOCOL_ORDER,
    REPLICATE_FIELDS,
    CalibrationThresholds,
    build_calibration_gate,
    choose_representative_replicate,
    load_calibration_inputs,
    render_calibration_svg,
    select_protocol,
    summarize_protocol_rows,
)


OUTPUT_NAMES = {
    "replicates": "protocol_replicate_metrics.csv",
    "per_residue": "wt_per_residue_energy.csv",
    "selection": "selected_scoring_protocol.json",
    "gate": "pyrosetta_scoring_calibration_gate.json",
    "figure": "pyrosetta_scoring_calibration_qc.svg",
    "representative": "selected_wt_prepared.pdb",
}
INIT_OPTIONS = (
    "-missing_density_to_jump true "
    "-detect_disulf true "
    "-ignore_unrecognized_res false "
    "-constant_seed -jran 8102026"
)
SCORE_FUNCTION = "ref2015"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage0-dir", type=Path, required=True)
    parser.add_argument("--structure-baseline-dir", type=Path, required=True)
    parser.add_argument("--import-gate-dir", type=Path, required=True)
    parser.add_argument("--experimental-structure", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=8)
    parser.add_argument("--base-seed", type=int, default=8102026)
    parser.add_argument("--interface-neighborhood-angstrom", type=float, default=8.0)
    parser.add_argument("--contact-cutoff-angstrom", type=float, default=4.0)
    parser.add_argument("--coordinate-constraint-sd-angstrom", type=float, default=0.25)
    parser.add_argument("--minimum-vhh-contact-retention", type=float, default=0.80)
    parser.add_argument(
        "--minimum-receptor-epitope-retention", type=float, default=0.90
    )
    parser.add_argument("--maximum-interface-ca-rmsd-angstrom", type=float, default=0.50)
    parser.add_argument("--maximum-dg-mad-reu", type=float, default=3.0)
    parser.add_argument("--generated-at")
    parser.add_argument("--expected-pyrosetta-version", default="2026.03")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    _validate_args(args)
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(
        timespec="seconds"
    )
    stage0_dir = _project_directory(args.stage0_dir)
    structure_dir = _project_directory(args.structure_baseline_dir)
    import_dir = _project_directory(args.import_gate_dir)
    structure_path = args.experimental_structure.expanduser().resolve(strict=True)
    if structure_path.is_symlink() or not structure_path.is_file():
        raise FileNotFoundError(f"Expected regular experimental structure: {structure_path}")

    calibration_inputs = load_calibration_inputs(
        stage0_dir=stage0_dir,
        import_gate_dir=import_dir,
    )
    structure_inputs = load_released_stage_inputs(
        stage0_dir=stage0_dir,
        structure_baseline_dir=structure_dir,
    )
    if (
        calibration_inputs["stage0_run_id"]
        != structure_inputs["stage0_run_id"]
    ):
        raise ValueError("Calibration and structure inputs refer to different stage-0 runs")
    contract = calibration_inputs["contract"]
    recorded_structure = Path(str(contract["inputs"]["experimental_structure"]["path"]))
    if (PROJECT_ROOT / recorded_structure).resolve(strict=True) != structure_path:
        raise ValueError("Experimental structure differs from the released stage-0 contract")

    output_dir = args.output_dir.expanduser().absolute()
    run_summary = args.run_summary.expanduser().absolute()
    requested_targets = [output_dir / name for name in OUTPUT_NAMES.values()]
    validated = validate_file_paths(
        project_root=PROJECT_ROOT,
        source_paths=[structure_path],
        target_paths=[*requested_targets, run_summary],
    )
    final_paths = dict(zip(OUTPUT_NAMES, validated.target_paths[:-1], strict=True))
    run_summary = validated.target_paths[-1]
    existing = [path for path in [*final_paths.values(), run_summary] if path.exists()]
    if existing:
        raise FileExistsError(
            "Refusing to overwrite existing outputs:\n" + "\n".join(map(str, existing))
        )
    for target in [*final_paths.values(), run_summary]:
        target.parent.mkdir(parents=True, exist_ok=True)

    pyrosetta = _import_and_initialize_pyrosetta()
    version = pyrosetta.version()
    if args.expected_pyrosetta_version not in version:
        raise RuntimeError(
            f"PyRosetta version does not contain {args.expected_pyrosetta_version!r}"
        )
    if platform.python_version_tuple()[:2] != ("3", "10"):
        raise RuntimeError("The pinned PyRosetta environment must use Python 3.10")

    raw_pose = pyrosetta.pose_from_file(str(structure_path))
    _assert_pose_safety(raw_pose, structure_inputs)
    scorefxn = pyrosetta.create_score_function(SCORE_FUNCTION)
    scorefxn(raw_pose)
    raw_contacts = _contact_sets(
        raw_pose,
        chain_a="C",
        chain_b="R",
        cutoff=args.contact_cutoff_angstrom,
    )
    expected_vhh = set(calibration_inputs["vhh_interface_auth_positions"])
    if not expected_vhh.issubset(raw_contacts["chain_a_auth_positions"]):
        raise RuntimeError(
            "Raw PyRosetta contact set loses one or more released VHH interface positions"
        )
    reference_contacts = {
        "chain_a_auth_positions": expected_vhh,
        "chain_b_auth_positions": raw_contacts["chain_b_auth_positions"],
    }
    local_indices = _interface_neighborhood(
        raw_pose,
        chain_a="C",
        chain_b="R",
        cutoff=args.interface_neighborhood_angstrom,
    )
    raw_ca = _ca_coordinates(raw_pose, local_indices)
    raw_cross = _cross_interface_energy(raw_pose, scorefxn, "C", "R")
    raw_total = float(scorefxn(raw_pose))
    raw_separated = _score_separated_without_repack(raw_pose, scorefxn, "C")
    raw_metrics = {
        "total_score": raw_total,
        "dG_separated_fixed_sidechains": raw_total - raw_separated,
        "cross_interface_energy": raw_cross["total"],
        "interface_fa_atr": raw_cross["fa_atr"],
        "interface_fa_rep": raw_cross["fa_rep"],
        "vhh_contact_count": len(raw_contacts["chain_a_auth_positions"]),
        "receptor_epitope_count": len(raw_contacts["chain_b_auth_positions"]),
        "local_repack_residue_count": len(local_indices),
    }

    replicate_rows: list[dict[str, object]] = []
    for protocol_number, protocol in enumerate(PROTOCOL_ORDER):
        for replicate in range(1, args.replicates + 1):
            seed = args.base_seed + protocol_number * 10000 + replicate
            print(
                f"Calibration progress: {protocol} replicate "
                f"{replicate}/{args.replicates} seed={seed}",
                flush=True,
            )
            prepared = _prepare_pose(
                raw_pose,
                scorefxn,
                local_indices=local_indices,
                protocol=protocol,
                seed=seed,
                coordinate_constraint_sd=args.coordinate_constraint_sd_angstrom,
            )
            replicate_rows.append(
                _measure_pose(
                    prepared,
                    scorefxn,
                    structure_inputs=structure_inputs,
                    local_indices=local_indices,
                    raw_ca=raw_ca,
                    raw_contacts=reference_contacts,
                    protocol=protocol,
                    replicate=replicate,
                    seed=seed,
                    contact_cutoff=args.contact_cutoff_angstrom,
                )
            )

    thresholds = CalibrationThresholds(
        minimum_vhh_contact_retention=args.minimum_vhh_contact_retention,
        minimum_receptor_epitope_retention=args.minimum_receptor_epitope_retention,
        maximum_interface_ca_rmsd_angstrom=args.maximum_interface_ca_rmsd_angstrom,
        maximum_dg_mad_reu=args.maximum_dg_mad_reu,
    )
    summaries = summarize_protocol_rows(
        replicate_rows,
        raw_interface_fa_rep=float(raw_metrics["interface_fa_rep"]),
        thresholds=thresholds,
    )
    selected_protocol, _ = select_protocol(summaries)
    representative_replicate = (
        choose_representative_replicate(replicate_rows, protocol=selected_protocol)
        if selected_protocol is not None
        else None
    )
    representative_pose = None
    representative_seed = None
    if selected_protocol is not None and representative_replicate is not None:
        representative_row = next(
            row
            for row in replicate_rows
            if row["protocol"] == selected_protocol
            and int(row["replicate"]) == representative_replicate
        )
        representative_seed = int(representative_row["seed"])
        representative_pose = _prepare_pose(
            raw_pose,
            scorefxn,
            local_indices=local_indices,
            protocol=selected_protocol,
            seed=representative_seed,
            coordinate_constraint_sd=args.coordinate_constraint_sd_angstrom,
        )
        _assert_pose_safety(representative_pose, structure_inputs)

    gate = build_calibration_gate(
        generated_at=generated_at,
        pyrosetta_version=version,
        score_function=SCORE_FUNCTION,
        thresholds=thresholds,
        raw_metrics=raw_metrics,
        protocol_summaries=summaries,
        selected_protocol=selected_protocol,
        representative_replicate=representative_replicate,
        stage0_run_id=str(calibration_inputs["stage0_run_id"]),
        import_gate_run_id=str(calibration_inputs["import_gate_run_id"]),
    )
    selection = {
        "schema_version": 1,
        "status": gate["status"],
        "selected_protocol": selected_protocol,
        "representative_replicate": representative_replicate,
        "representative_seed": representative_seed,
        "score_function": SCORE_FUNCTION,
        "protocol_order": list(PROTOCOL_ORDER),
        "selection_rule": "first_passing_protocol_prefer_interface_repack",
        "local_interface_definition": {
            "chains": "C_R",
            "heavy_atom_neighborhood_angstrom": args.interface_neighborhood_angstrom,
            "contact_retention_cutoff_angstrom": args.contact_cutoff_angstrom,
            "local_pose_indices": sorted(local_indices),
        },
        "protocol_parameters": {
            "interface_repack": {
                "backbone_movable": False,
                "sidechain_repack": "local_interface_neighborhood_only",
            },
            "interface_repack_constrained_min": {
                "sidechain_repack": "local_interface_neighborhood_only",
                "backbone_and_chi_minimize": "local_interface_neighborhood_only",
                "coordinate_constraint_atoms": "N_CA_C_O",
                "coordinate_constraint_sd_angstrom": (
                    args.coordinate_constraint_sd_angstrom
                ),
            },
        },
        "candidate_scoring_semantics": (
            "paired_relative_interface_signal_against_same_prepared_WT"
        ),
        "explicit_exclusions": [
            "global_relax",
            "missing_region_completion",
            "absolute_total_score_as_affinity",
            "membrane_protein_absolute_stability",
        ],
    }

    raw_per_residue = _per_residue_rows(
        raw_pose,
        scorefxn,
        structure_state="raw_import",
        protocol="raw_import",
        replicate=0,
        local_indices=local_indices,
    )
    representative_rows = (
        _per_residue_rows(
            representative_pose,
            scorefxn,
            structure_state="selected_prepared_WT",
            protocol=str(selected_protocol),
            replicate=int(representative_replicate),
            local_indices=local_indices,
        )
        if representative_pose is not None
        else []
    )

    with tempfile.TemporaryDirectory(prefix=".pyrosetta-calibration-", dir=PROJECT_ROOT) as temp:
        staging = Path(temp)
        staged = {key: staging / name for key, name in OUTPUT_NAMES.items()}
        staged_summary = staging / "run_summary.json"
        _write_csv(staged["replicates"], replicate_rows, REPLICATE_FIELDS)
        _write_csv(
            staged["per_residue"],
            [*raw_per_residue, *representative_rows],
            PER_RESIDUE_FIELDS,
        )
        _write_json(staged["selection"], selection)
        _write_json(staged["gate"], gate)
        render_calibration_svg(gate=gate, path=staged["figure"])
        if representative_pose is not None:
            representative_pose.dump_pdb(str(staged["representative"]))
        summary = {
            "schema_version": 1,
            "status": gate["status"],
            "generated_at": generated_at,
            "elapsed_seconds": round(time.perf_counter() - started, 6),
            "command_argv": [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
            "pyrosetta_version": version,
            "replicate_count_per_protocol": args.replicates,
            "selected_protocol": selected_protocol,
            "source_stage_ids": gate["source_stage_ids"],
            "outputs": {
                key: str(path)
                for key, path in final_paths.items()
                if key != "representative" or representative_pose is not None
            },
        }
        _write_json(staged_summary, summary)
        installs = {
            staged[key]: final_paths[key]
            for key in OUTPUT_NAMES
            if key != "representative" or representative_pose is not None
        }
        installs[staged_summary] = run_summary
        replace_staged_files(
            installs,
            project_root=PROJECT_ROOT,
            protected_source_paths=[structure_path],
        )
    return 0 if gate["status"] == "pass" else 2


def _validate_args(args: argparse.Namespace) -> None:
    if args.replicates < 3:
        raise ValueError("At least three replicates per protocol are required")
    if args.base_seed <= 0:
        raise ValueError("base seed must be positive")
    if args.interface_neighborhood_angstrom <= args.contact_cutoff_angstrom:
        raise ValueError("Interface neighborhood must exceed contact cutoff")
    if args.coordinate_constraint_sd_angstrom <= 0:
        raise ValueError("Coordinate-constraint SD must be positive")


def _import_and_initialize_pyrosetta():
    try:
        import pyrosetta
    except ImportError as exc:
        raise RuntimeError(
            "PyRosetta is required; use /data/software/env/luly25/multi_ligand"
        ) from exc
    pyrosetta.init(INIT_OPTIONS)
    return pyrosetta


def _prepare_pose(
    raw_pose,
    scorefxn,
    *,
    local_indices: set[int],
    protocol: str,
    seed: int,
    coordinate_constraint_sd: float,
):
    import pyrosetta

    pose = raw_pose.clone()
    pyrosetta.rosetta.numeric.random.rg().set_seed(seed)
    _repack_indices(pose, scorefxn, local_indices)
    if protocol == "interface_repack":
        return pose
    if protocol != "interface_repack_constrained_min":
        raise ValueError(f"Unknown calibration protocol: {protocol}")

    from pyrosetta.rosetta.core.id import AtomID
    from pyrosetta.rosetta.core.kinematics import MoveMap
    from pyrosetta.rosetta.core.scoring import ScoreType
    from pyrosetta.rosetta.core.scoring.constraints import CoordinateConstraint
    from pyrosetta.rosetta.core.scoring.func import HarmonicFunc
    from pyrosetta.rosetta.protocols.minimization_packing import MinMover

    anchor_index = next(index for index in range(1, pose.total_residue() + 1) if index not in local_indices)
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


def _repack_indices(pose, scorefxn, indices: set[int]) -> None:
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


def _measure_pose(
    pose,
    scorefxn,
    *,
    structure_inputs,
    local_indices: set[int],
    raw_ca: dict[int, tuple[float, float, float]],
    raw_contacts,
    protocol: str,
    replicate: int,
    seed: int,
    contact_cutoff: float,
) -> dict[str, object]:
    mapping_pass, breaks_pass, disulfide_pass = _pose_safety(pose, structure_inputs)
    total = float(scorefxn(pose))
    contacts = _contact_sets(pose, "C", "R", contact_cutoff)
    cross = _cross_interface_energy(pose, scorefxn, "C", "R")
    separated = _score_separated_with_repack(
        pose,
        scorefxn,
        moving_chain="C",
        repack_indices=local_indices,
        seed=seed + 1000000,
    )
    vhh_retention = _set_retention(
        raw_contacts["chain_a_auth_positions"], contacts["chain_a_auth_positions"]
    )
    receptor_retention = _set_retention(
        raw_contacts["chain_b_auth_positions"], contacts["chain_b_auth_positions"]
    )
    rmsd = _ca_rmsd(raw_ca, _ca_coordinates(pose, local_indices))
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
    status = (
        "pass"
        if mapping_pass and breaks_pass and disulfide_pass and finite
        else "blocked"
    )
    return {
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


def _pose_safety(pose, structure_inputs) -> tuple[bool, bool, bool]:
    from antibody_optimization.pyrosetta_import_gate import ResidueRecord

    pdb_info = pose.pdb_info()
    pose_records = [
        ResidueRecord(
            index=index,
            chain_id=str(pdb_info.chain(index)).strip(),
            auth_seq_id=int(pdb_info.number(index)),
            insertion_code=_normalize_icode(pdb_info.icode(index)),
            residue_name=str(pose.residue(index).name3()).strip().upper(),
        )
        for index in range(1, pose.total_residue() + 1)
    ]
    mapping_pass = not compare_pose_to_source(
        source_residues=structure_inputs["source_residues"],
        pose_residues=pose_records,
    )
    fold_tree = pose.fold_tree()
    fold_tree_cutpoints = {int(value) for value in fold_tree.cutpoints()}
    jump_cutpoints = {
        int(fold_tree.cutpoint_by_jump(number))
        for number in range(1, int(fold_tree.num_jump()) + 1)
    }
    bonded = _bonded_break_pairs(pose, structure_inputs["expected_breaks"])
    _, break_problems = evaluate_breaks(
        expected_breaks=structure_inputs["expected_breaks"],
        fold_tree_cutpoints=fold_tree_cutpoints,
        jump_cutpoints=jump_cutpoints,
        bonded_c_n_pairs=bonded,
    )
    disulfide_pass = _disulfide_is_bonded(pose, "C", {22, 95})
    return mapping_pass, not break_problems, disulfide_pass


def _assert_pose_safety(pose, structure_inputs) -> None:
    flags = _pose_safety(pose, structure_inputs)
    if flags != (True, True, True):
        raise RuntimeError(f"Pose failed mapping/break/disulfide safety: {flags}")


def _bonded_break_pairs(pose, expected_breaks) -> set[tuple[int, int]]:
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


def _disulfide_is_bonded(pose, chain_id: str, auth_positions: set[int]) -> bool:
    pdb_info = pose.pdb_info()
    indices = [
        index
        for index in range(1, pose.total_residue() + 1)
        if str(pdb_info.chain(index)).strip() == chain_id
        and int(pdb_info.number(index)) in auth_positions
    ]
    return len(indices) == 2 and bool(pose.residue(indices[0]).is_bonded(indices[1]))


def _interface_neighborhood(pose, chain_a: str, chain_b: str, cutoff: float) -> set[int]:
    pairs = _contact_pose_indices(pose, chain_a, chain_b, cutoff)
    return {index for pair in pairs for index in pair}


def _contact_sets(pose, chain_a: str, chain_b: str, cutoff: float) -> dict[str, object]:
    pairs = _contact_pose_indices(pose, chain_a, chain_b, cutoff)
    pdb_info = pose.pdb_info()
    minimum = _minimum_interchain_distance(pose, chain_a, chain_b)
    return {
        "chain_a_auth_positions": {int(pdb_info.number(left)) for left, _ in pairs},
        "chain_b_auth_positions": {int(pdb_info.number(right)) for _, right in pairs},
        "minimum_distance": minimum,
    }


def _contact_pose_indices(pose, chain_a: str, chain_b: str, cutoff: float) -> set[tuple[int, int]]:
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
    pairs: set[tuple[int, int]] = set()
    for left in left_indices:
        left_residue = pose.residue(left)
        for right in right_indices:
            right_residue = pose.residue(right)
            if _residues_within(left_residue, right_residue, cutoff_squared):
                pairs.add((left, right))
    return pairs


def _residues_within(left, right, cutoff_squared: float) -> bool:
    for left_atom in range(1, left.nheavyatoms() + 1):
        left_xyz = left.xyz(left_atom)
        for right_atom in range(1, right.nheavyatoms() + 1):
            delta = left_xyz - right.xyz(right_atom)
            if delta.length_squared() < cutoff_squared:
                return True
    return False


def _minimum_interchain_distance(pose, chain_a: str, chain_b: str) -> float:
    pdb_info = pose.pdb_info()
    left_indices = [i for i in range(1, pose.total_residue() + 1) if str(pdb_info.chain(i)).strip() == chain_a]
    right_indices = [i for i in range(1, pose.total_residue() + 1) if str(pdb_info.chain(i)).strip() == chain_b]
    minimum = math.inf
    for left in left_indices:
        left_residue = pose.residue(left)
        for right in right_indices:
            right_residue = pose.residue(right)
            for left_atom in range(1, left_residue.nheavyatoms() + 1):
                for right_atom in range(1, right_residue.nheavyatoms() + 1):
                    distance = (left_residue.xyz(left_atom) - right_residue.xyz(right_atom)).norm()
                    minimum = min(minimum, distance)
    return minimum


def _cross_interface_energy(pose, scorefxn, chain_a: str, chain_b: str) -> dict[str, float]:
    import pyrosetta

    scorefxn(pose)
    pdb_info = pose.pdb_info()
    left_indices = [i for i in range(1, pose.total_residue() + 1) if str(pdb_info.chain(i)).strip() == chain_a]
    right_indices = [i for i in range(1, pose.total_residue() + 1) if str(pdb_info.chain(i)).strip() == chain_b]
    graph = pose.energies().energy_graph()
    weights = scorefxn.weights()
    fa_atr_type = pyrosetta.rosetta.core.scoring.ScoreType.fa_atr
    fa_rep_type = pyrosetta.rosetta.core.scoring.ScoreType.fa_rep
    result = {"total": 0.0, "fa_atr": 0.0, "fa_rep": 0.0}
    for left in left_indices:
        for right in right_indices:
            edge = graph.find_energy_edge(left, right)
            if edge is None:
                continue
            result["total"] += float(edge.dot(weights))
            energy_map = pyrosetta.rosetta.core.scoring.EMapVector()
            edge.fill_energy_map(energy_map)
            result["fa_atr"] += float(energy_map[fa_atr_type] * weights[fa_atr_type])
            result["fa_rep"] += float(energy_map[fa_rep_type] * weights[fa_rep_type])
    return result


def _score_separated_without_repack(pose, scorefxn, moving_chain: str) -> float:
    separated = pose.clone()
    _translate_chain(separated, moving_chain, 1000.0)
    return float(scorefxn(separated))


def _score_separated_with_repack(
    pose,
    scorefxn,
    *,
    moving_chain: str,
    repack_indices: set[int],
    seed: int,
) -> float:
    import pyrosetta

    separated = pose.clone()
    _translate_chain(separated, moving_chain, 1000.0)
    pyrosetta.rosetta.numeric.random.rg().set_seed(seed)
    _repack_indices(separated, scorefxn, repack_indices)
    return float(scorefxn(separated))


def _translate_chain(pose, chain_id: str, distance: float) -> None:
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


def _ca_coordinates(pose, indices: set[int]) -> dict[int, tuple[float, float, float]]:
    result = {}
    for index in indices:
        xyz = pose.residue(index).xyz("CA")
        result[index] = (float(xyz.x), float(xyz.y), float(xyz.z))
    return result


def _ca_rmsd(
    reference: dict[int, tuple[float, float, float]],
    observed: dict[int, tuple[float, float, float]],
) -> float:
    if reference.keys() != observed.keys() or not reference:
        return math.nan
    squared = 0.0
    for index, ref in reference.items():
        obs = observed[index]
        squared += sum((a - b) ** 2 for a, b in zip(ref, obs, strict=True))
    return math.sqrt(squared / len(reference))


def _set_retention(reference: set[int], observed: set[int]) -> float:
    return len(reference & observed) / len(reference) if reference else math.nan


def _per_residue_rows(
    pose,
    scorefxn,
    *,
    structure_state: str,
    protocol: str,
    replicate: int,
    local_indices: set[int],
) -> list[dict[str, object]]:
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
        chain = str(pdb_info.chain(index)).strip()
        region = (
            "local_interface_neighborhood"
            if index in local_indices
            else ("Nb252_other" if chain == "C" else "NK2R_other")
        )
        weighted = {
            name: float(energies[score_type] * weights[score_type])
            for name, score_type in score_types.items()
        }
        total = sum(
            float(energies[score_type] * weights[score_type])
            for score_type in scorefxn.get_nonzero_weighted_scoretypes()
        )
        rows.append(
            {
                "structure_state": structure_state,
                "protocol": protocol,
                "replicate": replicate,
                "pose_index": index,
                "chain_id": chain,
                "auth_seq_id": int(pdb_info.number(index)),
                "residue_name": str(pose.residue(index).name3()).strip().upper(),
                "region": region,
                **weighted,
                "residue_total_score": total,
            }
        )
    return rows


def _project_directory(path: Path) -> Path:
    resolved = path.expanduser().resolve(strict=True)
    if not resolved.is_dir() or resolved.is_symlink():
        raise FileNotFoundError(f"Expected regular project directory: {resolved}")
    try:
        resolved.relative_to(PROJECT_ROOT.resolve(strict=True))
    except ValueError as exc:
        raise ValueError(f"Directory must be inside the project: {resolved}") from exc
    return resolved


def _normalize_icode(value: object) -> str:
    text = str(value or "").strip()
    return "" if text in {"?", "."} else text


def _write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, value: object) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )


if __name__ == "__main__":
    raise SystemExit(main())
