#!/usr/bin/env python3
"""Calibrate one reproducible local-interface PyRosetta scoring protocol."""

from __future__ import annotations

import argparse
import csv
import json
import math
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
    load_released_stage_inputs,
)
from antibody_optimization import pyrosetta_runtime as runtime  # noqa: E402
from antibody_optimization.pyrosetta_scoring_calibration import (  # noqa: E402
    CONTACT_CHANGE_FIELDS,
    PER_RESIDUE_FIELDS,
    PROTOCOL_ORDER,
    REPLICATE_FIELDS,
    CalibrationThresholds,
    audit_source_incomplete_sidechains,
    build_calibration_gate,
    build_contact_change_rows,
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
    "contact_changes": "selected_contact_changes.csv",
    "gate": "pyrosetta_scoring_calibration_gate.json",
    "figure": "pyrosetta_scoring_calibration_qc.svg",
    "representative": "selected_wt_prepared.pdb",
}
SCORE_FUNCTION = runtime.SCORE_FUNCTION


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
    source_incomplete_sidechains = audit_source_incomplete_sidechains(
        structure_baseline_dir=structure_dir,
        vhh_interface_auth_positions=calibration_inputs[
            "vhh_interface_auth_positions"
        ],
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

    pyrosetta = runtime.initialize_pyrosetta(
        expected_version=args.expected_pyrosetta_version
    )
    version = pyrosetta.version()

    raw_pose = pyrosetta.pose_from_file(str(structure_path))
    runtime.assert_pose_safety(raw_pose, structure_inputs)
    scorefxn = pyrosetta.create_score_function(SCORE_FUNCTION)
    scorefxn(raw_pose)
    raw_contacts = runtime.contact_sets(
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
    reference_contact_distances = {
        "C": runtime.auth_minimum_partner_distances(raw_pose, "C", "R"),
        "R": runtime.auth_minimum_partner_distances(raw_pose, "R", "C"),
    }
    local_indices = runtime.interface_neighborhood(
        raw_pose,
        chain_a="C",
        chain_b="R",
        cutoff=args.interface_neighborhood_angstrom,
    )
    source_incomplete_sidechains = [
        {
            **row,
            "local_interface_neighborhood": int(row["pose_index"]) in local_indices,
        }
        for row in source_incomplete_sidechains
    ]
    raw_ca = runtime.ca_coordinates(raw_pose, local_indices)
    raw_cross = runtime.cross_interface_energy(raw_pose, scorefxn, "C", "R")
    raw_total = float(scorefxn(raw_pose))
    raw_separated = runtime.score_separated_without_repack(raw_pose, scorefxn, "C")
    raw_metrics = {
        "total_score": raw_total,
        "dG_separated_fixed_sidechains": raw_total - raw_separated,
        "cross_interface_energy": raw_cross["total"],
        "interface_fa_atr": raw_cross["fa_atr"],
        "interface_fa_rep": raw_cross["fa_rep"],
        "vhh_contact_count": len(raw_contacts["chain_a_auth_positions"]),
        "receptor_epitope_count": len(raw_contacts["chain_b_auth_positions"]),
        "local_repack_residue_count": len(local_indices),
        "source_incomplete_sidechain_residue_count": len(
            source_incomplete_sidechains
        ),
        "source_missing_heavy_atom_count": sum(
            int(row["missing_heavy_atom_count"])
            for row in source_incomplete_sidechains
        ),
        "source_incomplete_sidechain_local_residue_count": sum(
            bool(row["local_interface_neighborhood"])
            for row in source_incomplete_sidechains
        ),
        "source_incomplete_vhh_interface_positions": [
            int(row["auth_seq_id"])
            for row in source_incomplete_sidechains
            if row["vhh_experimental_interface"]
        ],
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
            prepared = runtime.prepare_interface_pose(
                raw_pose,
                scorefxn,
                local_indices=local_indices,
                protocol=protocol,
                seed=seed,
                coordinate_constraint_sd=args.coordinate_constraint_sd_angstrom,
            )
            replicate_rows.append(
                runtime.measure_interface_pose(
                    prepared,
                    scorefxn,
                    structure_inputs=structure_inputs,
                    local_indices=local_indices,
                    reference_ca=raw_ca,
                    reference_contacts=reference_contacts,
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
    representative_row = None
    if selected_protocol is not None and representative_replicate is not None:
        representative_row = next(
            row
            for row in replicate_rows
            if row["protocol"] == selected_protocol
            and int(row["replicate"]) == representative_replicate
        )
        representative_seed = int(representative_row["seed"])
        representative_pose = runtime.prepare_interface_pose(
            raw_pose,
            scorefxn,
            local_indices=local_indices,
            protocol=selected_protocol,
            seed=representative_seed,
            coordinate_constraint_sd=args.coordinate_constraint_sd_angstrom,
        )
        runtime.assert_pose_safety(representative_pose, structure_inputs)

    contact_change_rows: list[dict[str, object]] = []
    if representative_pose is not None and representative_row is not None:
        prepared_contacts = runtime.contact_sets(
            representative_pose,
            chain_a="C",
            chain_b="R",
            cutoff=args.contact_cutoff_angstrom,
        )
        residue_names = runtime.auth_residue_names(representative_pose)
        prepared_contact_distances = {
            "C": runtime.auth_minimum_partner_distances(representative_pose, "C", "R"),
            "R": runtime.auth_minimum_partner_distances(representative_pose, "R", "C"),
        }
        contact_change_rows.extend(
            build_contact_change_rows(
                molecule_side="Nb252_VHH",
                chain_id="C",
                reference_positions=reference_contacts["chain_a_auth_positions"],
                prepared_positions=prepared_contacts["chain_a_auth_positions"],
                residue_names=residue_names["C"],
                reference_minimum_distances=reference_contact_distances["C"],
                prepared_minimum_distances=prepared_contact_distances["C"],
            )
        )
        contact_change_rows.extend(
            build_contact_change_rows(
                molecule_side="NK2R",
                chain_id="R",
                reference_positions=reference_contacts["chain_b_auth_positions"],
                prepared_positions=prepared_contacts["chain_b_auth_positions"],
                residue_names=residue_names["R"],
                reference_minimum_distances=reference_contact_distances["R"],
                prepared_minimum_distances=prepared_contact_distances["R"],
            )
        )
        _assert_representative_contact_metrics(
            representative_row=representative_row,
            contact_change_rows=contact_change_rows,
        )

    gate = build_calibration_gate(
        generated_at=generated_at,
        pyrosetta_version=version,
        score_function=SCORE_FUNCTION,
        thresholds=thresholds,
        raw_metrics=raw_metrics,
        protocol_summaries=summaries,
        selected_protocol=selected_protocol,
        representative_replicate=representative_replicate,
        contact_change_rows=contact_change_rows,
        stage0_run_id=str(calibration_inputs["stage0_run_id"]),
        import_gate_run_id=str(calibration_inputs["import_gate_run_id"]),
    )
    selection = {
        "schema_version": 2,
        "status": gate["status"],
        "selected_protocol": selected_protocol,
        "representative_replicate": representative_replicate,
        "representative_seed": representative_seed,
        "score_function": SCORE_FUNCTION,
        "protocol_order": list(PROTOCOL_ORDER),
        "selection_rule": "first_protocol_passing_all_structure_reproducibility_and_negative_binding_direction_gates",
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
        "source_incomplete_sidechains": source_incomplete_sidechains,
    }

    raw_per_residue = runtime.per_residue_rows(
        raw_pose,
        scorefxn,
        structure_state="raw_import",
        protocol="raw_import",
        replicate=0,
        local_indices=local_indices,
    )
    representative_rows = (
        runtime.per_residue_rows(
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
        _write_csv(
            staged["contact_changes"],
            contact_change_rows,
            CONTACT_CHANGE_FIELDS,
        )
        _write_json(staged["gate"], gate)
        render_calibration_svg(gate=gate, path=staged["figure"])
        if representative_pose is not None:
            representative_pose.dump_pdb(str(staged["representative"]))
        summary = {
            "schema_version": 2,
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


def _assert_representative_contact_metrics(
    *,
    representative_row: dict[str, object],
    contact_change_rows: list[dict[str, object]],
) -> None:
    for side, field in (
        ("Nb252_VHH", "vhh_contact_retention"),
        ("NK2R", "receptor_epitope_retention"),
    ):
        rows = [row for row in contact_change_rows if row["molecule_side"] == side]
        reference_count = sum(bool(row["reference_contact"]) for row in rows)
        retained_count = sum(row["contact_status"] == "retained" for row in rows)
        observed = retained_count / reference_count if reference_count else math.nan
        expected = float(representative_row[field])
        if not math.isclose(observed, expected, rel_tol=0.0, abs_tol=1e-12):
            raise RuntimeError(
                f"Representative {side} contact retention {observed} != {expected}"
            )


def _project_directory(path: Path) -> Path:
    resolved = path.expanduser().resolve(strict=True)
    if not resolved.is_dir() or resolved.is_symlink():
        raise FileNotFoundError(f"Expected regular project directory: {resolved}")
    try:
        resolved.relative_to(PROJECT_ROOT.resolve(strict=True))
    except ValueError as exc:
        raise ValueError(f"Directory must be inside the project: {resolved}") from exc
    return resolved


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
