#!/usr/bin/env python3
"""Run one paired-backbone Flex ddG pilot or production task."""

from __future__ import annotations

import argparse
import csv
import faulthandler
import json
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
from antibody_optimization.flex_ddg import (  # noqa: E402
    BACKRUB_NEIGHBORHOOD_ANGSTROM,
    BACKRUB_TEMPERATURE,
    BACKRUB_TRIALS,
)
from antibody_optimization.pyrosetta_import_gate import (  # noqa: E402
    load_released_stage_inputs,
)
from antibody_optimization import flex_ddg_runtime as flex_runtime  # noqa: E402
from antibody_optimization import pyrosetta_runtime as shared  # noqa: E402


OUTPUT_NAMES = {
    "result": "task_result.json",
    "energies": "energy_terms.csv",
    "contacts": "contact_qc.csv",
    "backbone": "backrub_backbone.pdb",
    "wt": "wt_final.pdb",
    "mutant": "mutant_final.pdb",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--run-kind", choices=("pilot", "production"), default="pilot")
    parser.add_argument("--task-index", type=int, required=True)
    parser.add_argument("--stage0-dir", type=Path, required=True)
    parser.add_argument("--structure-baseline-dir", type=Path, required=True)
    parser.add_argument("--calibration-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--backrub-trials", type=int, default=BACKRUB_TRIALS)
    parser.add_argument("--backrub-temperature", type=float, default=BACKRUB_TEMPERATURE)
    parser.add_argument(
        "--backrub-neighborhood-angstrom",
        type=float,
        default=BACKRUB_NEIGHBORHOOD_ANGSTROM,
    )
    parser.add_argument("--generated-at")
    parser.add_argument("--expected-pyrosetta-version", default="2026.03")
    parser.add_argument("--check_only", action="store_true")
    return parser.parse_args()


def main() -> int:
    faulthandler.enable()
    args = parse_args()
    if args.task_index < 0:
        raise ValueError("Task index must be nonnegative")
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(
        timespec="seconds"
    )
    plan_dir = _project_directory(args.plan_dir)
    stage0_dir = _project_directory(args.stage0_dir)
    structure_dir = _project_directory(args.structure_baseline_dir)
    calibration_dir = _project_directory(args.calibration_dir)
    prefix = f"flex_ddg_{args.run_kind}"
    plan_path = plan_dir / f"{prefix}_plan.json"
    manifest_path = plan_dir / f"{prefix}_manifest.csv"
    plan = _load_json(plan_path)
    manifest = _load_csv(manifest_path)
    expected_purpose = {
        "pilot": "production_parameter_runtime_and_protocol_feasibility_only",
        "production": "tier_1_2_plus_selected_tier_3_ensemble_affinity_review",
    }[args.run_kind]
    expected_tier_3_decision = args.run_kind == "production"
    if (
        plan.get("status") != "pass"
        or plan.get("purpose") != expected_purpose
        or plan.get("candidate_selection_performed") is not False
        or plan.get("tier_3_scope_decision_performed") is not expected_tier_3_decision
        or args.task_index >= len(manifest)
    ):
        raise ValueError(f"Flex ddG {args.run_kind} plan or task index is invalid")
    task = manifest[args.task_index]
    if int(task["task_index"]) != args.task_index:
        raise ValueError("Manifest task order does not match the task index")
    if args.backrub_trials != int(plan["backrub"]["trials"]):
        raise ValueError("Run must use the planned backrub trial count")
    if args.backrub_temperature != float(plan["backrub"]["temperature_kT"]):
        raise ValueError("Run must use the planned backrub temperature")
    if args.backrub_neighborhood_angstrom != float(
        plan["backrub"]["mutation_neighborhood_angstrom"]
    ):
        raise ValueError("Run must use the planned backrub neighborhood")

    pyrosetta = shared.initialize_pyrosetta(
        expected_version=args.expected_pyrosetta_version
    )
    api = flex_runtime.validate_backrub_api()

    selection = _load_json(calibration_dir / "selected_scoring_protocol.json")
    local_definition = selection["local_interface_definition"]
    local_indices = {
        int(value) for value in local_definition["local_pose_indices"]
    }
    contact_cutoff = float(local_definition["contact_retention_cutoff_angstrom"])
    structure_inputs = load_released_stage_inputs(
        stage0_dir=stage0_dir,
        structure_baseline_dir=structure_dir,
    )
    contact_rows = _load_csv(calibration_dir / "selected_contact_changes.csv")
    reference_contacts = {
        "chain_a_auth_positions": {
            int(row["auth_seq_id"])
            for row in contact_rows
            if row["chain_id"] == "C" and row["reference_contact"] == "True"
        },
        "chain_b_auth_positions": {
            int(row["auth_seq_id"])
            for row in contact_rows
            if row["chain_id"] == "R" and row["reference_contact"] == "True"
        },
    }
    starting_pdb = calibration_dir / "selected_wt_prepared.pdb"
    starting_pose = pyrosetta.pose_from_file(str(starting_pdb))
    shared.assert_pose_safety(starting_pose, structure_inputs)
    scorefxn = pyrosetta.create_score_function(shared.SCORE_FUNCTION)
    if args.check_only:
        smoke_results = []
        seen_candidates = set()
        representative_ids = set(plan.get("precheck_candidate_ids", ()))
        for candidate in manifest:
            candidate_id = candidate["candidate_id"]
            if candidate_id in seen_candidates or (
                representative_ids and candidate_id not in representative_ids
            ):
                continue
            seen_candidates.add(candidate_id)
            mutation_index = flex_runtime.locate_mutation_pose_index(
                starting_pose,
                chain_id=candidate["experimental_auth_asym_id"],
                auth_seq_id=int(candidate["experimental_auth_seq_id"]),
                insertion_code=candidate["experimental_insertion_code"],
            )
            neighborhood = flex_runtime.mutation_neighborhood(
                starting_pose,
                mutation_index,
                args.backrub_neighborhood_angstrom,
            )
            sampled = flex_runtime.sample_backrub(
                starting_pose,
                scorefxn,
                pivot_indices=neighborhood,
                seed=int(candidate["seed"]),
                trials=1,
                temperature=args.backrub_temperature,
            )
            shared.assert_pose_safety(sampled, structure_inputs)
            smoke_results.append(
                {
                    "candidate_id": candidate_id,
                    "neighborhood_residue_count": len(neighborhood),
                    "status": "pass",
                }
            )
        print(
            json.dumps(
                {
                    "status": "pass",
                    "api": api,
                    "real_pose_one_move_smoke": smoke_results,
                },
                sort_keys=True,
            )
        )
        return 0
    reference_ca = shared.ca_coordinates(starting_pose, local_indices)

    started = time.perf_counter()
    result = flex_runtime.run_paired_sample(
        starting_pose=starting_pose,
        scorefxn=scorefxn,
        structure_inputs=structure_inputs,
        candidate=task,
        seed=int(task["seed"]),
        backrub_trials=args.backrub_trials,
        backrub_temperature=args.backrub_temperature,
        backrub_neighborhood_angstrom=args.backrub_neighborhood_angstrom,
        local_interface_indices=local_indices,
        reference_ca=reference_ca,
        reference_contacts=reference_contacts,
        contact_cutoff=contact_cutoff,
    )
    wt = result["wt_metrics"]
    mutant = result["mutant_metrics"]
    status = "pass" if wt["status"] == "pass" and mutant["status"] == "pass" else "blocked"
    energy_rows = _energy_rows(wt, mutant)
    contact_rows_out = _contact_rows(wt, mutant)

    output_dir = args.output_dir.expanduser().absolute()
    validated = validate_file_paths(
        project_root=PROJECT_ROOT,
        source_paths=[
            plan_path,
            manifest_path,
            starting_pdb,
            calibration_dir / "selected_scoring_protocol.json",
            calibration_dir / "selected_contact_changes.csv",
        ],
        target_paths=[output_dir / name for name in OUTPUT_NAMES.values()],
    )
    final_paths = dict(zip(OUTPUT_NAMES, validated.target_paths, strict=True))
    existing = [path for path in final_paths.values() if path.exists()]
    if existing:
        raise FileExistsError(
            "Refusing to overwrite existing task outputs:\n" + "\n".join(map(str, existing))
        )
    for path in final_paths.values():
        path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".flex-ddg-task-", dir=PROJECT_ROOT) as temp:
        staging = Path(temp)
        staged = {key: staging / name for key, name in OUTPUT_NAMES.items()}
        result["backbone_pose"].dump_pdb(str(staged["backbone"]))
        result["wt_pose"].dump_pdb(str(staged["wt"]))
        result["mutant_pose"].dump_pdb(str(staged["mutant"]))
        _write_csv(staged["energies"], energy_rows)
        _write_csv(staged["contacts"], contact_rows_out)
        output_size = sum(path.stat().st_size for key, path in staged.items() if key != "result")
        task_result = _task_result(
            task=task,
            status=status,
            generated_at=generated_at,
            pyrosetta_version=pyrosetta.version(),
            result=result,
            wt=wt,
            mutant=mutant,
            output_size_bytes=output_size,
            total_elapsed_seconds=time.perf_counter() - started,
            backrub_trials=args.backrub_trials,
            run_kind=args.run_kind,
            tier_3_scope_decision_performed=expected_tier_3_decision,
        )
        _write_json(staged["result"], task_result)
        replace_staged_files(
            {staged[key]: final_paths[key] for key in OUTPUT_NAMES},
            project_root=PROJECT_ROOT,
            protected_source_paths=validated.source_paths,
        )
    return 0 if status == "pass" else 2


def _task_result(**kwargs) -> dict[str, object]:
    task = kwargs["task"]
    result = kwargs["result"]
    wt = kwargs["wt"]
    mutant = kwargs["mutant"]
    vhh_ref = set(wt["vhh_contact_auth_positions"])
    receptor_ref = set(wt["receptor_contact_auth_positions"])
    return {
        "schema_version": 1,
        "run_kind": kwargs["run_kind"],
        "task_index": int(task["task_index"]),
        "task_id": task["task_id"],
        "candidate_id": task["candidate_id"],
        "tier": task["tier"],
        "sample_index": int(task["sample_index"]),
        "seed": int(task["seed"]),
        "status": kwargs["status"],
        "generated_at": kwargs["generated_at"],
        "python": platform.python_version(),
        "pyrosetta_version": kwargs["pyrosetta_version"],
        "score_function": shared.SCORE_FUNCTION,
        "total_elapsed_seconds": kwargs["total_elapsed_seconds"],
        "initial_minimization_seconds": result["initial_minimization_seconds"],
        "backrub_seconds": result["backrub_seconds"],
        "wt_branch_seconds": result["wt_branch_seconds"],
        "mutant_branch_seconds": result["mutant_branch_seconds"],
        "measurement_seconds": result["measurement_seconds"],
        "peak_rss_mb": _peak_rss_mb(),
        "output_size_bytes": kwargs["output_size_bytes"],
        "backrub_trials": kwargs["backrub_trials"],
        "backrub_neighborhood_residue_count": len(
            result["backrub_neighborhood_pose_indices"]
        ),
        "backrub_neighborhood_pose_indices": result[
            "backrub_neighborhood_pose_indices"
        ],
        "delta_dG_separated": (
            float(mutant["dG_separated"]) - float(wt["dG_separated"])
        ),
        "delta_cross_interface_energy": (
            float(mutant["cross_interface_energy"])
            - float(wt["cross_interface_energy"])
        ),
        "delta_interface_fa_rep": (
            float(mutant["interface_fa_rep"]) - float(wt["interface_fa_rep"])
        ),
        "candidate_vs_paired_wt_vhh_contact_retention": shared.set_retention(
            vhh_ref, set(mutant["vhh_contact_auth_positions"])
        ),
        "candidate_vs_paired_wt_receptor_epitope_retention": shared.set_retention(
            receptor_ref, set(mutant["receptor_contact_auth_positions"])
        ),
        "wt_mapping_pass": wt["mapping_pass"],
        "wt_breaks_pass": wt["breaks_pass"],
        "wt_disulfide_pass": wt["disulfide_pass"],
        "mutant_mapping_pass": mutant["mapping_pass"],
        "mutant_breaks_pass": mutant["breaks_pass"],
        "mutant_disulfide_pass": mutant["disulfide_pass"],
        "candidate_selection_performed": False,
        "tier_3_scope_decision_performed": kwargs[
            "tier_3_scope_decision_performed"
        ],
    }


def _peak_rss_mb() -> float:
    """Return Linux peak resident memory in MiB for the remote-only runtime."""

    import resource

    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def _energy_rows(wt, mutant):
    fields = (
        "total_score",
        "dG_separated",
        "cross_interface_energy",
        "interface_fa_atr",
        "interface_fa_rep",
    )
    return [
        {"state": state, **{field: metrics[field] for field in fields}}
        for state, metrics in (("WT", wt), ("mutant", mutant))
    ]


def _contact_rows(wt, mutant):
    rows = []
    for side, field in (
        ("Nb252_VHH", "vhh_contact_auth_positions"),
        ("NK2R", "receptor_contact_auth_positions"),
    ):
        wt_set, mutant_set = set(wt[field]), set(mutant[field])
        for position in sorted(wt_set | mutant_set):
            rows.append(
                {
                    "molecule_side": side,
                    "auth_seq_id": position,
                    "wt_contact": position in wt_set,
                    "mutant_contact": position in mutant_set,
                }
            )
    return rows


def _project_directory(path: Path) -> Path:
    resolved = path.expanduser().resolve(strict=True)
    resolved.relative_to(PROJECT_ROOT.resolve(strict=True))
    if resolved.is_symlink() or not resolved.is_dir():
        raise ValueError(f"Expected regular project directory: {resolved}")
    return resolved


def _load_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _load_json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
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
