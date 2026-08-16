#!/usr/bin/env python3
"""Run nine paired WT/mutant AF3-VHH local reviews in PyRosetta."""

from __future__ import annotations

import argparse
import csv
import json
import platform
import sys
import time
from datetime import datetime
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from antibody_optimization import pyrosetta_runtime as runtime  # noqa: E402
from antibody_optimization.flex_ddg_runtime import locate_mutation_pose_index, mutation_neighborhood  # noqa: E402
from antibody_optimization.targeted_structure_review import (  # noqa: E402
    MUTATION_NEIGHBORHOOD_ANGSTROM,
    REPLICATES,
    TARGET_COUNT,
    build_runtime_gate,
    local_pose_energy_metrics,
)


FIELDS = [
    "candidate_id", "mutation", "review_group", "replicate", "seed",
    "af3_vhh_delta_total_score", "af3_vhh_delta_local_fa_rep",
    "af3_branch_pass", "status", "representative_af3_vhh_pdb",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--calibration-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=REPLICATES)
    parser.add_argument("--base-seed", type=int, default=8163000)
    parser.add_argument("--generated-at")
    parser.add_argument("--expected-pyrosetta-version", default="2026.03")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.replicates != REPLICATES or args.base_seed <= 0:
        raise ValueError("Targeted review requires exactly three replicates and a positive seed")
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    plan_dir = _project_dir(args.plan_dir)
    calibration_dir = _project_dir(args.calibration_dir)
    contract = _json(plan_dir / "targeted_structure_review_contract.json")
    if contract.get("status") != "pass" or contract.get("release") != "ready_for_remote_targeted_structure_review":
        raise ValueError("Targeted structure-review plan is not released")
    candidates = _csv(plan_dir / "targeted_structure_review_candidates.csv")
    if len(candidates) != TARGET_COUNT or any(row["mutant_residue"] == "P" for row in candidates):
        raise ValueError("Expected the exact nine non-Pro AF3 review candidates")
    output_dir = args.output_dir.expanduser().absolute()
    summary_path = args.run_summary.expanduser().absolute()
    if output_dir.exists() or summary_path.exists():
        raise FileExistsError("Targeted review outputs already exist; remove an interrupted run explicitly before rerunning")
    output_dir.mkdir(parents=True)
    structures_dir = output_dir / "representative_structures"
    structures_dir.mkdir()
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    selection = _json(calibration_dir / "selected_scoring_protocol.json")
    calibration_gate = _json(calibration_dir / "pyrosetta_scoring_calibration_gate.json")
    if calibration_gate.get("pyrosetta_affinity_scoring_release") != "pass":
        raise ValueError("PyRosetta calibration is not released")
    protocol = str(selection.get("selected_protocol"))
    if protocol != "interface_repack_constrained_min":
        raise ValueError("Unexpected calibrated protocol")
    coordinate_sd = float(selection["protocol_parameters"][protocol]["coordinate_constraint_sd_angstrom"])

    pyrosetta = runtime.initialize_pyrosetta(expected_version=args.expected_pyrosetta_version)
    scorefxn = pyrosetta.create_score_function(runtime.SCORE_FUNCTION)
    af3_start = pyrosetta.pose_from_file(str(plan_dir / "af3_nb252_parent_for_pyrosetta.pdb"))
    groups: dict[int, list[dict[str, str]]] = {}
    for candidate in candidates:
        groups.setdefault(int(candidate["sequence_index_1based"]), []).append(candidate)

    rows: list[dict[str, object]] = []
    completed = 0
    expected = len(candidates) * REPLICATES
    for position, position_candidates in sorted(groups.items()):
        pose_index = locate_mutation_pose_index(
            af3_start, chain_id="A", auth_seq_id=position, insertion_code=""
        )
        local_indices = set(
            mutation_neighborhood(af3_start, pose_index, MUTATION_NEIGHBORHOOD_ANGSTROM)
        )
        for replicate in range(1, REPLICATES + 1):
            seed = args.base_seed + position * 100 + replicate
            wt_pose = runtime.prepare_interface_pose(
                af3_start,
                scorefxn,
                local_indices=local_indices,
                protocol=protocol,
                seed=seed,
                coordinate_constraint_sd=coordinate_sd,
            )
            wt_metrics = local_pose_energy_metrics(wt_pose, scorefxn, local_indices)
            for candidate in position_candidates:
                completed += 1
                print(f"[{completed}/{expected}] {candidate['candidate_id']} replicate {replicate}/3", flush=True)
                mutant = af3_start.clone()
                runtime.mutate_pose_residue(
                    mutant,
                    chain_id="A",
                    auth_seq_id=position,
                    insertion_code="",
                    wt_residue=candidate["wt_residue"],
                    mutant_residue=candidate["mutant_residue"],
                )
                mutant = runtime.prepare_interface_pose(
                    mutant,
                    scorefxn,
                    local_indices=local_indices,
                    protocol=protocol,
                    seed=seed,
                    coordinate_constraint_sd=coordinate_sd,
                )
                metrics = local_pose_energy_metrics(mutant, scorefxn, local_indices)
                representative = ""
                if replicate == 1:
                    path = structures_dir / f"{candidate['candidate_id']}__af3_vhh.pdb"
                    mutant.dump_pdb(str(path))
                    representative = str(path)
                rows.append({
                    "candidate_id": candidate["candidate_id"],
                    "mutation": candidate["mutation"],
                    "review_group": candidate["review_group"],
                    "replicate": replicate,
                    "seed": seed,
                    "af3_vhh_delta_total_score": metrics["total_score"] - wt_metrics["total_score"],
                    "af3_vhh_delta_local_fa_rep": metrics["local_fa_rep"] - wt_metrics["local_fa_rep"],
                    "af3_branch_pass": True,
                    "status": "pass",
                    "representative_af3_vhh_pdb": representative,
                })

    gate = build_runtime_gate(candidates, rows)
    gate.update({
        "schema_version": 2,
        "generated_at": generated_at,
        "selected_protocol": protocol,
        "score_function": runtime.SCORE_FUNCTION,
        "mutation_neighborhood_angstrom": MUTATION_NEIGHBORHOOD_ANGSTROM,
        "new_computation_context": "AF3_VHH_alone",
        "existing_complex_evidence_reused_without_rerun": True,
    })
    _write_csv(output_dir / "targeted_structure_replicates.csv", rows, FIELDS)
    _write_json(output_dir / "targeted_structure_runtime_gate.json", gate)
    _write_json(summary_path, {
        "schema_version": 2,
        "status": gate["status"],
        "generated_at": generated_at,
        "elapsed_seconds": round(time.perf_counter() - started, 6),
        "python": platform.python_version(),
        "pyrosetta_version": pyrosetta.version(),
        "command_argv": [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
        "candidate_count": len(candidates),
        "replicate_count": len(rows),
        "candidate_filtering_applied_during_scoring": False,
        "combination_generated": False,
        "output_dir": str(output_dir),
    })
    return 0 if gate["status"] == "pass" else 2


def _project_dir(path: Path) -> Path:
    resolved = path.expanduser().resolve(strict=True)
    resolved.relative_to(PROJECT_ROOT.resolve(strict=True))
    if not resolved.is_dir() or resolved.is_symlink():
        raise ValueError(f"Expected regular project directory: {resolved}")
    return resolved


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


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
