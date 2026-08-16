#!/usr/bin/env python3
"""Run the shared pilot or full property-candidate PyRosetta protocol."""

from __future__ import annotations

import argparse
import csv
import json
import platform
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from antibody_optimization import pyrosetta_runtime as runtime  # noqa: E402
from antibody_optimization.affinity_scoring import (  # noqa: E402
    PAIRED_FIELDS, SUMMARY_FIELDS, WT_CONTROL_FIELDS,
    build_paired_row, build_wt_control_row, summarize_paired_rows,
)
from antibody_optimization.file_transaction import replace_staged_files, validate_file_paths  # noqa: E402
from antibody_optimization.flex_ddg_runtime import locate_mutation_pose_index, mutation_neighborhood  # noqa: E402
from antibody_optimization.property_affinity_review import (  # noqa: E402
    MUTATION_NEIGHBORHOOD_ANGSTROM, PILOT_CANDIDATES, POOL_SIZE, REPLICATES,
    PROPERTY_FIELDS, build_run_gate, combine_movable_indices,
)
from antibody_optimization.pyrosetta_import_gate import load_released_stage_inputs  # noqa: E402


EXTRA_RUNTIME_FIELDS = [
    "position_specific_wt_id", "mutation_pose_index",
    "mutation_neighborhood_pose_indices", "combined_movable_pose_indices",
    "combined_movable_residue_count",
]
OUTPUTS = {
    "wt": "property_affinity_wt_controls.csv",
    "paired": "property_affinity_candidate_replicates.csv",
    "summary": "property_affinity_candidate_summary.csv",
    "gate": "property_affinity_scoring_gate.json",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-kind", choices=("pilot", "full_scan"), required=True)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--pilot-result-dir", type=Path)
    parser.add_argument("--stage0-dir", type=Path, required=True)
    parser.add_argument("--structure-baseline-dir", type=Path, required=True)
    parser.add_argument("--calibration-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=REPLICATES)
    parser.add_argument("--base-seed", type=int, default=8162000)
    parser.add_argument("--generated-at")
    parser.add_argument("--expected-pyrosetta-version", default="2026.03")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.replicates != REPLICATES or args.base_seed <= 0:
        raise ValueError("This contract requires exactly 3 replicates and a positive base seed")
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    plan_dir = _project_dir(args.plan_dir)
    stage0_dir = _project_dir(args.stage0_dir)
    structure_dir = _project_dir(args.structure_baseline_dir)
    calibration_dir = _project_dir(args.calibration_dir)
    plan_contract = _json(plan_dir / "property_affinity_review_contract.json")
    if plan_contract.get("status") != "pass" or plan_contract.get("candidate_count") != POOL_SIZE:
        raise ValueError("Property-affinity plan is not released")
    all_candidates = _csv(plan_dir / "property_affinity_review_candidates.csv")
    if len(all_candidates) != POOL_SIZE:
        raise ValueError("Property-affinity candidate count mismatch")
    if args.run_kind == "pilot":
        if args.pilot_result_dir is not None:
            raise ValueError("Pilot mode does not accept --pilot-result-dir")
        candidates = [row for row in all_candidates if row["candidate_id"] in PILOT_CANDIDATES]
    else:
        if args.pilot_result_dir is None:
            raise ValueError("Full scan requires --pilot-result-dir")
        pilot_dir = _project_dir(args.pilot_result_dir)
        pilot_gate = _json(pilot_dir / OUTPUTS["gate"])
        if pilot_gate.get("status") != "pass" or pilot_gate.get("release") != "ready_for_full_property_affinity_scan":
            raise ValueError("Pilot gate does not release the full scan")
        candidates = all_candidates
    expected_count = 6 if args.run_kind == "pilot" else POOL_SIZE
    if len(candidates) != expected_count:
        raise ValueError("Run-kind candidate count mismatch")

    calibration_gate = _json(calibration_dir / "pyrosetta_scoring_calibration_gate.json")
    selection = _json(calibration_dir / "selected_scoring_protocol.json")
    if calibration_gate.get("pyrosetta_affinity_scoring_release") != "pass":
        raise ValueError("PyRosetta calibration is not released")
    protocol = str(selection.get("selected_protocol"))
    if protocol != "interface_repack_constrained_min":
        raise ValueError("Unexpected calibrated protocol")
    parameters = selection["protocol_parameters"][protocol]
    interface_definition = selection["local_interface_definition"]
    interface_indices = {int(value) for value in interface_definition["local_pose_indices"]}
    contact_cutoff = float(interface_definition["contact_retention_cutoff_angstrom"])
    coordinate_sd = float(parameters["coordinate_constraint_sd_angstrom"])

    structure_inputs = load_released_stage_inputs(stage0_dir=stage0_dir, structure_baseline_dir=structure_dir)
    contact_rows = _csv(calibration_dir / "selected_contact_changes.csv")
    reference_contacts = {
        "chain_a_auth_positions": {int(r["auth_seq_id"]) for r in contact_rows if r["chain_id"] == "C" and r["reference_contact"] == "True"},
        "chain_b_auth_positions": {int(r["auth_seq_id"]) for r in contact_rows if r["chain_id"] == "R" and r["reference_contact"] == "True"},
    }
    starting_pdb = calibration_dir / "selected_wt_prepared.pdb"
    sources = [plan_dir / "property_affinity_review_candidates.csv", plan_dir / "property_affinity_review_contract.json", calibration_dir / "pyrosetta_scoring_calibration_gate.json", calibration_dir / "selected_scoring_protocol.json", calibration_dir / "selected_contact_changes.csv", starting_pdb]
    if args.pilot_result_dir is not None:
        sources.append(args.pilot_result_dir / OUTPUTS["gate"])
    output_dir = args.output_dir.expanduser().absolute()
    summary_path = args.run_summary.expanduser().absolute()
    validated = validate_file_paths(project_root=PROJECT_ROOT, source_paths=sources, target_paths=[*[output_dir / name for name in OUTPUTS.values()], summary_path])
    finals = dict(zip(OUTPUTS, validated.target_paths[:-1], strict=True))
    summary_path = validated.target_paths[-1]
    existing = [p for p in [*finals.values(), summary_path] if p.exists()]
    if existing:
        raise FileExistsError("Refusing to overwrite existing outputs:\n" + "\n".join(map(str, existing)))
    for path in [*finals.values(), summary_path]:
        path.parent.mkdir(parents=True, exist_ok=True)

    pyrosetta = runtime.initialize_pyrosetta(expected_version=args.expected_pyrosetta_version)
    starting_pose = pyrosetta.pose_from_file(str(starting_pdb))
    runtime.assert_pose_safety(starting_pose, structure_inputs)
    scorefxn = pyrosetta.create_score_function(runtime.SCORE_FUNCTION)
    reference_ca = runtime.ca_coordinates(starting_pose, interface_indices)
    groups: dict[int, list[dict[str, str]]] = {}
    for candidate in candidates:
        groups.setdefault(int(candidate["sequence_index_1based"]), []).append(candidate)

    wt_rows: list[dict[str, object]] = []
    paired_rows: list[dict[str, object]] = []
    for position_number, (position, position_candidates) in enumerate(sorted(groups.items()), start=1):
        first = position_candidates[0]
        mutation_index = locate_mutation_pose_index(
            starting_pose,
            chain_id=first["experimental_auth_asym_id"],
            auth_seq_id=int(first["experimental_auth_seq_id"]),
            insertion_code=first["experimental_insertion_code"],
        )
        mutation_local = mutation_neighborhood(starting_pose, mutation_index, MUTATION_NEIGHBORHOOD_ANGSTROM)
        movable = set(combine_movable_indices(interface_indices, mutation_local))
        local_text = ";".join(map(str, sorted(mutation_local)))
        movable_text = ";".join(map(str, sorted(movable)))
        print(f"Position {position_number}/{len(groups)}: reported {position}, {len(position_candidates)} substitutions, movable={len(movable)}", flush=True)
        for replicate in range(1, args.replicates + 1):
            seed = args.base_seed + position * 100 + replicate
            wt_id = f"Nb252_WT_pos{position:03d}_rep{replicate:02d}_seed{seed}"
            wt_pose = runtime.prepare_interface_pose(starting_pose, scorefxn, local_indices=movable, protocol=protocol, seed=seed, coordinate_constraint_sd=coordinate_sd)
            wt_metrics = runtime.measure_interface_pose(wt_pose, scorefxn, structure_inputs=structure_inputs, local_indices=interface_indices, reference_ca=reference_ca, reference_contacts=reference_contacts, protocol=protocol, replicate=replicate, seed=seed, contact_cutoff=contact_cutoff, include_contact_sets=True)
            wt_row = build_wt_control_row(replicate=replicate, seed=seed, metrics=wt_metrics)
            wt_row.update({"wt_control_id": wt_id, "sequence_index_1based": position, "mutation_pose_index": mutation_index, "mutation_neighborhood_pose_indices": local_text, "combined_movable_pose_indices": movable_text, "combined_movable_residue_count": len(movable)})
            wt_rows.append(wt_row)
            for candidate_number, candidate in enumerate(position_candidates, start=1):
                print(f"  replicate {replicate}/3 candidate {candidate_number}/{len(position_candidates)} {candidate['candidate_id']}", flush=True)
                mutant_pose = starting_pose.clone()
                allowed = {(candidate["experimental_auth_asym_id"], int(candidate["experimental_auth_seq_id"]), candidate["experimental_insertion_code"]): candidate["mutant_residue"]}
                runtime.mutate_pose_residue(mutant_pose, chain_id=candidate["experimental_auth_asym_id"], auth_seq_id=int(candidate["experimental_auth_seq_id"]), insertion_code=candidate["experimental_insertion_code"], wt_residue=candidate["wt_residue"], mutant_residue=candidate["mutant_residue"])
                runtime.assert_pose_safety(mutant_pose, structure_inputs, allowed_mutations=allowed)
                mutant_pose = runtime.prepare_interface_pose(mutant_pose, scorefxn, local_indices=movable, protocol=protocol, seed=seed, coordinate_constraint_sd=coordinate_sd)
                mutant_metrics = runtime.measure_interface_pose(mutant_pose, scorefxn, structure_inputs=structure_inputs, local_indices=interface_indices, reference_ca=reference_ca, reference_contacts=reference_contacts, protocol=protocol, replicate=replicate, seed=seed, contact_cutoff=contact_cutoff, allowed_mutations=allowed, include_contact_sets=True)
                row = build_paired_row(candidate, replicate=replicate, seed=seed, wt_metrics=wt_metrics, mutant_metrics=mutant_metrics)
                row.update({"wt_control_id": wt_id, "position_specific_wt_id": wt_id, "mutation_pose_index": mutation_index, "mutation_neighborhood_pose_indices": local_text, "combined_movable_pose_indices": movable_text, "combined_movable_residue_count": len(movable)})
                paired_rows.append(row)

    summaries = summarize_paired_rows(paired_rows, expected_replicates=args.replicates)
    evidence_by_id = {row["candidate_id"]: row for row in candidates}
    for row in summaries:
        row.update({field: evidence_by_id[row["candidate_id"]][field] for field in PROPERTY_FIELDS})
        row["pilot_selected"] = evidence_by_id[row["candidate_id"]]["pilot_selected"]
    gate = build_run_gate(run_kind=args.run_kind, declared_candidate_ids=[r["candidate_id"] for r in candidates], declared_positions=list(groups), wt_controls=wt_rows, paired_rows=paired_rows, summaries=summaries, expected_replicates=args.replicates)
    gate.update({"generated_at": generated_at, "selected_protocol": protocol, "score_function": runtime.SCORE_FUNCTION, "mutation_neighborhood_angstrom": MUTATION_NEIGHBORHOOD_ANGSTROM})

    wt_fields = WT_CONTROL_FIELDS + ["sequence_index_1based", "mutation_pose_index", "mutation_neighborhood_pose_indices", "combined_movable_pose_indices", "combined_movable_residue_count"]
    paired_fields = PAIRED_FIELDS + EXTRA_RUNTIME_FIELDS
    summary_fields = SUMMARY_FIELDS + list(PROPERTY_FIELDS) + ["pilot_selected"]
    with tempfile.TemporaryDirectory(prefix=".property-affinity-", dir=PROJECT_ROOT) as temp:
        staging = Path(temp)
        staged = {key: staging / name for key, name in OUTPUTS.items()}
        staged_summary = staging / "run_summary.json"
        _write_csv(staged["wt"], wt_rows, wt_fields); _write_csv(staged["paired"], paired_rows, paired_fields); _write_csv(staged["summary"], summaries, summary_fields); _write_json(staged["gate"], gate)
        _write_json(staged_summary, {"schema_version": 1, "status": gate["status"], "generated_at": generated_at, "elapsed_seconds": round(time.perf_counter() - started, 6), "python": platform.python_version(), "pyrosetta_version": pyrosetta.version(), "command_argv": [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]], "run_kind": args.run_kind, "candidate_count": len(candidates), "position_count": len(groups), "wt_control_count": len(wt_rows), "mutant_evaluation_count": len(paired_rows), "replicate_count": args.replicates, "candidate_filtering_applied_during_scoring": False, "outputs": {key: str(path) for key, path in finals.items()}})
        replace_staged_files({**{staged[key]: finals[key] for key in OUTPUTS}, staged_summary: summary_path}, project_root=PROJECT_ROOT, protected_source_paths=validated.source_paths)
    return 0 if gate["status"] == "pass" else 2


def _project_dir(path: Path) -> Path:
    resolved = path.expanduser().resolve(strict=True); resolved.relative_to(PROJECT_ROOT.resolve(strict=True))
    if not resolved.is_dir() or resolved.is_symlink(): raise ValueError(f"Expected regular project directory: {resolved}")
    return resolved

def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle: return list(csv.DictReader(handle))

def _json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict): raise ValueError(f"Expected JSON object: {path}")
    return value

def _write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n"); writer.writeheader(); writer.writerows(rows)

def _write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")

if __name__ == "__main__":
    raise SystemExit(main())
