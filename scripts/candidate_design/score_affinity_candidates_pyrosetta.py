#!/usr/bin/env python3
"""Score every declared Nb252 candidate without applying candidate filters."""

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

from antibody_optimization.affinity_scoring import (  # noqa: E402
    PAIRED_FIELDS,
    SUMMARY_FIELDS,
    WT_CONTROL_FIELDS,
    build_paired_row,
    build_pilot_gate,
    build_wt_control_row,
    summarize_paired_rows,
)
from antibody_optimization.file_transaction import (  # noqa: E402
    replace_staged_files,
    validate_file_paths,
)
from antibody_optimization.pyrosetta_import_gate import (  # noqa: E402
    load_released_stage_inputs,
)
from antibody_optimization import pyrosetta_runtime as runtime  # noqa: E402


OUTPUT_NAMES = {
    "wt_controls": "wt_replicate_metrics.csv",
    "paired": "candidate_replicate_metrics.csv",
    "summary": "candidate_summary.csv",
    "gate": "affinity_scoring_pilot_gate.json",
    "figure": "affinity_scoring_pilot_qc.svg",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-dir", type=Path, required=True)
    parser.add_argument("--candidate-id-file", type=Path, required=True)
    parser.add_argument("--stage0-dir", type=Path, required=True)
    parser.add_argument("--structure-baseline-dir", type=Path, required=True)
    parser.add_argument("--calibration-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=3)
    parser.add_argument("--base-seed", type=int, default=8112100)
    parser.add_argument("--generated-at")
    parser.add_argument("--expected-pyrosetta-version", default="2026.03")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.replicates < 3 or args.base_seed <= 0:
        raise ValueError("At least three replicates and a positive base seed are required")
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(
        timespec="seconds"
    )
    candidate_dir = _project_directory(args.candidate_dir)
    stage0_dir = _project_directory(args.stage0_dir)
    structure_dir = _project_directory(args.structure_baseline_dir)
    calibration_dir = _project_directory(args.calibration_dir)
    candidate_id_file = args.candidate_id_file.expanduser().resolve(strict=True)
    if candidate_id_file.parent != candidate_dir:
        raise ValueError("Candidate ID file must belong to the declared candidate directory")

    candidate_gate = _load_json(candidate_dir / "affinity_candidate_gate.json")
    calibration_gate = _load_json(
        calibration_dir / "pyrosetta_scoring_calibration_gate.json"
    )
    selection = _load_json(calibration_dir / "selected_scoring_protocol.json")
    if candidate_gate.get("pyrosetta_pilot_release") != "ready_for_remote_pilot":
        raise ValueError("Candidate manifest does not release the remote pilot")
    if calibration_gate.get("pyrosetta_affinity_scoring_release") != "pass":
        raise ValueError("PyRosetta calibration does not release candidate scoring")
    protocol = str(selection.get("selected_protocol", ""))
    if protocol != "interface_repack_constrained_min":
        raise ValueError("Unexpected selected scoring protocol")
    parameters = selection["protocol_parameters"][protocol]
    local_definition = selection["local_interface_definition"]
    local_indices = {int(value) for value in local_definition["local_pose_indices"]}
    contact_cutoff = float(local_definition["contact_retention_cutoff_angstrom"])
    coordinate_sd = float(parameters["coordinate_constraint_sd_angstrom"])

    all_candidates = _load_csv(candidate_dir / "affinity_single_mutants.csv")
    by_id = {row["candidate_id"]: row for row in all_candidates}
    if len(by_id) != len(all_candidates):
        raise ValueError("Candidate manifest contains duplicate IDs")
    requested_ids = [
        line.strip()
        for line in candidate_id_file.read_text(encoding="utf-8-sig").splitlines()
        if line.strip()
    ]
    if not requested_ids or len(requested_ids) != len(set(requested_ids)):
        raise ValueError("Candidate ID file must contain unique nonempty IDs")
    missing = sorted(set(requested_ids) - set(by_id))
    if missing:
        raise ValueError(f"Unknown candidate IDs: {missing}")
    candidates = [by_id[candidate_id] for candidate_id in requested_ids]

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
    if len(reference_contacts["chain_a_auth_positions"]) != 24:
        raise ValueError("Calibration reference does not retain 24 VHH interface positions")

    starting_pdb = calibration_dir / "selected_wt_prepared.pdb"
    source_paths = [
        candidate_dir / "affinity_single_mutants.csv",
        candidate_dir / "affinity_candidate_gate.json",
        candidate_id_file,
        calibration_dir / "pyrosetta_scoring_calibration_gate.json",
        calibration_dir / "selected_scoring_protocol.json",
        calibration_dir / "selected_contact_changes.csv",
        starting_pdb,
    ]
    output_dir = args.output_dir.expanduser().absolute()
    run_summary = args.run_summary.expanduser().absolute()
    validated = validate_file_paths(
        project_root=PROJECT_ROOT,
        source_paths=source_paths,
        target_paths=[
            *[output_dir / name for name in OUTPUT_NAMES.values()],
            run_summary,
        ],
    )
    final_paths = dict(zip(OUTPUT_NAMES, validated.target_paths[:-1], strict=True))
    run_summary = validated.target_paths[-1]
    existing = [path for path in [*final_paths.values(), run_summary] if path.exists()]
    if existing:
        raise FileExistsError(
            "Refusing to overwrite existing outputs:\n" + "\n".join(map(str, existing))
        )
    for path in [*final_paths.values(), run_summary]:
        path.parent.mkdir(parents=True, exist_ok=True)

    pyrosetta = runtime.initialize_pyrosetta(
        expected_version=args.expected_pyrosetta_version
    )
    version = pyrosetta.version()
    starting_pose = pyrosetta.pose_from_file(str(starting_pdb))
    runtime.assert_pose_safety(starting_pose, structure_inputs)
    scorefxn = pyrosetta.create_score_function(runtime.SCORE_FUNCTION)
    reference_ca = runtime.ca_coordinates(starting_pose, local_indices)

    wt_control_rows = []
    paired_rows = []
    for replicate in range(1, args.replicates + 1):
        seed = args.base_seed + replicate
        print(
            f"Affinity pilot V2: replicate {replicate}/{args.replicates} seed={seed}",
            flush=True,
        )
        wt_pose = runtime.prepare_interface_pose(
            starting_pose,
            scorefxn,
            local_indices=local_indices,
            protocol=protocol,
            seed=seed,
            coordinate_constraint_sd=coordinate_sd,
        )
        wt_metrics = runtime.measure_interface_pose(
            wt_pose,
            scorefxn,
            structure_inputs=structure_inputs,
            local_indices=local_indices,
            reference_ca=reference_ca,
            reference_contacts=reference_contacts,
            protocol=protocol,
            replicate=replicate,
            seed=seed,
            contact_cutoff=contact_cutoff,
            include_contact_sets=True,
        )
        wt_control_rows.append(
            build_wt_control_row(replicate=replicate, seed=seed, metrics=wt_metrics)
        )
        for candidate_number, candidate in enumerate(candidates, start=1):
            print(
                f"  candidate {candidate_number}/{len(candidates)} "
                f"{candidate['candidate_id']}",
                flush=True,
            )
            mutant_pose = starting_pose.clone()
            chain_id = str(candidate["experimental_auth_asym_id"])
            auth_seq_id = int(candidate["experimental_auth_seq_id"])
            insertion_code = str(candidate["experimental_insertion_code"])
            mutant_residue = str(candidate["mutant_residue"])
            runtime.mutate_pose_residue(
                mutant_pose,
                chain_id=chain_id,
                auth_seq_id=auth_seq_id,
                insertion_code=insertion_code,
                wt_residue=str(candidate["wt_residue"]),
                mutant_residue=mutant_residue,
            )
            runtime.assert_pose_safety(
                mutant_pose,
                structure_inputs,
                allowed_mutations={(chain_id, auth_seq_id, insertion_code): mutant_residue},
            )
            mutant_pose = runtime.prepare_interface_pose(
                mutant_pose,
                scorefxn,
                local_indices=local_indices,
                protocol=protocol,
                seed=seed,
                coordinate_constraint_sd=coordinate_sd,
            )
            mutant_metrics = runtime.measure_interface_pose(
                mutant_pose,
                scorefxn,
                structure_inputs=structure_inputs,
                local_indices=local_indices,
                reference_ca=reference_ca,
                reference_contacts=reference_contacts,
                protocol=protocol,
                replicate=replicate,
                seed=seed,
                contact_cutoff=contact_cutoff,
                allowed_mutations={(chain_id, auth_seq_id, insertion_code): mutant_residue},
                include_contact_sets=True,
            )
            paired_rows.append(
                build_paired_row(
                    candidate,
                    replicate=replicate,
                    seed=seed,
                    wt_metrics=wt_metrics,
                    mutant_metrics=mutant_metrics,
                )
            )

    summaries = summarize_paired_rows(paired_rows, expected_replicates=args.replicates)
    gate = build_pilot_gate(
        wt_controls=wt_control_rows,
        paired_rows=paired_rows,
        summaries=summaries,
        expected_candidate_count=len(candidates),
        expected_replicates=args.replicates,
    )
    gate.update(
        {
            "generated_at": generated_at,
            "selected_protocol": protocol,
            "score_function": runtime.SCORE_FUNCTION,
            "candidate_id_file": str(candidate_id_file),
        }
    )

    with tempfile.TemporaryDirectory(prefix=".affinity-scoring-", dir=PROJECT_ROOT) as temp:
        staging = Path(temp)
        staged = {key: staging / name for key, name in OUTPUT_NAMES.items()}
        staged_summary = staging / "run_summary.json"
        _write_csv(staged["wt_controls"], wt_control_rows, WT_CONTROL_FIELDS)
        _write_csv(staged["paired"], paired_rows, PAIRED_FIELDS)
        _write_csv(staged["summary"], summaries, SUMMARY_FIELDS)
        _write_json(staged["gate"], gate)
        _render_svg(summaries, staged["figure"])
        _write_json(
            staged_summary,
            {
                "schema_version": 2,
                "status": gate["status"],
                "generated_at": generated_at,
                "elapsed_seconds": round(time.perf_counter() - started, 6),
                "python": platform.python_version(),
                "pyrosetta_version": version,
                "command_argv": [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
                "candidate_count": len(candidates),
                "wt_control_count": len(wt_control_rows),
                "mutant_evaluation_count": len(paired_rows),
                "replicate_count": args.replicates,
                "candidate_filtering_applied": False,
                "full_scan_contract": "score_all_declared_candidates_then_filter_once",
                "outputs": {key: str(path) for key, path in final_paths.items()},
            },
        )
        replace_staged_files(
            {
                **{staged[key]: final_paths[key] for key in OUTPUT_NAMES},
                staged_summary: run_summary,
            },
            project_root=PROJECT_ROOT,
            protected_source_paths=validated.source_paths,
        )
    return 0 if gate["status"] == "pass" else 2


def _render_svg(rows: list[dict[str, object]], path: Path) -> None:
    width, height = 900, max(360, 90 + len(rows) * 24)
    values = [float(row["delta_dG_separated_median"]) for row in rows] or [0.0]
    bound = max(1.0, max(abs(value) for value in values))
    center, scale = 520, 300 / (2 * bound)
    items = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="24" y="34" font-family="Arial" font-size="20">Paired PyRosetta affinity pilot V2</text>',
        '<text x="24" y="54" font-family="Arial" font-size="12">delta dG_separated = mutant - paired WT (REU; negative is favorable)</text>',
        '<text x="24" y="68" font-family="Arial" font-size="12">Unfiltered scan-stage results; no candidate selection is applied.</text>',
        f'<line x1="{center}" y1="76" x2="{center}" y2="{height - 24}" stroke="#555" stroke-width="1"/>',
    ]
    for index, row in enumerate(rows):
        y = 88 + index * 24
        value = float(row["delta_dG_separated_median"])
        x = center + value * scale
        color = "#2b8cbe" if row["status"] == "pass" else "#d7301f"
        items.extend(
            [
                f'<text x="24" y="{y + 4}" font-family="Arial" font-size="11">{row["candidate_id"]}</text>',
                f'<line x1="{center}" y1="{y}" x2="{x:.2f}" y2="{y}" stroke="{color}" stroke-width="5"/>',
                f'<circle cx="{x:.2f}" cy="{y}" r="4" fill="{color}"/>',
                f'<text x="830" y="{y + 4}" font-family="Arial" font-size="11" text-anchor="end">{value:.3f}</text>',
            ]
        )
    items.append("</svg>")
    path.write_text("\n".join(items) + "\n", encoding="utf-8", newline="\n")


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
