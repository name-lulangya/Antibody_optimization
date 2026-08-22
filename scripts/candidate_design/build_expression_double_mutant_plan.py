#!/usr/bin/env python3
"""Enumerate the complete 162-member double-mutant scoring plan."""

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

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from antibody_optimization.expression_double_mutants import (  # noqa: E402
    build_double_mutant_space,
    build_score_samples,
)
from antibody_optimization.file_transaction import replace_staged_files, validate_file_paths  # noqa: E402


NAMES = {
    "candidates": "expression_double_mutant_candidates.csv",
    "invalid": "expression_double_mutant_invalid_same_position_pairs.csv",
    "samples": "expression_double_mutant_score_samples.csv",
    "fasta": "expression_double_mutant_sequences.fasta",
    "word_changes": "expression_double_mutant_stable_word_changes.csv",
    "contract": "expression_double_mutant_plan_contract.json",
    "gate": "expression_double_mutant_plan_gate.json",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent19-dir", type=Path, required=True)
    parser.add_argument("--stable-word-library", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    parent_dir = args.parent19_dir.resolve(strict=True)
    parent_csv = parent_dir / "expression_single_mutant_parent19.csv"
    parent_gate_path = parent_dir / "expression_single_mutant_parent19_gate.json"
    word_path = args.stable_word_library.resolve(strict=True)
    output_dir = args.output_dir.absolute()
    output_dir.mkdir(parents=True, exist_ok=True)
    args.run_summary.parent.mkdir(parents=True, exist_ok=True)
    targets = [*(output_dir / name for name in NAMES.values()), args.run_summary.absolute()]
    validated = validate_file_paths(
        project_root=ROOT,
        source_paths=[parent_csv, parent_gate_path, word_path],
        target_paths=targets,
    )
    if existing := [path for path in validated.target_paths if path.exists()]:
        raise FileExistsError("Refusing to overwrite:\n" + "\n".join(map(str, existing)))
    words = [row["stable_word"] for row in _csv(word_path)]
    result = build_double_mutant_space(_csv(parent_csv), _json(parent_gate_path), words)
    samples = build_score_samples(str(result["parent_sequence"]), result["candidates"])
    facts = result["facts"]
    contract = {
        "schema_version": 1,
        "generated_at": generated_at,
        "optimization_target": "BL21 expression yield",
        "parent_release": "19 explicitly approved single mutants",
        "enumeration_rule": "all unordered pairs at distinct reported-sequence positions",
        "same_position_alternatives": "mutually_exclusive",
        "candidate_count": 162,
        "score_sample_count_including_wt": 163,
        "active_remote_predictors": ["NetSolP Distilled SU", "NanoMelt predicted apparent Tm"],
        "antifold_policy": (
            "reuse each constituent position's frozen structure-conditioned evidence; "
            "add deltas only when both positions use the same structural view; never interpret "
            "the sum as double-mutant epistasis"
        ),
        "sequence_risk_policy": "recompute on each complete double sequence",
        "stable_word_policy": "recompute exact overlapping degenerate-word occurrences on each complete double sequence",
        "candidate_selection_performed": False,
        "final_11_double_mutants_selected": False,
        "expected_remote_runtime": "under_5_hours",
        "slurm_route": "one_batch_job_one_gpu_12_cpus_sequential_netsolp_then_nanomelt",
    }
    gate = {
        "schema_version": 1,
        "generated_at": generated_at,
        "status": "pass" if not facts["hard_sequence_risk_count"] else "blocked",
        "release": (
            "ready_for_netsolp_nanomelt_double_scoring"
            if not facts["hard_sequence_risk_count"]
            else "blocked_by_combined_sequence_risk"
        ),
        **facts,
        "interpretation": "Complete unfiltered 162-double space; no final 11 selection has been performed.",
    }

    finals = dict(zip(NAMES, validated.target_paths[:-1], strict=True))
    run_summary = validated.target_paths[-1]
    with tempfile.TemporaryDirectory(prefix=".expression-double-plan-", dir=ROOT) as temp:
        stage = Path(temp)
        staged = {key: stage / name for key, name in NAMES.items()}
        staged_run = stage / "run_summary.json"
        _write_csv(staged["candidates"], result["candidates"])
        _write_csv(staged["invalid"], result["invalid_pairs"])
        _write_csv(staged["samples"], samples)
        _write_fasta(staged["fasta"], samples)
        _write_csv(staged["word_changes"], result["stable_word_changes"])
        _write_json(staged["contract"], contract)
        _write_json(staged["gate"], gate)
        _write_json(
            staged_run,
            {
                "schema_version": 1,
                "status": gate["status"],
                "generated_at": generated_at,
                "elapsed_seconds": round(time.perf_counter() - started, 6),
                "python": platform.python_version(),
                "command_argv": [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
                **facts,
            },
        )
        replace_staged_files(
            {**{staged[key]: finals[key] for key in NAMES}, staged_run: run_summary},
            project_root=ROOT,
            protected_source_paths=validated.source_paths,
        )
    if gate["status"] != "pass":
        raise RuntimeError("Double-mutant plan is blocked by a combined sequence risk")
    return 0


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8-sig", newline="")
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_fasta(path: Path, rows: list[dict[str, object]]) -> None:
    lines: list[str] = []
    for row in rows:
        lines.extend([f">{row['sample_uid']}", str(row["sequence_raw"])])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def _write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    raise SystemExit(main())
