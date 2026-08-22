#!/usr/bin/env python3
"""Score WT plus 162 active double mutants with NetSolP or NanoMelt."""

from __future__ import annotations

import argparse
import csv
import json
import platform
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

EXPECTED = 163


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tool", choices=("netsolp", "nanomelt"), required=True)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--tool-root", type=Path)
    parser.add_argument("--executable", type=Path)
    parser.add_argument("--generated-at")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.perf_counter()
    plan = args.plan_dir.resolve(strict=True)
    output = args.output_dir.absolute()
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite: {output}")
    samples = _csv(plan / "expression_double_mutant_score_samples.csv")
    gate = _json(plan / "expression_double_mutant_plan_gate.json")
    if (
        len(samples) != EXPECTED
        or len({row["sample_uid"] for row in samples}) != EXPECTED
        or gate.get("release") != "ready_for_netsolp_nanomelt_double_scoring"
    ):
        raise ValueError("Double-mutant scoring plan is incomplete or not released")
    output.mkdir(parents=True)
    if args.tool == "netsolp":
        rows, command = _score_netsolp(args, plan, samples, output)
    else:
        rows, command = _score_nanomelt(args, plan, samples, output)
    passed = sum(str(row["scoring_status"]) == "pass" for row in rows)
    status = "pass" if passed == EXPECTED else "failed"
    _write_csv(output / f"{args.tool}_sample_scores.csv", rows)
    _write_json(
        output / f"{args.tool}_model_run.json",
        {
            "schema_version": 1,
            "status": status,
            "generated_at": args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds"),
            "tool": args.tool,
            "python": platform.python_version(),
            "sample_count": EXPECTED,
            "pass_count": passed,
            "failure_count": EXPECTED - passed,
            "elapsed_seconds": round(time.perf_counter() - started, 6),
            "command": command,
            "candidate_filtering_applied": False,
        },
    )
    if status != "pass":
        raise RuntimeError(f"{args.tool} coverage gate failed: {passed}/{EXPECTED}")
    return 0


def _score_netsolp(args, plan, samples, output):
    from antibody_optimization.netsolp_yield import normalize_netsolp_scores

    if not args.tool_root:
        raise ValueError("NetSolP requires --tool-root")
    root = args.tool_root.resolve(strict=True)
    raw = output / "netsolp_raw_predictions.csv"
    command = [
        sys.executable,
        str(root / "predict.py"),
        "--FASTA_PATH",
        str(plan / "expression_double_mutant_sequences.fasta"),
        "--OUTPUT_PATH",
        str(raw),
        "--MODEL_TYPE",
        "Distilled",
        "--PREDICTION_TYPE",
        "SU",
        "--NUM_THREADS",
        "12",
    ]
    subprocess.run(command, cwd=root, check=True)
    normalized_samples = [
        {"sample_uid": row["sample_uid"], "sequence_raw": row["sequence_raw"]}
        for row in samples
    ]
    return normalize_netsolp_scores(normalized_samples, _csv(raw), expected_count=EXPECTED), command


def _score_nanomelt(args, plan, samples, output):
    from antibody_optimization.nanomelt_yield import normalize_nanomelt_scores

    if not args.executable:
        raise ValueError("NanoMelt requires --executable")
    import torch

    if not torch.cuda.is_available():
        raise RuntimeError("NanoMelt requires a visible CUDA GPU")
    raw = output / "nanomelt_raw_predictions.csv"
    command = [
        str(args.executable.resolve(strict=True)),
        "predict",
        "-i",
        str(plan / "expression_double_mutant_sequences.fasta"),
        "-o",
        str(raw),
        "-align",
        "-ncpu",
        "1",
        "-v",
    ]
    subprocess.run(command, check=True)
    normalized_samples = [
        {"sample_uid": row["sample_uid"], "sequence_raw": row["sequence_raw"]}
        for row in samples
    ]
    return normalize_nanomelt_scores(
        normalized_samples,
        _csv(raw),
        expected_pass_count=EXPECTED,
        expected_plan_count=EXPECTED,
    ), command


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    raise SystemExit(main())
