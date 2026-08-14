#!/usr/bin/env python3
"""Run the fixed external NetSolP entry point and normalize its 47 scores."""

from __future__ import annotations

import argparse
import csv
import json
import platform
import subprocess
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from antibody_optimization.file_transaction import replace_staged_files, validate_file_paths  # noqa: E402
from antibody_optimization.netsolp_yield import normalize_netsolp_scores  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--netsolp-workdir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--num-threads", type=int, default=12)
    parser.add_argument("--generated-at")
    args = parser.parse_args()
    if args.num_threads < 1:
        raise ValueError("--num-threads must be positive")
    started = time.perf_counter()
    generated = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    plan = args.plan_dir.resolve(strict=True)
    workdir = args.netsolp_workdir.resolve(strict=True)
    sources = [
        plan / "netsolp_validation_samples.csv",
        plan / "netsolp_validation_sequences.fasta",
        plan / "netsolp_yield_validation_contract.json",
        workdir / "predict.py",
    ]
    samples = _csv(sources[0])
    contract = _json(sources[2])
    if len(samples) != 47 or contract.get("release") != "ready_for_remote_netsolp_scoring":
        raise ValueError("NetSolP plan is incomplete or not released")
    if contract.get("model_type") != "Distilled" or contract.get("prediction_type") != "SU":
        raise ValueError("NetSolP plan does not match the fixed Distilled SU contract")
    output_dir = args.output_dir.absolute()
    targets = [
        output_dir / "netsolp_raw_predictions.csv",
        output_dir / "netsolp_sample_scores.csv",
        output_dir / "netsolp_model_run.json",
    ]
    valid = validate_file_paths(project_root=ROOT, source_paths=sources, target_paths=targets)
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite NetSolP score directory: {output_dir}")
    output_dir.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix=".netsolp-score-", dir=output_dir.parent) as temp:
        stage = Path(temp)
        raw = stage / targets[0].name
        command = [
            sys.executable,
            str(sources[3]),
            "--FASTA_PATH",
            str(sources[1]),
            "--OUTPUT_PATH",
            str(raw),
            "--MODEL_TYPE",
            "Distilled",
            "--PREDICTION_TYPE",
            "SU",
            "--NUM_THREADS",
            str(args.num_threads),
        ]
        completed = subprocess.run(command, cwd=workdir, check=True, text=True, capture_output=True)
        raw_rows = _csv(raw)
        score_rows = normalize_netsolp_scores(samples, raw_rows)
        scores = stage / targets[1].name
        model_run = stage / targets[2].name
        _write_csv(scores, score_rows)
        _write_json(
            model_run,
            {
                "schema_version": 1,
                "status": "pass",
                "generated_at": generated,
                "python": platform.python_version(),
                "elapsed_seconds": round(time.perf_counter() - started, 6),
                "command": command,
                "working_directory": str(workdir),
                "model_type": "Distilled",
                "prediction_type": "SU",
                "sample_count": 47,
                "stdout": completed.stdout.strip().splitlines(),
                "stderr": completed.stderr.strip().splitlines(),
            },
        )
        output_dir.mkdir(parents=False, exist_ok=False)
        replace_staged_files(
            {raw: targets[0], scores: targets[1], model_run: targets[2]},
            project_root=ROOT,
            protected_source_paths=valid.source_paths,
        )
    return 0


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
