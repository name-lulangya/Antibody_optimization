#!/usr/bin/env python3
"""Run one fixed NetSolP batch for a v2 reuse-validation or completion subset."""

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
    parser.add_argument("--sample-table", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--netsolp-workdir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--num-threads", type=int, default=12)
    parser.add_argument("--generated-at")
    args = parser.parse_args()
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    samples_path = args.sample_table.resolve(strict=True)
    fasta_path = args.fasta.resolve(strict=True)
    contract_path = args.contract.resolve(strict=True)
    workdir = args.netsolp_workdir.resolve(strict=True)
    predict = (workdir / "predict.py").resolve(strict=True)
    samples = _csv(samples_path)
    contract = _json(contract_path)
    spec = contract["netsolp"]
    if not samples or spec["model_type"] != "Distilled" or spec["prediction_type"] != "SU":
        raise ValueError("Invalid NetSolP subset contract")
    if len({row["score_id"] for row in samples}) != len(samples):
        raise ValueError("NetSolP score IDs are not unique")
    output_dir = args.output_dir.absolute()
    targets = [
        output_dir / "netsolp_raw_predictions.csv",
        output_dir / "netsolp_sample_scores.csv",
        output_dir / "netsolp_model_run.json",
    ]
    validated = validate_file_paths(
        project_root=ROOT,
        source_paths=[samples_path, fasta_path, contract_path, predict],
        target_paths=targets,
    )
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite NetSolP subset scores: {output_dir}")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    normalized_samples = [{"sample_uid": row["score_id"], "sequence_raw": row["sequence_raw"]} for row in samples]
    with tempfile.TemporaryDirectory(prefix=".expression-property-netsolp-", dir=output_dir.parent) as temp_name:
        staging = Path(temp_name)
        raw_path = staging / targets[0].name
        command = [
            sys.executable, str(predict), "--FASTA_PATH", str(fasta_path), "--OUTPUT_PATH", str(raw_path),
            "--MODEL_TYPE", "Distilled", "--PREDICTION_TYPE", "SU", "--NUM_THREADS", str(args.num_threads),
        ]
        completed = subprocess.run(command, cwd=workdir, check=True, text=True, capture_output=True)
        scores = normalize_netsolp_scores(normalized_samples, _csv(raw_path), expected_count=len(samples))
        score_path = staging / targets[1].name
        run_path = staging / targets[2].name
        _write_csv(score_path, scores)
        _write_json(run_path, {
            "schema_version": 1, "status": "pass", "generated_at": generated_at,
            "python": platform.python_version(), "elapsed_seconds": round(time.perf_counter() - started, 6),
            "command": command, "sample_count": len(scores),
            "pass_count": sum(row["scoring_status"] == "pass" for row in scores),
            "model_type": "Distilled", "prediction_type": "SU",
            "stdout": completed.stdout.strip().splitlines(), "stderr": completed.stderr.strip().splitlines(),
        })
        output_dir.mkdir(parents=False, exist_ok=False)
        replace_staged_files(
            {raw_path: targets[0], score_path: targets[1], run_path: targets[2]},
            project_root=ROOT, protected_source_paths=validated.source_paths,
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
        writer.writeheader(); writer.writerows(rows)


def _write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    raise SystemExit(main())
