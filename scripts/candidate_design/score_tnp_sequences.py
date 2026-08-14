#!/usr/bin/env python3
"""Score 47 planned VHHs sequentially with the fixed external TNP entry point."""

from __future__ import annotations

import argparse
import csv
import json
import os
import platform
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from antibody_optimization.tnp_yield import failed_tnp_result, normalize_tnp_result  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--tnp-source", type=Path, required=True)
    parser.add_argument("--tnp-executable", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--generated-at")
    args = parser.parse_args()
    started = time.perf_counter()
    generated = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    plan = args.plan_dir.resolve(strict=True)
    source = args.tnp_source.resolve(strict=True)
    executable = args.tnp_executable.resolve(strict=True)
    samples = _csv(plan / "tnp_validation_samples.csv")
    contract = _json(plan / "tnp_yield_validation_contract.json")
    if len(samples) != 47 or contract.get("release") != "ready_for_remote_single_process_tnp_scoring":
        raise ValueError("TNP plan is incomplete or not released")
    output_dir = args.output_dir.absolute()
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite TNP score directory: {output_dir}")
    output_dir.mkdir(parents=True)
    raw_dir = output_dir / "raw_tnp"
    raw_dir.mkdir()
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(source)
    rows = []
    for index, sample in enumerate(samples, 1):
        uid = str(sample["sample_uid"])
        print(f"[{index}/47] TNP scoring {uid}", flush=True)
        sample_dir = raw_dir / f"sample_{index:02d}"
        sample_started = time.perf_counter()
        command = [str(executable), "--seq", str(sample["sequence_raw"]), "--name", uid, "--output", str(sample_dir), "--hscale", "0", "--ncores", "1"]
        try:
            completed = subprocess.run(command, cwd=source, env=environment, text=True, capture_output=True, check=False)
            (sample_dir.parent / f"sample_{index:02d}.stdout.log").write_text(completed.stdout, encoding="utf-8", newline="\n")
            (sample_dir.parent / f"sample_{index:02d}.stderr.log").write_text(completed.stderr, encoding="utf-8", newline="\n")
            if completed.returncode:
                raise RuntimeError(f"TNP exit code {completed.returncode}")
            result_path = _single(sample_dir.glob("TNP_Results*.json"), f"result JSON for {uid}")
            details_path = _single(sample_dir.glob("Final_Models/*_Model_Details.json"), f"model details for {uid}")
            details = _json(details_path)
            modelled = str(details["sequences"]["H"])
            rows.append(normalize_tnp_result(sample, _json(result_path), modelled_sequence=modelled, elapsed_seconds=time.perf_counter() - sample_started))
        except Exception as error:  # preserve one explicit row per planned real sample
            rows.append(failed_tnp_result(sample, f"{type(error).__name__}: {error}", time.perf_counter() - sample_started))
        print(f"[{index}/47] {uid}: {rows[-1]['scoring_status']}", flush=True)
    _write_csv(output_dir / "tnp_sample_scores.csv", rows)
    passed = sum(row["scoring_status"] == "pass" for row in rows)
    _write_json(output_dir / "tnp_model_run.json", {"schema_version": 1, "status": "pass" if passed >= 46 else "failed", "generated_at": generated, "python": platform.python_version(), "elapsed_seconds": round(time.perf_counter() - started, 6), "execution": "single_process_sequential", "sample_count": 47, "pass_count": passed, "failure_count": 47 - passed, "tnp_executable": str(executable), "tnp_source": str(source), "hscale": 0, "ncores": 1})
    if passed < 46:
        raise RuntimeError(f"TNP coverage gate failed after scoring: {passed}/47 passed")
    return 0


def _single(paths, label):
    found = list(paths)
    if len(found) != 1: raise ValueError(f"Expected exactly one {label}; found {len(found)}")
    return found[0]
def _csv(path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle: return list(csv.DictReader(handle))
def _json(path): return json.loads(path.read_text(encoding="utf-8-sig"))
def _write_csv(path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n"); writer.writeheader(); writer.writerows(rows)
def _write_json(path, value): path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__": raise SystemExit(main())
