#!/usr/bin/env python3
"""Score the 43 TNP-applicable VHHs and preserve all 47 V2 status rows."""

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
from antibody_optimization.tnp_yield import (  # noqa: E402
    ELIGIBLE_COUNT,
    TNP_NOT_APPLICABLE,
    failed_tnp_result,
    normalize_tnp_result,
    not_applicable_tnp_result,
    verify_immune_builder_refine_patch,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--tnp-source", type=Path, required=True)
    parser.add_argument("--tnp-executable", type=Path, required=True)
    parser.add_argument("--immune-builder-refine", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--generated-at")
    args = parser.parse_args()
    started = time.perf_counter()
    generated = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    plan = args.plan_dir.resolve(strict=True)
    source = args.tnp_source.resolve(strict=True)
    executable = args.tnp_executable.resolve(strict=True)
    refine_file = args.immune_builder_refine.resolve(strict=True)
    samples = _csv(plan / "tnp_validation_samples.csv")
    contract = _json(plan / "tnp_yield_validation_contract.json")
    if len(samples) != 47 or contract.get("release") != "ready_for_remote_single_process_tnp_v2_scoring":
        raise ValueError("TNP plan is incomplete or not released")
    verify_immune_builder_refine_patch(refine_file.read_text(encoding="utf-8"))
    output_dir = args.output_dir.absolute()
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite TNP score directory: {output_dir}")
    output_dir.mkdir(parents=True)
    raw_dir = output_dir / "raw_tnp"
    raw_dir.mkdir()
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(source)
    rows = []
    eligible_index = 0
    for index, sample in enumerate(samples, 1):
        uid = str(sample["sample_uid"])
        if str(sample["tnp_applicability"]) == "not_applicable":
            rows.append(not_applicable_tnp_result(sample))
            print(f"[skip] {uid}: not_applicable ({sample['tnp_inapplicability_reason']})", flush=True)
            continue
        eligible_index += 1
        print(f"[{eligible_index}/{ELIGIBLE_COUNT}] TNP scoring {uid}", flush=True)
        sample_dir = raw_dir / f"sample_{index:02d}"
        sample_started = time.perf_counter()
        command = [str(executable), "--seq", str(sample["sequence_raw"]), "--name", uid, "--output", str(sample_dir), "--hscale", "0", "--ncores", "1"]
        try:
            completed = subprocess.run(command, cwd=source, env=environment, text=True, capture_output=True, check=False)
            (sample_dir.parent / f"sample_{index:02d}.stdout.log").write_text(completed.stdout, encoding="utf-8", newline="\n")
            (sample_dir.parent / f"sample_{index:02d}.stderr.log").write_text(completed.stderr, encoding="utf-8", newline="\n")
            if completed.returncode:
                raise RuntimeError(_tnp_failure_reason(completed, sample_dir))
            result_path = _single(sample_dir.glob("TNP_Results*.json"), f"result JSON for {uid}")
            details_path = _single(sample_dir.glob("Final_Models/*_Model_Details.json"), f"model details for {uid}")
            details = _json(details_path)
            modelled = str(details["sequences"]["H"])
            rows.append(normalize_tnp_result(sample, _json(result_path), modelled_sequence=modelled, elapsed_seconds=time.perf_counter() - sample_started))
        except Exception as error:  # preserve one explicit row per planned real sample
            rows.append(failed_tnp_result(sample, f"{type(error).__name__}: {error}", time.perf_counter() - sample_started))
        print(f"[{eligible_index}/{ELIGIBLE_COUNT}] {uid}: {rows[-1]['scoring_status']}", flush=True)
    _write_csv(output_dir / "tnp_sample_scores.csv", rows)
    passed = sum(row["scoring_status"] == "pass" for row in rows)
    not_applicable = sum(row["scoring_status"] == "not_applicable" for row in rows)
    eligible_failed = ELIGIBLE_COUNT - passed
    status = "pass" if passed == ELIGIBLE_COUNT and not_applicable == len(TNP_NOT_APPLICABLE) else "failed"
    _write_json(output_dir / "tnp_model_run.json", {"schema_version": 2, "status": status, "generated_at": generated, "python": platform.python_version(), "elapsed_seconds": round(time.perf_counter() - started, 6), "execution": "single_process_sequential_eligible_only", "planned_count": 47, "eligible_count": ELIGIBLE_COUNT, "pass_count": passed, "eligible_failure_count": eligible_failed, "not_applicable_count": not_applicable, "tnp_executable": str(executable), "tnp_source": str(source), "immune_builder_refine": str(refine_file), "immune_builder_patch_verified": True, "hscale": 0, "ncores": 1})
    if status != "pass":
        raise RuntimeError(f"TNP V2 coverage gate failed: {passed}/{ELIGIBLE_COUNT} eligible passed")
    return 0


def _single(paths, label):
    found = list(paths)
    if len(found) != 1: raise ValueError(f"Expected exactly one {label}; found {len(found)}")
    return found[0]
def _tnp_failure_reason(completed, sample_dir):
    internal = list(sample_dir.glob("*_TNP.log"))
    internal_text = internal[0].read_text(encoding="utf-8", errors="replace") if len(internal) == 1 else ""
    if "failed to number the sequence with ANARCI" in internal_text:
        return "TNP ANARCI numbering failed"
    if "unable to generate a model" in internal_text:
        stderr = completed.stderr.strip().splitlines()
        detail = next((line.strip() for line in stderr if "TypeError:" in line or "Error:" in line), "no detailed stderr")
        return f"NanoBodyBuilder2 failed: {detail}"
    return f"TNP exit code {completed.returncode}"
def _csv(path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle: return list(csv.DictReader(handle))
def _json(path): return json.loads(path.read_text(encoding="utf-8-sig"))
def _write_csv(path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n"); writer.writeheader(); writer.writerows(rows)
def _write_json(path, value): path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__": raise SystemExit(main())
