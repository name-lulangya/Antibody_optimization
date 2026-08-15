#!/usr/bin/env python3
"""Score the fixed 96-sequence unified Nb252 TNP review plan."""

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
    failed_tnp_result,
    normalize_tnp_result,
    verify_immune_builder_refine_patch,
)
from antibody_optimization.unified_tnp_review import SCORE_COUNT  # noqa: E402


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
    samples = _csv(plan / "unified_tnp_samples.csv")
    contract = _json(plan / "unified_tnp_review_contract.json")
    if (
        len(samples) != SCORE_COUNT
        or len({row["sample_uid"] for row in samples}) != SCORE_COUNT
        or contract.get("release") != "ready_for_remote_single_process_unified_tnp_review"
        or int(contract.get("planned_count", 0)) != SCORE_COUNT
    ):
        raise ValueError("Unified TNP review plan is incomplete or not released")
    verify_immune_builder_refine_patch(refine_file.read_text(encoding="utf-8"))

    output_dir = args.output_dir.absolute()
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite TNP score directory: {output_dir}")
    output_dir.mkdir(parents=True)
    raw_dir = output_dir / "raw_tnp"
    raw_dir.mkdir()
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(source)
    rows: list[dict[str, object]] = []
    for index, sample in enumerate(samples, 1):
        uid = str(sample["sample_uid"])
        print(f"[{index}/{SCORE_COUNT}] TNP scoring {uid}", flush=True)
        sample_dir = raw_dir / f"sample_{index:03d}"
        sample_started = time.perf_counter()
        command = [
            str(executable), "--seq", str(sample["sequence_raw"]), "--name", uid,
            "--output", str(sample_dir), "--hscale", "0", "--ncores", "1",
        ]
        try:
            completed = subprocess.run(
                command, cwd=source, env=environment, text=True,
                capture_output=True, check=False,
            )
            (raw_dir / f"sample_{index:03d}.stdout.log").write_text(
                completed.stdout, encoding="utf-8", newline="\n"
            )
            (raw_dir / f"sample_{index:03d}.stderr.log").write_text(
                completed.stderr, encoding="utf-8", newline="\n"
            )
            if completed.returncode:
                raise RuntimeError(_tnp_failure_reason(completed, sample_dir))
            result_path = _single(sample_dir.glob("TNP_Results*.json"), f"result JSON for {uid}")
            details_path = _single(
                sample_dir.glob("Final_Models/*_Model_Details.json"),
                f"model details for {uid}",
            )
            modelled = str(_json(details_path)["sequences"]["H"])
            row = normalize_tnp_result(
                sample,
                _json(result_path),
                modelled_sequence=modelled,
                elapsed_seconds=time.perf_counter() - sample_started,
            )
            if row["trimmed_n_terminal"] or row["trimmed_c_terminal"] != "GS":
                raise ValueError("TNP modeled an unexpected sequence domain")
            if int(row["modelled_length_aa"]) != 126:
                raise ValueError("TNP modeled length is not the contracted 126 aa")
            rows.append(row)
        except Exception as error:  # retain one explicit record per real input
            rows.append(
                failed_tnp_result(
                    sample,
                    f"{type(error).__name__}: {error}",
                    time.perf_counter() - sample_started,
                )
            )
        print(f"[{index}/{SCORE_COUNT}] {uid}: {rows[-1]['scoring_status']}", flush=True)

    _write_csv(output_dir / "tnp_sample_scores.csv", rows)
    passed = sum(row["scoring_status"] == "pass" for row in rows)
    status = "pass" if passed == SCORE_COUNT else "failed"
    _write_json(output_dir / "tnp_model_run.json", {
        "schema_version": 1,
        "status": status,
        "generated_at": generated,
        "python": platform.python_version(),
        "elapsed_seconds": round(time.perf_counter() - started, 6),
        "execution": "single_process_sequential",
        "planned_count": SCORE_COUNT,
        "pass_count": passed,
        "failure_count": SCORE_COUNT - passed,
        "tnp_executable": str(executable),
        "tnp_source": str(source),
        "immune_builder_refine": str(refine_file),
        "immune_builder_patch_verified": True,
        "hscale": 0,
        "ncores": 1,
    })
    if status != "pass":
        raise RuntimeError(f"Unified TNP coverage gate failed: {passed}/{SCORE_COUNT} passed")
    return 0


def _single(paths, label: str) -> Path:
    found = list(paths)
    if len(found) != 1:
        raise ValueError(f"Expected exactly one {label}; found {len(found)}")
    return found[0]


def _tnp_failure_reason(completed: subprocess.CompletedProcess[str], sample_dir: Path) -> str:
    internal = list(sample_dir.glob("*_TNP.log"))
    internal_text = internal[0].read_text(encoding="utf-8", errors="replace") if len(internal) == 1 else ""
    if "failed to number the sequence with ANARCI" in internal_text:
        return "TNP ANARCI numbering failed"
    if "unable to generate a model" in internal_text:
        stderr = completed.stderr.strip().splitlines()
        detail = next(
            (line.strip() for line in stderr if "TypeError:" in line or "Error:" in line),
            "no detailed stderr",
        )
        return f"NanoBodyBuilder2 failed: {detail}"
    return f"TNP exit code {completed.returncode}"


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
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8", newline="\n",
    )


if __name__ == "__main__":
    raise SystemExit(main())
