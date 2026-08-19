#!/usr/bin/env python3
"""Freeze AntiFold's valid role relative to the 47-sample BL21-yield set."""

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

from antibody_optimization.antifold_validation import build_antifold_yield_applicability_contract  # noqa: E402
from antibody_optimization.file_transaction import replace_staged_files, validate_file_paths  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--yield-samples", type=Path, required=True)
    parser.add_argument("--antifold-plan", type=Path, required=True)
    parser.add_argument("--structure-views", type=Path, required=True)
    parser.add_argument("--antifold-gate", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    args = parser.parse_args()
    started = time.perf_counter()
    generated = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    sources = [
        args.yield_samples.resolve(strict=True),
        args.antifold_plan.resolve(strict=True),
        args.structure_views.resolve(strict=True),
        args.antifold_gate.resolve(strict=True),
    ]
    contract = build_antifold_yield_applicability_contract(
        _csv(sources[0]), _json(sources[1]), _csv(sources[2]), _json(sources[3])
    )
    contract["generated_at"] = generated
    output = args.output_dir.absolute() / "antifold_yield_applicability_contract.json"
    summary = args.run_summary.absolute()
    valid = validate_file_paths(
        project_root=ROOT, source_paths=sources, target_paths=[output, summary]
    )
    if any(path.exists() for path in valid.target_paths):
        raise FileExistsError("Refusing to overwrite AntiFold applicability outputs")
    for path in valid.target_paths:
        path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".antifold-yield-role-", dir=ROOT) as temp:
        stage = Path(temp)
        staged_contract = stage / output.name
        staged_summary = stage / summary.name
        _write_json(staged_contract, contract)
        _write_json(
            staged_summary,
            {
                "schema_version": 1,
                "status": "pass",
                "generated_at": generated,
                "python": platform.python_version(),
                "elapsed_seconds": round(time.perf_counter() - started, 6),
                "classification_status": contract["classification_status"],
                "output": str(output),
            },
        )
        replace_staged_files(
            {staged_contract: output, staged_summary: summary},
            project_root=ROOT,
            protected_source_paths=valid.source_paths,
        )
    return 0


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _write_json(path: Path, value: object) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")


if __name__ == "__main__":
    raise SystemExit(main())
