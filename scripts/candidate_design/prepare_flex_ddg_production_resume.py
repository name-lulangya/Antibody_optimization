#!/usr/bin/env python3
"""Inspect production outputs and write pending-task array index files."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import datetime
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from antibody_optimization.flex_ddg_production import (  # noqa: E402
    DEFAULT_ARRAY_CHUNK_SIZE,
    DEFAULT_ARRAY_CONCURRENCY,
    assess_task_outputs,
    chunk_task_indices,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--task-root", type=Path, required=True)
    parser.add_argument("--submission-dir", type=Path, required=True)
    parser.add_argument("--concurrency", type=int, default=DEFAULT_ARRAY_CONCURRENCY)
    parser.add_argument("--chunk-size", type=int, default=DEFAULT_ARRAY_CHUNK_SIZE)
    parser.add_argument("--generated-at")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.concurrency <= 0:
        raise ValueError("Concurrency must be positive")
    plan_dir = _project_directory(args.plan_dir)
    task_root = _project_path(args.task_root)
    submission_dir = _project_path(args.submission_dir)
    if submission_dir.exists():
        raise FileExistsError(f"Submission directory already exists: {submission_dir}")
    manifest = _load_csv(plan_dir / "flex_ddg_production_manifest.csv")
    assessment = assess_task_outputs(manifest_rows=manifest, task_root=task_root)
    chunks = chunk_task_indices(assessment["pending_task_indices"], args.chunk_size)
    submission_dir.mkdir(parents=True)
    for chunk_number, indices in enumerate(chunks):
        (submission_dir / f"chunk_{chunk_number:03d}.txt").write_text(
            "".join(f"{index}\n" for index in indices), encoding="ascii", newline="\n"
        )
    state = {
        "schema_version": 1,
        "generated_at": args.generated_at
        or datetime.now().astimezone().isoformat(timespec="seconds"),
        "status": "blocked" if assessment["invalid_count"] else "pass",
        "concurrency": args.concurrency,
        "chunk_size": args.chunk_size,
        "chunk_count": len(chunks),
        **assessment,
    }
    (submission_dir / "resume_state.json").write_text(
        json.dumps(state, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(json.dumps({key: state[key] for key in ("status", "completed_count", "pending_count", "invalid_count", "chunk_count")}, sort_keys=True))
    return 2 if assessment["invalid_count"] else 0


def _project_directory(path: Path) -> Path:
    resolved = path.expanduser().resolve(strict=True)
    resolved.relative_to(PROJECT_ROOT.resolve(strict=True))
    if resolved.is_symlink() or not resolved.is_dir():
        raise ValueError(f"Expected regular project directory: {resolved}")
    return resolved


def _project_path(path: Path) -> Path:
    resolved = path.expanduser().absolute()
    resolved.relative_to(PROJECT_ROOT.absolute())
    return resolved


def _load_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


if __name__ == "__main__":
    raise SystemExit(main())
