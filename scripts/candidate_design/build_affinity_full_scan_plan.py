#!/usr/bin/env python3
"""Build the deterministic 12-shard plan for the unfiltered affinity scan."""

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

from antibody_optimization.affinity_full_scan import (  # noqa: E402
    CANDIDATES_PER_SHARD,
    EXPECTED_CANDIDATES,
    EXPECTED_REPLICATES,
    SHARD_COUNT,
    SHARD_MANIFEST_FIELDS,
    build_full_scan_shards,
)
from antibody_optimization.file_transaction import (  # noqa: E402
    replace_staged_files,
    validate_file_paths,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-dir", type=Path, required=True)
    parser.add_argument("--pilot-v2-review", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(
        timespec="seconds"
    )
    candidate_dir = _project_directory(args.candidate_dir)
    review_path = args.pilot_v2_review.expanduser().resolve(strict=True)
    review_path.relative_to(PROJECT_ROOT.resolve(strict=True))
    candidate_path = candidate_dir / "affinity_single_mutants.csv"
    gate_path = candidate_dir / "affinity_candidate_gate.json"
    candidates = _load_csv(candidate_path)
    candidate_gate = _load_json(gate_path)
    review = _load_json(review_path)
    if candidate_gate.get("candidate_count") != EXPECTED_CANDIDATES:
        raise ValueError("Candidate gate does not declare 456 candidates")
    scientific_release = review.get("scientific_release")
    if not isinstance(scientific_release, dict) or (
        scientific_release.get("status")
        != "released_for_full_456_scan_implementation"
        or scientific_release.get("full_affinity_scan_release") != "pass"
        or scientific_release.get("full_scan_contract")
        != "score_all_declared_candidates_then_filter_once"
    ):
        raise ValueError("Pilot V2 scientific review does not release the full scan")
    assignments, shard_ids = build_full_scan_shards(candidates)

    output_dir = args.output_dir.expanduser().absolute()
    run_summary = args.run_summary.expanduser().absolute()
    output_names = {
        "manifest": "full_scan_shard_manifest.csv",
        "plan": "full_scan_plan.json",
        **{
            shard_id: f"{shard_id}_candidate_ids.txt"
            for shard_id in sorted(shard_ids)
        },
    }
    validated = validate_file_paths(
        project_root=PROJECT_ROOT,
        source_paths=[candidate_path, gate_path, review_path],
        target_paths=[
            *[output_dir / name for name in output_names.values()],
            run_summary,
        ],
    )
    final_paths = dict(zip(output_names, validated.target_paths[:-1], strict=True))
    run_summary = validated.target_paths[-1]
    existing = [path for path in [*final_paths.values(), run_summary] if path.exists()]
    if existing:
        raise FileExistsError(
            "Refusing to overwrite existing outputs:\n" + "\n".join(map(str, existing))
        )
    for path in [*final_paths.values(), run_summary]:
        path.parent.mkdir(parents=True, exist_ok=True)

    plan = {
        "schema_version": 1,
        "status": "pass",
        "generated_at": generated_at,
        "plan_name": "nb252_affinity_pyrosetta_full_scan",
        "candidate_count": EXPECTED_CANDIDATES,
        "replicate_count_per_candidate": EXPECTED_REPLICATES,
        "mutant_evaluation_count": EXPECTED_CANDIDATES * EXPECTED_REPLICATES,
        "shard_count": SHARD_COUNT,
        "candidates_per_shard": CANDIDATES_PER_SHARD,
        "positions_per_shard": 2,
        "maximum_array_concurrency": 4,
        "base_seed": 8112100,
        "candidate_filtering_applied": False,
        "full_scan_contract": "score_all_declared_candidates_then_filter_once",
        "shards": [
            {
                "shard_id": shard_id,
                "shard_index": int(shard_id.rsplit("_", 1)[1]),
                "candidate_id_file": output_names[shard_id],
                "candidate_count": len(candidate_ids),
            }
            for shard_id, candidate_ids in sorted(shard_ids.items())
        ],
    }
    with tempfile.TemporaryDirectory(prefix=".affinity-full-plan-", dir=PROJECT_ROOT) as temp:
        staging = Path(temp)
        staged = {key: staging / name for key, name in output_names.items()}
        staged_summary = staging / "run_summary.json"
        _write_csv(staged["manifest"], assignments, SHARD_MANIFEST_FIELDS)
        _write_json(staged["plan"], plan)
        for shard_id, candidate_ids in shard_ids.items():
            staged[shard_id].write_text(
                "".join(f"{candidate_id}\n" for candidate_id in candidate_ids),
                encoding="utf-8",
                newline="\n",
            )
        _write_json(
            staged_summary,
            {
                "schema_version": 1,
                "status": "pass",
                "generated_at": generated_at,
                "elapsed_seconds": round(time.perf_counter() - started, 6),
                "python": platform.python_version(),
                "command_argv": [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
                "counts": {
                    "candidate_count": EXPECTED_CANDIDATES,
                    "mutant_evaluation_count": EXPECTED_CANDIDATES
                    * EXPECTED_REPLICATES,
                    "shard_count": SHARD_COUNT,
                    "candidates_per_shard": CANDIDATES_PER_SHARD,
                },
                "candidate_filtering_applied": False,
                "outputs": {key: str(path) for key, path in final_paths.items()},
            },
        )
        replace_staged_files(
            {
                **{staged[key]: final_paths[key] for key in output_names},
                staged_summary: run_summary,
            },
            project_root=PROJECT_ROOT,
            protected_source_paths=validated.source_paths,
        )
    return 0


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
