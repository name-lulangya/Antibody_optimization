#!/usr/bin/env python3
"""Uniformly tier the complete 456-candidate Nb252 affinity full scan."""

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

from antibody_optimization.affinity_post_scan import (  # noqa: E402
    tier_affinity_candidates,
)
from antibody_optimization.affinity_post_scan_plot import (  # noqa: E402
    render_affinity_post_scan_figure,
)
from antibody_optimization.file_transaction import (  # noqa: E402
    replace_staged_files,
    validate_file_paths,
)


OUTPUT_NAMES = {
    "candidates": "affinity_candidate_tiers.csv",
    "tier_summary": "affinity_tier_summary.csv",
    "review_pool": "affinity_review_pool.csv",
    "gate": "affinity_post_scan_gate.json",
    "figure_png": "affinity_post_scan_tiers.png",
    "figure_svg": "affinity_post_scan_tiers.svg",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full-scan-dir", type=Path, required=True)
    parser.add_argument("--candidate-dir", type=Path, required=True)
    parser.add_argument("--critical-residue-sets", type=Path, required=True)
    parser.add_argument("--calibration-gate", type=Path, required=True)
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
    full_scan_dir = _project_directory(args.full_scan_dir)
    candidate_dir = _project_directory(args.candidate_dir)
    sources = {
        "candidates": candidate_dir / "affinity_single_mutants.csv",
        "summaries": full_scan_dir / "candidate_summary.csv",
        "replicates": full_scan_dir / "candidate_replicate_metrics.csv",
        "merge_gate": full_scan_dir / "full_scan_merge_gate.json",
        "scientific_review": full_scan_dir / "affinity_full_scan_scientific_review.json",
        "critical_residue_sets": args.critical_residue_sets.expanduser().resolve(strict=True),
        "calibration_gate": args.calibration_gate.expanduser().resolve(strict=True),
    }
    result = tier_affinity_candidates(
        candidates=_load_csv(sources["candidates"]),
        summaries=_load_csv(sources["summaries"]),
        replicates=_load_csv(sources["replicates"]),
        merge_gate=_load_json(sources["merge_gate"]),
        scientific_review=_load_json(sources["scientific_review"]),
        critical_residue_sets=_load_json(sources["critical_residue_sets"]),
        calibration_gate=_load_json(sources["calibration_gate"]),
    )

    output_dir = args.output_dir.expanduser().absolute()
    run_summary = args.run_summary.expanduser().absolute()
    validated = validate_file_paths(
        project_root=PROJECT_ROOT,
        source_paths=list(sources.values()),
        target_paths=[
            *[output_dir / name for name in OUTPUT_NAMES.values()],
            run_summary,
        ],
    )
    final_paths = dict(zip(OUTPUT_NAMES, validated.target_paths[:-1], strict=True))
    run_summary = validated.target_paths[-1]
    existing = [path for path in [*final_paths.values(), run_summary] if path.exists()]
    if existing:
        raise FileExistsError(
            "Refusing to overwrite existing outputs:\n" + "\n".join(map(str, existing))
        )
    for path in [*final_paths.values(), run_summary]:
        path.parent.mkdir(parents=True, exist_ok=True)

    candidate_rows = result["candidate_rows"]
    tier_summary_rows = result["tier_summary_rows"]
    review_pool_rows = result["review_pool_rows"]
    gate = {
        "schema_version": 1,
        "gate_name": "nb252_affinity_post_scan_tiering",
        "status": "pass",
        "generated_at": generated_at,
        "candidate_count": len(candidate_rows),
        "tier_counts": result["tier_counts"],
        "strict_review_pool_count": len(review_pool_rows),
        "candidate_selection_performed": False,
        "combination_mutations_generated": False,
        "source_full_scan_complete": True,
        "all_candidates_assigned_exactly_one_tier": True,
        "thresholds": result["thresholds"],
        "release": "ready_for_strict_review_protocol_design",
        "interpretation": (
            "Evidence tiers organize paired-WT PyRosetta signals and structural QC; "
            "they are not measured affinity or a final experimental panel."
        ),
    }

    with tempfile.TemporaryDirectory(prefix=".affinity-post-scan-", dir=PROJECT_ROOT) as temp:
        staging = Path(temp)
        staged = {key: staging / name for key, name in OUTPUT_NAMES.items()}
        staged_summary = staging / "run_summary.json"
        _write_csv(staged["candidates"], candidate_rows)
        _write_csv(staged["tier_summary"], tier_summary_rows)
        _write_csv(staged["review_pool"], review_pool_rows)
        _write_json(staged["gate"], gate)
        render_affinity_post_scan_figure(
            rows=candidate_rows,
            png_path=staged["figure_png"],
            svg_path=staged["figure_svg"],
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
                    "candidate_count": len(candidate_rows),
                    "strict_review_pool_count": len(review_pool_rows),
                    "tier_counts": result["tier_counts"],
                },
                "candidate_selection_performed": False,
                "outputs": {key: str(path) for key, path in final_paths.items()},
            },
        )
        replace_staged_files(
            {
                **{staged[key]: final_paths[key] for key in OUTPUT_NAMES},
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


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write an empty CSV: {path}")
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
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
