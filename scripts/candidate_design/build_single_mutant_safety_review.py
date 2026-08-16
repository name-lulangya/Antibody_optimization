#!/usr/bin/env python3
"""Build the local Nb252 single-mutant safety and combination-qualification review."""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import sys
import tempfile
from datetime import datetime
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SRC = PROJECT_ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from antibody_optimization.single_mutant_safety import review_single_mutant_modules  # noqa: E402
from antibody_optimization.single_mutant_safety_plot import render_single_mutant_safety_review  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--affinity-evidence", type=Path, required=True)
    parser.add_argument("--property-review", type=Path, required=True)
    parser.add_argument("--unified-candidates", type=Path, required=True)
    parser.add_argument("--property-evidence", type=Path, required=True)
    parser.add_argument("--antifold-core-evidence", type=Path, required=True)
    parser.add_argument("--tnp-evidence", type=Path, required=True)
    parser.add_argument("--prepared-wt-pdb", type=Path, required=True)
    parser.add_argument("--critical-residue-sets", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at", default="")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output = args.output_dir.resolve()
    summary = args.run_summary.resolve()
    if output.exists() or summary.exists():
        raise FileExistsError("Output directory and run summary must not already exist")
    critical = _json(args.critical_residue_sets)
    missing = set(critical["experimental_missing_coordinates"]["reported_sequence_indices_1based"])
    result = review_single_mutant_modules(
        _csv(args.affinity_evidence),
        _csv(args.property_review),
        _csv(args.unified_candidates),
        _csv(args.property_evidence),
        _csv(args.antifold_core_evidence),
        _csv(args.tnp_evidence),
        args.prepared_wt_pdb.resolve(),
        missing,
    )
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    stage = Path(tempfile.mkdtemp(prefix=f".{output.name}-stage-", dir=output.parent))
    try:
        review_csv = stage / "single_mutant_safety_review.csv"
        summary_csv = stage / "single_mutant_safety_summary.csv"
        gate_json = stage / "single_mutant_safety_gate.json"
        png = stage / "single_mutant_safety_review.png"
        svg = stage / "single_mutant_safety_review.svg"
        _write_csv(review_csv, result["review_rows"])
        summary_rows = _summary_rows(result["facts"])
        _write_csv(summary_csv, summary_rows)
        gate = {
            "schema_version": 1,
            "status": "pass",
            "release": "ready_for_targeted_structure_review_not_combination_generation",
            "generated_at": generated_at,
            "facts": result["facts"],
            "contract": {
                "antifold_strong_negative_operational_threshold": -3.0,
                "exposed_hydrophobic_relative_sasa_threshold": 0.25,
                "targeted_affinity_alternative_support_floor": "13_of_20_for_each_metric_with_both_negative_medians",
                "threshold_scope": "project-specific conservative triage; not validated universal developability cutoffs",
                "property_affinity_requirement": "not directionally adverse; improvement is not required",
                "combination_mutations_generated": False,
            },
            "provenance": {
                "prepared_wt_pdb": str(args.prepared_wt_pdb),
                "affinity_evidence": str(args.affinity_evidence),
                "property_review": str(args.property_review),
                "unified_candidates": str(args.unified_candidates),
                "property_evidence": str(args.property_evidence),
                "antifold_core_evidence": str(args.antifold_core_evidence),
                "tnp_evidence": str(args.tnp_evidence),
                "critical_residue_sets": str(args.critical_residue_sets),
            },
            "interpretation": (
                "Combination qualification is a computational expert-review state. "
                "It is not measured affinity, expression, aggregation, glycosylation, or experimental validation."
            ),
        }
        _write_json(gate_json, gate)
        render_single_mutant_safety_review(result["review_rows"], png, svg)
        os.replace(stage, output)
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    summary.parent.mkdir(parents=True, exist_ok=True)
    _write_json(
        summary,
        {
            "schema_version": 1,
            "status": "pass",
            "generated_at": generated_at,
            "script": str(Path(__file__).relative_to(PROJECT_ROOT)),
            "candidate_count": result["facts"]["candidate_count"],
            "qualification_counts": result["facts"]["qualification_counts"],
            "combination_generated": False,
            "output_dir": str(args.output_dir),
        },
    )
    return 0


def _summary_rows(facts: dict[str, object]) -> list[dict[str, object]]:
    rows = [
        {"scope": "all", "status": status, "count": count}
        for status, count in sorted(facts["qualification_counts"].items())
    ]
    for track, counts in facts["qualification_counts_by_track"].items():
        rows.extend(
            {"scope": track, "status": status, "count": count}
            for status, count in sorted(counts.items())
        )
    return rows


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Cannot write empty CSV: {path}")
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, value: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    raise SystemExit(main())
