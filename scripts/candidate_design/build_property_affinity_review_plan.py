#!/usr/bin/env python3
"""Build the fixed 30-candidate property-affinity review pool locally."""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import tempfile
from datetime import datetime
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from antibody_optimization.property_affinity_plot import render_pool  # noqa: E402
from antibody_optimization.property_affinity_review import (  # noqa: E402
    MUTATION_NEIGHBORHOOD_ANGSTROM,
    PILOT_CANDIDATES,
    POOL_SIZE,
    REPLICATES,
    build_review_pool,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tnp-result-dir", type=Path, required=True)
    parser.add_argument("--unified-plan-dir", type=Path, required=True)
    parser.add_argument("--structure-baseline-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    for path in (args.tnp_result_dir, args.unified_plan_dir, args.structure_baseline_dir):
        path.resolve(strict=True).relative_to(PROJECT_ROOT.resolve(strict=True))
    args.output_dir = args.output_dir.expanduser().absolute()
    args.run_summary = args.run_summary.expanduser().absolute()
    if args.output_dir.exists():
        raise FileExistsError(f"Output directory already exists: {args.output_dir}")
    if args.run_summary.exists():
        raise FileExistsError(f"Run summary already exists: {args.run_summary}")
    review = _json(args.tnp_result_dir / "unified_tnp_scientific_review.json")
    gate = _json(args.tnp_result_dir / "unified_tnp_review_gate.json")
    decision = review.get("decision", {})
    if gate.get("status") != "pass" or not isinstance(decision, dict) or decision.get("status") != "pass" or decision.get("release") != "ready_for_multitool_shortlist_with_tnp_as_supporting_risk_only":
        raise ValueError("Unified TNP review is not scientifically released")
    rows = build_review_pool(
        _csv(args.tnp_result_dir / "unified_tnp_candidate_evidence.csv"),
        _csv(args.unified_plan_dir / "unified_single_mutant_candidates.csv"),
        _csv(args.structure_baseline_dir / "nb252_sequence_structure_mapping.csv"),
    )
    args.output_dir.mkdir(parents=True)
    fields = list(rows[0])
    _write_csv(args.output_dir / "property_affinity_review_candidates.csv", rows, fields)
    (args.output_dir / "pilot_candidate_ids.txt").write_text("\n".join(PILOT_CANDIDATES) + "\n", encoding="utf-8", newline="\n")
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    contract = {
        "schema_version": 1,
        "status": "pass",
        "release": "ready_for_property_affinity_pyrosetta_pilot",
        "generated_at": generated_at,
        "candidate_count": len(rows),
        "position_count": len({int(row["sequence_index_1based"]) for row in rows}),
        "pilot_candidate_count": len(PILOT_CANDIDATES),
        "replicates": REPLICATES,
        "run_kinds": ["pilot", "full_scan"],
        "mutation_neighborhood_angstrom": MUTATION_NEIGHBORHOOD_ANGSTROM,
        "movable_residue_definition": "calibrated_interface_neighborhood_union_mutation_8A_neighborhood",
        "paired_wt_semantics": "one_WT_per_position_replicate_shared_by_substitutions_at_that_position",
        "candidate_filtering_applied_during_scoring": False,
        "pilot_is_protocol_and_runtime_validation_only": True,
        "full_scan_requires_pilot_gate": "pass_ready_for_full_property_affinity_scan",
        "runtime_estimate": {"pilot": "10-25 minutes", "full_scan": "45-120 minutes", "exceeds_five_hours": False},
    }
    _write_json(args.output_dir / "property_affinity_review_contract.json", contract)
    with tempfile.TemporaryDirectory(prefix="nb252-mpl-") as mpl_config:
        os.environ.setdefault("MPLCONFIGDIR", mpl_config)
        render_pool(rows, args.output_dir / "property_affinity_review_pool.png", args.output_dir / "property_affinity_review_pool.svg")
    args.run_summary.parent.mkdir(parents=True, exist_ok=True)
    _write_json(args.run_summary, {**contract, "script": str(Path(__file__).relative_to(PROJECT_ROOT)), "outputs": [str(p.relative_to(PROJECT_ROOT)) for p in sorted(args.output_dir.iterdir())]})
    return 0


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def _write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def _write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    raise SystemExit(main())
