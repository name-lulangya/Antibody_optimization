#!/usr/bin/env python3
"""Independently review completed pilot and full property-affinity scans."""

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

from antibody_optimization.property_affinity_plot import render_scientific_review  # noqa: E402
from antibody_optimization.property_affinity_review import review_completed_scan  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot-result-dir", type=Path, required=True)
    parser.add_argument("--full-result-dir", type=Path, required=True)
    parser.add_argument("--pilot-run-summary", type=Path, required=True)
    parser.add_argument("--full-run-summary", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    pilot_dir = _project_dir(args.pilot_result_dir)
    full_dir = _project_dir(args.full_result_dir)
    pilot_run = _json(args.pilot_run_summary)
    full_run = _json(args.full_run_summary)
    pilot_gate = _json(pilot_dir / "property_affinity_scoring_gate.json")
    full_gate = _json(full_dir / "property_affinity_scoring_gate.json")
    if (
        pilot_gate.get("status") != "pass"
        or pilot_gate.get("release") != "ready_for_full_property_affinity_scan"
        or full_gate.get("status") != "pass"
        or full_gate.get("release") != "ready_for_property_affinity_noninferiority_review"
        or pilot_run.get("status") != "pass"
        or full_run.get("status") != "pass"
    ):
        raise ValueError("Pilot/full gates and run summaries must all be passed")

    review_rows, facts = review_completed_scan(
        summary_rows=_csv(full_dir / "property_affinity_candidate_summary.csv"),
        paired_rows=_csv(full_dir / "property_affinity_candidate_replicates.csv"),
        wt_rows=_csv(full_dir / "property_affinity_wt_controls.csv"),
        pilot_summary_rows=_csv(pilot_dir / "property_affinity_candidate_summary.csv"),
    )
    favorable = [row for row in review_rows if row["affinity_direction_class"] == "directionally_favorable"]
    intersections = [row for row in favorable if row["multitool_intersection_class"] == "rossetta_favorable_antifold_positive"]
    contact_changes = [row for row in review_rows if row["paired_contact_status"] != "preserved_all"]
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    scientific_review = {
        "schema_version": 1,
        "generated_at": generated_at,
        "status": "pass",
        "release": "ready_for_property_module_shortlist_review",
        "provenance": {
            "pilot_result_dir": str(pilot_dir.relative_to(PROJECT_ROOT)),
            "full_result_dir": str(full_dir.relative_to(PROJECT_ROOT)),
            "analysis_script": str(Path(__file__).relative_to(PROJECT_ROOT)),
            "score_function": full_gate["score_function"],
            "selected_protocol": full_gate["selected_protocol"],
            "score_semantics": full_gate["score_semantics"],
        },
        "runtime": {
            "pilot_elapsed_seconds": pilot_run["elapsed_seconds"],
            "full_scan_elapsed_seconds": full_run["elapsed_seconds"],
        },
        "validation": facts,
        "classification_rule": {
            "directionally_favorable": "both median deltas < 0 and each metric is < 0 in at least 2 of 3 replicates",
            "directionally_adverse": "both median deltas > 0 and each metric is < 0 in at most 1 of 3 replicates",
            "mixed": "all remaining sign patterns",
            "threshold_note": "direction-only descriptive classes; no REU effect-size or experimental non-inferiority threshold is claimed",
        },
        "results": {
            "directionally_favorable_candidates": [row["short_mutation"] for row in favorable],
            "all_three_replicates_both_energy_negative_candidates": [row["short_mutation"] for row in favorable if row["all_three_replicates_both_energy_negative"]],
            "rossetta_favorable_antifold_positive_candidates": [row["short_mutation"] for row in intersections],
            "contact_change_candidates": [
                {"mutation": row["short_mutation"], "lost_receptor_auth_positions": row["lost_receptor_auth_positions"]}
                for row in contact_changes
            ],
        },
        "decision": {
            "scientific_selection_performed": False,
            "candidate_filtering_performed": False,
            "allowed_use": "review the 9 directionally favorable property modules and their independent AntiFold/property/contact tradeoffs",
            "disallowed_use": "claim measured affinity, convert REU to KD, or treat the 9 candidates as final experimental selections",
        },
    }

    output_dir = args.output_dir.expanduser().absolute()
    run_summary = args.run_summary.expanduser().absolute()
    if output_dir.exists() or run_summary.exists():
        raise FileExistsError("Refusing to overwrite an existing scientific review")
    output_dir.mkdir(parents=True)
    run_summary.parent.mkdir(parents=True, exist_ok=True)
    review_csv = output_dir / "property_affinity_scientific_review.csv"
    summary_csv = output_dir / "property_affinity_scientific_review_summary.csv"
    review_json = output_dir / "property_affinity_scientific_review.json"
    png = output_dir / "property_affinity_scientific_review.png"
    svg = output_dir / "property_affinity_scientific_review.svg"
    _write_csv(review_csv, review_rows, list(review_rows[0]))
    summary_rows = [
        {"metric": "candidate_count", "value": facts["candidate_count"]},
        *[
            {"metric": f"direction_{key}", "value": value}
            for key, value in sorted(facts["direction_class_counts"].items())
        ],
        {"metric": "all_three_both_energy_negative_count", "value": facts["all_three_both_energy_negative_count"]},
        {"metric": "rossetta_favorable_antifold_positive_count", "value": len(intersections)},
        {"metric": "paired_contact_preserved_all_count", "value": facts["paired_contact_preserved_all_count"]},
        {"metric": "minimum_paired_receptor_contact_retention", "value": facts["minimum_paired_receptor_contact_retention"]},
        {"metric": "maximum_interface_ca_rmsd_angstrom", "value": facts["maximum_interface_ca_rmsd_angstrom"]},
        {"metric": "pilot_full_max_absolute_difference", "value": facts["pilot_full_max_absolute_difference"]},
        {"metric": "pilot_elapsed_seconds", "value": pilot_run["elapsed_seconds"]},
        {"metric": "full_scan_elapsed_seconds", "value": full_run["elapsed_seconds"]},
    ]
    _write_csv(summary_csv, summary_rows, ["metric", "value"])
    _write_json(review_json, scientific_review)
    with tempfile.TemporaryDirectory(prefix="nb252-review-mpl-") as mpl_config:
        os.environ.setdefault("MPLCONFIGDIR", mpl_config)
        render_scientific_review(review_rows, png, svg)
    _write_json(
        run_summary,
        {
            "schema_version": 1,
            "status": "pass",
            "generated_at": generated_at,
            "script": str(Path(__file__).relative_to(PROJECT_ROOT)),
            "candidate_count": facts["candidate_count"],
            "direction_class_counts": facts["direction_class_counts"],
            "scientific_selection_performed": False,
            "outputs": [str(path.relative_to(PROJECT_ROOT)) for path in (review_csv, summary_csv, review_json, png, svg)],
        },
    )
    return 0


def _project_dir(path: Path) -> Path:
    resolved = path.expanduser().resolve(strict=True)
    resolved.relative_to(PROJECT_ROOT.resolve(strict=True))
    if not resolved.is_dir() or resolved.is_symlink():
        raise ValueError(f"Expected regular project directory: {resolved}")
    return resolved


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path: Path) -> dict[str, object]:
    value = json.loads(path.expanduser().resolve(strict=True).read_text(encoding="utf-8-sig"))
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
