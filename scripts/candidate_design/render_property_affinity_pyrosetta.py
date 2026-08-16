#!/usr/bin/env python3
"""Render a completed pilot or full property-affinity PyRosetta result."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import datetime
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from antibody_optimization.property_affinity_plot import render_scoring  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-kind", choices=("pilot", "full_scan"), required=True)
    parser.add_argument("--result-dir", type=Path, required=True)
    args = parser.parse_args()
    result_dir = args.result_dir.resolve(strict=True)
    result_dir.relative_to(PROJECT_ROOT.resolve(strict=True))
    gate = json.loads((result_dir / "property_affinity_scoring_gate.json").read_text(encoding="utf-8-sig"))
    if gate.get("status") != "pass" or gate.get("run_kind") != args.run_kind:
        raise ValueError("Scoring gate is not a matching passed run")
    with (result_dir / "property_affinity_candidate_summary.csv").open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    png = result_dir / "property_affinity_scoring.png"
    svg = result_dir / "property_affinity_scoring.svg"
    manifest = result_dir / "property_affinity_scoring_plot_manifest.json"
    if png.exists() or svg.exists() or manifest.exists():
        raise FileExistsError("Refusing to overwrite existing scoring figures")
    render_scoring(rows, png, svg, args.run_kind)
    manifest.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
                "run_kind": args.run_kind,
                "plot_script": str(Path(__file__).relative_to(PROJECT_ROOT)),
                "exact_plot_data": str((result_dir / "property_affinity_candidate_summary.csv").relative_to(PROJECT_ROOT)),
                "upstream_gate": str((result_dir / "property_affinity_scoring_gate.json").relative_to(PROJECT_ROOT)),
                "outputs": [str(png.relative_to(PROJECT_ROOT)), str(svg.relative_to(PROJECT_ROOT))],
                "energy_semantics": "mutant_minus_position_specific_paired_WT_Rosetta_REU",
                "uncertainty_semantics": "three_replicates_summarized_by_median; contact_panel_shows_minimum",
                "candidate_filtering_applied_during_scoring": False,
            },
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        ) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
