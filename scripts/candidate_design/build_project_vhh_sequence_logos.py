#!/usr/bin/env python3
"""Build FR/CDR-annotated logos for natural, Nb252-neighbor, and project VHH sets."""

from __future__ import annotations

import argparse
import csv
import json
import os
import platform
import sys
import tempfile
import time
from collections import Counter
from datetime import datetime
from importlib import metadata as importlib_metadata
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))
os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "antibody_optimization_matplotlib")
)

from antibody_optimization.file_transaction import replace_staged_files, validate_file_paths
from antibody_optimization.vhh_conservation import (
    build_project_vhh_records,
    calculate_conservation,
)
from antibody_optimization.vhh_conservation_plot import render_frequency_logo


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sequence-review", type=Path, required=True)
    parser.add_argument("--sequence-positions", type=Path, required=True)
    parser.add_argument("--global-conservation", type=Path, required=True)
    parser.add_argument("--neighbor-conservation", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--expected-source-count", type=int, default=47)
    parser.add_argument("--expected-vhh-logo-count", type=int, default=45)
    parser.add_argument("--generated-at")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.perf_counter()
    generated = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    source_paths = (
        args.sequence_review,
        args.sequence_positions,
        args.global_conservation,
        args.neighbor_conservation,
    )
    for path in source_paths:
        path.resolve(strict=True)
    if args.output_dir.exists() or args.run_summary.exists():
        raise FileExistsError("Refusing to overwrite VHH sequence-logo outputs")

    review_rows = _csv(args.sequence_review)
    position_rows = _csv(args.sequence_positions)
    if len(review_rows) != args.expected_source_count:
        raise ValueError(
            f"Expected {args.expected_source_count} project sequences, found {len(review_rows)}"
        )
    project_records, project_audit = build_project_vhh_records(review_rows, position_rows)
    if len(project_records) != args.expected_vhh_logo_count:
        raise ValueError(
            f"Expected {args.expected_vhh_logo_count} eligible project VHHs, "
            f"found {len(project_records)}"
        )
    project_rows = calculate_conservation(
        project_records,
        subset_name="project_expression_numbered_heavy_chains",
    )
    global_rows = _csv(args.global_conservation)
    neighbor_rows = _csv(args.neighbor_conservation)
    exclusions = [row for row in project_audit if row["logo_status"] == "excluded"]
    counts = {
        "source_sequence_count": len(review_rows),
        "included_numbered_heavy_chain_count": len(project_records),
        "excluded_count": len(exclusions),
        "exclusion_reason_counts": dict(Counter(row["logo_reason"] for row in exclusions)),
        "global_logo_position_count": len(global_rows),
        "neighbor_logo_position_count": len(neighbor_rows),
        "project_logo_position_count": len(project_rows),
    }
    manifest = {
        "schema_version": 1,
        "status": "pass",
        "generated_at": generated,
        "imgt_region_definition": {
            "FR1": "1-26",
            "CDR1": "27-38",
            "FR2": "39-55",
            "CDR2": "56-65",
            "FR3": "66-104",
            "CDR3": "105-117",
            "FR4": "118-128",
        },
        "project_logo_scope": (
            "Equal-weight project sequences with frozen numbering_status=pass and chain_type=H. "
            "Yield values are not used as weights."
        ),
        "excluded_samples": exclusions,
        "counts": counts,
        "upstream_conservation_run": "vhh_conservation_20260818",
    }
    names = {
        "manifest": "sequence_logo_manifest.json",
        "audit": "project_sequence_logo_audit.csv",
        "project_frequencies": "project_expression_vhh_position_frequencies.csv",
        "global_png": "global_natural_vhh_sequence_logo_with_regions.png",
        "global_svg": "global_natural_vhh_sequence_logo_with_regions.svg",
        "neighbor_png": "nb252_neighbor_sequence_logo_with_regions.png",
        "neighbor_svg": "nb252_neighbor_sequence_logo_with_regions.svg",
        "project_png": "project_expression_vhh_sequence_logo_with_regions.png",
        "project_svg": "project_expression_vhh_sequence_logo_with_regions.svg",
    }
    output = args.output_dir.absolute()
    summary = args.run_summary.absolute()
    output.parent.mkdir(parents=True, exist_ok=True)
    summary.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".vhh-logos-", dir=PROJECT_ROOT) as tmp:
        stage = Path(tmp)
        staged = {key: stage / name for key, name in names.items()}
        _write_json(staged["manifest"], manifest)
        _write_csv(staged["audit"], project_audit)
        _write_csv(staged["project_frequencies"], project_rows)
        render_frequency_logo(
            global_rows,
            title="Global TNP natural VHH weighted residue frequencies",
            png_path=staged["global_png"],
            svg_path=staged["global_svg"],
        )
        render_frequency_logo(
            neighbor_rows,
            title="Nb252 framework-neighbor weighted residue frequencies",
            png_path=staged["neighbor_png"],
            svg_path=staged["neighbor_svg"],
        )
        render_frequency_logo(
            project_rows,
            title="Project expression-panel VHH residue frequencies (45 of 47 sequences)",
            png_path=staged["project_png"],
            svg_path=staged["project_svg"],
        )
        summary_stage = stage / "run_summary.json"
        _write_json(
            summary_stage,
            {
                "schema_version": 1,
                "status": "pass",
                "generated_at": generated,
                "elapsed_seconds": round(time.perf_counter() - started, 6),
                "python": platform.python_version(),
                "matplotlib": importlib_metadata.version("matplotlib"),
                "command_argv": [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
                "counts": counts,
                "outputs": {key: str(output / name) for key, name in names.items()},
            },
        )
        pairs = {staged[key]: output / name for key, name in names.items()}
        pairs[summary_stage] = summary
        for target in pairs.values():
            target.parent.mkdir(parents=True, exist_ok=True)
        validate_file_paths(
            project_root=PROJECT_ROOT,
            source_paths=source_paths,
            target_paths=pairs.values(),
        )
        replace_staged_files(
            pairs,
            project_root=PROJECT_ROOT,
            protected_source_paths=source_paths,
        )
    return 0


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
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
