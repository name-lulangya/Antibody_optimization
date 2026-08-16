#!/usr/bin/env python3
"""Narrow the released V2 Nb252 single-mutant pool using integrated risks."""

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

from antibody_optimization.file_transaction import replace_staged_files, validate_file_paths  # noqa: E402
from antibody_optimization.single_mutant_shortlist import build_single_mutant_shortlist  # noqa: E402
from antibody_optimization.single_mutant_shortlist_plot import render_single_mutant_shortlist  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v2-review", type=Path, required=True)
    parser.add_argument("--v2-gate", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    return parser.parse_args()


def main() -> int:
    args = parse_args(); started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    source_review = args.v2_review.expanduser().resolve(strict=True)
    source_gate = args.v2_gate.expanduser().resolve(strict=True)
    gate = _json(source_gate)
    if gate.get("status") != "pass" or gate.get("release") != "ready_for_combination_module_review":
        raise ValueError("V2 gate does not release shortlist review")
    output_dir = args.output_dir.expanduser().absolute(); run_summary = args.run_summary.expanduser().absolute()
    names = {
        "review": "single_mutant_shortlist_review.csv",
        "shortlist": "single_mutant_shortlist.csv",
        "summary": "single_mutant_shortlist_summary.csv",
        "plot_data": "single_mutant_shortlist_plot_data.csv",
        "gate": "single_mutant_shortlist_gate.json",
        "png": "single_mutant_shortlist.png",
        "svg": "single_mutant_shortlist.svg",
    }
    validated = validate_file_paths(
        project_root=PROJECT_ROOT,
        source_paths=[source_review, source_gate],
        target_paths=[*[output_dir / name for name in names.values()], run_summary],
    )
    finals = dict(zip(names, validated.target_paths[:-1], strict=True)); run_summary = validated.target_paths[-1]
    existing = [path for path in [*finals.values(), run_summary] if path.exists()]
    if existing: raise FileExistsError("Refusing to overwrite existing outputs:\n" + "\n".join(map(str, existing)))
    result = build_single_mutant_shortlist(_csv(source_review))
    output_gate = {
        "schema_version": 1, "status": "pass", "release": "ready_for_small_combination_contract",
        "generated_at": generated_at, **result["facts"],
        "contract": {
            "source_v2_gate": str(source_gate),
            "strong_negative_antifold_uses_existing_v2_flag": True,
            "affinity_candidate_states_changed": False,
            "historical_candidates_deleted": False,
            "combination_generated": False,
        },
        "interpretation": "The shortlist is computational prioritization, not measured affinity, expression, aggregation, or stability.",
    }
    summary_rows = [
        {"metric": "active_before", "count": result["facts"]["active_before_count"]},
        {"metric": "active_after", "count": result["facts"]["active_after_count"]},
        {"metric": "property_deprioritized", "count": result["facts"]["property_deprioritized_count"]},
    ]
    for path in [*finals.values(), run_summary]: path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".single-mutant-shortlist-", dir=PROJECT_ROOT) as temp:
        stage = Path(temp); staged = {key: stage / name for key, name in names.items()}; staged_run = stage / "run_summary.json"
        _write_csv(staged["review"], result["review_rows"], _fields(result["review_rows"]))
        _write_csv(staged["shortlist"], result["shortlist_rows"], _fields(result["shortlist_rows"]))
        _write_csv(staged["summary"], summary_rows, ["metric", "count"])
        plot_rows = render_single_mutant_shortlist(result["review_rows"], staged["png"], staged["svg"])
        _write_csv(staged["plot_data"], plot_rows, _fields(plot_rows))
        _write_json(staged["gate"], output_gate)
        _write_json(staged_run, {
            "schema_version": 1, "status": "pass", "generated_at": generated_at,
            "elapsed_seconds": round(time.perf_counter() - started, 6), "python": platform.python_version(),
            "command_argv": [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
            "active_before_count": result["facts"]["active_before_count"],
            "active_after_count": result["facts"]["active_after_count"],
            "property_deprioritized_count": result["facts"]["property_deprioritized_count"],
            "combination_generated": False,
            "outputs": {key: str(path) for key, path in finals.items()},
        })
        replace_staged_files(
            {**{staged[key]: finals[key] for key in names}, staged_run: run_summary},
            project_root=PROJECT_ROOT, protected_source_paths=validated.source_paths,
        )
    return 0


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle: return list(csv.DictReader(handle))


def _json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict): raise ValueError(f"Expected JSON object: {path}")
    return value


def _fields(rows: list[dict[str, object]]) -> list[str]:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields: fields.append(key)
    return fields


def _write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n", extrasaction="raise")
        writer.writeheader(); writer.writerows(rows)


def _write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    raise SystemExit(main())
