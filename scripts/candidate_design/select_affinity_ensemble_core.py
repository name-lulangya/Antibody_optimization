#!/usr/bin/env python3
"""Select the interpretable affinity core from the completed 20-sample ensemble."""

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

from antibody_optimization.affinity_ensemble import select_affinity_core_modules  # noqa: E402
from antibody_optimization.design_contract_plot import render_affinity_ensemble_figure  # noqa: E402
from antibody_optimization.file_transaction import replace_staged_files, validate_file_paths  # noqa: E402

OUTPUT_NAMES = {
    "evidence": "affinity_ensemble_evidence.csv",
    "cores": "affinity_core_modules.csv",
    "positions": "affinity_core_positions.csv",
    "gate": "affinity_ensemble_core_gate.json",
    "png": "affinity_ensemble_core.png",
    "svg": "affinity_ensemble_core.svg",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--production-result-dir", type=Path, required=True)
    parser.add_argument("--post-scan-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    parser.add_argument("--check_only", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    production = args.production_result_dir.resolve(strict=True)
    post_scan = args.post_scan_dir.resolve(strict=True)
    sources = {
        "ensemble": production / "flex_ddg_production_candidate_summary.csv",
        "ensemble_gate": production / "flex_ddg_production_gate.json",
        "ensemble_review": production / "flex_ddg_production_scientific_review.json",
        "post_scan": post_scan / "affinity_candidate_tiers.csv",
        "post_scan_gate": post_scan / "affinity_post_scan_gate.json",
    }
    inputs = {key: (_load_csv(path) if path.suffix == ".csv" else _load_json(path)) for key, path in sources.items()}
    _validate_upstream(inputs)
    result = select_affinity_core_modules(inputs["ensemble"], inputs["post_scan"])
    if args.check_only:
        print(json.dumps({"status": "pass", "counts": result["counts"]}, ensure_ascii=False))
        return 0
    output_dir = args.output_dir.absolute()
    targets = [output_dir / name for name in OUTPUT_NAMES.values()] + [args.run_summary.absolute()]
    validated = validate_file_paths(project_root=PROJECT_ROOT, source_paths=list(sources.values()), target_paths=targets)
    if any(path.exists() for path in validated.target_paths):
        raise FileExistsError("Refusing to overwrite affinity ensemble outputs")
    for path in validated.target_paths:
        path.parent.mkdir(parents=True, exist_ok=True)
    final = dict(zip((*OUTPUT_NAMES, "run_summary"), validated.target_paths, strict=True))
    gate = {
        "schema_version": 1, "gate_name": "nb252_affinity_ensemble_core", "status": "pass",
        "generated_at": generated_at, **result["counts"],
        "selection_definition": "both_metrics_negative_in_at_least_18_of_20_and_negative_medians",
        "weighted_composite_score_used": False, "candidate_selection_performed": True,
        "combination_mutations_generated": False, "same_position_variants_mutually_exclusive": True,
        "release": "ready_for_affinity_core_module_use",
        "interpretation": "PyRosetta ensemble ranking evidence, not measured affinity; retained risk columns must be reviewed during combination design.",
    }
    with tempfile.TemporaryDirectory(prefix=".affinity-ensemble-", dir=PROJECT_ROOT) as temp:
        staging = Path(temp)
        staged = {key: staging / Path(path).name for key, path in final.items()}
        _write_csv(staged["evidence"], result["evidence_rows"])
        _write_csv(staged["cores"], result["core_rows"])
        _write_csv(staged["positions"], result["position_rows"])
        _write_json(staged["gate"], gate)
        render_affinity_ensemble_figure(result["evidence_rows"], png_path=staged["png"], svg_path=staged["svg"])
        _write_json(staged["run_summary"], {
            "schema_version": 1, "status": "pass", "generated_at": generated_at,
            "elapsed_seconds": round(time.perf_counter() - started, 6), "python": platform.python_version(),
            "command_argv": [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
            "counts": result["counts"], "candidate_selection_performed": True,
            "outputs": {key: str(path) for key, path in final.items() if key != "run_summary"},
        })
        replace_staged_files({staged[key]: final[key] for key in staged}, project_root=PROJECT_ROOT, protected_source_paths=validated.source_paths)
    return 0


def _validate_upstream(inputs: dict[str, object]) -> None:
    if inputs["ensemble_gate"].get("status") != "pass" or inputs["ensemble_gate"].get("release") != "ready_for_ensemble_filter_design":
        raise ValueError("Flex ddG production is not released for ensemble filtering")
    if inputs["ensemble_review"].get("integrity_review", {}).get("status") != "pass" or inputs["post_scan_gate"].get("status") != "pass":
        raise ValueError("Upstream scientific review or post-scan gate did not pass")
    if inputs["post_scan_gate"].get("candidate_selection_performed") is not False:
        raise ValueError("Post-scan input must remain unselected")


def _load_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _load_json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def _write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    raise SystemExit(main())
