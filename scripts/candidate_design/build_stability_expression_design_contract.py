#!/usr/bin/env python3
"""Build the WT stability/expression discovery contract without mutations."""

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

from antibody_optimization.design_contract_plot import render_stability_expression_contract_figure  # noqa: E402
from antibody_optimization.file_transaction import replace_staged_files, validate_file_paths  # noqa: E402
from antibody_optimization.stability_expression_contract import build_stability_expression_contract  # noqa: E402

OUTPUT_NAMES = {
    "positions": "stability_expression_position_contract.csv",
    "contract": "stability_expression_design_contract.json",
    "gate": "stability_expression_contract_gate.json",
    "png": "stability_expression_contract.png",
    "svg": "stability_expression_contract.svg",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage0-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    parser.add_argument("--check_only", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    stage0 = args.stage0_dir.resolve(strict=True)
    sources = {
        "positions": stage0 / "mutable_position_inventory.csv",
        "stage2_contract": stage0 / "stage2_design_contract.json",
        "preflight": stage0 / "stage2_preflight.json",
    }
    positions = _load_csv(sources["positions"])
    stage2 = _load_json(sources["stage2_contract"])
    if _load_json(sources["preflight"]).get("status") != "pass":
        raise ValueError("Stage-2 preflight did not pass")
    result = build_stability_expression_contract(positions, stage2)
    if args.check_only:
        print(json.dumps({"status": "pass", "counts": result["counts"]}, ensure_ascii=False))
        return 0
    output_dir = args.output_dir.absolute()
    targets = [output_dir / name for name in OUTPUT_NAMES.values()] + [args.run_summary.absolute()]
    validated = validate_file_paths(project_root=PROJECT_ROOT, source_paths=list(sources.values()), target_paths=targets)
    if any(path.exists() for path in validated.target_paths):
        raise FileExistsError("Refusing to overwrite stability/expression contract outputs")
    for path in validated.target_paths:
        path.parent.mkdir(parents=True, exist_ok=True)
    final = dict(zip((*OUTPUT_NAMES, "run_summary"), validated.target_paths, strict=True))
    contract = {
        "schema_version": 1, "contract_name": "nb252_stability_expression_wt_discovery", "status": "pass",
        "generated_at": generated_at, "parent_contract": str(sources["stage2_contract"]),
        "position_counts": result["counts"], "module_rules": result["module_rules"],
        "frozen_scope": ["all_CDR_positions", "experimental_interface", "disulfide_cysteines", "terminal_SSGS"],
        "missing_framework_rule": "designable_cautiously_only_with_AF3_full_VHH_evidence",
        "conditional_background_rule": "install_each_affinity_core_then_freeze_it_and_recompute_allowed_framework_positions",
        "mutation_generation_performed": False,
    }
    gate = {
        "schema_version": 1, "gate_name": "nb252_stability_expression_design_contract", "status": "pass",
        "generated_at": generated_at, **result["counts"], "mutation_generation_performed": False,
        "release": "ready_for_nanobert_association_and_antifold_minimal_route_implementation",
        "interpretation": "This freezes design scope and evidence requirements; it does not predict stability, expression, or yield.",
    }
    with tempfile.TemporaryDirectory(prefix=".stability-expression-contract-", dir=PROJECT_ROOT) as temp:
        staging = Path(temp)
        staged = {key: staging / Path(path).name for key, path in final.items()}
        _write_csv(staged["positions"], result["position_rows"])
        _write_json(staged["contract"], contract); _write_json(staged["gate"], gate)
        render_stability_expression_contract_figure(result["position_rows"], png_path=staged["png"], svg_path=staged["svg"])
        _write_json(staged["run_summary"], {
            "schema_version": 1, "status": "pass", "generated_at": generated_at,
            "elapsed_seconds": round(time.perf_counter() - started, 6), "python": platform.python_version(),
            "command_argv": [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
            "counts": result["counts"], "mutation_generation_performed": False,
            "outputs": {key: str(path) for key, path in final.items() if key != "run_summary"},
        })
        replace_staged_files({staged[key]: final[key] for key in staged}, project_root=PROJECT_ROOT, protected_source_paths=validated.source_paths)
    return 0


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
