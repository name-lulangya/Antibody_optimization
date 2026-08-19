#!/usr/bin/env python3
"""Compare repeated NetSolP, NanoMelt, and AntiFold outputs with reusable legacy scores."""

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

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from antibody_optimization.antifold_validation import build_candidate_evidence, normalize_antifold_rows  # noqa: E402
from antibody_optimization.expression_property_completion import compare_repeat_scores  # noqa: E402
from antibody_optimization.file_transaction import replace_staged_files, validate_file_paths  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--antifold-plan-dir", type=Path, required=True)
    parser.add_argument("--netsolp-score-dir", type=Path, required=True)
    parser.add_argument("--nanomelt-score-dir", type=Path, required=True)
    parser.add_argument("--antifold-score-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    args = parser.parse_args()
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    plan_dir = args.plan_dir.resolve(strict=True)
    antifold_plan = args.antifold_plan_dir.resolve(strict=True)
    net_dir = args.netsolp_score_dir.resolve(strict=True)
    melt_dir = args.nanomelt_score_dir.resolve(strict=True)
    anti_dir = args.antifold_score_dir.resolve(strict=True)
    samples_path = plan_dir / "reuse_validation_samples.csv"
    expected_path = plan_dir / "reuse_validation_expected.csv"
    targets_path = plan_dir / "reuse_validation_antifold_targets.csv"
    contract_path = plan_dir / "expression_property_completion_contract.json"
    views_path = antifold_plan / "antifold_structure_views.csv"
    net_path = net_dir / "netsolp_sample_scores.csv"
    melt_path = melt_dir / "nanomelt_sample_scores.csv"
    anti_run_path = anti_dir / "antifold_model_run.json"
    sources = [samples_path, expected_path, targets_path, contract_path, views_path, net_path, melt_path, anti_run_path]
    views = {row["view_id"]: row for row in _csv(views_path)}
    indexed = {}
    for view_id, view in views.items():
        score_path = anti_dir / f"{view_id}.csv"
        score_path.resolve(strict=True)
        sources.append(score_path)
        indexed[view_id] = normalize_antifold_rows(_csv(score_path), view_id=view_id, vhh_chain=view["vhh_chain"])
    target_rows = _csv(targets_path)
    evidence, _ = build_candidate_evidence(target_rows, indexed)
    repeated_antifold = _wide_antifold(evidence)
    contract = _json(contract_path)
    tolerances = contract["repeat_validation"]
    comparisons, gate = compare_repeat_scores(
        _csv(samples_path), _csv(expected_path), _csv(net_path), _csv(melt_path), repeated_antifold,
        netsolp_tolerance=float(tolerances["netsolp_absolute_tolerance"]),
        nanomelt_tolerance=float(tolerances["nanomelt_tm_absolute_tolerance_c"]),
        antifold_tolerance=float(tolerances["antifold_absolute_tolerance"]),
    )
    gate = {**gate, "generated_at": generated_at}
    output_dir = args.output_dir.absolute()
    summary_path = args.run_summary.absolute()
    if output_dir.exists() or summary_path.exists():
        raise FileExistsError("Refusing to overwrite reuse-validation outputs")
    output_dir.parent.mkdir(parents=True, exist_ok=True); summary_path.parent.mkdir(parents=True, exist_ok=True)
    names = {"comparisons": "reuse_validation_comparisons.csv", "gate": "reuse_validation_gate.json"}
    with tempfile.TemporaryDirectory(prefix=".expression-property-validation-", dir=ROOT) as temp_name:
        staging = Path(temp_name)
        comparison_path = staging / names["comparisons"]
        gate_path = staging / names["gate"]
        summary_stage = staging / "run_summary.json"
        _write_csv(comparison_path, comparisons); _write_json(gate_path, gate)
        _write_json(summary_stage, {
            "schema_version": 1, "status": gate["status"], "generated_at": generated_at,
            "elapsed_seconds": round(time.perf_counter() - started, 6), "python": platform.python_version(),
            "counts": {"comparisons": len(comparisons), "failures": gate["failure_count"]},
            "release": gate["release"],
            "outputs": {key: str(output_dir / value) for key, value in names.items()},
        })
        pairs = {comparison_path: output_dir / names["comparisons"], gate_path: output_dir / names["gate"], summary_stage: summary_path}
        validate_file_paths(project_root=ROOT, source_paths=sources, target_paths=pairs.values())
        for target in pairs.values(): target.parent.mkdir(parents=True, exist_ok=True)
        replace_staged_files(pairs, project_root=ROOT, protected_source_paths=sources)
    if gate["status"] != "pass":
        raise RuntimeError("Historical score repeat-validation gate failed")
    return 0


def _wide_antifold(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    output: dict[str, dict[str, object]] = {}
    for row in rows:
        identifier = str(row["candidate_id"]); view = str(row["view_id"])
        wide = output.setdefault(identifier, {"candidate_id": identifier})
        for suffix in ("evaluation_status", "wt_log_probability", "mutant_log_probability", "delta_log_probability", "perplexity"):
            wide[f"{view}_{suffix}"] = row[suffix]
    return list(output.values())


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle: return list(csv.DictReader(handle))


def _json(path: Path) -> dict[str, object]: return json.loads(path.read_text(encoding="utf-8-sig"))


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n"); writer.writeheader(); writer.writerows(rows)


def _write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__": raise SystemExit(main())
