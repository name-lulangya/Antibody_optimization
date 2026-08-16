#!/usr/bin/env python3
"""Analyze targeted PyRosetta output and build qualification gate V2."""

from __future__ import annotations

import argparse
import csv
import json
import platform
import sys
import tempfile
import time
from collections import Counter
from datetime import datetime
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from antibody_optimization.file_transaction import replace_staged_files, validate_file_paths  # noqa: E402
from antibody_optimization.targeted_structure_review import qualify_v2  # noqa: E402
from antibody_optimization.targeted_structure_review_plot import render_targeted_review  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--safety-review-dir", type=Path, required=True)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--runtime-result-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    safety_dir = _project_dir(args.safety_review_dir)
    plan_dir = _project_dir(args.plan_dir)
    runtime_dir = _project_dir(args.runtime_result_dir)
    sources = [
        safety_dir / "single_mutant_safety_review.csv",
        safety_dir / "single_mutant_safety_gate.json",
        plan_dir / "targeted_structure_review_candidates.csv",
        plan_dir / "targeted_structure_review_contract.json",
        runtime_dir / "targeted_structure_replicates.csv",
        runtime_dir / "targeted_structure_runtime_gate.json",
    ]
    names = {
        "review": "single_mutant_safety_review_v2.csv",
        "summary": "single_mutant_safety_summary_v2.csv",
        "gate": "single_mutant_safety_gate_v2.json",
        "png": "targeted_structure_review_v2.png",
        "svg": "targeted_structure_review_v2.svg",
    }
    output_dir = args.output_dir.expanduser().absolute()
    run_summary = args.run_summary.expanduser().absolute()
    validated = validate_file_paths(project_root=PROJECT_ROOT, source_paths=sources, target_paths=[*[output_dir / value for value in names.values()], run_summary])
    finals = dict(zip(names, validated.target_paths[:-1], strict=True))
    run_summary = validated.target_paths[-1]
    existing = [path for path in [*finals.values(), run_summary] if path.exists()]
    if existing:
        raise FileExistsError("Refusing to overwrite existing outputs:\n" + "\n".join(map(str, existing)))
    runtime_gate = _json(sources[5])
    result = qualify_v2(_csv(sources[0]), _csv(sources[2]), _csv(sources[4]), runtime_gate)
    counts = Counter(str(row["v2_qualification_status"]) for row in result["review_rows"])
    summary_rows = [{"status": status, "count": count} for status, count in sorted(counts.items())]
    gate = {
        "schema_version": 2,
        "status": "pass",
        "release": "ready_for_combination_module_review",
        "generated_at": generated_at,
        **result["facts"],
        "contract": {
            "exact_contact_set_equality_required": False,
            "existing_complex_evidence_reused_without_rerun": True,
            "contact_requirement": "preserve the NK2R epitope and binding conformation while allowing local contact changes within the original epitope",
            "new_runtime_requirement": "three AF3 VHH replicates with non-increasing median total score and local fa_rep",
            "nonresolvable_flags_are_not_cleared_by_repack": True,
            "intrinsic_hard_risks_are_do_not_advance": True,
            "combination_generated": False,
        },
        "interpretation": "V2 qualification is computational triage, not measured affinity, expression, aggregation, or stability. Combination compatibility has not yet been evaluated.",
    }
    review_fields = _field_union(result["review_rows"])
    for path in [*finals.values(), run_summary]: path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".targeted-review-analysis-", dir=PROJECT_ROOT) as temp:
        staging = Path(temp)
        staged = {key: staging / value for key, value in names.items()}
        staged_run = staging / "run_summary.json"
        _write_csv(staged["review"], result["review_rows"], review_fields)
        _write_csv(staged["summary"], summary_rows, ["status", "count"])
        _write_json(staged["gate"], gate)
        render_targeted_review(result["review_rows"], staged["png"], staged["svg"])
        _write_json(staged_run, {
            "schema_version": 1, "status": "pass", "generated_at": generated_at,
            "elapsed_seconds": round(time.perf_counter() - started, 6),
            "python": platform.python_version(),
            "command_argv": [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
            "candidate_count": len(result["review_rows"]),
            "combination_ready_count": counts.get("combination_ready", 0),
            "combination_generated": False,
            "outputs": {key: str(path) for key, path in finals.items()},
        })
        replace_staged_files({**{staged[key]: finals[key] for key in names}, staged_run: run_summary}, project_root=PROJECT_ROOT, protected_source_paths=validated.source_paths)
    return 0


def _field_union(rows: list[dict[str, object]]) -> list[str]:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields: fields.append(key)
    return fields


def _project_dir(path: Path) -> Path:
    resolved = path.expanduser().resolve(strict=True); resolved.relative_to(PROJECT_ROOT.resolve(strict=True))
    if not resolved.is_dir() or resolved.is_symlink(): raise ValueError(f"Expected regular project directory: {resolved}")
    return resolved


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle: return list(csv.DictReader(handle))


def _json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict): raise ValueError(f"Expected JSON object: {path}")
    return value


def _write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n"); writer.writeheader(); writer.writerows(rows)


def _write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    raise SystemExit(main())
