#!/usr/bin/env python3
"""Merge validated legacy scores with the 126 newly completed scores into 847 rows."""

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

from antibody_optimization.expression_property_completion import build_complete_score_matrix  # noqa: E402
from antibody_optimization.file_transaction import replace_staged_files, validate_file_paths  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v2-contract-dir", type=Path, required=True)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--validation-dir", type=Path, required=True)
    parser.add_argument("--legacy-property-dir", type=Path, required=True)
    parser.add_argument("--legacy-antifold-dir", type=Path, required=True)
    parser.add_argument("--netsolp-score-dir", type=Path, required=True)
    parser.add_argument("--nanomelt-score-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    args = parser.parse_args()
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    v2_dir = args.v2_contract_dir.resolve(strict=True); plan_dir = args.plan_dir.resolve(strict=True)
    validation_dir = args.validation_dir.resolve(strict=True)
    legacy_property_dir = args.legacy_property_dir.resolve(strict=True)
    legacy_antifold_dir = args.legacy_antifold_dir.resolve(strict=True)
    net_dir = args.netsolp_score_dir.resolve(strict=True); melt_dir = args.nanomelt_score_dir.resolve(strict=True)
    candidates_path = v2_dir / "nb252_allowed_single_mutants.csv"
    audit_path = plan_dir / "reuse_mapping_audit.csv"
    contract_path = plan_dir / "expression_property_completion_contract.json"
    validation_gate_path = validation_dir / "reuse_validation_gate.json"
    old_property_path = legacy_property_dir / "unified_single_mutant_property_evidence.csv"
    old_antifold_path = legacy_antifold_dir / "unified_single_mutant_antifold_evidence.csv"
    net_path = net_dir / "netsolp_sample_scores.csv"
    melt_path = melt_dir / "nanomelt_sample_scores.csv"
    sources = [candidates_path, audit_path, contract_path, validation_gate_path, old_property_path, old_antifold_path, net_path, melt_path]
    current = _csv(candidates_path); contract = _json(contract_path)
    matrix, gate = build_complete_score_matrix(
        current, str(contract["authoritative_parent"]["sequence"]), _csv(audit_path),
        _csv(old_property_path), _csv(old_antifold_path), _csv(net_path), _csv(melt_path),
        _json(validation_gate_path),
    )
    gate = {**gate, "generated_at": generated_at}
    output_dir = args.output_dir.absolute(); summary_path = args.run_summary.absolute()
    if output_dir.exists() or summary_path.exists(): raise FileExistsError("Refusing to overwrite complete property matrix")
    output_dir.parent.mkdir(parents=True, exist_ok=True); summary_path.parent.mkdir(parents=True, exist_ok=True)
    names = {"matrix": "expression_single_mutant_property_matrix.csv", "gate": "expression_single_mutant_property_matrix_gate.json"}
    with tempfile.TemporaryDirectory(prefix=".expression-property-final-", dir=ROOT) as temp_name:
        staging = Path(temp_name); matrix_path = staging / names["matrix"]; gate_path = staging / names["gate"]
        summary_stage = staging / "run_summary.json"
        _write_csv(matrix_path, matrix); _write_json(gate_path, gate)
        _write_json(summary_stage, {
            "schema_version": 1, "status": gate["status"], "generated_at": generated_at,
            "elapsed_seconds": round(time.perf_counter() - started, 6), "python": platform.python_version(),
            "counts": {"candidates": len(matrix), **gate["antifold_scope_counts"], **gate["property_score_source_counts"]},
            "release": gate["release"],
            "candidate_selection_performed": False,
            "outputs": {key: str(output_dir / value) for key, value in names.items()},
        })
        pairs = {matrix_path: output_dir / names["matrix"], gate_path: output_dir / names["gate"], summary_stage: summary_path}
        validate_file_paths(project_root=ROOT, source_paths=sources, target_paths=pairs.values())
        for target in pairs.values(): target.parent.mkdir(parents=True, exist_ok=True)
        replace_staged_files(pairs, project_root=ROOT, protected_source_paths=sources)
    if gate["status"] != "pass": raise RuntimeError("Complete property matrix gate failed")
    return 0


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle: return list(csv.DictReader(handle))


def _json(path: Path) -> dict[str, object]: return json.loads(path.read_text(encoding="utf-8-sig"))


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n"); writer.writeheader(); writer.writerows(rows)


def _write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__": raise SystemExit(main())
