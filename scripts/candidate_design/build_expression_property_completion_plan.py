#!/usr/bin/env python3
"""Build the v2 reuse audit, repeat-validation panel, and 126-score completion plan."""

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

from antibody_optimization.expression_property_completion import build_reuse_plan  # noqa: E402
from antibody_optimization.file_transaction import replace_staged_files, validate_file_paths  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v2-contract-dir", type=Path, required=True)
    parser.add_argument("--preflight", type=Path, required=True)
    parser.add_argument("--legacy-property-dir", type=Path, required=True)
    parser.add_argument("--legacy-antifold-dir", type=Path, required=True)
    parser.add_argument("--legacy-property-plan-dir", type=Path, required=True)
    parser.add_argument("--antifold-plan-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    parser.add_argument("--check_only", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    contract_dir = args.v2_contract_dir.resolve(strict=True)
    preflight_path = args.preflight.resolve(strict=True)
    legacy_property_dir = args.legacy_property_dir.resolve(strict=True)
    legacy_antifold_dir = args.legacy_antifold_dir.resolve(strict=True)
    legacy_plan_dir = args.legacy_property_plan_dir.resolve(strict=True)
    antifold_plan_dir = args.antifold_plan_dir.resolve(strict=True)
    candidates_path = contract_dir / "nb252_allowed_single_mutants.csv"
    legacy_property_path = legacy_property_dir / "unified_single_mutant_property_evidence.csv"
    legacy_antifold_path = legacy_antifold_dir / "unified_single_mutant_antifold_evidence.csv"
    legacy_tool_contract_path = legacy_plan_dir / "unified_property_scoring_contract.json"
    antifold_contract_path = antifold_plan_dir / "antifold_environment_contract.json"
    views_path = antifold_plan_dir / "antifold_structure_views.csv"
    sources = [
        candidates_path, preflight_path, legacy_property_path, legacy_antifold_path,
        legacy_tool_contract_path, antifold_contract_path, views_path,
    ]
    for path in sources:
        path.resolve(strict=True)
    current = _csv(candidates_path)
    first = current[0]
    sequence = str(first["sequence"])
    index = int(first["reported_sequence_index_1based"])
    parent = sequence[: index - 1] + str(first["wt_residue"]) + sequence[index:]
    audit, validation_samples, validation_expected, anti_targets, completion_samples, gate = build_reuse_plan(
        current, _json(preflight_path), parent, _csv(legacy_property_path), _csv(legacy_antifold_path)
    )
    if args.check_only:
        print(json.dumps(gate, sort_keys=True))
        return 0

    output_dir = args.output_dir.absolute()
    summary_path = args.run_summary.absolute()
    if output_dir.exists() or summary_path.exists():
        raise FileExistsError("Refusing to overwrite expression-property completion plan")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    names = {
        "audit": "reuse_mapping_audit.csv",
        "validation_samples": "reuse_validation_samples.csv",
        "validation_fasta": "reuse_validation_sequences.fasta",
        "validation_expected": "reuse_validation_expected.csv",
        "antifold_targets": "reuse_validation_antifold_targets.csv",
        "completion_samples": "missing_property_samples.csv",
        "completion_fasta": "missing_property_sequences.fasta",
        "contract": "expression_property_completion_contract.json",
        "gate": "expression_property_completion_plan_gate.json",
    }
    legacy_tool_contract = _json(legacy_tool_contract_path)
    contract = {
        "schema_version": 1,
        "contract_name": "nb252_expression_property_completion_v2",
        "status": "pass",
        "generated_at": generated_at,
        "authoritative_parent": {"sample_uid": "LTT__Nb252", "length": 128, "sequence": parent},
        "candidate_count": 847,
        "reuse_join_key": [
            "reported_sequence_index_1based", "wt_residue", "mutant_residue", "full_128aa_sequence"
        ],
        "legacy_candidate_id_used_as_join_key": False,
        "netsolp": legacy_tool_contract["netsolp"],
        "nanomelt": legacy_tool_contract["nanomelt"],
        "antifold": {
            "environment_contract": str(antifold_contract_path.relative_to(ROOT)),
            "structure_views": str(views_path.relative_to(ROOT)),
            "repeat_inference_scope": "three_frozen_WT_structure_views_once_not_per_mutant",
            "reuse_semantics": "position-specific WT/mutant log-probability lookup",
        },
        "repeat_validation": {
            "candidate_count": 12,
            "score_row_count_including_wt": 13,
            "netsolp_absolute_tolerance": 5e-8,
            "nanomelt_tm_absolute_tolerance_c": 0.0050001,
            "antifold_absolute_tolerance": 5e-8,
            "must_pass_before_completion": True,
        },
        "runtime_estimate": {
            "validation": "approximately 5-30 minutes",
            "missing_126_scores": "approximately 15-90 minutes, NetSolP-dominated",
            "likely_over_one_hour": True,
            "likely_over_five_hours": False,
            "slurm_required": True,
            "resume_required": False,
        },
        "candidate_selection_performed": False,
        "ranking_performed": False,
    }
    gate = {**gate, "generated_at": generated_at}
    summary = {
        "schema_version": 1,
        "status": "pass",
        "generated_at": generated_at,
        "elapsed_seconds": round(time.perf_counter() - started, 6),
        "python": platform.python_version(),
        "counts": {
            "candidates": len(audit),
            "reused_property": sum(row["netsolp_reuse_status"] != "requires_new_score" for row in audit),
            "new_property": sum(row["netsolp_reuse_status"] == "requires_new_score" for row in audit),
            "antifold_reused": len(audit),
            "validation_score_rows": len(validation_samples),
            "completion_score_rows": len(completion_samples),
        },
        "release": gate["release"],
        "outputs": {key: str(output_dir / value) for key, value in names.items()},
    }
    with tempfile.TemporaryDirectory(prefix=".expression-property-plan-", dir=ROOT) as temp_name:
        staging = Path(temp_name)
        staged = {key: staging / value for key, value in names.items()}
        _write_csv(staged["audit"], audit)
        _write_csv(staged["validation_samples"], validation_samples)
        _write_fasta(staged["validation_fasta"], validation_samples)
        _write_csv(staged["validation_expected"], validation_expected)
        _write_csv(staged["antifold_targets"], anti_targets)
        _write_csv(staged["completion_samples"], completion_samples)
        _write_fasta(staged["completion_fasta"], completion_samples)
        _write_json(staged["contract"], contract)
        _write_json(staged["gate"], gate)
        staged_summary = staging / "run_summary.json"
        _write_json(staged_summary, summary)
        pairs = {staged[key]: output_dir / value for key, value in names.items()}
        pairs[staged_summary] = summary_path
        validate_file_paths(project_root=ROOT, source_paths=sources, target_paths=pairs.values())
        for target in pairs.values():
            target.parent.mkdir(parents=True, exist_ok=True)
        replace_staged_files(pairs, project_root=ROOT, protected_source_paths=sources)
    return 0


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fields = list(rows[0])
    fields.extend(key for row in rows[1:] for key in row if key not in fields)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_fasta(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(f">{row['score_id']}\n{row['sequence_raw']}\n")


def _write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    raise SystemExit(main())
