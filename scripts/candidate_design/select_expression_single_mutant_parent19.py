#!/usr/bin/env python3
"""Select the approved 19 single-mutant parents from the frozen trial30 panel."""

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

from antibody_optimization.expression_parent_panel import build_parent19_panel  # noqa: E402
from antibody_optimization.file_transaction import (  # noqa: E402
    replace_staged_files,
    validate_file_paths,
)


SELECTED_MUTATIONS = (
    "F30S", "F30A", "F30R",
    "Q1M", "Q1H", "Q1D",
    "T27F", "T27D", "T27S",
    "Q5V", "A23S", "F29T", "Y32L", "E44G",
    "S55G", "K86R", "P87T", "V97S", "T99F",
)

SELECTED_REASONS = {
    "F30S": "retain_f30_multi_family_candidate_with_strong_netsolp_s_support",
    "F30A": "retain_f30_multi_family_candidate_without_moderate_adverse_band",
    "F30R": "retain_f30_multi_family_candidate_with_complementary_charge_change",
    "Q1M": "retain_q1_strong_antifold_compatibility_hypothesis",
    "Q1H": "retain_q1_moderate_netsolp_s_and_weak_tm_support_without_adverse_band",
    "Q1D": "retain_q1_moderate_netsolp_s_support_with_distinct_acidic_substitution",
    "T27F": "retain_t27_strong_antifold_hypothesis_with_explicit_weak_s_tradeoff",
    "T27D": "retain_t27_moderate_antifold_hypothesis_without_adverse_band",
    "T27S": "retain_t27_moderate_antifold_hypothesis_without_adverse_band",
    "Q5V": "retain_only_q5_candidate_and_natural_consensus_reversion",
    "A23S": "retain_only_a23_candidate_under_one_per_nonfocal_position_rule",
    "F29T": "retain_f29_option_without_the_f29s_weak_tm_adverse_band",
    "Y32L": "retain_only_y32_candidate_under_one_per_nonfocal_position_rule",
    "E44G": "retain_e44_option_without_the_e44a_weak_u_adverse_band",
    "S55G": "retain_only_s55_candidate_under_one_per_nonfocal_position_rule",
    "K86R": "retain_only_k86_candidate_under_one_per_nonfocal_position_rule",
    "P87T": "retain_only_p87_candidate_under_one_per_nonfocal_position_rule",
    "V97S": "retain_only_v97_candidate_under_one_per_nonfocal_position_rule",
    "T99F": "retain_user_selected_stable_word_hypothesis_as_explicit_exploratory_parent",
}

NAMES = {
    "panel": "expression_single_mutant_parent19.csv",
    "audit": "expression_single_mutant_parent19_audit.csv",
    "fasta": "expression_single_mutant_parent19.fasta",
    "contract": "expression_single_mutant_parent19_contract.json",
    "gate": "expression_single_mutant_parent19_gate.json",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trial30", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    trial30 = args.trial30.expanduser().resolve(strict=True)
    output_dir = args.output_dir.expanduser().absolute()
    run_summary = args.run_summary.expanduser().absolute()
    output_dir.mkdir(parents=True, exist_ok=True)
    run_summary.parent.mkdir(parents=True, exist_ok=True)
    validated = validate_file_paths(
        project_root=PROJECT_ROOT,
        source_paths=[trial30],
        target_paths=[*[output_dir / name for name in NAMES.values()], run_summary],
    )
    finals = dict(zip(NAMES, validated.target_paths[:-1], strict=True))
    run_summary = validated.target_paths[-1]
    existing = [path for path in [*finals.values(), run_summary] if path.exists()]
    if existing:
        raise FileExistsError("Refusing to overwrite existing outputs:\n" + "\n".join(map(str, existing)))

    result = build_parent19_panel(_csv(trial30), SELECTED_MUTATIONS, SELECTED_REASONS)
    facts = result["facts"]
    contract = {
        "schema_version": 1,
        "generated_at": generated_at,
        "optimization_target": "BL21 expression yield",
        "source_panel_role": "frozen_existing_30_single_mutant_computational_trial_panel",
        "source_panel_modified": False,
        "current_final_panel_composition": {
            "single_mutants": 19,
            "double_mutants": 11,
            "total_candidates": 30,
        },
        "single_parent_selection_rule": (
            "retain three substitutions at reported F30, Q1, and T27; retain one "
            "substitution at each of the other ten positions represented in trial30"
        ),
        "selection_evidence_policy": (
            "reuse frozen categorical magnitude bands and risk annotations; do not rerun "
            "predictors and do not use within-band decimal differences as a global ranking"
        ),
        "selected_mutations_in_order": list(SELECTED_MUTATIONS),
        "selected_reasons": SELECTED_REASONS,
        "next_stage": (
            "enumerate all distinct-position pairs from these 19 parents; same-position "
            "alternative substitutions are mutually exclusive"
        ),
        "active_tools_for_future_double_scoring": ["NetSolP", "NanoMelt", "AntiFold"],
        "rosetta_used": False,
        "double_mutant_enumeration_performed": False,
        "final_30_candidate_panel_released": False,
    }
    gate = {
        "schema_version": 1,
        "generated_at": generated_at,
        "status": "pass",
        "release": "parent_19_single_mutants_ready_for_double_enumeration",
        **facts,
        "interpretation": (
            "The 19 rows are computationally selected single-mutant parents. Existing trial30 "
            "artifacts remain unchanged; the 11 double mutants have not yet been generated."
        ),
    }

    with tempfile.TemporaryDirectory(prefix=".expression-parent19-", dir=PROJECT_ROOT) as temp:
        stage = Path(temp)
        staged = {key: stage / name for key, name in NAMES.items()}
        staged_run = stage / "run_summary.json"
        _write_csv(staged["panel"], result["panel_rows"])
        _write_csv(staged["audit"], result["audit_rows"])
        _write_fasta(staged["fasta"], result["panel_rows"])
        _write_json(staged["contract"], contract)
        _write_json(staged["gate"], gate)
        _write_json(
            staged_run,
            {
                "schema_version": 1,
                "status": "pass",
                "generated_at": generated_at,
                "elapsed_seconds": round(time.perf_counter() - started, 6),
                "python": platform.python_version(),
                "command_argv": [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
                **facts,
                "outputs": {key: str(path) for key, path in finals.items()},
            },
        )
        replace_staged_files(
            {**{staged[key]: finals[key] for key in NAMES}, staged_run: run_summary},
            project_root=PROJECT_ROOT,
            protected_source_paths=validated.source_paths,
        )
    return 0


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _fields(rows: list[dict[str, object]]) -> list[str]:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    return fields


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=_fields(rows), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_fasta(path: Path, rows: list[dict[str, object]]) -> None:
    lines: list[str] = []
    for row in rows:
        code = str(row["mutation_reported_label"]).replace("Nb252 reported_seq ", "")
        lines.extend([f">{row['candidate_id']} mutation={code} parent19_order={row['parent19_selection_order']}", str(row["sequence"])])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def _write_json(path: Path, value: object) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )


if __name__ == "__main__":
    raise SystemExit(main())
