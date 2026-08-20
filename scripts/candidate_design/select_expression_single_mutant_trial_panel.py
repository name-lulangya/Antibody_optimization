#!/usr/bin/env python3
"""Select a magnitude-banded 30-member Nb252 expression single-mutant trial panel."""

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

from antibody_optimization.expression_panel_selection import (  # noqa: E402
    MAGNITUDE_THRESHOLDS,
    build_expression_trial_panel,
)
from antibody_optimization.expression_panel_selection_plot import (  # noqa: E402
    render_expression_trial_panel,
)
from antibody_optimization.file_transaction import (  # noqa: E402
    replace_staged_files,
    validate_file_paths,
)


NAMES = {
    "audit": "expression_single_mutant_selection_audit.csv",
    "shortlist": "expression_single_mutant_magnitude_shortlist.csv",
    "panel": "expression_single_mutant_trial30.csv",
    "reserve": "expression_single_mutant_reserve.csv",
    "panel_fasta": "expression_single_mutant_trial30.fasta",
    "reserve_fasta": "expression_single_mutant_reserve.fasta",
    "summary": "expression_single_mutant_selection_summary.csv",
    "contract": "expression_single_mutant_selection_contract.json",
    "gate": "expression_single_mutant_selection_gate.json",
    "png": "expression_single_mutant_selection.png",
    "svg": "expression_single_mutant_selection.svg",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--upstream-gate", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    matrix = args.matrix.expanduser().resolve(strict=True)
    upstream_gate_path = args.upstream_gate.expanduser().resolve(strict=True)
    output_dir = args.output_dir.expanduser().absolute()
    run_summary = args.run_summary.expanduser().absolute()
    output_dir.mkdir(parents=True, exist_ok=True)
    run_summary.parent.mkdir(parents=True, exist_ok=True)
    validated = validate_file_paths(
        project_root=PROJECT_ROOT,
        source_paths=[matrix, upstream_gate_path],
        target_paths=[*[output_dir / name for name in NAMES.values()], run_summary],
    )
    finals = dict(zip(NAMES, validated.target_paths[:-1], strict=True))
    run_summary = validated.target_paths[-1]
    existing = [path for path in [*finals.values(), run_summary] if path.exists()]
    if existing:
        raise FileExistsError("Refusing to overwrite existing outputs:\n" + "\n".join(map(str, existing)))

    rows = _csv(matrix)
    upstream_gate = _json(upstream_gate_path)
    result = build_expression_trial_panel(rows, upstream_gate)
    facts = result["facts"]
    contract = {
        "schema_version": 1,
        "generated_at": generated_at,
        "optimization_target": "BL21 expression yield",
        "candidate_kind": "single_substitution_only",
        "metric_direction": "higher_is_better_for_all_four_displayed_changes",
        "magnitude_thresholds": {
            key: {
                "negligible_to_weak_absolute_boundary": values[0],
                "weak_to_moderate_absolute_boundary": values[1],
                "moderate_to_strong_absolute_boundary": values[2],
            }
            for key, values in MAGNITUDE_THRESHOLDS.items()
        },
        "weak_and_negligible_selection_grade": 0,
        "netsolp_u_and_s_independent_family_count": 1,
        "magnitude_shortlist_rule": (
            "no strong-adverse metric and at least one moderate-or-strong favorable "
            "predictor family"
        ),
        "strict_core_rule": "magnitude shortlist, no hard sequence risk, no moderate-adverse metric",
        "controlled_tradeoff_rule": (
            "magnitude shortlist, no hard sequence risk, exactly one moderate-adverse metric; "
            "used only as an explicitly labelled diversity layer"
        ),
        "hard_sequence_risks": [
            "new_proline_backbone_constraint",
            "new_dense_local_hydrophobic_window",
            "new_extreme_local_charge_cluster",
        ],
        "selection_order": (
            "retain every strict-core candidate; fill remaining slots from the controlled-"
            "tradeoff layer by new-position coverage, categorical tier, strong favorable "
            "family count, favorable family count, AntiFold provenance, soft risk count, "
            "stable-word tie-break, candidate identifier; then apply the user-reviewed "
            "T99N-to-T99F exploratory substitution"
        ),
        "raw_within_band_values_used_for_ranking": False,
        "stable_word_role": (
            "last-order tie-break generally; T99F is one explicit hypothesis-testing "
            "exception and is not presented as predictor-supported optimization"
        ),
        "stable_word_exploratory_candidate": "Nb252 reported_seq T99F",
        "stable_word_exploratory_replaces": "Nb252 reported_seq T99N",
        "antifold_provenance_policy": (
            "experimental complex preferred; AF3 VHH-only fallback remains explicitly labelled"
        ),
        "final_experimental_panel_released": False,
    }
    gate = {
        "schema_version": 1,
        "status": "pass",
        "release": "trial_30_single_mutant_panel_ready_for_user_review",
        "generated_at": generated_at,
        **facts,
        "interpretation": (
            "Computational trial selection only; predictions are not measured expression, "
            "solubility, stability, or experimental structural validation."
        ),
    }
    summary_rows = [
        {"stage": "all_constrained_single_mutants", "count": facts["candidate_count"]},
        {"stage": "magnitude_shortlist", "count": facts["magnitude_shortlist_count"]},
        {"stage": "strict_core", "count": facts["strict_core_count"]},
        {"stage": "controlled_tradeoff", "count": facts["controlled_tradeoff_count"]},
        {
            "stage": "stable_word_exploratory_exception",
            "count": facts["trial_panel_stable_word_exploratory_count"],
        },
        {"stage": "blocked_sequence_risk", "count": facts["blocked_sequence_risk_count"]},
        {
            "stage": "blocked_multiple_moderate_adverse",
            "count": facts["blocked_multiple_moderate_adverse_count"],
        },
        {"stage": "trial_panel", "count": facts["trial_panel_count"]},
        {"stage": "reserve", "count": facts["reserve_count"]},
    ]

    with tempfile.TemporaryDirectory(prefix=".expression-panel-selection-", dir=PROJECT_ROOT) as temp:
        stage = Path(temp)
        staged = {key: stage / name for key, name in NAMES.items()}
        staged_run = stage / "run_summary.json"
        _write_csv(staged["audit"], result["audit_rows"])
        _write_csv(staged["shortlist"], result["shortlist_rows"])
        _write_csv(staged["panel"], result["panel_rows"])
        _write_csv(staged["reserve"], result["reserve_rows"])
        _write_fasta(staged["panel_fasta"], result["panel_rows"])
        _write_fasta(staged["reserve_fasta"], result["reserve_rows"])
        _write_csv(staged["summary"], summary_rows)
        _write_json(staged["contract"], contract)
        _write_json(staged["gate"], gate)
        render_expression_trial_panel(
            result["audit_rows"],
            result["panel_rows"],
            result["reserve_rows"],
            facts,
            staged["png"],
            staged["svg"],
        )
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


def _json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def _fields(rows: list[dict[str, object]]) -> list[str]:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    return fields


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=_fields(rows), lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def _write_fasta(path: Path, rows: list[dict[str, object]]) -> None:
    lines: list[str] = []
    for row in rows:
        lines.extend(
            [
                f">{row['candidate_id']} {row['mutation_reported_label']} "
                f"tier={row['selection_tier']} status={row['trial_selection_status']}",
                str(row["sequence"]),
            ]
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def _write_json(path: Path, value: object) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )


if __name__ == "__main__":
    raise SystemExit(main())
