#!/usr/bin/env python3
"""Select the V3 30-member Nb252 expression single-mutant panel locally."""

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

from antibody_optimization.expression_panel_selection import MAGNITUDE_THRESHOLDS  # noqa: E402
from antibody_optimization.expression_panel_selection_v3 import (  # noqa: E402
    ANTIFOLD_BOTTOM_COUNT,
    ANTIFOLD_DELTA_VETO,
    build_expression_single_mutant_panel_v3,
)
from antibody_optimization.expression_panel_selection_v3_plot import (  # noqa: E402
    render_expression_single_mutant_panel_v3,
)
from antibody_optimization.file_transaction import (  # noqa: E402
    replace_staged_files,
    validate_file_paths,
)


NAMES = {
    "audit": "expression_single_mutant_v3_audit.csv",
    "qualified": "expression_single_mutant_v3_qualified.csv",
    "panel": "expression_single_mutant_v3_final30.csv",
    "reserve": "expression_single_mutant_v3_reserve.csv",
    "panel_fasta": "expression_single_mutant_v3_final30.fasta",
    "reserve_fasta": "expression_single_mutant_v3_reserve.fasta",
    "summary": "expression_single_mutant_v3_summary.csv",
    "contract": "expression_single_mutant_v3_contract.json",
    "gate": "expression_single_mutant_v3_gate.json",
    "png": "expression_single_mutant_v3_selection.png",
    "svg": "expression_single_mutant_v3_selection.svg",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--upstream-gate", type=Path, required=True)
    parser.add_argument("--full-antifold-evidence", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    sources = [
        args.matrix.expanduser().resolve(strict=True),
        args.upstream_gate.expanduser().resolve(strict=True),
        args.full_antifold_evidence.expanduser().resolve(strict=True),
    ]
    output_dir = args.output_dir.expanduser().absolute()
    run_summary = args.run_summary.expanduser().absolute()
    output_dir.mkdir(parents=True, exist_ok=True)
    run_summary.parent.mkdir(parents=True, exist_ok=True)
    validated = validate_file_paths(
        project_root=PROJECT_ROOT,
        source_paths=sources,
        target_paths=[*[output_dir / name for name in NAMES.values()], run_summary],
    )
    finals = dict(zip(NAMES, validated.target_paths[:-1], strict=True))
    run_summary = validated.target_paths[-1]
    existing = [path for path in [*finals.values(), run_summary] if path.exists()]
    if existing:
        raise FileExistsError("Refusing to overwrite existing outputs:\n" + "\n".join(map(str, existing)))

    result = build_expression_single_mutant_panel_v3(
        _csv(validated.source_paths[0]),
        _json(validated.source_paths[1]),
        _csv(validated.source_paths[2]),
    )
    facts = result["facts"]
    contract = {
        "schema_version": 3,
        "generated_at": generated_at,
        "optimization_target": "BL21 expression yield",
        "deliverable": "30 single substitutions only",
        "positive_selection_metrics": ["NetSolP U", "NetSolP S", "NanoMelt predicted Tm"],
        "netsolp_u_and_s_counted_separately": True,
        "magnitude_thresholds": {
            key: {
                "negligible_to_weak_absolute_boundary": values[0],
                "weak_to_moderate_absolute_boundary": values[1],
                "moderate_to_strong_absolute_boundary": values[2],
            }
            for key, values in MAGNITUDE_THRESHOLDS.items()
            if key != "antifold_log_probability"
        },
        "antifold_role": "negative_veto_only_no_positive_selection_credit",
        "antifold_veto_rule": {
            "delta_log_probability_maximum": ANTIFOLD_DELTA_VETO,
            "within_position_rank_direction": "1_is_worst",
            "within_position_rank_maximum": ANTIFOLD_BOTTOM_COUNT,
            "amino_acid_state_count": 20,
            "rule": "delta_log_probability <= -3 and mutant rank among worst four of 20",
        },
        "qualification_layers": {
            "A_multi_metric": "at least two moderate/strong favorable metrics; no moderate/strong adverse metric",
            "B_single_metric_strong": "one strong favorable metric; other metrics weak/negligible",
            "C_single_metric_moderate": "one moderate favorable metric; other metrics weak/negligible",
            "D_controlled_tradeoff": "at least two moderate/strong favorable metrics and exactly one moderate adverse metric",
        },
        "selection_order": (
            "maximum three candidates per position; cover first, second, then third candidate per position; "
            "within each round use A>B>C>D, strong-positive count, positive count, fewer moderate adverse metrics, "
            "fewer soft sequence risks, stable-word gain, then candidate identifier"
        ),
        "raw_within_band_values_used_for_ranking": False,
        "stable_word_role": "last-order categorical tie-break only",
        "historical_19_plus_11_panel_used_for_ranking": False,
    }
    gate = {
        "schema_version": 3,
        "status": "pass",
        "release": "v3_final_30_single_mutants_ready_for_experimental_testing",
        "generated_at": generated_at,
        **facts,
        "interpretation": "Computationally prioritized single mutants; no expression improvement has yet been experimentally demonstrated.",
    }
    summary_rows = [
        {"stage": "all_constrained_single_mutants", "count": facts["candidate_count"]},
        {"stage": "antifold_combined_veto", "count": facts["antifold_veto_count"]},
        {"stage": "v3_qualified", "count": facts["qualified_count"]},
        *[
            {"stage": f"qualified_{tier}", "count": count}
            for tier, count in sorted(facts["qualified_tier_counts"].items())
        ],
        {"stage": "selected_final30_single", "count": facts["selected_count"]},
        {"stage": "qualified_reserve", "count": facts["reserve_count"]},
    ]

    with tempfile.TemporaryDirectory(prefix=".expression-panel-v3-", dir=PROJECT_ROOT) as temp:
        stage = Path(temp)
        staged = {key: stage / name for key, name in NAMES.items()}
        staged_run = stage / "run_summary.json"
        _write_csv(staged["audit"], result["audit_rows"])
        _write_csv(staged["qualified"], result["qualified_rows"])
        _write_csv(staged["panel"], result["panel_rows"])
        _write_csv(staged["reserve"], result["reserve_rows"])
        _write_fasta(staged["panel_fasta"], result["panel_rows"])
        _write_fasta(staged["reserve_fasta"], result["reserve_rows"])
        _write_csv(staged["summary"], summary_rows)
        _write_json(staged["contract"], contract)
        _write_json(staged["gate"], gate)
        render_expression_single_mutant_panel_v3(
            result["audit_rows"], result["panel_rows"], result["reserve_rows"], facts, staged["png"], staged["svg"]
        )
        _write_json(
            staged_run,
            {
                "schema_version": 3,
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
        mutation = str(row["mutation_reported_label"]).replace("Nb252 reported_seq ", "")
        lines.extend([f">{row['candidate_id']} mutation={mutation} tier={row['selection_tier_v3']}", str(row["sequence"])])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def _write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    raise SystemExit(main())
