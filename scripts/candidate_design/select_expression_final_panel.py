#!/usr/bin/env python3
"""Select 11 expression-oriented doubles and release the final 19+11 panel."""

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

from antibody_optimization.expression_final_panel import select_final_expression_panel  # noqa: E402
from antibody_optimization.expression_final_panel_plot import plot_final_expression_panel  # noqa: E402
from antibody_optimization.file_transaction import replace_staged_files, validate_file_paths  # noqa: E402


NAMES = {
    "audit": "expression_double_mutant_final_selection_audit.csv",
    "selected": "expression_double_mutant_selected11.csv",
    "reserves": "expression_double_mutant_reserves.csv",
    "final": "nb252_final_30_candidate_panel.csv",
    "fasta": "nb252_final_30_candidate_panel.fasta",
    "contract": "expression_final_panel_selection_contract.json",
    "gate": "expression_final_panel_gate.json",
    "plot_png": "expression_final_panel_overview.png",
    "plot_svg": "expression_final_panel_overview.svg",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--double-matrix", type=Path, required=True)
    parser.add_argument("--double-gate", type=Path, required=True)
    parser.add_argument("--parent19", type=Path, required=True)
    parser.add_argument("--parent19-contract", type=Path, required=True)
    parser.add_argument("--constraints", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    sources = [
        args.double_matrix,
        args.double_gate,
        args.parent19,
        args.parent19_contract,
        args.constraints,
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

    result = select_final_expression_panel(
        _csv(validated.source_paths[0]),
        _json(validated.source_paths[1]),
        _csv(validated.source_paths[2]),
        _json(validated.source_paths[3]),
        _json(validated.source_paths[4]),
    )
    facts = result["facts"]
    contract = {
        "schema_version": 1,
        "generated_at": generated_at,
        "optimization_target": "BL21 expression yield",
        "source_double_candidate_count": 162,
        "eligible_double_definition": {
            "hard_sequence_risk_count": 0,
            "moderate_or_strong_adverse_property_bands_allowed": False,
            "required_moderate_or_strong_favorable_predictor_family_count_minimum": 1,
            "stable_word_gain_can_rescue_adverse_property_band": False,
        },
        "predictor_families": {
            "NetSolP": "maximum categorical grade of U and S",
            "NanoMelt": "predicted apparent Tm categorical grade",
            "AntiFold": "worse categorical grade of the two constituent mutations; no double-structure rerun",
        },
        "evidence_layers": {
            "A_three_families": 3,
            "B_two_families": 2,
            "C_one_family": 1,
        },
        "diversity_constraints": {
            "maximum_selected_per_reported_position_pair": 1,
            "maximum_uses_per_exact_component_mutation": 2,
            "maximum_uses_per_reported_position": 3,
        },
        "optimizer": result["optimizer"],
        "active_predictors_rerun": False,
        "rosetta_used": False,
        "af3_batch_prediction_used": False,
        "interpretation": (
            "The selected sequences are computationally prioritized expression candidates. "
            "Categorical magnitude bands, not within-band decimals, drive selection."
        ),
    }
    gate = {
        "schema_version": 1,
        "generated_at": generated_at,
        "status": "pass",
        "release": "final_19_single_plus_11_double_panel_released_for_experimental_testing",
        **facts,
        "hard_constraint_validation": "pass",
        "diversity_constraint_validation": "pass",
        "predictor_rerun_performed": False,
        "interpretation": (
            "The panel contains 19 frozen approved singles and 11 doubles selected from the "
            "complete scored 162-row space. It is recommended for experimental testing, not "
            "reported as experimentally validated."
        ),
    }

    with tempfile.TemporaryDirectory(prefix=".expression-final-panel-", dir=PROJECT_ROOT) as temp:
        stage = Path(temp)
        staged = {key: stage / name for key, name in NAMES.items()}
        staged_run = stage / "run_summary.json"
        _write_csv(staged["audit"], result["audit_rows"])
        _write_csv(staged["selected"], result["selected_double_rows"])
        _write_csv(staged["reserves"], result["reserve_rows"])
        _write_csv(staged["final"], result["final_rows"])
        _write_fasta(staged["fasta"], result["final_rows"])
        _write_json(staged["contract"], contract)
        _write_json(staged["gate"], gate)
        plot_final_expression_panel(
            result["audit_rows"],
            result["selected_double_rows"],
            staged["plot_png"],
            staged["plot_svg"],
        )
        _write_json(
            staged_run,
            {
                "schema_version": 1,
                "status": "pass",
                "generated_at": generated_at,
                "elapsed_seconds": round(time.perf_counter() - started, 6),
                "python": platform.python_version(),
                "scipy_solver": result["optimizer"]["solver"],
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
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def _fields(rows: list[dict[str, object]]) -> list[str]:
    fields = []
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
    lines = []
    for row in rows:
        lines.extend(
            [
                f">{row['candidate_id']} order={row['final_panel_order']} kind={row['candidate_kind']} mutations={row['mutation_set']}",
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
