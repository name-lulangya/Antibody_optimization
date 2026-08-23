from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from antibody_optimization.expression_final_panel import (  # noqa: E402
    ExpressionFinalPanelError,
    select_final_expression_panel,
)


MATRIX_DIR = ROOT / "docs/result_artifacts/candidate_design/expression_double_mutant_property_matrix_20260822"
PARENT_DIR = ROOT / "docs/result_artifacts/candidate_design/expression_single_mutant_parent19_20260822"
CONSTRAINTS = ROOT / "docs/result_artifacts/input_baseline/vhh_conservation_consensus_v2_20260819/nb252_expression_design_constraints.json"
SCRIPT = ROOT / "scripts/candidate_design/select_expression_final_panel.py"


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _inputs() -> tuple[object, ...]:
    return (
        _csv(MATRIX_DIR / "expression_double_mutant_property_matrix.csv"),
        _json(MATRIX_DIR / "expression_double_mutant_property_matrix_gate.json"),
        _csv(PARENT_DIR / "expression_single_mutant_parent19.csv"),
        _json(PARENT_DIR / "expression_single_mutant_parent19_contract.json"),
        _json(CONSTRAINTS),
    )


def test_real_final_panel_is_exact_and_diverse() -> None:
    result = select_final_expression_panel(*_inputs())
    facts = result["facts"]
    assert facts["double_candidate_count"] == 162
    assert facts["eligible_double_count"] == 84
    assert facts["selected_double_count"] == 11
    assert facts["selected_evidence_layer_counts"] == {
        "A_three_families": 3,
        "B_two_families": 4,
        "C_one_family": 4,
    }
    assert facts["selected_unique_component_mutation_count"] == 13
    assert facts["selected_unique_reported_position_count"] == 10
    assert facts["final_candidate_count"] == 30
    assert Counter(row["candidate_kind"] for row in result["final_rows"]) == {
        "single_mutant": 19,
        "double_mutant": 11,
    }
    selected = result["selected_double_rows"]
    assert [row["mutation_set"] for row in selected] == [
        "F30R;Q5V",
        "F30S;S55G",
        "Q5V;E44G",
        "F30R;Q1H",
        "E44G;V97S",
        "S55G;P87T",
        "T27D;V97S",
        "Q1D;T27D",
        "Q1H;F29T",
        "F29T;K86R",
        "T27S;K86R",
    ]
    component_counts = Counter(
        row[key] for row in selected for key in ("mutation_a", "mutation_b")
    )
    position_counts = Counter(
        int(row[key])
        for row in selected
        for key in ("position_a_reported_1based", "position_b_reported_1based")
    )
    assert max(component_counts.values()) <= 2
    assert max(position_counts.values()) <= 3
    assert all(int(row["hard_sequence_risk_count"]) == 0 for row in selected)
    assert all(row["double_selection_eligibility"] == "eligible" for row in selected)


def test_selection_is_independent_of_input_row_order() -> None:
    inputs = list(_inputs())
    forward = select_final_expression_panel(*inputs)
    inputs[0] = list(reversed(inputs[0]))
    reverse = select_final_expression_panel(*inputs)
    assert [row["candidate_id"] for row in forward["selected_double_rows"]] == [
        row["candidate_id"] for row in reverse["selected_double_rows"]
    ]
    assert forward["optimizer"] == reverse["optimizer"]


def test_closed_matrix_gate_blocks_selection() -> None:
    inputs = list(_inputs())
    inputs[1] = {**inputs[1], "status": "blocked"}
    with pytest.raises(ExpressionFinalPanelError, match="gate is not open"):
        select_final_expression_panel(*inputs)


def test_cli_writes_complete_outputs_without_touching_sources() -> None:
    matrix = MATRIX_DIR / "expression_double_mutant_property_matrix.csv"
    parent = PARENT_DIR / "expression_single_mutant_parent19.csv"
    before = (matrix.read_bytes(), parent.read_bytes())
    with tempfile.TemporaryDirectory(prefix=".test-expression-final-", dir=ROOT) as temp:
        base = Path(temp)
        output = base / "result"
        run_summary = base / "run_summary.json"
        subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--double-matrix",
                str(matrix),
                "--double-gate",
                str(MATRIX_DIR / "expression_double_mutant_property_matrix_gate.json"),
                "--parent19",
                str(parent),
                "--parent19-contract",
                str(PARENT_DIR / "expression_single_mutant_parent19_contract.json"),
                "--constraints",
                str(CONSTRAINTS),
                "--output-dir",
                str(output),
                "--run-summary",
                str(run_summary),
                "--generated-at",
                "2026-08-22T20:30:00+08:00",
            ],
            cwd=ROOT,
            check=True,
        )
        gate = _json(output / "expression_final_panel_gate.json")
        assert gate["status"] == "pass"
        assert gate["final_candidate_count"] == 30
        assert len(_csv(output / "expression_double_mutant_final_selection_audit.csv")) == 162
        assert len(_csv(output / "expression_double_mutant_selected11.csv")) == 11
        assert len(_csv(output / "expression_double_mutant_reserves.csv")) == 73
        assert len(_csv(output / "nb252_final_30_candidate_panel.csv")) == 30
        assert (output / "expression_final_panel_overview.png").stat().st_size > 10_000
        assert (output / "expression_final_panel_overview.svg").stat().st_size > 10_000
        assert _json(run_summary)["status"] == "pass"
    assert (matrix.read_bytes(), parent.read_bytes()) == before
