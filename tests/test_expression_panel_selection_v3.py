import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from antibody_optimization.expression_panel_selection_v3 import (  # noqa: E402
    build_expression_single_mutant_panel_v3,
)


SOURCE = ROOT / "docs/result_artifacts/candidate_design/stable_word_single_mutant_v1_20260819"
ANTIFOLD = ROOT / "docs/result_artifacts/candidate_design/unified_single_mutant_antifold_20260815"
SCRIPT = ROOT / "scripts/candidate_design/select_expression_single_mutant_panel_v3.py"


def _csv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _result():
    return build_expression_single_mutant_panel_v3(
        _csv(SOURCE / "expression_single_mutant_property_stable_word_matrix.csv"),
        json.loads((SOURCE / "stable_word_single_mutant_gate.json").read_text(encoding="utf-8")),
        _csv(ANTIFOLD / "unified_single_mutant_antifold_evidence.csv"),
    )


def test_real_v3_selection_uses_three_metrics_and_combined_antifold_veto():
    result = _result()
    facts = result["facts"]
    assert facts["candidate_count"] == 847
    assert facts["antifold_veto_count"] == 151
    assert facts["antifold_veto_experimental_complex_count"] == 130
    assert facts["antifold_veto_af3_fallback_count"] == 21
    assert facts["qualified_count"] == 61
    assert facts["qualified_tier_counts"] == {
        "A_multi_metric": 8,
        "B_single_metric_strong": 11,
        "C_single_metric_moderate": 42,
    }
    assert facts["selected_count"] == 30
    assert facts["reserve_count"] == 31
    assert facts["selected_unique_position_count"] == 23
    assert facts["selected_maximum_per_position"] == 2
    assert facts["antifold_positive_credit_used"] is False
    assert facts["raw_within_band_values_used_for_ranking"] is False
    panel = result["panel_rows"]
    assert all(row["antifold_veto_status"] == "pass" for row in panel)
    assert all(int(row["positive_metric_count_v3"]) >= 1 for row in panel)
    assert all(int(row["hard_sequence_risk_count_v3"]) == 0 for row in panel)
    assert all(int(row["strong_adverse_property_count_v3"]) == 0 for row in panel)
    assert len({row["sequence"] for row in panel}) == 30


def test_antifold_rank_is_complete_and_veto_requires_both_conditions():
    audit = _result()["audit_rows"]
    assert {int(row["antifold_position_state_count"]) for row in audit} == {20}
    for row in audit:
        expected = (
            float(row["antifold_selection_delta_log_probability"]) <= -3.0
            and int(row["antifold_mutant_rank_worst_first"]) <= 4
        )
        assert (row["antifold_veto_status"] == "veto") is expected
    assert any(
        float(row["antifold_selection_delta_log_probability"]) <= -3.0
        and int(row["antifold_mutant_rank_worst_first"]) > 4
        and row["antifold_veto_status"] == "pass"
        for row in audit
    )


def test_v3_cli_writes_final30_without_changing_sources():
    matrix = SOURCE / "expression_single_mutant_property_stable_word_matrix.csv"
    full_antifold = ANTIFOLD / "unified_single_mutant_antifold_evidence.csv"
    matrix_bytes = matrix.read_bytes()
    antifold_bytes = full_antifold.read_bytes()
    with tempfile.TemporaryDirectory(prefix=".test-expression-v3-", dir=ROOT) as temp:
        base = Path(temp)
        output = base / "selection"
        summary = base / "run.json"
        subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--matrix",
                str(matrix),
                "--upstream-gate",
                str(SOURCE / "stable_word_single_mutant_gate.json"),
                "--full-antifold-evidence",
                str(full_antifold),
                "--output-dir",
                str(output),
                "--run-summary",
                str(summary),
                "--generated-at",
                "2026-08-25T12:00:00+08:00",
            ],
            check=True,
            cwd=ROOT,
        )
        gate = json.loads((output / "expression_single_mutant_v3_gate.json").read_text(encoding="utf-8"))
        contract = json.loads((output / "expression_single_mutant_v3_contract.json").read_text(encoding="utf-8"))
        assert gate["status"] == "pass"
        assert gate["selected_count"] == 30
        assert gate["release"] == "v3_final_30_single_mutants_ready_for_experimental_testing"
        assert contract["netsolp_u_and_s_counted_separately"] is True
        assert contract["antifold_role"] == "negative_veto_only_no_positive_selection_credit"
        assert len(_csv(output / "expression_single_mutant_v3_final30.csv")) == 30
        assert len(_csv(output / "expression_single_mutant_v3_reserve.csv")) == 31
        assert (output / "expression_single_mutant_v3_selection.png").stat().st_size > 1000
        assert (output / "expression_single_mutant_v3_selection.svg").stat().st_size > 1000
    assert matrix.read_bytes() == matrix_bytes
    assert full_antifold.read_bytes() == antifold_bytes
