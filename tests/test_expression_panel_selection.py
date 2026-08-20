import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from antibody_optimization.expression_panel_selection import (
    build_expression_trial_panel,
    classify_change,
)

SOURCE = ROOT / "docs/result_artifacts/candidate_design/stable_word_single_mutant_v1_20260819"


def _csv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def test_magnitude_boundaries_do_not_promote_micro_changes():
    thresholds = (0.005, 0.010, 0.015)
    assert classify_change(0.004999, thresholds) == "negligible"
    assert classify_change(0.005, thresholds) == "weak_favorable"
    assert classify_change(-0.009999, thresholds) == "weak_adverse"
    assert classify_change(0.010, thresholds) == "moderate_favorable"
    assert classify_change(-0.015, thresholds) == "strong_adverse"


def test_real_847_space_yields_traceable_trial_30():
    rows = _csv(SOURCE / "expression_single_mutant_property_stable_word_matrix.csv")
    gate = json.loads((SOURCE / "stable_word_single_mutant_gate.json").read_text(encoding="utf-8"))
    result = build_expression_trial_panel(rows, gate)
    facts = result["facts"]
    assert facts["candidate_count"] == 847
    assert facts["magnitude_shortlist_count"] == 40
    assert facts["magnitude_shortlist_stable_word_gain_count"] == 0
    assert facts["strict_core_count"] == 26
    assert facts["controlled_tradeoff_count"] == 10
    assert facts["blocked_sequence_risk_count"] == 2
    assert facts["blocked_multiple_moderate_adverse_count"] == 2
    assert facts["trial_panel_count"] == 30
    assert facts["reserve_count"] == 7
    assert facts["trial_panel_unique_position_count"] >= 12
    assert facts["trial_panel_tier_counts"]["D_controlled_tradeoff"] == 3
    assert facts["trial_panel_tier_counts"]["E_stable_word_exploratory"] == 1
    assert facts["trial_panel_stable_word_gain_count"] == 1
    panel = result["panel_rows"]
    assert len({row["sequence"] for row in panel}) == 30
    assert sum(int(row["favorable_family_count"]) < 1 for row in panel) == 1
    assert all(int(row["strong_adverse_metric_count"]) == 0 for row in panel)
    assert all(int(row["hard_sequence_risk_count"]) == 0 for row in panel)
    assert sum(row["selection_eligibility_class"] == "strict_core" for row in panel) == 26
    by_mutation = {
        row["mutation_reported_label"].replace("Nb252 reported_seq ", ""): row
        for row in result["audit_rows"]
    }
    assert by_mutation["A23P"]["selection_eligibility_class"] == "blocked_sequence_risk"
    assert by_mutation["F30P"]["selection_eligibility_class"] == "blocked_sequence_risk"
    panel_mutations = {
        row["mutation_reported_label"].replace("Nb252 reported_seq ", "") for row in panel
    }
    assert "T99F" in panel_mutations
    assert "T99N" not in panel_mutations
    assert by_mutation["T99F"]["selection_tier"] == "E_stable_word_exploratory"


def test_cli_writes_trial_panel_and_visualization():
    with tempfile.TemporaryDirectory(prefix=".test-expression-panel-", dir=ROOT) as temp:
        base = Path(temp)
        output = base / "selection"
        summary = base / "run.json"
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "scripts/candidate_design/select_expression_single_mutant_trial_panel.py"),
                "--matrix",
                str(SOURCE / "expression_single_mutant_property_stable_word_matrix.csv"),
                "--upstream-gate",
                str(SOURCE / "stable_word_single_mutant_gate.json"),
                "--output-dir",
                str(output),
                "--run-summary",
                str(summary),
                "--generated-at",
                "2026-08-20T03:00:00+08:00",
            ],
            check=True,
            cwd=ROOT,
        )
        gate = json.loads((output / "expression_single_mutant_selection_gate.json").read_text(encoding="utf-8"))
        assert gate["status"] == "pass"
        assert gate["release"] == "trial_30_single_mutant_panel_ready_for_user_review"
        assert gate["final_experimental_panel_released"] is False
        assert len(_csv(output / "expression_single_mutant_trial30.csv")) == 30
        assert len(_csv(output / "expression_single_mutant_reserve.csv")) == 7
        assert (output / "expression_single_mutant_selection.png").stat().st_size > 1000
        assert (output / "expression_single_mutant_selection.svg").stat().st_size > 1000
