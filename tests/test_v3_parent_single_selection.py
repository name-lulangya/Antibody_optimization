import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from antibody_optimization.v3_parent_single_selection import (  # noqa: E402
    HARD_EXPERT_EXCLUSION_IDS,
    SELECTED_PARENT_IDS,
    V3ParentSingleSelectionError,
    V3_T99F_ID,
    build_v3_parent_single_selection,
)
from antibody_optimization.v3_parent_single_selection_plot import (  # noqa: E402
    build_v3_parent_selection_plot_rows,
)


REVIEW = (
    ROOT
    / "docs/result_artifacts/candidate_design"
    / "v3_parent_single_expert_review_20260825"
    / "v3_parent_single_expert_review.csv"
)
COMPLETE_AUDIT = (
    ROOT
    / "docs/result_artifacts/candidate_design"
    / "expression_single_mutant_selection_v3_20260825"
    / "expression_single_mutant_v3_audit.csv"
)
SCRIPT = ROOT / "scripts/candidate_design/select_v3_parent_single_panel.py"


def _csv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _result():
    return build_v3_parent_single_selection(_csv(REVIEW), _csv(COMPLETE_AUDIT))


def test_real_review_pool_yields_exact_15_parents_and_102_valid_pairs():
    result = _result()
    facts = result["facts"]
    selected = result["selected_rows"]
    assert {row["candidate_id"] for row in selected} == set(SELECTED_PARENT_IDS)
    assert facts == {
        "review_candidate_count": 31,
        "selected_parent_single_count": 15,
        "not_selected_candidate_count": 16,
        "high_confidence_expert_risk_exclusion_count": 5,
        "competitive_not_selected_count": 11,
        "selected_with_strong_favorable_metric_count": 9,
        "selected_user_directed_exploration_count": 1,
        "selected_unique_position_count": 12,
        "selected_same_position_groups": {"11": 2, "30": 2, "75": 2},
        "theoretical_unordered_pair_count": 105,
        "invalid_same_position_pair_count": 3,
        "valid_unordered_double_mutant_count": 102,
        "double_mutant_enumeration_performed": False,
    }
    assert len({row["sequence"] for row in selected}) == 15
    assert all(len(row["sequence"]) == 128 for row in selected)
    assert all(row["sequence"].endswith("SSGS") for row in selected)
    assert all(row["sequence"].count("C") == 2 for row in selected)
    assert all(
        any(
            row[field] in {"moderate_favorable", "strong_favorable"}
            for field in (
                "netsolp_u_band_v3",
                "netsolp_s_band_v3",
                "nanomelt_tm_band_v3",
            )
        )
        for row in selected
        if row["candidate_id"] != V3_T99F_ID
    )


def test_every_decision_preserves_expert_review_evidence():
    review = {row["candidate_id"]: row for row in _csv(REVIEW)}
    result = _result()
    fields = {
        "candidate_id",
        "mutation_reported_label",
        "reported_sequence_index_1based",
        "wt_residue",
        "mutant_residue",
        "sequence",
        "netsolp_delta_u",
        "netsolp_u_band_v3",
        "netsolp_delta_s",
        "netsolp_s_band_v3",
        "nanomelt_delta_tm_c",
        "nanomelt_tm_band_v3",
        "antifold_selection_source",
        "antifold_delta_logp",
        "antifold_mutant_rank_worst_first",
        "antifold_veto_status",
        "stable_word_effect",
        "expert_structural_assessment",
        "expert_solubility_expectation",
        "expert_thermal_stability_expectation",
        "expert_confidence",
        "expert_primary_concern",
        "expert_rationale_cn",
    }
    assert len(result["audit_rows"]) == 31
    for row in result["audit_rows"]:
        source = review[row["candidate_id"]]
        assert {field: row[field] for field in fields} == {
            field: source[field] for field in fields
        }
        assert row["v3_parent_decision_reason_cn"]


def test_t99f_remains_an_explicit_upstream_ineligible_exploration_exception():
    row = next(
        row for row in _result()["audit_rows"] if row["candidate_id"] == V3_T99F_ID
    )
    assert row["v3_parent_selection_status"] == "selected"
    assert row["v3_parent_decision_class"] == (
        "selected_user_directed_stable_word_exploration"
    )
    assert row["stable_word_effect"] == "gain_only"
    assert row["upstream_selection_tier_v3"] == "not_eligible"
    assert row["upstream_selection_status_v3"] == "not_selected"
    assert all(
        row[field] not in {"moderate_favorable", "strong_favorable"}
        for field in (
            "netsolp_u_band_v3",
            "netsolp_s_band_v3",
            "nanomelt_tm_band_v3",
        )
    )


def test_only_five_high_confidence_concrete_expert_risks_are_hard_excluded():
    rows = _result()["audit_rows"]
    hard = {
        row["candidate_id"]
        for row in rows
        if row["v3_parent_high_confidence_expert_risk_exclusion"] is True
    }
    assert hard == HARD_EXPERT_EXCLUSION_IDS
    assert all(
        row["expert_structural_assessment"] == "structurally_concerning"
        and row["expert_confidence"] == "high"
        and row["v3_parent_selection_status"] == "not_selected"
        for row in rows
        if row["candidate_id"] in hard
    )
    competitive = {
        row["candidate_id"]
        for row in rows
        if row["v3_parent_selection_status"] == "not_selected"
        and row["candidate_id"] not in hard
    }
    assert len(competitive) == 11
    assert all(
        row["v3_parent_decision_class"]
        != "rejected_high_confidence_expert_risk"
        for row in rows
        if row["candidate_id"] in competitive
    )


def test_compact_plot_data_preserves_decisions_and_neutralizes_weak_changes():
    rows = build_v3_parent_selection_plot_rows(_result()["audit_rows"])
    assert len(rows) == 31
    assert sum(row["plot_disposition"] == "selected_parent" for row in rows) == 15
    assert sum(
        row["plot_disposition"] == "high_confidence_expert_risk_exclusion"
        for row in rows
    ) == 5
    assert sum(
        row["plot_disposition"] == "competitive_not_selected" for row in rows
    ) == 11
    t99f = next(row for row in rows if row["candidate_id"] == V3_T99F_ID)
    assert t99f["stable_word_effect"] == "gain_only"
    assert t99f["netsolp_u_display_grade"] == 0
    assert t99f["netsolp_s_display_grade"] == 0
    assert t99f["nanomelt_tm_display_grade"] == 0


def test_selection_rejects_missing_review_candidate_and_sequence_mismatch():
    review = _csv(REVIEW)
    complete = _csv(COMPLETE_AUDIT)
    with pytest.raises(V3ParentSingleSelectionError, match="31 expert-review"):
        build_v3_parent_single_selection(review[:-1], complete)
    mutated = [dict(row) for row in complete]
    source = next(row for row in mutated if row["candidate_id"] == V3_T99F_ID)
    source["sequence"] = source["sequence"][:-1] + "A"
    with pytest.raises(V3ParentSingleSelectionError, match="sequence mismatch"):
        build_v3_parent_single_selection(review, mutated)


def test_cli_writes_new_artifacts_without_modifying_upstream_reviews():
    original_review = REVIEW.read_bytes()
    original_audit = COMPLETE_AUDIT.read_bytes()
    with tempfile.TemporaryDirectory(prefix=".test-v3-parent15-", dir=ROOT) as temp:
        output = Path(temp) / "result"
        run_summary = Path(temp) / "run_summary.json"
        subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--expert-review-csv",
                str(REVIEW),
                "--complete-v3-audit-csv",
                str(COMPLETE_AUDIT),
                "--output-dir",
                str(output),
                "--run-summary",
                str(run_summary),
                "--generated-at",
                "2026-08-25T16:00:00+08:00",
            ],
            check=True,
            cwd=ROOT,
        )
        assert len(_csv(output / "v3_parent_single_selection_audit.csv")) == 31
        assert len(_csv(output / "v3_parent_single_selected15.csv")) == 15
        assert len(_csv(output / "v3_parent_single_selection_plot_data.csv")) == 31
        assert (output / "v3_parent_single_selection_overview.png").stat().st_size > 0
        assert (output / "v3_parent_single_selection_overview.svg").stat().st_size > 0
        fasta = (output / "v3_parent_single_selected15.fasta").read_text(
            encoding="utf-8"
        )
        assert fasta.count(">") == 15
        manifest = json.loads(
            (output / "v3_parent_single_selection_manifest.json").read_text(
                encoding="utf-8"
            )
        )
        assert manifest["status"] == "pass"
        assert manifest["gate"] == {
            "v3_parent_single_selection": "pass",
            "selected_parent_count": 15,
            "detailed_decision_audit_count": 31,
            "valid_double_mutant_space_ready_for_generation": True,
            "double_mutant_enumeration": "not_performed",
            "final_15_double_mutant_selection": "not_performed",
        }
        assert manifest["facts"]["valid_unordered_double_mutant_count"] == 102
        summary = json.loads(run_summary.read_text(encoding="utf-8"))
        assert summary["status"] == "pass"
        assert summary["facts"]["selected_parent_single_count"] == 15
        assert summary["verification"]["plot_source_rows"] == 31
    assert REVIEW.read_bytes() == original_review
    assert COMPLETE_AUDIT.read_bytes() == original_audit
