from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
from collections import Counter
from functools import lru_cache
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from antibody_optimization.v3_double_mutant_selection import (  # noqa: E402
    build_v3_double_mutant_selection,
)


MATRIX_DIR = (
    ROOT
    / "docs/result_artifacts/candidate_design"
    / "v3_double_mutant_property_matrix_20260825"
)
PARENT_DIR = (
    ROOT
    / "docs/result_artifacts/candidate_design"
    / "v3_parent_single_selection_20260825"
)
POST_SYNC_REVIEW = (
    ROOT
    / "docs/result_artifacts/candidate_design"
    / "v3_double_mutant_post_sync_review_20260825"
    / "v3_double_mutant_post_sync_review.json"
)

EXPECTED_SELECTED_MUTATION_SETS = (
    "F30S;Q5V",
    "S55G;K43A",
    "K86S;Q5V",
    "L11Y;K86S",
    "S55G;F30N",
    "N76G;L11M",
    "F30S;K75E",
    "K86S;K43A",
    "A23R;S55G",
    "K43A;N76G",
    "K75E;Q1D",
    "L11Y;K75A",
    "Q5V;N76G",
    "L11Y;Q1D",
    "F30S;Q1D",
)

EXPECTED_OUTPUT_NAMES = {
    "audit": "v3_double_mutant_final_selection_audit102.csv",
    "selected": "v3_double_mutant_selected15.csv",
    "final_panel": "v3_final_panel30.csv",
    "fasta": "v3_final_panel30.fasta",
    "plot_data": "v3_final_panel_plot_data.csv",
    "png": "v3_final_panel_overview.png",
    "svg": "v3_final_panel_overview.svg",
    "manifest": "v3_final_panel_manifest.json",
}


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    assert isinstance(value, dict)
    return value


@lru_cache(maxsize=1)
def _result() -> dict[str, object]:
    return build_v3_double_mutant_selection(
        _csv(MATRIX_DIR / "v3_double_mutant_property_matrix102.csv"),
        _csv(PARENT_DIR / "v3_parent_single_selected15.csv"),
        _csv(PARENT_DIR / "v3_parent_single_selection_audit.csv"),
        _json(POST_SYNC_REVIEW),
    )


def _expected_enhanced_review(row: dict[str, object]) -> bool:
    favorable = int(row["moderate_or_strong_favorable_metric_count"])
    moderate_adverse = int(row["moderate_adverse_metric_count"])
    strong_adverse = int(row["strong_adverse_metric_count"])
    two_metric_support_without_adverse = (
        favorable >= 2 and moderate_adverse == 0 and strong_adverse == 0
    )
    adverse_tradeoff = moderate_adverse > 0 or strong_adverse > 0
    local_pair = (
        row["pair_spatial_class"]
        != "spatially_separated_ca_at_least_10A"
    )
    return two_metric_support_without_adverse or adverse_tradeoff or local_pair


def _bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() == "true"


def test_all_102_doubles_receive_common_expert_review_with_58_44_depths():
    audit = _result()["audit_rows"]
    assert len(audit) == 102
    assert len({str(row["double_candidate_id"]) for row in audit}) == 102
    assert Counter(str(row["expert_review_depth"]) for row in audit) == {
        "enhanced": 58,
        "standard": 44,
    }
    assert sum(
        int(row["moderate_or_strong_favorable_metric_count"]) >= 2
        and int(row["moderate_adverse_metric_count"]) == 0
        and int(row["strong_adverse_metric_count"]) == 0
        for row in audit
    ) == 42
    assert sum(
        int(row["moderate_adverse_metric_count"]) > 0
        or int(row["strong_adverse_metric_count"]) > 0
        for row in audit
    ) == 7
    assert sum(
        row["pair_spatial_class"] != "spatially_separated_ca_at_least_10A"
        for row in audit
    ) == 12
    for row in audit:
        expected = "enhanced" if _expected_enhanced_review(row) else "standard"
        assert row["expert_review_depth"] == expected
        assert _bool(row["expert_review_completed"])
        assert str(row["expert_review_triggers"]).strip()
        assert str(row["double_expert_assessment"]).strip()
        assert str(row["double_expert_confidence"]).strip()
        assert str(row["double_expert_rationale_cn"]).strip()
        assert str(row["final_double_decision_reason_cn"]).strip()


def test_t99f_has_no_mutation_specific_review_or_selection_rule():
    result = _result()
    audit = result["audit_rows"]
    t99f_rows = [
        row
        for row in audit
        if _bool(row["contains_t99f_stable_word_exploration_parent"])
    ]
    assert len(t99f_rows) == 14
    assert Counter(str(row["expert_review_depth"]) for row in t99f_rows) == {
        "enhanced": 2,
        "standard": 12,
    }
    for row in t99f_rows:
        expected = "enhanced" if _expected_enhanced_review(row) else "standard"
        assert row["expert_review_depth"] == expected
    assert result["selection_policy"]["mutation_specific_quota_or_exception"] is False
    assert result["selection_policy"]["stable_word_role"] == (
        "uniform_soft_evidence_only_not_a_selection_requirement"
    )


def test_legacy_review_flags_t99f_labels_and_within_band_decimals_do_not_drive_selection():
    matrix = _csv(MATRIX_DIR / "v3_double_mutant_property_matrix102.csv")
    for row in matrix:
        row["contains_t99f_stable_word_exploration_parent"] = "False"
        row["detailed_expert_review_required"] = "False"
        row["machine_structure_triage_status"] = "routine_context_recorded"
        row["machine_structure_triage_triggers"] = ""
        row["netsolp_u_delta_vs_wt"] = str(
            float(row["netsolp_u_delta_vs_wt"]) + 1e-8
        )
        row["netsolp_s_delta_vs_wt"] = str(
            float(row["netsolp_s_delta_vs_wt"]) - 1e-8
        )
        row["nanomelt_tm_c_delta_vs_wt"] = str(
            float(row["nanomelt_tm_c_delta_vs_wt"]) + 1e-8
        )
    perturbed = build_v3_double_mutant_selection(
        matrix,
        _csv(PARENT_DIR / "v3_parent_single_selected15.csv"),
        _csv(PARENT_DIR / "v3_parent_single_selection_audit.csv"),
        _json(POST_SYNC_REVIEW),
    )
    assert [row["mutation_set"] for row in perturbed["selected_double_rows"]] == list(
        EXPECTED_SELECTED_MUTATION_SETS
    )
    assert Counter(row["expert_review_depth"] for row in perturbed["audit_rows"]) == {
        "enhanced": 58,
        "standard": 44,
    }


def test_post_sync_deamidation_erratum_is_applied_without_rewriting_source_matrix():
    matrix = MATRIX_DIR / "v3_double_mutant_property_matrix102.csv"
    before = matrix.read_bytes()
    audit = _result()["audit_rows"]
    corrected = next(row for row in audit if row["mutation_set"] == "N76G;F30N")
    assert corrected["deamidation_motif_delta"] in {0, "0"}
    assert "new_deamidation_motif" in str(
        corrected["soft_sequence_risk_flags"]
    ).split("|")
    assert int(corrected["soft_sequence_risk_count"]) >= 1
    assert _bool(corrected["post_sync_annotation_erratum_applied"])
    assert matrix.read_bytes() == before


def test_selection_rejects_an_incomplete_or_duplicate_102_row_identity_space():
    matrix = _csv(MATRIX_DIR / "v3_double_mutant_property_matrix102.csv")
    inputs = (
        _csv(PARENT_DIR / "v3_parent_single_selected15.csv"),
        _csv(PARENT_DIR / "v3_parent_single_selection_audit.csv"),
        _json(POST_SYNC_REVIEW),
    )
    with pytest.raises(ValueError):
        build_v3_double_mutant_selection(matrix[:-1], *inputs)
    duplicate = [dict(row) for row in matrix]
    duplicate[-1] = dict(duplicate[0])
    with pytest.raises(ValueError):
        build_v3_double_mutant_selection(duplicate, *inputs)


def test_exact_15_doubles_are_frozen_in_non_efficacy_display_order():
    result = _result()
    selected = result["selected_double_rows"]
    audit = result["audit_rows"]
    assert [row["mutation_set"] for row in selected] == list(
        EXPECTED_SELECTED_MUTATION_SETS
    )
    assert [int(row["final_double_panel_order_not_efficacy_rank"]) for row in selected] == list(
        range(1, 16)
    )
    assert Counter(str(row["final_double_selection_status"]) for row in audit) == {
        "not_selected": 87,
        "selected": 15,
    }
    assert all(
        int(row["moderate_or_strong_favorable_metric_count"]) >= 2
        and int(row["moderate_adverse_metric_count"]) == 0
        and int(row["strong_adverse_metric_count"]) == 0
        and row["expert_review_depth"] == "enhanced"
        for row in selected
    )
    assert all(
        row["pair_spatial_class"] == "spatially_separated_ca_at_least_10A"
        for row in selected
    )
    assert all(row["antifold_constituent_gate"] == "pass" for row in selected)
    assert all(not _bool(row["antifold_double_mutant_scored"]) for row in selected)
    assert all(
        not _bool(row["antifold_component_values_combined"]) for row in selected
    )


def test_selected_15_satisfy_all_frozen_diversity_and_risk_caps():
    selected = _result()["selected_double_rows"]
    component_use = Counter(
        mutation
        for row in selected
        for mutation in (str(row["mutation_a"]), str(row["mutation_b"]))
    )
    position_use = Counter(
        int(position)
        for row in selected
        for position in (
            row["position_a_reported_1based"],
            row["position_b_reported_1based"],
        )
    )
    position_pairs = Counter(
        tuple(
            sorted(
                (
                    int(row["position_a_reported_1based"]),
                    int(row["position_b_reported_1based"]),
                )
            )
        )
        for row in selected
    )
    assert len(component_use) == 13
    assert max(component_use.values()) <= 3
    assert len(position_use) == 10
    assert max(position_use.values()) <= 4
    assert max(position_pairs.values()) == 1
    assert sum(
        row["pair_structure_distance_source"]
        == "af3_vhh_only_due_missing_experimental_coordinate"
        for row in selected
    ) == 4
    assert sum(int(row["soft_sequence_risk_count"]) > 0 for row in selected) == 2
    assert sum(
        row["pair_spatial_class"] != "spatially_separated_ca_at_least_10A"
        for row in selected
    ) == 0
    assert all(int(row["hard_sequence_risk_count"]) == 0 for row in selected)


def test_final_panel_contains_authoritative_15_singles_plus_selected_15_doubles():
    result = _result()
    final_rows = result["final_panel_rows"]
    parent_rows = _csv(PARENT_DIR / "v3_parent_single_selected15.csv")
    parent_audit = {
        row["candidate_id"]: row
        for row in _csv(PARENT_DIR / "v3_parent_single_selection_audit.csv")
        if row["v3_parent_selection_status"] == "selected"
    }
    selected = result["selected_double_rows"]
    assert len(final_rows) == 30
    assert Counter(str(row["candidate_kind"]) for row in final_rows) == {
        "single_mutant": 15,
        "double_mutant": 15,
    }
    assert [int(row["final_panel_order_not_efficacy_rank"]) for row in final_rows] == list(
        range(1, 31)
    )
    assert [row["candidate_id"] for row in final_rows[:15]] == [
        row["candidate_id"] for row in parent_rows
    ]
    assert [row["candidate_id"] for row in final_rows[15:]] == [
        row["double_candidate_id"] for row in selected
    ]
    sequences = [str(row["sequence"]) for row in final_rows]
    assert len(set(sequences)) == 30
    assert all(
        len(sequence) == 128
        and sequence.endswith("SSGS")
        and sequence.count("C") == 2
        for sequence in sequences
    )
    reconstructed_parents = set()
    for row in parent_rows:
        sequence = str(row["sequence"])
        source = parent_audit[row["candidate_id"]]
        position = int(source["reported_sequence_index_1based"])
        reconstructed_parents.add(
            sequence[: position - 1] + str(source["wt_residue"]) + sequence[position:]
        )
    assert len(reconstructed_parents) == 1
    parent = reconstructed_parents.pop()
    for row in final_rows:
        difference_count = sum(
            wt != mutant for wt, mutant in zip(parent, str(row["sequence"]), strict=True)
        )
        expected = 1 if row["candidate_kind"] == "single_mutant" else 2
        assert difference_count == expected


def test_selection_cli_writes_auditable_release_and_refuses_implicit_overwrite():
    protected = (
        MATRIX_DIR / "v3_double_mutant_property_matrix102.csv",
        MATRIX_DIR / "v3_double_mutant_property_matrix_manifest.json",
        PARENT_DIR / "v3_parent_single_selected15.csv",
        PARENT_DIR / "v3_parent_single_selection_audit.csv",
        POST_SYNC_REVIEW,
    )
    source_bytes = {path: path.read_bytes() for path in protected}
    with tempfile.TemporaryDirectory(prefix=".test-v3-final-panel-", dir=ROOT) as temp:
        work = Path(temp)
        output = work / "v3_final_15plus15_panel_20260825"
        summary = work / "run_summary.json"
        command = [
            sys.executable,
            str(ROOT / "scripts/candidate_design/select_v3_double_mutant_panel.py"),
            "--output-dir",
            str(output),
            "--run-summary",
            str(summary),
            "--generated-at",
            "2026-08-25T20:00:00+08:00",
        ]
        subprocess.run(command, cwd=ROOT, check=True)
        manifest = _json(output / EXPECTED_OUTPUT_NAMES["manifest"])
        run_summary = _json(summary)
        assert manifest["status"] == "pass"
        assert manifest["workflow"] == "v3_double_mutant_final_selection"
        assert manifest["gate"]["v3_double_expert_review"] == "pass"
        assert manifest["gate"]["final_15_double_mutant_selection"] == "pass"
        assert manifest["gate"]["final_30_panel_release"] == "pass"
        assert manifest["facts"]["source_double_candidate_count"] == 102
        assert manifest["facts"]["enhanced_expert_review_count"] == 58
        assert manifest["facts"]["standard_expert_review_count"] == 44
        assert manifest["facts"]["selected_double_mutant_count"] == 15
        assert manifest["facts"]["final_panel_candidate_count"] == 30
        assert manifest["selection_policy"]["mutation_specific_quota_or_exception"] is False
        assert manifest["selection_policy"]["within_band_raw_decimals_used_as_rank"] is False
        assert manifest["selection_policy"]["predictors_rerun"] is False
        assert manifest["selection_policy"]["antifold_role"] == (
            "constituent_negative_veto_only_no_double_score_no_positive_rank"
        )
        assert run_summary["status"] == "pass"
        assert run_summary["workflow"] == "v3_double_mutant_final_selection"
        assert len(_csv(output / EXPECTED_OUTPUT_NAMES["audit"])) == 102
        assert len(_csv(output / EXPECTED_OUTPUT_NAMES["selected"])) == 15
        assert len(_csv(output / EXPECTED_OUTPUT_NAMES["final_panel"])) == 30
        assert len(_csv(output / EXPECTED_OUTPUT_NAMES["plot_data"])) > 0
        assert (output / EXPECTED_OUTPUT_NAMES["png"]).stat().st_size > 1000
        assert (output / EXPECTED_OUTPUT_NAMES["svg"]).stat().st_size > 1000
        for key, record in manifest["outputs"].items():
            if key == "manifest":
                continue
            path = ROOT / record["path"]
            if not path.is_absolute():
                path = (ROOT / path).resolve()
            assert path.exists()
            assert record["sha256"]
        with pytest.raises(subprocess.CalledProcessError):
            subprocess.run(command, cwd=ROOT, check=True)
    assert {path: path.read_bytes() for path in protected} == source_bytes
