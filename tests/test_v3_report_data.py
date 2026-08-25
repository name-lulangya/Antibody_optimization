import copy
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from antibody_optimization.v3_report_data import (  # noqa: E402
    ANTIFOLD_REPORT_ROLE,
    AUTHORITATIVE_PARENT_SHA256,
    V3ReportDataError,
    load_v3_report_data,
    validate_v3_report_data,
)


EXPECTED_PARENT_MUTATIONS = [
    "L11Y",
    "F30S",
    "K86S",
    "A23R",
    "Q5V",
    "S55G",
    "K75A",
    "F29Q",
    "K43A",
    "N76G",
    "F30N",
    "K75E",
    "L11M",
    "Q1D",
    "T99F",
]

EXPECTED_DOUBLE_MUTATIONS = [
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
]


@pytest.fixture(scope="module")
def report_data():
    return load_v3_report_data(ROOT)


def test_loader_exposes_stable_report_api_and_authoritative_parent(report_data):
    required = {
        "parent_sequence",
        "parent_hash",
        "counts",
        "constraints",
        "conservation",
        "tool_validation",
        "magnitude_thresholds",
        "single_audit847",
        "parent_selected15",
        "double_audit102",
        "selected_doubles15",
        "final_panel30",
        "audit_evidence",
        "antifold_policy",
    }
    assert required <= report_data.keys()
    assert report_data["parent_hash"] == AUTHORITATIVE_PARENT_SHA256
    assert len(report_data["parent_sequence"]) == 128
    assert report_data["parent_sequence"].endswith("SSGS")
    assert report_data["parent_sequence"].count("C") == 2


def test_report_chain_preserves_released_v3_counts_and_identities(report_data):
    assert report_data["counts"] == {
        "allowed_single_mutants": 847,
        "upstream_single_shortlist": 30,
        "parent_expert_review_pool": 31,
        "selected_parent_singles": 15,
        "selected_parent_reported_positions": 12,
        "theoretical_parent_pairs": 105,
        "invalid_same_position_pairs": 3,
        "valid_double_mutants": 102,
        "selected_double_mutants": 15,
        "final_panel_sequences": 30,
        "hard_frozen_positions": 80,
        "experimental_interface_positions": 24,
        "antifold_vetoed_singles": 151,
        "antifold_veto_pass_singles": 696,
        "v3_qualified_singles": 61,
    }
    assert len(report_data["single_audit847"]) == 847
    assert len(report_data["double_audit102"]) == 102
    assert [
        row["mutation_reported_label"].split()[-1]
        for row in report_data["parent_selected15"]
    ] == EXPECTED_PARENT_MUTATIONS
    assert [row["mutation_set"] for row in report_data["selected_doubles15"]] == (
        EXPECTED_DOUBLE_MUTATIONS
    )
    assert len(report_data["final_panel30"]) == 30
    assert len({row["candidate_id"] for row in report_data["final_panel30"]}) == 30
    assert len({row["sequence"] for row in report_data["final_panel30"]}) == 30


def test_antifold_is_only_a_negative_single_state_veto(report_data):
    policy = report_data["antifold_policy"]
    assert policy == ANTIFOLD_REPORT_ROLE
    assert policy["role"] == "negative_risk_exclusion_only"
    assert policy["positive_candidate_credit"] is False
    assert policy["proposes_candidates"] is False
    assert policy["ranks_candidates"] is False
    assert policy["double_mutant_scoring"] == "not_performed"
    assert policy["double_mutant_use"] == (
        "constituent_single_mutant_veto_evidence_only"
    )
    assert policy["component_values_combined"] is False
    assert all(
        row["antifold_constituent_gate"] == "pass"
        and row["antifold_double_mutant_scored"] == "False"
        and row["antifold_component_values_combined"] == "False"
        and row["antifold_double_mutant_score"] == ""
        for row in report_data["double_audit102"]
    )


def test_tool_validation_summary_preserves_declared_scope(report_data):
    tools = report_data["tool_validation"]
    assert tools["netsolp"]["sample_count"] == 47
    assert tools["netsolp"]["evidence_level"] == "compatibility_filter_only"
    assert tools["nanomelt"]["planned_count"] == 47
    assert tools["nanomelt"]["scored_count"] == 43
    assert tools["nanomelt"]["evidence_level"] == "no_supported_use"
    assert tools["rp3net"]["candidate_role"] == "not_used"
    assert tools["plm_sol"]["candidate_role"] == "not_used"
    assert tools["plm_sol"]["increment_over_netsolp_s"] == pytest.approx(
        -0.1404500629533707
    )


def test_conservation_and_audit_evidence_are_report_ready(report_data):
    conservation = report_data["conservation"]
    assert conservation["source_record_count"] == 4059
    assert conservation["eligible_sequence_count"] == 4057
    assert conservation["redundancy_cluster_count"] == 3784
    assert conservation["neighbor_sequence_count"] == 1564
    assert conservation["neighbor_cluster_count"] == 1532
    assert report_data["audit_evidence"]["machine_check_summary"] == {
        "total": 18,
        "passed": 18,
        "failed": 0,
        "failed_names": [],
    }
    assert report_data["audit_evidence"]["status"] == "pass_with_material_caveats"


def test_renderer_boundary_rejects_any_positive_or_double_antifold_role(report_data):
    modified = copy.deepcopy(report_data)
    modified["tool_validation"]["antifold"]["proposes_candidates"] = True
    with pytest.raises(V3ReportDataError, match="AntiFold report role was altered"):
        validate_v3_report_data(modified)

    modified = copy.deepcopy(report_data)
    modified["tool_validation"]["antifold"]["double_mutant_scoring"] = "performed"
    with pytest.raises(V3ReportDataError, match="AntiFold report role was altered"):
        validate_v3_report_data(modified)
