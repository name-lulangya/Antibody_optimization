from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/reporting/audit_v3_release.py"


def _module():
    spec = importlib.util.spec_from_file_location("audit_v3_release", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_v3_release_audit_reconstructs_the_complete_identity_chain() -> None:
    result = _module().build_audit(
        ROOT,
        "2026-08-25T21:10:00+08:00",
        "862394d78d229618d13048229457f9be1ed2f759",
    )

    assert result["status"] == "pass_with_material_caveats"
    assert result["machine_check_summary"] == {
        "total": 18,
        "passed": 18,
        "failed": 0,
        "failed_names": [],
    }
    assert result["source_identity"]["authoritative_parent_length"] == 128
    assert (
        result["source_identity"]["authoritative_parent_sha256"]
        == "df5b83ddde8a3486383c12afe45e22af6a358f507eab5503d5dbd4430710288d"
    )
    assert result["stage_counts"] == {
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
    }


def test_v3_release_audit_records_material_caveats_without_hiding_them() -> None:
    result = _module().build_audit(
        ROOT,
        "2026-08-25T21:10:00+08:00",
        "862394d78d229618d13048229457f9be1ed2f759",
    )
    risk = result["final_panel_risk_and_evidence"]

    assert risk["final_constructs_using_AF3_only_position_evidence"] == 7
    assert risk["final_constructs_with_recorded_soft_sequence_liability"] == 4
    assert risk["final_constructs_with_reported_position_30_mutation"] == 6
    assert len(risk["selected_parents_with_antifold_delta_logp_le_minus3_but_gate_pass"]) == 8
    assert risk["final_constructs_with_at_least_one_such_antifold_component"] == 20
    assert risk["double_sidechain_modeling_performed"] is False
    assert risk["selected_t99f_double_count"] == 0

    provenance = result["provenance"]
    assert provenance["manifest_hash_mismatch_count"] == 0
    assert len(provenance["missing_remote_raw_bindings"]) == 4
    assert provenance["unexpected_missing_bindings"] == []

    validity = result["predictor_validity_and_dependence"]
    assert validity["netsolp_yield_evidence_level"] == "compatibility_filter_only"
    assert validity["nanomelt_yield_evidence_level"] == "no_supported_use"
    assert validity["antifold_yield_classification_status"] == "not_applicable"
    assert 0.26 < validity["netsolp_u_s_single847_spearman"] < 0.28
    assert 0.29 < validity["netsolp_u_s_double102_spearman"] < 0.31

    conflicts = result["machine_readable_semantic_conflicts"]
    assert conflicts["upstream30_gate_still_claims_final_experimental_release"] is True
    assert conflicts["old_critical_interface_mutation_semantics"] == "cautious_not_forbidden"
    assert conflicts["current_constraint_freezes_all_24_interface_positions"] is True

    readiness = result["report_readiness"]
    assert readiness["existing_report_is_historical_v2"] is True
    assert (
        readiness["historical_v2_materials_status"]
        == "expected_read_only_provenance_not_a_defect_or_blocker"
    )
    assert readiness["existing_report_active_route"] == "BL21_expression_19_single_plus_11_double"
    assert readiness["existing_report_single_mutant_count"] == 19
    assert readiness["existing_report_double_mutant_count"] == 11
    assert readiness["existing_report_source_double_candidate_count"] == 162
    assert readiness["existing_delivery_contains_v2_parent19_and_selected11_files"] is True
    assert readiness["existing_report_builder_points_to_historical_v2_inputs"] is True
    assert "current_external_materials_release_blocked" not in readiness
    assert readiness["v3_report_generation_status"] == "not_started"
    assert (
        readiness["v3_report_directory"]
        == "docs/result_artifacts/weekly_report_result/Nb252_V3_expression_report"
    )
    assert readiness["v3_audit_report"].endswith("/Nb252_V3_audit_report.md")
    assert readiness["v3_audit_evidence"].endswith("/Nb252_V3_audit_evidence.json")
    assert (ROOT / readiness["v3_audit_report"]).is_file()
    assert readiness["new_v3_report_drafting_allowed"] is True
    assert readiness["new_v3_report_finalization_allowed"] is False
    assert (
        readiness["v3_report_finalization_status"]
        == "pending_generation_identity_check_and_render_QA"
    )
