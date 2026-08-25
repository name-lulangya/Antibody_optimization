import csv
import importlib.util
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from antibody_optimization.vhh_expert_review import (  # noqa: E402
    VHHExpertReviewError,
    build_expert_review_rows,
    derive_position_contexts,
)
from antibody_optimization.vhh_expert_review_assessments import (  # noqa: E402
    get_all_v3_expert_assessments,
    validate_v3_expert_assessments,
)


SHORTLIST = (
    ROOT
    / "docs/result_artifacts/candidate_design"
    / "expression_single_mutant_selection_v3_20260825"
    / "expression_single_mutant_v3_final30.csv"
)
MAPPING = (
    ROOT
    / "docs/result_artifacts/input_baseline/structure_released_20260810"
    / "nb252_sequence_structure_mapping.csv"
)
EXPERIMENTAL_CIF = ROOT / "data/structures/cxs_exports/NK2R-252__native.cif"
AF3_CIF = ROOT / "data/structures/cxs_exports/fold_2r_252_nomg_model_0__native.cif"
RENDER_SCRIPT = (
    ROOT
    / "scripts/candidate_design/render_v3_parent_single_structure_review_chimerax.py"
)
REVIEW_DIR = (
    ROOT
    / "docs/result_artifacts/candidate_design"
    / "v3_parent_single_expert_review_20260825"
)
VISUAL_MANIFEST = REVIEW_DIR / "structure_views/structure_review_views.csv"
REVIEW_CSV = REVIEW_DIR / "v3_parent_single_expert_review.csv"


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _real_inputs():
    candidates = _csv(SHORTLIST)
    positions = sorted(
        {int(row["reported_sequence_index_1based"]) for row in candidates}
    )
    contexts = derive_position_contexts(
        _csv(MAPPING),
        positions,
        EXPERIMENTAL_CIF,
        AF3_CIF,
    )
    return candidates, contexts


def _test_assessments(
    candidates: list[dict[str, str]],
    contexts: dict[int, dict[str, object]],
) -> dict[str, dict[str, str]]:
    """Supply complete test-only judgements without encoding a real selection."""

    assessments: dict[str, dict[str, str]] = {}
    for row in candidates:
        position = int(row["reported_sequence_index_1based"])
        af3_only = (
            contexts[position]["primary_structure_source"]
            == "af3_only_due_missing_experimental_coordinates"
        )
        assessments[row["candidate_id"]] = {
            "structural_facts_cn": "测试专用结构事实。",
            "chimerax_single_rotamer_observation_cn": "测试专用可视化观察。",
            "vhh_expert_inference_cn": "测试专用专家推断。",
            "expert_structural_assessment": "reasonable_with_caution",
            "expert_solubility_expectation": "neutral_or_uncertain",
            "expert_thermal_stability_expectation": "neutral_or_uncertain",
            "expert_confidence": "low" if af3_only else "medium",
            "expert_primary_concern": "test_only_structural_context",
            "expert_rationale_cn": "仅用于测试表连接和证据来源，不代表正式专家结论。",
            "expert_uncertainty_cn": "未在测试中选择父单突。",
            "expert_rule_flags": ("test_only",),
        }
    return assessments


@pytest.fixture(scope="module")
def real_review():
    candidates, contexts = _real_inputs()
    rows = build_expert_review_rows(
        candidates,
        contexts,
        _test_assessments(candidates, contexts),
    )
    return candidates, contexts, rows


def test_real_review_has_30_candidates_23_positions_and_source_coverage(real_review):
    candidates, contexts, rows = real_review

    assert len(candidates) == len(rows) == 30
    assert len(contexts) == 23
    assert len({row["candidate_id"] for row in rows}) == 30
    assert len({int(row["reported_sequence_index_1based"]) for row in rows}) == 23

    sources = [row["primary_structure_source"] for row in rows]
    assert sources.count("experimental_complex") == 26
    assert sources.count("af3_only_due_missing_experimental_coordinates") == 4

    af3_only = {
        row["candidate_id"]
        for row in rows
        if row["primary_structure_source"]
        == "af3_only_due_missing_experimental_coordinates"
    }
    assert af3_only == {
        "Nb252_expr_seq011_L11Y",
        "Nb252_expr_seq011_L11M",
        "Nb252_expr_seq028_I28Y",
        "Nb252_expr_seq029_F29Q",
    }
    assert all(
        row["experimental_coordinate_status"] == "missing_coordinates"
        and row["structure_evidence_limit"]
        == "predicted_context_only_not_experimentally_observed"
        and row["expert_confidence"] == "low"
        for row in rows
        if row["candidate_id"] in af3_only
    )


def test_review_does_not_select_parents_or_change_upstream_evidence(real_review):
    candidates, _, rows = real_review
    source_by_id = {row["candidate_id"]: row for row in candidates}

    assert {row["parent_single_selection_status"] for row in rows} == {
        "not_performed"
    }
    for row in rows:
        source = source_by_id[row["candidate_id"]]
        assert row["sequence"] == source["sequence"]
        assert row["wt_residue"] == source["wt_residue"]
        assert row["mutant_residue"] == source["mutant_residue"]
        assert row["netsolp_delta_u"] == pytest.approx(
            float(source["netsolp_delta_usability_vs_current_wt"])
        )
        assert row["netsolp_u_band_v3"] == source["netsolp_u_magnitude_band_v3"]
        assert row["netsolp_delta_s"] == pytest.approx(
            float(source["netsolp_delta_solubility_vs_current_wt"])
        )
        assert row["netsolp_s_band_v3"] == source["netsolp_s_magnitude_band_v3"]
        assert row["nanomelt_delta_tm_c"] == pytest.approx(
            float(source["nanomelt_delta_predicted_apparent_tm_c_vs_current_wt"])
        )
        assert row["nanomelt_tm_band_v3"] == source["nanomelt_tm_magnitude_band_v3"]
        assert row["antifold_selection_source"] == source["antifold_selection_source"]
        assert row["antifold_delta_logp"] == pytest.approx(
            float(source["antifold_selection_delta_log_probability"])
        )
        assert row["antifold_mutant_rank_worst_first"] == int(
            source["antifold_mutant_rank_worst_first"]
        )
        assert row["antifold_veto_status"] == source["antifold_veto_status"]
        assert row["upstream_hard_sequence_risk_flags"] == source[
            "hard_sequence_risk_flags_v3"
        ]
        assert row["upstream_soft_sequence_risk_flags"] == source[
            "soft_sequence_risk_flags_v3"
        ]
        assert row["stable_word_effect"] == source["stable_word_effect"]


def test_real_parent_structures_recover_decision_relevant_environments(real_review):
    _, contexts, _ = real_review

    assert contexts[49]["primary_structure_source"] == "experimental_complex"
    assert contexts[49]["primary_exposure_class"] == "buried"
    assert float(contexts[49]["primary_relative_sasa"]) < 0.15

    assert contexts[96]["primary_structure_source"] == "experimental_complex"
    assert contexts[96]["primary_exposure_class"] == "buried"
    assert float(contexts[96]["primary_relative_sasa"]) < 0.15

    assert float(contexts[76]["primary_phi_degrees"]) > 0.0
    assert contexts[76]["primary_backbone_class_heuristic"] == (
        "loop_or_turn_like_phi_psi"
    )

    assert contexts[50]["primary_structure_source"] == "experimental_complex"
    assert float(contexts[50]["experimental_minimum_receptor_distance_a"]) < 4.5
    assert contexts[50]["experimental_nearest_receptor_residue"]


def test_review_rejects_assessment_identity_mismatch(real_review):
    candidates, contexts, _ = real_review
    assessments = _test_assessments(candidates, contexts)
    assessments.pop("Nb252_expr_seq001_Q1A")
    assessments["not_a_v3_candidate"] = {
        "structural_facts_cn": "测试用错误身份。",
        "chimerax_single_rotamer_observation_cn": "测试用错误身份。",
        "vhh_expert_inference_cn": "测试用错误身份。",
        "expert_structural_assessment": "indeterminate",
        "expert_solubility_expectation": "neutral_or_uncertain",
        "expert_thermal_stability_expectation": "neutral_or_uncertain",
        "expert_confidence": "low",
        "expert_primary_concern": "identity_mismatch",
        "expert_rationale_cn": "测试用错误身份。",
        "expert_uncertainty_cn": "不适用。",
        "expert_rule_flags": ("test_only",),
    }

    with pytest.raises(VHHExpertReviewError, match="Assessment identity mismatch"):
        build_expert_review_rows(candidates, contexts, assessments)


def test_review_rejects_candidate_sequence_identity_mismatch(real_review):
    candidates, contexts, _ = real_review
    changed = [dict(row) for row in candidates]
    changed[0]["sequence"] = changed[0]["sequence"][:-1]

    with pytest.raises(VHHExpertReviewError, match="Candidate sequence mismatch"):
        build_expert_review_rows(
            changed,
            contexts,
            _test_assessments(changed, contexts),
        )


def test_curated_assessments_cover_all_candidates_and_visual_observations():
    candidates = _csv(SHORTLIST)
    identifiers = [row["candidate_id"] for row in candidates]
    validate_v3_expert_assessments(identifiers)
    assessments = get_all_v3_expert_assessments()

    assert set(assessments) == set(identifiers)
    assert all(
        assessment["chimerax_single_rotamer_observation_cn"].strip()
        for assessment in assessments.values()
    )
    assert all("selected" not in assessment for assessment in assessments.values())
    assert all("rank" not in assessment for assessment in assessments.values())


def test_render_plan_and_released_visual_artifacts_cover_all_30_candidates():
    spec = importlib.util.spec_from_file_location("v3_expert_review_render", RENDER_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    candidates = _csv(SHORTLIST)
    mappings = module._mapping(_csv(MAPPING))
    plan = module.build_review_plan(candidates, mappings)

    assert len(plan["site_views"]) == 23
    assert plan["primary_candidate_view_count"] == 30
    assert plan["af3_sensitivity_view_count"] == 6
    primary_plan = [
        row for row in plan["candidate_views"] if row["view_kind"] == "candidate_primary"
    ]
    assert sum(row["source_model_name"] == module.EXPERIMENTAL_MODEL for row in primary_plan) == 26
    assert sum(row["source_model_name"] == module.AF3_MODEL for row in primary_plan) == 4
    assert len({row["view_id"] for row in plan["candidate_views"]}) == 36
    assert len({row["image_path"] for row in plan["candidate_views"]}) == 36

    visual_rows = _csv(VISUAL_MANIFEST)
    assert len(visual_rows) == 61
    primary_rows = [row for row in visual_rows if row["view_kind"] == "candidate_primary"]
    assert len(primary_rows) == 30
    assert {row["candidate_id"] for row in primary_rows} == {
        row["candidate_id"] for row in candidates
    }
    for row in visual_rows:
        image = VISUAL_MANIFEST.parent / row["image_path"]
        assert image.is_file() and image.stat().st_size > 0
        assert row["molecular_structure_saved"] == "False"
        assert row["candidate_selection_performed"] == "False"


def test_released_expert_review_is_non_ranking_and_fully_visualized():
    rows = _csv(REVIEW_CSV)

    assert len(rows) == 30
    assert {row["parent_single_selection_status"] for row in rows} == {
        "not_performed"
    }
    assert {row["manual_visual_review_status"] for row in rows} == {
        "reviewed_in_chimerax_1_12_single_rotamer_view"
    }
    assert all(row["chimerax_single_rotamer_observation_cn"].strip() for row in rows)
    assert sum(
        row["primary_structure_source"] == "experimental_complex" for row in rows
    ) == 26
    assert sum(
        row["primary_structure_source"]
        == "af3_only_due_missing_experimental_coordinates"
        for row in rows
    ) == 4
