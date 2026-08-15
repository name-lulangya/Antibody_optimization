from __future__ import annotations

import csv
import json
from pathlib import Path

from antibody_optimization.unified_tnp_plot import render_unified_tnp_review
from antibody_optimization.unified_tnp_review import (
    BLOCKED_PRODUCTION_CYS_IDS,
    SCORE_COUNT,
    TNP_FLAGS,
    TNP_METRICS,
    analyze_unified_tnp_scores,
    build_unified_tnp_samples,
    magnitude_label,
)


ROOT = Path(__file__).resolve().parents[1]
PROPERTY_RESULT = ROOT / "docs/result_artifacts/candidate_design/unified_property_scoring_result_20260815"
PROPERTY_PLAN = ROOT / "docs/result_artifacts/candidate_design/unified_property_scoring_plan_20260815"
FLEX_RESULT = ROOT / "docs/result_artifacts/candidate_design/flex_ddg_production_result_20260812"
UNIFIED_PLAN = ROOT / "docs/result_artifacts/candidate_design/unified_single_mutant_plan_20260815"


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _real_samples():
    return build_unified_tnp_samples(
        _csv(PROPERTY_RESULT / "unified_single_mutant_property_evidence.csv"),
        _csv(PROPERTY_PLAN / "unified_property_samples.csv"),
        _csv(FLEX_RESULT / "flex_ddg_production_candidate_summary.csv"),
        _csv(UNIFIED_PLAN / "unified_single_mutant_candidates.csv"),
    )


def test_real_unified_tnp_plan_has_exact_released_pools():
    samples, blocked = _real_samples()
    assert len(samples) == SCORE_COUNT
    assert samples[0]["sample_uid"] == "LTT__Nb252__WT"
    assert sum(row["candidate_source"] == "property_pareto_front_1" for row in samples) == 49
    assert sum(row["candidate_source"] == "affinity_flex_ddg_20_sample_pool" for row in samples) == 46
    assert {row["candidate_id"] for row in blocked} == BLOCKED_PRODUCTION_CYS_IDS
    parent = samples[0]["sequence_raw"]
    for row in samples:
        assert len(row["sequence_raw"]) == 128
        assert row["sequence_raw"].endswith("SSGS")
        assert row["sequence_raw"][21] == row["sequence_raw"][94] == "C"
        if row is not samples[0]:
            assert sum(a != b for a, b in zip(parent, row["sequence_raw"], strict=True)) == 1


def test_magnitude_labels_use_inclusive_neutral_boundaries():
    assert magnitude_label(-0.010001, 0.01) == "adverse"
    assert magnitude_label(-0.01, 0.01) == "neutral"
    assert magnitude_label(0.01, 0.01) == "neutral"
    assert magnitude_label(0.010001, 0.01) == "favorable"


def test_unified_tnp_analysis_preserves_all_metrics_and_risk_semantics(tmp_path):
    samples, _ = _real_samples()
    scores = []
    for index, sample in enumerate(samples):
        score = {
            "sample_uid": sample["sample_uid"],
            "sequence_raw": sample["sequence_raw"],
            "modelled_length_aa": 126,
            "trimmed_n_terminal": "",
            "trimmed_c_terminal": "GS",
            "scoring_status": "pass",
        }
        for metric_index, metric in enumerate(TNP_METRICS):
            score[metric] = 20.0 + metric_index + index / 100.0
        for flag in TNP_FLAGS:
            score[flag] = "green"
        score["tnp_flag_psh"] = "amber"
        if index == 1:
            score["tnp_flag_psh"] = "green"
        elif index == 2:
            score["tnp_flag_ppc"] = "amber"
        elif index == 3:
            score["tnp_flag_pnc"] = "red"
        scores.append(score)

    evidence, summary, gate = analyze_unified_tnp_scores(samples, scores)
    assert len(evidence) == 95
    assert sum(row["tnp_developability_review"] == "flag_improvement" for row in evidence) == 1
    assert sum(row["tnp_developability_review"] == "flag_regression" for row in evidence) == 1
    assert sum(row["tnp_developability_review"] == "new_red_flag" for row in evidence) == 1
    assert gate["status"] == "pass"
    assert gate["candidate_selection_performed"] is False
    assert sum(row["candidate_count"] for row in summary) == 95
    render_unified_tnp_review(
        evidence,
        png_path=tmp_path / "unified_tnp.png",
        svg_path=tmp_path / "unified_tnp.svg",
    )
    assert (tmp_path / "unified_tnp.png").stat().st_size > 0
    assert (tmp_path / "unified_tnp.svg").stat().st_size > 0


def test_unified_tnp_slurm_contract_is_single_process_and_non_resumable():
    slurm = (ROOT / "scripts/candidate_design/submit_unified_tnp_review.slurm").read_text(encoding="utf-8")
    wrapper = (ROOT / "scripts/candidate_design/submit_unified_tnp_review.sh").read_text(encoding="utf-8")
    assert "#SBATCH --partition=batch" in slurm
    assert "#SBATCH --gres=gpu:1" in slurm
    assert "#SBATCH --cpus-per-task=12" in slurm
    assert "#SBATCH --time=02:00:00" in slurm
    assert "logs/unified_tnp_review/run-%j.log" in slurm
    assert "--array" not in slurm
    assert slurm.index("conda activate /data/software/env/luly25/tnp") < slurm.index("set -u")
    assert "Score directory already exists" in wrapper
    assert "resume" not in wrapper.lower()
