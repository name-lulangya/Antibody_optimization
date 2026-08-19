from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

from antibody_optimization.expression_property_completion import (
    WT_SCORE_ID,
    build_complete_score_matrix,
    build_reuse_plan,
    compare_repeat_scores,
)


ROOT = Path(__file__).resolve().parents[1]
V2 = ROOT / "docs/result_artifacts/input_baseline/vhh_conservation_consensus_v2_20260819"
PREFLIGHT = ROOT / "docs/result_artifacts/candidate_design/expression_single_mutant_contract_v2_20260819/expression_single_mutant_contract_preflight.json"
OLD_PROPERTY = ROOT / "docs/result_artifacts/candidate_design/unified_property_scoring_result_20260815/unified_single_mutant_property_evidence.csv"
OLD_ANTIFOLD = ROOT / "docs/result_artifacts/candidate_design/unified_single_mutant_antifold_20260815/unified_single_mutant_antifold_evidence.csv"


def _csv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _real_plan():
    current = _csv(V2 / "nb252_allowed_single_mutants.csv")
    first = current[0]
    index = int(first["reported_sequence_index_1based"])
    parent = first["sequence"][: index - 1] + first["wt_residue"] + first["sequence"][index:]
    result = build_reuse_plan(current, _json(PREFLIGHT), parent, _csv(OLD_PROPERTY), _csv(OLD_ANTIFOLD))
    return current, parent, result


def test_real_v2_reuse_plan_has_exact_complete_split_and_validation_panel():
    _, _, (audit, validation, expected, targets, completion, gate) = _real_plan()
    assert gate["status"] == "pass"
    assert len(audit) == 847 and len(validation) == len(expected) == 13
    assert len(targets) == 12 and len(completion) == 127
    assert Counter(row["netsolp_reuse_status"] for row in audit) == {
        "reused_pending_repeat_validation": 721,
        "requires_new_score": 126,
    }
    assert Counter(row["antifold_evaluation_scope"] for row in audit) == {"three_views": 721, "af3_only": 126}
    assert Counter(int(row["reported_sequence_index_1based"]) for row in completion[1:]) == {
        11: 18, 14: 18, 24: 18, 26: 18, 27: 18, 28: 18, 29: 18,
    }
    assert any(row["candidate_id"] == "Nb252_expr_seq005_Q5V" for row in validation)
    assert all(row["exact_sequence_match"] == "True" or row["exact_sequence_match"] is True for row in audit)


def test_repeat_comparison_accepts_recorded_precision_and_blocks_a_changed_score():
    _, _, (_, samples, expected, _, _, _) = _real_plan()
    expected_by_id = {row["score_id"]: row for row in expected}
    net = [{
        "sample_uid": row["score_id"], "sequence_raw": row["sequence_raw"],
        "predicted_usability": expected_by_id[row["score_id"]]["netsolp_predicted_usability"],
        "predicted_solubility": expected_by_id[row["score_id"]]["netsolp_predicted_solubility"],
    } for row in samples]
    melt = [{
        "sample_uid": row["score_id"], "sequence_raw": row["sequence_raw"],
        "nanomelt_predicted_apparent_tm_c": expected_by_id[row["score_id"]]["nanomelt_predicted_apparent_tm_c"],
        "scored_length_aa": 126, "trimmed_c_terminal": "GS",
    } for row in samples]
    anti = []
    for row in expected[1:]:
        candidate = {"candidate_id": row["candidate_id"]}
        for view in ("experimental_vhh_only", "experimental_complex_context", "af3_vhh_only"):
            for suffix in ("evaluation_status", "wt_log_probability", "mutant_log_probability", "delta_log_probability", "perplexity"):
                candidate[f"{view}_{suffix}"] = row[f"{view}_{suffix}"]
        anti.append(candidate)
    comparisons, gate = compare_repeat_scores(samples, expected, net, melt, anti)
    assert gate["status"] == "pass" and comparisons
    net[1]["predicted_solubility"] = float(net[1]["predicted_solubility"]) + 0.01
    _, blocked = compare_repeat_scores(samples, expected, net, melt, anti)
    assert blocked["status"] == "blocked" and blocked["failure_count"] == 1


def test_complete_matrix_reuses_721_and_fills_only_126():
    current, parent, (audit, _, expected, _, completion, _) = _real_plan()
    wt = expected[0]
    net = []
    melt = []
    for number, sample in enumerate(completion):
        net.append({
            "sample_uid": sample["score_id"], "sequence_raw": sample["sequence_raw"],
            "predicted_usability": float(wt["netsolp_predicted_usability"]) + number / 100000,
            "predicted_solubility": float(wt["netsolp_predicted_solubility"]) + number / 100000,
        })
        melt.append({
            "sample_uid": sample["score_id"], "sequence_raw": sample["sequence_raw"],
            "nanomelt_predicted_apparent_tm_c": float(wt["nanomelt_predicted_apparent_tm_c"]) + number / 100,
            "scored_length_aa": 126, "trimmed_c_terminal": "GS",
        })
    matrix, gate = build_complete_score_matrix(
        current, parent, audit, _csv(OLD_PROPERTY), _csv(OLD_ANTIFOLD), net, melt,
        {"status": "pass", "release": "legacy_scores_validated_for_exact_reuse"},
    )
    assert gate["status"] == "pass" and len(matrix) == 847
    assert gate["property_score_source_counts"] == {
        "new_missing_position_score": 126,
        "reused_20260815_repeat_validated": 721,
    }
    assert gate["antifold_scope_counts"] == {"af3_only": 126, "three_views": 721}
    assert all(not row["candidate_selection_performed"] for row in matrix)


def test_completion_slurm_is_sequential_gated_and_has_no_resume_or_array():
    slurm = (ROOT / "scripts/candidate_design/submit_expression_property_completion_v2.slurm").read_text(encoding="utf-8")
    wrapper = (ROOT / "scripts/candidate_design/submit_expression_property_completion_v2.sh").read_text(encoding="utf-8")
    assert "#SBATCH --partition=batch" in slurm and "#SBATCH --gres=gpu:1" in slurm
    assert "logs/expression_property_completion_v2/run-%j.log" in slurm
    assert slurm.index("validate_expression_property_reuse.py") < slurm.index("missing_property_samples.csv")
    assert "--array" not in slurm and "resume" not in slurm.lower()
    assert "sbatch --parsable" in wrapper
