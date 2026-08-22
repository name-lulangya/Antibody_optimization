from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from antibody_optimization.expression_double_mutants import (
    WT_SCORE_ID,
    build_double_mutant_space,
    build_score_samples,
    merge_property_scores,
)
from antibody_optimization.stable_words import compare_stable_word_occurrences


PARENT_DIR = ROOT / "docs/result_artifacts/candidate_design/expression_single_mutant_parent19_20260822"
WORDS = ROOT / "docs/result_artifacts/candidate_design/stable_word_single_mutant_v1_20260819/stable_word_library.csv"


def _csv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _result():
    parents = _csv(PARENT_DIR / "expression_single_mutant_parent19.csv")
    gate = json.loads((PARENT_DIR / "expression_single_mutant_parent19_gate.json").read_text(encoding="utf-8"))
    words = [row["stable_word"] for row in _csv(WORDS)]
    return parents, build_double_mutant_space(parents, gate, words)


def test_real_parent19_enumerates_complete_distinct_position_space():
    parents, result = _result()
    candidates = result["candidates"]
    assert len(candidates) == 162
    assert len(result["invalid_pairs"]) == 9
    assert len(build_score_samples(result["parent_sequence"], candidates)) == 163
    assert all(row["position_a_reported_1based"] != row["position_b_reported_1based"] for row in candidates)
    assert all(row["sequence"].endswith("SSGS") for row in candidates)
    assert all(row["sequence"].count("C") == result["parent_sequence"].count("C") for row in candidates)
    assert result["facts"]["antifold_same_view_additive_evaluable_count"] == 102
    assert result["facts"]["stable_word_gain_candidate_count"] == 19
    assert result["facts"]["hard_sequence_risk_count"] == 0
    assert len(parents) == 19


def test_stable_word_comparison_recomputes_jointly_created_window():
    comparison = compare_stable_word_occurrences("DAD", "AGA", ("sss",))
    assert comparison["created_stable_word_occurrence_count"] == 1
    assert comparison["lost_stable_word_occurrence_count"] == 0
    assert comparison["stable_word_effect"] == "gain_only"


def test_complete_property_merge_calculates_deltas_and_nonadditivity():
    parents, result = _result()
    samples = build_score_samples(result["parent_sequence"], result["candidates"])
    net = []
    melt = []
    for index, sample in enumerate(samples):
        net.append(
            {
                "sample_uid": sample["sample_uid"],
                "sequence_raw": sample["sequence_raw"],
                "predicted_usability": 0.4 + index / 10000,
                "predicted_solubility": 0.5 + index / 10000,
                "scoring_status": "pass",
            }
        )
        melt.append(
            {
                "sample_uid": sample["sample_uid"],
                "sequence_raw": sample["sequence_raw"],
                "nanomelt_predicted_apparent_tm_c": 65 + index / 100,
                "scored_length_aa": 126,
                "trimmed_c_terminal": "GS",
                "scoring_status": "pass",
            }
        )
    merged = merge_property_scores(result["candidates"], parents, net, melt)
    assert len(merged) == 162
    assert merged[0]["netsolp_u_delta_vs_wt"] == pytest.approx(0.0001)
    assert merged[0]["nanomelt_tm_c_delta_vs_wt"] == pytest.approx(0.01)
    assert "netsolp_u_interaction_residual" in merged[0]
    assert merged[0]["candidate_selection_performed"] is False


def test_plan_cli_writes_new_artifacts_without_touching_parent():
    parent_bytes = (PARENT_DIR / "expression_single_mutant_parent19.csv").read_bytes()
    with tempfile.TemporaryDirectory(prefix=".test-double-plan-", dir=ROOT) as temp:
        output = Path(temp) / "plan"
        summary = Path(temp) / "summary.json"
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "scripts/candidate_design/build_expression_double_mutant_plan.py"),
                "--parent19-dir",
                str(PARENT_DIR),
                "--stable-word-library",
                str(WORDS),
                "--output-dir",
                str(output),
                "--run-summary",
                str(summary),
                "--generated-at",
                "2026-08-22T18:00:00+08:00",
            ],
            cwd=ROOT,
            check=True,
        )
        gate = json.loads((output / "expression_double_mutant_plan_gate.json").read_text(encoding="utf-8"))
        assert gate["valid_double_mutant_count"] == 162
        assert len(_csv(output / "expression_double_mutant_score_samples.csv")) == 163
    assert (PARENT_DIR / "expression_single_mutant_parent19.csv").read_bytes() == parent_bytes


def test_active_slurm_contract_is_expression_only_and_non_array():
    slurm = (ROOT / "scripts/candidate_design/submit_expression_double_mutant_scan.slurm").read_text(encoding="utf-8")
    wrapper = (ROOT / "scripts/candidate_design/run_expression_double_mutant_scan.sh").read_text(encoding="utf-8")
    assert "#SBATCH --partition=batch" in slurm
    assert "#SBATCH --gres=gpu:1" in slurm
    assert "#SBATCH --cpus-per-task=12" in slurm
    assert "logs/expression_double_mutant_scan/run-%j.log" in slurm
    assert "--array" not in slurm
    assert "pyrosetta" not in slurm.lower()
    assert "tnp" not in slurm.lower()
    assert "score_expression_double_mutant_properties.py" in slurm
    assert "sbatch scripts/candidate_design/submit_expression_double_mutant_scan.slurm" in wrapper


def test_finalizer_cli_accepts_complete_synthetic_tool_outputs():
    plan = ROOT / "docs/result_artifacts/candidate_design/expression_double_mutant_plan_20260822"
    samples = _csv(plan / "expression_double_mutant_score_samples.csv")
    with tempfile.TemporaryDirectory(prefix=".test-double-final-", dir=ROOT) as temp:
        work = Path(temp)
        net_dir, melt_dir = work / "net", work / "melt"
        net_dir.mkdir()
        melt_dir.mkdir()
        net_rows = []
        melt_rows = []
        for index, sample in enumerate(samples):
            net_rows.append(
                {
                    "sample_uid": sample["sample_uid"],
                    "sequence_raw": sample["sequence_raw"],
                    "predicted_usability": 0.4 + index / 10000,
                    "predicted_solubility": 0.5 + index / 10000,
                    "scoring_status": "pass",
                }
            )
            melt_rows.append(
                {
                    "sample_uid": sample["sample_uid"],
                    "sequence_raw": sample["sequence_raw"],
                    "nanomelt_predicted_apparent_tm_c": 65 + index / 100,
                    "scored_length_aa": 126,
                    "trimmed_c_terminal": "GS",
                    "scoring_status": "pass",
                }
            )
        for path, rows in (
            (net_dir / "netsolp_sample_scores.csv", net_rows),
            (melt_dir / "nanomelt_sample_scores.csv", melt_rows),
        ):
            with path.open("w", encoding="utf-8-sig", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
        (net_dir / "netsolp_model_run.json").write_text('{"status":"pass"}\n', encoding="utf-8")
        (melt_dir / "nanomelt_model_run.json").write_text('{"status":"pass"}\n', encoding="utf-8")
        output, summary = work / "output", work / "summary.json"
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "scripts/candidate_design/finalize_expression_double_mutant_matrix.py"),
                "--plan-dir",
                str(plan),
                "--parent19-dir",
                str(PARENT_DIR),
                "--netsolp-score-dir",
                str(net_dir),
                "--nanomelt-score-dir",
                str(melt_dir),
                "--output-dir",
                str(output),
                "--run-summary",
                str(summary),
                "--generated-at",
                "2026-08-22T18:00:00+08:00",
            ],
            cwd=ROOT,
            check=True,
        )
        gate = json.loads((output / "expression_double_mutant_property_matrix_gate.json").read_text(encoding="utf-8"))
        assert gate["candidate_count"] == 162
        assert gate["final_11_double_mutants_selected"] is False
        assert (output / "expression_double_mutant_property_overview.png").stat().st_size > 1000
