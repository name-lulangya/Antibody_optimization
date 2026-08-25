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

from antibody_optimization.stable_words import parse_stable_words
from antibody_optimization.v3_double_mutants import (
    V3DoubleMutantError,
    build_v3_double_mutant_space,
    build_v3_score_samples,
    merge_v3_property_scores,
)


PARENT_DIR = (
    ROOT
    / "docs/result_artifacts/candidate_design"
    / "v3_parent_single_selection_20260825"
)
PLAN_DIR = (
    ROOT
    / "docs/result_artifacts/candidate_design"
    / "v3_double_mutant_plan_20260825"
)
MAPPING = (
    ROOT
    / "docs/result_artifacts/input_baseline/structure_released_20260810"
    / "nb252_sequence_structure_mapping.csv"
)
EXPERIMENTAL = ROOT / "data/structures/cxs_exports/NK2R-252__native.cif"
AF3 = ROOT / "data/structures/cxs_exports/fold_2r_252_nomg_model_0__native.cif"
WORDS = ROOT / "data/Stable_word_SS_3D_1336 (1).txt"


def _csv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


@lru_cache(maxsize=1)
def _result():
    words = parse_stable_words(
        WORDS.read_text(encoding="utf-8-sig").splitlines(keepends=True)
    )
    return build_v3_double_mutant_space(
        _csv(PARENT_DIR / "v3_parent_single_selected15.csv"),
        _csv(PARENT_DIR / "v3_parent_single_selection_audit.csv"),
        words,
        _csv(MAPPING),
        EXPERIMENTAL,
        AF3,
    )


def test_real_parent15_enumerates_exact_v3_pair_space():
    result = _result()
    assert len(result["parents"]) == 15
    assert result["facts"]["parent_unique_position_count"] == 12
    assert result["facts"]["theoretical_unordered_pair_count"] == 105
    assert len(result["invalid_pairs"]) == 3
    assert len(result["candidates"]) == 102
    assert len(build_v3_score_samples(result["parent_sequence"], result["candidates"])) == 103
    invalid = {
        frozenset((row["mutation_a"], row["mutation_b"]))
        for row in result["invalid_pairs"]
    }
    assert invalid == {
        frozenset(("L11Y", "L11M")),
        frozenset(("F30S", "F30N")),
        frozenset(("K75A", "K75E")),
    }


def test_real_v3_double_sequences_are_unique_exact_constructs():
    result = _result()
    parent = result["parent_sequence"]
    candidates = result["candidates"]
    assert len({row["double_candidate_id"] for row in candidates}) == 102
    assert len({row["sequence"] for row in candidates}) == 102
    assert sum(row["contains_t99f_stable_word_exploration_parent"] for row in candidates) == 14
    for row in candidates:
        sequence = row["sequence"]
        assert len(sequence) == 128
        assert sequence.endswith("SSGS")
        assert sequence[21] == "C" and sequence[94] == "C" and sequence.count("C") == 2
        differences = [
            index
            for index, (wt, mutant) in enumerate(zip(parent, sequence, strict=True), 1)
            if wt != mutant
        ]
        assert differences == sorted(
            (row["position_a_reported_1based"], row["position_b_reported_1based"])
        )
        assert row["candidate_prefiltering_applied"] is False


def test_real_v3_double_stable_words_and_sequence_risks_are_recomputed():
    rows = _result()["candidates"]
    stable = Counter(row["stable_word_effect"] for row in rows)
    assert stable == {"unchanged": 87, "gain_only": 15}
    gains = [row for row in rows if row["stable_word_effect"] == "gain_only"]
    assert sum(row["contains_t99f_stable_word_exploration_parent"] for row in gains) == 14
    assert [row["mutation_set"] for row in gains if not row["contains_t99f_stable_word_exploration_parent"]] == [
        "N76G;K75E"
    ]
    assert all(row["hard_sequence_risk_count"] == 0 for row in rows)
    assert Counter(row["soft_sequence_risk_count"] for row in rows) == {0: 78, 1: 23, 2: 1}
    flags = Counter()
    for row in rows:
        flags.update(value for value in row["soft_sequence_risk_flags"].split("|") if value)
    assert flags == {"more_M_or_W": 13, "new_deamidation_motif": 12}


def test_antifold_is_constituent_only_and_never_added():
    for row in _result()["candidates"]:
        assert row["antifold_component_a_veto_status"] == "pass"
        assert row["antifold_component_b_veto_status"] == "pass"
        assert row["antifold_constituent_gate"] == "pass"
        assert row["antifold_double_mutant_scored"] is False
        assert row["antifold_component_values_combined"] is False
        assert row["antifold_double_mutant_score"] == ""
        assert not any("additive" in key or "sum" in key for key in row)


def test_real_pair_geometry_preserves_experimental_and_af3_evidence_boundaries():
    rows = _result()["candidates"]
    assert Counter(row["pair_structure_distance_source"] for row in rows) == {
        "experimental_complex_vhh": 64,
        "af3_vhh_only_due_missing_experimental_coordinate": 38,
    }
    assert Counter(row["pair_spatial_class"] for row in rows) == {
        "direct_local_neighborhood_under_4p5A": 6,
        "nearby_ca_under_10A": 6,
        "spatially_separated_ca_at_least_10A": 90,
    }
    assert Counter(row["machine_structure_triage_status"] for row in rows) == {
        "detailed_review_triggered": 53,
        "routine_context_recorded": 49,
    }
    missing = [row for row in rows if row["pair_experimental_coordinate_status"] != "both_observed"]
    assert len(missing) == 38
    assert all(row["experimental_pair_ca_distance_a"] == "" for row in missing)
    assert all(float(row["af3_pair_ca_distance_a"]) > 0 for row in missing)


def test_complete_score_merge_uses_same_wt_and_labels_model_nonadditivity():
    result = _result()
    samples = build_v3_score_samples(result["parent_sequence"], result["candidates"])
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
                "nanomelt_predicted_apparent_tm_c": 65.0 + index / 100,
                "scored_length_aa": 126,
                "trimmed_c_terminal": "GS",
                "scoring_status": "pass",
            }
        )
    merged = merge_v3_property_scores(result["candidates"], net_rows, melt_rows)
    assert len(merged) == 102
    assert merged[0]["netsolp_u_delta_vs_wt"] == pytest.approx(0.0001)
    assert merged[0]["nanomelt_tm_c_delta_vs_wt"] == pytest.approx(0.01)
    assert "netsolp_u_model_nonadditivity_residual" in merged[0]
    assert merged[0]["model_nonadditivity_interpretation"] == (
        "predictor_output_residual_not_physical_epistasis"
    )
    assert merged[0]["final_double_selection_status"] == "not_performed"
    with pytest.raises(V3DoubleMutantError):
        merge_v3_property_scores(result["candidates"], net_rows[:-1], melt_rows)


def test_plan_cli_preserves_parent_artifacts_and_writes_complete_outputs():
    selected_path = PARENT_DIR / "v3_parent_single_selected15.csv"
    selected_bytes = selected_path.read_bytes()
    with tempfile.TemporaryDirectory(prefix=".test-v3-double-plan-", dir=ROOT) as temp:
        output = Path(temp) / "plan"
        summary = Path(temp) / "run_summary.json"
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "scripts/candidate_design/build_v3_double_mutant_plan.py"),
                "--output-dir",
                str(output),
                "--run-summary",
                str(summary),
                "--generated-at",
                "2026-08-25T18:30:00+08:00",
            ],
            cwd=ROOT,
            check=True,
        )
        manifest = json.loads(
            (output / "v3_double_mutant_plan_manifest.json").read_text(encoding="utf-8")
        )
        assert manifest["gate"]["v3_double_mutant_plan"] == "pass"
        assert manifest["facts"]["valid_double_mutant_count"] == 102
        assert len(_csv(output / "v3_double_mutant_score_samples103.csv")) == 103
        assert (output / "v3_double_mutant_plan_overview.png").stat().st_size > 1000
    assert selected_path.read_bytes() == selected_bytes


def test_v3_remote_slurm_contract_is_sequential_expression_only():
    slurm = (
        ROOT / "scripts/candidate_design/submit_v3_double_mutant_scan.slurm"
    ).read_text(encoding="utf-8")
    wrapper = (
        ROOT / "scripts/candidate_design/run_v3_double_mutant_scan.sh"
    ).read_text(encoding="utf-8")
    assert "#SBATCH --partition=batch" in slurm
    assert "#SBATCH --gres=gpu:1" in slurm
    assert "#SBATCH --cpus-per-task=12" in slurm
    assert "#SBATCH --time=02:00:00" in slurm
    assert "logs/v3_double_mutant_scan/run-%j.log" in slurm
    assert "#SBATCH --array" not in slurm
    assert "#SBATCH --mem" not in slurm
    assert "pyrosetta" not in slurm.lower()
    assert "tnp" not in slurm.lower()
    assert "score_v3_double_mutant_properties.py" in slurm
    assert "finalize_v3_double_mutant_matrix.py" in slurm
    assert "sbatch scripts/candidate_design/submit_v3_double_mutant_scan.slurm" in wrapper


def test_v3_finalizer_accepts_complete_test_only_tool_outputs():
    samples = _csv(PLAN_DIR / "v3_double_mutant_score_samples103.csv")
    with tempfile.TemporaryDirectory(prefix=".test-v3-double-final-", dir=ROOT) as temp:
        work = Path(temp)
        net_dir = work / "netsolp"
        melt_dir = work / "nanomelt"
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
                    "aligned_sequence": "",
                    "scored_ungapped_sequence": sample["sequence_raw"][:-2],
                    "scored_length_aa": 126,
                    "trimmed_n_terminal": "",
                    "trimmed_c_terminal": "GS",
                    "nanomelt_predicted_apparent_tm_c": 65.0 + index / 100,
                    "scoring_status": "pass",
                    "not_scored_reason": "",
                }
            )
        _write_csv(net_dir / "netsolp_sample_scores.csv", net_rows)
        _write_csv(melt_dir / "nanomelt_sample_scores.csv", melt_rows)
        (net_dir / "netsolp_model_run.json").write_text(
            '{"status":"pass","tool":"netsolp","pass_count":103}\n',
            encoding="utf-8",
        )
        (melt_dir / "nanomelt_model_run.json").write_text(
            '{"status":"pass","tool":"nanomelt","pass_count":103}\n',
            encoding="utf-8",
        )
        output = work / "output"
        summary = work / "run_summary.json"
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "scripts/candidate_design/finalize_v3_double_mutant_matrix.py"),
                "--plan-dir",
                str(PLAN_DIR),
                "--netsolp-score-dir",
                str(net_dir),
                "--nanomelt-score-dir",
                str(melt_dir),
                "--output-dir",
                str(output),
                "--run-summary",
                str(summary),
                "--generated-at",
                "2026-08-25T19:00:00+08:00",
            ],
            cwd=ROOT,
            check=True,
        )
        manifest = json.loads(
            (output / "v3_double_mutant_property_matrix_manifest.json").read_text(
                encoding="utf-8"
            )
        )
        assert manifest["gate"]["v3_double_complete_property_matrix"] == "pass"
        assert manifest["gate"]["final_15_double_mutant_selection"] == "not_performed"
        assert len(_csv(output / "v3_double_mutant_property_matrix102.csv")) == 102
        assert len(_csv(output / "v3_double_mutant_property_plot_data.csv")) == 306
        assert (output / "v3_double_mutant_property_overview.png").stat().st_size > 1000
