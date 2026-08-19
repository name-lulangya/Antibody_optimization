from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import pytest

from antibody_optimization.rp3net_yield import (
    RP3NetYieldError,
    analyze_rp3net_associations,
    build_rp3net_validation_inputs,
    normalize_rp3net_scores,
)
from antibody_optimization.yield_classification import (
    fixed_yield_apparent_classification,
    fixed_yield_nested_classification,
    nested_yield_classification,
)


ROOT = Path(__file__).resolve().parents[1]


def test_rp3net_scoring_entry_import_does_not_require_analysis_packages() -> None:
    completed = subprocess.run(
        [
            sys.executable,
            "-S",
            str(ROOT / "scripts/candidate_design/score_rp3net_sequences.py"),
            "--help",
        ],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _real_samples() -> list[dict[str, object]]:
    return build_rp3net_validation_inputs(
        _csv(ROOT / "docs/result_artifacts/nb_expression/nb_expression_records.csv"),
        _csv(ROOT / "docs/result_artifacts/input_baseline/sequence/sequence_numbering_review.csv"),
        _csv(ROOT / "docs/result_artifacts/input_baseline/sequence/sequence_numbering_positions.csv"),
    )["sample_rows"]


def test_rp3net_plan_and_normalization_cover_all_47_sequences() -> None:
    samples = _real_samples()
    assert len(samples) == 47
    raw = [{"id": row["sample_uid"], "score": str((index + 1) / 50)} for index, row in enumerate(samples)]
    normalized = normalize_rp3net_scores(samples, raw)
    assert len(normalized) == 47
    assert normalized[0]["sequence_raw"] == samples[0]["sequence_raw"]
    raw[0]["score"] = "1.1"
    with pytest.raises(RP3NetYieldError, match=r"outside \[0, 1\]"):
        normalize_rp3net_scores(samples, raw)


def test_nested_classification_fits_thresholds_without_held_out_sample() -> None:
    rows = [
        {
            "sample_uid": f"x{index}", "provider_code": "LTT" if index < 6 else "WCC",
            "numeric_yield_value": float(index + 1), "score": float(index + 1) / 12,
            "sequence_cluster_90": f"c{index // 2}",
        }
        for index in range(12)
    ]
    result = nested_yield_classification(rows, "score", outer_scheme="leave_one_out")
    assert result["summary"]["n"] == 12
    assert len(result["prediction_rows"]) == 12
    held = next(row for row in result["prediction_rows"] if row["sample_uid"] == "x0")
    assert held["provider_training_yield_threshold"] == 4.0
    assert 0 <= result["summary"]["roc_auc"] <= 1


def test_fixed_yield_classification_keeps_outcome_cutoff_constant() -> None:
    rows = [
        {
            "sample_uid": f"x{index}", "provider_code": "LTT" if index < 6 else "WCC",
            "numeric_yield_value": float(index + 1), "score": float(index + 1) / 12,
            "sequence_cluster_90": f"c{index // 2}",
        }
        for index in range(12)
    ]
    nested = fixed_yield_nested_classification(
        rows, "score", outer_scheme="leave_one_cluster_out", yield_threshold=7.0,
    )
    assert nested["summary"]["positive_count"] == 6
    assert {row["fixed_yield_threshold"] for row in nested["prediction_rows"]} == {7.0}
    assert {row["observed_high_yield"] for row in nested["prediction_rows"] if row["numeric_yield_value"] == 7.0} == {1}
    apparent = fixed_yield_apparent_classification(rows, "score", yield_threshold=7.0)
    assert apparent["summary"]["threshold_fit_scope"] == "all_numeric_rows_apparent"
    assert apparent["summary"]["mcc"] == 1.0


def test_rp3net_analysis_preserves_llj_semantics_and_reports_both_outer_schemes() -> None:
    samples = _real_samples()
    raw = []
    for row in samples:
        score = (
            min(float(row["numeric_yield_value"]) / 30.0, 0.99)
            if row["observation_semantics"] == "individual_approximate"
            else int(row["llj_ordinal_level"]) / 4.0
        )
        raw.append({"id": row["sample_uid"], "score": score})
    result = analyze_rp3net_associations(samples, normalize_rp3net_scores(samples, raw))
    assert result["primary"]["llj_ordinal_n"] == 16
    assert {row["outer_scheme"] for row in result["classification_rows"]} == {
        "leave_one_out", "leave_one_cluster_out"
    }
    assert len(result["classification_prediction_rows"]) == 62
