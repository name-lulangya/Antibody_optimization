from __future__ import annotations

import ast
import csv
import math
from pathlib import Path

import pytest

from antibody_optimization.plm_sol_yield import (
    PLMSolYieldError,
    analyze_plm_sol_associations,
    normalize_plm_sol_scores,
)
from antibody_optimization.plm_sol_yield_plot import (
    render_plm_sol_fixed5_figure,
    render_plm_sol_yield_figure,
)


ROOT = Path(__file__).resolve().parents[1]
RP3_PLAN = ROOT / "docs/result_artifacts/candidate_design/rp3net_yield_validation_plan_20260818"
NET_RESULTS = ROOT / "docs/result_artifacts/candidate_design/netsolp_yield_validation_result_20260814"
RP3_RESULTS = ROOT / "docs/result_artifacts/candidate_design/rp3net_yield_validation_result_20260818"


def test_plm_sol_tool_scripts_avoid_python311_path_write_text_api() -> None:
    scripts = (
        ROOT / "scripts/candidate_design/score_plm_sol_embeddings.py",
        ROOT / "scripts/candidate_design/score_plm_sol_classifier.py",
    )
    for script in scripts:
        tree = ast.parse(script.read_text(encoding="utf-8"), filename=str(script), feature_version=(3, 8))
        unsupported = [
            node.lineno
            for node in ast.walk(tree)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "write_text"
            and any(keyword.arg == "newline" for keyword in node.keywords)
        ]
        assert unsupported == [], f"Python 3.8-incompatible Path.write_text(newline=...) at {script}:{unsupported}"


def test_normalize_plm_sol_scores_maps_hashes_by_exact_sequence() -> None:
    samples = _csv(RP3_PLAN / "rp3net_validation_samples.csv")
    raw = [
        {
            "Unnamed: 0": index,
            "protein_ID": f"hash_{index:03d}",
            "sequence": row["sequence_raw"],
            "predict_result": 0.2 + index / 100,
        }
        for index, row in enumerate(reversed(samples))
    ]
    result = normalize_plm_sol_scores(samples, raw)
    assert len(result) == 47
    assert [row["sample_uid"] for row in result] == [row["sample_uid"] for row in samples]
    assert result[0]["embedding_key"] == "hash_046"
    assert result[0]["sequence_length_aa"] == len(samples[0]["sequence_raw"])
    assert "Unnamed: 0" not in result[0]


def test_normalize_plm_sol_scores_rejects_sequence_mismatch() -> None:
    samples = _csv(RP3_PLAN / "rp3net_validation_samples.csv")
    raw = [
        {"protein_ID": f"hash_{index}", "sequence": row["sequence_raw"], "predict_result": 0.5}
        for index, row in enumerate(samples)
    ]
    raw[0]["sequence"] = raw[0]["sequence"][:-1]
    with pytest.raises(PLMSolYieldError, match="do not match"):
        normalize_plm_sol_scores(samples, raw)


def test_plm_sol_analysis_preserves_all_semantics_and_outputs(tmp_path: Path) -> None:
    samples = _csv(RP3_PLAN / "rp3net_validation_samples.csv")
    netsolp = _csv(NET_RESULTS / "netsolp_yield_sample_evidence.csv")
    rp3net = _csv(RP3_RESULTS / "rp3net_yield_sample_evidence.csv")
    net_by_id = {row["sample_uid"]: row for row in netsolp}
    scores = [
        {
            "sample_uid": row["sample_uid"],
            "sequence_raw": row["sequence_raw"],
            "embedding_key": f"hash_{index:03d}",
            "plm_sol_solubility_score": net_by_id[row["sample_uid"]]["predicted_solubility"],
            "scoring_status": "pass",
        }
        for index, row in enumerate(samples)
    ]
    result = analyze_plm_sol_associations(samples, scores, netsolp, rp3net)
    assert len(result["sample_rows"]) == 47
    assert sum(row["observation_semantics"] == "individual_approximate" for row in result["sample_rows"]) == 31
    assert len(result["classification_rows"]) == 2
    assert len(result["classification_prediction_rows"]) == 62
    assert len(result["fixed5_metric_rows"]) == 3
    assert len(result["comparison_rows"]) == 5
    assert result["comparison_rows"][0]["comparison"] == "NetSolP S"
    assert float(result["comparison_rows"][0]["all_47_spearman"]) == pytest.approx(1.0)
    assert math.isfinite(float(result["comparison_rows"][3]["partial_spearman_provider_length_netsolp_s"]))
    render_plm_sol_yield_figure(
        result["sample_rows"], result["metric_rows"], result["classification_rows"], result["comparison_rows"],
        png_path=tmp_path / "main.png", svg_path=tmp_path / "main.svg",
    )
    render_plm_sol_fixed5_figure(
        result["sample_rows"], result["fixed5_metric_rows"],
        png_path=tmp_path / "fixed.png", svg_path=tmp_path / "fixed.svg",
    )
    assert all((tmp_path / name).stat().st_size > 0 for name in ("main.png", "main.svg", "fixed.png", "fixed.svg"))


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))
