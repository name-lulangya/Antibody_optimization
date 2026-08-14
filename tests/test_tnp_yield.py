from __future__ import annotations

import csv
from pathlib import Path

import pytest

from antibody_optimization.tnp_yield import TNPYieldError, analyze_tnp_associations, build_tnp_validation_inputs, normalize_tnp_result
from antibody_optimization.tnp_yield_plot import render_tnp_yield_figure


ROOT = Path(__file__).resolve().parents[1]


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _real_samples() -> list[dict[str, object]]:
    return build_tnp_validation_inputs(
        _csv(ROOT / "docs/result_artifacts/nb_expression/nb_expression_records.csv"),
        _csv(ROOT / "docs/result_artifacts/input_baseline/sequence/sequence_numbering_review.csv"),
        _csv(ROOT / "docs/result_artifacts/input_baseline/sequence/sequence_numbering_positions.csv"),
    )["sample_rows"]


def _result(uid: str, psh: float) -> dict[str, object]:
    return {uid: {"name": uid, "Total CDR Length": 35, "CDR3 Length": 20, "CDR3 Compactness": 1.2, "PSH": psh, "PPC": 0.1, "PNC": 0.2, "Flags": {"L": "green", "L3": "green", "C": "green", "PSH": "amber", "PPC": "green", "PNC": "green"}}}


def test_real_47_tnp_plan_reuses_frozen_clusters() -> None:
    samples = _real_samples()
    assert len(samples) == 47
    assert len({row["sequence_cluster_90"] for row in samples}) == 40
    assert next(row for row in samples if row["sample_uid"] == "WCC__4-28")["numbering_status"] == "failed"


def test_normalize_tnp_result_records_unique_terminal_trim_and_flags() -> None:
    sample = _real_samples()[0]
    sequence = str(sample["sequence_raw"])
    row = normalize_tnp_result(sample, _result(str(sample["sample_uid"]), 141.2), modelled_sequence=sequence[:-2], elapsed_seconds=2.5)
    assert row["trimmed_n_terminal"] == ""
    assert row["trimmed_c_terminal"] == sequence[-2:]
    assert row["tnp_flag_psh"] == "amber"
    with pytest.raises(TNPYieldError, match="unique subsequence"):
        normalize_tnp_result(sample, _result(str(sample["sample_uid"]), 141.2), modelled_sequence="AAAA", elapsed_seconds=2.5)


def test_tnp_analysis_keeps_psh_primary_compares_netsolp_and_renders(tmp_path: Path) -> None:
    samples = _real_samples()
    scores, netsolp = [], []
    for index, sample in enumerate(samples):
        psh = 200.0 - float(sample["numeric_yield_value"]) if sample["observation_semantics"] == "individual_approximate" else 150.0 - 10.0 * int(sample["llj_ordinal_level"])
        scores.append(normalize_tnp_result(sample, _result(str(sample["sample_uid"]), psh), modelled_sequence=str(sample["sequence_raw"]), elapsed_seconds=1.0))
        netsolp.append({"sample_uid": sample["sample_uid"], "predicted_usability": 0.2 + index / 100.0})
    result = analyze_tnp_associations(samples, scores, netsolp)
    assert result["primary"]["feature"] == "tnp_psh"
    assert result["coverage"] == {"planned": 47, "passed": 47, "numeric_passed": 31, "llj_passed": 16}
    assert [row["model"] for row in result["cv_rows"]] == ["Provider only", "NetSolP U", "TNP PSH", "NetSolP U + TNP PSH"]
    assert all("pooled_spearman_bh_fdr" in row for row in result["metric_rows"])
    png, svg = tmp_path / "tnp.png", tmp_path / "tnp.svg"
    render_tnp_yield_figure(result["sample_rows"], result["metric_rows"], result["cv_rows"], png_path=png, svg_path=svg)
    assert png.stat().st_size > 1000
    assert svg.stat().st_size > 1000
