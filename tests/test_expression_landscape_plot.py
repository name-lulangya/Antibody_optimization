from __future__ import annotations

import csv
from pathlib import Path

import pytest

from antibody_optimization.expression_landscape_plot import (
    ExpressionLandscapeError,
    build_expression_landscape_rows,
    render_expression_landscape,
    render_expression_scatter,
)


ROOT = Path(__file__).resolve().parents[1]
MATRIX = ROOT / "docs/result_artifacts/candidate_design/stable_word_single_mutant_v1_20260819/expression_single_mutant_property_stable_word_matrix.csv"


def _rows() -> list[dict[str, str]]:
    with MATRIX.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def test_real_landscape_contract_counts() -> None:
    rows, facts = build_expression_landscape_rows(_rows())
    assert len(rows) == 847
    assert facts == {
        "candidate_count": 847,
        "reported_position_count": 48,
        "stable_word_gain_candidate_count": 22,
        "experimental_complex_antifold_pass_count": 721,
        "experimental_complex_antifold_not_evaluable_count": 126,
        "af3_antifold_fallback_count": 126,
        "candidate_selection_performed": False,
    }
    assert all(row["antifold_landscape_delta_log_probability"] != "" for row in rows)
    assert sum(row["antifold_landscape_source"].startswith("af3_vhh_only") for row in rows) == 126


def test_antifold_status_value_mismatch_is_rejected() -> None:
    rows = _rows()
    rows[0] = dict(rows[0])
    rows[0]["experimental_complex_context_evaluation_status"] = "not_evaluable"
    with pytest.raises(ExpressionLandscapeError, match="value/status mismatch"):
        build_expression_landscape_rows(rows)


def test_real_landscape_renders_png_and_svg(tmp_path: Path) -> None:
    rows, _ = build_expression_landscape_rows(_rows())
    png, svg = tmp_path / "landscape.png", tmp_path / "landscape.svg"
    scatter_png, scatter_svg = tmp_path / "scatter.png", tmp_path / "scatter.svg"
    render_expression_landscape(rows, png_path=png, svg_path=svg)
    render_expression_scatter(rows, png_path=scatter_png, svg_path=scatter_svg)
    assert png.stat().st_size > 100_000
    assert svg.stat().st_size > 10_000
    assert scatter_png.stat().st_size > 100_000
    assert scatter_svg.stat().st_size > 10_000
