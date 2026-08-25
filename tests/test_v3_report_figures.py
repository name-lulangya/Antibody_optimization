import csv
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from antibody_optimization.v3_report_figures import (
    load_v3_report_figure_data,
    render_v3_report_figures,
)


def test_v3_report_figure_inputs_are_bound_to_released_v3_panel() -> None:
    data = load_v3_report_figure_data(PROJECT_ROOT)

    selected_parents = [
        row for row in data.parent_rows if row["v3_parent_selection_status"] == "selected"
    ]
    selected_doubles = [
        row for row in data.double_rows if row["final_double_selection_status"] == "selected"
    ]

    assert len(data.parent_rows) == 31
    assert len(selected_parents) == 15
    assert len(data.double_rows) == 102
    assert len(selected_doubles) == 15
    assert data.single_gate["candidate_count"] - data.single_gate["antifold_veto_count"] == 696
    assert data.single_contract["antifold_role"] == "negative_veto_only_no_positive_selection_credit"
    assert data.final_manifest["facts"]["selected_three_metric_positive_count"] == 6
    assert data.final_manifest["facts"]["selected_two_metric_positive_count"] == 9


def test_v3_report_figures_render_png_and_svg(tmp_path: Path) -> None:
    outputs = render_v3_report_figures(PROJECT_ROOT, tmp_path)

    assert set(outputs) == {
        "single_selection_flow",
        "parent15_property_heatmap",
        "double_selection_flow",
        "double15_property_heatmap",
        "tool_validation_summary",
        "source_data_csv",
    }
    for key in (
        "single_selection_flow",
        "parent15_property_heatmap",
        "double_selection_flow",
        "double15_property_heatmap",
        "tool_validation_summary",
    ):
        png_path, svg_path = outputs[key]
        assert png_path.is_file() and png_path.stat().st_size > 10_000
        assert svg_path.is_file() and svg_path.stat().st_size > 10_000
    source_csv = outputs["source_data_csv"]
    assert source_csv.is_file() and source_csv.stat().st_size > 1_000
    with source_csv.open(encoding="utf-8", newline="") as handle:
        source_rows = list(csv.DictReader(handle))
    antifold_role = [
        row
        for row in source_rows
        if row["item_id"] == "AntiFold" and row["metric"] == "V3 role"
    ]
    assert antifold_role[0]["value"] == (
        "AntiFold only excludes high-risk substitutions, never proposes, "
        "rewards, or ranks candidates."
    )
