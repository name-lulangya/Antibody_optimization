"""Build report-sized figures from the authoritative Nb252 V3 artifacts.

The module is deliberately presentation-only.  It validates and visualizes
the already released 15-single plus 15-double panel without recalculating or
re-ranking candidates.  NetSolP U, NetSolP S, and NanoMelt predicted Tm are
shown as separate ordinal evidence bands.  AntiFold has one role throughout:
negative risk exclusion of single substitutions; it never proposes, rewards,
or ranks a candidate, and no double-mutant AntiFold score is inferred.

Inputs are the tracked V3 CSV/JSON artifacts beneath ``project_root``.  The
return value of :func:`load_v3_report_figure_data` is a validated in-memory
snapshot for report rendering.  The module does not write report documents
and does not read historical V1/V2 result artifacts.
"""

from __future__ import annotations

import csv
import json
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(tempfile.gettempdir()) / "antibody_optimization_matplotlib"),
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap


REPORT_SINGLE_DIR = Path(
    "docs/result_artifacts/candidate_design/v3_parent_single_selection_20260825"
)
UPSTREAM_SINGLE_DIR = Path(
    "docs/result_artifacts/candidate_design/"
    "expression_single_mutant_selection_v3_20260825"
)
REPORT_FINAL_DIR = Path(
    "docs/result_artifacts/candidate_design/v3_final_15plus15_panel_20260825"
)
VALIDATION_PATHS = {
    "NetSolP": Path(
        "docs/result_artifacts/candidate_design/"
        "netsolp_yield_validation_result_20260814/netsolp_yield_associations.csv"
    ),
    "NanoMelt": Path(
        "docs/result_artifacts/candidate_design/"
        "nanomelt_yield_validation_result_20260815/nanomelt_yield_associations.csv"
    ),
    "RP3Net": Path(
        "docs/result_artifacts/candidate_design/"
        "rp3net_yield_validation_result_20260818/rp3net_yield_associations.csv"
    ),
    "PLM_Sol": Path(
        "docs/result_artifacts/candidate_design/"
        "plm_sol_yield_validation_result_20260819/plm_sol_yield_associations.csv"
    ),
}
FIXED5_PATH = Path(
    "docs/result_artifacts/candidate_design/"
    "fixed5mg_predictor_classification_20260819/fixed5mg_classification_metrics.csv"
)
PLM_FIXED5_PATH = Path(
    "docs/result_artifacts/candidate_design/"
    "plm_sol_yield_validation_result_20260819/plm_sol_fixed5mg_metrics.csv"
)

BAND_TO_GRADE = {
    "strong_adverse": -2,
    "moderate_adverse": -1,
    "weak_adverse": 0,
    "negligible": 0,
    "weak_favorable": 0,
    "moderate_favorable": 1,
    "strong_favorable": 2,
}

_HEATMAP_CMAP = ListedColormap(
    ["#9b2d30", "#d98b78", "#f3f2ef", "#7db6d6", "#194f7a"]
)
_HEATMAP_NORM = BoundaryNorm(
    [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5], _HEATMAP_CMAP.N
)


@dataclass(frozen=True)
class V3ReportFigureData:
    """Validated source rows and counts required by the four report figures."""

    single_gate: Mapping[str, Any]
    single_contract: Mapping[str, Any]
    parent_manifest: Mapping[str, Any]
    final_manifest: Mapping[str, Any]
    parent_rows: tuple[Mapping[str, str], ...]
    double_rows: tuple[Mapping[str, str], ...]
    validation_associations: tuple[Mapping[str, str], ...]
    fixed5_metrics: tuple[Mapping[str, str], ...]


def _read_json(path: Path) -> Mapping[str, Any]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def _read_csv(path: Path) -> tuple[Mapping[str, str], ...]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = tuple(dict(row) for row in csv.DictReader(handle))
    if not rows:
        raise ValueError(f"Expected non-empty CSV: {path}")
    return rows


def _as_int(value: object) -> int:
    return int(float(str(value)))


def _selected_parent_rows(
    data: V3ReportFigureData,
) -> list[Mapping[str, str]]:
    rows = [row for row in data.parent_rows if row["v3_parent_selection_status"] == "selected"]
    return sorted(rows, key=lambda row: _as_int(row["v3_parent_panel_order_not_efficacy_rank"]))


def _selected_double_rows(
    data: V3ReportFigureData,
) -> list[Mapping[str, str]]:
    rows = [row for row in data.double_rows if row["final_double_selection_status"] == "selected"]
    return sorted(rows, key=lambda row: _as_int(row["final_double_panel_order_not_efficacy_rank"]))


def load_v3_report_figure_data(project_root: Path) -> V3ReportFigureData:
    """Load and cross-check the released V3 report inputs.

    The checks protect figure semantics, not candidate selection.  In
    particular they bind the 31-row expert-review pool, selected 15 parents,
    102 reviewed doubles, selected 15 doubles, and the frozen AntiFold role.
    """

    project_root = Path(project_root)
    single_gate = _read_json(
        project_root / UPSTREAM_SINGLE_DIR / "expression_single_mutant_v3_gate.json"
    )
    single_contract = _read_json(
        project_root / UPSTREAM_SINGLE_DIR / "expression_single_mutant_v3_contract.json"
    )
    parent_manifest = _read_json(
        project_root / REPORT_SINGLE_DIR / "v3_parent_single_selection_manifest.json"
    )
    final_manifest = _read_json(
        project_root / REPORT_FINAL_DIR / "v3_final_panel_manifest.json"
    )
    parent_rows = _read_csv(
        project_root / REPORT_SINGLE_DIR / "v3_parent_single_selection_audit.csv"
    )
    double_rows = _read_csv(
        project_root / REPORT_FINAL_DIR / "v3_double_mutant_final_selection_audit102.csv"
    )
    validation_associations: list[Mapping[str, str]] = []
    for tool, relative_path in VALIDATION_PATHS.items():
        for source_row in _read_csv(project_root / relative_path):
            row = dict(source_row)
            row["tool"] = tool
            row["source_artifact"] = relative_path.as_posix()
            validation_associations.append(row)
    fixed5_metrics = [
        dict(row, source_artifact=FIXED5_PATH.as_posix())
        for row in _read_csv(project_root / FIXED5_PATH)
        if row["outer_scheme"] == "leave_one_cluster_out"
    ]
    for source_row in _read_csv(project_root / PLM_FIXED5_PATH):
        if source_row["outer_scheme"] == "leave_one_cluster_out":
            fixed5_metrics.append(
                dict(
                    source_row,
                    feature="plm_sol_solubility_score",
                    predictor="PLM_Sol",
                    source_artifact=PLM_FIXED5_PATH.as_posix(),
                )
            )
    data = V3ReportFigureData(
        single_gate=single_gate,
        single_contract=single_contract,
        parent_manifest=parent_manifest,
        final_manifest=final_manifest,
        parent_rows=parent_rows,
        double_rows=double_rows,
        validation_associations=tuple(validation_associations),
        fixed5_metrics=tuple(fixed5_metrics),
    )

    if single_gate.get("candidate_count") != 847:
        raise ValueError("V3 report expects the frozen 847-single search space")
    if single_gate.get("antifold_veto_count") != 151:
        raise ValueError("V3 report expects 151 AntiFold negative-veto exclusions")
    if single_gate.get("qualified_count") != 61 or single_gate.get("selected_count") != 30:
        raise ValueError("Unexpected upstream V3 single-selection counts")
    if single_contract.get("antifold_role") != "negative_veto_only_no_positive_selection_credit":
        raise ValueError("AntiFold role must remain negative exclusion only")
    if len(parent_rows) != 31 or len({row["candidate_id"] for row in parent_rows}) != 31:
        raise ValueError("Expected 31 unique parent-review rows")
    if len(_selected_parent_rows(data)) != 15:
        raise ValueError("Expected 15 selected parent singles")
    if len(double_rows) != 102 or len({row["double_candidate_id"] for row in double_rows}) != 102:
        raise ValueError("Expected 102 unique reviewed double mutants")
    if len(_selected_double_rows(data)) != 15:
        raise ValueError("Expected 15 selected double mutants")
    policy = final_manifest.get("selection_policy", {})
    if policy.get("antifold_role") != "constituent_negative_veto_only_no_double_score_no_positive_rank":
        raise ValueError("Double-mutant report must not imply an AntiFold score or positive rank")
    if any(str(row["antifold_double_mutant_scored"]).lower() != "false" for row in double_rows):
        raise ValueError("The released V3 double space must not contain AntiFold double scores")
    expected_association_features = {
        "predicted_usability",
        "predicted_solubility",
        "nanomelt_predicted_apparent_tm_c",
        "rp3net_expression_probability",
        "plm_sol_solubility_score",
    }
    available_association_features = {
        row["feature"] for row in validation_associations
    }
    if not expected_association_features <= available_association_features:
        raise ValueError("Missing one or more report validation association rows")
    if {row["predictor"] for row in fixed5_metrics} != {
        "NetSolP U",
        "NetSolP S",
        "RP3Net",
        "PLM_Sol",
    }:
        raise ValueError("Expected four fixed-5-mg/L display predictors")
    return data


def _configure_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [
                "Microsoft YaHei",
                "Noto Sans CJK SC",
                "SimHei",
                "Arial Unicode MS",
                "DejaVu Sans",
            ],
            "axes.unicode_minus": False,
            "axes.edgecolor": "#333333",
            "axes.labelcolor": "#222222",
            "text.color": "#222222",
            "xtick.color": "#333333",
            "ytick.color": "#333333",
        }
    )


def _save(fig: plt.Figure, stem: Path) -> tuple[Path, Path]:
    stem.parent.mkdir(parents=True, exist_ok=True)
    png_path = stem.with_suffix(".png")
    svg_path = stem.with_suffix(".svg")
    fig.savefig(png_path, dpi=400, bbox_inches="tight", facecolor="white")
    fig.savefig(svg_path, bbox_inches="tight", facecolor="white", metadata={"Date": None})
    plt.close(fig)
    svg_path.write_text(
        "\n".join(line.rstrip() for line in svg_path.read_text(encoding="utf-8").splitlines())
        + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return png_path, svg_path


def render_single_selection_flow(data: V3ReportFigureData, stem: Path) -> tuple[Path, Path]:
    """Render the exact 847→696→61→30→31→15 single-selection flow."""

    values = [847, 696, 61, 30, 31, 15]
    labels = [
        "完整约束单突空间",
        "AntiFold风险排除后",
        "性质规则合格",
        "上游短名单",
        "专家审查池\n（+T99F）",
        "父单突15条",
    ]
    fig, ax = plt.subplots(figsize=(7.35, 2.35))
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    x = np.linspace(0.075, 0.925, len(values))
    ax.plot(x, [0.58] * len(x), color="#9aa3aa", linewidth=1.4, zorder=1)
    for index, (xx, value, label) in enumerate(zip(x, values, labels)):
        color = "#194f7a" if index in {0, 5} else "#5b8fb5"
        ax.scatter(xx, 0.58, s=500, color=color, edgecolor="white", linewidth=1.2, zorder=3, clip_on=False)
        ax.text(xx, 0.58, str(value), ha="center", va="center", color="white", weight="bold", fontsize=9.6)
        ax.text(xx, 0.31, label, ha="center", va="top", fontsize=8.4, linespacing=1.25)
        if index < len(values) - 1:
            ax.annotate("", xy=(x[index + 1] - 0.035, 0.58), xytext=(xx + 0.035, 0.58),
                        arrowprops={"arrowstyle": "->", "color": "#59636b", "lw": 1.0})
    ax.text(
        0.5,
        0.02,
        "AntiFold仅用于负向风险排除，不提供正向加分、不排序、也不提议候选。",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=8.5,
        color="#333333",
    )
    ax.set_title("V3父单突筛选路径", loc="left", fontsize=12, weight="bold", pad=8)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.83, bottom=0.02)
    return _save(fig, stem)


def _heatmap(
    ax: plt.Axes,
    rows: Sequence[Mapping[str, str]],
    band_keys: Sequence[str],
    row_labels: Sequence[str],
    metric_labels: Sequence[str],
) -> plt.cm.ScalarMappable:
    matrix = np.asarray([[BAND_TO_GRADE[row[key]] for key in band_keys] for row in rows], dtype=float)
    image = ax.imshow(matrix, cmap=_HEATMAP_CMAP, norm=_HEATMAP_NORM, aspect="auto")
    ax.set_yticks(np.arange(len(row_labels)), row_labels, fontsize=8.4)
    ax.set_xticks(np.arange(len(metric_labels)), metric_labels, fontsize=9)
    ax.tick_params(length=0)
    for y in range(len(rows) + 1):
        ax.axhline(y - 0.5, color="white", lw=0.8)
    for x in range(len(metric_labels) + 1):
        ax.axvline(x - 0.5, color="white", lw=0.8)
    return image


def _add_band_colorbar(
    fig: plt.Figure,
    image: plt.cm.ScalarMappable,
    ax: plt.Axes | None = None,
    cax: plt.Axes | None = None,
    orientation: str = "vertical",
) -> None:
    colorbar = fig.colorbar(
        image,
        ax=ax,
        cax=cax,
        fraction=0.05,
        pad=0.035,
        ticks=[-2, -1, 0, 1, 2],
        orientation=orientation,
    )
    if orientation == "horizontal":
        colorbar.ax.set_xticklabels(
            ["强\n不利", "中等\n不利", "弱/中性", "中等\n有利", "强\n有利"]
        )
    else:
        colorbar.ax.set_yticklabels(
            ["强不利", "中等不利", "微弱/中性", "中等有利", "强有利"]
        )
    colorbar.ax.tick_params(labelsize=7.8, length=0)


def render_parent15_heatmap(data: V3ReportFigureData, stem: Path) -> tuple[Path, Path]:
    """Render parent-single U/S/Tm evidence with the T99F exploration marker."""

    rows = _selected_parent_rows(data)
    labels = [row["mutation_reported_label"].replace("Nb252 reported_seq ", "") for row in rows]
    fig = plt.figure(figsize=(7.35, 5.75))
    ax = fig.add_axes([0.15, 0.15, 0.52, 0.74])
    ax_notes = fig.add_axes([0.69, 0.15, 0.10, 0.74])
    ax_cbar = fig.add_axes([0.82, 0.25, 0.025, 0.54])
    image = _heatmap(
        ax,
        rows,
        ("netsolp_u_band_v3", "netsolp_s_band_v3", "nanomelt_tm_band_v3"),
        labels,
        ("NetSolP U", "NetSolP S", "NanoMelt预测Tm"),
    )
    fig.text(0.15, 0.94, "15条父单突的分档性质证据", ha="left", va="top", fontsize=12, weight="bold")
    ax.set_xlabel("弱变化统一归入中性显示，避免由微小数值差异驱动排序", fontsize=8.2, labelpad=9)
    ax_notes.set_ylim(ax.get_ylim())
    ax_notes.set_xlim(0, 1)
    ax_notes.set_xticks([])
    ax_notes.set_yticks([])
    ax_notes.set_title("注记", fontsize=8.3, pad=5)
    for index, row in enumerate(rows):
        notes: list[str] = []
        if labels[index] == "T99F":
            notes.append("稳定词探索")
        if notes:
            ax_notes.text(0.02, index, "\n".join(notes), ha="left", va="center", fontsize=7.0, color="#444444", linespacing=1.15)
    ax_notes.spines[:].set_visible(False)
    _add_band_colorbar(fig, image, cax=ax_cbar)
    fig.text(
        0.13,
        0.015,
        "T99F：稳定词探索项；AntiFold仅执行联合负向否决，不因预测改善而推荐任何单突。",
        ha="left",
        va="bottom",
        fontsize=7.8,
        linespacing=1.35,
    )
    return _save(fig, stem)


def render_double_selection_flow(data: V3ReportFigureData, stem: Path) -> tuple[Path, Path]:
    """Render the 102→42→15 double-selection flow and 6/9 evidence split."""

    facts = data.final_manifest["facts"]
    values = [facts["source_double_candidate_count"], facts["multi_metric_review_trigger_count"], facts["final_double_count"]]
    labels = ["全部有效双突", "最高性质档\n（2–3项中/强改善且无中/强不利）", "最终双突15条"]
    fig, ax = plt.subplots(figsize=(7.35, 2.55))
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    x = [0.12, 0.50, 0.85]
    for index, (xx, value, label) in enumerate(zip(x, values, labels)):
        ax.scatter(xx, 0.60, s=620, color="#194f7a" if index == 2 else "#5b8fb5", edgecolor="white", lw=1.2)
        ax.text(xx, 0.60, str(value), ha="center", va="center", color="white", weight="bold", fontsize=13)
        ax.text(xx, 0.29, label, ha="center", va="top", fontsize=8.6, linespacing=1.25)
        if index < 2:
            ax.annotate("", xy=(x[index + 1] - 0.075, 0.60), xytext=(xx + 0.075, 0.60),
                        arrowprops={"arrowstyle": "->", "color": "#59636b", "lw": 1.2})
    ax.text(0.85, 0.91, "6条：3项改善\n9条：2项改善", ha="center", va="top", fontsize=8.8, color="#194f7a")
    ax.text(
        0.5,
        0.02,
        "42条同属最高性质档，最终15条再结合性质幅度与组合多样性选择；AntiFold未对双突评分。",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=7.2,
    )
    ax.set_title("V3双突筛选路径", loc="left", fontsize=12, weight="bold", pad=8)
    fig.subplots_adjust(left=0.03, right=0.97, top=0.84, bottom=0.04)
    return _save(fig, stem)


def render_double15_heatmap(data: V3ReportFigureData, stem: Path) -> tuple[Path, Path]:
    """Render selected-double property bands and parent/position coverage."""

    rows = _selected_double_rows(data)
    labels = [row["mutation_set"].replace(";", "+") for row in rows]
    facts = data.final_manifest["facts"]
    fig = plt.figure(figsize=(7.35, 6.05))
    ax = fig.add_axes([0.14, 0.16, 0.46, 0.73])
    ax_cbar = fig.add_axes([0.65, 0.82, 0.28, 0.022])
    ax_summary = fig.add_axes([0.72, 0.27, 0.22, 0.47])
    image = _heatmap(
        ax,
        rows,
        ("netsolp_u_magnitude_band", "netsolp_s_magnitude_band", "nanomelt_tm_c_magnitude_band"),
        labels,
        ("NetSolP U", "NetSolP S", "NanoMelt预测Tm"),
    )
    fig.text(0.14, 0.94, "最终15条双突的分档性质证据", ha="left", va="top", fontsize=12, weight="bold")
    _add_band_colorbar(fig, image, cax=ax_cbar, orientation="horizontal")

    parent_coverage = facts["selected_parent_component_count"]
    position_coverage = facts["selected_reported_position_count"]
    summary_labels = ["父单突覆盖", "位置覆盖", "3项改善", "2项改善"]
    numerators = [parent_coverage, position_coverage, facts["selected_three_metric_positive_count"], facts["selected_two_metric_positive_count"]]
    denominators = [15, 12, 15, 15]
    y = np.arange(4)
    ax_summary.barh(y, denominators, color="#e6e8ea", height=0.45)
    ax_summary.barh(y, numerators, color="#3977a5", height=0.45)
    ax_summary.set_yticks(y, summary_labels, fontsize=8.2)
    ax_summary.set_xlim(0, 15)
    ax_summary.invert_yaxis()
    ax_summary.set_xticks([])
    ax_summary.set_title("组合覆盖与证据组成", fontsize=9.3, loc="left")
    for yy, numerator, denominator in zip(y, numerators, denominators):
        ax_summary.text(denominator + 0.25, yy, f"{numerator}/{denominator}", va="center", fontsize=8.2)
    ax_summary.spines[:].set_visible(False)
    fig.text(
        0.15,
        0.045,
        "15条均有2或3项中等/强改善，且无中等/强不利；展示顺序不是效力排序。",
        fontsize=7.8,
        ha="left",
        va="bottom",
        linespacing=1.35,
    )
    return _save(fig, stem)


def _validation_rows(data: V3ReportFigureData) -> list[Mapping[str, str]]:
    wanted = [
        ("NetSolP U*", "predicted_usability"),
        ("NetSolP S*", "predicted_solubility"),
        ("NanoMelt", "nanomelt_predicted_apparent_tm_c"),
        ("RP3Net", "rp3net_expression_probability"),
        ("PLM_Sol", "plm_sol_solubility_score"),
    ]
    indexed = {row["feature"]: row for row in data.validation_associations}
    return [dict(indexed[feature], display_label=label) for label, feature in wanted]


def render_tool_validation_summary(data: V3ReportFigureData, stem: Path) -> tuple[Path, Path]:
    """Render continuous, fixed-threshold, and tool-role evidence panels."""

    association_rows = _validation_rows(data)
    fixed5_order = ["NetSolP U", "NetSolP S", "RP3Net", "PLM_Sol"]
    fixed5_index = {row["predictor"]: row for row in data.fixed5_metrics}
    fixed5_rows = [fixed5_index[label] for label in fixed5_order]

    fig, axes = plt.subplots(2, 2, figsize=(7.35, 5.65))
    ax_assoc, ax_auc, ax_quality, ax_roles = axes.flat
    labels = [row["display_label"] for row in association_rows]
    rho = [float(row["stratified_spearman_rho"]) for row in association_rows]
    y = np.arange(len(labels))
    colors = ["#3977a5", "#3977a5", "#7aa6c4", "#9da4aa", "#9da4aa"]
    bars = ax_assoc.barh(y, rho, color=colors, height=0.58)
    ax_assoc.axvline(0, color="#555555", lw=0.8)
    ax_assoc.set_yticks(y, labels, fontsize=8)
    ax_assoc.invert_yaxis()
    ax_assoc.set_xlim(-0.1, 0.6)
    ax_assoc.set_xlabel("来源分层 Spearman ρ", fontsize=8.4)
    ax_assoc.set_title("A  连续关联", loc="left", fontsize=10, weight="bold")
    ax_assoc.bar_label(bars, labels=[f"{value:.3f}" for value in rho], padding=3, fontsize=7.5)
    ax_assoc.spines[["top", "right"]].set_visible(False)

    x = np.arange(len(fixed5_rows))
    width = 0.36
    roc = [float(row["roc_auc"]) for row in fixed5_rows]
    pr = [float(row["pr_auc_average_precision"]) for row in fixed5_rows]
    ax_auc.bar(x - width / 2, roc, width, label="ROC-AUC", color="#194f7a")
    ax_auc.bar(x + width / 2, pr, width, label="PR-AUC", color="#7db6d6")
    ax_auc.set_xticks(x, fixed5_order, rotation=20, ha="right", fontsize=7.5)
    ax_auc.set_ylim(0, 1)
    ax_auc.set_ylabel("0–1", fontsize=8)
    ax_auc.set_title("B  固定5 mg/L展示：区分度", loc="left", fontsize=10, weight="bold")
    ax_auc.legend(frameon=False, fontsize=7.2, ncols=2, loc="upper left")
    ax_auc.spines[["top", "right"]].set_visible(False)

    mcc = [float(row["mcc"]) for row in fixed5_rows]
    ba = [float(row["balanced_accuracy"]) for row in fixed5_rows]
    ax_quality.bar(x - width / 2, mcc, width, label="MCC", color="#3977a5")
    ax_quality.bar(x + width / 2, ba, width, label="Balanced accuracy", color="#9fc4dc")
    ax_quality.set_xticks(x, fixed5_order, rotation=20, ha="right", fontsize=7.5)
    ax_quality.set_ylim(0, 1)
    ax_quality.set_ylabel("0–1", fontsize=8)
    ax_quality.set_title("C  固定5 mg/L展示：阈值分类", loc="left", fontsize=10, weight="bold")
    ax_quality.legend(frameon=False, fontsize=7.2, ncols=2, loc="upper left")
    ax_quality.spines[["top", "right"]].set_visible(False)

    ax_roles.axis("off")
    ax_roles.set_title("D  V3中的工具角色", loc="left", fontsize=10, weight="bold")
    role_lines = [
        "AntiFold only excludes high-risk substitutions,",
        "never proposes, rewards, or ranks candidates.",
        "",
        "NetSolP U/S：同一模型的两个输出，分别保留分档证据。",
        "NanoMelt：预测稳定性约束，不作为已验证的产量排序器。",
        "RP3Net / PLM_Sol：当前证据不足，不进入候选筛选。",
        "5 mg/L分类仅用于展示，不作为工具准入或候选阈值。",
    ]
    ax_roles.text(0.01, 0.96, "\n".join(role_lines), ha="left", va="top", fontsize=7.8, linespacing=1.45)
    fig.suptitle("表达相关预测工具的项目内验证与角色", fontsize=12.5, weight="bold", y=0.985)
    fig.text(
        0.08,
        0.008,
        "* NetSolP U和S是同一模型的两个输出，不代表两个独立模型。固定5 mg/L结果使用序列簇留一，且仅纳入有个体数值产量的31条记录；NanoMelt未纳入该固定阈值展示。",
        fontsize=7.1,
        ha="left",
        va="bottom",
    )
    fig.subplots_adjust(left=0.12, right=0.97, top=0.91, bottom=0.12, hspace=0.48, wspace=0.47)
    return _save(fig, stem)


def build_compact_source_rows(data: V3ReportFigureData) -> list[dict[str, object]]:
    """Return a compact long table containing every plotted report value."""

    rows: list[dict[str, object]] = []

    def add(figure: str, item: str, metric: str, value: object, source: str, note: str = "") -> None:
        rows.append(
            {
                "figure_id": figure,
                "item_id": item,
                "metric": metric,
                "value": value,
                "source_artifact": source,
                "note": note,
            }
        )

    single_gate_source = (UPSTREAM_SINGLE_DIR / "expression_single_mutant_v3_gate.json").as_posix()
    for item, value in zip(
        ("complete_space", "after_antifold_veto", "qualified", "upstream_shortlist", "expert_review_pool", "selected_parents"),
        (847, 696, 61, 30, 31, 15),
    ):
        add("single_selection_flow", item, "candidate_count", value, single_gate_source)
    parent_source = (REPORT_SINGLE_DIR / "v3_parent_single_selection_audit.csv").as_posix()
    for row in _selected_parent_rows(data):
        label = row["mutation_reported_label"].replace("Nb252 reported_seq ", "")
        for metric, key in (
            ("NetSolP U band", "netsolp_u_band_v3"),
            ("NetSolP S band", "netsolp_s_band_v3"),
            ("NanoMelt predicted Tm band", "nanomelt_tm_band_v3"),
        ):
            add("parent15_property_heatmap", label, metric, row[key], parent_source)

    final_source = (REPORT_FINAL_DIR / "v3_double_mutant_final_selection_audit102.csv").as_posix()
    for item, value in (
        ("valid_doubles", 102),
        ("top_property_tier", 42),
        ("selected_doubles", 15),
        ("three_metric_selected", 6),
        ("two_metric_selected", 9),
    ):
        add("double_selection_flow", item, "candidate_count", value, final_source)
    for row in _selected_double_rows(data):
        label = row["mutation_set"].replace(";", "+")
        for metric, key in (
            ("NetSolP U band", "netsolp_u_magnitude_band"),
            ("NetSolP S band", "netsolp_s_magnitude_band"),
            ("NanoMelt predicted Tm band", "nanomelt_tm_c_magnitude_band"),
        ):
            add("double15_property_heatmap", label, metric, row[key], final_source)
    facts = data.final_manifest["facts"]
    selected_parents = _selected_parent_rows(data)
    parent_position_count = len(
        {int(row["reported_sequence_index_1based"]) for row in selected_parents}
    )
    for metric, value in (
        ("selected parent coverage", f"{facts['selected_parent_component_count']}/{len(selected_parents)}"),
        ("selected position coverage", f"{facts['selected_reported_position_count']}/{parent_position_count}"),
    ):
        add("double15_property_heatmap", "panel_summary", metric, value, final_source)

    for row in _validation_rows(data):
        add(
            "tool_validation_summary",
            row["display_label"],
            "source-stratified Spearman rho",
            row["stratified_spearman_rho"],
            row["source_artifact"],
        )
    for row in data.fixed5_metrics:
        for metric in ("roc_auc", "pr_auc_average_precision", "mcc", "balanced_accuracy"):
            add(
                "tool_validation_summary",
                row["predictor"],
                metric,
                row[metric],
                row["source_artifact"],
                "fixed 5 mg/L; leave-one-cluster-out; display only",
            )
    add(
        "tool_validation_summary",
        "AntiFold",
        "V3 role",
        "AntiFold only excludes high-risk substitutions, never proposes, rewards, or ranks candidates.",
        (UPSTREAM_SINGLE_DIR / "expression_single_mutant_v3_contract.json").as_posix(),
    )
    return rows


def write_compact_source_csv(data: V3ReportFigureData, path: Path) -> Path:
    """Write the compact source table for the report figures."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = build_compact_source_rows(data)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return path


def render_v3_report_figures(project_root: Path, output_dir: Path) -> Mapping[str, object]:
    """Render five V3 report figures and their compact source-data CSV."""

    _configure_style()
    data = load_v3_report_figure_data(project_root)
    output_dir = Path(output_dir)
    outputs: dict[str, object] = {
        "single_selection_flow": render_single_selection_flow(data, output_dir / "v3_single_selection_flow"),
        "parent15_property_heatmap": render_parent15_heatmap(data, output_dir / "v3_parent15_property_heatmap"),
        "double_selection_flow": render_double_selection_flow(data, output_dir / "v3_double_selection_flow"),
        "double15_property_heatmap": render_double15_heatmap(data, output_dir / "v3_double15_property_heatmap"),
        "tool_validation_summary": render_tool_validation_summary(data, output_dir / "v3_tool_validation_summary"),
    }
    outputs["source_data_csv"] = write_compact_source_csv(
        data, output_dir / "v3_report_figure_source_data.csv"
    )
    return outputs


__all__ = [
    "BAND_TO_GRADE",
    "V3ReportFigureData",
    "build_compact_source_rows",
    "load_v3_report_figure_data",
    "render_double15_heatmap",
    "render_double_selection_flow",
    "render_parent15_heatmap",
    "render_single_selection_flow",
    "render_tool_validation_summary",
    "render_v3_report_figures",
    "write_compact_source_csv",
]
