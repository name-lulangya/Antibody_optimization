"""Compact PLM_Sol BL21-yield validation figure."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np


def render_plm_sol_yield_figure(
    sample_rows: Sequence[Mapping[str, object]],
    metric_rows: Sequence[Mapping[str, object]],
    classification_rows: Sequence[Mapping[str, object]],
    comparison_rows: Sequence[Mapping[str, object]],
    *,
    png_path: Path,
    svg_path: Path,
) -> None:
    """Render continuous, classification, and comparator evidence."""

    numeric = [row for row in sample_rows if row["observation_semantics"] == "individual_approximate"]
    primary = metric_rows[0]
    colors = {"LTT": "#0072B2", "WCC": "#D55E00"}
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 8.2), constrained_layout=True)
    for provider in ("LTT", "WCC"):
        selected = [row for row in numeric if row["provider_code"] == provider]
        axes[0, 0].scatter(
            [float(row["plm_sol_solubility_score"]) for row in selected],
            [float(row["numeric_yield_value"]) for row in selected],
            label=provider, color=colors[provider], edgecolor="white", linewidth=0.7, s=48,
        )
    axes[0, 0].set_yscale("log")
    axes[0, 0].set(xlabel="PLM_Sol solubility score", ylabel="Reported yield (mg/1 L)", title="A  Numeric observations")
    axes[0, 0].legend(frameon=False)

    labels = ["Pooled", "Within-provider", "Length-adjusted", "Cluster-CV increment"]
    values = [float(primary[key]) for key in (
        "pooled_spearman_rho", "stratified_spearman_rho",
        "length_adjusted_partial_spearman", "cluster_cv_increment_over_provider",
    )]
    axes[0, 1].barh(labels, values, color=["#009E73" if value >= 0 else "#CC79A7" for value in values])
    axes[0, 1].axvline(0, color="#555555", linewidth=1)
    axes[0, 1].set(xlim=(-1, 1), xlabel="Association / incremental Spearman rho", title="B  Continuous evidence")

    schemes = ["leave_one_out", "leave_one_cluster_out"]
    metrics = ["roc_auc", "pr_auc_average_precision", "balanced_accuracy", "mcc"]
    metric_labels = ["ROC-AUC", "PR-AUC", "Balanced acc.", "MCC"]
    x = np.arange(len(metrics))
    width = 0.36
    for index, scheme in enumerate(schemes):
        row = next(item for item in classification_rows if item["outer_scheme"] == scheme)
        axes[1, 0].bar(x + (index - 0.5) * width, [float(row[key]) for key in metrics], width, label=scheme.replace("_", " "))
    axes[1, 0].axhline(0, color="#777777", linewidth=0.8)
    axes[1, 0].set_xticks(x, metric_labels, rotation=20)
    axes[1, 0].set_ylim(-1, 1)
    axes[1, 0].set(title="C  Nested classification", ylabel="Held-out metric")
    axes[1, 0].legend(frameon=False, fontsize=8)

    correlations = comparison_rows[:3]
    names = [str(row["comparison"]) for row in correlations]
    all_values = [float(row["all_47_spearman"]) for row in correlations]
    numeric_values = [float(row["numeric_31_spearman"]) for row in correlations]
    x = np.arange(len(names))
    axes[1, 1].bar(x - 0.18, all_values, 0.36, label="All 47")
    axes[1, 1].bar(x + 0.18, numeric_values, 0.36, label="Numeric 31")
    increment = float(comparison_rows[3]["increment_over_netsolp_s"])
    axes[1, 1].set_xticks(x, names, rotation=18)
    axes[1, 1].set_ylim(-1, 1)
    axes[1, 1].set(title=f"D  Predictor overlap; NetSolP-S cluster-CV increment={increment:.3f}", ylabel="Spearman rho")
    axes[1, 1].legend(frameon=False)
    fig.suptitle("PLM_Sol validation against collaborator-reported BL21 yield", fontsize=14)
    fig.savefig(png_path, dpi=600)
    fig.savefig(svg_path)
    plt.close(fig)


def render_plm_sol_fixed5_figure(
    sample_rows: Sequence[Mapping[str, object]],
    metric_rows: Sequence[Mapping[str, object]],
    *,
    png_path: Path,
    svg_path: Path,
) -> None:
    """Render the explicitly display-only fixed-5-mg PLM_Sol view."""

    numeric = [row for row in sample_rows if row["observation_semantics"] == "individual_approximate"]
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.2), constrained_layout=True)
    low = [row for row in numeric if float(row["numeric_yield_value"]) < 5.0]
    high = [row for row in numeric if float(row["numeric_yield_value"]) >= 5.0]
    axes[0].scatter(
        [float(row["plm_sol_solubility_score"]) for row in low],
        [float(row["numeric_yield_value"]) for row in low], label="<5 mg", color="#D55E00",
    )
    axes[0].scatter(
        [float(row["plm_sol_solubility_score"]) for row in high],
        [float(row["numeric_yield_value"]) for row in high], label=">=5 mg", color="#0072B2",
    )
    axes[0].axhline(5.0, color="#777777", linestyle="--", linewidth=1)
    axes[0].set(xlabel="PLM_Sol score", ylabel="Reported yield (mg/1 L)", title="A  Fixed 5 mg labels")
    axes[0].legend(frameon=False)
    schemes = ["leave_one_out", "leave_one_cluster_out"]
    labels = ["ROC-AUC", "PR-AUC", "Balanced acc.", "MCC"]
    keys = ["roc_auc", "pr_auc_average_precision", "balanced_accuracy", "mcc"]
    x = np.arange(len(keys)); width = 0.36
    for index, scheme in enumerate(schemes):
        row = next(item for item in metric_rows if item["outer_scheme"] == scheme)
        axes[1].bar(x + (index - 0.5) * width, [float(row[key]) for key in keys], width, label=scheme.replace("_", " "))
    axes[1].set_xticks(x, labels, rotation=20)
    axes[1].set_ylim(-1, 1)
    axes[1].set(title="B  Held-out display metrics", ylabel="Metric")
    axes[1].legend(frameon=False, fontsize=8)
    fig.suptitle("PLM_Sol fixed-5-mg exploratory display (not a predictor gate)")
    fig.savefig(png_path, dpi=600)
    fig.savefig(svg_path)
    plt.close(fig)
