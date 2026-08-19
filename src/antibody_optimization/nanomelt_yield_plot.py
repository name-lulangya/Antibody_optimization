"""Plot NanoMelt predicted-apparent-Tm association with reported BL21 yield."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np


def render_nanomelt_yield_figure(
    sample_rows: Sequence[Mapping[str, object]],
    metric: Mapping[str, object],
    cv_rows: Sequence[Mapping[str, object]],
    leave_one_out_rows: Sequence[Mapping[str, object]],
    classification_rows: Sequence[Mapping[str, object]],
    classification_prediction_rows: Sequence[Mapping[str, object]],
    *,
    png_path: Path,
    svg_path: Path,
) -> None:
    """Render the exact numeric, stratified, ordinal, and transfer evidence."""

    numeric = [
        row for row in sample_rows
        if row["scoring_status"] == "pass" and row["observation_semantics"] == "individual_approximate"
    ]
    colors = {"LTT": "#0072B2", "WCC": "#D55E00"}
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 8.2), constrained_layout=True)

    for provider in ("LTT", "WCC"):
        selected = [row for row in numeric if row["provider_code"] == provider]
        axes[0, 0].scatter(
            [float(row["nanomelt_predicted_apparent_tm_c"]) for row in selected],
            [float(row["numeric_yield_value"]) for row in selected],
            label=provider,
            color=colors[provider],
            edgecolor="white",
            linewidth=0.7,
            s=48,
        )
    nb252 = next(row for row in numeric if row["sample_uid"] == "LTT__Nb252")
    axes[0, 0].scatter(
        [float(nb252["nanomelt_predicted_apparent_tm_c"])],
        [float(nb252["numeric_yield_value"])],
        marker="*",
        s=180,
        color="#000000",
        label="Nb252",
        zorder=4,
    )
    axes[0, 0].set(
        xlabel="NanoMelt predicted apparent Tm (°C; aligned VHH domain)",
        ylabel="Reported yield (source numeric value; log scale)",
        title="A  Individual numeric observations",
    )
    axes[0, 0].set_yscale("log")
    axes[0, 0].legend(frameon=False)

    effects = [
        float(metric["stratified_spearman_rho"]),
        float(metric["ltt_spearman_rho"]),
        float(metric["wcc_spearman_rho"]),
        float(metric["without_nb252_stratified_spearman_rho"]),
    ]
    labels = ["Stratified", "LTT", "WCC", "Without Nb252"]
    axes[0, 1].barh(labels, effects, color=["#009E73" if value >= 0 else "#CC79A7" for value in effects])
    axes[0, 1].axvline(0, color="#555555", linewidth=1)
    axes[0, 1].set(xlabel="Spearman ρ", xlim=(-1, 1), title="B  Direction and Nb252 sensitivity")

    schemes = ["leave_one_out", "leave_one_cluster_out"]
    metric_keys = ["roc_auc", "pr_auc_average_precision", "balanced_accuracy", "mcc"]
    metric_labels = ["ROC-AUC", "PR-AUC", "Balanced acc.", "MCC"]
    x = np.arange(len(metric_keys))
    width = 0.36
    for index, scheme in enumerate(schemes):
        row = next(item for item in classification_rows if item["outer_scheme"] == scheme)
        axes[1, 0].bar(
            x + (index - 0.5) * width,
            [float(row[key]) for key in metric_keys],
            width,
            label=scheme.replace("_", " "),
        )
    axes[1, 0].axhline(0, color="#777777", linewidth=0.8)
    axes[1, 0].set_xticks(x, metric_labels, rotation=20)
    axes[1, 0].set_ylim(-1, 1)
    axes[1, 0].set(title="C  Nested classification", ylabel="Held-out metric")
    axes[1, 0].legend(frameon=False, fontsize=8)

    loo_predictions = [
        row for row in classification_prediction_rows if row["outer_scheme"] == "leave_one_out"
    ]
    thresholds = [float(row["training_score_threshold"]) for row in loo_predictions]
    axes[1, 1].hist(
        thresholds,
        bins=min(10, max(4, len(set(thresholds)))),
        color="#999999",
        edgecolor="white",
    )
    axes[1, 1].axvline(float(np.median(thresholds)), color="#D55E00", linewidth=2, label="Median")
    axes[1, 1].set(
        xlabel="Training-fold NanoMelt Tm threshold (°C)",
        ylabel="Held-out samples",
        title="D  Classification-threshold stability",
    )
    axes[1, 1].legend(frameon=False)

    for axis in axes.flat:
        axis.title.set_fontweight("bold")
        axis.title.set_ha("left")
        axis.title.set_position((0, 1.0))
    for path in (png_path, svg_path):
        path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=600, bbox_inches="tight", facecolor="white")
    with plt.rc_context({"svg.hashsalt": "nb252-nanomelt-yield-v1"}):
        fig.savefig(svg_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
