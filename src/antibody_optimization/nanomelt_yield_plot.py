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
    *,
    png_path: Path,
    svg_path: Path,
) -> None:
    """Render the exact numeric, stratified, ordinal, and transfer evidence."""

    numeric = [row for row in sample_rows if row["observation_semantics"] == "individual_approximate"]
    llj = [row for row in sample_rows if row["provider_code"] == "LLJ"]
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

    rng = np.random.default_rng(252)
    for level in (1, 2, 3):
        selected = [row for row in llj if int(row["llj_ordinal_level"]) == level]
        axes[1, 0].scatter(
            level + rng.uniform(-0.07, 0.07, len(selected)),
            [float(row["nanomelt_predicted_apparent_tm_c"]) for row in selected],
            color="#6A51A3",
            edgecolor="white",
            linewidth=0.6,
            s=42,
        )
    axes[1, 0].set_xticks([1, 2, 3], ["~2", "~10", ">20"])
    axes[1, 0].set(
        xlabel="LLJ reported yield group",
        ylabel="NanoMelt predicted apparent Tm (°C)",
        title="C  LLJ ordinal/censored context",
    )

    model = next(row for row in cv_rows if row["model"] == "provider_plus_nanomelt_tm")
    increments = [float(model["increment_over_provider"]), float(model["cluster_increment_over_provider"])]
    axes[1, 1].bar(["LOOCV", "Leave-cluster-out"], increments, color=["#56B4E9", "#E69F00"])
    axes[1, 1].axhline(0, color="#555555", linewidth=1)
    loo_values = [float(row["stratified_spearman_rho"]) for row in leave_one_out_rows]
    axes[1, 1].text(
        0.02,
        0.98,
        f"Leave-one-out ρ range: {min(loo_values):.2f} to {max(loo_values):.2f}",
        transform=axes[1, 1].transAxes,
        va="top",
        fontsize=9,
    )
    axes[1, 1].set(
        ylabel="Spearman increment over provider-only baseline",
        title="D  Out-of-sample incremental value",
    )

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
