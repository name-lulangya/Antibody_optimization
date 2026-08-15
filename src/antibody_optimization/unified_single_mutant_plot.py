"""Figures for the unified Nb252 single-mutant plan and AntiFold landscape."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Mapping, Sequence


def render_unified_space(
    positions: Sequence[Mapping[str, object]],
    candidates: Sequence[Mapping[str, object]],
    *,
    png_path: Path,
    svg_path: Path,
) -> None:
    """Render the exact position scope and candidate-status counts."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    plt.rcParams["svg.hashsalt"] = "nb252-unified-single-mutant"

    colors = {"eligible": "#3B82A0", "missing": "#E6A15A", "hard": "#A6A6A6"}
    fig, axes = plt.subplots(2, 1, figsize=(11, 6.2), gridspec_kw={"height_ratios": [1.1, 1]})
    ax = axes[0]
    for row in positions:
        x = int(row["sequence_index_1based"])
        category = "hard" if _bool(row["hard_immutable"]) else "missing" if _bool(row["experimental_missing_coordinates"]) else "eligible"
        ax.scatter(x, 0, marker="s", s=34, color=colors[category], edgecolor="none")
        if _bool(row["experimental_interface"]):
            ax.scatter(x, 0.18, marker="|", s=80, linewidth=1.5, color="#C44E52")
    ax.set_xlim(0, 129); ax.set_ylim(-0.25, 0.35); ax.set_yticks([])
    ax.set_xlabel("Nb252 reported-sequence position (1-based)")
    ax.set_title("A  Unified position scope (interface is cautious, not forbidden)", loc="left", fontsize=11)
    ax.legend(handles=[
        Line2D([0], [0], marker="s", color="none", markerfacecolor=colors["eligible"], label="Current-round position", markersize=8),
        Line2D([0], [0], marker="s", color="none", markerfacecolor=colors["missing"], label="Experimental coordinates missing; deferred", markersize=8),
        Line2D([0], [0], marker="s", color="none", markerfacecolor=colors["hard"], label="Hard immutable", markersize=8),
        Line2D([0], [0], marker="|", color="#C44E52", label="Experimental interface", markersize=12),
    ], ncol=2, frameon=False, loc="upper left")

    counts = Counter(str(row["design_status"]) for row in candidates)
    order = ["eligible_current_round", "deferred_missing_experimental_coordinates", "blocked_new_unpaired_cys"]
    labels = ["Eligible current round", "Deferred: missing coordinates", "Blocked: new unpaired Cys"]
    bars = axes[1].barh(labels, [counts[key] for key in order], color=[colors["eligible"], colors["missing"], "#C44E52"])
    axes[1].set_xlabel("Single substitutions")
    axes[1].set_title("B  Complete 19-substitution enumeration before scientific filtering", loc="left", fontsize=11)
    for bar, key in zip(bars, order, strict=True):
        axes[1].text(bar.get_width() + 18, bar.get_y() + bar.get_height() / 2, f"{counts[key]:,}", va="center")
    axes[1].set_xlim(0, max(counts.values()) * 1.13)
    fig.text(0.01, 0.01, "No affinity, AntiFold, stability, expression, or yield threshold applied.", fontsize=8, color="#444444")
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    fig.savefig(png_path, dpi=600, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight", metadata={"Date": None})
    plt.close(fig)
    _normalize_svg(svg_path)


def render_antifold_landscape(rows: Sequence[Mapping[str, object]], *, png_path: Path, svg_path: Path) -> None:
    """Render per-position complex-context compatibility and view sensitivity."""
    import matplotlib.pyplot as plt
    import numpy as np

    plt.rcParams["svg.hashsalt"] = "nb252-unified-antifold"

    eligible = [row for row in rows if row["design_status"] == "eligible_current_round"]
    x = np.asarray([int(row["sequence_index_1based"]) for row in eligible])
    y = np.asarray([float(row["experimental_complex_context_delta_log_probability"]) for row in eligible])
    af3 = np.asarray([float(row["af3_vhh_only_delta_log_probability"]) for row in eligible])
    core = np.asarray([_bool(row["affinity_core_module"]) for row in eligible])
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.7))
    axes[0].scatter(x[~core], y[~core], s=8, alpha=.35, color="#3B82A0", rasterized=True)
    axes[0].scatter(x[core], y[core], s=34, color="#C44E52", edgecolor="white", linewidth=.5, label="Existing affinity core")
    axes[0].axhline(0, color="#555555", lw=.8); axes[0].set_xlabel("Reported-sequence position")
    axes[0].set_ylabel("AntiFold Δlog probability\nexperimental complex context")
    axes[0].set_title("A  Primary complex-context landscape", loc="left", fontsize=11); axes[0].legend(frameon=False)
    axes[1].scatter(y[~core], af3[~core], s=8, alpha=.35, color="#3B82A0", rasterized=True)
    axes[1].scatter(y[core], af3[core], s=34, color="#C44E52", edgecolor="white", linewidth=.5)
    lo = min(y.min(), af3.min()); hi = max(y.max(), af3.max())
    axes[1].plot([lo, hi], [lo, hi], "--", color="#777777", lw=.8)
    axes[1].axhline(0, color="#999999", lw=.6); axes[1].axvline(0, color="#999999", lw=.6)
    axes[1].set_xlabel("Δlog probability: experimental complex")
    axes[1].set_ylabel("Δlog probability: AF3 VHH")
    axes[1].set_title("B  Experimental-versus-predicted sensitivity", loc="left", fontsize=11)
    fig.text(0.01, 0.01, "AntiFold is structure-conditioned sequence compatibility, not affinity, stability, expression, yield, or experiment.", fontsize=8)
    fig.tight_layout(rect=(0, .035, 1, 1)); fig.savefig(png_path, dpi=600, bbox_inches="tight"); fig.savefig(svg_path, bbox_inches="tight", metadata={"Date": None}); plt.close(fig)
    _normalize_svg(svg_path)


def _bool(value: object) -> bool:
    return value is True or str(value).lower() == "true"


def _normalize_svg(path: Path) -> None:
    """Remove generator whitespace and normalize SVG line endings."""
    lines = path.read_text(encoding="utf-8").splitlines()
    path.write_text("\n".join(line.rstrip() for line in lines) + "\n", encoding="utf-8", newline="\n")
