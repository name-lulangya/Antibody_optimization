"""Rendering for the compact Nb252 AntiFold validation result."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np


VIEW_FIELDS = (
    ("experimental_vhh_only_delta_log_probability", "Experimental VHH only"),
    ("experimental_complex_context_delta_log_probability", "Experimental complex"),
    ("af3_vhh_only_delta_log_probability", "AF3 VHH"),
)


def render_antifold_validation(
    rows: Sequence[Mapping[str, object]], *, png_path: Path, svg_path: Path
) -> None:
    """Render candidate deltas and experimental-versus-AF3 sensitivity."""

    labels = [f"{row['wt_residue']}{row['sequence_index_1based']}{row['mutant_residue']}" for row in rows]
    x = np.arange(len(rows), dtype=float)
    colors = ("#4C78A8", "#F58518", "#54A24B")
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 5.0), constrained_layout=True)
    width = 0.24
    for offset, ((field, label), color) in enumerate(zip(VIEW_FIELDS, colors, strict=True)):
        values = [float(row[field]) for row in rows]
        axes[0].bar(x + (offset - 1) * width, values, width=width, label=label, color=color)
    axes[0].axhline(0, color="#444444", linewidth=0.8)
    axes[0].set_xticks(x, labels, rotation=45, ha="right")
    axes[0].set_ylabel("Δ log P(mutant − WT), AntiFold")
    axes[0].set_title("A  Structure-conditioned substitution compatibility")
    axes[0].legend(frameon=False, fontsize=8)

    exp = np.asarray([float(row["experimental_complex_context_delta_log_probability"]) for row in rows])
    af3 = np.asarray([float(row["af3_vhh_only_delta_log_probability"]) for row in rows])
    perplexity = np.asarray([float(row["experimental_complex_context_perplexity"]) for row in rows])
    sizes = 35 + 8 * np.sqrt(perplexity)
    edgecolors = ["#B22222" if str(row.get("risk_flags", "")) else "#333333" for row in rows]
    axes[1].scatter(exp, af3, s=sizes, c="#72B7B2", edgecolors=edgecolors, linewidths=1.2)
    for index, label in enumerate(labels):
        axes[1].annotate(label, (exp[index], af3[index]), xytext=(4, 4), textcoords="offset points", fontsize=8)
    lo = float(min(exp.min(), af3.min(), 0.0))
    hi = float(max(exp.max(), af3.max(), 0.0))
    pad = max(0.2, (hi - lo) * 0.08)
    axes[1].plot([lo - pad, hi + pad], [lo - pad, hi + pad], linestyle="--", color="#777777", linewidth=0.8)
    axes[1].axhline(0, color="#BBBBBB", linewidth=0.7)
    axes[1].axvline(0, color="#BBBBBB", linewidth=0.7)
    axes[1].set_xlim(lo - pad, hi + pad)
    axes[1].set_ylim(lo - pad, hi + pad)
    axes[1].set_xlabel("Experimental complex Δ log P")
    axes[1].set_ylabel("AF3 VHH Δ log P")
    axes[1].set_title("B  Conformation/context sensitivity")
    axes[1].text(
        0.02, 0.02, "Marker size: experimental-complex perplexity\nRed edge: retained upstream risk flag",
        transform=axes[1].transAxes, fontsize=8, va="bottom",
    )
    fig.suptitle(
        "Nb252 AntiFold minimal validation\nCompatibility signal only; not affinity, expression, yield, Tm, or experimental evidence",
        fontsize=12,
    )
    fig.savefig(png_path, dpi=600)
    fig.savefig(svg_path)
    plt.close(fig)
