"""Plot exact targeted-structure qualification results."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Mapping, Sequence

import matplotlib.pyplot as plt


def render_targeted_review(rows: Sequence[Mapping[str, object]], png: Path, svg: Path) -> None:
    """Render qualification counts and paired structural diagnostics."""

    reviewed = [row for row in rows if str(row.get("targeted_runtime_reviewed")).lower() == "true" or row.get("targeted_runtime_reviewed") is True]
    counts = Counter(str(row["v2_qualification_status"]) for row in rows)
    order = ["combination_ready", "single_mutant_test_only", "blocked_pending_structure", "targeted_alternative_review", "not_prioritized", "blocked", "do_not_advance"]
    colors = {"combination_ready": "#009E73", "single_mutant_test_only": "#E69F00", "blocked_pending_structure": "#CC79A7", "targeted_alternative_review": "#56B4E9", "not_prioritized": "#999999", "blocked": "#D55E00", "do_not_advance": "#7A0019"}
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5), constrained_layout=True)
    present = [status for status in order if counts.get(status)]
    bars = axes[0].barh(present, [counts[status] for status in present], color=[colors[status] for status in present])
    axes[0].bar_label(bars, padding=3)
    axes[0].invert_yaxis(); axes[0].set_xlabel("Candidates"); axes[0].set_title("A  V2 qualification")
    axes[0].spines[["top", "right"]].set_visible(False)
    selected = [row for row in reviewed if row.get("review_group") == "gap_boundary_nonproline"]
    if selected:
        axes[1].scatter(
            [float(row["median_af3_vhh_delta_total_score"]) for row in selected],
            [float(row["median_af3_vhh_delta_local_fa_rep"]) for row in selected],
            label="AF3 gap-boundary non-Pro", color="#0072B2", marker="o", s=48, alpha=0.85,
        )
    axes[1].axvline(0, color="#444444", lw=1, ls="--")
    axes[1].axhline(0, color="#444444", lw=1, ls=":")
    axes[1].set_xlabel("Median Δ total score, AF3 VHH (REU)")
    axes[1].set_ylabel("Median Δ local fa_rep, AF3 VHH (REU)")
    axes[1].set_title("B  Nonredundant AF3 local review")
    axes[1].legend(frameon=False, fontsize=8)
    axes[1].spines[["top", "right"]].set_visible(False)
    fig.text(0.01, 0.01, "Computational triage only; prior complex evidence is reused and exact contact-set equality is not a hard gate.", fontsize=8)
    fig.savefig(png, dpi=600)
    fig.savefig(svg)
    plt.close(fig)
    lines = svg.read_text(encoding="utf-8").splitlines()
    svg.write_text("\n".join(line.rstrip() for line in lines) + "\n", encoding="utf-8", newline="\n")
