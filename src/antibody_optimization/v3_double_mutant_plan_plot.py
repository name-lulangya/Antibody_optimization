"""Render the unfiltered V3 15-parent to 102-double enumeration overview."""

from __future__ import annotations

from collections import Counter
from typing import Mapping, Sequence

import numpy as np


def build_v3_double_plan_plot_rows(
    parents: Sequence[Mapping[str, object]],
    candidates: Sequence[Mapping[str, object]],
    invalid_pairs: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    """Return one exact plotting row for each of the 105 theoretical pairs."""

    order = {
        str(row["candidate_id"]): int(
            row["v3_parent_panel_order_not_efficacy_rank"]
        )
        for row in parents
    }
    mutation = {
        str(row["candidate_id"]): (
            f"{row['wt_residue']}{int(row['reported_sequence_index_1based'])}"
            f"{row['mutant_residue']}"
        )
        for row in parents
    }
    rows: list[dict[str, object]] = []
    for source in candidates:
        first = str(source["parent_a_candidate_id"])
        second = str(source["parent_b_candidate_id"])
        rows.append(
            {
                "theoretical_pair_order": int(source["theoretical_pair_order"]),
                "parent_a_order": order[first],
                "parent_b_order": order[second],
                "parent_a_mutation": mutation[first],
                "parent_b_mutation": mutation[second],
                "pair_status": "valid_released_for_complete_scoring",
                "structure_triage_status": source[
                    "machine_structure_triage_status"
                ],
                "pair_structure_distance_source": source[
                    "pair_structure_distance_source"
                ],
                "pair_spatial_class": source["pair_spatial_class"],
                "stable_word_effect": source["stable_word_effect"],
                "hard_sequence_risk_count": int(
                    source["hard_sequence_risk_count"]
                ),
                "contains_t99f": bool(
                    source["contains_t99f_stable_word_exploration_parent"]
                ),
            }
        )
    for source in invalid_pairs:
        first = str(source["parent_a_candidate_id"])
        second = str(source["parent_b_candidate_id"])
        rows.append(
            {
                "theoretical_pair_order": int(source["theoretical_pair_order"]),
                "parent_a_order": order[first],
                "parent_b_order": order[second],
                "parent_a_mutation": mutation[first],
                "parent_b_mutation": mutation[second],
                "pair_status": "invalid_same_position",
                "structure_triage_status": "not_applicable",
                "pair_structure_distance_source": "not_applicable",
                "pair_spatial_class": "not_applicable",
                "stable_word_effect": "not_applicable",
                "hard_sequence_risk_count": 0,
                "contains_t99f": False,
            }
        )
    rows.sort(key=lambda row: int(row["theoretical_pair_order"]))
    if len(rows) != 105:
        raise ValueError("V3 double-plan plot requires 105 theoretical pairs")
    return rows


def render_v3_double_mutant_plan(
    plot_rows: Sequence[Mapping[str, object]],
    png_path,
    svg_path,
) -> None:
    """Render legal-pair topology, counts, and pre-score triage summaries."""

    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch

    labels = [""] * 15
    matrix = np.full((15, 15), np.nan)
    np.fill_diagonal(matrix, 0)
    for row in plot_rows:
        first = int(row["parent_a_order"]) - 1
        second = int(row["parent_b_order"]) - 1
        labels[first] = str(row["parent_a_mutation"])
        labels[second] = str(row["parent_b_mutation"])
        value = 1 if row["pair_status"] == "valid_released_for_complete_scoring" else -1
        matrix[first, second] = value
        matrix[second, first] = value
    masked = np.ma.masked_invalid(matrix)
    figure = plt.figure(figsize=(12.6, 4.9))
    grid = figure.add_gridspec(1, 3, width_ratios=(1.7, 0.72, 1.08), wspace=0.38)

    axis_a = figure.add_subplot(grid[0, 0])
    axis_a.imshow(
        masked,
        cmap=ListedColormap(["#c95b52", "#eeeeee", "#4f92bd"]),
        vmin=-1,
        vmax=1,
        interpolation="none",
    )
    axis_a.set_xticks(range(15), labels, rotation=65, ha="right", fontsize=7.5)
    axis_a.set_yticks(range(15), labels, fontsize=7.5)
    axis_a.set_title("A  Complete unordered-pair map", loc="left", fontsize=10.5)
    axis_a.set_xlabel("Parent single (display order; not efficacy rank)")
    axis_a.set_ylabel("Parent single")
    axis_a.legend(
        handles=(
            Patch(facecolor="#4f92bd", label="valid double"),
            Patch(facecolor="#c95b52", label="same-position invalid"),
            Patch(facecolor="#eeeeee", label="self pair"),
        ),
        frameon=False,
        fontsize=7.5,
        loc="lower right",
    )

    axis_b = figure.add_subplot(grid[0, 1])
    valid = sum(
        row["pair_status"] == "valid_released_for_complete_scoring"
        for row in plot_rows
    )
    invalid = len(plot_rows) - valid
    bars = axis_b.bar(
        ["theoretical", "valid", "invalid"],
        [len(plot_rows), valid, invalid],
        color=["#7b8a97", "#4f92bd", "#c95b52"],
        width=0.68,
    )
    axis_b.bar_label(bars, padding=3, fontsize=9)
    axis_b.set_ylim(0, 116)
    axis_b.set_ylabel("Pair count")
    axis_b.set_title("B  Enumeration", loc="left", fontsize=10.5)
    axis_b.tick_params(axis="x", rotation=45, labelsize=8)
    axis_b.spines[["top", "right"]].set_visible(False)

    axis_c = figure.add_subplot(grid[0, 2])
    valid_rows = [
        row
        for row in plot_rows
        if row["pair_status"] == "valid_released_for_complete_scoring"
    ]
    counts = Counter(str(row["structure_triage_status"]) for row in valid_rows)
    categories = ["routine_context_recorded", "detailed_review_triggered"]
    values = [counts.get(category, 0) for category in categories]
    bars = axis_c.barh(
        ["routine context", "detailed-review trigger"],
        values,
        color=["#8bb8d3", "#efaa62"],
    )
    axis_c.bar_label(bars, padding=3, fontsize=9)
    axis_c.set_xlim(0, max(values + [1]) * 1.18)
    axis_c.set_xlabel("Valid doubles")
    axis_c.set_title("C  WT-structure triage", loc="left", fontsize=10.5)
    axis_c.spines[["top", "right"]].set_visible(False)
    t99f = sum(bool(row["contains_t99f"]) for row in valid_rows)
    stable_gain = sum(
        str(row["stable_word_effect"]) in {"gain_only", "net_gain"}
        for row in valid_rows
    )
    hard_risk = sum(int(row["hard_sequence_risk_count"]) > 0 for row in valid_rows)
    axis_c.text(
        0.98,
        0.04,
        (
            f"T99F-containing: {t99f}\n"
            f"stable-word net gain: {stable_gain}\n"
            f"hard sequence-risk flags: {hard_risk}"
        ),
        transform=axis_c.transAxes,
        ha="right",
        va="bottom",
        fontsize=8.5,
        linespacing=1.4,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.88},
    )

    figure.suptitle(
        "Nb252 V3 double-mutant plan: 15 released parents, no property prefilter",
        fontsize=12,
        y=1.02,
    )
    figure.savefig(png_path, dpi=600, bbox_inches="tight")
    figure.savefig(svg_path, bbox_inches="tight")
    plt.close(figure)
    _normalize_svg_line_endings(svg_path)


def _normalize_svg_line_endings(svg_path) -> None:
    """Remove Matplotlib path-line padding so Git whitespace checks stay clean."""

    text = svg_path.read_text(encoding="utf-8")
    svg_path.write_text(
        "\n".join(line.rstrip() for line in text.splitlines()) + "\n",
        encoding="utf-8",
        newline="\n",
    )


__all__ = ["build_v3_double_plan_plot_rows", "render_v3_double_mutant_plan"]
