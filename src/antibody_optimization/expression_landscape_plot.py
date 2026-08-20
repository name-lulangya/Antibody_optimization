"""Four-metric landscape rendering for the released Nb252 single-mutant space."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Mapping, Sequence


AMINO_ACIDS = tuple("ACDEFGHIKLMNPQRSTVWY")
METRICS = (
    (
        "netsolp_delta_usability_vs_current_wt",
        "A  NetSolP usability",
        "ΔU vs WT",
    ),
    (
        "netsolp_delta_solubility_vs_current_wt",
        "B  NetSolP solubility",
        "ΔS vs WT",
    ),
    (
        "nanomelt_delta_predicted_apparent_tm_c_vs_current_wt",
        "C  NanoMelt apparent Tm",
        "ΔTm vs WT (°C)",
    ),
    (
        "antifold_landscape_delta_log_probability",
        "D  AntiFold: experimental complex with AF3 fallback",
        "Δlog probability vs WT",
    ),
)


class ExpressionLandscapeError(ValueError):
    """Raised when the released single-mutant matrix violates the plot contract."""


def build_expression_landscape_rows(
    matrix_rows: Sequence[Mapping[str, object]],
) -> tuple[list[dict[str, object]], dict[str, object]]:
    """Validate and normalize the 847-row matrix for plotting.

    The function preserves one row per released candidate. AntiFold uses the
    experimental-complex value when it is evaluable and the AF3 VHH-only value
    only for the 126 candidates whose experimental coordinates are missing.
    """

    if len(matrix_rows) != 847:
        raise ExpressionLandscapeError("Landscape requires exactly 847 candidate rows")
    seen_ids: set[str] = set()
    seen_sequences: set[str] = set()
    position_identity: dict[int, tuple[str, str]] = {}
    normalized: list[dict[str, object]] = []
    for source in matrix_rows:
        candidate_id = str(source["candidate_id"])
        sequence = str(source["sequence"])
        if candidate_id in seen_ids or sequence in seen_sequences:
            raise ExpressionLandscapeError("Candidate IDs and sequences must be unique")
        seen_ids.add(candidate_id)
        seen_sequences.add(sequence)
        position = int(source["reported_sequence_index_1based"])
        wt, mutant = str(source["wt_residue"]), str(source["mutant_residue"])
        region = str(source["region"])
        if mutant not in AMINO_ACIDS or wt not in AMINO_ACIDS or mutant == wt:
            raise ExpressionLandscapeError(f"Invalid substitution for {candidate_id}")
        prior = position_identity.setdefault(position, (wt, region))
        if prior != (wt, region):
            raise ExpressionLandscapeError(f"Inconsistent WT or region at reported position {position}")
        anti_status = str(source["experimental_complex_context_evaluation_status"])
        if anti_status not in {"pass", "not_evaluable"}:
            raise ExpressionLandscapeError(f"Unexpected AntiFold status for {candidate_id}: {anti_status}")
        anti_raw = str(source["experimental_complex_context_delta_log_probability"])
        if (anti_status == "pass") != bool(anti_raw):
            raise ExpressionLandscapeError(f"AntiFold value/status mismatch for {candidate_id}")
        af3_status = str(source["af3_vhh_only_evaluation_status"])
        af3_raw = str(source["af3_vhh_only_delta_log_probability"])
        if af3_status != "pass" or not af3_raw:
            raise ExpressionLandscapeError(f"AF3 AntiFold fallback unavailable for {candidate_id}")
        af3_value = _finite_float(af3_raw, candidate_id, "af3_vhh_only_delta_log_probability")
        experimental_value = (
            _finite_float(
                anti_raw, candidate_id, "experimental_complex_context_delta_log_probability"
            )
            if anti_status == "pass" else ""
        )
        landscape_value = experimental_value if anti_status == "pass" else af3_value
        landscape_source = (
            "experimental_complex_context"
            if anti_status == "pass"
            else "af3_vhh_only_fallback_for_missing_experimental_coordinates"
        )
        row: dict[str, object] = {
            "candidate_id": candidate_id,
            "reported_sequence_index_1based": position,
            "wt_residue": wt,
            "mutant_residue": mutant,
            "mutation_reported_label": str(source["mutation_reported_label"]),
            "imgt_position_label": str(source["imgt_position_label"]),
            "region": region,
            "stable_word_effect": str(source["stable_word_effect"]),
            "created_stable_word_occurrence_count": int(
                source["created_stable_word_occurrence_count"]
            ),
            "experimental_complex_context_evaluation_status": anti_status,
            "experimental_complex_context_delta_log_probability": experimental_value,
            "af3_vhh_only_evaluation_status": af3_status,
            "af3_vhh_only_delta_log_probability": af3_value,
            "antifold_landscape_source": landscape_source,
            "antifold_landscape_delta_log_probability": landscape_value,
        }
        for metric, _, _ in METRICS[:3]:
            row[metric] = _finite_float(source[metric], candidate_id, metric)
        normalized.append(row)
    normalized.sort(
        key=lambda row: (
            int(row["reported_sequence_index_1based"]),
            AMINO_ACIDS.index(str(row["mutant_residue"])),
        )
    )
    status_counts = Counter(
        str(row["experimental_complex_context_evaluation_status"]) for row in normalized
    )
    facts = {
        "candidate_count": len(normalized),
        "reported_position_count": len(position_identity),
        "stable_word_gain_candidate_count": sum(
            str(row["stable_word_effect"]) in {"gain_only", "net_gain"} for row in normalized
        ),
        "experimental_complex_antifold_pass_count": status_counts["pass"],
        "experimental_complex_antifold_not_evaluable_count": status_counts["not_evaluable"],
        "af3_antifold_fallback_count": sum(
            str(row["antifold_landscape_source"]).startswith("af3_vhh_only")
            for row in normalized
        ),
        "candidate_selection_performed": False,
    }
    if facts != {
        "candidate_count": 847,
        "reported_position_count": 48,
        "stable_word_gain_candidate_count": 22,
        "experimental_complex_antifold_pass_count": 721,
        "experimental_complex_antifold_not_evaluable_count": 126,
        "af3_antifold_fallback_count": 126,
        "candidate_selection_performed": False,
    }:
        raise ExpressionLandscapeError(f"Unexpected released landscape counts: {facts}")
    return normalized, facts


def render_expression_landscape(
    rows: Sequence[Mapping[str, object]],
    *,
    png_path: Path,
    svg_path: Path,
) -> None:
    """Render four aligned substitution heatmaps from normalized plot rows."""

    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    positions = sorted({int(row["reported_sequence_index_1based"]) for row in rows})
    position_index = {position: index for index, position in enumerate(positions)}
    aa_index = {amino_acid: index for index, amino_acid in enumerate(AMINO_ACIDS)}
    identity = {
        int(row["reported_sequence_index_1based"]): (str(row["wt_residue"]), str(row["region"]))
        for row in rows
    }
    matrices: dict[str, object] = {}
    for metric, _, _ in METRICS:
        matrix = np.full((len(AMINO_ACIDS), len(positions)), np.nan)
        for row in rows:
            value = row[metric]
            if value != "":
                matrix[
                    aa_index[str(row["mutant_residue"])],
                    position_index[int(row["reported_sequence_index_1based"])],
                ] = float(value)
        matrices[metric] = matrix
    stable_cells = [
        (
            position_index[int(row["reported_sequence_index_1based"])] + 0.5,
            aa_index[str(row["mutant_residue"])] + 0.5,
        )
        for row in rows
        if str(row["stable_word_effect"]) in {"gain_only", "net_gain"}
    ]
    mpl.rcParams["svg.hashsalt"] = "nb252-expression-landscape-v1"
    cmap = mpl.colormaps["RdBu"].with_extremes(bad="#E8E8E8")
    fig, axes = plt.subplots(2, 2, figsize=(19.2, 12.8), constrained_layout=False)
    region_colors = {"FR1": "#B8CCE4", "FR2": "#B8CCE4", "FR3": "#B8CCE4",
                     "CDR1": "#F2C879", "CDR2": "#F2C879", "CDR3": "#F2C879"}
    for panel_index, (axis, (metric, title, colorbar_label)) in enumerate(
        zip(axes.flat, METRICS, strict=True)
    ):
        matrix = matrices[metric]
        finite = matrix[np.isfinite(matrix)]
        limit = float(max(abs(finite.min()), abs(finite.max())))
        mesh = axis.pcolormesh(
            np.arange(len(positions) + 1),
            np.arange(len(AMINO_ACIDS) + 1),
            matrix,
            cmap=cmap,
            norm=mpl.colors.TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit),
            edgecolors="#FFFFFF",
            linewidth=0.22,
            rasterized=True,
        )
        axis.set_xlim(0, len(positions)); axis.set_ylim(0, len(AMINO_ACIDS))
        axis.set_xticks(np.arange(len(positions)) + 0.5)
        axis.set_xticklabels(
            [f"{identity[position][0]}{position}" for position in positions],
            rotation=90,
            fontsize=7,
        )
        axis.set_yticks(np.arange(len(AMINO_ACIDS)) + 0.5)
        axis.set_yticklabels(AMINO_ACIDS, fontsize=8)
        axis.set_xlabel("WT residue and reported-sequence position")
        axis.set_ylabel("Mutant residue")
        axis.set_title(title, loc="left", fontweight="normal")
        axis.tick_params(length=0)
        for spine in axis.spines.values():
            spine.set_linewidth(0.6); spine.set_color("#666666")
        _draw_region_strip(axis, positions, identity, region_colors)
        if stable_cells:
            axis.scatter(
                [cell[0] for cell in stable_cells],
                [cell[1] for cell in stable_cells],
                marker="*", s=34, c="#FFD43B", edgecolors="#111111", linewidths=0.45,
                zorder=4,
            )
        colorbar = fig.colorbar(mesh, ax=axis, fraction=0.027, pad=0.012)
        colorbar.set_label(colorbar_label)
        colorbar.ax.tick_params(labelsize=8)
    legend = [
        Line2D([0], [0], marker="*", linestyle="none", markersize=9,
               markerfacecolor="#FFD43B", markeredgecolor="#111111",
               label="Creates ≥1 stable-word occurrence"),
        Patch(facecolor="#D9D9D9", edgecolor="#666666",
              label="AntiFold panel: 721 experimental-complex + 126 AF3 fallback"),
        Patch(facecolor="#E8E8E8", edgecolor="none", label="Not in released candidate space"),
    ]
    fig.legend(handles=legend, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.018))
    fig.suptitle("Nb252 expression-design single-mutant landscape (847 released candidates)", y=0.992)
    fig.text(
        0.5, 0.002,
        "Blue is favorable relative to WT; red is unfavorable. AntiFold uses the experimental complex where evaluable and AF3 VHH-only values for 126 missing-coordinate candidates. Predictions are not measured expression or stability.",
        ha="center", fontsize=9,
    )
    fig.subplots_adjust(left=0.045, right=0.978, top=0.935, bottom=0.085, wspace=0.12, hspace=0.29)
    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=300, bbox_inches="tight", metadata={"Date": None})
    fig.savefig(svg_path, bbox_inches="tight", metadata={"Date": None})
    plt.close(fig)
    _normalize_svg(svg_path)


def render_expression_scatter(
    rows: Sequence[Mapping[str, object]],
    *,
    png_path: Path,
    svg_path: Path,
) -> None:
    """Render the current-space analogue of the historical two-panel scatter."""

    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    plt.rcParams["svg.hashsalt"] = "nb252-expression-scatter-v1"
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.4), constrained_layout=False)
    region_groups = {
        "Framework": ({"FR1", "FR2", "FR3", "FR4"}, "#4C78A8"),
        "CDR": ({"CDR1", "CDR2", "CDR3"}, "#E08E3C"),
    }
    for label, (regions, color) in region_groups.items():
        subset = [row for row in rows if str(row["region"]) in regions]
        axes[0].scatter(
            [float(row[METRICS[0][0]]) for row in subset],
            [float(row[METRICS[1][0]]) for row in subset],
            s=18, alpha=0.52, c=color, edgecolors="none", label=label,
            rasterized=True,
        )
        for source, marker in (
            ("experimental_complex_context", "o"),
            ("af3_vhh_only_fallback_for_missing_experimental_coordinates", "^"),
        ):
            selected = [row for row in subset if row["antifold_landscape_source"] == source]
            axes[1].scatter(
                [float(row[METRICS[3][0]]) for row in selected],
                [float(row[METRICS[2][0]]) for row in selected],
                s=20 if marker == "o" else 26, alpha=0.52, c=color, marker=marker,
                edgecolors="none", rasterized=True,
            )
    stable = [
        row for row in rows if str(row["stable_word_effect"]) in {"gain_only", "net_gain"}
    ]
    axes[0].scatter(
        [float(row[METRICS[0][0]]) for row in stable],
        [float(row[METRICS[1][0]]) for row in stable],
        s=72, marker="*", c="#FFD43B", edgecolors="#111111", linewidths=0.55, zorder=4,
    )
    axes[1].scatter(
        [float(row[METRICS[3][0]]) for row in stable],
        [float(row[METRICS[2][0]]) for row in stable],
        s=72, marker="*", c="#FFD43B", edgecolors="#111111", linewidths=0.55, zorder=4,
    )
    for axis in axes:
        axis.axhline(0, color="#666666", linewidth=0.8)
        axis.axvline(0, color="#666666", linewidth=0.8)
        axis.grid(color="#D8D8D8", linewidth=0.45, alpha=0.7)
        axis.spines[["top", "right"]].set_visible(False)
    axes[0].set(
        xlabel="NetSolP ΔU vs WT", ylabel="NetSolP ΔS vs WT",
        title="A  NetSolP relative predictions",
    )
    axes[1].set(
        xlabel="AntiFold Δlog probability vs WT",
        ylabel="NanoMelt predicted ΔTm vs WT (°C)",
        title="B  AntiFold compatibility and predicted thermal stability",
    )
    legend = [
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="#4C78A8",
               markeredgecolor="none", label="Framework position"),
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="#E08E3C",
               markeredgecolor="none", label="CDR position"),
        Line2D([0], [0], marker="o", linestyle="none", color="#666666",
               label="Experimental-complex AntiFold"),
        Line2D([0], [0], marker="^", linestyle="none", color="#666666",
               label="AF3 AntiFold fallback (126)"),
        Line2D([0], [0], marker="*", linestyle="none", markersize=10,
               markerfacecolor="#FFD43B", markeredgecolor="#111111",
               label="Creates ≥1 stable-word occurrence"),
    ]
    fig.legend(handles=legend, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.015))
    fig.suptitle("Nb252 expression-design single-mutant property scatter (847 candidates)", y=0.985)
    fig.text(
        0.5, 0.002,
        "All values are predictions relative to WT. AntiFold prioritizes the experimental complex and uses AF3 only where experimental coordinates are unavailable.",
        ha="center", fontsize=8.5,
    )
    fig.subplots_adjust(left=0.075, right=0.985, top=0.88, bottom=0.20, wspace=0.23)
    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=300, bbox_inches="tight", metadata={"Date": None})
    fig.savefig(svg_path, bbox_inches="tight", metadata={"Date": None})
    plt.close(fig)
    _normalize_svg(svg_path)


def _draw_region_strip(axis, positions, identity, colors) -> None:
    from matplotlib.patches import Rectangle

    groups: list[tuple[int, int, str]] = []
    start = 0
    current = identity[positions[0]][1]
    for index, position in enumerate(positions[1:], 1):
        region = identity[position][1]
        if region != current:
            groups.append((start, index, current))
            start, current = index, region
    groups.append((start, len(positions), current))
    for start, end, region in groups:
        axis.add_patch(
            Rectangle(
                (start, 1.012), end - start, 0.028,
                transform=axis.get_xaxis_transform(), clip_on=False,
                facecolor=colors[region], edgecolor="none",
            )
        )
        axis.text(
            (start + end) / 2, 1.047, region, transform=axis.get_xaxis_transform(),
            ha="center", va="bottom", fontsize=7, clip_on=False,
        )


def _normalize_svg(path: Path) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    path.write_text("\n".join(line.rstrip() for line in lines) + "\n", encoding="utf-8", newline="\n")


def _finite_float(value: object, candidate_id: str, field: str) -> float:
    import math

    try:
        parsed = float(value)
    except (TypeError, ValueError) as error:
        raise ExpressionLandscapeError(f"Non-numeric {field} for {candidate_id}") from error
    if not math.isfinite(parsed):
        raise ExpressionLandscapeError(f"Non-finite {field} for {candidate_id}")
    return parsed
