"""Plot complete NetSolP/NanoMelt evidence for 102 V3 double mutants."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Mapping, Sequence


BAND_GRADE = {
    "strong_adverse": -2,
    "moderate_adverse": -1,
    "weak_adverse": 0,
    "negligible": 0,
    "weak_favorable": 0,
    "moderate_favorable": 1,
    "strong_favorable": 2,
}


def build_v3_double_result_plot_rows(
    rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    """Return three long-form metric rows per scored double."""

    specifications = (
        (
            "NetSolP U",
            "netsolp_u_delta_vs_wt",
            "netsolp_u_magnitude_band",
            "netsolp_u_model_nonadditivity_residual",
            0.015,
        ),
        (
            "NetSolP S",
            "netsolp_s_delta_vs_wt",
            "netsolp_s_magnitude_band",
            "netsolp_s_model_nonadditivity_residual",
            0.030,
        ),
        (
            "NanoMelt predicted Tm",
            "nanomelt_tm_c_delta_vs_wt",
            "nanomelt_tm_c_magnitude_band",
            "nanomelt_tm_c_model_nonadditivity_residual",
            1.5,
        ),
    )
    output: list[dict[str, object]] = []
    for source in rows:
        for metric, delta_field, band_field, residual_field, strong_threshold in specifications:
            band = str(source[band_field])
            output.append(
                {
                    "v3_double_plan_order_not_efficacy_rank": source[
                        "v3_double_plan_order_not_efficacy_rank"
                    ],
                    "double_candidate_id": source["double_candidate_id"],
                    "mutation_set": source["mutation_set"],
                    "metric": metric,
                    "delta_vs_wt": float(source[delta_field]),
                    "magnitude_band": band,
                    "magnitude_band_grade": BAND_GRADE[band],
                    "model_nonadditivity_residual": float(source[residual_field]),
                    "residual_scaled_by_strong_threshold": float(source[residual_field])
                    / strong_threshold,
                    "property_review_class": source["property_review_class"],
                    "machine_structure_triage_status": source[
                        "machine_structure_triage_status"
                    ],
                    "contains_t99f": source[
                        "contains_t99f_stable_word_exploration_parent"
                    ],
                }
            )
    if len(output) != 306:
        raise ValueError("V3 double-result plot requires 102 x 3 rows")
    return output


def render_v3_double_mutant_results(plot_rows, png_path, svg_path) -> None:
    """Render magnitude, property, non-additivity, and review-priority views."""

    import matplotlib.pyplot as plt
    import numpy as np

    metric_order = ["NetSolP U", "NetSolP S", "NanoMelt predicted Tm"]
    band_order = [
        "strong_adverse",
        "moderate_adverse",
        "weak_adverse",
        "negligible",
        "weak_favorable",
        "moderate_favorable",
        "strong_favorable",
    ]
    colors = {
        "strong_adverse": "#7f1d1d",
        "moderate_adverse": "#df7a68",
        "weak_adverse": "#f2d0c9",
        "negligible": "#eeeeee",
        "weak_favorable": "#d3e7f1",
        "moderate_favorable": "#6baed6",
        "strong_favorable": "#08306b",
    }
    band_labels = {value: value.replace("_", " ") for value in band_order}
    figure, axes = plt.subplots(2, 2, figsize=(11.4, 8.0))

    axis = axes[0, 0]
    bottoms = np.zeros(3)
    for band in band_order:
        values = [
            sum(
                row["metric"] == metric and row["magnitude_band"] == band
                for row in plot_rows
            )
            for metric in metric_order
        ]
        axis.bar(
            metric_order,
            values,
            bottom=bottoms,
            color=colors[band],
            label=band_labels[band],
        )
        bottoms += np.asarray(values)
    axis.set_ylabel("Double mutants")
    axis.set_title("A  Frozen magnitude bands", loc="left")
    axis.tick_params(axis="x", rotation=18)
    axis.legend(frameon=False, fontsize=7, ncol=2, loc="upper right")

    by_candidate: dict[str, dict[str, Mapping[str, object]]] = {}
    for row in plot_rows:
        by_candidate.setdefault(str(row["double_candidate_id"]), {})[
            str(row["metric"])
        ] = row
    axis = axes[0, 1]
    for identifier, metrics in by_candidate.items():
        u = metrics["NetSolP U"]
        s = metrics["NetSolP S"]
        t = metrics["NanoMelt predicted Tm"]
        marker = "*" if str(u["contains_t99f"]).lower() == "true" else "o"
        axis.scatter(
            float(u["delta_vs_wt"]),
            float(s["delta_vs_wt"]),
            c=float(t["magnitude_band_grade"]),
            cmap="RdBu",
            vmin=-2,
            vmax=2,
            marker=marker,
            s=48 if marker == "*" else 24,
            edgecolor="#333333",
            linewidth=0.25,
            alpha=0.85,
        )
    axis.axhline(0, color="#777777", linewidth=0.7)
    axis.axvline(0, color="#777777", linewidth=0.7)
    axis.set_xlabel("NetSolP ΔU")
    axis.set_ylabel("NetSolP ΔS")
    axis.set_title("B  Complete-sequence property landscape", loc="left")
    axis.text(0.02, 0.02, "★ contains T99F", transform=axis.transAxes, fontsize=8)

    axis = axes[1, 0]
    residuals = [
        [
            float(row["residual_scaled_by_strong_threshold"])
            for row in plot_rows
            if row["metric"] == metric
        ]
        for metric in metric_order
    ]
    boxes = axis.boxplot(residuals, tick_labels=metric_order, patch_artist=True)
    for patch, color in zip(boxes["boxes"], ("#9ecae1", "#6baed6", "#fdae6b"), strict=True):
        patch.set_facecolor(color)
    axis.axhline(0, color="#555555", linewidth=0.8)
    axis.set_ylabel("Residual / strong-band threshold")
    axis.set_title("C  Predictor non-additivity diagnostic", loc="left")
    axis.tick_params(axis="x", rotation=18)
    axis.text(
        0.02,
        0.02,
        "Model-output residual; not physical epistasis",
        transform=axis.transAxes,
        fontsize=8,
    )

    axis = axes[1, 1]
    candidate_rows = [metrics["NetSolP U"] for metrics in by_candidate.values()]
    table = Counter(
        (
            str(row["property_review_class"]),
            str(row["machine_structure_triage_status"]),
        )
        for row in candidate_rows
    )
    property_classes = [
        value
        for value in (
            "strong_adverse_property_requires_review",
            "property_supported_with_moderate_adverse_tradeoff",
            "property_supported_no_moderate_or_strong_adverse",
            "no_moderate_or_strong_positive_metric",
        )
        if any(key[0] == value for key in table)
    ]
    display_labels = {
        "strong_adverse_property_requires_review": "Strong adverse\nproperty",
        "property_supported_with_moderate_adverse_tradeoff": (
            "Supported with\nmoderate trade-off"
        ),
        "property_supported_no_moderate_or_strong_adverse": (
            "Supported; no\nmoderate/strong adverse"
        ),
        "no_moderate_or_strong_positive_metric": "No moderate/strong\npositive",
    }
    routine = [table[(value, "routine_context_recorded")] for value in property_classes]
    detailed = [table[(value, "detailed_review_triggered")] for value in property_classes]
    y = np.arange(len(property_classes))
    axis.barh(y, routine, color="#8bb8d3", label="routine structure context")
    axis.barh(y, detailed, left=routine, color="#efaa62", label="detailed-review trigger")
    axis.set_yticks(
        y,
        [display_labels[value] for value in property_classes],
        fontsize=7.5,
    )
    axis.set_xlabel("Double mutants")
    axis.set_title("D  Post-score review strata", loc="left")
    axis.legend(frameon=False, fontsize=8)

    figure.suptitle(
        "Nb252 V3 complete 102-double property matrix (no final selection)",
        fontsize=12,
    )
    figure.tight_layout()
    figure.savefig(png_path, dpi=600, bbox_inches="tight")
    figure.savefig(svg_path, bbox_inches="tight")
    plt.close(figure)
    svg_output = Path(svg_path)
    svg_text = svg_output.read_text(encoding="utf-8")
    svg_output.write_text(
        "\n".join(line.rstrip() for line in svg_text.splitlines()) + "\n",
        encoding="utf-8",
        newline="\n",
    )


__all__ = ["build_v3_double_result_plot_rows", "render_v3_double_mutant_results"]
