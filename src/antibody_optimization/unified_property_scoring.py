"""Unified NetSolP/NanoMelt evidence for released Nb252 single mutants.

This module freezes the WT plus 1,962 currently released single mutants and
joins relative-to-WT property predictions to the existing AntiFold and
PyRosetta evidence.  Pareto layers are computed separately for the affinity
and stability/developability discovery tracks.  No weighted score, yield
prediction, experimental claim, or final candidate selection is produced.
"""

from __future__ import annotations

from collections import Counter
from typing import Mapping, Sequence

import numpy as np


WT_SCORE_ID = "LTT__Nb252__WT"
EXPECTED_CANDIDATES = 1962
EXPECTED_SCORE_ROWS = EXPECTED_CANDIDATES + 1


class UnifiedPropertyScoringError(ValueError):
    """Raised when unified property-scoring evidence is inconsistent."""


def build_property_samples(
    candidate_rows: Sequence[Mapping[str, object]],
    plan_gate: Mapping[str, object],
) -> list[dict[str, object]]:
    """Return WT plus the released candidates in deterministic input order."""

    if plan_gate.get("status") != "pass" or plan_gate.get("release") != "ready_for_unified_property_scoring":
        raise UnifiedPropertyScoringError("Unified AntiFold landscape is not released")
    eligible = [row for row in candidate_rows if str(row["design_status"]) == "eligible_current_round"]
    if len(eligible) != EXPECTED_CANDIDATES:
        raise UnifiedPropertyScoringError("Expected 1,962 released candidates")
    ids = [str(row["candidate_id"]) for row in eligible]
    if len(set(ids)) != EXPECTED_CANDIDATES:
        raise UnifiedPropertyScoringError("Candidate IDs are not unique")
    parent_sequences = set()
    for row in eligible:
        sequence = str(row["sequence"])
        index = int(row["sequence_index_1based"])
        wt = str(row["wt_residue"]); mutant = str(row["mutant_residue"])
        if len(sequence) != 128 or sequence[index - 1] != mutant:
            raise UnifiedPropertyScoringError(f"Invalid candidate sequence: {row['candidate_id']}")
        parent_sequences.add(sequence[: index - 1] + wt + sequence[index:])
    if len(parent_sequences) != 1:
        raise UnifiedPropertyScoringError("Candidates do not reconstruct one parent sequence")
    parent = parent_sequences.pop()
    rows = [{
        "score_id": WT_SCORE_ID,
        "candidate_id": "WT",
        "sequence_raw": parent,
        "design_track": "wild_type_control",
        "sequence_index_1based": "",
        "wt_residue": "",
        "mutant_residue": "",
        "is_wt_control": True,
    }]
    rows.extend({
        "score_id": str(row["candidate_id"]),
        "candidate_id": str(row["candidate_id"]),
        "sequence_raw": str(row["sequence"]),
        "design_track": str(row["design_track"]),
        "sequence_index_1based": int(row["sequence_index_1based"]),
        "wt_residue": str(row["wt_residue"]),
        "mutant_residue": str(row["mutant_residue"]),
        "is_wt_control": False,
    } for row in eligible)
    return rows


def build_property_evidence(
    candidate_rows: Sequence[Mapping[str, object]],
    sample_rows: Sequence[Mapping[str, object]],
    netsolp_rows: Sequence[Mapping[str, object]],
    nanomelt_rows: Sequence[Mapping[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, object]]:
    """Join predictions, calculate WT deltas, and assign track-specific Pareto layers."""

    if len(sample_rows) != EXPECTED_SCORE_ROWS:
        raise UnifiedPropertyScoringError("Property plan must contain WT plus 1,962 candidates")
    net = _unique(netsolp_rows, "sample_uid", EXPECTED_SCORE_ROWS, "NetSolP")
    melt = _unique(nanomelt_rows, "sample_uid", EXPECTED_SCORE_ROWS, "NanoMelt")
    samples = {str(row["score_id"]): row for row in sample_rows}
    sample_ids = set(samples)
    if len(samples) != EXPECTED_SCORE_ROWS:
        raise UnifiedPropertyScoringError("Property plan score IDs are not unique")
    if set(net) != sample_ids or set(melt) != sample_ids:
        raise UnifiedPropertyScoringError("Property score IDs do not match the plan")
    for identifier in sample_ids:
        if str(net[identifier]["scoring_status"]) != "pass" or str(melt[identifier]["scoring_status"]) != "pass":
            raise UnifiedPropertyScoringError(f"Property scoring did not pass for {identifier}")
        expected_sequence = str(samples[identifier]["sequence_raw"])
        if str(net[identifier]["sequence_raw"]) != expected_sequence or str(melt[identifier]["sequence_raw"]) != expected_sequence:
            raise UnifiedPropertyScoringError(f"Property sequence mismatch for {identifier}")
        if str(melt[identifier]["trimmed_n_terminal"]) or str(melt[identifier]["trimmed_c_terminal"]) != "GS":
            raise UnifiedPropertyScoringError(f"Unexpected NanoMelt scoring domain for {identifier}")
        if str(melt[identifier]["scored_ungapped_sequence"]) != expected_sequence[:-2]:
            raise UnifiedPropertyScoringError(f"NanoMelt did not preserve the expected 126-aa domain for {identifier}")

    wt_net = net[WT_SCORE_ID]; wt_melt = melt[WT_SCORE_ID]
    wt_u = float(wt_net["predicted_usability"]); wt_s = float(wt_net["predicted_solubility"])
    wt_tm = float(wt_melt["nanomelt_predicted_apparent_tm_c"])
    candidates = {str(row["candidate_id"]): row for row in candidate_rows if str(row["design_status"]) == "eligible_current_round"}
    if len(candidates) != EXPECTED_CANDIDATES:
        raise UnifiedPropertyScoringError("Candidate evidence does not contain 1,962 released rows")
    output: list[dict[str, object]] = []
    for sample in sample_rows[1:]:
        identifier = str(sample["score_id"]); source = candidates[identifier]
        u = float(net[identifier]["predicted_usability"]); s = float(net[identifier]["predicted_solubility"])
        tm = float(melt[identifier]["nanomelt_predicted_apparent_tm_c"])
        flags = [flag for flag in str(source.get("new_liability_flags", "")).split("|") if flag and flag.lower() != "nan"]
        output.append({
            **dict(source),
            "netsolp_predicted_usability": u,
            "netsolp_delta_usability_vs_wt": u - wt_u,
            "netsolp_predicted_solubility": s,
            "netsolp_delta_solubility_vs_wt": s - wt_s,
            "nanomelt_predicted_apparent_tm_c": tm,
            "nanomelt_delta_predicted_apparent_tm_c_vs_wt": tm - wt_tm,
            "nanomelt_scored_length_aa": int(melt[identifier]["scored_length_aa"]),
            "nanomelt_trimmed_c_terminal": str(melt[identifier]["trimmed_c_terminal"]),
            "chemical_risk_count": len(flags),
            "chemical_risk_review": "review" if flags else "none_detected_by_current_heuristics",
            "hard_constraint_status": "pass",
            "yield_prediction_performed": False,
            "candidate_selection_performed": False,
        })

    objectives = {
        "affinity_existing_interface_scan": [
            ("netsolp_delta_usability_vs_wt", 1),
            ("netsolp_delta_solubility_vs_wt", 1),
            ("nanomelt_delta_predicted_apparent_tm_c_vs_wt", 1),
            ("experimental_complex_context_delta_log_probability", 1),
            ("pyrosetta_delta_dG_separated_median", -1),
            ("pyrosetta_delta_cross_interface_energy_median", -1),
        ],
        "stability_developability_discovery": [
            ("netsolp_delta_usability_vs_wt", 1),
            ("netsolp_delta_solubility_vs_wt", 1),
            ("nanomelt_delta_predicted_apparent_tm_c_vs_wt", 1),
            ("experimental_complex_context_delta_log_probability", 1),
        ],
    }
    for track, definitions in objectives.items():
        track_rows = [row for row in output if row["design_track"] == track]
        layers = first_two_pareto_layers(track_rows, definitions)
        for row, layer in zip(track_rows, layers, strict=True):
            row["property_pareto_layer"] = layer
            row["preliminary_property_tier"] = "pareto_front_1" if layer == 1 else "pareto_front_2" if layer == 2 else "background"
            row["pareto_objectives"] = "|".join(f"{'max' if direction == 1 else 'min'}:{name}" for name, direction in definitions)

    output.sort(key=lambda row: (str(row["design_track"]), int(row["sequence_index_1based"]), str(row["mutant_residue"])))
    tier_counts = Counter((str(row["design_track"]), str(row["preliminary_property_tier"])) for row in output)
    summaries = [{"design_track": track, "preliminary_property_tier": tier, "candidate_count": count}
                 for (track, tier), count in sorted(tier_counts.items())]
    gate = {
        "schema_version": 1,
        "gate_name": "nb252_unified_single_mutant_property_scoring",
        "status": "pass",
        "candidate_count": len(output),
        "wt_control_count": 1,
        "netsolp_pass_count": len(net),
        "nanomelt_pass_count": len(melt),
        "track_counts": dict(sorted(Counter(str(row["design_track"]) for row in output).items())),
        "candidate_selection_performed": False,
        "weighted_composite_score_used": False,
        "yield_prediction_performed": False,
        "release": "ready_for_preliminary_property_pool_review",
        "interpretation": "Predicted relative property and compatibility evidence only; no measured affinity, stability, expression, Tm, or yield.",
    }
    return output, summaries, gate


def first_two_pareto_layers(
    rows: Sequence[Mapping[str, object]], objectives: Sequence[tuple[str, int]]
) -> list[int]:
    """Return Pareto layer 1, 2, or 3 (meaning third-or-later)."""

    if not rows or not objectives:
        raise UnifiedPropertyScoringError("Pareto inputs must be non-empty")
    values = np.asarray([[float(row[name]) * direction for name, direction in objectives] for row in rows], dtype=float)
    if not np.isfinite(values).all():
        raise UnifiedPropertyScoringError("Pareto objectives contain non-finite values")
    layers = np.full(len(rows), 3, dtype=int)
    remaining = np.arange(len(rows))
    for layer in (1, 2):
        current = values[remaining]
        front_mask = np.ones(len(remaining), dtype=bool)
        for index, point in enumerate(current):
            dominates = np.all(current >= point, axis=1) & np.any(current > point, axis=1)
            if np.any(dominates):
                front_mask[index] = False
        front = remaining[front_mask]
        layers[front] = layer
        remaining = remaining[~front_mask]
        if not len(remaining):
            break
    return layers.tolist()


def _unique(rows: Sequence[Mapping[str, object]], key: str, expected: int, label: str) -> dict[str, Mapping[str, object]]:
    indexed = {str(row[key]): row for row in rows}
    if len(rows) != expected or len(indexed) != expected:
        raise UnifiedPropertyScoringError(f"{label} must contain {expected} unique rows")
    return indexed
