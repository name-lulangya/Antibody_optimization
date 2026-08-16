"""Deterministic Nb252 double-mutant design and scoring contracts.

The design space is the unordered pairwise product of the released 14-single
shortlist. Pairs that attempt two substitutions at the same reported sequence
position are invalid. Every other pair is retained without score-based
filtering. AntiFold evidence is carried forward only as the exact sum of
fixed-backbone per-position log-probability deltas; it is explicitly not a
double-mutant or epistasis prediction.
"""

from __future__ import annotations
from collections import Counter
from itertools import combinations
from statistics import median
from typing import Mapping, Sequence
from .unified_single_mutants import sequence_liability_deltas

PARENT_LENGTH = 128
EXPECTED_SINGLES = 14
EXPECTED_DOUBLES = 86
EXPECTED_INVALID = 5
REPLICATES = 3
WT_SCORE_ID = "LTT__Nb252__WT"

class DoubleMutantDesignError(ValueError):
    """Raised when released inputs or double-mutant results disagree."""

def build_double_mutant_space(shortlist_rows: Sequence[Mapping[str, object]], shortlist_gate: Mapping[str, object], mapping_rows: Sequence[Mapping[str, object]]) -> dict[str, object]:
    """Build the complete valid pairwise space and same-position audit."""
    if shortlist_gate.get("status") != "pass" or shortlist_gate.get("release") != "ready_for_small_combination_contract":
        raise DoubleMutantDesignError("Single-mutant shortlist is not released")
    if len(shortlist_rows) != EXPECTED_SINGLES or any(str(row.get("shortlist_decision")) != "retain_active" for row in shortlist_rows):
        raise DoubleMutantDesignError("Expected exactly 14 released active singles")
    parent_sequences = set(); normalized = []
    for source in shortlist_rows:
        row = dict(source); position = int(row["sequence_index_1based"]); wt = str(row["wt_residue"]); mutant = str(row["mutant_residue"]); sequence = str(row["sequence"])
        if len(sequence) != PARENT_LENGTH or sequence[position - 1] != mutant:
            raise DoubleMutantDesignError(f"Invalid single-mutant sequence: {row['mutation']}")
        parent_sequences.add(sequence[: position - 1] + wt + sequence[position:]); normalized.append((position, str(row["mutation"]), row))
    if len(parent_sequences) != 1: raise DoubleMutantDesignError("Active singles do not reconstruct one parent")
    parent = parent_sequences.pop()
    if parent[124:128] != "SSGS" or parent[21] != "C" or parent[94] != "C": raise DoubleMutantDesignError("Authoritative SSGS/disulfide constraints are not preserved")
    mapping = {int(row["sequence_index_1based"]): row for row in mapping_rows if str(row.get("source_model_role")) == "experimental_nk2r_nb252"}
    normalized.sort(key=lambda item: (item[0], item[1])); candidates = []; invalid = []
    for (_, mutation_a, a), (_, mutation_b, b) in combinations(normalized, 2):
        pos_a = int(a["sequence_index_1based"]); pos_b = int(b["sequence_index_1based"])
        if pos_a == pos_b:
            invalid.append({"mutation_a": mutation_a, "mutation_b": mutation_b, "sequence_index_1based": pos_a, "exclusion_reason": "same_position_mutually_exclusive_substitutions"}); continue
        mutant_sequence = list(parent); mutant_sequence[pos_a - 1] = str(a["mutant_residue"]); mutant_sequence[pos_b - 1] = str(b["mutant_residue"]); sequence = "".join(mutant_sequence)
        differences = [index + 1 for index, (x, y) in enumerate(zip(parent, sequence, strict=True)) if x != y]
        if differences != [pos_a, pos_b]: raise DoubleMutantDesignError("Double-mutant sequence does not contain exactly the declared substitutions")
        if any(position in {22, 95, 125, 126, 127, 128} for position in differences): raise DoubleMutantDesignError("A double mutant changes an immutable position")
        map_a = mapping.get(pos_a); map_b = mapping.get(pos_b)
        if not map_a or not map_b or map_a["coordinate_status"] != "observed" or map_b["coordinate_status"] != "observed": raise DoubleMutantDesignError("All positions must be observed in the experimental structure")
        track_a = str(a["design_track"]); track_b = str(b["design_track"]); pair_track = "_x_".join(sorted((track_a, track_b)))
        delta_a = float(a["antifold_complex_delta_log_probability"]); delta_b = float(b["antifold_complex_delta_log_probability"])
        liabilities = sequence_liability_deltas(parent, sequence)
        candidates.append({
            "candidate_id": f"Nb252_double_{mutation_a}__{mutation_b}", "mutation_set": f"{mutation_a};{mutation_b}", "combination_track": pair_track,
            "mutation_a": mutation_a, "mutation_b": mutation_b, "position_a": pos_a, "position_b": pos_b,
            "wt_a": str(a["wt_residue"]), "mutant_a": str(a["mutant_residue"]), "wt_b": str(b["wt_residue"]), "mutant_b": str(b["mutant_residue"]),
            "region_a": str(a["region"]), "region_b": str(b["region"]), "source_candidate_a": str(a["candidate_id"]), "source_candidate_b": str(b["candidate_id"]),
            "source_role_a": str(a["shortlist_role"]), "source_role_b": str(b["shortlist_role"]),
            "experimental_chain_a": str(map_a["auth_asym_id"]), "experimental_auth_seq_a": int(map_a["auth_seq_id"]), "experimental_insertion_a": str(map_a["insertion_code"]),
            "experimental_chain_b": str(map_b["auth_asym_id"]), "experimental_auth_seq_b": int(map_b["auth_seq_id"]), "experimental_insertion_b": str(map_b["insertion_code"]),
            "sequence": sequence, "antifold_additive_fixed_backbone_delta_log_probability": delta_a + delta_b,
            "antifold_additive_component_a": delta_a, "antifold_additive_component_b": delta_b,
            "antifold_double_model_rerun": False, "antifold_epistasis_evaluated": False,
            "hard_constraint_status": "pass", "candidate_filtering_applied": False, "scan_status": "planned_unfiltered",
            **liabilities,
        })
    counts = Counter(str(row["combination_track"]) for row in candidates); expected = {"affinity_x_affinity": 26, "affinity_x_property": 48, "property_x_property": 12}
    if len(candidates) != EXPECTED_DOUBLES or len(invalid) != EXPECTED_INVALID or dict(counts) != expected: raise DoubleMutantDesignError(f"Unexpected pair space: {len(candidates)}, {len(invalid)}, {dict(counts)}")
    return {"parent_sequence": parent, "candidates": candidates, "invalid_pairs": invalid, "facts": {"active_single_count": 14, "mathematical_pair_count": 91, "invalid_same_position_pair_count": 5, "valid_double_count": 86, "combination_track_counts": expected, "candidate_filtering_applied": False}}

def build_score_samples(parent_sequence: str, candidates: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    """Return WT plus all double-mutant sequences."""
    rows = [{"sample_uid": WT_SCORE_ID, "candidate_id": "WT", "sequence_raw": parent_sequence, "is_wt_control": True}]
    rows.extend({"sample_uid": str(row["candidate_id"]), "candidate_id": str(row["candidate_id"]), "sequence_raw": str(row["sequence"]), "is_wt_control": False} for row in candidates)
    return rows

def summarize_double_replicates(rows: Sequence[Mapping[str, object]], metric_fields: Sequence[str]) -> list[dict[str, object]]:
    """Summarize three paired repeats per double mutant."""
    grouped = {}
    for row in rows: grouped.setdefault(str(row["candidate_id"]), []).append(row)
    summaries = []
    for candidate_id, group in sorted(grouped.items()):
        if len(group) != REPLICATES or len({int(row["replicate"]) for row in group}) != REPLICATES: raise DoubleMutantDesignError(f"Expected three unique repeats for {candidate_id}")
        summary = {"candidate_id": candidate_id, "replicate_count": REPLICATES}
        for field in metric_fields:
            values = [float(row[field]) for row in group]; summary[f"{field}_median"] = median(values); summary[f"{field}_negative_count"] = sum(value < 0 for value in values)
        summaries.append(summary)
    if len(summaries) != EXPECTED_DOUBLES: raise DoubleMutantDesignError("Expected summaries for all doubles")
    return summaries
