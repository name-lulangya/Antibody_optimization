"""Degenerate stable-word features for Nb252 single-mutant evaluation."""

from __future__ import annotations

from collections import defaultdict
from typing import Mapping, Sequence


AA_TO_DEGENERATE = {
    "D": "a", "E": "a",
    "K": "b", "R": "b",
    "L": "h", "I": "h", "V": "h",
    "N": "n", "Q": "n",
    "S": "o", "T": "o",
    "A": "s", "G": "s",
    "F": "@", "W": "@",
    "C": "C", "H": "H", "M": "M", "P": "P", "Y": "Y",
}
STABLE_WORD_ALPHABET = frozenset("abhnos@CHMPY")
PRIMARY_YIELD_FEATURE = "stable_word_window_normalized_density"
YIELD_FEATURES = (
    "stable_word_occurrence_count",
    "stable_word_occurrences_per_residue",
    PRIMARY_YIELD_FEATURE,
    "stable_word_unique_fraction",
)


class StableWordError(ValueError):
    """Raised when a stable-word input or sequence violates the fixed contract."""


def parse_stable_words(lines: Sequence[str]) -> tuple[str, ...]:
    """Validate a case-sensitive, one-word-per-line stable-word library.

    Blank lines, surrounding whitespace, unsupported symbols, and exact duplicate
    words are rejected rather than silently normalized.
    """

    words: list[str] = []
    for line_number, raw in enumerate(lines, 1):
        word = raw.rstrip("\r\n")
        if not word:
            raise StableWordError(f"Stable-word line {line_number} is blank")
        if word != word.strip():
            raise StableWordError(f"Stable-word line {line_number} has surrounding whitespace")
        unsupported = sorted(set(word) - STABLE_WORD_ALPHABET)
        if unsupported:
            raise StableWordError(
                f"Stable-word line {line_number} has unsupported symbols: {unsupported}"
            )
        words.append(word)
    if not words:
        raise StableWordError("Stable-word library is empty")
    if len(set(words)) != len(words):
        raise StableWordError("Stable-word library contains exact duplicate words")
    return tuple(words)


def encode_degenerate_sequence(sequence: str) -> str:
    """Map an uppercase standard amino-acid sequence to the 12-symbol alphabet."""

    if not sequence:
        raise StableWordError("Cannot encode an empty amino-acid sequence")
    unsupported = sorted(set(sequence) - set(AA_TO_DEGENERATE))
    if unsupported:
        raise StableWordError(f"Sequence contains unsupported residues: {unsupported}")
    return "".join(AA_TO_DEGENERATE[residue] for residue in sequence)


def stable_word_occurrences(
    sequence: str,
    stable_words: Sequence[str],
) -> list[dict[str, object]]:
    """Return all overlapping exact word occurrences with 1-based coordinates."""

    encoded = encode_degenerate_sequence(sequence)
    by_length: dict[int, set[str]] = defaultdict(set)
    for word in stable_words:
        by_length[len(word)].add(word)
    occurrences: list[dict[str, object]] = []
    for length in sorted(by_length):
        if length > len(encoded):
            continue
        words = by_length[length]
        for start_0based in range(len(encoded) - length + 1):
            segment = encoded[start_0based : start_0based + length]
            if segment in words:
                occurrences.append(
                    {
                        "stable_word": segment,
                        "stable_word_length": length,
                        "start_reported_1based": start_0based + 1,
                        "end_reported_1based": start_0based + length,
                        "amino_acid_segment": sequence[start_0based : start_0based + length],
                        "degenerate_segment": segment,
                    }
                )
    return sorted(
        occurrences,
        key=lambda row: (
            int(row["start_reported_1based"]),
            int(row["stable_word_length"]),
            str(row["stable_word"]),
        ),
    )


def stable_word_sequence_features(
    sequence: str,
    stable_words: Sequence[str],
) -> dict[str, object]:
    """Calculate raw and length-normalized stable-word sequence descriptors."""

    occurrences = stable_word_occurrences(sequence, stable_words)
    potential_windows = sum(max(0, len(sequence) - len(word) + 1) for word in stable_words)
    if potential_windows <= 0:
        raise StableWordError("Sequence is shorter than every stable word")
    unique_words = {str(row["stable_word"]) for row in occurrences}
    return {
        "stable_word_occurrence_count": len(occurrences),
        "stable_word_unique_count": len(unique_words),
        "stable_word_potential_window_count": potential_windows,
        "stable_word_occurrences_per_residue": len(occurrences) / len(sequence),
        "stable_word_window_normalized_density": len(occurrences) / potential_windows,
        "stable_word_unique_fraction": len(unique_words) / len(stable_words),
    }


def compare_stable_word_occurrences(
    parent_sequence: str,
    mutant_sequence: str,
    stable_words: Sequence[str],
) -> dict[str, object]:
    """Compare exact degenerate-word occurrences for any same-length mutant.

    Occurrences are keyed by word and coordinates, so overlapping matches are
    retained.  This generic comparison is suitable for single or multiple
    substitutions; it deliberately does not infer that two single-mutant word
    effects add when the combined sequence can create a new joint window.
    """

    if len(parent_sequence) != len(mutant_sequence):
        raise StableWordError("Parent and mutant sequences must have the same length")
    parent = stable_word_occurrences(parent_sequence, stable_words)
    mutant = stable_word_occurrences(mutant_sequence, stable_words)
    parent_by_key = {_occurrence_key(row): row for row in parent}
    mutant_by_key = {_occurrence_key(row): row for row in mutant}
    if len(parent_by_key) != len(parent) or len(mutant_by_key) != len(mutant):
        raise StableWordError("Stable-word occurrences are not uniquely keyed")
    created_keys = sorted(set(mutant_by_key) - set(parent_by_key), key=_key_order)
    lost_keys = sorted(set(parent_by_key) - set(mutant_by_key), key=_key_order)
    created = [dict(mutant_by_key[key]) for key in created_keys]
    lost = [dict(parent_by_key[key]) for key in lost_keys]
    created_count, lost_count = len(created), len(lost)
    if len(mutant) - len(parent) != created_count - lost_count:
        raise StableWordError("Stable-word count reconciliation failed")
    return {
        "wt_stable_word_occurrence_count": len(parent),
        "mutant_stable_word_occurrence_count": len(mutant),
        "created_stable_word_occurrence_count": created_count,
        "lost_stable_word_occurrence_count": lost_count,
        "net_stable_word_occurrence_delta": created_count - lost_count,
        "created_unique_stable_word_count": len({row["stable_word"] for row in created}),
        "lost_unique_stable_word_count": len({row["stable_word"] for row in lost}),
        "longest_created_stable_word_length": max(
            (int(row["stable_word_length"]) for row in created), default=""
        ),
        "longest_lost_stable_word_length": max(
            (int(row["stable_word_length"]) for row in lost), default=""
        ),
        "stable_word_created": bool(created_count),
        "stable_word_lost": bool(lost_count),
        "stable_word_effect": classify_stable_word_effect(created_count, lost_count),
        "stable_word_selection_role": "soft_preference_not_hard_filter",
        "created_occurrences": created,
        "lost_occurrences": lost,
    }


def evaluate_single_mutants(
    parent_sequence: str,
    candidate_rows: Sequence[Mapping[str, object]],
    stable_words: Sequence[str],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Compare stable-word occurrences for validated single mutants against WT.

    Returns one summary row per candidate and one long-form row per created or
    lost occurrence. Candidate identity, reported position, WT residue, mutant
    residue, and the full sequence are all verified before comparison.
    """

    parent_occurrences = stable_word_occurrences(parent_sequence, stable_words)
    parent_by_key = {_occurrence_key(row): row for row in parent_occurrences}
    if len(parent_by_key) != len(parent_occurrences):
        raise StableWordError("WT stable-word occurrences are not uniquely keyed")
    seen_ids: set[str] = set()
    summaries: list[dict[str, object]] = []
    changes: list[dict[str, object]] = []
    for candidate in candidate_rows:
        identifier = str(candidate["candidate_id"])
        if identifier in seen_ids:
            raise StableWordError(f"Duplicate candidate ID: {identifier}")
        seen_ids.add(identifier)
        position = int(candidate["reported_sequence_index_1based"])
        wt = str(candidate["wt_residue"])
        mutant = str(candidate["mutant_residue"])
        sequence = str(candidate["sequence"])
        _validate_single_mutant(parent_sequence, sequence, position, wt, mutant, identifier)
        mutant_occurrences = stable_word_occurrences(sequence, stable_words)
        mutant_by_key = {_occurrence_key(row): row for row in mutant_occurrences}
        created_keys = sorted(set(mutant_by_key) - set(parent_by_key), key=_key_order)
        lost_keys = sorted(set(parent_by_key) - set(mutant_by_key), key=_key_order)
        for change_type, keys, selected in (
            ("created", created_keys, mutant_by_key),
            ("lost", lost_keys, parent_by_key),
        ):
            for key in keys:
                occurrence = selected[key]
                start = int(occurrence["start_reported_1based"])
                end = int(occurrence["end_reported_1based"])
                if not start <= position <= end:
                    raise StableWordError(
                        f"Stable-word change outside mutation window for {identifier}: {key}"
                    )
                offset = slice(start - 1, end)
                changes.append(
                    {
                        "candidate_id": identifier,
                        "reported_sequence_index_1based": position,
                        "mutation_reported_label": str(candidate["mutation_reported_label"]),
                        "change_type": change_type,
                        "stable_word": str(occurrence["stable_word"]),
                        "stable_word_length": int(occurrence["stable_word_length"]),
                        "start_reported_1based": start,
                        "end_reported_1based": end,
                        "overlaps_mutation": True,
                        "wt_amino_acid_segment": parent_sequence[offset],
                        "mutant_amino_acid_segment": sequence[offset],
                        "wt_degenerate_segment": encode_degenerate_sequence(parent_sequence[offset]),
                        "mutant_degenerate_segment": encode_degenerate_sequence(sequence[offset]),
                    }
                )
        created_count, lost_count = len(created_keys), len(lost_keys)
        if len(mutant_occurrences) - len(parent_occurrences) != created_count - lost_count:
            raise StableWordError(f"Stable-word count reconciliation failed for {identifier}")
        summaries.append(
            {
                "candidate_id": identifier,
                "reported_sequence_index_1based": position,
                "wt_residue": wt,
                "mutant_residue": mutant,
                "mutation_reported_label": str(candidate["mutation_reported_label"]),
                "wt_degenerate_symbol": AA_TO_DEGENERATE[wt],
                "mutant_degenerate_symbol": AA_TO_DEGENERATE[mutant],
                "degenerate_symbol_changed": AA_TO_DEGENERATE[wt] != AA_TO_DEGENERATE[mutant],
                "wt_stable_word_occurrence_count": len(parent_occurrences),
                "mutant_stable_word_occurrence_count": len(mutant_occurrences),
                "created_stable_word_occurrence_count": created_count,
                "lost_stable_word_occurrence_count": lost_count,
                "net_stable_word_occurrence_delta": created_count - lost_count,
                "created_unique_stable_word_count": len({key[0] for key in created_keys}),
                "lost_unique_stable_word_count": len({key[0] for key in lost_keys}),
                "longest_created_stable_word_length": max((len(key[0]) for key in created_keys), default=""),
                "longest_lost_stable_word_length": max((len(key[0]) for key in lost_keys), default=""),
                "stable_word_created": bool(created_count),
                "stable_word_lost": bool(lost_count),
                "stable_word_effect": classify_stable_word_effect(created_count, lost_count),
                "stable_word_selection_role": "soft_preference_not_hard_filter",
            }
        )
    return summaries, changes


def classify_stable_word_effect(created_count: int, lost_count: int) -> str:
    """Classify created/lost occurrence counts without hiding exchanges."""

    if created_count < 0 or lost_count < 0:
        raise StableWordError("Stable-word change counts cannot be negative")
    if created_count == 0 and lost_count == 0:
        return "unchanged"
    if created_count > 0 and lost_count == 0:
        return "gain_only"
    if created_count == 0 and lost_count > 0:
        return "loss_only"
    if created_count > lost_count:
        return "net_gain"
    if created_count < lost_count:
        return "net_loss"
    return "balanced_exchange"


def analyze_stable_word_yield(
    sample_rows: Sequence[Mapping[str, object]],
    stable_words: Sequence[str],
) -> dict[str, object]:
    """Evaluate stable-word descriptors against the frozen 47-sequence yield panel."""

    from .nanobert_yield import (
        classify_primary_evidence,
        stratified_bootstrap_ci,
        stratified_permutation_p,
    )
    from .netsolp_yield import yield_metric_row
    from .yield_classification import nested_yield_classification

    if len(sample_rows) != 47 or len({str(row["sample_uid"]) for row in sample_rows}) != 47:
        raise StableWordError("Stable-word yield validation requires 47 unique samples")
    combined: list[dict[str, object]] = []
    for sample in sample_rows:
        sequence = str(sample["sequence_raw"])
        if int(sample["sequence_length_aa"]) != len(sequence):
            raise StableWordError(f"Sequence length mismatch for {sample['sample_uid']}")
        row = dict(sample)
        row.update(stable_word_sequence_features(sequence, stable_words))
        combined.append(row)
    numeric = [row for row in combined if row["observation_semantics"] == "individual_approximate"]
    llj = [row for row in combined if row["provider_code"] == "LLJ"]
    if len(numeric) != 31 or len(llj) != 16:
        raise StableWordError("Frozen yield semantics must contain 31 numeric and 16 LLJ rows")
    metric_rows = [yield_metric_row(numeric, llj, feature) for feature in YIELD_FEATURES]
    primary = next(row for row in metric_rows if row["feature"] == PRIMARY_YIELD_FEATURE)
    low, high = stratified_bootstrap_ci(numeric, PRIMARY_YIELD_FEATURE)
    primary["bootstrap_95ci_low"] = low
    primary["bootstrap_95ci_high"] = high
    primary["stratified_permutation_p"] = stratified_permutation_p(
        numeric, PRIMARY_YIELD_FEATURE, float(primary["stratified_spearman_rho"])
    )
    for metric in metric_rows:
        metric.setdefault("bootstrap_95ci_low", "")
        metric.setdefault("bootstrap_95ci_high", "")
        metric.setdefault("stratified_permutation_p", "")
    evidence_level, reasons = classify_primary_evidence(primary)
    loo = nested_yield_classification(numeric, PRIMARY_YIELD_FEATURE, outer_scheme="leave_one_out")
    cluster = nested_yield_classification(
        numeric, PRIMARY_YIELD_FEATURE, outer_scheme="leave_one_cluster_out"
    )
    return {
        "sample_rows": combined,
        "metric_rows": metric_rows,
        "classification_rows": [
            {"outer_scheme": "leave_one_out", **loo["summary"]},
            {"outer_scheme": "leave_one_cluster_out", **cluster["summary"]},
        ],
        "classification_prediction_rows": loo["prediction_rows"] + cluster["prediction_rows"],
        "primary": primary,
        "empirical_yield_evidence_level": evidence_level,
        "decision_reasons": reasons,
    }


def _occurrence_key(row: Mapping[str, object]) -> tuple[str, int, int]:
    return (
        str(row["stable_word"]),
        int(row["start_reported_1based"]),
        int(row["end_reported_1based"]),
    )


def _key_order(key: tuple[str, int, int]) -> tuple[int, int, str]:
    return key[1], key[2], key[0]


def _validate_single_mutant(
    parent: str,
    sequence: str,
    position: int,
    wt: str,
    mutant: str,
    identifier: str,
) -> None:
    if len(sequence) != len(parent):
        raise StableWordError(f"Candidate length mismatch for {identifier}")
    if not 1 <= position <= len(parent):
        raise StableWordError(f"Candidate position out of range for {identifier}")
    if parent[position - 1] != wt or sequence[position - 1] != mutant or wt == mutant:
        raise StableWordError(f"Candidate mutation identity mismatch for {identifier}")
    differences = [index for index, pair in enumerate(zip(parent, sequence, strict=True), 1) if pair[0] != pair[1]]
    if differences != [position]:
        raise StableWordError(f"Candidate is not the declared single mutant: {identifier}")
