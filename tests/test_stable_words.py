from __future__ import annotations

import csv
from pathlib import Path

import pytest

from antibody_optimization.stable_words import (
    StableWordError,
    analyze_stable_word_yield,
    encode_degenerate_sequence,
    evaluate_single_mutants,
    parse_stable_words,
    stable_word_occurrences,
)


ROOT = Path(__file__).resolve().parents[1]
WORDS = ROOT / "data/Stable_word_SS_3D_1336 (1).txt"
MATRIX = (
    ROOT / "docs/result_artifacts/candidate_design/expression_property_complete_matrix_v2_20260819"
    / "expression_single_mutant_property_matrix.csv"
)
YIELD_SAMPLES = (
    ROOT / "docs/result_artifacts/candidate_design/netsolp_yield_validation_plan_20260814"
    / "netsolp_validation_samples.csv"
)
PARENT = "QVQLQESGGGLVQAGGSLRLSCAASGTIFFGYDMGWYRQAPGKEREFVASITTGSNTNYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAVDTIDYIIEWNVYYYIFSYWGQGTQVTVSSGS"


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _words() -> tuple[str, ...]:
    return parse_stable_words(WORDS.read_text(encoding="utf-8-sig").splitlines(keepends=True))


def test_degenerate_encoding_uses_exact_case_sensitive_twelve_symbol_contract():
    assert encode_degenerate_sequence("DEKRLIVNQSTAGFWCHMPY") == "aabbhhhnnooss@@CHMPY"
    assert parse_stable_words(["h\n", "H\n"]) == ("h", "H")
    with pytest.raises(StableWordError, match="unsupported symbols"):
        parse_stable_words(["Aaa\n"])


def test_occurrence_scan_counts_overlaps_and_nested_words():
    rows = stable_word_occurrences("AAA", ("ss", "sss"))
    assert [(row["stable_word"], row["start_reported_1based"]) for row in rows] == [
        ("ss", 1), ("sss", 1), ("ss", 2),
    ]


def test_same_degenerate_class_is_unchanged_and_exchange_is_not_hidden():
    same = [{
        "candidate_id": "A1G", "reported_sequence_index_1based": 1,
        "wt_residue": "A", "mutant_residue": "G", "mutation_reported_label": "A1G",
        "sequence": "GAAAA",
    }]
    summary, changes = evaluate_single_mutants("AAAAA", same, ("sssss",))
    assert summary[0]["stable_word_effect"] == "unchanged"
    assert changes == []

    exchange = [{
        "candidate_id": "A1F", "reported_sequence_index_1based": 1,
        "wt_residue": "A", "mutant_residue": "F", "mutation_reported_label": "A1F",
        "sequence": "FAAAA",
    }]
    summary, changes = evaluate_single_mutants("AAAAA", exchange, ("sssss", "@ssss"))
    assert summary[0]["stable_word_effect"] == "balanced_exchange"
    assert {row["change_type"] for row in changes} == {"created", "lost"}
    assert all(row["overlaps_mutation"] for row in changes)


def test_real_library_and_847_single_mutants_reconcile():
    words = _words()
    assert len(words) == len(set(words)) == 1336
    summaries, changes = evaluate_single_mutants(PARENT, _csv(MATRIX), words)
    assert len(summaries) == 847
    assert len({row["candidate_id"] for row in summaries}) == 847
    assert all(row["overlaps_mutation"] for row in changes)
    assert all(
        int(row["mutant_stable_word_occurrence_count"])
        - int(row["wt_stable_word_occurrence_count"])
        == int(row["net_stable_word_occurrence_delta"])
        for row in summaries
    )
    assert all(
        row["stable_word_effect"] == "unchanged"
        for row in summaries if not row["degenerate_symbol_changed"]
    )


def test_real_47_sequence_validation_preserves_frozen_semantics():
    result = analyze_stable_word_yield(_csv(YIELD_SAMPLES), _words())
    assert len(result["sample_rows"]) == 47
    assert len(result["metric_rows"]) == 4
    assert len({tuple(row) for row in result["metric_rows"]}) == 1
    assert {row["outer_scheme"] for row in result["classification_rows"]} == {
        "leave_one_out", "leave_one_cluster_out",
    }
    assert len(result["classification_prediction_rows"]) == 62
    assert result["empirical_yield_evidence_level"] in {
        "weak_ranking_evidence", "compatibility_filter_only", "no_supported_use",
    }
