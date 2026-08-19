#!/usr/bin/env python3
"""Add stable-word occurrence changes to the complete 847-mutant property matrix."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import sys
import tempfile
import time
from collections import Counter
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from antibody_optimization.file_transaction import replace_staged_files, validate_file_paths  # noqa: E402
from antibody_optimization.stable_words import (  # noqa: E402
    AA_TO_DEGENERATE,
    STABLE_WORD_ALPHABET,
    evaluate_single_mutants,
    parse_stable_words,
)


NAMES = {
    "library": "stable_word_library.csv",
    "summary": "stable_word_single_mutant_summary.csv",
    "changes": "stable_word_occurrence_changes.csv",
    "matrix": "expression_single_mutant_property_stable_word_matrix.csv",
    "contract": "stable_word_evaluation_contract.json",
    "gate": "stable_word_single_mutant_gate.json",
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stable-words", type=Path, required=True)
    parser.add_argument("--property-matrix-dir", type=Path, required=True)
    parser.add_argument("--completion-plan-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    args = parser.parse_args()
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    stable_path = args.stable_words.resolve(strict=True)
    matrix_dir = args.property_matrix_dir.resolve(strict=True)
    plan_dir = args.completion_plan_dir.resolve(strict=True)
    matrix_path = matrix_dir / "expression_single_mutant_property_matrix.csv"
    matrix_gate_path = matrix_dir / "expression_single_mutant_property_matrix_gate.json"
    completion_contract_path = plan_dir / "expression_property_completion_contract.json"
    sources = [stable_path, matrix_path, matrix_gate_path, completion_contract_path]
    matrix_gate = _json(matrix_gate_path)
    completion_contract = _json(completion_contract_path)
    if matrix_gate.get("status") != "pass" or int(matrix_gate.get("candidate_count", 0)) != 847:
        raise ValueError("Complete property matrix gate is not released for 847 candidates")
    stable_bytes = stable_path.read_bytes()
    stable_words = parse_stable_words(stable_bytes.decode("utf-8-sig").splitlines(keepends=True))
    matrix_rows = _csv(matrix_path)
    if len(matrix_rows) != 847:
        raise ValueError("Stable-word evaluation requires the complete 847-row property matrix")
    parent = str(completion_contract["authoritative_parent"]["sequence"])
    summaries, changes = evaluate_single_mutants(parent, matrix_rows, stable_words)
    summary_by_id = {str(row["candidate_id"]): row for row in summaries}
    augmented = []
    identity_fields = {
        "candidate_id", "reported_sequence_index_1based", "wt_residue",
        "mutant_residue", "mutation_reported_label",
    }
    for row in matrix_rows:
        combined = dict(row)
        combined.update({
            key: value for key, value in summary_by_id[str(row["candidate_id"])].items()
            if key not in identity_fields
        })
        augmented.append(combined)
    library_rows = [
        {"stable_word_id": f"stable_word_{index:04d}", "stable_word": word, "stable_word_length": len(word)}
        for index, word in enumerate(stable_words, 1)
    ]
    effect_counts = Counter(str(row["stable_word_effect"]) for row in summaries)
    contract = {
        "schema_version": 1,
        "status": "pass",
        "generated_at": generated_at,
        "contract_name": "nb252_degenerate_stable_word_evaluation_v1",
        "source_path": str(stable_path.relative_to(ROOT)),
        "source_sha256": hashlib.sha256(stable_bytes).hexdigest().upper(),
        "source_line_count": len(stable_words),
        "case_sensitive": True,
        "alphabet": sorted(STABLE_WORD_ALPHABET),
        "amino_acid_to_symbol": AA_TO_DEGENERATE,
        "match_semantics": {
            "algorithm": "full_exact_contiguous_substring_scan_after_deterministic_degenerate_encoding",
            "overlapping_occurrences_counted": True,
            "nested_and_different_length_words_kept_separate": True,
            "occurrence_key": ["stable_word", "start_reported_1based", "end_reported_1based"],
            "created": "mutant_occurrences_minus_WT_occurrences",
            "lost": "WT_occurrences_minus_mutant_occurrences",
            "net_delta": "created_occurrence_count_minus_lost_occurrence_count",
        },
        "selection_semantics": {
            "role": "soft_preference_not_hard_filter",
            "preference_order": ["gain_only", "net_gain", "unchanged", "balanced_exchange", "net_loss", "loss_only"],
            "cannot_override": ["immutable_positions", "new_cysteine_exclusion", "clear_property_regression"],
            "causal_stability_or_expression_claim": False,
        },
        "candidate_count": len(summaries),
        "candidate_selection_performed": False,
        "ranking_performed": False,
    }
    gate = {
        "schema_version": 1,
        "status": "pass",
        "generated_at": generated_at,
        "gate_name": "nb252_stable_word_single_mutant_features_v1",
        "candidate_count": len(summaries),
        "occurrence_change_row_count": len(changes),
        "stable_word_effect_counts": dict(sorted(effect_counts.items())),
        "same_degenerate_symbol_count": sum(not bool(row["degenerate_symbol_changed"]) for row in summaries),
        "all_occurrence_changes_overlap_mutation": all(bool(row["overlaps_mutation"]) for row in changes),
        "candidate_selection_performed": False,
        "ranking_performed": False,
        "release": "stable_word_features_ready_for_yield_validation_and_later_soft_preference",
    }
    output_dir = args.output_dir.absolute()
    summary_path = args.run_summary.absolute()
    targets = [output_dir / name for name in NAMES.values()] + [summary_path]
    valid = validate_file_paths(project_root=ROOT, source_paths=sources, target_paths=targets)
    if any(path.exists() for path in valid.target_paths):
        raise FileExistsError("Refusing to overwrite stable-word single-mutant outputs")
    for path in valid.target_paths:
        path.parent.mkdir(parents=True, exist_ok=True)
    final = dict(zip((*NAMES, "run_summary"), valid.target_paths, strict=True))
    with tempfile.TemporaryDirectory(prefix=".stable-word-mutants-", dir=ROOT) as temp_name:
        stage = Path(temp_name)
        staged = {key: stage / path.name for key, path in final.items()}
        _write_csv(staged["library"], library_rows)
        _write_csv(staged["summary"], summaries)
        _write_csv(staged["changes"], changes)
        _write_csv(staged["matrix"], augmented)
        _write_json(staged["contract"], contract)
        _write_json(staged["gate"], gate)
        _write_json(staged["run_summary"], {
            "schema_version": 1,
            "status": "pass",
            "generated_at": generated_at,
            "python": platform.python_version(),
            "elapsed_seconds": round(time.perf_counter() - started, 6),
            "counts": {
                "stable_words": len(stable_words),
                "candidates": len(summaries),
                "occurrence_changes": len(changes),
            },
            "candidate_selection_performed": False,
            "outputs": {key: str(value) for key, value in final.items() if key != "run_summary"},
        })
        replace_staged_files(
            {staged[key]: final[key] for key in staged},
            project_root=ROOT,
            protected_source_paths=valid.source_paths,
        )
    return 0


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Cannot write empty CSV: {path.name}")
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def _write_json(path: Path, value: object) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8", newline="\n",
    )


if __name__ == "__main__":
    raise SystemExit(main())
