#!/usr/bin/env python3
"""Validate stable-word descriptors against the frozen 47-sequence BL21 panel."""

from __future__ import annotations

import argparse
import csv
import json
import platform
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from antibody_optimization.file_transaction import replace_staged_files, validate_file_paths  # noqa: E402
from antibody_optimization.stable_words import (  # noqa: E402
    PRIMARY_YIELD_FEATURE,
    analyze_stable_word_yield,
    parse_stable_words,
)


NAMES = {
    "samples": "stable_word_yield_sample_evidence.csv",
    "metrics": "stable_word_yield_associations.csv",
    "classification": "stable_word_yield_classification.csv",
    "predictions": "stable_word_yield_classification_predictions.csv",
    "gate": "stable_word_yield_validation_gate.json",
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stable-word-result-dir", type=Path, required=True)
    parser.add_argument("--yield-samples", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    args = parser.parse_args()
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    stable_dir = args.stable_word_result_dir.resolve(strict=True)
    sample_path = args.yield_samples.resolve(strict=True)
    library_path = stable_dir / "stable_word_library.csv"
    contract_path = stable_dir / "stable_word_evaluation_contract.json"
    feature_gate_path = stable_dir / "stable_word_single_mutant_gate.json"
    sources = [library_path, contract_path, feature_gate_path, sample_path]
    contract, feature_gate = _json(contract_path), _json(feature_gate_path)
    if contract.get("status") != "pass" or feature_gate.get("status") != "pass":
        raise ValueError("Stable-word feature contract is not released")
    library_rows = _csv(library_path)
    stable_words = parse_stable_words([str(row["stable_word"]) + "\n" for row in library_rows])
    if len(stable_words) != int(contract["source_line_count"]):
        raise ValueError("Normalized stable-word library does not match its frozen contract")
    result = analyze_stable_word_yield(_csv(sample_path), stable_words)
    gate = {
        "schema_version": 1,
        "status": "pass",
        "generated_at": generated_at,
        "gate_name": "nb252_stable_word_bl21_reported_yield_validation_v1",
        "sample_count": 47,
        "numeric_individual_count": 31,
        "llj_ordinal_censored_count": 16,
        "primary_feature": PRIMARY_YIELD_FEATURE,
        "higher_is_predeclared_as_preferred": True,
        "continuous_statistics": result["primary"],
        "classification_statistics": result["classification_rows"],
        "empirical_yield_evidence_level": result["empirical_yield_evidence_level"],
        "decision_reasons": result["decision_reasons"],
        "candidate_selection_role": "user_directed_soft_tie_breaker_not_hard_filter",
        "candidate_selection_performed": False,
        "high_capacity_model_trained": False,
        "release": "stable_word_feature_ready_as_interpretable_soft_preference_only",
        "interpretation": (
            "Stable-word descriptors are external sequence-pattern priors, not measured stability, "
            "BL21 yield, or a causal expression model."
        ),
    }
    output_dir = args.output_dir.absolute()
    summary_path = args.run_summary.absolute()
    targets = [output_dir / name for name in NAMES.values()] + [summary_path]
    valid = validate_file_paths(project_root=ROOT, source_paths=sources, target_paths=targets)
    if any(path.exists() for path in valid.target_paths):
        raise FileExistsError("Refusing to overwrite stable-word yield-validation outputs")
    for path in valid.target_paths:
        path.parent.mkdir(parents=True, exist_ok=True)
    final = dict(zip((*NAMES, "run_summary"), valid.target_paths, strict=True))
    with tempfile.TemporaryDirectory(prefix=".stable-word-yield-", dir=ROOT) as temp_name:
        stage = Path(temp_name)
        staged = {key: stage / path.name for key, path in final.items()}
        _write_csv(staged["samples"], result["sample_rows"])
        _write_csv(staged["metrics"], result["metric_rows"])
        _write_csv(staged["classification"], result["classification_rows"])
        _write_csv(staged["predictions"], result["classification_prediction_rows"])
        _write_json(staged["gate"], gate)
        _write_json(staged["run_summary"], {
            "schema_version": 1,
            "status": "pass",
            "generated_at": generated_at,
            "python": platform.python_version(),
            "elapsed_seconds": round(time.perf_counter() - started, 6),
            "empirical_yield_evidence_level": result["empirical_yield_evidence_level"],
            "counts": {
                "samples": len(result["sample_rows"]),
                "metrics": len(result["metric_rows"]),
                "classification_predictions": len(result["classification_prediction_rows"]),
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
