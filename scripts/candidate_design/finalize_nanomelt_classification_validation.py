#!/usr/bin/env python3
"""Add leakage-controlled classification to the released NanoMelt evidence."""

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
from antibody_optimization.nanomelt_yield import analyze_nanomelt_associations  # noqa: E402
from antibody_optimization.nanomelt_yield_plot import render_nanomelt_yield_figure  # noqa: E402


NAMES = {
    "classification": "nanomelt_yield_classification.csv",
    "predictions": "nanomelt_yield_classification_predictions.csv",
    "gate": "nanomelt_yield_classification_gate.json",
    "png": "nanomelt_yield_classification.png",
    "svg": "nanomelt_yield_classification.svg",
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--yield-samples", type=Path, required=True)
    parser.add_argument("--prior-sample-evidence", type=Path, required=True)
    parser.add_argument("--prior-gate", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    args = parser.parse_args()
    started = time.perf_counter()
    generated = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    sources = [
        args.yield_samples.resolve(strict=True),
        args.prior_sample_evidence.resolve(strict=True),
        args.prior_gate.resolve(strict=True),
    ]
    prior_gate = _json(sources[2])
    if prior_gate.get("status") != "pass" or prior_gate.get("evidence_level") != "no_supported_use":
        raise ValueError("Prior NanoMelt evidence gate is not the released no-supported-use result")
    result = analyze_nanomelt_associations(_csv(sources[0]), _csv(sources[1]))
    gate = {
        "schema_version": 1,
        "gate_name": "nanomelt_bl21_yield_continuous_and_classification_decision",
        "status": "pass",
        "generated_at": generated,
        "planned_sample_count": 47,
        "scored_sample_count": int(result["primary"]["scored_n"]),
        "numeric_classification_sample_count": int(result["primary"]["numeric_n"]),
        "not_scored_sample_uids": result["not_scored_uids"],
        "feature": "nanomelt_predicted_apparent_tm_c",
        "expected_direction": "higher_predicted_apparent_tm_higher_reported_yield",
        "classification_label_definition": (
            "high_vs_low_relative_to_matching_provider_median_fitted_inside_each_outer_training_fold"
        ),
        "score_threshold_selection": (
            "maximize_training_MCC_then_balanced_accuracy_then_higher_threshold"
        ),
        "classification_results": result["classification_rows"],
        "continuous_evidence_level": result["primary"]["continuous_evidence_level"],
        "combined_evidence_level": result["evidence_level"],
        "decision_reasons": result["decision_reasons"],
        "yield_ranking_supported": False,
        "selection_role": "predicted_stability_compatibility_constraint_only",
        "release": "nanomelt_not_supported_for_yield_ranking",
        "interpretation": (
            "Predicted apparent Tm is retained as a stability constraint; the current data do not "
            "support treating it as BL21-yield ranker or classifier."
        ),
    }
    targets = [args.output_dir.absolute() / name for name in NAMES.values()] + [
        args.run_summary.absolute()
    ]
    valid = validate_file_paths(
        project_root=ROOT, source_paths=sources, target_paths=targets
    )
    if any(path.exists() for path in valid.target_paths):
        raise FileExistsError("Refusing to overwrite NanoMelt classification outputs")
    for path in valid.target_paths:
        path.parent.mkdir(parents=True, exist_ok=True)
    final = dict(zip((*NAMES, "summary"), valid.target_paths, strict=True))
    with tempfile.TemporaryDirectory(prefix=".nanomelt-classification-", dir=ROOT) as temp:
        stage = Path(temp)
        staged = {key: stage / Path(value).name for key, value in final.items()}
        _write_csv(staged["classification"], result["classification_rows"])
        _write_csv(staged["predictions"], result["classification_prediction_rows"])
        _write_json(staged["gate"], gate)
        render_nanomelt_yield_figure(
            result["sample_rows"],
            result["primary"],
            result["cv_rows"],
            result["leave_one_out_rows"],
            result["classification_rows"],
            result["classification_prediction_rows"],
            png_path=staged["png"],
            svg_path=staged["svg"],
        )
        _write_json(
            staged["summary"],
            {
                "schema_version": 1,
                "status": "pass",
                "generated_at": generated,
                "python": platform.python_version(),
                "elapsed_seconds": round(time.perf_counter() - started, 6),
                "combined_evidence_level": result["evidence_level"],
                "outputs": {key: str(value) for key, value in final.items() if key != "summary"},
            },
        )
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
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, value: object) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")


if __name__ == "__main__":
    raise SystemExit(main())
