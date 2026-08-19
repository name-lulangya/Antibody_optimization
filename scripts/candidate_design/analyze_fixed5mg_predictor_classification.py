#!/usr/bin/env python3
"""Render a display-only RP3Net/NetSolP comparison at a fixed 5 mg cutoff."""

from __future__ import annotations

import argparse
import csv
import json
import sys
import tempfile
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from antibody_optimization.file_transaction import replace_staged_files, validate_file_paths  # noqa: E402
from antibody_optimization.fixed_yield_plot import render_fixed_yield_classification_figure  # noqa: E402
from antibody_optimization.yield_classification import (  # noqa: E402
    fixed_yield_apparent_classification,
    fixed_yield_nested_classification,
)


FEATURES = {
    "rp3net_expression_probability": "RP3Net",
    "predicted_usability": "NetSolP U",
    "predicted_solubility": "NetSolP S",
}
NAMES = {
    "metrics": "fixed5mg_classification_metrics.csv",
    "predictions": "fixed5mg_classification_predictions.csv",
    "contract": "fixed5mg_classification_contract.json",
    "png": "fixed5mg_classification.png",
    "svg": "fixed5mg_classification.svg",
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rp3net-samples", type=Path, required=True)
    parser.add_argument("--netsolp-samples", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--yield-threshold", type=float, default=5.0)
    parser.add_argument("--generated-at")
    args = parser.parse_args()
    generated = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    rp3net_path = args.rp3net_samples.resolve(strict=True)
    netsolp_path = args.netsolp_samples.resolve(strict=True)
    combined = _join_numeric_rows(_csv(rp3net_path), _csv(netsolp_path))

    metric_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    for feature, display in FEATURES.items():
        for scheme in ("leave_one_out", "leave_one_cluster_out"):
            result = fixed_yield_nested_classification(
                combined, feature, outer_scheme=scheme, yield_threshold=args.yield_threshold,
            )
            metric_rows.append({"feature": feature, "predictor": display, "outer_scheme": scheme, **result["summary"]})
            prediction_rows.extend({"predictor": display, **row} for row in result["prediction_rows"])
        apparent = fixed_yield_apparent_classification(combined, feature, yield_threshold=args.yield_threshold)
        metric_rows.append({"feature": feature, "predictor": display, "outer_scheme": "apparent_full_sample", **apparent["summary"]})
        prediction_rows.extend({"predictor": display, **row} for row in apparent["prediction_rows"])

    output_dir = args.output_dir.absolute()
    output_dir.mkdir(parents=True, exist_ok=True)
    targets = [output_dir / name for name in NAMES.values()] + [args.run_summary.absolute()]
    valid = validate_file_paths(
        project_root=ROOT, source_paths=[rp3net_path, netsolp_path], target_paths=targets,
    )
    if any(path.exists() for path in valid.target_paths):
        raise FileExistsError("Refusing to overwrite fixed-10-mg classification outputs")
    for path in valid.target_paths:
        path.parent.mkdir(parents=True, exist_ok=True)
    final = dict(zip((*NAMES, "summary"), valid.target_paths, strict=True))
    high_count = sum(float(row["numeric_yield_value"]) >= args.yield_threshold for row in combined)
    contract = {
        "schema_version": 1, "status": "pass", "generated_at": generated,
        "numeric_sample_count": len(combined), "high_yield_count": high_count,
        "low_yield_count": len(combined) - high_count,
        "yield_label_rule": f"numeric_yield_value >= {args.yield_threshold:g}",
        "yield_threshold": args.yield_threshold, "yield_unit": "reported_source_value_mg_per_1L",
        "included_semantics": "individual_approximate", "excluded_semantics": "LLJ_group_level_censored",
        "predictors": FEATURES, "predictor_direction": "higher_is_better",
        "score_threshold_rule": "maximize_training_MCC_then_balanced_accuracy_then_higher_threshold",
        "decision_use": "display_only_not_a_predictor_gate_or_candidate_filter",
        "primary_evaluation": "none_display_only",
        "apparent_full_sample_warning": "Threshold fitted and evaluated on the same 31 rows; descriptive only.",
    }
    with tempfile.TemporaryDirectory(prefix=".fixed5mg-", dir=ROOT) as temp:
        stage = Path(temp)
        staged = {key: stage / Path(path).name for key, path in final.items()}
        _write_csv(staged["metrics"], metric_rows)
        _write_csv(staged["predictions"], prediction_rows)
        _write_json(staged["contract"], contract)
        render_fixed_yield_classification_figure(
            combined, metric_rows, yield_threshold=args.yield_threshold,
            png_path=staged["png"], svg_path=staged["svg"],
        )
        _write_json(staged["summary"], {
            "schema_version": 1, "status": "pass", "generated_at": generated,
            "counts": {"numeric_samples": len(combined), "high": high_count, "low": len(combined) - high_count},
            "outputs": {key: str(path) for key, path in final.items() if key != "summary"},
        })
        replace_staged_files(
            {staged[key]: final[key] for key in staged}, project_root=ROOT,
            protected_source_paths=valid.source_paths,
        )
    return 0


def _join_numeric_rows(rp3net_rows: list[dict[str, str]], netsolp_rows: list[dict[str, str]]) -> list[dict[str, object]]:
    netsolp = {row["sample_uid"]: row for row in netsolp_rows}
    if len(netsolp) != 47 or len(rp3net_rows) != 47 or {row["sample_uid"] for row in rp3net_rows} != set(netsolp):
        raise ValueError("RP3Net and NetSolP inputs must contain the same 47 unique samples")
    combined = []
    for row in rp3net_rows:
        other = netsolp[row["sample_uid"]]
        for key in ("sequence_raw", "provider_code", "observation_semantics", "numeric_yield_value", "sequence_cluster_90"):
            if row[key] != other[key]:
                raise ValueError(f"Input mismatch for {row['sample_uid']}: {key}")
        if row["observation_semantics"] != "individual_approximate":
            continue
        item: dict[str, object] = dict(row)
        item["predicted_usability"] = float(other["predicted_usability"])
        item["predicted_solubility"] = float(other["predicted_solubility"])
        combined.append(item)
    if len(combined) != 31:
        raise ValueError("Expected exactly 31 individual approximate numeric observations")
    return combined


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def _write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    raise SystemExit(main())
