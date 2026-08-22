#!/usr/bin/env python3
"""Merge double-mutant property scores without selecting the final 11."""

from __future__ import annotations

import argparse
import csv
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

from antibody_optimization.expression_double_mutants import merge_property_scores  # noqa: E402
from antibody_optimization.file_transaction import replace_staged_files, validate_file_paths  # noqa: E402

NAMES = {
    "matrix": "expression_double_mutant_property_matrix.csv",
    "bands": "expression_double_mutant_magnitude_band_counts.csv",
    "plot": "expression_double_mutant_property_overview.png",
    "gate": "expression_double_mutant_property_matrix_gate.json",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--parent19-dir", type=Path, required=True)
    parser.add_argument("--netsolp-score-dir", type=Path, required=True)
    parser.add_argument("--nanomelt-score-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    plan = args.plan_dir.resolve(strict=True)
    parent = args.parent19_dir.resolve(strict=True)
    netsolp = args.netsolp_score_dir.resolve(strict=True)
    nanomelt = args.nanomelt_score_dir.resolve(strict=True)
    sources = [
        plan / "expression_double_mutant_candidates.csv",
        plan / "expression_double_mutant_plan_gate.json",
        parent / "expression_single_mutant_parent19.csv",
        netsolp / "netsolp_sample_scores.csv",
        netsolp / "netsolp_model_run.json",
        nanomelt / "nanomelt_sample_scores.csv",
        nanomelt / "nanomelt_model_run.json",
    ]
    for run_path in (sources[4], sources[6]):
        if _json(run_path).get("status") != "pass":
            raise ValueError(f"Tool run did not pass: {run_path}")
    if _json(sources[1]).get("release") != "ready_for_netsolp_nanomelt_double_scoring":
        raise ValueError("Double-mutant plan gate is not released")
    rows = merge_property_scores(_csv(sources[0]), _csv(sources[2]), _csv(sources[3]), _csv(sources[5]))
    band_rows = _band_counts(rows)
    output = args.output_dir.absolute()
    output.mkdir(parents=True, exist_ok=True)
    args.run_summary.parent.mkdir(parents=True, exist_ok=True)
    targets = [*(output / name for name in NAMES.values()), args.run_summary.absolute()]
    validated = validate_file_paths(project_root=ROOT, source_paths=sources, target_paths=targets)
    if existing := [path for path in validated.target_paths if path.exists()]:
        raise FileExistsError("Refusing to overwrite:\n" + "\n".join(map(str, existing)))
    gate = {
        "schema_version": 1,
        "generated_at": generated_at,
        "status": "pass",
        "release": "ready_for_magnitude_aware_double_mutant_selection",
        "candidate_count": len(rows),
        "netsolp_complete_count": len(rows),
        "nanomelt_complete_count": sum(row["nanomelt_scoring_status"] == "pass" for row in rows),
        "antifold_same_view_additive_evaluable_count": sum(
            str(row["antifold_same_view_additive_evaluable"]).lower() == "true" for row in rows
        ),
        "antifold_mixed_view_component_only_count": sum(
            str(row["antifold_same_view_additive_evaluable"]).lower() != "true" for row in rows
        ),
        "candidate_selection_performed": False,
        "final_11_double_mutants_selected": False,
        "interpretation": (
            "All 162 doubles have complete NetSolP/NanoMelt scores and retained constituent "
            "AntiFold evidence. Interaction residuals describe predictor non-additivity, not physical epistasis."
        ),
    }
    finals = dict(zip(NAMES, validated.target_paths[:-1], strict=True))
    run_summary = validated.target_paths[-1]
    with tempfile.TemporaryDirectory(prefix=".expression-double-final-", dir=ROOT) as temp:
        stage = Path(temp)
        staged = {key: stage / name for key, name in NAMES.items()}
        staged_run = stage / "run_summary.json"
        _write_csv(staged["matrix"], rows)
        _write_csv(staged["bands"], band_rows)
        _plot(staged["plot"], rows)
        _write_json(staged["gate"], gate)
        _write_json(
            staged_run,
            {
                "schema_version": 1,
                "status": "pass",
                "generated_at": generated_at,
                "elapsed_seconds": round(time.perf_counter() - started, 6),
                "python": platform.python_version(),
                "command_argv": [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
                "candidate_count": len(rows),
                "candidate_selection_performed": False,
            },
        )
        replace_staged_files(
            {**{staged[key]: finals[key] for key in NAMES}, staged_run: run_summary},
            project_root=ROOT,
            protected_source_paths=validated.source_paths,
        )
    return 0


def _band_counts(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    fields = (
        "netsolp_u_magnitude_band",
        "netsolp_s_magnitude_band",
        "nanomelt_tm_c_magnitude_band",
        "antifold_worst_component_band",
    )
    output = []
    for field in fields:
        for band, count in sorted(Counter(str(row[field]) for row in rows).items()):
            output.append({"metric": field, "magnitude_band": band, "candidate_count": count})
    return output


def _plot(path: Path, rows: list[dict[str, object]]) -> None:
    import matplotlib.pyplot as plt

    metrics = (
        ("netsolp_u_delta_vs_wt", "NetSolP ΔU"),
        ("netsolp_s_delta_vs_wt", "NetSolP ΔS"),
        ("nanomelt_tm_c_delta_vs_wt", "NanoMelt ΔTm (°C)"),
    )
    figure, axes = plt.subplots(1, 3, figsize=(10.5, 3.4))
    for axis, (field, label) in zip(axes, metrics, strict=True):
        values = [float(row[field]) for row in rows]
        axis.hist(values, bins=18, color="#4f92bd", edgecolor="white")
        axis.axvline(0, color="#333333", linewidth=1)
        axis.set_xlabel(label)
        axis.set_ylabel("Double mutants")
    figure.suptitle("Nb252 162-double property landscape (no selection)")
    figure.tight_layout()
    figure.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    raise SystemExit(main())
