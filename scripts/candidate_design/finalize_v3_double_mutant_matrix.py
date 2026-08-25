#!/usr/bin/env python3
"""Join complete V3 double NetSolP/NanoMelt scores without selecting 15.

The workflow requires 103/103 passing rows from both tools, calculates
complete-sequence changes relative to the same WT, preserves frozen magnitude
bands and model non-additivity diagnostics, and releases a 102-row matrix for
expert review and a later explicit final-double selection stage.
"""

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

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from antibody_optimization.file_transaction import (  # noqa: E402
    replace_staged_files,
    validate_file_paths,
)
from antibody_optimization.v3_double_mutant_result_plot import (  # noqa: E402
    build_v3_double_result_plot_rows,
    render_v3_double_mutant_results,
)
from antibody_optimization.v3_double_mutants import (  # noqa: E402
    EXPECTED_SCORE_SAMPLE_COUNT,
    EXPECTED_VALID_DOUBLE_COUNT,
    merge_v3_property_scores,
)


DEFAULT_PLAN_DIR = (
    ROOT
    / "docs/result_artifacts/candidate_design"
    / "v3_double_mutant_plan_20260825"
)
DEFAULT_SCORE_ROOT = (
    ROOT / "results/candidate_design/v3_double_mutant_scan_20260825"
)
DEFAULT_OUTPUT_DIR = (
    ROOT
    / "docs/result_artifacts/candidate_design"
    / "v3_double_mutant_property_matrix_20260825"
)
DEFAULT_RUN_SUMMARY = (
    ROOT
    / "docs/run_summaries/candidate_design"
    / "v3_double_mutant_property_matrix_20260825"
    / "run_summary.json"
)
OUTPUT_NAMES = {
    "matrix": "v3_double_mutant_property_matrix102.csv",
    "bands": "v3_double_mutant_magnitude_band_counts.csv",
    "nonadditivity": "v3_double_mutant_model_nonadditivity_summary.csv",
    "review": "v3_double_mutant_review_priority.csv",
    "plot_data": "v3_double_mutant_property_plot_data.csv",
    "png": "v3_double_mutant_property_overview.png",
    "svg": "v3_double_mutant_property_overview.svg",
    "manifest": "v3_double_mutant_property_matrix_manifest.json",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN_DIR)
    parser.add_argument("--netsolp-score-dir", type=Path)
    parser.add_argument("--nanomelt-score-dir", type=Path)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--run-summary", type=Path, default=DEFAULT_RUN_SUMMARY)
    parser.add_argument("--generated-at", default="")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(
        timespec="seconds"
    )
    plan_dir = args.plan_dir.resolve(strict=True)
    score_root = DEFAULT_SCORE_ROOT
    netsolp_dir = (
        args.netsolp_score_dir or score_root / "netsolp"
    ).resolve(strict=True)
    nanomelt_dir = (
        args.nanomelt_score_dir or score_root / "nanomelt"
    ).resolve(strict=True)
    candidate_path = plan_dir / "v3_double_mutant_candidates102.csv"
    plan_manifest_path = plan_dir / "v3_double_mutant_plan_manifest.json"
    net_score_path = netsolp_dir / "netsolp_sample_scores.csv"
    net_run_path = netsolp_dir / "netsolp_model_run.json"
    melt_score_path = nanomelt_dir / "nanomelt_sample_scores.csv"
    melt_run_path = nanomelt_dir / "nanomelt_model_run.json"
    sources = (
        candidate_path,
        plan_manifest_path,
        net_score_path,
        net_run_path,
        melt_score_path,
        melt_run_path,
    )
    plan_manifest = _json(plan_manifest_path)
    _verify_tool_runs(plan_manifest, _json(net_run_path), _json(melt_run_path))
    net_rows = _csv(net_score_path)
    melt_rows = _csv(melt_score_path)
    invalid_nanomelt_domain = [
        row.get("sample_uid", "")
        for row in melt_rows
        if not (
            row.get("scoring_status") == "pass"
            and int(row.get("scored_length_aa", 0)) == 126
            and row.get("trimmed_n_terminal", "") == ""
            and row.get("trimmed_c_terminal", "") == "GS"
        )
    ]
    if invalid_nanomelt_domain:
        raise ValueError(
            "NanoMelt scored-domain contract failed for: "
            + ",".join(invalid_nanomelt_domain)
        )
    rows = merge_v3_property_scores(
        _csv(candidate_path), net_rows, melt_rows
    )
    rows = [_annotate_review_class(row) for row in rows]
    band_rows = _band_counts(rows)
    nonadditivity_rows = _nonadditivity_summary(rows)
    review_rows = _review_rows(rows)
    plot_rows = build_v3_double_result_plot_rows(rows)

    output_dir = args.output_dir.resolve()
    run_summary = args.run_summary.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    run_summary.parent.mkdir(parents=True, exist_ok=True)
    finals = {key: output_dir / name for key, name in OUTPUT_NAMES.items()}
    existing = [path for path in (*finals.values(), run_summary) if path.exists()]
    if existing and not args.overwrite:
        raise FileExistsError(
            "Refusing to overwrite V3 double property outputs:\n"
            + "\n".join(str(path) for path in existing)
        )
    validated = validate_file_paths(
        project_root=ROOT,
        source_paths=sources,
        target_paths=(*finals.values(), run_summary),
    )
    facts = _facts(rows)
    with tempfile.TemporaryDirectory(prefix=".v3-double-matrix-", dir=ROOT) as temp:
        stage = Path(temp)
        staged = {key: stage / name for key, name in OUTPUT_NAMES.items()}
        staged_run = stage / "run_summary.json"
        _write_csv(staged["matrix"], rows)
        _write_csv(staged["bands"], band_rows)
        _write_csv(staged["nonadditivity"], nonadditivity_rows)
        _write_csv(staged["review"], review_rows)
        _write_csv(staged["plot_data"], plot_rows)
        render_v3_double_mutant_results(plot_rows, staged["png"], staged["svg"])
        elapsed_seconds = round(time.perf_counter() - started, 6)
        inputs = {
            "plan_manifest": {
                "path": _relative(plan_manifest_path),
                "sha256": _sha256(plan_manifest_path),
            },
            "candidate_plan": {
                "path": _relative(candidate_path),
                "sha256": _sha256(candidate_path),
            },
            "netsolp_scores": {
                "path": _relative(net_score_path),
                "sha256": _sha256(net_score_path),
            },
            "netsolp_run": {
                "path": _relative(net_run_path),
                "sha256": _sha256(net_run_path),
            },
            "nanomelt_scores": {
                "path": _relative(melt_score_path),
                "sha256": _sha256(melt_score_path),
            },
            "nanomelt_run": {
                "path": _relative(melt_run_path),
                "sha256": _sha256(melt_run_path),
            },
        }
        output_keys = tuple(key for key in OUTPUT_NAMES if key != "manifest")
        gate = {
            "v3_double_complete_property_matrix": "pass",
            "release": "ready_for_expert_review_and_explicit_final15_double_selection",
            "candidate_count": 102,
            "netsolp_coverage": "103_of_103_pass",
            "nanomelt_coverage": "103_of_103_pass",
            "candidate_prefiltering_applied": False,
            "final_15_double_mutant_selection": "not_performed",
            "final_30_panel_release": "not_performed",
        }
        manifest = {
            "schema_version": 1,
            "status": "pass",
            "workflow": "v3_double_mutant_complete_property_matrix",
            "generated_at": generated_at,
            "optimization_target": "Nb252 BL21 expression yield",
            "scientific_scope": (
                "Join complete NetSolP and NanoMelt predictions for all 102 "
                "released doubles, retain constituent AntiFold evidence and "
                "pre-score structure/risk annotations, and perform no final selection."
            ),
            "runtime_contract": {
                "execution": "project_environment_after_sequential_remote_scoring",
                "expected_runtime": "under_five_minutes_after_tool_scores_exist",
                "checkpoint_or_resume": False,
                "elapsed_seconds": elapsed_seconds,
                "python": platform.python_version(),
                "command_argv": [
                    sys.executable,
                    str(Path(__file__).resolve()),
                    *sys.argv[1:],
                ],
            },
            "metric_contract": {
                "positive_metrics": [
                    "NetSolP U",
                    "NetSolP S",
                    "NanoMelt predicted apparent Tm",
                ],
                "magnitude_thresholds": {
                    "netsolp_u": [0.005, 0.010, 0.015],
                    "netsolp_s": [0.010, 0.020, 0.030],
                    "nanomelt_tm_c": [0.5, 1.0, 1.5],
                },
                "within_band_decimal_ranking_performed": False,
                "model_nonadditivity": (
                    "double delta minus the two parent-single deltas; diagnostic "
                    "of predictor output only, not physical epistasis"
                ),
                "antifold": (
                    "two constituent pass records only; no addition and no "
                    "double-mutant AntiFold score"
                ),
            },
            "facts": facts,
            "inputs": inputs,
            "outputs": {
                key: {"path": _relative(finals[key]), "sha256": _sha256(staged[key])}
                for key in output_keys
            },
            "run_summary": _relative(run_summary),
            "verification": {
                "complete_double_candidate_rows": len(rows),
                "plot_rows": len(plot_rows),
                "both_tools_use_identical_103_id_sequence_table": True,
                "all_nanomelt_rows_score_126aa_and_trim_terminal_gs": True,
                "all_antifold_evidence_is_constituent_only": True,
                "model_nonadditivity_not_labeled_physical_epistasis": True,
                "candidate_selection_performed": False,
            },
            "gate": gate,
        }
        _write_json(staged["manifest"], manifest)
        _write_json(
            staged_run,
            {
                "schema_version": 1,
                "status": "pass",
                "workflow": manifest["workflow"],
                "generated_at": generated_at,
                "elapsed_seconds": elapsed_seconds,
                "python": platform.python_version(),
                "command_argv": manifest["runtime_contract"]["command_argv"],
                "inputs": inputs,
                "facts": facts,
                "verification": manifest["verification"],
                "outputs": {
                    **{key: _relative(finals[key]) for key in OUTPUT_NAMES},
                    "run_summary": _relative(run_summary),
                },
                "gate": gate,
            },
        )
        replace_staged_files(
            {
                **{staged[key]: finals[key] for key in OUTPUT_NAMES},
                staged_run: run_summary,
            },
            project_root=ROOT,
            protected_source_paths=validated.source_paths,
        )
    print("Released complete V3 102-double property matrix; final15 selection not performed")
    return 0


def _annotate_review_class(source):
    row = dict(source)
    bands = [
        str(row["netsolp_u_magnitude_band"]),
        str(row["netsolp_s_magnitude_band"]),
        str(row["nanomelt_tm_c_magnitude_band"]),
    ]
    favorable = sum(value in {"moderate_favorable", "strong_favorable"} for value in bands)
    moderate_adverse = sum(value == "moderate_adverse" for value in bands)
    strong_adverse = sum(value == "strong_adverse" for value in bands)
    if strong_adverse:
        review_class = "strong_adverse_property_requires_review"
    elif favorable and not moderate_adverse:
        review_class = "property_supported_no_moderate_or_strong_adverse"
    elif favorable:
        review_class = "property_supported_with_moderate_adverse_tradeoff"
    else:
        review_class = "no_moderate_or_strong_positive_metric"
    detailed_reasons = [
        value
        for value in str(row["machine_structure_triage_triggers"]).split("|")
        if value
    ]
    if favorable:
        detailed_reasons.append("property_supported_candidate_requires_expert_review")
    row.update(
        {
            "moderate_or_strong_favorable_metric_count": favorable,
            "moderate_adverse_metric_count": moderate_adverse,
            "strong_adverse_metric_count": strong_adverse,
            "property_review_class": review_class,
            "detailed_expert_review_required": bool(detailed_reasons),
            "detailed_expert_review_reasons": "|".join(dict.fromkeys(detailed_reasons)),
            "candidate_selection_performed": False,
        }
    )
    return row


def _band_counts(rows):
    output = []
    for metric, field in (
        ("NetSolP U", "netsolp_u_magnitude_band"),
        ("NetSolP S", "netsolp_s_magnitude_band"),
        ("NanoMelt predicted Tm", "nanomelt_tm_c_magnitude_band"),
    ):
        for band, count in sorted(Counter(str(row[field]) for row in rows).items()):
            output.append({"metric": metric, "magnitude_band": band, "candidate_count": count})
    return output


def _nonadditivity_summary(rows):
    output = []
    for metric, field in (
        ("NetSolP U", "netsolp_u_model_nonadditivity_residual"),
        ("NetSolP S", "netsolp_s_model_nonadditivity_residual"),
        ("NanoMelt predicted Tm", "nanomelt_tm_c_model_nonadditivity_residual"),
    ):
        values = np.asarray([float(row[field]) for row in rows], dtype=float)
        output.append(
            {
                "metric": metric,
                "candidate_count": len(values),
                "minimum": float(values.min()),
                "q25": float(np.quantile(values, 0.25)),
                "median": float(np.median(values)),
                "mean": float(values.mean()),
                "q75": float(np.quantile(values, 0.75)),
                "maximum": float(values.max()),
                "interpretation": "predictor_output_residual_not_physical_epistasis",
            }
        )
    return output


def _review_rows(rows):
    fields = (
        "v3_double_plan_order_not_efficacy_rank",
        "double_candidate_id",
        "mutation_set",
        "property_review_class",
        "moderate_or_strong_favorable_metric_count",
        "moderate_adverse_metric_count",
        "strong_adverse_metric_count",
        "machine_structure_triage_status",
        "machine_structure_triage_triggers",
        "detailed_expert_review_required",
        "detailed_expert_review_reasons",
        "hard_sequence_risk_flags",
        "soft_sequence_risk_flags",
        "stable_word_effect",
        "contains_t99f_stable_word_exploration_parent",
        "final_double_selection_status",
    )
    return [{field: row[field] for field in fields} for row in rows]


def _facts(rows):
    classes = Counter(str(row["property_review_class"]) for row in rows)
    return {
        "candidate_count": len(rows),
        "property_review_class_counts": dict(sorted(classes.items())),
        "moderate_or_strong_positive_candidate_count": sum(
            int(row["moderate_or_strong_favorable_metric_count"]) > 0 for row in rows
        ),
        "strong_adverse_candidate_count": sum(
            int(row["strong_adverse_metric_count"]) > 0 for row in rows
        ),
        "detailed_expert_review_required_count": sum(
            str(row["detailed_expert_review_required"]).lower() == "true"
            or row["detailed_expert_review_required"] is True
            for row in rows
        ),
        "stable_word_gain_candidate_count": sum(
            int(row["net_stable_word_occurrence_delta"]) > 0 for row in rows
        ),
        "hard_sequence_risk_candidate_count": sum(
            int(row["hard_sequence_risk_count"]) > 0 for row in rows
        ),
        "candidate_selection_performed": False,
    }


def _verify_tool_runs(plan_manifest, net_run, melt_run):
    gate = plan_manifest.get("gate", {})
    if not (
        plan_manifest.get("status") == "pass"
        and gate.get("release") == "ready_for_complete_netsolp_nanomelt_scoring"
        and net_run.get("status") == "pass"
        and net_run.get("tool") == "netsolp"
        and int(net_run.get("pass_count", -1)) == EXPECTED_SCORE_SAMPLE_COUNT
        and melt_run.get("status") == "pass"
        and melt_run.get("tool") == "nanomelt"
        and int(melt_run.get("pass_count", -1)) == EXPECTED_SCORE_SAMPLE_COUNT
    ):
        raise ValueError("V3 plan or one of the 103/103 tool runs did not pass")


def _csv(path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _fields(rows):
    fields = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    return fields


def _write_csv(path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=_fields(rows), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _write_json(path, value):
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def _sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _relative(path):
    resolved = path.resolve()
    try:
        return resolved.relative_to(ROOT).as_posix()
    except ValueError:
        return str(resolved)


if __name__ == "__main__":
    raise SystemExit(main())
