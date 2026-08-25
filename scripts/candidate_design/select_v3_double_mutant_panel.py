#!/usr/bin/env python3
"""Release the reviewed V3 final 15 doubles and 15+15 experimental panel.

This local workflow is expected to finish within one minute.  It reuses the
released 102-row property matrix and parent expert evidence; it does not rerun
predictors, model double-mutant structures, require Slurm, or create restart
state.
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
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from antibody_optimization.file_transaction import (  # noqa: E402
    replace_staged_files,
    validate_file_paths,
)
from antibody_optimization.v3_double_mutant_selection import (  # noqa: E402
    SELECTED_DOUBLE_MUTATION_SETS,
    build_v3_double_mutant_selection,
    selected_double_export_rows,
)
from antibody_optimization.v3_double_mutant_selection_plot import (  # noqa: E402
    build_v3_final_panel_plot_rows,
    render_v3_final_panel_selection,
)


DEFAULT_OUTPUT_DIR = (
    ROOT
    / "docs/result_artifacts/candidate_design"
    / "v3_final_15plus15_panel_20260825"
)
DEFAULT_RUN_SUMMARY = (
    ROOT
    / "docs/run_summaries/candidate_design"
    / "v3_final_15plus15_panel_20260825"
    / "run_summary.json"
)
OUTPUT_NAMES = {
    "audit": "v3_double_mutant_final_selection_audit102.csv",
    "selected": "v3_double_mutant_selected15.csv",
    "final_panel": "v3_final_panel30.csv",
    "fasta": "v3_final_panel30.fasta",
    "plot_data": "v3_final_panel_plot_data.csv",
    "png": "v3_final_panel_overview.png",
    "svg": "v3_final_panel_overview.svg",
    "manifest": "v3_final_panel_manifest.json",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--matrix-dir",
        type=Path,
        default=(
            ROOT
            / "docs/result_artifacts/candidate_design"
            / "v3_double_mutant_property_matrix_20260825"
        ),
    )
    parser.add_argument(
        "--parent-selection-dir",
        type=Path,
        default=(
            ROOT
            / "docs/result_artifacts/candidate_design"
            / "v3_parent_single_selection_20260825"
        ),
    )
    parser.add_argument(
        "--post-sync-review",
        type=Path,
        default=(
            ROOT
            / "docs/result_artifacts/candidate_design"
            / "v3_double_mutant_post_sync_review_20260825"
            / "v3_double_mutant_post_sync_review.json"
        ),
    )
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
    matrix_dir = args.matrix_dir.resolve(strict=True)
    parent_dir = args.parent_selection_dir.resolve(strict=True)
    post_sync_review = args.post_sync_review.resolve(strict=True)
    matrix_csv = (matrix_dir / "v3_double_mutant_property_matrix102.csv").resolve(
        strict=True
    )
    matrix_manifest_path = (
        matrix_dir / "v3_double_mutant_property_matrix_manifest.json"
    ).resolve(strict=True)
    parent_selected_csv = (
        parent_dir / "v3_parent_single_selected15.csv"
    ).resolve(strict=True)
    parent_audit_csv = (
        parent_dir / "v3_parent_single_selection_audit.csv"
    ).resolve(strict=True)
    parent_manifest_path = (
        parent_dir / "v3_parent_single_selection_manifest.json"
    ).resolve(strict=True)
    output_dir = args.output_dir.resolve()
    run_summary = args.run_summary.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    run_summary.parent.mkdir(parents=True, exist_ok=True)
    finals = {key: output_dir / name for key, name in OUTPUT_NAMES.items()}
    existing = [path for path in (*finals.values(), run_summary) if path.exists()]
    if existing and not args.overwrite:
        raise FileExistsError(
            "Refusing to overwrite existing V3 final-panel outputs:\n"
            + "\n".join(str(path) for path in existing)
        )
    validated = validate_file_paths(
        project_root=ROOT,
        source_paths=(
            matrix_csv,
            matrix_manifest_path,
            parent_selected_csv,
            parent_audit_csv,
            parent_manifest_path,
            post_sync_review,
        ),
        target_paths=(*finals.values(), run_summary),
    )

    matrix_manifest = _json(matrix_manifest_path)
    parent_manifest = _json(parent_manifest_path)
    post_sync = _json(post_sync_review)
    _validate_upstream_gates(
        matrix_manifest,
        parent_manifest,
        matrix_csv=matrix_csv,
        parent_selected_csv=parent_selected_csv,
        parent_audit_csv=parent_audit_csv,
    )
    result = build_v3_double_mutant_selection(
        _csv(matrix_csv),
        _csv(parent_selected_csv),
        _csv(parent_audit_csv),
        post_sync,
    )
    selected_export = selected_double_export_rows(result["selected_double_rows"])
    plot_rows = build_v3_final_panel_plot_rows(result["audit_rows"])
    parent_mutations = [
        str(row["mutation_reported_label"]).replace("Nb252 reported_seq ", "")
        for row in _csv(parent_selected_csv)
    ]

    with tempfile.TemporaryDirectory(prefix=".v3-final-panel-", dir=ROOT) as temp:
        stage = Path(temp)
        staged = {key: stage / name for key, name in OUTPUT_NAMES.items()}
        staged_summary = stage / "run_summary.json"
        _write_csv(staged["audit"], result["audit_rows"])
        _write_csv(staged["selected"], selected_export)
        _write_csv(staged["final_panel"], result["final_panel_rows"])
        _write_fasta(staged["fasta"], result["final_panel_rows"])
        _write_csv(staged["plot_data"], plot_rows)
        render_v3_final_panel_selection(
            plot_rows,
            parent_mutations,
            staged["png"],
            staged["svg"],
        )
        elapsed_seconds = round(time.perf_counter() - started, 6)
        output_keys = ("audit", "selected", "final_panel", "fasta", "plot_data", "png", "svg")
        command_argv = [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]]
        inputs = {
            "double_property_matrix": _record(matrix_csv),
            "double_property_matrix_manifest": _record(matrix_manifest_path),
            "parent_selected15": _record(parent_selected_csv),
            "parent_selection_audit31": _record(parent_audit_csv),
            "parent_selection_manifest": _record(parent_manifest_path),
            "post_sync_annotation_review": _record(post_sync_review),
        }
        manifest = {
            "schema_version": 1,
            "status": "pass",
            "generated_at": generated_at,
            "workflow": "v3_double_mutant_final_selection",
            "optimization_target": "Nb252 BL21 expression yield",
            "scientific_scope": (
                "Complete a common expert review of all 102 valid V3 doubles, use "
                "enhanced versus standard only as review depth, explicitly select 15 "
                "doubles, and release the final 15-single plus 15-double panel."
            ),
            "runtime_contract": {
                "execution": "local_cpu",
                "expected_runtime": "under_one_minute",
                "exceeds_one_hour": False,
                "exceeds_five_hours": False,
                "slurm_required": False,
                "checkpoint_or_resume": False,
                "elapsed_seconds": elapsed_seconds,
                "python": platform.python_version(),
                "command_argv": command_argv,
            },
            "selection_policy": result["selection_policy"],
            "review_depth_contract": {
                "enhanced_definition": (
                    "at least two moderate/strong favorable metrics with no moderate/strong "
                    "adverse, or any moderate/strong adverse, or a non-separated WT site pair"
                ),
                "standard_definition": "all remaining candidates",
                "interpretation": "review detail only; not eligibility, rank, or elimination",
            },
            "selected_double_mutation_sets_in_display_order_not_efficacy_rank": list(
                SELECTED_DOUBLE_MUTATION_SETS
            ),
            "selected_double_ids_in_display_order_not_efficacy_rank": [
                row["double_candidate_id"] for row in selected_export
            ],
            "facts": result["facts"],
            "inputs": inputs,
            "outputs": {
                key: {"path": _relative(finals[key]), "sha256": _sha256(staged[key])}
                for key in output_keys
            },
            "run_summary": _relative(run_summary),
            "verification": {
                "all_102_doubles_share_common_review_fields": True,
                "review_depth_does_not_control_selection": True,
                "t99f_has_no_mutation_specific_rule": True,
                "post_sync_deamidation_erratum_applied_as_overlay": True,
                "source_matrix_and_parent_artifacts_modified": False,
                "selected_doubles_have_two_or_three_positive_bands": True,
                "selected_doubles_have_no_moderate_or_strong_adverse_band": True,
                "antifold_is_constituent_negative_veto_only": True,
                "double_sidechain_modeling_performed": False,
                "final_sequences_are_unique_128aa_and_preserve_terminal_ssgs": True,
            },
            "gate": {
                "v3_double_expert_review": "pass",
                "final_15_double_mutant_selection": "pass",
                "final_30_panel_release": "pass",
                "report_and_presentation_sync": "not_performed",
            },
        }
        _write_json(staged["manifest"], manifest)
        _write_json(
            staged_summary,
            {
                "schema_version": 1,
                "status": "pass",
                "workflow": manifest["workflow"],
                "generated_at": generated_at,
                "elapsed_seconds": elapsed_seconds,
                "python": platform.python_version(),
                "command_argv": command_argv,
                "inputs": inputs,
                "facts": result["facts"],
                "selection_policy": result["selection_policy"],
                "verification": manifest["verification"],
                "outputs": {
                    **{key: _relative(finals[key]) for key in OUTPUT_NAMES},
                    "run_summary": _relative(run_summary),
                },
                "gate": manifest["gate"],
            },
        )
        replace_staged_files(
            {
                **{staged[key]: finals[key] for key in OUTPUT_NAMES},
                staged_summary: run_summary,
            },
            project_root=ROOT,
            protected_source_paths=validated.source_paths,
        )
    print(
        "Released V3 final panel: "
        f"{result['facts']['source_double_candidate_count']} doubles reviewed, "
        f"{result['facts']['selected_double_mutant_count']} doubles selected, "
        f"{result['facts']['final_panel_candidate_count']} total candidates"
    )
    return 0


def _validate_upstream_gates(
    matrix_manifest: dict[str, object],
    parent_manifest: dict[str, object],
    *,
    matrix_csv: Path,
    parent_selected_csv: Path,
    parent_audit_csv: Path,
) -> None:
    if matrix_manifest.get("status") != "pass" or matrix_manifest.get("gate", {}).get(
        "v3_double_complete_property_matrix"
    ) != "pass":
        raise ValueError("V3 double-property matrix gate has not passed")
    if parent_manifest.get("status") != "pass" or parent_manifest.get("gate", {}).get(
        "v3_parent_single_selection"
    ) != "pass":
        raise ValueError("V3 parent-single selection gate has not passed")
    expected = (
        (matrix_manifest["outputs"]["matrix"]["sha256"], matrix_csv),
        (parent_manifest["outputs"]["selected"]["sha256"], parent_selected_csv),
        (parent_manifest["outputs"]["audit"]["sha256"], parent_audit_csv),
    )
    for expected_sha, path in expected:
        if _sha256(path) != expected_sha:
            raise ValueError(f"Released input identity mismatch: {path}")


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise TypeError(f"Expected JSON object: {path}")
    return value


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Cannot write empty CSV: {path}")
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_fasta(path: Path, rows: list[dict[str, object]]) -> None:
    lines: list[str] = []
    for row in rows:
        lines.append(
            f">{row['candidate_id']} kind={row['candidate_kind']} mutations={row['mutation_set']}"
        )
        sequence = str(row["sequence"])
        lines.extend(sequence[index : index + 80] for index in range(0, len(sequence), 80))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def _write_json(path: Path, value: object) -> None:
    path.write_text(
        json.dumps(value, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _record(path: Path) -> dict[str, str]:
    return {"path": _relative(path), "sha256": _sha256(path)}


def _relative(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT).as_posix()
    except ValueError:
        return str(path.resolve())


if __name__ == "__main__":
    raise SystemExit(main())
