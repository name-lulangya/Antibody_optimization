#!/usr/bin/env python3
"""Release the approved 15 V3 parent single mutants and 31-row audit.

This local workflow is expected to finish within one minute.  It reuses the
released predictor and expert-review evidence, writes no structures, does not
run any predictor, and does not enumerate or select double mutants.
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
from antibody_optimization.v3_parent_single_selection import (  # noqa: E402
    HARD_EXPERT_EXCLUSION_IDS,
    SELECTED_PARENT_IDS,
    V3_T99F_ID,
    build_v3_parent_single_selection,
    selected_parent_export_rows,
)
from antibody_optimization.v3_parent_single_selection_plot import (  # noqa: E402
    build_v3_parent_selection_plot_rows,
    render_v3_parent_single_selection,
)


DEFAULT_OUTPUT_DIR = (
    ROOT
    / "docs/result_artifacts/candidate_design"
    / "v3_parent_single_selection_20260825"
)
DEFAULT_RUN_SUMMARY = (
    ROOT
    / "docs/run_summaries/candidate_design"
    / "v3_parent_single_selection_20260825"
    / "run_summary.json"
)
OUTPUT_NAMES = {
    "audit": "v3_parent_single_selection_audit.csv",
    "selected": "v3_parent_single_selected15.csv",
    "fasta": "v3_parent_single_selected15.fasta",
    "plot_data": "v3_parent_single_selection_plot_data.csv",
    "png": "v3_parent_single_selection_overview.png",
    "svg": "v3_parent_single_selection_overview.svg",
    "manifest": "v3_parent_single_selection_manifest.json",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--expert-review-csv",
        type=Path,
        default=(
            ROOT
            / "docs/result_artifacts/candidate_design"
            / "v3_parent_single_expert_review_20260825"
            / "v3_parent_single_expert_review.csv"
        ),
    )
    parser.add_argument(
        "--complete-v3-audit-csv",
        type=Path,
        default=(
            ROOT
            / "docs/result_artifacts/candidate_design"
            / "expression_single_mutant_selection_v3_20260825"
            / "expression_single_mutant_v3_audit.csv"
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
    expert_review = args.expert_review_csv.resolve(strict=True)
    complete_audit = args.complete_v3_audit_csv.resolve(strict=True)
    output_dir = args.output_dir.resolve()
    run_summary = args.run_summary.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    run_summary.parent.mkdir(parents=True, exist_ok=True)
    finals = {key: output_dir / name for key, name in OUTPUT_NAMES.items()}
    existing = [path for path in (*finals.values(), run_summary) if path.exists()]
    if existing and not args.overwrite:
        raise FileExistsError(
            "Refusing to overwrite existing V3 parent-selection outputs:\n"
            + "\n".join(str(path) for path in existing)
        )
    validated = validate_file_paths(
        project_root=ROOT,
        source_paths=(expert_review, complete_audit),
        target_paths=(*finals.values(), run_summary),
    )

    result = build_v3_parent_single_selection(
        _csv(expert_review),
        _csv(complete_audit),
    )
    selected_export = selected_parent_export_rows(result["selected_rows"])
    plot_rows = build_v3_parent_selection_plot_rows(result["audit_rows"])
    facts = result["facts"]

    with tempfile.TemporaryDirectory(prefix=".v3-parent-selection-", dir=ROOT) as temp:
        stage = Path(temp)
        staged = {key: stage / name for key, name in OUTPUT_NAMES.items()}
        staged_run_summary = stage / "run_summary.json"
        _write_csv(staged["audit"], result["audit_rows"])
        _write_csv(staged["selected"], selected_export)
        _write_fasta(staged["fasta"], selected_export)
        _write_csv(staged["plot_data"], plot_rows)
        render_v3_parent_single_selection(plot_rows, staged["png"], staged["svg"])
        elapsed_seconds = round(time.perf_counter() - started, 6)
        artifact_output_keys = ("audit", "selected", "fasta", "plot_data", "png", "svg")
        manifest = {
            "schema_version": 1,
            "status": "pass",
            "generated_at": generated_at,
            "workflow": "v3_parent_single_selection",
            "optimization_target": "Nb252 BL21 expression yield",
            "scientific_scope": (
                "Freeze the user-approved 15 parent single mutants from the released "
                "31-candidate expert-review pool and preserve a detailed decision for "
                "all selected and non-selected candidates."
            ),
            "runtime_contract": {
                "execution": "local_cpu",
                "expected_runtime": "under_one_minute",
                "slurm_required": False,
                "checkpoint_or_resume": False,
                "elapsed_seconds": elapsed_seconds,
                "python": platform.python_version(),
                "command_argv": [
                    sys.executable,
                    str(Path(__file__).resolve()),
                    *sys.argv[1:],
                ],
            },
            "selection_policy": {
                "positive_metrics": [
                    "NetSolP U magnitude band",
                    "NetSolP S magnitude band",
                    "NanoMelt predicted Tm magnitude band",
                ],
                "within_band_raw_decimals_used_as_global_rank": False,
                "strong_improvement_priority": (
                    "Select every reviewed candidate with at least one strong-favorable "
                    "property band unless it has a high-confidence concrete expert risk."
                ),
                "expert_hard_exclusion_rule": (
                    "Hard exclusion requires structurally_concerning plus high expert "
                    "confidence and a concrete physical risk documented in the review."
                ),
                "antifold_role": "negative veto only; all 31 reviewed candidates passed",
                "stable_word_role": (
                    "soft hypothesis only; T99F is one user-directed exploratory exception"
                ),
                "same_position_rule": (
                    "Alternatives may both be parent singles but their mutual pair is invalid"
                ),
            },
            "selected_parent_ids_in_display_order_not_efficacy_rank": list(
                SELECTED_PARENT_IDS
            ),
            "selected_mutations_in_display_order_not_efficacy_rank": [
                row["mutation_reported_label"].replace("Nb252 reported_seq ", "")
                for row in selected_export
            ],
            "high_confidence_expert_risk_exclusion_ids": sorted(
                HARD_EXPERT_EXCLUSION_IDS
            ),
            "user_directed_exploration_id": V3_T99F_ID,
            "facts": facts,
            "inputs": {
                "expert_review": {
                    "path": _relative(expert_review),
                    "sha256": _sha256(expert_review),
                },
                "complete_v3_audit": {
                    "path": _relative(complete_audit),
                    "sha256": _sha256(complete_audit),
                },
                "upstream_files_modified": False,
            },
            "outputs": {
                key: {
                    "path": _relative(finals[key]),
                    "sha256": _sha256(staged[key]),
                }
                for key in artifact_output_keys
            },
            "run_summary": _relative(run_summary),
            "verification": {
                "decision_rows": len(result["audit_rows"]),
                "selected_rows": len(selected_export),
                "plot_source_rows": len(plot_rows),
                "selected_sequences_are_unique_128aa_single_mutants": True,
                "parent_ssgs_and_two_cysteines_preserved": True,
                "upstream_expert_review_and_v3_audit_are_read_only": True,
            },
            "gate": {
                "v3_parent_single_selection": "pass",
                "selected_parent_count": 15,
                "detailed_decision_audit_count": 31,
                "valid_double_mutant_space_ready_for_generation": True,
                "double_mutant_enumeration": "not_performed",
                "final_15_double_mutant_selection": "not_performed",
            },
        }
        _write_json(staged["manifest"], manifest)
        _write_json(
            staged_run_summary,
            {
                "schema_version": 1,
                "status": "pass",
                "workflow": "v3_parent_single_selection",
                "generated_at": generated_at,
                "elapsed_seconds": elapsed_seconds,
                "python": platform.python_version(),
                "command_argv": manifest["runtime_contract"]["command_argv"],
                "inputs": manifest["inputs"],
                "facts": facts,
                "selection_policy": manifest["selection_policy"],
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
                staged_run_summary: run_summary,
            },
            project_root=ROOT,
            protected_source_paths=validated.source_paths,
        )
    print(
        "Released 15 V3 parent singles and 31 detailed decisions; "
        f"future valid unordered double-mutant space={facts['valid_unordered_double_mutant_count']}"
    )
    return 0


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _fields(rows: list[dict[str, object]]) -> list[str]:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    return fields


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=_fields(rows), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_fasta(path: Path, rows: list[dict[str, object]]) -> None:
    lines: list[str] = []
    for row in rows:
        mutation = str(row["mutation_reported_label"]).replace(
            "Nb252 reported_seq ", ""
        )
        lines.extend(
            [
                (
                    f">{row['candidate_id']} mutation={mutation} "
                    f"parent_panel_order_not_efficacy_rank="
                    f"{row['v3_parent_panel_order_not_efficacy_rank']}"
                ),
                str(row["sequence"]),
            ]
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def _write_json(path: Path, value: object) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _relative(path: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(ROOT).as_posix()
    except ValueError:
        return str(resolved)


if __name__ == "__main__":
    raise SystemExit(main())
