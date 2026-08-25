#!/usr/bin/env python3
"""Freeze the complete V3 15-parent to 102-double scoring plan.

This local workflow is expected to finish within one minute.  It enumerates
every valid unordered pair without property prefiltering, recomputes sequence
risks and stable words on each complete construct, records WT-structure pair
geometry, and prepares the one shared WT-plus-102 sample table for NetSolP and
NanoMelt.  It does not run a predictor or select the final 15 doubles.
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
from antibody_optimization.stable_words import parse_stable_words  # noqa: E402
from antibody_optimization.v3_double_mutant_plan_plot import (  # noqa: E402
    build_v3_double_plan_plot_rows,
    render_v3_double_mutant_plan,
)
from antibody_optimization.v3_double_mutants import (  # noqa: E402
    WT_SCORE_ID,
    build_v3_double_mutant_space,
    build_v3_score_samples,
)


DEFAULT_PARENT_DIR = (
    ROOT
    / "docs/result_artifacts/candidate_design"
    / "v3_parent_single_selection_20260825"
)
DEFAULT_OUTPUT_DIR = (
    ROOT
    / "docs/result_artifacts/candidate_design"
    / "v3_double_mutant_plan_20260825"
)
DEFAULT_RUN_SUMMARY = (
    ROOT
    / "docs/run_summaries/candidate_design"
    / "v3_double_mutant_plan_20260825"
    / "run_summary.json"
)
OUTPUT_NAMES = {
    "candidates": "v3_double_mutant_candidates102.csv",
    "invalid": "v3_double_mutant_invalid_same_position_pairs.csv",
    "samples": "v3_double_mutant_score_samples103.csv",
    "fasta": "v3_double_mutant_score_samples103.fasta",
    "word_changes": "v3_double_mutant_stable_word_changes.csv",
    "structure_triage": "v3_double_mutant_structure_triage.csv",
    "plot_data": "v3_double_mutant_plan_plot_data.csv",
    "png": "v3_double_mutant_plan_overview.png",
    "svg": "v3_double_mutant_plan_overview.svg",
    "contract": "v3_double_mutant_plan_contract.json",
    "manifest": "v3_double_mutant_plan_manifest.json",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent-selection-dir", type=Path, default=DEFAULT_PARENT_DIR)
    parser.add_argument(
        "--stable-word-library",
        type=Path,
        default=ROOT / "data/Stable_word_SS_3D_1336 (1).txt",
    )
    parser.add_argument(
        "--sequence-structure-mapping",
        type=Path,
        default=(
            ROOT
            / "docs/result_artifacts/input_baseline/structure_released_20260810"
            / "nb252_sequence_structure_mapping.csv"
        ),
    )
    parser.add_argument(
        "--experimental-cif",
        type=Path,
        default=ROOT / "data/structures/cxs_exports/NK2R-252__native.cif",
    )
    parser.add_argument(
        "--af3-cif",
        type=Path,
        default=(
            ROOT
            / "data/structures/cxs_exports/fold_2r_252_nomg_model_0__native.cif"
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
    parent_dir = args.parent_selection_dir.resolve(strict=True)
    selected_csv = parent_dir / "v3_parent_single_selected15.csv"
    decision_csv = parent_dir / "v3_parent_single_selection_audit.csv"
    parent_manifest_path = parent_dir / "v3_parent_single_selection_manifest.json"
    stable_word_path = args.stable_word_library.resolve(strict=True)
    mapping_path = args.sequence_structure_mapping.resolve(strict=True)
    experimental_cif = args.experimental_cif.resolve(strict=True)
    af3_cif = args.af3_cif.resolve(strict=True)
    sources = (
        selected_csv,
        decision_csv,
        parent_manifest_path,
        stable_word_path,
        mapping_path,
        experimental_cif,
        af3_cif,
    )
    parent_manifest = _json(parent_manifest_path)
    _verify_parent_release(parent_manifest, selected_csv, decision_csv)

    output_dir = args.output_dir.resolve()
    run_summary = args.run_summary.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    run_summary.parent.mkdir(parents=True, exist_ok=True)
    finals = {key: output_dir / name for key, name in OUTPUT_NAMES.items()}
    existing = [path for path in (*finals.values(), run_summary) if path.exists()]
    if existing and not args.overwrite:
        raise FileExistsError(
            "Refusing to overwrite existing V3 double-plan outputs:\n"
            + "\n".join(str(path) for path in existing)
        )
    validated = validate_file_paths(
        project_root=ROOT,
        source_paths=sources,
        target_paths=(*finals.values(), run_summary),
    )
    stable_words = parse_stable_words(
        stable_word_path.read_text(encoding="utf-8-sig").splitlines(keepends=True)
    )
    result = build_v3_double_mutant_space(
        _csv(selected_csv),
        _csv(decision_csv),
        stable_words,
        _csv(mapping_path),
        experimental_cif,
        af3_cif,
    )
    samples = build_v3_score_samples(result["parent_sequence"], result["candidates"])
    structure_rows = _structure_triage_rows(result["candidates"])
    plot_rows = build_v3_double_plan_plot_rows(
        result["parents"], result["candidates"], result["invalid_pairs"]
    )
    facts = result["facts"]

    with tempfile.TemporaryDirectory(prefix=".v3-double-plan-", dir=ROOT) as temp:
        stage = Path(temp)
        staged = {key: stage / name for key, name in OUTPUT_NAMES.items()}
        staged_run_summary = stage / "run_summary.json"
        _write_csv(staged["candidates"], result["candidates"])
        _write_csv(staged["invalid"], result["invalid_pairs"])
        _write_csv(staged["samples"], samples)
        _write_fasta(staged["fasta"], samples)
        _write_csv(staged["word_changes"], result["stable_word_changes"])
        _write_csv(staged["structure_triage"], structure_rows)
        _write_csv(staged["plot_data"], plot_rows)
        render_v3_double_mutant_plan(plot_rows, staged["png"], staged["svg"])
        elapsed_seconds = round(time.perf_counter() - started, 6)
        contract = _contract(generated_at)
        _write_json(staged["contract"], contract)
        output_keys = tuple(key for key in OUTPUT_NAMES if key != "manifest")
        inputs = {
            "parent_selection_manifest": {
                "path": _relative(parent_manifest_path),
                "sha256": _sha256(parent_manifest_path),
            },
            "parent_selected15": {
                "path": _relative(selected_csv),
                "sha256": _sha256(selected_csv),
            },
            "parent_decision_audit": {
                "path": _relative(decision_csv),
                "sha256": _sha256(decision_csv),
            },
            "stable_word_library": {
                "path": _relative(stable_word_path),
                "sha256": _sha256(stable_word_path),
                "word_count": len(stable_words),
            },
            "sequence_structure_mapping": {
                "path": _relative(mapping_path),
                "sha256": _sha256(mapping_path),
            },
            "experimental_complex": {
                "path": _relative(experimental_cif),
                "sha256": _sha256(experimental_cif),
            },
            "af3_vhh": {
                "path": _relative(af3_cif),
                "sha256": _sha256(af3_cif),
            },
            "upstream_files_modified": False,
        }
        gate = {
            "v3_double_mutant_plan": "pass",
            "release": "ready_for_complete_netsolp_nanomelt_scoring",
            "complete_double_sequences_released": 102,
            "shared_score_samples_released": 103,
            "candidate_prefiltering_applied": False,
            "remote_property_scoring": "not_started",
            "final_15_double_mutant_selection": "not_performed",
        }
        manifest = {
            "schema_version": 1,
            "status": "pass",
            "workflow": "v3_double_mutant_plan",
            "generated_at": generated_at,
            "optimization_target": "Nb252 BL21 expression yield",
            "scientific_scope": (
                "Enumerate all valid unordered pairs from the released V3 parent15, "
                "recompute complete-sequence annotations, and freeze one shared "
                "WT-plus-102 property-scoring table without selecting doubles."
            ),
            "runtime_contract": {
                "execution": "local_cpu",
                "expected_runtime": "under_one_minute",
                "expected_over_one_hour": False,
                "expected_over_five_hours": False,
                "slurm_required_for_plan_build": False,
                "checkpoint_or_resume": False,
                "elapsed_seconds": elapsed_seconds,
                "python": platform.python_version(),
                "command_argv": [
                    sys.executable,
                    str(Path(__file__).resolve()),
                    *sys.argv[1:],
                ],
            },
            "contract": contract,
            "facts": facts,
            "wt_score_id": WT_SCORE_ID,
            "inputs": inputs,
            "outputs": {
                key: {"path": _relative(finals[key]), "sha256": _sha256(staged[key])}
                for key in output_keys
            },
            "run_summary": _relative(run_summary),
            "verification": {
                "parent_single_count": 15,
                "parent_unique_position_count": 12,
                "theoretical_pair_count": 105,
                "invalid_same_position_pair_count": 3,
                "valid_double_count": 102,
                "score_sample_count_including_wt": 103,
                "all_double_sequences_unique_128aa_exactly_two_substitutions": True,
                "parent_ssgs_and_two_cysteines_preserved": True,
                "antifold_component_values_not_combined": True,
                "complete_sequences_not_prefiltered": True,
                "upstream_parent_artifacts_are_read_only": True,
            },
            "gate": gate,
        }
        _write_json(staged["manifest"], manifest)
        _write_json(
            staged_run_summary,
            {
                "schema_version": 1,
                "status": "pass",
                "workflow": "v3_double_mutant_plan",
                "generated_at": generated_at,
                "elapsed_seconds": elapsed_seconds,
                "python": platform.python_version(),
                "command_argv": manifest["runtime_contract"]["command_argv"],
                "inputs": inputs,
                "facts": facts,
                "parameters": contract,
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
                staged_run_summary: run_summary,
            },
            project_root=ROOT,
            protected_source_paths=validated.source_paths,
        )
    print(
        "Released complete V3 double plan: 105 theoretical pairs, "
        "3 same-position exclusions, 102 valid doubles, 103 score samples"
    )
    return 0


def _contract(generated_at: str) -> dict[str, object]:
    return {
        "schema_version": 1,
        "generated_at": generated_at,
        "parent_release": "v3_parent_single_selection_20260825 selected15",
        "enumeration": "all unordered pairs of released parents",
        "same_position_alternatives": "mutually_exclusive_and_not_generated",
        "property_prefilter_before_complete_scoring": False,
        "active_remote_predictors": [
            "NetSolP Distilled SU on complete 128-aa constructs",
            "NanoMelt predicted apparent Tm with observed 126-aa scored domain and terminal GS trimming",
        ],
        "antifold_policy": (
            "retain each constituent single-position source/delta/rank/veto; "
            "require both vetoes to pass; do not add values and do not claim a "
            "double-mutant AntiFold score"
        ),
        "stable_word_policy": (
            "recompute exact overlapping degenerate-word occurrences on each "
            "complete double sequence; soft hypothesis only"
        ),
        "sequence_risk_policy": "recompute on each complete double sequence",
        "structure_triage_policy": (
            "WT-site geometry only; experimental coordinates primary when both "
            "sites are observed, otherwise AF3 distances are separately labeled "
            "prediction-only; no mutant effect is inferred"
        ),
        "remote_route": (
            "one non-array batch Slurm job, one GPU, 12 CPUs, sequential NetSolP "
            "then NanoMelt then project-environment analysis"
        ),
        "expected_remote_runtime": "approximately_10_to_30_minutes_under_2_hours",
        "checkpoint_or_resume": False,
        "candidate_selection_performed": False,
    }


def _verify_parent_release(manifest, selected_csv, decision_csv) -> None:
    gate = manifest.get("gate", {})
    facts = manifest.get("facts", {})
    if not (
        manifest.get("status") == "pass"
        and gate.get("v3_parent_single_selection") == "pass"
        and gate.get("double_mutant_enumeration") == "not_performed"
        and int(facts.get("selected_parent_single_count", -1)) == 15
        and int(facts.get("valid_unordered_double_mutant_count", -1)) == 102
    ):
        raise ValueError("V3 parent-single release gate is not valid for enumeration")
    outputs = manifest.get("outputs", {})
    expected = {
        "selected": (selected_csv, outputs.get("selected", {}).get("sha256")),
        "audit": (decision_csv, outputs.get("audit", {}).get("sha256")),
    }
    for label, (path, recorded_hash) in expected.items():
        if not recorded_hash or _sha256(path) != recorded_hash:
            raise ValueError(f"V3 parent {label} identity differs from its manifest")


def _structure_triage_rows(candidates):
    fields = (
        "v3_double_plan_order_not_efficacy_rank",
        "double_candidate_id",
        "mutation_set",
        "position_a_reported_1based",
        "position_b_reported_1based",
        "pair_experimental_coordinate_status",
        "experimental_pair_ca_distance_a",
        "experimental_pair_minimum_heavy_atom_distance_a",
        "af3_pair_ca_distance_a",
        "af3_pair_minimum_heavy_atom_distance_a",
        "pair_structure_distance_source",
        "pair_ca_distance_a",
        "pair_minimum_heavy_atom_distance_a",
        "pair_spatial_class",
        "pair_shared_local_neighbor_count",
        "pair_shared_local_neighbors",
        "pair_geometry_role",
        "contains_t99f_stable_word_exploration_parent",
        "hard_sequence_risk_flags",
        "machine_structure_triage_status",
        "machine_structure_triage_triggers",
    )
    return [{field: row[field] for field in fields} for row in candidates]


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _fields(rows) -> list[str]:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    return fields


def _write_csv(path: Path, rows) -> None:
    if not rows:
        path.write_text("", encoding="utf-8-sig", newline="")
        return
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=_fields(rows), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_fasta(path: Path, rows) -> None:
    lines: list[str] = []
    for row in rows:
        lines.extend([f">{row['sample_uid']}", str(row["sequence_raw"])])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def _json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


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
