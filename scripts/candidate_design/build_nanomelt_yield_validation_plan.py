#!/usr/bin/env python3
"""Build the fixed 47-sequence NanoMelt–BL21 yield validation plan."""

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
from antibody_optimization.nanomelt_yield import (  # noqa: E402
    BOOTSTRAP_SEED,
    RESAMPLING_REPLICATES,
    build_nanomelt_validation_inputs,
)


NAMES = {
    "samples": "nanomelt_validation_samples.csv",
    "fasta": "nanomelt_validation_sequences.fasta",
    "contract": "nanomelt_yield_validation_contract.json",
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expression-records", type=Path, required=True)
    parser.add_argument("--numbering-review", type=Path, required=True)
    parser.add_argument("--numbering-positions", type=Path, required=True)
    parser.add_argument("--allowed-use-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    parser.add_argument("--check_only", action="store_true")
    args = parser.parse_args()
    started = time.perf_counter()
    generated = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    sources = [
        args.expression_records.resolve(strict=True),
        args.numbering_review.resolve(strict=True),
        args.numbering_positions.resolve(strict=True),
        args.allowed_use_manifest.resolve(strict=True),
    ]
    allowed = _json(sources[3])
    if allowed["gates"]["cross_assay_pooling_gate"] != "pass" or allowed["counts"]["samples"] != 47:
        raise ValueError("Expression audit does not release the fixed 47-sample exploratory analysis")
    samples = build_nanomelt_validation_inputs(_csv(sources[0]), _csv(sources[1]), _csv(sources[2]))["sample_rows"]
    if args.check_only:
        print(json.dumps({"status": "pass", "samples": len(samples), "numeric": 31, "llj": 16}, ensure_ascii=False))
        return 0

    targets = [args.output_dir.absolute() / name for name in NAMES.values()] + [args.run_summary.absolute()]
    valid = validate_file_paths(project_root=ROOT, source_paths=sources, target_paths=targets)
    if any(path.exists() for path in valid.target_paths):
        raise FileExistsError("Refusing to overwrite NanoMelt validation plan")
    for path in valid.target_paths:
        path.parent.mkdir(parents=True, exist_ok=True)
    final = dict(zip((*NAMES, "summary"), valid.target_paths, strict=True))
    contract = {
        "schema_version": 1,
        "status": "pass",
        "generated_at": generated,
        "expression_system": "BL21_E_coli",
        "remote_environment": "/data/software/env/luly25/nanomelt",
        "remote_entry_point": "/data/software/env/luly25/nanomelt/bin/nanomelt",
        "software": {
            "nanomelt": "1.4.0",
            "python": "3.10.20",
            "torch": "2.7.1+cu126",
            "transformers": "4.56.1",
            "immune_builder": "1.2",
            "anarci_bioconda": "2024.05.21",
            "openmm": "8.5.2",
            "pdbfixer": "1.12.0",
        },
        "required_openmm_platforms": ["Reference", "CPU"],
        "anarci_runtime_contract": {
            "required_api": ["anarci", "run_anarci", "validate_sequence", "scheme_short_to_long"],
            "required_hmm_database": "dat/HMMs/ALL.hmm_with_pressed_indexes",
            "distribution_metadata_is_not_a_version_gate": True,
        },
        "immune_builder_local_patch": {
            "patch_id": "openmm_threads_mapping_fix_20260815",
            "file": "/data/software/env/luly25/nanomelt/lib/python3.10/site-packages/ImmuneBuilder/refine.py",
            "old_expression": "{'Threads', str(n_threads)}",
            "new_expression": "{'Threads': str(n_threads)}",
        },
        "esm_cache": "/homes/Tianlab/luly25/.cache/torch/hub/checkpoints",
        "scoring": {"do_align": True, "ncpu": 1, "single_batch": True, "expected_rows": 47},
        "smoke_result": {
            "sample_uid": "LTT__Nb252",
            "input_length_aa": 128,
            "scored_length_aa": 126,
            "trimmed_reported_positions": "127-128",
            "trimmed_sequence": "GS",
            "predicted_apparent_tm_c": 65.18,
        },
        "primary_feature": "nanomelt_predicted_apparent_tm_c",
        "expected_direction": "higher_predicted_apparent_tm_higher_reported_yield",
        "primary_numeric_analysis": "31_LTT_WCC_individual_approximate",
        "llj_use": "16_group_ordinal_or_censored_observations_only",
        "sequence_cluster_identity": 0.9,
        "resampling": {"replicates": RESAMPLING_REPLICATES, "bootstrap_seed": BOOTSTRAP_SEED, "permutation_seed": BOOTSTRAP_SEED + 1},
        "influence_check": "leave_each_numeric_sample_out_and_report_without_Nb252",
        "high_capacity_model_training": False,
        "coverage_gate": {"planned": 47, "scoring_pass_required": 47, "numeric_required": 31, "llj_required": 16},
        "weak_ranking_gate": {
            "stratified_spearman_min": 0.30,
            "bootstrap_95ci_low_strictly_positive": True,
            "stratified_permutation_p_max": 0.05,
            "ltt_and_wcc_directions_positive": True,
            "scored_length_adjusted_direction_positive": True,
            "cluster_cv_increment_positive": True,
            "without_nb252_direction_positive": True,
        },
        "evidence_levels": ["weak_ranking_evidence", "compatibility_filter_only", "no_supported_use"],
        "release": "ready_for_remote_single_batch_nanomelt_scoring",
    }
    with tempfile.TemporaryDirectory(prefix=".nanomelt-plan-", dir=ROOT) as temp:
        stage = Path(temp)
        staged = {key: stage / Path(path).name for key, path in final.items()}
        _write_csv(staged["samples"], samples)
        staged["fasta"].write_text(
            "".join(f'>{row["sample_uid"]}\n{row["sequence_raw"]}\n' for row in samples),
            encoding="utf-8",
            newline="\n",
        )
        _write_json(staged["contract"], contract)
        _write_json(
            staged["summary"],
            {
                "schema_version": 1,
                "status": "pass",
                "generated_at": generated,
                "python": platform.python_version(),
                "elapsed_seconds": round(time.perf_counter() - started, 6),
                "sample_count": 47,
                "numeric_count": 31,
                "llj_ordinal_count": 16,
                "sequence_cluster_count_90": len({row["sequence_cluster_90"] for row in samples}),
                "outputs": {key: str(value) for key, value in final.items() if key != "summary"},
            },
        )
        replace_staged_files(
            {staged[key]: final[key] for key in staged},
            project_root=ROOT,
            protected_source_paths=valid.source_paths,
        )
    return 0


def _csv(path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _write_csv(path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    raise SystemExit(main())
