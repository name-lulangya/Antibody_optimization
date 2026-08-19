#!/usr/bin/env python3
"""Build the fixed 47-sequence PLM_Sol--BL21 yield validation plan."""

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


NAMES = {
    "samples": "plm_sol_validation_samples.csv",
    "fasta": "plm_sol_validation_sequences.fasta",
    "contract": "plm_sol_yield_validation_contract.json",
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-plan-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    parser.add_argument("--check_only", action="store_true")
    args = parser.parse_args()
    started = time.perf_counter()
    generated = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    reference = args.reference_plan_dir.resolve(strict=True)
    source_samples = reference / "rp3net_validation_samples.csv"
    source_fasta = reference / "rp3net_validation_sequences.fasta"
    source_contract = reference / "rp3net_yield_validation_contract.json"
    samples = _csv(source_samples)
    contract = _json(source_contract)
    if len(samples) != 47 or len({row["sample_uid"] for row in samples}) != 47:
        raise ValueError("Reference validation plan must contain 47 unique samples")
    if len({row["sequence_raw"] for row in samples}) != 47:
        raise ValueError("PLM_Sol sequence-based output mapping requires 47 unique sequences")
    expected_fasta = "".join(f'>{row["sample_uid"]}\n{row["sequence_raw"]}\n' for row in samples)
    if source_fasta.read_text(encoding="utf-8") != expected_fasta:
        raise ValueError("Reference FASTA does not match the fixed sample table")
    if contract.get("status") != "pass" or contract.get("expression_system") != "BL21_E_coli":
        raise ValueError("Reference yield-validation contract is not released")
    if args.check_only:
        print(json.dumps({"status": "pass", "samples": 47}, ensure_ascii=False))
        return 0

    targets = [args.output_dir.absolute() / name for name in NAMES.values()] + [args.run_summary.absolute()]
    sources = [source_samples, source_fasta, source_contract]
    valid = validate_file_paths(project_root=ROOT, source_paths=sources, target_paths=targets)
    if any(path.exists() for path in valid.target_paths):
        raise FileExistsError("Refusing to overwrite PLM_Sol validation plan")
    for path in valid.target_paths:
        path.parent.mkdir(parents=True, exist_ok=True)
    final = dict(zip((*NAMES, "summary"), valid.target_paths, strict=True))
    plm_contract = {
        "schema_version": 1,
        "status": "pass",
        "generated_at": generated,
        "expression_system": "BL21_E_coli",
        "sample_semantics": {"total": 47, "numeric_individual": 31, "llj_ordinal_censored": 16},
        "software": {
            "name": "PLM_Sol",
            "release": "V1.0",
            "source_root": "/homes/Tianlab/luly25/software/PLM_Sol/V1.0/Violet969-PLM_Sol-5156de9",
            "classifier_checkpoint": "/homes/Tianlab/luly25/software/PLM_Sol/V1.0/Violet969-PLM_Sol-5156de9/model_param/model_param.t7",
            "embedding_model": "/homes/Tianlab/luly25/software/PLM_Sol_models/prottrans_t5_xl_u50",
        },
        "environments": {
            "embedding": "/data/software/env/luly25/plm_sol_embed",
            "classifier": "/data/software/env/luly25/plm_sol",
            "analysis": "/data/software/env/luly25/ab_optim",
        },
        "embedding": {"protocol": "prottrans_t5_xl_u50", "per_residue_dimension": 1024, "storage_dtype": "float16"},
        "classifier": {"primary_score": "predict_result", "direction": "higher_is_better", "range": [0, 1]},
        "vendor_fixes": [
            "disable_missing_models_legacy_import",
            "resolve_train_arguments_yml_from_checkpoint_directory",
            "call_solver_predict_evaluation_with_dataset_only",
            "import_torch_pad_sequence",
            "install_pyaml_21_10_1",
        ],
        "classification": {
            "primary_outcome": "matching_provider_outer_training_fold_median",
            "score_threshold": "maximize_training_MCC_then_balanced_accuracy_then_higher_threshold",
            "outer_schemes": ["leave_one_out", "leave_one_cluster_out"],
            "fixed_5mg": "display_only_not_a_predictor_gate",
        },
        "comparators": ["NetSolP_S", "NetSolP_U", "RP3Net_historical"],
        "runtime_estimate": "10_to_30_minutes_and_less_than_1_hour",
        "release": "ready_for_remote_plm_sol_scoring",
    }
    with tempfile.TemporaryDirectory(prefix=".plm-sol-plan-", dir=ROOT) as temp:
        stage = Path(temp)
        staged = {key: stage / Path(path).name for key, path in final.items()}
        _write_csv(staged["samples"], samples)
        staged["fasta"].write_text(expected_fasta, encoding="utf-8", newline="\n")
        _write_json(staged["contract"], plm_contract)
        _write_json(staged["summary"], {
            "schema_version": 1, "status": "pass", "generated_at": generated,
            "python": platform.python_version(), "elapsed_seconds": round(time.perf_counter() - started, 6),
            "counts": {"samples": 47, "numeric": 31, "llj_ordinal_censored": 16},
            "outputs": {key: str(value) for key, value in final.items() if key != "summary"},
        })
        replace_staged_files({staged[key]: final[key] for key in staged}, project_root=ROOT, protected_source_paths=valid.source_paths)
    return 0


def _csv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _write_csv(path: Path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def _write_json(path: Path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    raise SystemExit(main())
