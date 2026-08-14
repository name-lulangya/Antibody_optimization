#!/usr/bin/env python3
"""Build the fixed single-process 47-sequence TNP–BL21 validation plan."""

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
from antibody_optimization.tnp_yield import BOOTSTRAP_SEED, RESAMPLING_REPLICATES, build_tnp_validation_inputs  # noqa: E402


NAMES = {"samples": "tnp_validation_samples.csv", "fasta": "tnp_validation_sequences.fasta", "contract": "tnp_yield_validation_contract.json"}


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
    sources = [args.expression_records.resolve(strict=True), args.numbering_review.resolve(strict=True), args.numbering_positions.resolve(strict=True), args.allowed_use_manifest.resolve(strict=True)]
    allowed = _json(sources[3])
    if allowed["gates"]["cross_assay_pooling_gate"] != "pass" or allowed["counts"]["samples"] != 47:
        raise ValueError("Expression audit does not release the fixed 47-sample exploratory analysis")
    samples = build_tnp_validation_inputs(_csv(sources[0]), _csv(sources[1]), _csv(sources[2]))["sample_rows"]
    if args.check_only:
        print(json.dumps({"status": "pass", "samples": len(samples), "sequence_clusters_90": len({row["sequence_cluster_90"] for row in samples})}, ensure_ascii=False))
        return 0
    targets = [args.output_dir.absolute() / name for name in NAMES.values()] + [args.run_summary.absolute()]
    valid = validate_file_paths(project_root=ROOT, source_paths=sources, target_paths=targets)
    if any(path.exists() for path in valid.target_paths):
        raise FileExistsError("Refusing to overwrite TNP validation plan")
    for path in valid.target_paths:
        path.parent.mkdir(parents=True, exist_ok=True)
    final = dict(zip((*NAMES, "summary"), valid.target_paths, strict=True))
    contract = {
        "schema_version": 1, "status": "pass", "generated_at": generated, "expression_system": "BL21_E_coli",
        "software": {"tnp": "0.0.1", "tnp_commit": "29dcac72f1380e8538e8870f45a699d3c6156162", "immune_builder": "1.2", "anarci": "2024.05.21", "biopython": "1.77", "openmm": "8.5.2", "dssp": "4.6.1", "torch": "2.7.1+cu126"},
        "remote_environment": "/data/software/env/luly25/tnp", "remote_source": "/homes/Tianlab/luly25/software/TNP",
        "remote_entry_point": "/data/software/env/luly25/tnp/bin/TNP", "hydrophobicity_scale": {"argument": 0, "name": "Kyte-Doolittle"},
        "required_pythonpath": "/homes/Tianlab/luly25/software/TNP",
        "nanobodybuilder2_weight_md5": {"nanobody_model_1": "4591075f467ca9f76a37a5d1d3cfe591", "nanobody_model_2": "620fc916720bc7068cd18c7afa1aea8d", "nanobody_model_3": "f1ef0a66a54efb9d14ddbb97fdc30785", "nanobody_model_4": "a7e14a33f4c00e96df896120e9cf8522"},
        "execution": {"slurm_jobs": 1, "processes": 1, "samples_sequential": 47, "tnp_ncores": 1, "array": False},
        "primary_score": "tnp_psh", "expected_primary_direction": "higher_PSH_associated_with_lower_yield",
        "secondary_scores": ["tnp_total_cdr_length", "tnp_cdr3_length", "tnp_cdr3_compactness", "tnp_ppc", "tnp_pnc"],
        "coverage_gate": {"total_pass_min": 46, "numeric_pass_min": 30, "require_ltt_and_wcc": True},
        "sequence_cluster_identity": 0.9, "resampling": {"replicates": RESAMPLING_REPLICATES, "seed": BOOTSTRAP_SEED},
        "high_capacity_model_training": False, "release": "ready_for_remote_single_process_tnp_scoring",
    }
    with tempfile.TemporaryDirectory(prefix=".tnp-plan-", dir=ROOT) as temp:
        stage = Path(temp)
        staged = {key: stage / Path(path).name for key, path in final.items()}
        _write_csv(staged["samples"], samples)
        staged["fasta"].write_text("".join(f'>{row["sample_uid"]}\n{row["sequence_raw"]}\n' for row in samples), encoding="utf-8", newline="\n")
        _write_json(staged["contract"], contract)
        _write_json(staged["summary"], {"schema_version": 1, "status": "pass", "generated_at": generated, "python": platform.python_version(), "elapsed_seconds": round(time.perf_counter() - started, 6), "sample_count": 47, "numeric_count": 31, "llj_ordinal_count": 16, "sequence_cluster_count_90": len({row["sequence_cluster_90"] for row in samples}), "execution": "one_slurm_job_one_process_sequential_47", "outputs": {key: str(value) for key, value in final.items() if key != "summary"}})
        replace_staged_files({staged[key]: final[key] for key in staged}, project_root=ROOT, protected_source_paths=valid.source_paths)
    return 0


def _csv(path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle: return list(csv.DictReader(handle))
def _json(path): return json.loads(path.read_text(encoding="utf-8-sig"))
def _write_csv(path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n"); writer.writeheader(); writer.writerows(rows)
def _write_json(path, value): path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__": raise SystemExit(main())
