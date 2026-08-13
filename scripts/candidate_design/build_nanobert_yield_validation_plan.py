#!/usr/bin/env python3
"""Build the frozen 47-sequence nanoBERT–yield validation plan."""

from __future__ import annotations

import argparse, csv, json, platform, sys, tempfile, time
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from antibody_optimization.file_transaction import replace_staged_files, validate_file_paths  # noqa: E402
from antibody_optimization.nanobert_yield import (  # noqa: E402
    BOOTSTRAP_REPLICATES, BOOTSTRAP_SEED, PERMUTATION_REPLICATES, build_validation_inputs,
)

NAMES = {"samples": "nanobert_validation_samples.csv", "regions": "nanobert_validation_regions.csv", "contract": "nanobert_yield_validation_contract.json"}


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
    started = time.perf_counter(); generated = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    sources = [args.expression_records.resolve(strict=True), args.numbering_review.resolve(strict=True), args.numbering_positions.resolve(strict=True), args.allowed_use_manifest.resolve(strict=True)]
    allowed = _json(sources[3])
    if allowed["gates"]["cross_assay_pooling_gate"] != "pass" or allowed["counts"]["samples"] != 47:
        raise ValueError("Expression audit does not release the fixed 47-sample exploratory analysis")
    result = build_validation_inputs(_csv(sources[0]), _csv(sources[1]), _csv(sources[2]))
    if args.check_only:
        print(json.dumps({"status": "pass", "samples": len(result["sample_rows"]), "region_rows": len(result["region_rows"])}, ensure_ascii=False)); return 0
    targets = [args.output_dir.absolute() / name for name in NAMES.values()] + [args.run_summary.absolute()]
    valid = validate_file_paths(project_root=ROOT, source_paths=sources, target_paths=targets)
    if any(path.exists() for path in valid.target_paths): raise FileExistsError("Refusing to overwrite nanoBERT validation plan")
    for path in valid.target_paths: path.parent.mkdir(parents=True, exist_ok=True)
    final = dict(zip((*NAMES, "summary"), valid.target_paths, strict=True))
    contract = {
        "schema_version": 1, "status": "pass", "generated_at": generated,
        "model_id": "NaturalAntibody/nanoBERT", "model_revision": "edc8182ad89a827f8737fa572c6b5fac6197e6b0",
        "model_license": "cc-by-nc-sa-4.0",
        "score_definition": "per-residue single-mask pseudo-log-likelihood; natural log of WT residue probability",
        "primary_score": "mean_pll_over_complete_reported_sequence", "primary_numeric_analysis": "31_LTT_WCC_individual_approximate",
        "source_control": "within-provider_stratified_ranks_and_provider_intercept_LOOCV",
        "llj_use": "16_group_ordinal_or_censored_observations_only", "high_capacity_model_training": False,
        "bootstrap": {"replicates": BOOTSTRAP_REPLICATES, "seed": BOOTSTRAP_SEED},
        "permutation": {"replicates": PERMUTATION_REPLICATES, "seed": BOOTSTRAP_SEED + 1},
        "evidence_levels": ["weak_ranking_evidence", "compatibility_filter_only", "no_supported_use"],
        "release": "ready_for_remote_nanobert_scoring",
    }
    with tempfile.TemporaryDirectory(prefix=".nanobert-plan-", dir=ROOT) as temp:
        stage = Path(temp); staged = {key: stage / Path(path).name for key, path in final.items()}
        _write_csv(staged["samples"], result["sample_rows"]); _write_csv(staged["regions"], result["region_rows"]); _write_json(staged["contract"], contract)
        _write_json(staged["summary"], {"schema_version": 1, "status": "pass", "generated_at": generated, "python": platform.python_version(), "elapsed_seconds": round(time.perf_counter()-started, 6), "sample_count": 47, "numeric_count": 31, "llj_ordinal_count": 16, "outputs": {k: str(v) for k,v in final.items() if k != "summary"}})
        replace_staged_files({staged[k]: final[k] for k in staged}, project_root=ROOT, protected_source_paths=valid.source_paths)
    return 0


def _csv(path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle: return list(csv.DictReader(handle))
def _json(path): return json.loads(path.read_text(encoding="utf-8-sig"))
def _write_csv(path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer=csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n"); writer.writeheader(); writer.writerows(rows)
def _write_json(path, value): path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True)+"\n", encoding="utf-8", newline="\n")
if __name__ == "__main__": raise SystemExit(main())
