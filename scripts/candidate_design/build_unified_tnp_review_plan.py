#!/usr/bin/env python3
"""Freeze WT plus 95 released Nb252 candidates for one TNP review run."""

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
from antibody_optimization.file_transaction import replace_staged_files, validate_file_paths  # noqa: E402
from antibody_optimization.unified_tnp_review import (  # noqa: E402
    AFFINITY_COUNT,
    BLOCKED_PRODUCTION_CYS_IDS,
    MAGNITUDE_THRESHOLDS,
    PROPERTY_COUNT,
    SCORE_COUNT,
    build_unified_tnp_samples,
)


NAMES = {
    "samples": "unified_tnp_samples.csv",
    "fasta": "unified_tnp_sequences.fasta",
    "blocked": "unified_tnp_blocked_audit.csv",
    "contract": "unified_tnp_review_contract.json",
    "gate": "unified_tnp_review_plan_gate.json",
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--property-result-dir", type=Path, required=True)
    parser.add_argument("--property-plan-dir", type=Path, required=True)
    parser.add_argument("--flex-ddg-result-dir", type=Path, required=True)
    parser.add_argument("--unified-plan-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    parser.add_argument("--check_only", action="store_true")
    args = parser.parse_args()
    started = time.perf_counter()
    generated = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")

    property_result = args.property_result_dir.resolve(strict=True)
    property_plan = args.property_plan_dir.resolve(strict=True)
    flex_result = args.flex_ddg_result_dir.resolve(strict=True)
    unified_plan = args.unified_plan_dir.resolve(strict=True)
    sources = [
        property_result / "unified_single_mutant_property_evidence.csv",
        property_result / "unified_property_scoring_gate.json",
        property_result / "unified_property_scoring_scientific_review.json",
        property_plan / "unified_property_samples.csv",
        flex_result / "flex_ddg_production_candidate_summary.csv",
        flex_result / "flex_ddg_production_gate.json",
        flex_result / "flex_ddg_production_scientific_review.json",
        unified_plan / "unified_single_mutant_candidates.csv",
        unified_plan / "unified_single_mutant_plan_gate.json",
    ]
    for path in sources:
        path.resolve(strict=True)
    if _json(sources[1]).get("status") != "pass":
        raise ValueError("Unified property result gate is not passed")
    if _json(sources[2]).get("release") != "ready_for_preliminary_property_pool_definition":
        raise ValueError("Unified property scientific review is not released")
    if _json(sources[5]).get("status") != "pass" or int(_json(sources[5]).get("candidate_count", 0)) != 50:
        raise ValueError("Flex ddG production result is not the fixed 50-candidate review")
    if _json(sources[8]).get("status") != "pass":
        raise ValueError("Unified single-mutant plan is not passed")

    samples, blocked = build_unified_tnp_samples(
        _csv(sources[0]), _csv(sources[3]), _csv(sources[4]), _csv(sources[7])
    )
    source_counts = Counter(str(row["candidate_source"]) for row in samples[1:])
    if args.check_only:
        print(json.dumps({
            "status": "pass",
            "score_count": len(samples),
            "candidate_count": len(samples) - 1,
            "source_counts": dict(sorted(source_counts.items())),
            "blocked_audit_count": len(blocked),
        }, sort_keys=True))
        return 0

    output_dir = args.output_dir.absolute()
    run_summary = args.run_summary.absolute()
    targets = [output_dir / name for name in NAMES.values()] + [run_summary]
    valid = validate_file_paths(project_root=ROOT, source_paths=sources, target_paths=targets)
    if any(path.exists() for path in valid.target_paths):
        raise FileExistsError("Refusing to overwrite unified TNP review plan")
    for path in valid.target_paths:
        path.parent.mkdir(parents=True, exist_ok=True)
    final = dict(zip((*NAMES, "summary"), valid.target_paths, strict=True))
    contract = {
        "schema_version": 1,
        "contract_name": "nb252_unified_tnp_candidate_review",
        "status": "pass",
        "generated_at": generated,
        "planned_count": SCORE_COUNT,
        "candidate_count": SCORE_COUNT - 1,
        "wt_control_count": 1,
        "candidate_sources": {
            "property_pareto_front_1": PROPERTY_COUNT,
            "affinity_flex_ddg_20_sample_pool": AFFINITY_COUNT,
        },
        "blocked_audit_ids": sorted(BLOCKED_PRODUCTION_CYS_IDS),
        "magnitude_thresholds": {
            "delta_netsolp_usability": {"neutral_inclusive": [-0.01, 0.01]},
            "delta_netsolp_solubility": {"neutral_inclusive": [-0.02, 0.02]},
            "delta_nanomelt_predicted_apparent_tm_c": {"neutral_inclusive": [-1.0, 1.0]},
        },
        "magnitude_threshold_role": "operational_relative_to_scanned_library_not_calibrated_uncertainty",
        "primary_property_axes": ["netsolp_usability", "nanomelt_predicted_apparent_tm_c"],
        "supporting_property_axis": "netsolp_solubility",
        "compatibility_constraint": "experimental_complex_context_antifold_delta_log_probability",
        "risk_constraints": ["TNP_flags", "chemical_liability_flags"],
        "software": {
            "tnp": "0.0.1",
            "tnp_commit": "29dcac72f1380e8538e8870f45a699d3c6156162",
            "immune_builder": "1.2",
            "anarci": "2024.05.21",
            "biopython": "1.77",
            "openmm": "8.5.2",
            "dssp": "4.6.1",
            "torch": "2.7.1+cu126",
        },
        "remote_environment": "/data/software/env/luly25/tnp",
        "remote_source": "/homes/Tianlab/luly25/software/TNP",
        "remote_entry_point": "/data/software/env/luly25/tnp/bin/TNP",
        "required_pythonpath": "/homes/Tianlab/luly25/software/TNP",
        "immune_builder_refine": "/data/software/env/luly25/tnp/lib/python3.10/site-packages/ImmuneBuilder/refine.py",
        "hydrophobicity_scale": {"argument": 0, "name": "Kyte-Doolittle"},
        "input_length_aa": 128,
        "expected_modelled_length_aa": 126,
        "required_trimmed_c_terminal": "GS",
        "execution": {
            "slurm_jobs": 1,
            "processes": 1,
            "samples_sequential": SCORE_COUNT,
            "tnp_ncores": 1,
            "array": False,
            "resume": False,
        },
        "runtime_estimate": {
            "expected": "approximately 50-60 minutes from the prior 43-sequence real run",
            "likely_over_one_hour": False,
            "likely_over_five_hours": False,
        },
        "yield_prediction": False,
        "candidate_selection": False,
        "release": "ready_for_remote_single_process_unified_tnp_review",
    }
    if set(MAGNITUDE_THRESHOLDS.values()) != {0.01, 0.02, 1.0}:
        raise RuntimeError("Magnitude threshold contract drift")
    gate = {
        "schema_version": 1,
        "gate_name": "nb252_unified_tnp_review_plan",
        "status": "pass",
        "generated_at": generated,
        "planned_count": len(samples),
        "candidate_count": len(samples) - 1,
        "wt_control_count": 1,
        "source_counts": dict(sorted(source_counts.items())),
        "blocked_audit_count": len(blocked),
        "release": contract["release"],
    }

    with tempfile.TemporaryDirectory(prefix=".unified-tnp-plan-", dir=ROOT) as temp:
        stage = Path(temp)
        staged = {key: stage / Path(path).name for key, path in final.items()}
        _write_csv(staged["samples"], samples)
        _write_csv(staged["blocked"], blocked)
        staged["fasta"].write_text(
            "".join(f'>{row["sample_uid"]}\n{row["sequence_raw"]}\n' for row in samples),
            encoding="utf-8",
            newline="\n",
        )
        _write_json(staged["contract"], contract)
        _write_json(staged["gate"], gate)
        _write_json(staged["summary"], {
            "schema_version": 1,
            "status": "pass",
            "generated_at": generated,
            "python": platform.python_version(),
            "elapsed_seconds": round(time.perf_counter() - started, 6),
            "planned_count": len(samples),
            "candidate_count": len(samples) - 1,
            "source_counts": dict(sorted(source_counts.items())),
            "blocked_audit_count": len(blocked),
            "outputs": {key: str(path) for key, path in final.items() if key != "summary"},
        })
        replace_staged_files(
            {staged[key]: final[key] for key in staged},
            project_root=ROOT,
            protected_source_paths=valid.source_paths,
        )
    return 0


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    raise SystemExit(main())
