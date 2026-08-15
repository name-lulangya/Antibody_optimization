#!/usr/bin/env python3
"""Build IMGT-numbered structures and the eight-core AntiFold validation plan."""

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

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from antibody_optimization.antifold_validation import (  # noqa: E402
    build_core_candidate_panel,
    prepare_imgt_structure,
)
from antibody_optimization.file_transaction import replace_staged_files, validate_file_paths  # noqa: E402


OUTPUTS = {
    "environment": "antifold_environment_contract.json",
    "candidates": "antifold_candidate_panel.csv",
    "views": "antifold_structure_views.csv",
    "mapping": "antifold_structure_mapping.csv",
    "plan": "antifold_validation_plan.json",
    "gate": "antifold_validation_plan_gate.json",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage0-dir", type=Path, required=True)
    parser.add_argument("--structure-baseline-dir", type=Path, required=True)
    parser.add_argument("--critical-facts", type=Path, required=True)
    parser.add_argument("--affinity-core-dir", type=Path, required=True)
    parser.add_argument("--experimental-cif", type=Path, required=True)
    parser.add_argument("--af3-cif", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    parser.add_argument("--check_only", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    sources = {
        "stage2": args.stage0_dir / "stage2_design_contract.json",
        "stage0_preflight": args.stage0_dir / "stage2_preflight.json",
        "mapping": args.structure_baseline_dir / "nb252_sequence_structure_mapping.csv",
        "structure_manifest": args.structure_baseline_dir / "structure_baseline_manifest.json",
        "critical_facts": args.critical_facts,
        "core_modules": args.affinity_core_dir / "affinity_core_modules.csv",
        "core_gate": args.affinity_core_dir / "affinity_ensemble_core_gate.json",
        "experimental_cif": args.experimental_cif,
        "af3_cif": args.af3_cif,
    }
    for path in sources.values():
        path.resolve(strict=True)
    if _json(sources["stage0_preflight"]).get("status") != "pass":
        raise ValueError("Stage-0 preflight is not passed")
    if _json(sources["structure_manifest"]).get("status") != "pass":
        raise ValueError("Released structure baseline is not passed")
    stage2 = _json(sources["stage2"])
    critical = _json(sources["critical_facts"])
    candidates = build_core_candidate_panel(
        _csv(sources["core_modules"]), _json(sources["core_gate"]), stage2, critical
    )
    mapping_rows = _csv(sources["mapping"])
    _validate_candidate_mapping(candidates, mapping_rows)
    if args.check_only:
        print(json.dumps({"status": "pass", "core_candidate_count": len(candidates), "view_count": 3}))
        return 0

    output_dir = args.output_dir.absolute()
    run_summary = args.run_summary.absolute()
    if output_dir.exists() or run_summary.exists():
        raise FileExistsError("Refusing to overwrite AntiFold validation plan outputs")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    run_summary.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix=".antifold-plan-", dir=PROJECT_ROOT) as temp_name:
        staging = Path(temp_name)
        structures = staging / "structures"
        exp_vhh = structures / "nb252_experimental_vhh_imgt.pdb"
        exp_complex = structures / "nb252_experimental_complex_imgt.pdb"
        af3_vhh = structures / "nb252_af3_vhh_imgt.pdb"
        view_specs = [
            ("experimental_vhh_only", sources["experimental_cif"], "NK2R-252.pdb", "C", ["C"], exp_vhh, "experimental_observed_only"),
            ("experimental_complex_context", sources["experimental_cif"], "NK2R-252.pdb", "C", ["C", "R"], exp_complex, "experimental_complex_context"),
            ("af3_vhh_only", sources["af3_cif"], "fold_2r_252_nomg_model_0.cif", "A", ["A"], af3_vhh, "predicted_full_vhh_sensitivity"),
        ]
        views: list[dict[str, object]] = []
        for view_id, source, model_name, vhh_chain, chains, staged_path, evidence_scope in view_specs:
            meta = prepare_imgt_structure(
                source_path=source,
                source_model_name=model_name,
                vhh_chain=vhh_chain,
                retained_chains=chains,
                mapping_rows=mapping_rows,
                output_path=staged_path,
            )
            views.append({
                "view_id": view_id,
                "structure_path": (output_dir / "structures" / staged_path.name).relative_to(PROJECT_ROOT).as_posix(),
                "source_model_name": model_name,
                "vhh_chain": vhh_chain,
                "antigen_chain": "R" if view_id == "experimental_complex_context" else "",
                "evidence_scope": evidence_scope,
                "vhh_observed_residue_count": meta["vhh_observed_residue_count"],
                "unnumbered_terminal_residue_count_removed": meta[
                    "unnumbered_terminal_residue_count_removed"
                ],
                "missing_coordinates_completed": False,
                "experimental_predicted_coordinates_mixed": False,
                "structure_sha256": _sha256(staged_path),
            })

        mapping_out = _mapping_output(mapping_rows, views)
        environment = {
            "schema_version": 1,
            "status": "validated_by_remote_smoke",
            "environment_path": "/data/software/env/luly25/antifold",
            "python": "3.10.20",
            "packages": {
                "antifold": "0.3.1", "torch": "2.2.0", "torch_cuda_build": "12.1",
                "torch_geometric": "2.4.0", "torch_scatter": "2.1.2", "biopython": "1.83",
                "biotite": "0.38.0", "pygam": "0.9.1", "numpy": "1.26.4",
                "pandas": "2.3.3", "fsspec": "2026.7.0",
            },
            "model_path": "/homes/Tianlab/luly25/software/AntiFold/models/model.pt",
            "model_sha256": "d5c442fa0372c28f4d0026d2f551b6f8ba7e7a127cb6837813a88093ed233e9e",
            "loader_expected_path": "/data/software/env/luly25/antifold/lib/python3.10/site-packages/models/model.pt",
            "loader_expected_path_resolves_to": "/homes/Tianlab/luly25/software/AntiFold/models/model.pt",
            "num_threads": 0,
            "gpu_smoke": {"device": "NVIDIA A100-PCIE-40GB", "model_type": "GVPTransformerModel", "status": "pass"},
            "source_version": "0.3.1",
            "source_commit": None,
            "source_commit_note": "Unavailable because the installed source tree has no .git metadata.",
            "known_loader_behavior": "AntiFold 0.3.1 ignores ordinary checkpoint_path values and loads loader_expected_path; the deployed symlink is required.",
        }
        plan = {
            "schema_version": 1,
            "plan_name": "nb252_antifold_minimal_validation",
            "status": "pass",
            "generated_at": generated_at,
            "authoritative_parent": stage2["authoritative_parent"],
            "candidate_count": len(candidates),
            "evaluation_entity_count": len(candidates) + 1,
            "view_count": len(views),
            "views": [row["view_id"] for row in views],
            "scoring_definition": "mutant_log_probability_minus_wt_log_probability_at_same_imgt_position",
            "interpretation": "supportive structure-conditioned sequence compatibility; not affinity, stability, expression, yield, or experimental validation",
            "candidate_generation_performed": False,
            "candidate_filtering_applied": False,
            "inputs": {key: str(path) for key, path in sources.items()},
        }
        gate = {
            "schema_version": 1,
            "gate_name": "nb252_antifold_validation_plan",
            "status": "pass",
            "generated_at": generated_at,
            "core_candidate_count": len(candidates),
            "view_count": len(views),
            "all_core_positions_uniquely_mapped": True,
            "missing_coordinates_completed": False,
            "experimental_predicted_coordinates_mixed": False,
            "release": "ready_for_remote_antifold_wt_backbone_scoring",
        }
        staged_outputs = {key: staging / name for key, name in OUTPUTS.items()}
        _write_json(staged_outputs["environment"], environment)
        _write_csv(staged_outputs["candidates"], candidates)
        _write_csv(staged_outputs["views"], views)
        _write_csv(staged_outputs["mapping"], mapping_out)
        _write_json(staged_outputs["plan"], plan)
        _write_json(staged_outputs["gate"], gate)

        final_pairs: dict[Path, Path] = {}
        for key, staged_path in staged_outputs.items():
            final_pairs[staged_path] = output_dir / OUTPUTS[key]
        for staged_path in structures.iterdir():
            final_pairs[staged_path] = output_dir / "structures" / staged_path.name
        summary_stage = staging / "run_summary.json"
        _write_json(summary_stage, {
            "schema_version": 1, "status": "pass", "generated_at": generated_at,
            "elapsed_seconds": round(time.perf_counter() - started, 6), "python": platform.python_version(),
            "command_argv": [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
            "counts": {"core_candidates": len(candidates), "views": len(views)},
            "outputs": {path.name: str(target) for path, target in final_pairs.items()},
        })
        final_pairs[summary_stage] = run_summary
        targets = list(final_pairs.values())
        validated = validate_file_paths(
            project_root=PROJECT_ROOT, source_paths=list(sources.values()), target_paths=targets
        )
        for path in validated.target_paths:
            path.parent.mkdir(parents=True, exist_ok=True)
        replace_staged_files(final_pairs, project_root=PROJECT_ROOT, protected_source_paths=validated.source_paths)
    return 0


def _validate_candidate_mapping(candidates: list[dict[str, object]], mapping: list[dict[str, str]]) -> None:
    for candidate in candidates:
        index = int(candidate["sequence_index_1based"])
        rows = [row for row in mapping if int(row["sequence_index_1based"]) == index]
        for model_name in ("NK2R-252.pdb", "fold_2r_252_nomg_model_0.cif"):
            matched = [row for row in rows if row["source_model_name"] == model_name and row["coordinate_status"] == "observed"]
            if len(matched) != 1 or matched[0]["residue_aa"] != candidate["wt_residue"]:
                raise ValueError(f"Candidate mapping is not unique for {candidate['candidate_id']} in {model_name}")


def _mapping_output(mapping: list[dict[str, str]], views: list[dict[str, object]]) -> list[dict[str, object]]:
    output = []
    for view in views:
        for row in mapping:
            if row["source_model_name"] != view["source_model_name"]:
                continue
            output.append({
                "view_id": view["view_id"], "source_model_name": row["source_model_name"],
                "source_auth_asym_id": row["auth_asym_id"], "source_auth_seq_id": row["auth_seq_id"],
                "source_insertion_code": row["insertion_code"], "derived_vhh_chain": view["vhh_chain"],
                "reported_sequence_index_1based": row["sequence_index_1based"], "residue_aa": row["residue_aa"],
                "imgt_position_label": row["numbering_position_label"], "coordinate_status": row["coordinate_status"],
                "evidence_scope": view["evidence_scope"],
            })
    return output


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def _write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


if __name__ == "__main__":
    raise SystemExit(main())
