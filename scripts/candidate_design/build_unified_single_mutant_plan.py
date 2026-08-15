#!/usr/bin/env python3
"""Build the unified Nb252 contract and complete single-mutant universe."""

from __future__ import annotations

import argparse, csv, json, platform, sys, tempfile, time
from collections import Counter
from datetime import datetime
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from antibody_optimization.file_transaction import replace_staged_files, validate_file_paths
from antibody_optimization.unified_single_mutant_plot import render_unified_space
from antibody_optimization.unified_single_mutants import build_unified_space


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage0-dir", type=Path, required=True)
    parser.add_argument("--critical-facts", type=Path, required=True)
    parser.add_argument("--affinity-core-dir", type=Path, required=True)
    parser.add_argument("--affinity-single-mutant-dir", type=Path, required=True)
    parser.add_argument("--affinity-full-scan-dir", type=Path, required=True)
    parser.add_argument("--antifold-plan-dir", type=Path, required=True)
    parser.add_argument("--antifold-result-dir", type=Path, required=True)
    parser.add_argument("--historical-stability-contract-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at"); parser.add_argument("--check_only", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args(); started = time.perf_counter()
    generated = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    source = {
        "stage2": args.stage0_dir / "stage2_design_contract.json",
        "positions": args.stage0_dir / "mutable_position_inventory.csv",
        "critical": args.critical_facts,
        "core": args.affinity_core_dir / "affinity_core_modules.csv",
        "core_gate": args.affinity_core_dir / "affinity_ensemble_core_gate.json",
        "affinity_candidates": args.affinity_single_mutant_dir / "affinity_single_mutants.csv",
        "affinity_candidate_gate": args.affinity_single_mutant_dir / "affinity_candidate_gate.json",
        "affinity_summary": args.affinity_full_scan_dir / "candidate_summary.csv",
        "affinity_merge_gate": args.affinity_full_scan_dir / "full_scan_merge_gate.json",
        "antifold_plan_gate": args.antifold_plan_dir / "antifold_validation_plan_gate.json",
        "antifold_result_gate": args.antifold_result_dir / "antifold_validation_gate.json",
        "historical_stability_gate": args.historical_stability_contract_dir / "stability_expression_contract_gate.json",
    }
    for path in source.values(): path.resolve(strict=True)
    for key in ("core_gate", "affinity_candidate_gate", "affinity_merge_gate", "antifold_plan_gate", "antifold_result_gate", "historical_stability_gate"):
        if _json(source[key]).get("status") != "pass": raise ValueError(f"Upstream gate not passed: {key}")
    positions, candidates = build_unified_space(_json(source["stage2"]), _csv(source["positions"]), _json(source["critical"]), _csv(source["core"]), _csv(source["affinity_candidates"]), _csv(source["affinity_summary"]))
    counts = Counter(str(row["design_status"]) for row in candidates)
    if args.check_only:
        print(json.dumps({"status":"pass","position_count":len(positions),"candidate_count":len(candidates),"design_status_counts":dict(counts)}, sort_keys=True)); return 0
    output = args.output_dir.absolute(); summary = args.run_summary.absolute()
    if output.exists() or summary.exists(): raise FileExistsError("Refusing to overwrite unified single-mutant outputs")
    output.parent.mkdir(parents=True, exist_ok=True); summary.parent.mkdir(parents=True, exist_ok=True)
    names = {"contract":"unified_single_mutant_contract.json","positions":"unified_position_space.csv","candidates":"unified_single_mutant_candidates.csv","fasta":"unified_single_mutant_candidates.fasta","counts":"unified_candidate_status_counts.csv","gate":"unified_single_mutant_plan_gate.json","png":"unified_single_mutant_space.png","svg":"unified_single_mutant_space.svg"}
    contract = {
        "schema_version":1,"contract_name":"nb252_unified_single_mutant_design","status":"pass","generated_at":generated,
        "authoritative_parent":_json(source["stage2"])["authoritative_parent"],
        "scope":"all_19_non_wt_substitutions_at_each_non_hard_immutable_reported_position_with_separate_release_and_evidence_routes",
        "hard_constraints":["do_not_mutate_reported_positions_22_95_125_126_127_128","do_not_introduce_extra_unpaired_cysteine","validate_parent_wt_identity"],
        "policies":{"fr_cdr_partition_is_not_a_design_boundary":True,"experimental_interface":"reuse_existing_456_member_pyrosetta_scan_and_do_not_rescore","experimental_missing_coordinates":"audit_only_not_released_current_round","noninterface_observed":"stability_developability_discovery_then_shortlist_for_affinity_noninferiority","candidate_selection_performed":False},
        "antifold_reuse":{"primary_view":"experimental_complex_context","sensitivity_views":["experimental_vhh_only","af3_vhh_only"],"raw_score_source":"existing_antifold_validation_model_scores","model_rerun_required":False,"interpretation":"structure-conditioned sequence compatibility only"},
        "chemistry_features":{"formal_charge":"D/E=-1; K/R=+1; all others=0","hydrophobic_set":"AVILMFWY","motifs":["N-X-[ST], X != P","N-[GST]","D-[GST]","M/W count"],"semantics":"sequence heuristics; no scientific filtering applied"},
        "historical_stability_contract_role":"provenance_only_superseded_scope_not_active_filter",
        "pyrosetta_reuse":{"interface_candidate_count":456,"source_candidate_table":str(source["affinity_candidates"]),"source_result_table":str(source["affinity_summary"]),"model_rerun_required":False},
    }
    gate={"schema_version":1,"gate_name":"nb252_unified_single_mutant_plan","status":"pass","generated_at":generated,"position_count":128,"enumerated_position_count":122,"candidate_count":2318,"design_status_counts":dict(sorted(counts.items())),"design_track_counts":dict(sorted(Counter(str(r["design_track"]) for r in candidates).items())),"existing_affinity_scan_candidate_count":sum(_bool(r["existing_affinity_scan_candidate"]) for r in candidates),"affinity_core_count":sum(_bool(r["affinity_core_module"]) for r in candidates),"candidate_selection_performed":False,"scientific_score_threshold_applied":False,"pyrosetta_rescoring_performed":False,"antifold_model_rerun_required":False,"release":"ready_for_existing_antifold_landscape_join"}
    count_rows=[{"design_status":key,"candidate_count":counts[key]} for key in ("eligible_current_round","deferred_missing_experimental_coordinates","blocked_new_unpaired_cys")]
    with tempfile.TemporaryDirectory(prefix=".unified-plan-", dir=PROJECT_ROOT) as tmp:
        stage=Path(tmp); staged={k:stage/v for k,v in names.items()}
        _write_json(staged["contract"],contract); _write_csv(staged["positions"],positions); _write_csv(staged["candidates"],candidates); _write_csv(staged["counts"],count_rows); _write_json(staged["gate"],gate)
        with staged["fasta"].open("w",encoding="utf-8",newline="\n") as handle:
            for row in candidates: handle.write(f">{row['candidate_id']} design_status={row['design_status']}\n{row['sequence']}\n")
        render_unified_space(positions,candidates,png_path=staged["png"],svg_path=staged["svg"])
        summary_stage=stage/"run_summary.json"; _write_json(summary_stage,{"schema_version":1,"status":"pass","generated_at":generated,"elapsed_seconds":round(time.perf_counter()-started,6),"python":platform.python_version(),"command_argv":[sys.executable,str(Path(__file__).resolve()),*sys.argv[1:]],"counts":{"positions":128,"candidates":2318,**dict(counts)},"release":gate["release"],"outputs":{k:str(output/v) for k,v in names.items()}})
        pairs={staged[k]:output/v for k,v in names.items()}; pairs[summary_stage]=summary
        validate_file_paths(project_root=PROJECT_ROOT,source_paths=source.values(),target_paths=pairs.values())
        for path in pairs.values(): path.parent.mkdir(parents=True,exist_ok=True)
        replace_staged_files(pairs,project_root=PROJECT_ROOT,protected_source_paths=source.values())
    return 0


def _csv(path):
    with path.open("r",encoding="utf-8-sig",newline="") as h:return list(csv.DictReader(h))
def _json(path): return json.loads(path.read_text(encoding="utf-8-sig"))
def _write_json(path,value): path.write_text(json.dumps(value,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8",newline="\n")
def _write_csv(path,rows):
    with path.open("w",encoding="utf-8-sig",newline="") as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0]),lineterminator="\n");w.writeheader();w.writerows(rows)
def _bool(value): return value is True or str(value).lower()=="true"
if __name__=="__main__": raise SystemExit(main())
