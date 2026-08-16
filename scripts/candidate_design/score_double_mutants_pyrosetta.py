#!/usr/bin/env python3
"""Score all 86 Nb252 double mutants with three paired PyRosetta repeats."""
from __future__ import annotations
import argparse,csv,json,platform,sys,time
from collections import Counter
from datetime import datetime
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/"src"))
from antibody_optimization import pyrosetta_runtime as runtime
from antibody_optimization.affinity_scoring import PAIRED_FIELDS,SUMMARY_FIELDS,WT_CONTROL_FIELDS,build_paired_row,build_wt_control_row,summarize_paired_rows
from antibody_optimization.double_mutant_design import EXPECTED_DOUBLES,REPLICATES
from antibody_optimization.flex_ddg_runtime import locate_mutation_pose_index,mutation_neighborhood
from antibody_optimization.pyrosetta_import_gate import load_released_stage_inputs

NEIGHBORHOOD=8.0
EXTRA=["position_pair","mutation_pose_indices","mutation_neighborhood_pose_indices","combined_movable_pose_indices","combined_movable_residue_count"]
PLAN_FIELDS=["mutation_set","combination_track","mutation_a","mutation_b","position_a","position_b","source_role_a","source_role_b","antifold_additive_fixed_backbone_delta_log_probability"]
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument("--plan-dir",type=Path,required=True);p.add_argument("--stage0-dir",type=Path,required=True);p.add_argument("--structure-baseline-dir",type=Path,required=True);p.add_argument("--calibration-dir",type=Path,required=True);p.add_argument("--output-dir",type=Path,required=True);p.add_argument("--run-summary",type=Path,required=True);p.add_argument("--replicates",type=int,default=3);p.add_argument("--base-seed",type=int,default=8163000);p.add_argument("--expected-pyrosetta-version",default="2026.03");p.add_argument("--generated-at");a=p.parse_args();started=time.perf_counter()
 if a.replicates!=REPLICATES:raise ValueError("Double-mutant contract requires three repeats")
 plan=a.plan_dir.resolve(strict=True);cal=a.calibration_dir.resolve(strict=True);out=a.output_dir.absolute();summary_path=a.run_summary.absolute();gate=_json(plan/"double_mutant_plan_gate.json");candidates=_csv(plan/"double_mutant_candidates.csv")
 if gate.get("release")!="ready_for_unfiltered_double_mutant_scoring" or len(candidates)!=EXPECTED_DOUBLES:raise ValueError("Double-mutant plan is not released")
 if out.exists() or summary_path.exists():raise FileExistsError("Refusing to overwrite double-mutant PyRosetta outputs")
 selection=_json(cal/"selected_scoring_protocol.json");cal_gate=_json(cal/"pyrosetta_scoring_calibration_gate.json")
 if cal_gate.get("pyrosetta_affinity_scoring_release")!="pass":raise ValueError("PyRosetta calibration is not released")
 protocol=str(selection["selected_protocol"]);parameters=selection["protocol_parameters"][protocol];interface_definition=selection["local_interface_definition"];interface_indices={int(x) for x in interface_definition["local_pose_indices"]};contact_cutoff=float(interface_definition["contact_retention_cutoff_angstrom"]);coordinate_sd=float(parameters["coordinate_constraint_sd_angstrom"])
 structure_inputs=load_released_stage_inputs(stage0_dir=a.stage0_dir.resolve(strict=True),structure_baseline_dir=a.structure_baseline_dir.resolve(strict=True));contact_rows=_csv(cal/"selected_contact_changes.csv");reference_contacts={"chain_a_auth_positions":{int(r["auth_seq_id"]) for r in contact_rows if r["chain_id"]=="C" and r["reference_contact"]=="True"},"chain_b_auth_positions":{int(r["auth_seq_id"]) for r in contact_rows if r["chain_id"]=="R" and r["reference_contact"]=="True"}}
 pyrosetta=runtime.initialize_pyrosetta(expected_version=a.expected_pyrosetta_version);starting_pose=pyrosetta.pose_from_file(str(cal/"selected_wt_prepared.pdb"));runtime.assert_pose_safety(starting_pose,structure_inputs);scorefxn=pyrosetta.create_score_function(runtime.SCORE_FUNCTION);reference_ca=runtime.ca_coordinates(starting_pose,interface_indices)
 groups={}
 for row in candidates:groups.setdefault((int(row["position_a"]),int(row["position_b"])),[]).append(row)
 wt_rows=[];paired=[];group_total=len(groups)
 for group_number,(position_pair,group) in enumerate(sorted(groups.items()),1):
  first=group[0];targets=[];locals_union=set()
  for suffix in ("a","b"):
   index=locate_mutation_pose_index(starting_pose,chain_id=first[f"experimental_chain_{suffix}"],auth_seq_id=int(first[f"experimental_auth_seq_{suffix}"]),insertion_code=first[f"experimental_insertion_{suffix}"]);targets.append(index);locals_union.update(mutation_neighborhood(starting_pose,index,NEIGHBORHOOD))
  movable=interface_indices|locals_union;position_text=f"{position_pair[0]};{position_pair[1]}";target_text=";".join(map(str,targets));local_text=";".join(map(str,sorted(locals_union)));movable_text=";".join(map(str,sorted(movable)));print(f"Position pair {group_number}/{group_total}: {position_text}, variants={len(group)}, movable={len(movable)}",flush=True)
  for replicate in range(1,REPLICATES+1):
   seed=a.base_seed+position_pair[0]*10000+position_pair[1]*10+replicate;wt_id=f"Nb252_WT_pair{position_pair[0]:03d}_{position_pair[1]:03d}_rep{replicate:02d}_seed{seed}";wt_pose=runtime.prepare_interface_pose(starting_pose,scorefxn,local_indices=movable,protocol=protocol,seed=seed,coordinate_constraint_sd=coordinate_sd);wt_metrics=runtime.measure_interface_pose(wt_pose,scorefxn,structure_inputs=structure_inputs,local_indices=interface_indices,reference_ca=reference_ca,reference_contacts=reference_contacts,protocol=protocol,replicate=replicate,seed=seed,contact_cutoff=contact_cutoff,include_contact_sets=True);wt=build_wt_control_row(replicate=replicate,seed=seed,metrics=wt_metrics);wt.update({"wt_control_id":wt_id,"position_pair":position_text,"mutation_pose_indices":target_text,"mutation_neighborhood_pose_indices":local_text,"combined_movable_pose_indices":movable_text,"combined_movable_residue_count":len(movable)});wt_rows.append(wt)
   for variant_number,candidate in enumerate(group,1):
    print(f"  repeat {replicate}/3 variant {variant_number}/{len(group)} {candidate['mutation_set']}",flush=True);pose=starting_pose.clone();allowed={}
    for suffix in ("a","b"):
     key=(candidate[f"experimental_chain_{suffix}"],int(candidate[f"experimental_auth_seq_{suffix}"]),candidate[f"experimental_insertion_{suffix}"]);allowed[key]=candidate[f"mutant_{suffix}"];runtime.mutate_pose_residue(pose,chain_id=key[0],auth_seq_id=key[1],insertion_code=key[2],wt_residue=candidate[f"wt_{suffix}"],mutant_residue=candidate[f"mutant_{suffix}"])
    runtime.assert_pose_safety(pose,structure_inputs,allowed_mutations=allowed);pose=runtime.prepare_interface_pose(pose,scorefxn,local_indices=movable,protocol=protocol,seed=seed,coordinate_constraint_sd=coordinate_sd);metrics=runtime.measure_interface_pose(pose,scorefxn,structure_inputs=structure_inputs,local_indices=interface_indices,reference_ca=reference_ca,reference_contacts=reference_contacts,protocol=protocol,replicate=replicate,seed=seed,contact_cutoff=contact_cutoff,allowed_mutations=allowed,include_contact_sets=True)
    adapted={"candidate_id":candidate["candidate_id"],"mutation_reported_label":candidate["mutation_set"],"mutation_numbering_label":candidate["mutation_set"],"mutation_source_auth_label":f"C:{candidate['experimental_auth_seq_a']}{candidate['mutant_a']};C:{candidate['experimental_auth_seq_b']}{candidate['mutant_b']}","sequence_index_1based":position_text,"wt_residue":f"{candidate['wt_a']};{candidate['wt_b']}","mutant_residue":f"{candidate['mutant_a']};{candidate['mutant_b']}","region":f"{candidate['region_a']};{candidate['region_b']}","prepared_contact_sensitive":True};row=build_paired_row(adapted,replicate=replicate,seed=seed,wt_metrics=wt_metrics,mutant_metrics=metrics);row.update({"wt_control_id":wt_id,"position_pair":position_text,"mutation_pose_indices":target_text,"mutation_neighborhood_pose_indices":local_text,"combined_movable_pose_indices":movable_text,"combined_movable_residue_count":len(movable)});paired.append(row)
 summaries=summarize_paired_rows(paired,expected_replicates=REPLICATES);by_id={r["candidate_id"]:r for r in candidates}
 for row in summaries:row.update({field:by_id[row["candidate_id"]][field] for field in PLAN_FIELDS})
 blockers=[]
 if len(wt_rows)!=len(groups)*3: blockers.append("wt_control_count_mismatch")
 if len(paired)!=258 or len(summaries)!=86: blockers.append("candidate_coverage_mismatch")
 if any(r["status"]!="pass" for r in paired): blockers.append("runtime_failure")
 run_gate={"schema_version":1,"gate_name":"nb252_double_mutant_pyrosetta_full_scan","status":"pass" if not blockers else "blocked","release":"ready_for_post_scan_joint_analysis" if not blockers else "blocked","blockers":blockers,"candidate_count":len(summaries),"position_pair_count":len(groups),"wt_control_count":len(wt_rows),"replicate_count":3,"mutant_evaluation_count":len(paired),"candidate_filtering_applied":False,"score_semantics":"double_mutant_minus_position_pair_matched_WT_Rosetta_ranking_signal","selected_protocol":protocol,"score_function":runtime.SCORE_FUNCTION,"mutation_neighborhood_angstrom":NEIGHBORHOOD}
 out.mkdir(parents=True);_write_csv(out/"double_mutant_wt_controls.csv",wt_rows,WT_CONTROL_FIELDS+EXTRA);_write_csv(out/"double_mutant_candidate_replicates.csv",paired,PAIRED_FIELDS+EXTRA);_write_csv(out/"double_mutant_candidate_summary.csv",summaries,SUMMARY_FIELDS+PLAN_FIELDS);_write_json(out/"double_mutant_scoring_gate.json",run_gate);summary_path.parent.mkdir(parents=True,exist_ok=True);_write_json(summary_path,{"schema_version":1,"status":run_gate["status"],"generated_at":a.generated_at or datetime.now().astimezone().isoformat(timespec="seconds"),"elapsed_seconds":round(time.perf_counter()-started,6),"python":platform.python_version(),"pyrosetta_version":pyrosetta.version(),"candidate_count":86,"position_pair_count":len(groups),"mutant_evaluation_count":258,"candidate_filtering_applied":False,"output_dir":str(out)});return 0 if not blockers else 2
def _csv(path):
 with path.open("r",encoding="utf-8-sig",newline="") as h:return list(csv.DictReader(h))
def _json(path):return json.loads(path.read_text(encoding="utf-8-sig"))
def _write_csv(path,rows,fields):
 with path.open("w",encoding="utf-8-sig",newline="") as h:w=csv.DictWriter(h,fieldnames=fields,lineterminator="\n");w.writeheader();w.writerows(rows)
def _write_json(path,value):path.write_text(json.dumps(value,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8",newline="\n")
if __name__=="__main__":raise SystemExit(main())
