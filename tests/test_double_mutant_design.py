import csv,json,subprocess,sys,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/"src"))
from antibody_optimization.double_mutant_design import build_double_mutant_space,build_score_samples
from antibody_optimization.double_mutant_analysis import build_joint_evidence
SHORT=ROOT/"docs/result_artifacts/candidate_design/single_mutant_shortlist_20260816";MAPPING=ROOT/"docs/result_artifacts/input_baseline/structure_released_20260810/nb252_sequence_structure_mapping.csv"
def _csv(p):
 with p.open("r",encoding="utf-8-sig",newline="") as h:return list(csv.DictReader(h))
def test_real_shortlist_builds_all_valid_pairs_without_filtering():
 result=build_double_mutant_space(_csv(SHORT/"single_mutant_shortlist.csv"),json.loads((SHORT/"single_mutant_shortlist_gate.json").read_text(encoding="utf-8")),_csv(MAPPING));assert result["facts"]["combination_track_counts"]=={"affinity_x_affinity":26,"affinity_x_property":48,"property_x_property":12};assert len(result["candidates"])==86 and len(result["invalid_pairs"])==5;assert len(build_score_samples(result["parent_sequence"],result["candidates"]))==87;assert all(row["candidate_filtering_applied"] is False for row in result["candidates"]);assert {row["mutation_set"] for row in result["candidates"]}>={"Q1D;S55G","F30A;D101W","D101W;I103W"}
def test_cli_writes_plan_and_figure():
 with tempfile.TemporaryDirectory(prefix=".test-double-",dir=ROOT) as temp:
  out=Path(temp)/"out";summary=Path(temp)/"run.json";subprocess.run([sys.executable,str(ROOT/"scripts/candidate_design/build_double_mutant_plan.py"),"--shortlist-dir",str(SHORT),"--mapping",str(MAPPING),"--output-dir",str(out),"--run-summary",str(summary),"--generated-at","2026-08-16T20:00:00+08:00"],check=True,cwd=ROOT);gate=json.loads((out/"double_mutant_plan_gate.json").read_text(encoding="utf-8"));assert gate["valid_double_count"]==86 and gate["status"]=="pass";assert len(_csv(out/"double_mutant_candidates.csv"))==86 and (out/"double_mutant_plan.png").stat().st_size>1000 and summary.is_file()
def test_joint_analysis_requires_complete_sets_and_uses_existing_magnitude_rules():
 plan=ROOT/"docs/result_artifacts/candidate_design/double_mutant_plan_20260816";candidates=_csv(plan/"double_mutant_candidates.csv");samples=_csv(plan/"double_mutant_score_samples.csv");net=[];melt=[];tnp=[]
 for sample in samples:
  uid=sample["sample_uid"];is_wt=uid=="LTT__Nb252__WT";net.append({**sample,"predicted_usability":0.5 if is_wt else 0.515,"predicted_solubility":0.5,"scoring_status":"pass"});melt.append({**sample,"nanomelt_predicted_apparent_tm_c":65.0 if is_wt else 65.5,"scoring_status":"pass"});tnp.append({**sample,"tnp_psh":140.0,"tnp_flag_total_cdr_length":"green","tnp_flag_cdr3_length":"green","tnp_flag_cdr3_compactness":"green","tnp_flag_psh":"amber","tnp_flag_ppc":"green","tnp_flag_pnc":"green","scoring_status":"pass"})
 rose=[{"candidate_id":row["candidate_id"],"delta_dG_separated_median":-1.0,"delta_cross_interface_energy_median":-0.5} for row in candidates];rows,gate=build_joint_evidence(candidates,net,melt,tnp,rose);assert len(rows)==86 and gate["status"]=="pass";assert "balanced_supported" in {row["joint_evidence_class"] for row in rows};assert all(row["netsolp_usability_magnitude"]=="favorable" for row in rows);assert all(row["joint_evidence_class"]=="tradeoff_or_no_clear_joint_support" for row in rows if row["new_liability_flags"])
