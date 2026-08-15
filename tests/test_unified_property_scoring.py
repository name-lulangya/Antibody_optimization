from __future__ import annotations
import csv,json
from pathlib import Path
from antibody_optimization.nanomelt_yield import normalize_nanomelt_scores
from antibody_optimization.netsolp_yield import normalize_netsolp_scores
from antibody_optimization.unified_property_plot import render_property_result
from antibody_optimization.unified_property_scoring import build_property_evidence,build_property_samples

ROOT=Path(__file__).resolve().parents[1]
UNIFIED=ROOT/"docs/result_artifacts/candidate_design/unified_single_mutant_antifold_20260815"
def _csv(path):
 with path.open("r",encoding="utf-8-sig",newline="") as h:return list(csv.DictReader(h))
def _json(path):return json.loads(path.read_text(encoding="utf-8-sig"))

def _samples_and_candidates():
 candidates=_csv(UNIFIED/"unified_single_mutant_antifold_evidence.csv")
 return candidates,build_property_samples(candidates,_json(UNIFIED/"unified_single_mutant_antifold_gate.json"))

def test_real_property_plan_contains_wt_and_released_tracks():
 candidates,samples=_samples_and_candidates()
 assert len(samples)==1963 and samples[0]["score_id"]=="LTT__Nb252__WT"
 assert sum(row["design_track"]=="affinity_existing_interface_scan" for row in samples)==432
 assert sum(row["design_track"]=="stability_developability_discovery" for row in samples)==1530
 parent=samples[0]["sequence_raw"]
 for row in samples[1:]:
  assert sum(a!=b for a,b in zip(parent,row["sequence_raw"],strict=True))==1

def test_generalized_property_normalizers_preserve_arbitrary_plan_size():
 samples=[{"sample_uid":"wt","sequence_raw":"AAAA"},{"sample_uid":"m1","sequence_raw":"AATA"}]
 net=normalize_netsolp_scores(samples,[{"sid":"wt","fasta":"AAAA","predicted_usability":.5,"predicted_solubility":.6},{"sid":"m1","fasta":"AATA","predicted_usability":.7,"predicted_solubility":.4}],expected_count=2)
 assert [row["sample_uid"] for row in net]==["wt","m1"]
 melt=normalize_nanomelt_scores(samples,[{"ID":"wt","Aligned Sequence":"AAAA","Sequence":"AAAA","NanoMelt Tm (C)":60},{"ID":"m1","Aligned Sequence":"AATA","Sequence":"AATA","NanoMelt Tm (C)":61}],expected_pass_count=2,expected_plan_count=2)
 assert [row["scoring_status"] for row in melt]==["pass","pass"]

def test_real_unified_property_join_and_pareto_contract(tmp_path):
 candidates,samples=_samples_and_candidates();net=[];melt=[]
 for number,row in enumerate(samples):
  sequence=row["sequence_raw"]
  net.append({"sample_uid":row["score_id"],"sequence_raw":sequence,"predicted_usability":.5+(number%17)/1000,"predicted_solubility":.6+(number%13)/1000,"scoring_status":"pass"})
  melt.append({"sample_uid":row["score_id"],"sequence_raw":sequence,"scored_ungapped_sequence":sequence[:-2],"scored_length_aa":126,"trimmed_n_terminal":"","trimmed_c_terminal":"GS","nanomelt_predicted_apparent_tm_c":65+(number%19)/10,"scoring_status":"pass"})
 rows,summaries,gate=build_property_evidence(candidates,samples,net,melt)
 assert gate["status"]=="pass" and len(rows)==1962
 assert gate["track_counts"]=={"affinity_existing_interface_scan":432,"stability_developability_discovery":1530}
 assert {row["preliminary_property_tier"] for row in rows}.issubset({"pareto_front_1","pareto_front_2","background"})
 assert sum(row["candidate_selection_performed"] for row in rows)==0
 assert sum(row["candidate_count"] for row in summaries)==1962
 render_property_result(rows,png_path=tmp_path/"result.png",svg_path=tmp_path/"result.svg")
 assert (tmp_path/"result.png").stat().st_size>0
 assert (tmp_path/"result.svg").stat().st_size>0

def test_slurm_contract_is_parallel_without_arrays_or_resume():
 wrapper=(ROOT/"scripts/candidate_design/submit_unified_property_scoring.sh").read_text(encoding="utf-8")
 assert "afterok:${NET_JOB}:${MELT_JOB}" in wrapper
 assert "--array" not in wrapper and "resume" not in wrapper.lower()
 for name in ("submit_unified_property_netsolp.slurm","submit_unified_property_nanomelt.slurm","submit_unified_property_analysis.slurm"):
  text=(ROOT/"scripts/candidate_design"/name).read_text(encoding="utf-8")
  assert "#SBATCH --partition=batch" in text and "#SBATCH --gres=gpu:1" in text
  assert text.index("conda activate")<text.index("set -u")
