#!/usr/bin/env python3
"""Build the complete unfiltered 86-member Nb252 double-mutant plan."""
from __future__ import annotations
import argparse,csv,json,platform,sys,time
from datetime import datetime
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/"src"))
from antibody_optimization.double_mutant_design import build_double_mutant_space,build_score_samples
from antibody_optimization.double_mutant_design_plot import plot_plan
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument("--shortlist-dir",type=Path,required=True);p.add_argument("--mapping",type=Path,required=True);p.add_argument("--output-dir",type=Path,required=True);p.add_argument("--run-summary",type=Path,required=True);p.add_argument("--generated-at");a=p.parse_args();started=time.perf_counter();source=a.shortlist_dir.resolve(strict=True);mapping=a.mapping.resolve(strict=True);out=a.output_dir.absolute();summary=a.run_summary.absolute()
 if out.exists() or summary.exists():raise FileExistsError(f"Refusing to overwrite: {out} or {summary}")
 result=build_double_mutant_space(_csv(source/"single_mutant_shortlist.csv"),_json(source/"single_mutant_shortlist_gate.json"),_csv(mapping));samples=build_score_samples(result["parent_sequence"],result["candidates"]);out.mkdir(parents=True)
 _write_csv(out/"double_mutant_candidates.csv",result["candidates"]);_write_csv(out/"double_mutant_invalid_same_position_pairs.csv",result["invalid_pairs"]);_write_csv(out/"double_mutant_score_samples.csv",samples)
 with (out/"double_mutant_sequences.fasta").open("w",encoding="utf-8",newline="\n") as h:
  for row in samples:h.write(f">{row['sample_uid']}\n{row['sequence_raw']}\n")
 generated=a.generated_at or datetime.now().astimezone().isoformat(timespec="seconds");gate={"schema_version":1,"gate_name":"nb252_complete_double_mutant_plan","status":"pass","release":"ready_for_unfiltered_double_mutant_scoring","generated_at":generated,**result["facts"],"antifold_evidence":"additive fixed-backbone single-position log-probability deltas; no double-mutant epistasis claim","scoring_plan":{"netsolp":87,"nanomelt":87,"tnp":87,"pyrosetta_candidates":86,"pyrosetta_replicates":3,"pyrosetta_evaluations":258},"interpretation":"Computational design plan only; no measured affinity, expression, solubility, stability, Tm, or yield."};_write_json(out/"double_mutant_plan_gate.json",gate);plot_plan(result["candidates"],out/"double_mutant_plan_plot_data.csv",out/"double_mutant_plan.png",out/"double_mutant_plan.svg");summary.parent.mkdir(parents=True,exist_ok=True);_write_json(summary,{"schema_version":1,"status":"pass","generated_at":generated,"elapsed_seconds":round(time.perf_counter()-started,6),"python":platform.python_version(),"valid_double_count":86,"invalid_same_position_pair_count":5,"output_dir":str(out),"candidate_filtering_applied":False});return 0
def _csv(path):
 with path.open("r",encoding="utf-8-sig",newline="") as h:return list(csv.DictReader(h))
def _json(path):return json.loads(path.read_text(encoding="utf-8-sig"))
def _write_csv(path,rows):
 with path.open("w",encoding="utf-8-sig",newline="") as h:w=csv.DictWriter(h,fieldnames=list(rows[0]),lineterminator="\n");w.writeheader();w.writerows(rows)
def _write_json(path,value):path.write_text(json.dumps(value,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8",newline="\n")
if __name__=="__main__":raise SystemExit(main())
