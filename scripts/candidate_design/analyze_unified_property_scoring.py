#!/usr/bin/env python3
"""Join unified NetSolP/NanoMelt scores and assign track-specific Pareto layers."""

from __future__ import annotations
import argparse,csv,json,platform,sys,tempfile,time
from datetime import datetime
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/"src"))
from antibody_optimization.file_transaction import replace_staged_files,validate_file_paths
from antibody_optimization.unified_property_plot import render_property_result
from antibody_optimization.unified_property_scoring import build_property_evidence

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument("--plan-dir",type=Path,required=True);p.add_argument("--unified-antifold-dir",type=Path,required=True);p.add_argument("--netsolp-score-dir",type=Path,required=True);p.add_argument("--nanomelt-score-dir",type=Path,required=True);p.add_argument("--output-dir",type=Path,required=True);p.add_argument("--run-summary",type=Path,required=True);p.add_argument("--generated-at");a=p.parse_args();started=time.perf_counter();generated=a.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
 plan=a.plan_dir.resolve(strict=True);unified=a.unified_antifold_dir.resolve(strict=True);net=a.netsolp_score_dir.resolve(strict=True);melt=a.nanomelt_score_dir.resolve(strict=True)
 sources=[plan/"unified_property_samples.csv",plan/"unified_property_scoring_plan_gate.json",unified/"unified_single_mutant_antifold_evidence.csv",unified/"unified_single_mutant_antifold_gate.json",net/"netsolp_sample_scores.csv",net/"netsolp_model_run.json",melt/"nanomelt_sample_scores.csv",melt/"nanomelt_model_run.json"]
 for path in sources:path.resolve(strict=True)
 if _json(sources[1]).get("status")!="pass" or _json(sources[3]).get("status")!="pass" or _json(sources[5]).get("status")!="pass" or _json(sources[7]).get("status")!="pass":raise ValueError("One or more unified property inputs are not passed")
 rows,summaries,gate=build_property_evidence(_csv(sources[2]),_csv(sources[0]),_csv(sources[4]),_csv(sources[6]));gate["generated_at"]=generated
 output=a.output_dir.absolute();summary=a.run_summary.absolute()
 if output.exists() or summary.exists():raise FileExistsError("Refusing to overwrite unified property results")
 output.parent.mkdir(parents=True,exist_ok=True);summary.parent.mkdir(parents=True,exist_ok=True);names={"evidence":"unified_single_mutant_property_evidence.csv","summary":"unified_property_tier_counts.csv","gate":"unified_property_scoring_gate.json","png":"unified_property_scoring.png","svg":"unified_property_scoring.svg"}
 with tempfile.TemporaryDirectory(prefix=".property-analysis-",dir=ROOT) as tmp:
  stage=Path(tmp);staged={k:stage/v for k,v in names.items()};_write_csv(staged["evidence"],rows);_write_csv(staged["summary"],summaries);_write_json(staged["gate"],gate);render_property_result(rows,png_path=staged["png"],svg_path=staged["svg"]);rs=stage/"run_summary.json";_write_json(rs,{"schema_version":1,"status":gate["status"],"generated_at":generated,"elapsed_seconds":round(time.perf_counter()-started,6),"python":platform.python_version(),"counts":{"candidates":len(rows),"wt_controls":1,"tier_rows":len(summaries)},"release":gate["release"],"outputs":{k:str(output/v) for k,v in names.items()}});pairs={staged[k]:output/v for k,v in names.items()};pairs[rs]=summary;validate_file_paths(project_root=ROOT,source_paths=sources,target_paths=pairs.values())
  for path in pairs.values():path.parent.mkdir(parents=True,exist_ok=True)
  replace_staged_files(pairs,project_root=ROOT,protected_source_paths=sources)
 return 0
def _csv(path):
 with path.open("r",encoding="utf-8-sig",newline="") as h:return list(csv.DictReader(h))
def _json(path):return json.loads(path.read_text(encoding="utf-8-sig"))
def _write_csv(path,rows):
 with path.open("w",encoding="utf-8-sig",newline="") as h:w=csv.DictWriter(h,fieldnames=list(rows[0]),lineterminator="\n");w.writeheader();w.writerows(rows)
def _write_json(path,value):path.write_text(json.dumps(value,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8",newline="\n")
if __name__=="__main__":raise SystemExit(main())
