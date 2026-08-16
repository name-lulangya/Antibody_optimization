#!/usr/bin/env python3
"""Score WT plus 86 double mutants with one selected sequence-property tool."""
from __future__ import annotations
import argparse,csv,json,os,platform,subprocess,sys,time
from datetime import datetime
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/"src"))

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument("--tool",choices=("netsolp","nanomelt","tnp"),required=True);p.add_argument("--plan-dir",type=Path,required=True);p.add_argument("--output-dir",type=Path,required=True);p.add_argument("--tool-root",type=Path);p.add_argument("--executable",type=Path);p.add_argument("--immune-builder-refine",type=Path);p.add_argument("--generated-at");a=p.parse_args();started=time.perf_counter();plan=a.plan_dir.resolve(strict=True);out=a.output_dir.absolute()
 samples=_csv(plan/"double_mutant_score_samples.csv");gate=_json(plan/"double_mutant_plan_gate.json")
 if len(samples)!=87 or len({row["sample_uid"] for row in samples})!=87 or gate.get("release")!="ready_for_unfiltered_double_mutant_scoring":raise ValueError("Double-mutant plan is incomplete")
 if out.exists():raise FileExistsError(f"Refusing to overwrite: {out}")
 out.mkdir(parents=True);generated=a.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
 if a.tool=="netsolp":rows,command=_netsolp(a,plan,samples,out)
 elif a.tool=="nanomelt":rows,command=_nanomelt(a,plan,samples,out)
 else:rows,command=_tnp(a,samples,out)
 passed=sum(row["scoring_status"]=="pass" for row in rows);status="pass" if passed==87 else "failed";_write_csv(out/f"{a.tool}_sample_scores.csv",rows);_write_json(out/f"{a.tool}_model_run.json",{"schema_version":1,"status":status,"generated_at":generated,"tool":a.tool,"python":platform.python_version(),"sample_count":87,"pass_count":passed,"failure_count":87-passed,"elapsed_seconds":round(time.perf_counter()-started,6),"command":command,"candidate_filtering_applied":False})
 if status!="pass":raise RuntimeError(f"{a.tool} coverage gate failed: {passed}/87")
 return 0

def _netsolp(a,plan,samples,out):
 from antibody_optimization.netsolp_yield import normalize_netsolp_scores
 if not a.tool_root:raise ValueError("NetSolP requires --tool-root")
 root=a.tool_root.resolve(strict=True);predict=root/"predict.py";raw=out/"netsolp_raw_predictions.csv";command=[sys.executable,str(predict),"--FASTA_PATH",str(plan/"double_mutant_sequences.fasta"),"--OUTPUT_PATH",str(raw),"--MODEL_TYPE","Distilled","--PREDICTION_TYPE","SU","--NUM_THREADS","12"];subprocess.run(command,cwd=root,check=True);normalized=[{"sample_uid":r["sample_uid"],"sequence_raw":r["sequence_raw"]} for r in samples];return normalize_netsolp_scores(normalized,_csv(raw),expected_count=87),command
def _nanomelt(a,plan,samples,out):
 from antibody_optimization.nanomelt_yield import normalize_nanomelt_scores
 if not a.executable:raise ValueError("NanoMelt requires --executable")
 import torch
 if not torch.cuda.is_available():raise RuntimeError("NanoMelt requires a visible CUDA GPU")
 raw=out/"nanomelt_raw_predictions.csv";command=[str(a.executable.resolve(strict=True)),"predict","-i",str(plan/"double_mutant_sequences.fasta"),"-o",str(raw),"-align","-ncpu","1","-v"];subprocess.run(command,check=True);normalized=[{"sample_uid":r["sample_uid"],"sequence_raw":r["sequence_raw"]} for r in samples];return normalize_nanomelt_scores(normalized,_csv(raw),expected_pass_count=87,expected_plan_count=87),command
def _tnp(a,samples,out):
 from antibody_optimization.tnp_yield import failed_tnp_result,normalize_tnp_result,verify_immune_builder_refine_patch
 if not a.tool_root or not a.executable or not a.immune_builder_refine:raise ValueError("TNP requires --tool-root, --executable and --immune-builder-refine")
 source=a.tool_root.resolve(strict=True);executable=a.executable.resolve(strict=True);verify_immune_builder_refine_patch(a.immune_builder_refine.resolve(strict=True).read_text(encoding="utf-8"));raw=out/"raw_tnp";raw.mkdir();env=os.environ.copy();env["PYTHONPATH"]=str(source);rows=[]
 for index,sample in enumerate(samples,1):
  uid=sample["sample_uid"];sample_dir=raw/f"sample_{index:03d}";started=time.perf_counter();command=[str(executable),"--seq",sample["sequence_raw"],"--name",uid,"--output",str(sample_dir),"--hscale","0","--ncores","1"];print(f"[{index}/87] TNP {uid}",flush=True)
  try:
   completed=subprocess.run(command,cwd=source,env=env,text=True,capture_output=True,check=False);(raw/f"sample_{index:03d}.stdout.log").write_text(completed.stdout,encoding="utf-8");(raw/f"sample_{index:03d}.stderr.log").write_text(completed.stderr,encoding="utf-8")
   if completed.returncode:raise RuntimeError(f"TNP exit code {completed.returncode}")
   result=_single(sample_dir.glob("TNP_Results*.json"));details=_single(sample_dir.glob("Final_Models/*_Model_Details.json"));modelled=str(_json(details)["sequences"]["H"]);row=normalize_tnp_result(sample,_json(result),modelled_sequence=modelled,elapsed_seconds=time.perf_counter()-started)
   if row["trimmed_n_terminal"] or row["trimmed_c_terminal"]!="GS" or int(row["modelled_length_aa"])!=126:raise ValueError("TNP modeled an unexpected sequence domain")
   rows.append(row)
  except Exception as error:rows.append(failed_tnp_result(sample,f"{type(error).__name__}: {error}",time.perf_counter()-started))
 return rows,[str(executable),"--seq","<87 planned sequences>"]
def _single(paths):
 found=list(paths)
 if len(found)!=1:raise ValueError(f"Expected one result; found {len(found)}")
 return found[0]
def _csv(path):
 with path.open("r",encoding="utf-8-sig",newline="") as h:return list(csv.DictReader(h))
def _json(path):return json.loads(path.read_text(encoding="utf-8-sig"))
def _write_csv(path,rows):
 with path.open("w",encoding="utf-8-sig",newline="") as h:w=csv.DictWriter(h,fieldnames=list(rows[0]),lineterminator="\n");w.writeheader();w.writerows(rows)
def _write_json(path,value):path.write_text(json.dumps(value,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8",newline="\n")
if __name__=="__main__":raise SystemExit(main())
