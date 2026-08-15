#!/usr/bin/env python3
"""Run one official NetSolP batch for WT plus 1,962 released single mutants."""

from __future__ import annotations
import argparse,csv,json,platform,subprocess,sys,tempfile,time
from datetime import datetime
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/"src"))
from antibody_optimization.file_transaction import replace_staged_files,validate_file_paths
from antibody_optimization.netsolp_yield import normalize_netsolp_scores

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument("--plan-dir",type=Path,required=True);p.add_argument("--netsolp-workdir",type=Path,required=True);p.add_argument("--output-dir",type=Path,required=True);p.add_argument("--num-threads",type=int,default=12);p.add_argument("--generated-at");a=p.parse_args();started=time.perf_counter();generated=a.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
 plan=a.plan_dir.resolve(strict=True);work=a.netsolp_workdir.resolve(strict=True);samples_path=plan/"unified_property_samples.csv";fasta=plan/"unified_property_sequences.fasta";contract_path=plan/"unified_property_scoring_contract.json";predict=work/"predict.py";samples=_csv(samples_path);contract=_json(contract_path)
 if len(samples)!=1963 or contract.get("status")!="pass":raise ValueError("Unified property plan is incomplete")
 spec=contract["netsolp"]
 if spec["model_type"]!="Distilled" or spec["prediction_type"]!="SU":raise ValueError("Unexpected NetSolP contract")
 output=a.output_dir.absolute();targets=[output/"netsolp_raw_predictions.csv",output/"netsolp_sample_scores.csv",output/"netsolp_model_run.json"]
 valid=validate_file_paths(project_root=ROOT,source_paths=[samples_path,fasta,contract_path,predict],target_paths=targets)
 if output.exists():raise FileExistsError("Refusing to overwrite unified NetSolP scores")
 output.parent.mkdir(parents=True,exist_ok=True)
 normalized_samples=[{"sample_uid":row["score_id"],"sequence_raw":row["sequence_raw"]} for row in samples]
 with tempfile.TemporaryDirectory(prefix=".property-netsolp-",dir=output.parent) as tmp:
  stage=Path(tmp);raw=stage/targets[0].name;command=[sys.executable,str(predict),"--FASTA_PATH",str(fasta),"--OUTPUT_PATH",str(raw),"--MODEL_TYPE","Distilled","--PREDICTION_TYPE","SU","--NUM_THREADS",str(a.num_threads)]
  completed=subprocess.run(command,cwd=work,check=True,text=True,capture_output=True);scores=normalize_netsolp_scores(normalized_samples,_csv(raw),expected_count=1963);score_path=stage/targets[1].name;run=stage/targets[2].name;_write_csv(score_path,scores);_write_json(run,{"schema_version":1,"status":"pass","generated_at":generated,"python":platform.python_version(),"elapsed_seconds":round(time.perf_counter()-started,6),"command":command,"sample_count":len(scores),"pass_count":sum(row["scoring_status"]=="pass" for row in scores),"model_type":"Distilled","prediction_type":"SU","stdout":completed.stdout.strip().splitlines(),"stderr":completed.stderr.strip().splitlines()});output.mkdir(parents=False,exist_ok=False);replace_staged_files({raw:targets[0],score_path:targets[1],run:targets[2]},project_root=ROOT,protected_source_paths=valid.source_paths)
 return 0
def _csv(path):
 with path.open("r",encoding="utf-8-sig",newline="") as h:return list(csv.DictReader(h))
def _json(path):return json.loads(path.read_text(encoding="utf-8-sig"))
def _write_csv(path,rows):
 with path.open("w",encoding="utf-8-sig",newline="") as h:w=csv.DictWriter(h,fieldnames=list(rows[0]),lineterminator="\n");w.writeheader();w.writerows(rows)
def _write_json(path,value):path.write_text(json.dumps(value,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8",newline="\n")
if __name__=="__main__":raise SystemExit(main())
