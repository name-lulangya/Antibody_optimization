#!/usr/bin/env python3
"""Run one NanoMelt batch for WT plus 1,962 released single mutants."""

from __future__ import annotations
import argparse,csv,importlib.metadata,json,platform,subprocess,sys,tempfile,time
from datetime import datetime
from pathlib import Path
import anarci as anarci_module
import torch
from openmm import Platform
from packaging.version import Version
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/"src"))
from antibody_optimization.file_transaction import replace_staged_files,validate_file_paths
from antibody_optimization.nanomelt_yield import normalize_nanomelt_scores,verify_anarci_runtime,verify_required_openmm_platforms

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument("--plan-dir",type=Path,required=True);p.add_argument("--nanomelt-executable",type=Path,required=True);p.add_argument("--immune-builder-refine",type=Path,required=True);p.add_argument("--output-dir",type=Path,required=True);p.add_argument("--generated-at");a=p.parse_args();started=time.perf_counter();generated=a.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
 plan=a.plan_dir.resolve(strict=True);executable=a.nanomelt_executable.resolve(strict=True);refine=a.immune_builder_refine.resolve(strict=True);samples_path=plan/"unified_property_samples.csv";fasta=plan/"unified_property_sequences.fasta";contract_path=plan/"unified_property_scoring_contract.json";samples=_csv(samples_path);contract=_json(contract_path);spec=contract["nanomelt"]
 if len(samples)!=1963 or Path(sys.prefix)!=Path(spec["environment"]):raise ValueError("NanoMelt plan/environment mismatch")
 for key,dist in (("nanomelt","nanomelt"),("torch","torch"),("transformers","transformers"),("immune_builder","ImmuneBuilder"),("openmm","OpenMM"),("pdbfixer","pdbfixer")):
  actual=torch.__version__ if key=="torch" else importlib.metadata.version(dist)
  if Version(actual)!=Version(str(spec["software"][key]).split("_")[0]):raise ValueError(f"NanoMelt environment mismatch for {key}: {actual}")
 if platform.python_version()!=spec["software"]["python"]:raise ValueError("NanoMelt Python mismatch")
 anarci_runtime=verify_anarci_runtime(anarci_module,Path(sys.prefix),expected_conda_version=spec["software"]["anarci_bioconda"])
 text=refine.read_text(encoding="utf-8")
 if "platform, {'Threads', str(n_threads)})" in text or "platform, {'Threads': str(n_threads)})" not in text:raise ValueError("ImmuneBuilder Threads patch is missing")
 if not torch.cuda.is_available():raise RuntimeError("NanoMelt requires a visible CUDA GPU")
 openmm_platforms=verify_required_openmm_platforms([Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())],spec["required_openmm_platforms"])
 output=a.output_dir.absolute();targets=[output/"nanomelt_raw_predictions.csv",output/"nanomelt_sample_scores.csv",output/"nanomelt_model_run.json"];valid=validate_file_paths(project_root=ROOT,source_paths=[samples_path,fasta,contract_path,executable,refine],target_paths=targets)
 if output.exists():raise FileExistsError("Refusing to overwrite unified NanoMelt scores")
 output.parent.mkdir(parents=True,exist_ok=True);normalized_samples=[{"sample_uid":row["score_id"],"sequence_raw":row["sequence_raw"]} for row in samples]
 with tempfile.TemporaryDirectory(prefix=".property-nanomelt-",dir=output.parent) as tmp:
  stage=Path(tmp);raw=stage/targets[0].name;command=[str(executable),"predict","-i",str(fasta),"-o",str(raw),"-align","-ncpu","1","-v"];subprocess.run(command,check=True,text=True);scores=normalize_nanomelt_scores(normalized_samples,_csv(raw),expected_pass_count=1963,expected_plan_count=1963);score_path=stage/targets[1].name;run=stage/targets[2].name;_write_csv(score_path,scores);_write_json(run,{"schema_version":1,"status":"pass","generated_at":generated,"python":platform.python_version(),"elapsed_seconds":round(time.perf_counter()-started,6),"command":command,"sample_count":len(scores),"pass_count":sum(row["scoring_status"]=="pass" for row in scores),"cuda_device":torch.cuda.get_device_name(0),"torch":torch.__version__,"torch_cuda":torch.version.cuda,"anarci_runtime":anarci_runtime,"openmm_platforms":openmm_platforms});output.mkdir(parents=False,exist_ok=False);replace_staged_files({raw:targets[0],score_path:targets[1],run:targets[2]},project_root=ROOT,protected_source_paths=valid.source_paths)
 return 0
def _csv(path):
 with path.open("r",encoding="utf-8-sig",newline="") as h:return list(csv.DictReader(h))
def _json(path):return json.loads(path.read_text(encoding="utf-8-sig"))
def _write_csv(path,rows):
 with path.open("w",encoding="utf-8-sig",newline="") as h:w=csv.DictWriter(h,fieldnames=list(rows[0]),lineterminator="\n");w.writeheader();w.writerows(rows)
def _write_json(path,value):path.write_text(json.dumps(value,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8",newline="\n")
if __name__=="__main__":raise SystemExit(main())
