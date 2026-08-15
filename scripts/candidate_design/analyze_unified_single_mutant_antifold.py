#!/usr/bin/env python3
"""Join the complete single-mutant universe to existing AntiFold raw scores."""

from __future__ import annotations

import argparse, csv, json, platform, sys, tempfile, time
from datetime import datetime
from pathlib import Path

PROJECT_ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(PROJECT_ROOT/"src"))
from antibody_optimization.antifold_validation import normalize_antifold_rows
from antibody_optimization.file_transaction import replace_staged_files, validate_file_paths
from antibody_optimization.unified_single_mutant_plot import render_antifold_landscape
from antibody_optimization.unified_single_mutants import evaluate_antifold_landscape

VIEWS=("experimental_vhh_only","experimental_complex_context","af3_vhh_only")

def parse_args():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument("--plan-dir",type=Path,required=True); p.add_argument("--antifold-plan-dir",type=Path,required=True); p.add_argument("--score-dir",type=Path,required=True); p.add_argument("--output-dir",type=Path,required=True); p.add_argument("--run-summary",type=Path,required=True); p.add_argument("--generated-at"); return p.parse_args()

def main():
    a=parse_args(); started=time.perf_counter(); generated=a.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    plan=a.plan_dir.resolve(strict=True); afplan=a.antifold_plan_dir.resolve(strict=True); scores=a.score_dir.resolve(strict=True)
    plan_gate=_json(plan/"unified_single_mutant_plan_gate.json"); model_run=_json(scores/"antifold_model_run.json")
    if plan_gate.get("status")!="pass" or plan_gate.get("release")!="ready_for_existing_antifold_landscape_join": raise ValueError("Unified plan is not released")
    if model_run.get("status")!="pass": raise ValueError("Existing AntiFold model run is not passed")
    views={row["view_id"]:row for row in _csv(afplan/"antifold_structure_views.csv")}
    if tuple(view for view in VIEWS if view in views)!=VIEWS or len(views)!=3: raise ValueError("Expected exactly three AntiFold structure views")
    indexed={view:normalize_antifold_rows(_csv(scores/f"{view}.csv"),view_id=view,vhh_chain=views[view]["vhh_chain"]) for view in VIEWS}
    rows,gate=evaluate_antifold_landscape(_csv(plan/"unified_single_mutant_candidates.csv"),indexed); gate["generated_at"]=generated
    output=a.output_dir.absolute(); summary=a.run_summary.absolute()
    if output.exists() or summary.exists(): raise FileExistsError("Refusing to overwrite unified AntiFold outputs")
    output.parent.mkdir(parents=True,exist_ok=True); summary.parent.mkdir(parents=True,exist_ok=True)
    names={"evidence":"unified_single_mutant_antifold_evidence.csv","gate":"unified_single_mutant_antifold_gate.json","png":"unified_single_mutant_antifold.png","svg":"unified_single_mutant_antifold.svg"}
    with tempfile.TemporaryDirectory(prefix=".unified-antifold-",dir=PROJECT_ROOT) as tmp:
        stage=Path(tmp); staged={k:stage/v for k,v in names.items()}; _write_csv(staged["evidence"],rows); _write_json(staged["gate"],gate); render_antifold_landscape(rows,png_path=staged["png"],svg_path=staged["svg"])
        rs=stage/"run_summary.json"; _write_json(rs,{"schema_version":1,"status":gate["status"],"generated_at":generated,"elapsed_seconds":round(time.perf_counter()-started,6),"python":platform.python_version(),"command_argv":[sys.executable,str(Path(__file__).resolve()),*sys.argv[1:]],"counts":{"candidates":len(rows),**gate["evaluation_scope_counts"]},"release":gate["release"],"model_rerun_performed":False,"source_model_run":str(scores/"antifold_model_run.json"),"outputs":{k:str(output/v) for k,v in names.items()}})
        pairs={staged[k]:output/v for k,v in names.items()};pairs[rs]=summary
        sources=[plan/"unified_single_mutant_candidates.csv",plan/"unified_single_mutant_plan_gate.json",afplan/"antifold_structure_views.csv",scores/"antifold_model_run.json",*[scores/f"{v}.csv" for v in VIEWS]]
        validate_file_paths(project_root=PROJECT_ROOT,source_paths=sources,target_paths=pairs.values())
        for p in pairs.values():p.parent.mkdir(parents=True,exist_ok=True)
        replace_staged_files(pairs,project_root=PROJECT_ROOT,protected_source_paths=sources)
    return 0 if gate["status"]=="pass" else 2

def _csv(path):
    with path.open("r",encoding="utf-8-sig",newline="") as h:return list(csv.DictReader(h))
def _json(path):return json.loads(path.read_text(encoding="utf-8-sig"))
def _write_json(path,value):path.write_text(json.dumps(value,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8",newline="\n")
def _write_csv(path,rows):
    with path.open("w",encoding="utf-8-sig",newline="") as h:w=csv.DictWriter(h,fieldnames=list(rows[0]),lineterminator="\n");w.writeheader();w.writerows(rows)
if __name__=="__main__":raise SystemExit(main())
