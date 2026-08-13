#!/usr/bin/env python3
"""Analyze nanoBERT scores against reported yield and simple sequence baselines."""

from __future__ import annotations

import argparse, csv, json, platform, sys, tempfile, time
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from antibody_optimization.file_transaction import replace_staged_files, validate_file_paths  # noqa: E402
from antibody_optimization.nanobert_yield import analyze_associations  # noqa: E402
from antibody_optimization.nanobert_yield_plot import render_nanobert_yield_figure  # noqa: E402

NAMES={"samples":"nanobert_yield_sample_evidence.csv","metrics":"nanobert_yield_associations.csv","gate":"nanobert_yield_validation_gate.json","png":"nanobert_yield_validation.png","svg":"nanobert_yield_validation.svg"}


def main()->int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir",type=Path,required=True);parser.add_argument("--score-dir",type=Path,required=True)
    parser.add_argument("--output-dir",type=Path,required=True);parser.add_argument("--run-summary",type=Path,required=True);parser.add_argument("--generated-at")
    args=parser.parse_args();started=time.perf_counter();generated=args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    plan=args.plan_dir.resolve(strict=True);scores=args.score_dir.resolve(strict=True)
    sources=[plan/"nanobert_validation_samples.csv",plan/"nanobert_yield_validation_contract.json",scores/"nanobert_sample_scores.csv",scores/"nanobert_model_run.json"]
    contract=_json(sources[1]);model_run=_json(sources[3])
    if contract.get("status")!="pass" or model_run.get("status")!="pass" or model_run.get("model_revision")!=contract.get("model_revision"):raise ValueError("Plan/model run identity or status mismatch")
    result=analyze_associations(_csv(sources[0]),_csv(sources[2]));targets=[args.output_dir.absolute()/name for name in NAMES.values()]+[args.run_summary.absolute()]
    valid=validate_file_paths(project_root=ROOT,source_paths=sources,target_paths=targets)
    if any(path.exists() for path in valid.target_paths):raise FileExistsError("Refusing to overwrite nanoBERT validation outputs")
    for path in valid.target_paths:path.parent.mkdir(parents=True,exist_ok=True)
    final=dict(zip((*NAMES,"summary"),valid.target_paths,strict=True));primary=result["primary"]
    gate={"schema_version":1,"gate_name":"nb252_nanobert_reported_yield_validation","status":"pass","generated_at":generated,"sample_count":47,"numeric_individual_count":31,"llj_ordinal_censored_count":16,"primary_feature":"nanobert_mean_pll_raw","evidence_level":result["evidence_level"],"decision_reasons":result["decision_reasons"],"primary_statistics":primary,"high_capacity_model_trained":False,"nb252_expression_prediction_validated":False,"release":("ready_for_weak_nanobert_ranking_use" if result["evidence_level"]=="weak_ranking_evidence" else "nanobert_compatibility_filter_only" if result["evidence_level"]=="compatibility_filter_only" else "nanobert_not_supported_for_candidate_use"),"interpretation":"Association with collaborator-reported yield; not an expression-rate or mg/L predictor."}
    with tempfile.TemporaryDirectory(prefix=".nanobert-analysis-",dir=ROOT) as temp:
        stage=Path(temp);staged={k:stage/Path(v).name for k,v in final.items()}
        _write_csv(staged["samples"],result["sample_rows"]);_write_csv(staged["metrics"],result["metric_rows"]);_write_json(staged["gate"],gate)
        render_nanobert_yield_figure(result["sample_rows"],result["metric_rows"],png_path=staged["png"],svg_path=staged["svg"])
        _write_json(staged["summary"],{"schema_version":1,"status":"pass","generated_at":generated,"python":platform.python_version(),"elapsed_seconds":round(time.perf_counter()-started,6),"evidence_level":result["evidence_level"],"counts":{"samples":47,"numeric":31,"llj_ordinal_censored":16,"features":len(result["metric_rows"])},"outputs":{k:str(v) for k,v in final.items() if k!="summary"}})
        replace_staged_files({staged[k]:final[k] for k in staged},project_root=ROOT,protected_source_paths=valid.source_paths)
    return 0


def _csv(path):
    with path.open("r",encoding="utf-8-sig",newline="") as handle:return list(csv.DictReader(handle))
def _json(path):return json.loads(path.read_text(encoding="utf-8-sig"))
def _write_csv(path,rows):
    with path.open("w",encoding="utf-8-sig",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]),lineterminator="\n");writer.writeheader();writer.writerows(rows)
def _write_json(path,value):path.write_text(json.dumps(value,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8",newline="\n")
if __name__=="__main__":raise SystemExit(main())
