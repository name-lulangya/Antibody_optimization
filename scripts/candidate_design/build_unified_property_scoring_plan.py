#!/usr/bin/env python3
"""Freeze WT plus 1,962 released Nb252 single mutants for property scoring."""

from __future__ import annotations

import argparse, csv, json, platform, sys, tempfile, time
from collections import Counter
from datetime import datetime
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/"src"))
from antibody_optimization.file_transaction import replace_staged_files,validate_file_paths
from antibody_optimization.unified_property_scoring import build_property_samples

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--unified-antifold-dir",type=Path,required=True);p.add_argument("--output-dir",type=Path,required=True);p.add_argument("--run-summary",type=Path,required=True);p.add_argument("--generated-at");p.add_argument("--check_only",action="store_true");a=p.parse_args()
    started=time.perf_counter();generated=a.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    source=a.unified_antifold_dir.resolve(strict=True);evidence=source/"unified_single_mutant_antifold_evidence.csv";source_gate=source/"unified_single_mutant_antifold_gate.json";review=source/"unified_single_mutant_antifold_scientific_review.json"
    for path in (evidence,source_gate,review):path.resolve(strict=True)
    if _json(review).get("release")!="ready_for_multitool_property_scoring_without_antifold_only_selection":raise ValueError("Unified AntiFold scientific review is not released")
    samples=build_property_samples(_csv(evidence),_json(source_gate));tracks=Counter(row["design_track"] for row in samples[1:])
    if a.check_only:
        print(json.dumps({"status":"pass","score_rows":len(samples),"candidate_rows":len(samples)-1,"track_counts":dict(tracks)},sort_keys=True));return 0
    output=a.output_dir.absolute();summary=a.run_summary.absolute()
    if output.exists() or summary.exists():raise FileExistsError("Refusing to overwrite unified property plan")
    output.parent.mkdir(parents=True,exist_ok=True);summary.parent.mkdir(parents=True,exist_ok=True)
    names={"samples":"unified_property_samples.csv","fasta":"unified_property_sequences.fasta","contract":"unified_property_scoring_contract.json","gate":"unified_property_scoring_plan_gate.json"}
    contract={
      "schema_version":1,"contract_name":"nb252_unified_single_mutant_property_scoring","status":"pass","generated_at":generated,
      "score_row_count":1963,"candidate_count":1962,"wt_score_id":"LTT__Nb252__WT","candidate_selection_performed":False,
      "netsolp":{"environment":"/data/software/env/luly25/netsolp","working_directory":"/homes/Tianlab/luly25/software/netsolp","model_type":"Distilled","prediction_type":"SU","outputs":["predicted_usability","predicted_solubility"],"semantics":"compatibility signals relative to WT; not yield"},
      "nanomelt":{"environment":"/data/software/env/luly25/nanomelt","entry_point":"/data/software/env/luly25/nanomelt/bin/nanomelt","immune_builder_refine":"/data/software/env/luly25/nanomelt/lib/python3.10/site-packages/ImmuneBuilder/refine.py","input_length_aa":128,"expected_scored_length_aa":126,"required_trimmed_c_terminal":"GS","expected_pass_count":1963,"required_openmm_platforms":["Reference","CPU"],"software":{"python":"3.10.20","nanomelt":"1.4.0","torch":"2.7.1+cu126","transformers":"4.56.1","immune_builder":"1.2","openmm":"8.5.2","pdbfixer":"1.12.0","anarci_bioconda":"2024.05.21"},"semantics":"predicted apparent Tm relative to WT; not experimental Tm or yield"},
      "analysis":{"tracks":["affinity_existing_interface_scan","stability_developability_discovery"],"pareto_layers":"track_specific_first_two_fronts","weighted_composite_score":False,"yield_prediction":False,"chemical_risks_are_separate_review_flags":True},
      "runtime_estimate":{"netsolp":"approximately 15-60 minutes","nanomelt":"approximately 5-20 minutes","each_likely_over_one_hour":False,"each_likely_over_five_hours":False,"resume_required":False}
    }
    gate={"schema_version":1,"gate_name":"nb252_unified_property_scoring_plan","status":"pass","generated_at":generated,"score_row_count":len(samples),"candidate_count":len(samples)-1,"track_counts":dict(sorted(tracks.items())),"wt_control_count":1,"release":"ready_for_parallel_netsolp_and_nanomelt_scoring"}
    with tempfile.TemporaryDirectory(prefix=".property-plan-",dir=ROOT) as tmp:
        stage=Path(tmp);staged={k:stage/v for k,v in names.items()};_write_csv(staged["samples"],samples)
        with staged["fasta"].open("w",encoding="utf-8",newline="\n") as h:
            for row in samples:h.write(f">{row['score_id']}\n{row['sequence_raw']}\n")
        _write_json(staged["contract"],contract);_write_json(staged["gate"],gate);rs=stage/"run_summary.json";_write_json(rs,{"schema_version":1,"status":"pass","generated_at":generated,"elapsed_seconds":round(time.perf_counter()-started,6),"python":platform.python_version(),"counts":{"score_rows":1963,"candidate_rows":1962,**dict(tracks)},"release":gate["release"],"outputs":{k:str(output/v) for k,v in names.items()}})
        pairs={staged[k]:output/v for k,v in names.items()};pairs[rs]=summary;validate_file_paths(project_root=ROOT,source_paths=[evidence,source_gate,review],target_paths=pairs.values())
        for path in pairs.values():path.parent.mkdir(parents=True,exist_ok=True)
        replace_staged_files(pairs,project_root=ROOT,protected_source_paths=[evidence,source_gate,review])
    return 0

def _csv(path):
    with path.open("r",encoding="utf-8-sig",newline="") as h:return list(csv.DictReader(h))
def _json(path):return json.loads(path.read_text(encoding="utf-8-sig"))
def _write_csv(path,rows):
    with path.open("w",encoding="utf-8-sig",newline="") as h:w=csv.DictWriter(h,fieldnames=list(rows[0]),lineterminator="\n");w.writeheader();w.writerows(rows)
def _write_json(path,value):path.write_text(json.dumps(value,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8",newline="\n")
if __name__=="__main__":raise SystemExit(main())
