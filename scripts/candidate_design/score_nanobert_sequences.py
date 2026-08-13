#!/usr/bin/env python3
"""Score the fixed 47 sequences by single-mask nanoBERT pseudo-log-likelihood."""

from __future__ import annotations

import argparse, csv, importlib.metadata, json, math, platform, time
from pathlib import Path

import torch
from transformers import AutoModelForMaskedLM, AutoTokenizer

REVISION = "edc8182ad89a827f8737fa572c6b5fac6197e6b0"


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, required=True); parser.add_argument("--model-snapshot", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True); parser.add_argument("--mask-batch-size", type=int, default=64); parser.add_argument("--check_only", action="store_true")
    args=parser.parse_args(); plan=args.plan_dir.resolve(strict=True); snapshot=args.model_snapshot.resolve(strict=True)
    if snapshot.name != REVISION: raise ValueError("Model snapshot revision path does not match the pinned revision")
    samples=_csv(plan/"nanobert_validation_samples.csv"); regions=_csv(plan/"nanobert_validation_regions.csv"); contract=_json(plan/"nanobert_yield_validation_contract.json")
    if contract.get("status") != "pass" or len(samples) != 47: raise ValueError("Validation plan is not released")
    tokenizer=AutoTokenizer.from_pretrained(snapshot, local_files_only=True); model=AutoModelForMaskedLM.from_pretrained(snapshot, local_files_only=True)
    if tokenizer.mask_token_id is None: raise ValueError("Tokenizer has no mask token")
    device=torch.device("cuda" if torch.cuda.is_available() else "cpu"); model.to(device); model.eval()
    _validate_tokenizer(tokenizer, str(samples[0]["sequence_raw"]))
    if args.check_only:
        _score_one(str(samples[0]["sequence_raw"]), tokenizer, model, device, min(args.mask_batch_size, 16)); print(json.dumps({"status":"pass","device":str(device),"samples":47})); return 0
    if args.output_dir.exists(): raise FileExistsError("Refusing to overwrite nanoBERT score directory")
    args.output_dir.mkdir(parents=True); started=time.perf_counter(); position_output=[]; sample_output=[]
    region_map={(str(row["sample_uid"]),int(row["sequence_index_1based"])):str(row["region"]) for row in regions}
    for index, sample in enumerate(samples, start=1):
        sequence=str(sample["sequence_raw"]); scores=_score_one(sequence, tokenizer, model, device, args.mask_batch_size)
        for position,(aa,logp,prob) in enumerate(zip(sequence,scores[0],scores[1]), start=1):
            position_output.append({"sample_uid":sample["sample_uid"],"sequence_index_1based":position,"residue_aa":aa,"region":region_map.get((str(sample["sample_uid"]),position),"unassigned"),"wt_probability":prob,"wt_log_probability":logp})
        numbered=[row for row in position_output if row["sample_uid"]==sample["sample_uid"] and row["region"]!="unassigned"]
        fr=[row for row in numbered if str(row["region"]).startswith("FR")]; cdr=[row for row in numbered if str(row["region"]).startswith("CDR")]
        sample_output.append({"sample_uid":sample["sample_uid"],"scoring_status":"pass","sequence_length_aa":len(sequence),"nanobert_sum_pll_raw":sum(scores[0]),"nanobert_mean_pll_raw":sum(scores[0])/len(sequence),"nanobert_mean_pll_numbered":_mean(numbered),"nanobert_mean_pll_fr":_mean(fr),"nanobert_mean_pll_cdr":_mean(cdr)})
        print(f"Scored {index}/47 {sample['sample_uid']}", flush=True)
    _write_csv(args.output_dir/"nanobert_sample_scores.csv",sample_output); _write_csv(args.output_dir/"nanobert_position_scores.csv",position_output)
    _write_json(args.output_dir/"nanobert_model_run.json",{"schema_version":1,"status":"pass","model_id":"NaturalAntibody/nanoBERT","model_revision":REVISION,"model_snapshot":str(snapshot),"device":str(device),"torch":torch.__version__,"transformers":importlib.metadata.version("transformers"),"python":platform.python_version(),"mask_batch_size":args.mask_batch_size,"sample_count":47,"position_count":len(position_output),"elapsed_seconds":round(time.perf_counter()-started,6),"score_definition":"single-mask WT natural-log probability"})
    return 0


def _validate_tokenizer(tokenizer, sequence):
    ids=tokenizer(sequence, add_special_tokens=True)["input_ids"]
    if len(ids)!=len(sequence)+2: raise ValueError("nanoBERT tokenizer is not one token per residue")
    decoded=[tokenizer.convert_ids_to_tokens(value) for value in ids[1:-1]]
    if decoded!=list(sequence): raise ValueError("nanoBERT tokenizer does not preserve residue identities")
def _score_one(sequence,tokenizer,model,device,batch_size):
    encoded=tokenizer(sequence,add_special_tokens=True,return_tensors="pt"); ids=encoded["input_ids"][0]; logs=[]; probs=[]
    for start in range(0,len(sequence),batch_size):
        positions=list(range(start+1,min(start+batch_size,len(sequence))+1)); batch=ids.repeat(len(positions),1)
        for row,pos in enumerate(positions): batch[row,pos]=tokenizer.mask_token_id
        with torch.inference_mode(): logits=model(input_ids=batch.to(device),attention_mask=torch.ones_like(batch,device=device)).logits
        log_probs=torch.log_softmax(logits,dim=-1)
        for row,pos in enumerate(positions):
            value=float(log_probs[row,pos,int(ids[pos])].cpu()); logs.append(value); probs.append(math.exp(value))
    return logs,probs
def _mean(rows): return sum(float(row["wt_log_probability"]) for row in rows)/len(rows) if rows else ""
def _csv(path):
    with path.open("r",encoding="utf-8-sig",newline="") as handle:return list(csv.DictReader(handle))
def _json(path):return json.loads(path.read_text(encoding="utf-8-sig"))
def _write_csv(path,rows):
    with path.open("w",encoding="utf-8-sig",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]),lineterminator="\n");writer.writeheader();writer.writerows(rows)
def _write_json(path,value):path.write_text(json.dumps(value,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8",newline="\n")
if __name__=="__main__":raise SystemExit(main())
