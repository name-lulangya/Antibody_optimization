#!/usr/bin/env python3
"""Join complete double-mutant scans and render decision-facing evidence."""
from __future__ import annotations
import argparse,csv,json,platform,sys,time
from datetime import datetime
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/"src"))
from antibody_optimization.double_mutant_analysis import build_joint_evidence
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument("--plan-dir",type=Path,required=True);p.add_argument("--score-root",type=Path,required=True);p.add_argument("--output-dir",type=Path,required=True);p.add_argument("--run-summary",type=Path,required=True);p.add_argument("--generated-at");a=p.parse_args();started=time.perf_counter();plan=a.plan_dir.resolve(strict=True);scores=a.score_root.resolve(strict=True);out=a.output_dir.absolute();summary=a.run_summary.absolute()
 if out.exists() or summary.exists():raise FileExistsError(f"Refusing to overwrite: {out} or {summary}")
 rows,gate=build_joint_evidence(_csv(plan/"double_mutant_candidates.csv"),_csv(scores/"netsolp/netsolp_sample_scores.csv"),_csv(scores/"nanomelt/nanomelt_sample_scores.csv"),_csv(scores/"tnp/tnp_sample_scores.csv"),_csv(scores/"pyrosetta/double_mutant_candidate_summary.csv"));out.mkdir(parents=True);_write_csv(out/"double_mutant_joint_evidence.csv",rows);_write_json(out/"double_mutant_joint_evidence_gate.json",gate);_plot(rows,out/"double_mutant_joint_evidence.png",out/"double_mutant_joint_evidence.svg");summary.parent.mkdir(parents=True,exist_ok=True);_write_json(summary,{"schema_version":1,"status":"pass","generated_at":a.generated_at or datetime.now().astimezone().isoformat(timespec="seconds"),"python":platform.python_version(),"elapsed_seconds":round(time.perf_counter()-started,6),"candidate_count":86,"release":gate["release"],"output_dir":str(out),"final_candidate_selection_performed":False});return 0
def _plot(rows,png,svg):
 import matplotlib.pyplot as plt
 colors={"balanced_supported":"#228833","affinity_supported_property_nonadverse":"#4477AA","property_supported_affinity_nonadverse":"#CCBB44","tradeoff_or_no_clear_joint_support":"#BBBBBB"};fig,axes=plt.subplots(1,3,figsize=(13.2,4.1))
 for label,color in colors.items():
  selected=[r for r in rows if r["joint_evidence_class"]==label];axes[0].scatter([r["pyrosetta_delta_dG_separated_median"] for r in selected],[r["pyrosetta_delta_cross_interface_energy_median"] for r in selected],s=18,label=label.replace("_"," "),color=color,alpha=.8);axes[1].scatter([r["netsolp_delta_usability_vs_wt"] for r in selected],[r["nanomelt_delta_predicted_apparent_tm_c_vs_wt"] for r in selected],s=18,color=color,alpha=.8)
 axes[0].axhline(0,color="#555555",lw=.7);axes[0].axvline(0,color="#555555",lw=.7);axes[0].set(xlabel="PyRosetta ΔdG separated (REU)",ylabel="PyRosetta Δcross-interface energy (REU)");axes[1].axhline(0,color="#555555",lw=.7);axes[1].axvline(0,color="#555555",lw=.7);axes[1].set(xlabel="NetSolP ΔU vs WT",ylabel="NanoMelt Δpredicted Tm vs WT (°C)")
 counts=[sum(r["joint_evidence_class"]==label for r in rows) for label in colors];axes[2].bar(range(4),counts,color=list(colors.values()));axes[2].set_xticks(range(4),["Balanced","Affinity +\nproperty nonadverse","Property +\naffinity nonadverse","Tradeoff /\nunclear"],rotation=12);axes[2].set_ylabel("Candidate count")
 for ax in axes:ax.spines[["top","right"]].set_visible(False)
 handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc="upper center",ncol=2,frameon=False);fig.text(.5,.01,"Predicted ranking/risk evidence only; no final experimental selection.",ha="center",fontsize=8);fig.tight_layout(rect=(0,.05,1,.87));fig.savefig(png,dpi=600);fig.savefig(svg);plt.close(fig)
def _csv(path):
 with path.open("r",encoding="utf-8-sig",newline="") as h:return list(csv.DictReader(h))
def _write_csv(path,rows):
 with path.open("w",encoding="utf-8-sig",newline="") as h:w=csv.DictWriter(h,fieldnames=list(rows[0]),lineterminator="\n");w.writeheader();w.writerows(rows)
def _write_json(path,value):path.write_text(json.dumps(value,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8",newline="\n")
if __name__=="__main__":raise SystemExit(main())
