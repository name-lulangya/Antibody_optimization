"""Magnitude-aware joint review of complete Nb252 double-mutant scans.

The review joins four model outputs only after all 86 candidates have been
scored. It preserves raw WT-relative values and assigns evidence classes; it
does not select the final 30 experimental sequences or interpret predictions
as measured affinity, expression, solubility, stability, Tm, or yield.
"""
from __future__ import annotations
from collections import Counter
from typing import Mapping,Sequence
from .unified_tnp_review import MAGNITUDE_THRESHOLDS,magnitude_label

EXPECTED=86;WT_ID="LTT__Nb252__WT"
FLAGS=("tnp_flag_total_cdr_length","tnp_flag_cdr3_length","tnp_flag_cdr3_compactness","tnp_flag_psh","tnp_flag_ppc","tnp_flag_pnc")
FLAG_RANK={"green":0,"amber":1,"red":2}
class DoubleMutantAnalysisError(ValueError):pass

def build_joint_evidence(candidates:Sequence[Mapping[str,object]],netsolp:Sequence[Mapping[str,object]],nanomelt:Sequence[Mapping[str,object]],tnp:Sequence[Mapping[str,object]],pyrosetta:Sequence[Mapping[str,object]])->tuple[list[dict[str,object]],dict[str,object]]:
    """Join complete score sets and assign non-final evidence classes."""
    cand=_unique(candidates,"candidate_id",EXPECTED,"candidates");net=_unique(netsolp,"sample_uid",EXPECTED+1,"NetSolP");melt=_unique(nanomelt,"sample_uid",EXPECTED+1,"NanoMelt");tnp_map=_unique(tnp,"sample_uid",EXPECTED+1,"TNP");rose=_unique(pyrosetta,"candidate_id",EXPECTED,"PyRosetta")
    expected={WT_ID,*cand}
    if set(net)!=expected or set(melt)!=expected or set(tnp_map)!=expected:raise DoubleMutantAnalysisError("Sequence-tool identities do not equal WT plus 86 candidates")
    for label,table in (("NetSolP",net),("NanoMelt",melt),("TNP",tnp_map)):
        if any(str(row["scoring_status"])!="pass" for row in table.values()):raise DoubleMutantAnalysisError(f"{label} does not have complete pass coverage")
    wt_net=net[WT_ID];wt_melt=melt[WT_ID];wt_tnp=tnp_map[WT_ID];rows=[]
    for identifier,source in cand.items():
        n=net[identifier];m=melt[identifier];t=tnp_map[identifier];r=rose[identifier]
        if str(n["sequence_raw"])!=str(source["sequence"]) or str(m["sequence_raw"])!=str(source["sequence"]) or str(t["sequence_raw"])!=str(source["sequence"]):raise DoubleMutantAnalysisError(f"Sequence mismatch for {identifier}")
        du=float(n["predicted_usability"])-float(wt_net["predicted_usability"]);ds=float(n["predicted_solubility"])-float(wt_net["predicted_solubility"]);dt=float(m["nanomelt_predicted_apparent_tm_c"])-float(wt_melt["nanomelt_predicted_apparent_tm_c"])
        labels=(magnitude_label(du,MAGNITUDE_THRESHOLDS["netsolp_delta_usability_vs_wt"]),magnitude_label(ds,MAGNITUDE_THRESHOLDS["netsolp_delta_solubility_vs_wt"]),magnitude_label(dt,MAGNITUDE_THRESHOLDS["nanomelt_delta_predicted_apparent_tm_c_vs_wt"]))
        favorable=labels.count("favorable");adverse=labels.count("adverse");regressions=sum(FLAG_RANK[str(t[field]).lower()]>FLAG_RANK[str(wt_tnp[field]).lower()] for field in FLAGS);improvements=sum(FLAG_RANK[str(t[field]).lower()]<FLAG_RANK[str(wt_tnp[field]).lower()] for field in FLAGS)
        dg=float(r["delta_dG_separated_median"]);cross=float(r["delta_cross_interface_energy_median"]);affinity_supported=dg<0 and cross<0
        chemical=bool(str(source.get("new_liability_flags","")).strip());property_nonadverse=adverse==0 and regressions==0 and not chemical
        if affinity_supported and property_nonadverse and favorable:classification="balanced_supported"
        elif affinity_supported and property_nonadverse:classification="affinity_supported_property_nonadverse"
        elif dg<=0 and cross<=0 and favorable and regressions==0 and not chemical:classification="property_supported_affinity_nonadverse"
        else:classification="tradeoff_or_no_clear_joint_support"
        rows.append({**dict(source),"netsolp_predicted_usability":float(n["predicted_usability"]),"netsolp_delta_usability_vs_wt":du,"netsolp_usability_magnitude":labels[0],"netsolp_predicted_solubility":float(n["predicted_solubility"]),"netsolp_delta_solubility_vs_wt":ds,"netsolp_solubility_magnitude":labels[1],"nanomelt_predicted_apparent_tm_c":float(m["nanomelt_predicted_apparent_tm_c"]),"nanomelt_delta_predicted_apparent_tm_c_vs_wt":dt,"nanomelt_tm_magnitude":labels[2],"property_material_favorable_count":favorable,"property_material_adverse_count":adverse,"tnp_psh":float(t["tnp_psh"]),"tnp_psh_delta_vs_wt":float(t["tnp_psh"])-float(wt_tnp["tnp_psh"]),"tnp_flag_regression_count":regressions,"tnp_flag_improvement_count":improvements,"pyrosetta_delta_dG_separated_median":dg,"pyrosetta_delta_cross_interface_energy_median":cross,"pyrosetta_affinity_direction_supported":affinity_supported,"joint_evidence_class":classification,"final_candidate_selection_performed":False})
    counts=dict(Counter(str(row["joint_evidence_class"]) for row in rows));gate={"schema_version":1,"gate_name":"nb252_double_mutant_joint_evidence","status":"pass","release":"ready_for_post_scan_scientific_review","candidate_count":len(rows),"joint_evidence_class_counts":counts,"magnitude_thresholds":{"netsolp_usability":0.01,"netsolp_solubility":0.02,"nanomelt_predicted_tm_c":1.0},"candidate_filtering_applied_during_scoring":False,"final_candidate_selection_performed":False,"interpretation":"Predicted ranking and risk evidence only; final 30-member experimental panel not selected."}
    return sorted(rows,key=lambda row:(str(row["joint_evidence_class"]),str(row["candidate_id"]))),gate
def _unique(rows,key,expected,label):
    result={str(row[key]):row for row in rows}
    if len(rows)!=expected or len(result)!=expected:raise DoubleMutantAnalysisError(f"{label} must contain {expected} unique rows")
    return result
