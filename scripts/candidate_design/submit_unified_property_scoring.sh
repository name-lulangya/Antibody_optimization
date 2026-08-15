#!/usr/bin/env bash
set -euo pipefail
PLAN=docs/result_artifacts/candidate_design/unified_property_scoring_plan_20260815
SCORE_ROOT=results/candidate_design/unified_property_scoring_20260815
RESULT=docs/result_artifacts/candidate_design/unified_property_scoring_result_20260815
SUMMARY=docs/run_summaries/candidate_design/unified_property_scoring_result_20260815.json
LOG_DIR=logs/unified_property_scoring
test -s "$PLAN/unified_property_scoring_plan_gate.json"
if [[ -e "$SCORE_ROOT/netsolp_scores" || -e "$SCORE_ROOT/nanomelt_scores" || -e "$RESULT" || -e "$SUMMARY" ]]; then
  echo "Unified property scoring output already exists; refusing to overwrite." >&2
  exit 1
fi
mkdir -p "$LOG_DIR"
NET_RAW=$(sbatch --parsable scripts/candidate_design/submit_unified_property_netsolp.slurm)
MELT_RAW=$(sbatch --parsable scripts/candidate_design/submit_unified_property_nanomelt.slurm)
NET_JOB=${NET_RAW%%;*}
MELT_JOB=${MELT_RAW%%;*}
MERGE_RAW=$(sbatch --parsable --dependency="afterok:${NET_JOB}:${MELT_JOB}" scripts/candidate_design/submit_unified_property_analysis.slurm)
MERGE_JOB=${MERGE_RAW%%;*}
echo "NetSolP job: ${NET_JOB}"
echo "NanoMelt job: ${MELT_JOB}"
echo "Analysis job: ${MERGE_JOB} (afterok:${NET_JOB}:${MELT_JOB})"
