#!/usr/bin/env bash
set -euo pipefail

LOG_DIR="logs/unified_tnp_review"
PLAN_DIR="docs/result_artifacts/candidate_design/unified_tnp_review_plan_20260815"
SCORE_DIR="results/candidate_design/unified_tnp_review_20260815/model_scores"
ARTIFACT_DIR="docs/result_artifacts/candidate_design/unified_tnp_review_result_20260815"
SUMMARY="docs/run_summaries/candidate_design/unified_tnp_review_result_20260815.json"

[[ -d "${PLAN_DIR}" ]] || { echo "Missing unified TNP plan: ${PLAN_DIR}" >&2; exit 2; }
[[ ! -e "${SCORE_DIR}" ]] || { echo "Score directory already exists: ${SCORE_DIR}" >&2; exit 2; }
[[ ! -e "${ARTIFACT_DIR}" ]] || { echo "Artifact directory already exists: ${ARTIFACT_DIR}" >&2; exit 2; }
[[ ! -e "${SUMMARY}" ]] || { echo "Run summary already exists: ${SUMMARY}" >&2; exit 2; }
mkdir -p "${LOG_DIR}" "$(dirname "${SCORE_DIR}")" "$(dirname "${SUMMARY}")"
JOB_ID=$(sbatch --parsable scripts/candidate_design/submit_unified_tnp_review.slurm)
echo "Submitted single-process unified TNP review job ${JOB_ID%%;*}"
