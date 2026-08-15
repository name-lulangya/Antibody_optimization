#!/usr/bin/env bash
set -euo pipefail

LOG_DIR="logs/tnp_yield_validation_v2"
PLAN_DIR="docs/result_artifacts/candidate_design/tnp_yield_validation_plan_v2_20260814"
SCORE_DIR="results/candidate_design/tnp_yield_validation_v2_20260814/model_scores"
ARTIFACT_DIR="docs/result_artifacts/candidate_design/tnp_yield_validation_result_v2_20260814"

[[ -d "${PLAN_DIR}" ]] || { echo "Missing validation plan: ${PLAN_DIR}" >&2; exit 2; }
[[ ! -e "${ARTIFACT_DIR}" || -d "${ARTIFACT_DIR}" ]] || { echo "Artifact path is not a directory: ${ARTIFACT_DIR}" >&2; exit 2; }
if [[ -e "${SCORE_DIR}" ]]; then
  [[ -d "${SCORE_DIR}" ]] || { echo "Score path is not a directory: ${SCORE_DIR}" >&2; exit 2; }
  [[ -s "${SCORE_DIR}/tnp_sample_scores.csv" ]] || { echo "Incomplete score directory: missing tnp_sample_scores.csv" >&2; exit 2; }
  [[ -s "${SCORE_DIR}/tnp_model_run.json" ]] || { echo "Incomplete score directory: missing tnp_model_run.json" >&2; exit 2; }
  echo "Found completed-score files; the Slurm job will validate and reuse them."
fi
mkdir -p "${LOG_DIR}" "$(dirname "${SCORE_DIR}")"
JOB_ID=$(sbatch --parsable scripts/candidate_design/submit_tnp_yield_validation.slurm)
echo "Submitted single-process TNP-yield validation V2 job ${JOB_ID%%;*}"
