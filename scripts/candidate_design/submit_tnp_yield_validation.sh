#!/usr/bin/env bash
set -euo pipefail

LOG_DIR="logs/tnp_yield_validation"
PLAN_DIR="docs/result_artifacts/candidate_design/tnp_yield_validation_plan_20260814"
SCORE_DIR="results/candidate_design/tnp_yield_validation_20260814/model_scores"
ARTIFACT_DIR="docs/result_artifacts/candidate_design/tnp_yield_validation_result_20260814"

[[ -d "${PLAN_DIR}" ]] || { echo "Missing validation plan: ${PLAN_DIR}" >&2; exit 2; }
[[ ! -e "${SCORE_DIR}" ]] || { echo "Score directory already exists: ${SCORE_DIR}" >&2; exit 2; }
[[ ! -e "${ARTIFACT_DIR}" ]] || { echo "Artifact directory already exists: ${ARTIFACT_DIR}" >&2; exit 2; }
mkdir -p "${LOG_DIR}" "$(dirname "${SCORE_DIR}")"
JOB_ID=$(sbatch --parsable scripts/candidate_design/submit_tnp_yield_validation.slurm)
echo "Submitted single-process TNP-yield validation job ${JOB_ID%%;*}"
