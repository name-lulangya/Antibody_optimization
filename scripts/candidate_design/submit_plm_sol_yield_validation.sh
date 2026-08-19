#!/usr/bin/env bash
set -euo pipefail

LOG_DIR="logs/plm_sol_yield_validation"
PLAN_DIR="docs/result_artifacts/candidate_design/plm_sol_yield_validation_plan_20260819"
RUN_ROOT="results/candidate_design/plm_sol_yield_validation_20260819"
ARTIFACT_DIR="docs/result_artifacts/candidate_design/plm_sol_yield_validation_result_20260819"

[[ -d "${PLAN_DIR}" ]] || { echo "Missing validation plan: ${PLAN_DIR}" >&2; exit 2; }
[[ ! -e "${RUN_ROOT}" ]] || { echo "PLM_Sol run directory already exists: ${RUN_ROOT}" >&2; exit 2; }
[[ ! -e "${ARTIFACT_DIR}" ]] || { echo "Artifact directory already exists: ${ARTIFACT_DIR}" >&2; exit 2; }
mkdir -p "${LOG_DIR}" "$(dirname "${RUN_ROOT}")"
JOB_ID=$(sbatch --parsable scripts/candidate_design/submit_plm_sol_yield_validation.slurm)
echo "Submitted PLM_Sol-yield validation job ${JOB_ID%%;*}"
