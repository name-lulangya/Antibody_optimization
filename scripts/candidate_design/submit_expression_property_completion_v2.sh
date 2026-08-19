#!/usr/bin/env bash
set -euo pipefail

PLAN_DIR="docs/result_artifacts/candidate_design/expression_property_completion_plan_v2_20260819"
SCORE_ROOT="results/candidate_design/expression_property_completion_v2_20260819"
VALIDATION_RESULT="docs/result_artifacts/candidate_design/expression_property_reuse_validation_v2_20260819"
FINAL_RESULT="docs/result_artifacts/candidate_design/expression_property_complete_matrix_v2_20260819"
LOG_DIR="logs/expression_property_completion_v2"

[[ -s "${PLAN_DIR}/expression_property_completion_plan_gate.json" ]] || {
  echo "Missing completion plan gate: ${PLAN_DIR}" >&2
  exit 2
}
[[ ! -e "${SCORE_ROOT}" ]] || { echo "Score root already exists: ${SCORE_ROOT}" >&2; exit 2; }
[[ ! -e "${VALIDATION_RESULT}" ]] || { echo "Validation result already exists: ${VALIDATION_RESULT}" >&2; exit 2; }
[[ ! -e "${FINAL_RESULT}" ]] || { echo "Final result already exists: ${FINAL_RESULT}" >&2; exit 2; }

mkdir -p "${LOG_DIR}" "$(dirname "${SCORE_ROOT}")"
JOB_ID=$(sbatch --parsable scripts/candidate_design/submit_expression_property_completion_v2.slurm)
echo "Submitted expression-property completion job ${JOB_ID%%;*}"
