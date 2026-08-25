#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/../.."

PLAN_MANIFEST="docs/result_artifacts/candidate_design/v3_double_mutant_plan_20260825/v3_double_mutant_plan_manifest.json"
SCORE_ROOT="results/candidate_design/v3_double_mutant_scan_20260825"
FINAL_DIR="docs/result_artifacts/candidate_design/v3_double_mutant_property_matrix_20260825"
RUN_SUMMARY="docs/run_summaries/candidate_design/v3_double_mutant_property_matrix_20260825/run_summary.json"

test -s "${PLAN_MANIFEST}"
if [[ -e "${SCORE_ROOT}" || -e "${FINAL_DIR}" || -e "${RUN_SUMMARY}" ]]; then
  echo "A formal output path already exists; use a new explicitly named run after resolving the cause." >&2
  exit 1
fi

mkdir -p logs/v3_double_mutant_scan
sbatch scripts/candidate_design/submit_v3_double_mutant_scan.slurm
