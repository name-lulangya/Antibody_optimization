#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/../.."

PLAN_GATE="docs/result_artifacts/candidate_design/expression_double_mutant_plan_20260822/expression_double_mutant_plan_gate.json"
SCORE_ROOT="results/candidate_design/expression_double_mutant_scan_20260822"
FINAL_DIR="docs/result_artifacts/candidate_design/expression_double_mutant_property_matrix_20260822"

test -s "${PLAN_GATE}"
if [[ -e "${SCORE_ROOT}" || -e "${FINAL_DIR}" ]]; then
  echo "Output already exists; this under-five-hour workflow must be rerun from a clean new output path." >&2
  exit 1
fi

mkdir -p logs/expression_double_mutant_scan
sbatch scripts/candidate_design/submit_expression_double_mutant_scan.slurm
