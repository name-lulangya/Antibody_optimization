#!/usr/bin/env bash
set -euo pipefail

PLAN_DIR=docs/result_artifacts/candidate_design/unified_single_mutant_plan_20260815
SCORE_DIR=results/candidate_design/antifold_validation_20260815/model_scores
OUTPUT_DIR=docs/result_artifacts/candidate_design/unified_single_mutant_antifold_20260815
SUMMARY=docs/run_summaries/candidate_design/unified_single_mutant_antifold_20260815.json
LOG_DIR=logs/unified_single_mutant_antifold

test -s "$PLAN_DIR/unified_single_mutant_plan_gate.json"
test -s "$SCORE_DIR/antifold_model_run.json"
for view in experimental_vhh_only experimental_complex_context af3_vhh_only; do
  test -s "$SCORE_DIR/${view}.csv"
done
if [[ -e "$OUTPUT_DIR" || -e "$SUMMARY" ]]; then
  echo "Unified AntiFold output already exists; refusing to overwrite." >&2
  exit 1
fi
mkdir -p "$LOG_DIR"
sbatch scripts/candidate_design/submit_unified_single_mutant_antifold.slurm
