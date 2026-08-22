#!/usr/bin/env bash
#SBATCH --job-name=nb252-expr-double
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00
#SBATCH --output=logs/expression_double_mutant_scan/run-%j.log
#SBATCH --error=logs/expression_double_mutant_scan/run-%j.log

set -eo pipefail

PLAN_DIR="docs/result_artifacts/candidate_design/expression_double_mutant_plan_20260822"
PARENT_DIR="docs/result_artifacts/candidate_design/expression_single_mutant_parent19_20260822"
SCORE_ROOT="results/candidate_design/expression_double_mutant_scan_20260822"
FINAL_DIR="docs/result_artifacts/candidate_design/expression_double_mutant_property_matrix_20260822"

source /data/software/anaconda3/etc/profile.d/conda.sh
export PYTHONUNBUFFERED=1

echo "[1/3] NetSolP: WT plus 162 complete double mutants"
conda activate /data/software/env/luly25/netsolp
set -u
python scripts/candidate_design/score_expression_double_mutant_properties.py \
  --tool netsolp \
  --plan-dir "${PLAN_DIR}" \
  --output-dir "${SCORE_ROOT}/netsolp" \
  --tool-root /homes/Tianlab/luly25/software/netsolp

echo "[2/3] NanoMelt: the identical WT-plus-162 sample table"
set +u
conda activate /data/software/env/luly25/nanomelt
set -u
python scripts/candidate_design/score_expression_double_mutant_properties.py \
  --tool nanomelt \
  --plan-dir "${PLAN_DIR}" \
  --output-dir "${SCORE_ROOT}/nanomelt" \
  --executable /data/software/env/luly25/nanomelt/bin/nanomelt

echo "[3/3] Merge complete matrix without selecting the final 11"
set +u
conda activate /data/software/env/luly25/ab_optim
set -u
python scripts/candidate_design/finalize_expression_double_mutant_matrix.py \
  --plan-dir "${PLAN_DIR}" \
  --parent19-dir "${PARENT_DIR}" \
  --netsolp-score-dir "${SCORE_ROOT}/netsolp" \
  --nanomelt-score-dir "${SCORE_ROOT}/nanomelt" \
  --output-dir "${FINAL_DIR}" \
  --run-summary docs/run_summaries/candidate_design/expression_double_mutant_property_matrix_20260822.json

echo "Expression double-mutant scan finished"
