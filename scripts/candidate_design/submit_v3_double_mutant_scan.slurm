#!/usr/bin/env bash
#SBATCH --job-name=nb252-v3-double
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00
#SBATCH --output=logs/v3_double_mutant_scan/run-%j.log
#SBATCH --error=logs/v3_double_mutant_scan/run-%j.log

set -eo pipefail

PLAN_DIR="docs/result_artifacts/candidate_design/v3_double_mutant_plan_20260825"
SCORE_ROOT="results/candidate_design/v3_double_mutant_scan_20260825"
FINAL_DIR="docs/result_artifacts/candidate_design/v3_double_mutant_property_matrix_20260825"

source /data/software/anaconda3/etc/profile.d/conda.sh
export PYTHONUNBUFFERED=1

echo "[1/3] NetSolP: WT plus all 102 complete V3 doubles"
conda activate /data/software/env/luly25/netsolp
set -u
python scripts/candidate_design/score_v3_double_mutant_properties.py \
  --tool netsolp \
  --plan-dir "${PLAN_DIR}" \
  --output-dir "${SCORE_ROOT}/netsolp" \
  --netsolp-workdir /homes/Tianlab/luly25/software/netsolp

echo "[2/3] NanoMelt: the identical WT-plus-102 ID/sequence table"
set +u
conda activate /data/software/env/luly25/nanomelt
set -u
python scripts/candidate_design/score_v3_double_mutant_properties.py \
  --tool nanomelt \
  --plan-dir "${PLAN_DIR}" \
  --output-dir "${SCORE_ROOT}/nanomelt" \
  --nanomelt-executable /data/software/env/luly25/nanomelt/bin/nanomelt \
  --immune-builder-refine /data/software/env/luly25/nanomelt/lib/python3.10/site-packages/ImmuneBuilder/refine.py

echo "[3/3] Join the complete matrix without selecting the final 15 doubles"
set +u
conda activate /data/software/env/luly25/ab_optim
set -u
python scripts/candidate_design/finalize_v3_double_mutant_matrix.py \
  --plan-dir "${PLAN_DIR}" \
  --netsolp-score-dir "${SCORE_ROOT}/netsolp" \
  --nanomelt-score-dir "${SCORE_ROOT}/nanomelt" \
  --output-dir "${FINAL_DIR}" \
  --run-summary docs/run_summaries/candidate_design/v3_double_mutant_property_matrix_20260825/run_summary.json

echo "V3 complete double-mutant scoring and matrix join finished"
