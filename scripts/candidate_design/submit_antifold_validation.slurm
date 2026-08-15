#!/usr/bin/env bash
#SBATCH --job-name=nb252-antifold
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --time=00:30:00
#SBATCH --output=logs/antifold_validation/run-%j.log
#SBATCH --error=logs/antifold_validation/run-%j.log

set -euo pipefail

PLAN_DIR="docs/result_artifacts/candidate_design/antifold_validation_plan_20260815"
SCORE_DIR="results/candidate_design/antifold_validation_20260815/model_scores"
ARTIFACT_DIR="docs/result_artifacts/candidate_design/antifold_validation_result_20260815"
SUMMARY="docs/run_summaries/candidate_design/antifold_validation_result_20260815.json"

source /data/software/anaconda3/etc/profile.d/conda.sh
conda activate /data/software/env/luly25/antifold
export PYTHONUNBUFFERED=1

python scripts/candidate_design/score_antifold_validation.py \
  --plan-dir "${PLAN_DIR}" \
  --output-dir "${SCORE_DIR}"

conda activate /data/software/env/luly25/ab_optim
export MPLCONFIGDIR="${TMPDIR:-/tmp}/nb252-matplotlib-${SLURM_JOB_ID}"

python scripts/candidate_design/analyze_antifold_validation.py \
  --plan-dir "${PLAN_DIR}" \
  --score-dir "${SCORE_DIR}" \
  --output-dir "${ARTIFACT_DIR}" \
  --run-summary "${SUMMARY}"
