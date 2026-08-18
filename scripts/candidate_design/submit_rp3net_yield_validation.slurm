#!/usr/bin/env bash
#SBATCH --job-name=nb252-rp3net-yield
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --time=01:00:00
#SBATCH --output=logs/rp3net_yield_validation/run-%j.log
#SBATCH --error=logs/rp3net_yield_validation/run-%j.log

set -eo pipefail

PLAN_DIR="docs/result_artifacts/candidate_design/rp3net_yield_validation_plan_20260818"
SCORE_DIR="results/candidate_design/rp3net_yield_validation_20260818/model_scores"
ARTIFACT_DIR="docs/result_artifacts/candidate_design/rp3net_yield_validation_result_20260818"
SUMMARY="docs/run_summaries/candidate_design/rp3net_yield_validation_result_20260818.json"

source /data/software/anaconda3/etc/profile.d/conda.sh
conda activate /data/software/env/luly25/rp3net

python scripts/candidate_design/score_rp3net_sequences.py \
  --plan-dir "${PLAN_DIR}" \
  --output-dir "${SCORE_DIR}" \
  --device cuda:0 \
  --batch-size 8

conda activate /data/software/env/luly25/ab_optim
export MPLCONFIGDIR="${TMPDIR:-/tmp}/nb252-rp3net-matplotlib-${SLURM_JOB_ID}"
python scripts/candidate_design/analyze_rp3net_yield_validation.py \
  --plan-dir "${PLAN_DIR}" \
  --score-dir "${SCORE_DIR}" \
  --output-dir "${ARTIFACT_DIR}" \
  --run-summary "${SUMMARY}"
