#!/usr/bin/env bash
#SBATCH --job-name=nb252-nanobert-yield
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --time=01:00:00
#SBATCH --output=logs/nanobert_yield_validation/run-%j.log
#SBATCH --error=logs/nanobert_yield_validation/run-%j.log

set -euo pipefail

PLAN_DIR="docs/result_artifacts/candidate_design/nanobert_yield_validation_plan_20260813"
SCORE_DIR="results/candidate_design/nanobert_yield_validation_20260813/model_scores"
ARTIFACT_DIR="docs/result_artifacts/candidate_design/nanobert_yield_validation_result_20260813"
SUMMARY="docs/run_summaries/candidate_design/nanobert_yield_validation_result_20260813.json"
MODEL_SNAPSHOT="/homes/Tianlab/luly25/software/huggingface_cache/vhh-lm/hub/models--NaturalAntibody--nanoBERT/snapshots/edc8182ad89a827f8737fa572c6b5fac6197e6b0"

source /data/software/anaconda3/etc/profile.d/conda.sh
conda activate /data/software/env/luly25/vhh-lm
export HF_HOME=/homes/Tianlab/luly25/software/huggingface_cache/vhh-lm
export TRANSFORMERS_OFFLINE=1
export HF_HUB_OFFLINE=1

python scripts/candidate_design/score_nanobert_sequences.py \
  --plan-dir "${PLAN_DIR}" \
  --model-snapshot "${MODEL_SNAPSHOT}" \
  --output-dir "${SCORE_DIR}" \
  --mask-batch-size 64

conda activate /data/software/env/luly25/ab_optim
export MPLCONFIGDIR="${TMPDIR:-/tmp}/nb252-matplotlib-${SLURM_JOB_ID}"
python scripts/candidate_design/analyze_nanobert_yield_validation.py \
  --plan-dir "${PLAN_DIR}" \
  --score-dir "${SCORE_DIR}" \
  --output-dir "${ARTIFACT_DIR}" \
  --run-summary "${SUMMARY}"
