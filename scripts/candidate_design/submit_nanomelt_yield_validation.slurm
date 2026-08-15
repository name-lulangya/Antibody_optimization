#!/usr/bin/env bash
#SBATCH --job-name=nb252-nanomelt-yield
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00
#SBATCH --output=logs/nanomelt_yield_validation/run-%j.log
#SBATCH --error=logs/nanomelt_yield_validation/run-%j.log

set -euo pipefail

PLAN_DIR="docs/result_artifacts/candidate_design/nanomelt_yield_validation_plan_20260815"
SCORE_DIR="results/candidate_design/nanomelt_yield_validation_20260815/model_scores"
ARTIFACT_DIR="docs/result_artifacts/candidate_design/nanomelt_yield_validation_result_20260815"
SUMMARY="docs/run_summaries/candidate_design/nanomelt_yield_validation_result_20260815.json"
NANOMELT_EXECUTABLE="/data/software/env/luly25/nanomelt/bin/nanomelt"
IMMUNE_BUILDER_REFINE="/data/software/env/luly25/nanomelt/lib/python3.10/site-packages/ImmuneBuilder/refine.py"

source /data/software/anaconda3/etc/profile.d/conda.sh
conda activate /data/software/env/luly25/nanomelt
export TORCH_HOME="/homes/Tianlab/luly25/.cache/torch"
export PYTHONUNBUFFERED=1

python scripts/candidate_design/score_nanomelt_sequences.py \
  --plan-dir "${PLAN_DIR}" \
  --nanomelt-executable "${NANOMELT_EXECUTABLE}" \
  --immune-builder-refine "${IMMUNE_BUILDER_REFINE}" \
  --output-dir "${SCORE_DIR}"

conda activate /data/software/env/luly25/ab_optim
export MPLCONFIGDIR="${TMPDIR:-/tmp}/nb252-matplotlib-${SLURM_JOB_ID}"
python scripts/candidate_design/analyze_nanomelt_yield_validation.py \
  --plan-dir "${PLAN_DIR}" \
  --score-dir "${SCORE_DIR}" \
  --output-dir "${ARTIFACT_DIR}" \
  --run-summary "${SUMMARY}"
