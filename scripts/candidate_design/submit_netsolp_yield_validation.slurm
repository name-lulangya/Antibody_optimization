#!/usr/bin/env bash
#SBATCH --job-name=nb252-netsolp-yield
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --time=00:30:00
#SBATCH --output=logs/netsolp_yield_validation/run-%j.log
#SBATCH --error=logs/netsolp_yield_validation/run-%j.log

set -euo pipefail

PLAN_DIR="docs/result_artifacts/candidate_design/netsolp_yield_validation_plan_20260814"
SCORE_DIR="results/candidate_design/netsolp_yield_validation_20260814/model_scores"
ARTIFACT_DIR="docs/result_artifacts/candidate_design/netsolp_yield_validation_result_20260814"
SUMMARY="docs/run_summaries/candidate_design/netsolp_yield_validation_result_20260814.json"
NETSOLP_WORKDIR="/homes/Tianlab/luly25/software/netsolp"

source /data/software/anaconda3/etc/profile.d/conda.sh
conda activate /data/software/env/luly25/netsolp

python scripts/candidate_design/score_netsolp_sequences.py \
  --plan-dir "${PLAN_DIR}" \
  --netsolp-workdir "${NETSOLP_WORKDIR}" \
  --output-dir "${SCORE_DIR}" \
  --num-threads "${SLURM_CPUS_PER_TASK}"

conda activate /data/software/env/luly25/ab_optim
export MPLCONFIGDIR="${TMPDIR:-/tmp}/nb252-matplotlib-${SLURM_JOB_ID}"
python scripts/candidate_design/analyze_netsolp_yield_validation.py \
  --plan-dir "${PLAN_DIR}" \
  --score-dir "${SCORE_DIR}" \
  --output-dir "${ARTIFACT_DIR}" \
  --run-summary "${SUMMARY}"
