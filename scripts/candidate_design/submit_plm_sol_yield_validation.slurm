#!/usr/bin/env bash
#SBATCH --job-name=nb252-plm-sol-yield
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --time=01:00:00
#SBATCH --output=logs/plm_sol_yield_validation/run-%j.log
#SBATCH --error=logs/plm_sol_yield_validation/run-%j.log

set -eo pipefail

PLAN_DIR="docs/result_artifacts/candidate_design/plm_sol_yield_validation_plan_20260819"
EMBEDDING_DIR="results/candidate_design/plm_sol_yield_validation_20260819/embeddings"
SCORE_DIR="results/candidate_design/plm_sol_yield_validation_20260819/model_scores"
ARTIFACT_DIR="docs/result_artifacts/candidate_design/plm_sol_yield_validation_result_20260819"
SUMMARY="docs/run_summaries/candidate_design/plm_sol_yield_validation_result_20260819.json"

source /data/software/anaconda3/etc/profile.d/conda.sh
conda activate /data/software/env/luly25/plm_sol_embed
python scripts/candidate_design/score_plm_sol_embeddings.py \
  --plan-dir "${PLAN_DIR}" \
  --output-dir "${EMBEDDING_DIR}"

conda activate /data/software/env/luly25/plm_sol
python scripts/candidate_design/score_plm_sol_classifier.py \
  --plan-dir "${PLAN_DIR}" \
  --embedding-dir "${EMBEDDING_DIR}" \
  --output-dir "${SCORE_DIR}"

conda activate /data/software/env/luly25/ab_optim
export MPLCONFIGDIR="${TMPDIR:-/tmp}/nb252-plm-sol-matplotlib-${SLURM_JOB_ID}"
python scripts/candidate_design/analyze_plm_sol_yield_validation.py \
  --plan-dir "${PLAN_DIR}" \
  --score-dir "${SCORE_DIR}" \
  --netsolp-samples docs/result_artifacts/candidate_design/netsolp_yield_validation_result_20260814/netsolp_yield_sample_evidence.csv \
  --rp3net-samples docs/result_artifacts/candidate_design/rp3net_yield_validation_result_20260818/rp3net_yield_sample_evidence.csv \
  --output-dir "${ARTIFACT_DIR}" \
  --run-summary "${SUMMARY}"
