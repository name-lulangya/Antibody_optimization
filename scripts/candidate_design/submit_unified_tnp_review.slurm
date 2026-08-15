#!/usr/bin/env bash
#SBATCH --job-name=nb252-unified-tnp
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00
#SBATCH --output=logs/unified_tnp_review/run-%j.log
#SBATCH --error=logs/unified_tnp_review/run-%j.log

set -eo pipefail

PLAN_DIR="docs/result_artifacts/candidate_design/unified_tnp_review_plan_20260815"
SCORE_DIR="results/candidate_design/unified_tnp_review_20260815/model_scores"
ARTIFACT_DIR="docs/result_artifacts/candidate_design/unified_tnp_review_result_20260815"
SUMMARY="docs/run_summaries/candidate_design/unified_tnp_review_result_20260815.json"
TNP_SOURCE="/homes/Tianlab/luly25/software/TNP"
TNP_EXECUTABLE="/data/software/env/luly25/tnp/bin/TNP"
IMMUNE_BUILDER_REFINE="/data/software/env/luly25/tnp/lib/python3.10/site-packages/ImmuneBuilder/refine.py"

source /data/software/anaconda3/etc/profile.d/conda.sh
conda activate /data/software/env/luly25/tnp
export PATH="/data/software/env/luly25/tnp/bin:${PATH}"
export PYTHONPATH="${TNP_SOURCE}"
set -u

python scripts/candidate_design/score_unified_tnp_candidates.py \
  --plan-dir "${PLAN_DIR}" \
  --tnp-source "${TNP_SOURCE}" \
  --tnp-executable "${TNP_EXECUTABLE}" \
  --immune-builder-refine "${IMMUNE_BUILDER_REFINE}" \
  --output-dir "${SCORE_DIR}"

set +u
conda activate /data/software/env/luly25/ab_optim
set -u
export MPLCONFIGDIR="${TMPDIR:-/tmp}/nb252-matplotlib-${SLURM_JOB_ID}"
python scripts/candidate_design/analyze_unified_tnp_review.py \
  --plan-dir "${PLAN_DIR}" \
  --score-dir "${SCORE_DIR}" \
  --output-dir "${ARTIFACT_DIR}" \
  --run-summary "${SUMMARY}"
