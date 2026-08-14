#!/usr/bin/env bash
#SBATCH --job-name=nb252-tnp-v2
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --time=01:00:00
#SBATCH --output=logs/tnp_yield_validation_v2/run-%j.log
#SBATCH --error=logs/tnp_yield_validation_v2/run-%j.log

set -euo pipefail

PLAN_DIR="docs/result_artifacts/candidate_design/tnp_yield_validation_plan_v2_20260814"
SCORE_DIR="results/candidate_design/tnp_yield_validation_v2_20260814/model_scores"
ARTIFACT_DIR="docs/result_artifacts/candidate_design/tnp_yield_validation_result_v2_20260814"
SUMMARY="docs/run_summaries/candidate_design/tnp_yield_validation_result_v2_20260814.json"
NETSOLP_EVIDENCE="docs/result_artifacts/candidate_design/netsolp_yield_validation_result_20260814/netsolp_yield_sample_evidence.csv"
TNP_SOURCE="/homes/Tianlab/luly25/software/TNP"
TNP_EXECUTABLE="/data/software/env/luly25/tnp/bin/TNP"
IMMUNE_BUILDER_REFINE="/data/software/env/luly25/tnp/lib/python3.10/site-packages/ImmuneBuilder/refine.py"

source /data/software/anaconda3/etc/profile.d/conda.sh
conda activate /data/software/env/luly25/tnp
export PATH="/data/software/env/luly25/tnp/bin:${PATH}"
export PYTHONPATH="${TNP_SOURCE}"

python scripts/candidate_design/score_tnp_sequences.py \
  --plan-dir "${PLAN_DIR}" \
  --tnp-source "${TNP_SOURCE}" \
  --tnp-executable "${TNP_EXECUTABLE}" \
  --immune-builder-refine "${IMMUNE_BUILDER_REFINE}" \
  --output-dir "${SCORE_DIR}"

conda activate /data/software/env/luly25/ab_optim
export MPLCONFIGDIR="${TMPDIR:-/tmp}/nb252-matplotlib-${SLURM_JOB_ID}"
python scripts/candidate_design/analyze_tnp_yield_validation.py \
  --plan-dir "${PLAN_DIR}" \
  --score-dir "${SCORE_DIR}" \
  --netsolp-evidence "${NETSOLP_EVIDENCE}" \
  --output-dir "${ARTIFACT_DIR}" \
  --run-summary "${SUMMARY}"
