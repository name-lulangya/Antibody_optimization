#!/bin/bash
#SBATCH --job-name=nb252_flex_summary
#SBATCH --partition=batch
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=12
#SBATCH --time=00:30:00
#SBATCH --output=nb252-flex-summary-%j.log
#SBATCH --error=nb252-flex-summary-%j.log

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"

source /data/software/anaconda3/etc/profile.d/conda.sh
conda activate /data/software/env/luly25/ab_optim

python scripts/candidate_design/summarize_flex_ddg_pilot.py \
  --plan-dir docs/result_artifacts/candidate_design/flex_ddg_pilot_plan_20260812 \
  --task-root results/candidate_design/flex_ddg_pilot_20260812/tasks \
  --output-dir docs/result_artifacts/candidate_design/flex_ddg_pilot_result_20260812 \
  --run-summary docs/run_summaries/candidate_design/flex_ddg_pilot_result_20260812.json
