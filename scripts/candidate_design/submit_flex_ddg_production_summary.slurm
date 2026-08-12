#!/usr/bin/env bash
#SBATCH --job-name=nb252-flex-prod-summary
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --time=00:30:00
#SBATCH --output=logs/flex_ddg_production/summary-%j.log
#SBATCH --error=logs/flex_ddg_production/summary-%j.log

set -euo pipefail
source /data/software/anaconda3/etc/profile.d/conda.sh
conda activate /data/software/env/luly25/ab_optim

python scripts/candidate_design/summarize_flex_ddg_production.py \
  --plan-dir docs/result_artifacts/candidate_design/flex_ddg_production_plan_20260812 \
  --task-root results/candidate_design/flex_ddg_production_20260812/tasks \
  --output-dir docs/result_artifacts/candidate_design/flex_ddg_production_result_20260812 \
  --run-summary docs/run_summaries/candidate_design/flex_ddg_production_result_20260812.json
