#!/usr/bin/env bash
#SBATCH --job-name=nb252-uni-af
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --time=00:30:00
#SBATCH --output=logs/unified_single_mutant_antifold/run-%j.log
#SBATCH --error=logs/unified_single_mutant_antifold/run-%j.log

set -eo pipefail
source /data/software/anaconda3/etc/profile.d/conda.sh
conda activate /data/software/env/luly25/ab_optim
set -u

python scripts/candidate_design/analyze_unified_single_mutant_antifold.py \
  --plan-dir docs/result_artifacts/candidate_design/unified_single_mutant_plan_20260815 \
  --antifold-plan-dir docs/result_artifacts/candidate_design/antifold_validation_plan_20260815 \
  --score-dir results/candidate_design/antifold_validation_20260815/model_scores \
  --output-dir docs/result_artifacts/candidate_design/unified_single_mutant_antifold_20260815 \
  --run-summary docs/run_summaries/candidate_design/unified_single_mutant_antifold_20260815.json
