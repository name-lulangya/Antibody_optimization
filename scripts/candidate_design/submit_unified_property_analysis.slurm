#!/usr/bin/env bash
#SBATCH --job-name=nb252-prop-merge
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --time=00:30:00
#SBATCH --output=logs/unified_property_scoring/analysis-%j.log
#SBATCH --error=logs/unified_property_scoring/analysis-%j.log
set -eo pipefail
source /data/software/anaconda3/etc/profile.d/conda.sh
conda activate /data/software/env/luly25/ab_optim
set -u
python scripts/candidate_design/analyze_unified_property_scoring.py \
  --plan-dir docs/result_artifacts/candidate_design/unified_property_scoring_plan_20260815 \
  --unified-antifold-dir docs/result_artifacts/candidate_design/unified_single_mutant_antifold_20260815 \
  --netsolp-score-dir results/candidate_design/unified_property_scoring_20260815/netsolp_scores \
  --nanomelt-score-dir results/candidate_design/unified_property_scoring_20260815/nanomelt_scores \
  --output-dir docs/result_artifacts/candidate_design/unified_property_scoring_result_20260815 \
  --run-summary docs/run_summaries/candidate_design/unified_property_scoring_result_20260815.json
