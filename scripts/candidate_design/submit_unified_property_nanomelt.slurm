#!/usr/bin/env bash
#SBATCH --job-name=nb252-prop-melt
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --time=01:00:00
#SBATCH --output=logs/unified_property_scoring/nanomelt-%j.log
#SBATCH --error=logs/unified_property_scoring/nanomelt-%j.log
set -eo pipefail
source /data/software/anaconda3/etc/profile.d/conda.sh
conda activate /data/software/env/luly25/nanomelt
set -u
python scripts/candidate_design/score_unified_property_nanomelt.py \
  --plan-dir docs/result_artifacts/candidate_design/unified_property_scoring_plan_20260815 \
  --nanomelt-executable /data/software/env/luly25/nanomelt/bin/nanomelt \
  --immune-builder-refine /data/software/env/luly25/nanomelt/lib/python3.10/site-packages/ImmuneBuilder/refine.py \
  --output-dir results/candidate_design/unified_property_scoring_20260815/nanomelt_scores
