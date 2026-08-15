#!/usr/bin/env bash
#SBATCH --job-name=nb252-prop-net
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00
#SBATCH --output=logs/unified_property_scoring/netsolp-%j.log
#SBATCH --error=logs/unified_property_scoring/netsolp-%j.log
set -eo pipefail
source /data/software/anaconda3/etc/profile.d/conda.sh
conda activate /data/software/env/luly25/netsolp
set -u
python scripts/candidate_design/score_unified_property_netsolp.py \
  --plan-dir docs/result_artifacts/candidate_design/unified_property_scoring_plan_20260815 \
  --netsolp-workdir /homes/Tianlab/luly25/software/netsolp \
  --output-dir results/candidate_design/unified_property_scoring_20260815/netsolp_scores \
  --num-threads 12
