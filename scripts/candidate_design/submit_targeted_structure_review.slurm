#!/bin/bash
#SBATCH --job-name=nb252-targeted-structure
#SBATCH --partition=batch
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=12
#SBATCH --time=01:00:00
#SBATCH --output=logs/targeted_structure_review/run-%j.log
#SBATCH --error=logs/targeted_structure_review/run-%j.log

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"

source /data/software/anaconda3/etc/profile.d/conda.sh
set +u
conda activate /data/software/env/luly25/multi_ligand
set -u

python scripts/candidate_design/score_targeted_structure_review_pyrosetta.py \
  --plan-dir docs/result_artifacts/candidate_design/targeted_structure_review_plan_20260816 \
  --calibration-dir docs/result_artifacts/structure_preparation/pyrosetta_scoring_calibration_v2_20260811 \
  --output-dir results/candidate_design/targeted_structure_review_runtime_20260816 \
  --run-summary docs/run_summaries/candidate_design/targeted_structure_review_runtime_20260816.json \
  --replicates 3 \
  --base-seed 8163000

set +u
conda activate /data/software/env/luly25/ab_optim
set -u

python scripts/candidate_design/analyze_targeted_structure_review.py \
  --safety-review-dir docs/result_artifacts/candidate_design/single_mutant_safety_review_20260816 \
  --plan-dir docs/result_artifacts/candidate_design/targeted_structure_review_plan_20260816 \
  --runtime-result-dir results/candidate_design/targeted_structure_review_runtime_20260816 \
  --output-dir docs/result_artifacts/candidate_design/targeted_structure_review_result_20260816 \
  --run-summary docs/run_summaries/candidate_design/targeted_structure_review_result_20260816.json
