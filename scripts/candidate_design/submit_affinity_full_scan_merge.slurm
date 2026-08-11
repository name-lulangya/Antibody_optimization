#!/bin/bash
#SBATCH --job-name=nb252_aff_merge
#SBATCH --partition=batch
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=12
#SBATCH --time=00:30:00
#SBATCH --output=nb252-affinity-merge-%j.log
#SBATCH --error=nb252-affinity-merge-%j.log

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"

source /data/software/anaconda3/etc/profile.d/conda.sh
conda activate /data/software/env/luly25/ab_optim

python scripts/candidate_design/merge_affinity_full_scan.py \
  --candidate-dir docs/result_artifacts/candidate_design/affinity_single_mutants_20260811 \
  --scan-plan-dir docs/result_artifacts/candidate_design/affinity_pyrosetta_full_scan_plan_20260811 \
  --shard-root results/candidate_design/affinity_pyrosetta_full_scan_20260811/shards \
  --output-dir docs/result_artifacts/candidate_design/affinity_pyrosetta_full_scan_20260811 \
  --run-summary docs/run_summaries/candidate_design/affinity_pyrosetta_full_scan_20260811.json
