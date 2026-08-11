#!/bin/bash
#SBATCH --job-name=nb252_aff_full
#SBATCH --partition=batch
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=12
#SBATCH --array=0-11%4
#SBATCH --time=02:00:00
#SBATCH --output=nb252-affinity-full-%A_%a.log
#SBATCH --error=nb252-affinity-full-%A_%a.log

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"

source /data/software/anaconda3/etc/profile.d/conda.sh
conda activate /data/software/env/luly25/multi_ligand

printf -v SHARD_ID "shard_%02d" "${SLURM_ARRAY_TASK_ID}"
PLAN_DIR="docs/result_artifacts/candidate_design/affinity_pyrosetta_full_scan_plan_20260811"
SHARD_ROOT="results/candidate_design/affinity_pyrosetta_full_scan_20260811/shards"
OUTPUT_DIR="${SHARD_ROOT}/${SHARD_ID}"
mkdir -p "${OUTPUT_DIR}"

python scripts/candidate_design/score_affinity_candidates_pyrosetta.py \
  --run-kind full_scan_shard \
  --shard-id "${SHARD_ID}" \
  --candidate-dir docs/result_artifacts/candidate_design/affinity_single_mutants_20260811 \
  --candidate-id-file "${PLAN_DIR}/${SHARD_ID}_candidate_ids.txt" \
  --scan-plan-dir "${PLAN_DIR}" \
  --pilot-v2-review docs/result_artifacts/candidate_design/affinity_pyrosetta_pilot_v2_20260811/affinity_pilot_v2_scientific_review.json \
  --stage0-dir docs/result_artifacts/candidate_design/stage0_contract_20260810 \
  --structure-baseline-dir docs/result_artifacts/input_baseline/structure_released_20260810 \
  --calibration-dir docs/result_artifacts/structure_preparation/pyrosetta_scoring_calibration_v2_20260811 \
  --output-dir "${OUTPUT_DIR}" \
  --run-summary "${OUTPUT_DIR}/shard_run_summary.json" \
  --replicates 3 \
  --base-seed 8112100
