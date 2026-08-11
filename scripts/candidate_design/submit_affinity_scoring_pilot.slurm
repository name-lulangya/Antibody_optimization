#!/bin/bash
#SBATCH --job-name=nb252_aff_pilot
#SBATCH --partition=batch
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=12
#SBATCH --time=02:00:00
#SBATCH --output=nb252-affinity-pilot-%j.log
#SBATCH --error=nb252-affinity-pilot-%j.log

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"

source /data/software/anaconda3/etc/profile.d/conda.sh
conda activate /data/software/env/luly25/multi_ligand

python scripts/candidate_design/score_affinity_candidates_pyrosetta.py \
  --candidate-dir docs/result_artifacts/candidate_design/affinity_single_mutants_20260811 \
  --candidate-id-file docs/result_artifacts/candidate_design/affinity_single_mutants_20260811/pilot_candidate_ids.txt \
  --stage0-dir docs/result_artifacts/candidate_design/stage0_contract_20260810 \
  --structure-baseline-dir docs/result_artifacts/input_baseline/structure_released_20260810 \
  --calibration-dir docs/result_artifacts/structure_preparation/pyrosetta_scoring_calibration_v2_20260811 \
  --output-dir docs/result_artifacts/candidate_design/affinity_pyrosetta_pilot_20260811 \
  --run-summary docs/run_summaries/candidate_design/affinity_pyrosetta_pilot_20260811.json \
  --replicates 3 \
  --base-seed 8112100
