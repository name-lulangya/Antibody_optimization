#!/bin/bash
#SBATCH --job-name=nb252-property-affinity
#SBATCH --partition=batch
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=12
#SBATCH --time=02:30:00
#SBATCH --output=logs/property_affinity_pyrosetta/run-%j.log
#SBATCH --error=logs/property_affinity_pyrosetta/run-%j.log

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"

: "${RUN_KIND:?RUN_KIND must be pilot or full_scan}"
case "${RUN_KIND}" in
  pilot)
    OUTPUT_DIR="docs/result_artifacts/candidate_design/property_affinity_pyrosetta_pilot_20260816"
    RUN_SUMMARY="docs/run_summaries/candidate_design/property_affinity_pyrosetta_pilot_20260816.json"
    PILOT_ARGS=()
    ;;
  full_scan)
    OUTPUT_DIR="docs/result_artifacts/candidate_design/property_affinity_pyrosetta_full_scan_20260816"
    RUN_SUMMARY="docs/run_summaries/candidate_design/property_affinity_pyrosetta_full_scan_20260816.json"
    PILOT_ARGS=(--pilot-result-dir docs/result_artifacts/candidate_design/property_affinity_pyrosetta_pilot_20260816)
    ;;
  *)
    echo "RUN_KIND must be pilot or full_scan" >&2
    exit 2
    ;;
esac

source /data/software/anaconda3/etc/profile.d/conda.sh
set +u
conda activate /data/software/env/luly25/multi_ligand
set -u

python scripts/candidate_design/score_property_affinity_pyrosetta.py \
  --run-kind "${RUN_KIND}" \
  --plan-dir docs/result_artifacts/candidate_design/property_affinity_review_plan_20260816 \
  "${PILOT_ARGS[@]}" \
  --stage0-dir docs/result_artifacts/candidate_design/stage0_contract_20260810 \
  --structure-baseline-dir docs/result_artifacts/input_baseline/structure_released_20260810 \
  --calibration-dir docs/result_artifacts/structure_preparation/pyrosetta_scoring_calibration_v2_20260811 \
  --output-dir "${OUTPUT_DIR}" \
  --run-summary "${RUN_SUMMARY}" \
  --replicates 3 \
  --base-seed 8162000

set +u
conda activate /data/software/env/luly25/ab_optim
set -u
python scripts/candidate_design/render_property_affinity_pyrosetta.py \
  --run-kind "${RUN_KIND}" \
  --result-dir "${OUTPUT_DIR}"
