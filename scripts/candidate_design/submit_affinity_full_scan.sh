#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")/../.."
mkdir -p results/candidate_design/affinity_pyrosetta_full_scan_20260811/shards

/data/software/env/luly25/ab_optim/bin/python -c "import matplotlib, numpy"

ARRAY_JOB_ID=$(sbatch --parsable scripts/candidate_design/submit_affinity_full_scan_array.slurm)
ARRAY_JOB_ID=${ARRAY_JOB_ID%%;*}
MERGE_JOB_ID=$(
  sbatch --parsable \
    --dependency="afterok:${ARRAY_JOB_ID}" \
    scripts/candidate_design/submit_affinity_full_scan_merge.slurm
)
MERGE_JOB_ID=${MERGE_JOB_ID%%;*}

echo "Affinity full-scan array job: ${ARRAY_JOB_ID}"
echo "Dependent merge job: ${MERGE_JOB_ID}"
