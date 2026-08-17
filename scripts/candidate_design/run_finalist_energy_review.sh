#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$PROJECT_ROOT"

source /data/software/anaconda3/etc/profile.d/conda.sh
conda activate /data/software/env/luly25/ab_optim

python scripts/candidate_design/build_finalist_energy_review.py \
  --preliminary-dir docs/result_artifacts/candidate_design/preliminary_panel_20260817 \
  --affinity-result-dir docs/result_artifacts/candidate_design/affinity_pyrosetta_full_scan_20260811 \
  --property-result-dir docs/result_artifacts/candidate_design/property_affinity_pyrosetta_full_scan_20260816 \
  --double-pyrosetta-dir results/candidate_design/double_mutant_scan_20260816/pyrosetta \
  --output-dir docs/result_artifacts/candidate_design/finalist_energy_review_20260817 \
  --run-summary docs/run_summaries/candidate_design/finalist_energy_review_20260817.json
