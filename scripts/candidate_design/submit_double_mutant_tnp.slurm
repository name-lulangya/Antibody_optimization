#!/usr/bin/env bash
#SBATCH --job-name=nb252-double-tnp
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00
#SBATCH --output=logs/double_mutant_scan/tnp-%j.log
#SBATCH --error=logs/double_mutant_scan/tnp-%j.log
set -eo pipefail
cd "${SLURM_SUBMIT_DIR}"
source /data/software/anaconda3/etc/profile.d/conda.sh
conda activate /data/software/env/luly25/tnp
export PATH="/data/software/env/luly25/tnp/bin:${PATH}"
set -u
python scripts/candidate_design/score_double_mutant_properties.py --tool tnp --plan-dir docs/result_artifacts/candidate_design/double_mutant_plan_20260816 --tool-root /homes/Tianlab/luly25/software/TNP --executable /data/software/env/luly25/tnp/bin/TNP --immune-builder-refine /data/software/env/luly25/tnp/lib/python3.10/site-packages/ImmuneBuilder/refine.py --output-dir results/candidate_design/double_mutant_scan_20260816/tnp
