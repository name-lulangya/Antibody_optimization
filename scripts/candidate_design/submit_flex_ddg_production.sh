#!/usr/bin/env bash
set -euo pipefail

CONCURRENCY=12
CHUNK_SIZE=900
while [[ $# -gt 0 ]]; do
  case "$1" in
    --concurrency) CONCURRENCY="$2"; shift 2 ;;
    --chunk-size) CHUNK_SIZE="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done
if [[ ! "${CONCURRENCY}" =~ ^[1-9][0-9]*$ ]]; then
  echo "--concurrency must be a positive integer" >&2
  exit 2
fi
if [[ ! "${CHUNK_SIZE}" =~ ^[1-9][0-9]*$ ]]; then
  echo "--chunk-size must be a positive integer" >&2
  exit 2
fi

PLAN_DIR="docs/result_artifacts/candidate_design/flex_ddg_production_plan_20260812"
TASK_ROOT="results/candidate_design/flex_ddg_production_20260812/tasks"
SUBMISSION_ROOT="results/candidate_design/flex_ddg_production_20260812/submissions"
LOG_DIR="logs/flex_ddg_production"
PRECHECK_DIR="results/candidate_design/flex_ddg_production_20260812/prechecks"
STAMP=$(date -u +%Y%m%dT%H%M%SZ)
SUBMISSION_DIR="${SUBMISSION_ROOT}/submission_${STAMP}"
CODE_REVISION=$(git rev-parse --verify HEAD)
PRECHECK_MARKER="${PRECHECK_DIR}/${CODE_REVISION}.pass"

mkdir -p "${TASK_ROOT}" "${SUBMISSION_ROOT}" "${LOG_DIR}" "${PRECHECK_DIR}"

if [[ ! -f "${PRECHECK_MARKER}" ]]; then
  source /data/software/anaconda3/etc/profile.d/conda.sh
  conda activate /data/software/env/luly25/multi_ligand
  python scripts/candidate_design/run_flex_ddg_task_pyrosetta.py \
    --run-kind production \
    --plan-dir "${PLAN_DIR}" \
    --task-index 0 \
    --stage0-dir docs/result_artifacts/candidate_design/stage0_contract_20260810 \
    --structure-baseline-dir docs/result_artifacts/input_baseline/structure_released_20260810 \
    --calibration-dir docs/result_artifacts/structure_preparation/pyrosetta_scoring_calibration_v2_20260811 \
    --output-dir results/candidate_design/flex_ddg_production_20260812/precheck-not-written \
    --check_only
  printf 'pass\t%s\n' "${CODE_REVISION}" > "${PRECHECK_MARKER}.tmp"
  mv "${PRECHECK_MARKER}.tmp" "${PRECHECK_MARKER}"
else
  echo "Reusing successful real-data precheck for revision ${CODE_REVISION}."
fi

source /data/software/anaconda3/etc/profile.d/conda.sh
conda activate /data/software/env/luly25/ab_optim
python scripts/candidate_design/prepare_flex_ddg_production_resume.py \
  --plan-dir "${PLAN_DIR}" \
  --task-root "${TASK_ROOT}" \
  --submission-dir "${SUBMISSION_DIR}" \
  --concurrency "${CONCURRENCY}" \
  --chunk-size "${CHUNK_SIZE}"

JOBS_FILE="${SUBMISSION_DIR}/submission_jobs.tsv"
: > "${JOBS_FILE}"
PREVIOUS_JOB_ID=""
for INDEX_FILE in "${SUBMISSION_DIR}"/chunk_*.txt; do
  [[ -f "${INDEX_FILE}" ]] || continue
  TASK_COUNT=$(wc -l < "${INDEX_FILE}")
  ARRAY_SPEC="0-$((TASK_COUNT - 1))%${CONCURRENCY}"
  DEPENDENCY=()
  if [[ -n "${PREVIOUS_JOB_ID}" ]]; then
    DEPENDENCY=(--dependency="afterok:${PREVIOUS_JOB_ID}")
  fi
  JOB_ID=$(sbatch --parsable \
    --kill-on-invalid-dep=yes \
    --array="${ARRAY_SPEC}" \
    --export="ALL,FLEX_DDG_INDEX_FILE=${INDEX_FILE}" \
    "${DEPENDENCY[@]}" \
    scripts/candidate_design/submit_flex_ddg_production_array.slurm)
  JOB_ID="${JOB_ID%%;*}"
  printf 'array\t%s\n' "${JOB_ID}" >> "${JOBS_FILE}"
  PREVIOUS_JOB_ID="${JOB_ID}"
done

if [[ -z "${PREVIOUS_JOB_ID}" ]]; then
  echo "All 1000 tasks are already complete; no arrays submitted."
  SUMMARY_JOB_ID=$(sbatch --parsable scripts/candidate_design/submit_flex_ddg_production_summary.slurm)
  SUMMARY_JOB_ID="${SUMMARY_JOB_ID%%;*}"
  printf 'summary\t%s\n' "${SUMMARY_JOB_ID}" >> "${JOBS_FILE}"
  echo "Submitted summary job ${SUMMARY_JOB_ID}."
  exit 0
fi
SUMMARY_JOB_ID=$(sbatch --parsable \
  --kill-on-invalid-dep=yes \
  --dependency="afterok:${PREVIOUS_JOB_ID}" \
  scripts/candidate_design/submit_flex_ddg_production_summary.slurm)
SUMMARY_JOB_ID="${SUMMARY_JOB_ID%%;*}"
printf 'summary\t%s\n' "${SUMMARY_JOB_ID}" >> "${JOBS_FILE}"
echo "Submitted production arrays listed in ${JOBS_FILE}"
echo "Change all submitted array throttles with:"
echo "bash scripts/candidate_design/set_flex_ddg_production_concurrency.sh --submission-dir ${SUBMISSION_DIR} --concurrency N"
echo "Summary job ${SUMMARY_JOB_ID} will run after the final array succeeds."
