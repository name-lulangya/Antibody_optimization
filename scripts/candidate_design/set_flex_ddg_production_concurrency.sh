#!/usr/bin/env bash
set -euo pipefail

SUBMISSION_DIR=""
CONCURRENCY=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --submission-dir) SUBMISSION_DIR="$2"; shift 2 ;;
    --concurrency) CONCURRENCY="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done
if [[ -z "${SUBMISSION_DIR}" || ! -f "${SUBMISSION_DIR}/submission_jobs.tsv" ]]; then
  echo "A valid --submission-dir is required" >&2
  exit 2
fi
if [[ ! "${CONCURRENCY}" =~ ^[1-9][0-9]*$ ]]; then
  echo "--concurrency must be a positive integer" >&2
  exit 2
fi

while IFS=$'\t' read -r kind job_id; do
  [[ "${kind}" == "array" ]] || continue
  scontrol update JobId="${job_id}" ArrayTaskThrottle="${CONCURRENCY}"
  echo "Updated array job ${job_id} throttle to ${CONCURRENCY}"
done < "${SUBMISSION_DIR}/submission_jobs.tsv"
