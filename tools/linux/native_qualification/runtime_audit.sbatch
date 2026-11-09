#!/bin/bash
#SBATCH --job-name=native-p2-runtime-audit
#SBATCH --partition=phyall
#SBATCH --nodelist=phyn004
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=4G
#SBATCH --time=00:02:00
# Read-only dependency audit after job 19171; no assembly or evolution call.
set -euo pipefail
: "${RF_SOURCE:?}" "${RF_OUTPUT:?}" "${RF_SNAPSHOTS:?}" "${TC:?}"
cd "$RF_SOURCE"
scontrol show job "$SLURM_JOB_ID" > "$RF_OUTPUT/runtime-audit-allocation.txt"
"$TC/python/bin/python" -I -B tools/linux/native_qualification/runtime_audit.py \
  "$RF_SNAPSHOTS" > "$RF_OUTPUT/compute-runtime-audit.json"
