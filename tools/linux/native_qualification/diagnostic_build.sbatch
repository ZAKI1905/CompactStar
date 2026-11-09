#!/bin/bash
#SBATCH --job-name=native-checkpoint-report-build
#SBATCH --partition=phyall
#SBATCH --nodelist=phyn004
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
set -euo pipefail
: "${RF_SOURCE:?}" "${RF_OUTPUT:?}" "${RF_RETAINED:?}" "${TC:?}"
cd "$RF_SOURCE"
export TMPDIR="$RF_OUTPUT/tmp" PYTHONDONTWRITEBYTECODE=1
mkdir -p "$TMPDIR"
scontrol show job "$SLURM_JOB_ID" > "$RF_OUTPUT/build-allocation.txt"
"$TC/python/bin/python" -I -B tools/linux/measure.py "$TC/python/bin/python" -B \
 tools/linux/native_qualification/build_diagnostic.py "$RF_RETAINED" "$RF_OUTPUT/build" > "$RF_OUTPUT/build.log" 2>&1
