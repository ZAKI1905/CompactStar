#!/bin/bash
#SBATCH --job-name=native-p2-bounded
#SBATCH --partition=phyall
#SBATCH --nodelist=phyn004
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
# Submission of this file requires the separately requested owner trajectory authorization.
set -euo pipefail
: "${RF_SOURCE:?}" "${RF_OUTPUT:?}" "${RF_SNAPSHOTS:?}" "${RF_YEARS:?}" "${TC:?}"
cd "$RF_SOURCE"
mkdir -p "$RF_OUTPUT"
scontrol show job "$SLURM_JOB_ID" > "$RF_OUTPUT/allocation.txt"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$TC/python/bin/python" -I -B tools/linux/measure.py \
 "$TC/python/bin/python" -B tools/linux/native_pair/snapshot.py stage \
 "$RF_SNAPSHOTS" "$RF_YEARS" "$RF_OUTPUT/run" > "$RF_OUTPUT/stage.log" 2>&1
