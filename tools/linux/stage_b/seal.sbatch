#!/bin/bash
set -euo pipefail
: "${RF_SOURCE:?}" "${RF_SEAL:?}" "${RF_SEAL_RUN:?}" "${RF_R1:?}" "${RF_R2:?}" "${RF_R3:?}"
TC=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain
export PATH="$TC/bin:$PATH" RF_SEAL_WORKERS=24
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$RF_SOURCE"
# Logs and allocation receipts must be outside all roots being sealed.
scontrol show job "$SLURM_JOB_ID" > "$RF_SEAL_RUN/allocation-$SLURM_JOB_ID.txt"
srun --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
 "$TC/python/bin/python" -B tools/linux/stage_b/seal.py "$RF_SEAL" "$RF_R1" "$RF_R2" "$RF_R3"
