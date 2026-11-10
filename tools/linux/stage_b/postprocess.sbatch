#!/bin/bash
set -euo pipefail
: "${RF_SOURCE:?}" "${RF_OUTPUT:?}" "${RF_PHASE:?}"
TC=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain
export PATH="$TC/bin:$PATH" RF_STAGE_WORKERS=24
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$RF_SOURCE"
scontrol show job "$SLURM_JOB_ID" > "$RF_OUTPUT/allocation-$SLURM_JOB_ID.txt"
run() { srun --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores "$TC/python/bin/python" -B "$@"; }
case "$RF_PHASE" in
 build)
  run tools/linux/stage_b/supplement.py "$RF_OUTPUT" build
  run tools/linux/stage_b/observation_trace.py "$RF_OUTPUT" build
  ;;
 source|control)
  run tools/linux/stage_b/supplement.py "$RF_OUTPUT" "$RF_PHASE"
  run tools/linux/stage_b/observation_trace.py "$RF_OUTPUT" "$RF_PHASE"
  ;;
 assess)
  run tools/linux/stage_b/residuals.py "$RF_OUTPUT" 24
  run tools/linux/stage_b/observation_residuals.py "$RF_OUTPUT" 24
  run tools/linux/stage_b/assess.py "$RF_OUTPUT" 24
  run tools/linux/stage_b/chemical_response.py "$RF_OUTPUT"
  run tools/linux/stage_b/events.py "$RF_OUTPUT"
  run tools/linux/stage_b/ledger.py "$RF_OUTPUT" 24
  run tools/linux/stage_b/local_reference.py "$RF_OUTPUT" 24
  run tools/linux/stage_b/accuracy.py "$RF_OUTPUT" 24
  ;;
 *) exit 2 ;;
esac
