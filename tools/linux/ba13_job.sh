#!/bin/bash
# Static BA13 only. Source eku_environment.sh before sbatch --export=ALL.
set -euo pipefail
SRC=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/sources/CompactStar
OUT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/ba13-recovery
CHAR="$TC/corrections/static-14487/characterization/run-_kg8zl8c/result.json"
export PYTHONDONTWRITEBYTECODE=1
stage=${1:?stage}
run=${2:?run}
args=()
case "$stage" in
 stars) dest="$OUT/$run/stars";;
 adapters) dest="$OUT/$run/adapters"; args=(--input "$OUT/$run/stars");;
 native-adapters) stage=adapters; dest="$OUT/$run/native-adapters"; args=(--input "$OUT/$run/stars" --native-tail);;
 coefficient) index=${SLURM_ARRAY_TASK_ID:?}; dest="$OUT/$run/coefficient-$index"; args=(--input "$OUT/$run/native-adapters" --index "$index");;
 *) exit 2;;
esac
cd "$SRC"
exec "$TC/python/bin/python" -B tools/linux/ba13_stage.py "$stage" --source "$SRC" --executable "$TC/build/CompactStar-Release/tests/phase6a1_frozen_sensitivity" --characterization "$CHAR" --output "$dest" "${args[@]}"
