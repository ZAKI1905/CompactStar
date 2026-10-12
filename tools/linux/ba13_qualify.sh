#!/bin/bash
# Static-only recipe. Source eku_environment.sh on login before sbatch --export=ALL.
# A complete candidate is not an uncertainty-qualified or promoted certificate.
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1 PYTHONOPTIMIZE=0
SRC=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/sources/CompactStar
RAW=${1:?native raw directory}
REPEAT=${2:?independent repeat raw directory}
ARCHIVE=${3:?verified extracted historical package directory}
OUT=${4:?fresh qualification output directory}
mkdir "$OUT"
cd "$SRC"
PY="$TC/python/bin/python"
MAC="$ARCHIVE/original/phase6a1-frozen-certificate-final"
"$PY" -B tools/linux/measure.py "$PY" -B tools/linux/ba13_fit.py --work "$RAW" --output "$OUT/native" > "$OUT/native-fit.log" 2>&1
"$PY" -B tools/linux/measure.py "$PY" -B tools/linux/ba13_fit.py --work "$REPEAT" --output "$OUT/repeat" > "$OUT/repeat-fit.log" 2>&1
cmp "$OUT/native/certificate.json" "$OUT/repeat/certificate.json"
cmp "$OUT/native/certificate.tsv" "$OUT/repeat/certificate.tsv"
"$PY" -B tools/linux/measure.py "$PY" -B tools/linux/ba13_fit.py --work "$MAC" --historical-package "$ARCHIVE" --output "$OUT/mac-replay" > "$OUT/mac-replay.log" 2>&1
"$PY" -B tools/linux/measure.py "$PY" -B tools/linux/ba13_refit.py --raw "$RAW" --certificate "$OUT/native/certificate.json" --output "$OUT/native-refit.json" > "$OUT/native-refit.log" 2>&1
"$PY" -B tools/linux/measure.py "$PY" -B tools/linux/ba13_refit.py --raw "$MAC" --certificate "$MAC/certificate.json" --output "$OUT/mac-refit.json" > "$OUT/mac-refit.log" 2>&1
"$PY" -B tools/linux/measure.py "$PY" -B tools/linux/ba13_platform_compare.py --native "$RAW" --mac "$MAC" --native-certificate "$OUT/native/certificate.json" --mac-certificate "$MAC/certificate.json" --output "$OUT/platform-comparison.json" > "$OUT/platform-comparison.log" 2>&1
"$PY" -B tools/linux/measure.py "$TC/build/CompactStar-Release/tests/phase6a1_frozen_monitor" "$OUT/native/certificate.tsv" > "$OUT/monitor.log" 2>&1
"$PY" -B tools/linux/measure.py "$TC/build/CompactStar-Release/tests/phase6a1_ba13_native_certificate_binding" "$RAW/chemical-profiles/target-0" "$RAW/chemical-profiles/target-0/certificate.txt" "$OUT/native/certificate.tsv" "$OUT/native-binding-work" > "$OUT/native-binding.log" 2>&1
echo 'STATIC_CANDIDATE_COMPLETE; PLATFORM_AGREEMENT_UNQUALIFIED; NO PROMOTION OR TRAJECTORY'
