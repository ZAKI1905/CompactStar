#!/bin/bash
set -euo pipefail
export GSL_ROOT_DIR="$TC/install/gsl-2.7.1"
export CMAKE_PREFIX_PATH="$TC/sysroot/usr:$GSL_ROOT_DIR:$CMAKE_PREFIX_PATH"
cd "$TC/build/gsl"
set +e
bash "$TC/source/gsl-2.7.1/pkgconfig.test" > "$TC/logs/gsl-pkgconfig-rerun.log" 2>&1
rc=$?
echo "$rc" > "$TC/logs/gsl-pkgconfig-rerun.exit"
for mode in Debug Release; do
 ctest --test-dir "$TC/build/CONFIND-$mode" -R '^installed_package_consumer$' --output-on-failure > "$TC/logs/CONFIND-$mode-consumer-rerun.log" 2>&1
 rc=$?
 echo "$rc" > "$TC/logs/CONFIND-$mode-consumer-rerun.exit"
done
