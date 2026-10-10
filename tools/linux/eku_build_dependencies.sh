#!/bin/bash
set -euo pipefail
: "${TC:?}" "${CC:?}" "${CXX:?}"
SRC=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/sources
HANDOFF=/mnt/sdd/zaki/CompactStar/incoming/CompactStar-EKU-handoff-20261010
mkdir -p "$TC/source" "$TC/build/gsl" "$TC/install"
"$TC/bin/git" -C "$SRC/Zaki" rev-parse HEAD
"$TC/bin/git" -C "$SRC/Zaki" status --porcelain
if [ ! -d "$TC/source/gsl-2.7.1" ]; then tar -xzf "$HANDOFF/upstream/gsl-2.7.1.tar.gz" -C "$TC/source"; fi
cd "$TC/build/gsl"
"$TC/source/gsl-2.7.1/configure" --prefix="$TC/install/gsl-2.7.1" --disable-shared --enable-static CC="$CC" CFLAGS='-O2 -g -ffp-contract=off' > "$TC/logs/gsl-configure.log" 2>&1
make -j4 > "$TC/logs/gsl-build.log" 2>&1
set +e
make -j4 check > "$TC/logs/gsl-check.log" 2>&1
rc=$?
set -e
echo "$rc" > "$TC/logs/gsl-check.exit"
make install > "$TC/logs/gsl-install.log" 2>&1
GSL="$TC/install/gsl-2.7.1"
common=(-DCMAKE_C_COMPILER="$CC" -DCMAKE_CXX_COMPILER="$CXX" -DCMAKE_INSTALL_LIBDIR=lib -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_C_FLAGS=-ffp-contract=off -DCMAKE_CXX_FLAGS=-ffp-contract=off -DGIT_EXECUTABLE="$TC/bin/git" -DGSL_ROOT_DIR="$GSL" -DGSL_INCLUDE_DIR="$GSL/include" -DGSL_LIBRARY="$GSL/lib/libgsl.a" -DGSL_CBLAS_LIBRARY="$GSL/lib/libgslcblas.a" -DZLIB_INCLUDE_DIR="$TC/sysroot/usr/include" -DZLIB_LIBRARY="$TC/sysroot/usr/lib64/libz.so" -DCMAKE_FIND_USE_PACKAGE_REGISTRY=OFF -DCMAKE_FIND_USE_SYSTEM_PACKAGE_REGISTRY=OFF)
for mode in Debug Release; do
 for dep in Zaki CONFIND; do
  extra=(-DZAKI_BUILD_TESTS=ON -DZAKI_BUILD_EXAMPLES=OFF)
  if [ "$dep" = CONFIND ]; then extra=(-DCONFIND_BUILD_TESTS=ON -DZaki_DIR="$TC/install/$mode/Zaki/lib/cmake/Zaki"); fi
  cmake -S "$SRC/$dep" -B "$TC/build/$dep-$mode" -DCMAKE_BUILD_TYPE="$mode" -DCMAKE_INSTALL_PREFIX="$TC/install/$mode/$dep" "${common[@]}" "${extra[@]}" > "$TC/logs/$dep-$mode-configure.log" 2>&1
  cmake --build "$TC/build/$dep-$mode" --parallel 4 > "$TC/logs/$dep-$mode-build.log" 2>&1
  cmake --install "$TC/build/$dep-$mode" > "$TC/logs/$dep-$mode-install.log" 2>&1
  set +e
  ctest --test-dir "$TC/build/$dep-$mode" --output-on-failure --output-junit "$TC/logs/$dep-$mode-tests.xml" -j1 > "$TC/logs/$dep-$mode-tests.log" 2>&1
  rc=$?
  set -e
  echo "$rc" > "$TC/logs/$dep-$mode-tests.exit"
 done
done
