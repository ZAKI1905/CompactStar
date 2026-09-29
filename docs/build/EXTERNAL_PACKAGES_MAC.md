# Qualified external packages on the local Mac

The migration candidate uses explicit, authenticated Zaki 2.0.1 and CONFIND
2.0.0 packages. Historical `dependencies/` archives and headers are retained
as scientific oracles and are not build or install inputs. See ADR-0019 and
[the migration qualification report](../validation/COMPACTSTAR_EXTERNAL_DEPENDENCY_MIGRATION_QUALIFICATION.md)
for authority and exact artifact hashes.

CMake 3.18 or later is required. The qualified machine uses AppleClang 21,
CMake 4.2.1, `/opt/local` GSL 2.7.1 and the macOS SDK's zlib. Python 3.12.10
with NumPy/SciPy/mpmath runs the existing Python tests; it is not linked into
any numerical executable. Packages must already exist; configuration never
downloads or substitutes dependencies.

For the current candidate, use these out-of-source roots:

```sh
migration_packages=/Users/keeper/Documents/CompactStar/external/qualification/compactstar-migration/e263a6e-b0cbd510
cmake -S . -B /private/tmp/compactstar-debug \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCOMPACTSTAR_ZAKI_PREFIX="$migration_packages/Debug/Zaki" \
  -DCOMPACTSTAR_CONFIND_PREFIX="$migration_packages/Debug/CONFIND" \
  -DGSL_ROOT_DIR=/opt/local \
  -DGSL_CONFIG_EXECUTABLE=/opt/local/bin/gsl-config \
  -DGSL_INCLUDE_DIR=/opt/local/include \
  -DGSL_LIBRARY=/opt/local/lib/libgsl.dylib \
  -DGSL_CBLAS_LIBRARY=/opt/local/lib/libgslcblas.dylib \
  -DPython3_EXECUTABLE=/Users/keeper/miniforge3/bin/python3 \
  -DCOMPACTSTAR_EOS_DATA_ROOT=/Users/keeper/Documents/CompactStar/data/compose
cmake --build /private/tmp/compactstar-debug -j6
ctest --test-dir /private/tmp/compactstar-debug --output-on-failure
```

For Release, change the build directory, build type and both package paths
to Release. Also supply the Debug packages used by nested governed Phase-5D
qualification:

```sh
-DCOMPACTSTAR_QUALIFICATION_DEBUG_ZAKI_PREFIX="$migration_packages/Debug/Zaki"
-DCOMPACTSTAR_QUALIFICATION_DEBUG_CONFIND_PREFIX="$migration_packages/Debug/CONFIND"
```

Debug and Release packages cannot be mixed. `dependency-provenance.json` in
the build directory records what was authenticated. `cmake/dependency-lock.json`
is the immutable identity list for this qualification, not a promise that
arbitrary rebuilds have identical archive bytes. Updating it requires explicit
package provenance and renewed same-mode numerical qualification.

Three existing tests have Debug-only reference authority: the Phase-5B
structural-response regression, Phase-5C coefficient regression, and
Phase-5D coupled-oracle certificate. They run in Debug. The Phase-5D1
regression constructs its own fresh Debug build even when invoked from a
Release build. The migration report additionally compares same-mode Release
artifacts exactly; it never treats Release as required to equal Debug. The
existing Hartle-monopole and baryon-number tests use separately frozen OLD
Release migration references in Release; their Debug baselines and explicit
ADR-0012 reference overrides remain unchanged.

Use `ctest -R '^taskmanager_'` for the fast sequence/observer/replay checks
and the real stellar T2 test. T2 requires the authenticated external EOS.
The temporary working root is short and TaskManager uses one thread.

The two previously undefined spin diagnostic entry points now throw explicit
`logic_error` exceptions rather than relying on dynamic lookup. They remain
scientifically unimplemented. The dipole normalization needs separate
scientific qualification; this migration does not select it.
