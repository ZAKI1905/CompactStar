# ---------------------------------------------------------------------------
#  Canonical automated-test root for CompactStar.
#
#  This directory holds tests that run unattended under CTest. It is distinct
#  from main/Test/, which holds manual demo and debugging programs and is NOT
#  the automated suite.
#
#  Reached only when BUILD_TESTING is ON (standard CMake/CTest semantics).
#
#  Phase 1 scope: build/link/run infrastructure only. No scientific reference
#  values belong here; those are roadmap Phase 2A/2B work.
# ---------------------------------------------------------------------------

add_executable(compactstar_library_smoke smoke/compactstar_library_smoke.cpp)

target_link_libraries(compactstar_library_smoke PRIVATE CompactStar)

# Test binaries are build artifacts only — never installed.
add_test(NAME compactstar_library_smoke COMMAND compactstar_library_smoke)

# ---------------------------------------------------------------------------
#  Thermal verification (roadmap Phase 2A-1, ADR-0002 §V1 Tier A).
#
#  Self-contained: builds its own synthetic CompOSE fixture in a temp directory
#  and requires NO external EOS assets, so it runs in a clean checkout.
#  It is verification evidence, not a scientific baseline.
# ---------------------------------------------------------------------------
add_executable(heat_capacity_v1 thermal/heat_capacity_v1.cpp)
target_link_libraries(heat_capacity_v1 PRIVATE CompactStar)
add_test(NAME heat_capacity_v1 COMMAND heat_capacity_v1)
