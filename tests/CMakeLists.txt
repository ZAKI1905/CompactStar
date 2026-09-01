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
