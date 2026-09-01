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

# ---------------------------------------------------------------------------
#  ADR-0002 PhotonCooling source-conformance (Phase 2A). Self-contained.
# ---------------------------------------------------------------------------
add_executable(photon_cooling_conformance thermal/photon_cooling_conformance.cpp)
target_link_libraries(photon_cooling_conformance PRIVATE CompactStar)
add_test(NAME photon_cooling_conformance COMMAND photon_cooling_conformance)
set_tests_properties(photon_cooling_conformance PROPERTIES LABELS "thermal")

# ---------------------------------------------------------------------------
#  Tier B — real canonical-star heat-capacity verification (ADR-0002 §V1).
#
#  Requires authenticated official CompOSE tables. Registered ONLY when
#  COMPACTSTAR_EOS_DATA_ROOT points at them, so a clean checkout still
#  configures and Tier A still runs with no external data. The test itself
#  exits non-zero if the data are missing — it never silently skips.
#
#  Expected layout under the root:
#     DS-CMF-1-with-crust/DS(CMF)-1_with_crust.eos
#     DNS-CMF-Hadronic-with-electrons/eos.{t,nb,yq,thermo}
# ---------------------------------------------------------------------------
set(COMPACTSTAR_EOS_DATA_ROOT "" CACHE PATH
    "Path to authenticated official CompOSE tables; enables the Tier-B thermal test.")

if(COMPACTSTAR_EOS_DATA_ROOT AND EXISTS "${COMPACTSTAR_EOS_DATA_ROOT}")
    add_executable(heat_capacity_real_star thermal/heat_capacity_real_star.cpp)
    target_link_libraries(heat_capacity_real_star PRIVATE CompactStar)
    add_test(NAME heat_capacity_real_star
             COMMAND heat_capacity_real_star "${COMPACTSTAR_EOS_DATA_ROOT}")
    set_tests_properties(heat_capacity_real_star PROPERTIES
        LABELS "thermal;external-data;scientific")
    # Phase-2B passive-cooling regression harness (GOVERNANCE.md §3.1 condition 7).
    #
    # BUILT BUT DELIBERATELY NOT REGISTERED AS A TEST. The authenticated live thermal
    # configuration cannot currently run: Config::stepper defaults to MSBDF, an implicit
    # BDF method, while GSLIntegrator::Integrate hard-codes sys.jacobian = nullptr, so GSL
    # calls a null function pointer and the process dies with SIGSEGV. See
    # docs/validation/PASSIVE_COOLING_BASELINE.md. Registering a segfaulting test would be
    # dishonest; the add_test() line belongs in the increment that repairs the integrator
    # and freezes the golden values in the same change.
    add_executable(passive_cooling_regression thermal/passive_cooling_regression.cpp)
    target_link_libraries(passive_cooling_regression PRIVATE CompactStar)

    message(STATUS "CompactStar: Tier-B thermal test enabled (EOS data root supplied)")
else()
    message(STATUS "CompactStar: Tier-B thermal test NOT registered "
                   "(set COMPACTSTAR_EOS_DATA_ROOT to enable)")
endif()
