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
#  TOV equation-level reference against the exact Schwarzschild interior. Self-contained.
# ---------------------------------------------------------------------------
add_executable(tov_reference_analytic core/tov_reference_analytic.cpp)
target_link_libraries(tov_reference_analytic PRIVATE CompactStar)
add_test(NAME tov_reference_analytic COMMAND tov_reference_analytic)
set_tests_properties(tov_reference_analytic PROPERTIES LABELS "core;tov;reference")

# ---------------------------------------------------------------------------
#  INV-12 cache-correctness contracts. Self-contained.
#  The registered test runs ONLY the supported-contract cases. The known-hazard
#  reproductions live behind --audit-known-hazards and are deliberately NOT a
#  CTest: they document defects, they do not assert that defects are correct.
# ---------------------------------------------------------------------------
add_executable(cache_contract evolution/cache_contract.cpp)
target_link_libraries(cache_contract PRIVATE CompactStar)
add_test(NAME cache_contract COMMAND cache_contract)
set_tests_properties(cache_contract PROPERTIES LABELS "evolution;cache;contract")

# ---------------------------------------------------------------------------
#  INV-12 cache contracts on the thermal path (Y_q, NeutrinoCooling payload).
#  Self-contained. Same two-mode split as cache_contract: only the supported
#  contracts are registered; --audit-known-hazards is documentation, not a test.
# ---------------------------------------------------------------------------
add_executable(cache_thermal_contract thermal/cache_thermal_contract.cpp)
target_link_libraries(cache_thermal_contract PRIVATE CompactStar)
add_test(NAME cache_thermal_contract COMMAND cache_thermal_contract)
set_tests_properties(cache_thermal_contract PROPERTIES LABELS "thermal;cache;contract")

# ---------------------------------------------------------------------------
#  Evolution stepper contract + thermal-only execution. Self-contained.
# ---------------------------------------------------------------------------
add_executable(evolution_stepper_contract thermal/evolution_stepper_contract.cpp)
target_link_libraries(evolution_stepper_contract PRIVATE CompactStar)
add_test(NAME evolution_stepper_contract COMMAND evolution_stepper_contract)
set_tests_properties(evolution_stepper_contract PROPERTIES LABELS "thermal")

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
    # Phase-2B passive-cooling regression baseline (GOVERNANCE.md §3.1 condition 7).
    # Requires the authenticated external CMF tables; ~11 s runtime.
    add_executable(passive_cooling_regression thermal/passive_cooling_regression.cpp)
    target_link_libraries(passive_cooling_regression PRIVATE CompactStar)
    add_test(NAME passive_cooling_regression
             COMMAND passive_cooling_regression "${COMPACTSTAR_EOS_DATA_ROOT}"
                     --compare
                     "${CMAKE_CURRENT_SOURCE_DIR}/baselines/passive_cooling_cmf_1p6_debug.tsv")
    set_tests_properties(passive_cooling_regression PROPERTIES
        LABELS "thermal;scientific;external-data;regression")

    # Phase-2B-2 TOV reference validation against the official CompOSE eos.mr.
    # Reference data are third-party and live beside the .eos table in the
    # authenticated distribution; the test SKIPs (77) if eos.mr is absent.
    add_executable(tov_reference_cmf core/tov_reference_cmf.cpp)
    target_link_libraries(tov_reference_cmf PRIVATE CompactStar)
    add_test(NAME tov_reference_cmf
             COMMAND tov_reference_cmf "${COMPACTSTAR_EOS_DATA_ROOT}")
    set_tests_properties(tov_reference_cmf PROPERTIES
        LABELS "core;tov;reference;external-data;scientific"
        SKIP_RETURN_CODE 77)

    message(STATUS "CompactStar: Tier-B thermal test enabled (EOS data root supplied)")
else()
    message(STATUS "CompactStar: Tier-B thermal test NOT registered "
                   "(set COMPACTSTAR_EOS_DATA_ROOT to enable)")
endif()
