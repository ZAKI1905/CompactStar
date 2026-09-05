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
#  Phase 4C-I0 EOS thermodynamic-derivative authority (ADR-0007 P5). Self-contained:
#  builds its own synthetic EOS tables in a temp directory. Contract only — it exercises
#  no O(Omega^2) physics.
# ---------------------------------------------------------------------------
add_executable(eos_derivative_contract eos/eos_derivative_contract.cpp)
target_link_libraries(eos_derivative_contract PRIVATE CompactStar)
add_test(NAME eos_derivative_contract COMMAND eos_derivative_contract)
set_tests_properties(eos_derivative_contract PROPERTIES LABELS "eos;derivative;contract")

# ---------------------------------------------------------------------------
#  Phase 5A-2 accepted ADR-0010 local cold charge-neutral thermodynamic
#  contract. Self-contained analytic free leptons plus independent toy and
#  projection fixtures; no external EOS data and no scientific baseline.
# ---------------------------------------------------------------------------
add_executable(rotochemical_local_thermodynamics
    eos/rotochemical_local_thermodynamics.cpp)
target_link_libraries(rotochemical_local_thermodynamics PRIVATE CompactStar)
add_test(NAME rotochemical_local_thermodynamics
         COMMAND rotochemical_local_thermodynamics)
set_tests_properties(rotochemical_local_thermodynamics PROPERTIES
    LABELS "rotochemical;thermodynamics;self-contained;scientific;contract")

# ---------------------------------------------------------------------------
#  Phase 5A-3 source-faithful Track-R cold ideal npe-mu local provider.
#  Self-contained analytic/numerical validation only: no star, external EOS
#  data, electrostatic projection, or scientific baseline.
# ---------------------------------------------------------------------------
add_executable(rotochemical_trackr_freegas_local
    eos/rotochemical_trackr_freegas_local.cpp)
target_link_libraries(rotochemical_trackr_freegas_local PRIVATE CompactStar)
add_test(NAME rotochemical_trackr_freegas_local
         COMMAND rotochemical_trackr_freegas_local)
set_tests_properties(rotochemical_trackr_freegas_local PROPERTIES
    LABELS "rotochemical;thermodynamics;track-r;free-gas;self-contained;scientific;contract")

# Phase 5A-4 R1: active npe chart and species boundaries; independent local oracles.
add_executable(rotochemical_trackr_npe eos/rotochemical_trackr_npe.cpp)
target_link_libraries(rotochemical_trackr_npe PRIVATE CompactStar)
add_test(NAME rotochemical_trackr_npe COMMAND rotochemical_trackr_npe)
set_tests_properties(rotochemical_trackr_npe PROPERTIES
    LABELS "rotochemical;thermodynamics;track-r;free-gas;npe;self-contained;scientific;contract")

# Phase 5A-5: vacuum, active pe chart, neutron threshold and pe/npe limit.
add_executable(rotochemical_trackr_pe eos/rotochemical_trackr_pe.cpp)
target_link_libraries(rotochemical_trackr_pe PRIVATE CompactStar)
add_test(NAME rotochemical_trackr_pe COMMAND rotochemical_trackr_pe)
set_tests_properties(rotochemical_trackr_pe PROPERTIES
    LABELS "rotochemical;thermodynamics;track-r;free-gas;pe;self-contained;scientific;contract")

# ---------------------------------------------------------------------------
#  Phase 4C-I1 governed O(Omega^2) monopole IMPLEMENTATION contract (ADR-0007).
#  Self-contained. Protects contracts — fail-closed inputs, regular-centre start,
#  provenance/staleness, seed invariance, materialization algebra, solve accounting.
#  It is NOT independent physical validation of the ODE; that is Phase 4D.
# ---------------------------------------------------------------------------
add_executable(hartle_monopole_contract rotation/hartle_monopole_contract.cpp)
target_include_directories(hartle_monopole_contract PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/rotation)
target_link_libraries(hartle_monopole_contract PRIVATE CompactStar)
add_test(NAME hartle_monopole_contract COMMAND hartle_monopole_contract)
set_tests_properties(hartle_monopole_contract PROPERTIES
    LABELS "rotation;hartle;monopole;contract")

# ---------------------------------------------------------------------------
#  Phase 4D INDEPENDENT physical validation of the governed O(Omega^2) monopole response on
#  the exact constant-density background (ADR-0007 §7). Self-contained. The oracles are a
#  test-only (m0,h0) solver and a continuum solver in hartle_monopole_reference.hpp.
# ---------------------------------------------------------------------------
add_executable(hartle_monopole_physics_analytic rotation/hartle_monopole_physics_analytic.cpp)
target_include_directories(hartle_monopole_physics_analytic PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/rotation)
target_link_libraries(hartle_monopole_physics_analytic PRIVATE CompactStar)
add_test(NAME hartle_monopole_physics_analytic COMMAND hartle_monopole_physics_analytic)
set_tests_properties(hartle_monopole_physics_analytic PROPERTIES
    LABELS "rotation;hartle;monopole;scientific;reference")

# Phase 4D-RI implementation contract of the MEASURE-COMPLETE EOS energy-density source
# (ADR-0008, ACCEPTED 2026-09-03): constant-density interior measure identically zero, smooth-EOS
# equivalence with the superseded differential form, per-segment measure identity and
# same-partition source accounting, terminal surface atom applied exactly once. Self-contained.
add_executable(hartle_monopole_measure_contract rotation/hartle_monopole_measure_contract.cpp)
target_include_directories(hartle_monopole_measure_contract PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/rotation)
target_link_libraries(hartle_monopole_measure_contract PRIVATE CompactStar)
add_test(NAME hartle_monopole_measure_contract COMMAND hartle_monopole_measure_contract)
set_tests_properties(hartle_monopole_measure_contract PROPERTIES
    LABELS "rotation;hartle;monopole;contract;measure"
    TIMEOUT 900)

# Phase 4D Experiment K: the governed O(Omega^2) monopole response against PUBLISHED
# second-order results — Chandrasekhar & Miller 1974 Table I (homogeneous stars, analytic
# background) and Hartle & Thorne 1968 Tables 3/5 on the printed Harrison-Wheeler EOS
# (historical fixture embedded in tests/rotation/hartle_thorne_1968_hw_eos.hpp). Self-contained.
add_executable(hartle_monopole_published rotation/hartle_monopole_published.cpp)
target_include_directories(hartle_monopole_published PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/rotation)
target_link_libraries(hartle_monopole_published PRIVATE CompactStar)
add_test(NAME hartle_monopole_published COMMAND hartle_monopole_published)
set_tests_properties(hartle_monopole_published PROPERTIES
    LABELS "rotation;hartle;monopole;scientific;published"
    TIMEOUT 900)

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
#  Phase-2B-4B scale-free first-order Hartle moment of inertia. Self-contained.
# ---------------------------------------------------------------------------
add_executable(hartle_moment_inertia_analytic rotation/hartle_moment_inertia_analytic.cpp)
target_include_directories(hartle_moment_inertia_analytic PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/rotation)
target_link_libraries(hartle_moment_inertia_analytic PRIVATE CompactStar)
add_test(NAME hartle_moment_inertia_analytic COMMAND hartle_moment_inertia_analytic)
set_tests_properties(hartle_moment_inertia_analytic PROPERTIES
    LABELS "rotation;hartle;scientific")

# ---------------------------------------------------------------------------
#  Phase-4A ADR-0006 first-order physical-normalization contract. Self-contained.
#  Tolerances are predeclared in ADR-0006 SS7 and must not be widened here.
# ---------------------------------------------------------------------------
add_executable(hartle_normalization_contract rotation/hartle_normalization_contract.cpp)
target_include_directories(hartle_normalization_contract PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/rotation)
target_link_libraries(hartle_normalization_contract PRIVATE CompactStar)
add_test(NAME hartle_normalization_contract COMMAND hartle_normalization_contract)
set_tests_properties(hartle_normalization_contract PROPERTIES
    LABELS "rotation;hartle;scientific;contract")

# ---------------------------------------------------------------------------
#  Phase-4B independent PHYSICAL validation of the normalized first-order Hartle
#  response (profile shape, exterior identities, volume identity, weak-field
#  limit). Self-contained. Bounds are predeclared in the source from the Phase-2B
#  record and must not be widened here.
# ---------------------------------------------------------------------------
add_executable(hartle_first_order_physics_analytic
    rotation/hartle_first_order_physics_analytic.cpp)
target_include_directories(hartle_first_order_physics_analytic PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/rotation)
target_link_libraries(hartle_first_order_physics_analytic PRIVATE CompactStar)
add_test(NAME hartle_first_order_physics_analytic
         COMMAND hartle_first_order_physics_analytic)
set_tests_properties(hartle_first_order_physics_analytic PROPERTIES
    LABELS "rotation;hartle;scientific;reference")

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
#  ADR-0004 canonical proper-volume / metric primitive. Self-contained.
#  Exercises CompactStar/Geometry.hpp directly: the regular-centre limit, the
#  fail-closed domain, the metric factor against two independent algebraic
#  forms, and the w_V composition. No EOS assets, no TOV solve.
# ---------------------------------------------------------------------------
add_executable(proper_volume_contract core/proper_volume_contract.cpp)
target_link_libraries(proper_volume_contract PRIVATE CompactStar)
add_test(NAME proper_volume_contract COMMAND proper_volume_contract)
set_tests_properties(proper_volume_contract PROPERTIES
    LABELS "core;geometry;contract;scientific")

# ---------------------------------------------------------------------------
#  ADR-0004 GeometryCache conformance. Self-contained.
#  Pins the stored-Lambda and derived-Lambda routes to the same cached arrays,
#  ties the cached WV to the canonical primitive, and proves WVExpNu/WVExp2Nu
#  are composed from the ONE WV array rather than independently re-derived.
# ---------------------------------------------------------------------------
add_executable(geometry_cache_measure_contract evolution/geometry_cache_measure_contract.cpp)
target_link_libraries(geometry_cache_measure_contract PRIVATE CompactStar)
add_test(NAME geometry_cache_measure_contract COMMAND geometry_cache_measure_contract)
set_tests_properties(geometry_cache_measure_contract PROPERTIES
    LABELS "evolution;geometry;cache;contract")

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

    # Phase-2B-4A nonrotating TOV -> cooling grid-convergence harness.
    # Registered only with the authenticated external data root.
    add_executable(grid_convergence_cmf thermal/grid_convergence_cmf.cpp)
    target_link_libraries(grid_convergence_cmf PRIVATE CompactStar)
    add_test(NAME grid_convergence_cmf
             COMMAND grid_convergence_cmf "${COMPACTSTAR_EOS_DATA_ROOT}")
    set_tests_properties(grid_convergence_cmf PROPERTIES
        LABELS "thermal;scientific;external-data;convergence"
        TIMEOUT 1800)

    # Phase 4D INDEPENDENT physical validation of the governed O(Omega^2) monopole response on
    # the authenticated CMF sequence (ADR-0007 §7): (m0,h0) reference solver, exterior identity,
    # radial convergence, EOS-derivative sensitivity, homogeneous sequence-derivative identity.
    add_executable(hartle_monopole_physics_cmf rotation/hartle_monopole_physics_cmf.cpp)
    target_include_directories(hartle_monopole_physics_cmf PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/rotation)
    target_link_libraries(hartle_monopole_physics_cmf PRIVATE CompactStar)
    add_test(NAME hartle_monopole_physics_cmf
             COMMAND hartle_monopole_physics_cmf "${COMPACTSTAR_EOS_DATA_ROOT}")
    set_tests_properties(hartle_monopole_physics_cmf PROPERTIES
        LABELS "rotation;hartle;monopole;scientific;external-data;reference"
        TIMEOUT 1800)

    # Phase 4D-BL production regression, established AFTER independent verification.
    # This baseline is not an independent physics oracle.
    add_executable(hartle_monopole_regression rotation/hartle_monopole_regression.cpp)
    target_link_libraries(hartle_monopole_regression PRIVATE CompactStar)
    target_compile_definitions(hartle_monopole_regression PRIVATE
        COMPACTSTAR_BASELINE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/baselines")
    add_test(NAME hartle_monopole_regression
             COMMAND hartle_monopole_regression "${COMPACTSTAR_EOS_DATA_ROOT}")
    set_tests_properties(hartle_monopole_regression PROPERTIES
        LABELS "rotation;hartle;monopole;scientific;regression;external-data"
        TIMEOUT 900)

    # Phase 4C-I1 governed O(Omega^2) monopole conformance on the authenticated CMF
    # sequence (ADR-0007). Structural conformance and diagnostics only — the independent
    # scientific validation is Phase 4D, and no golden is produced here.
    add_executable(hartle_monopole_cmf rotation/hartle_monopole_cmf.cpp)
    target_include_directories(hartle_monopole_cmf PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/rotation)
    target_link_libraries(hartle_monopole_cmf PRIVATE CompactStar)
    add_test(NAME hartle_monopole_cmf
             COMMAND hartle_monopole_cmf "${COMPACTSTAR_EOS_DATA_ROOT}")
    set_tests_properties(hartle_monopole_cmf PROPERTIES
        LABELS "rotation;hartle;monopole;external-data"
        TIMEOUT 900)

    # Phase 4C-I0 EOS d(eps)/dp on the authenticated CMF sequence (ADR-0007 P5).
    # Asserts completeness, finiteness and positivity; prints the retired profile
    # finite difference as a diagnostic only. No O(Omega^2) physics.
    add_executable(eos_derivative_cmf eos/eos_derivative_cmf.cpp)
    target_link_libraries(eos_derivative_cmf PRIVATE CompactStar)
    add_test(NAME eos_derivative_cmf
             COMMAND eos_derivative_cmf "${COMPACTSTAR_EOS_DATA_ROOT}")
    set_tests_properties(eos_derivative_cmf PROPERTIES
        LABELS "eos;derivative;scientific;external-data"
        TIMEOUT 900)

    # Phase-2B-4B scale-free Hartle I on the authenticated CMF sequence.
    add_executable(hartle_moment_inertia_cmf rotation/hartle_moment_inertia_cmf.cpp)
    target_include_directories(hartle_moment_inertia_cmf PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/rotation)
    target_link_libraries(hartle_moment_inertia_cmf PRIVATE CompactStar)
    add_test(NAME hartle_moment_inertia_cmf
             COMMAND hartle_moment_inertia_cmf "${COMPACTSTAR_EOS_DATA_ROOT}")
    set_tests_properties(hartle_moment_inertia_cmf PROPERTIES
        LABELS "rotation;hartle;scientific;external-data"
        TIMEOUT 900)

    # Phase-4A ADR-0006 normalization contract on the authenticated CMF sequence,
    # plus the V9 slow-rotation regime diagnostic (reported, never enforced).
    add_executable(hartle_normalization_cmf rotation/hartle_normalization_cmf.cpp)
    target_include_directories(hartle_normalization_cmf PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/rotation)
    target_link_libraries(hartle_normalization_cmf PRIVATE CompactStar)
    add_test(NAME hartle_normalization_cmf
             COMMAND hartle_normalization_cmf "${COMPACTSTAR_EOS_DATA_ROOT}")
    set_tests_properties(hartle_normalization_cmf PROPERTIES
        LABELS "rotation;hartle;scientific;external-data;contract"
        TIMEOUT 900)

    # Phase-4B independent physical validation of the normalized first-order
    # response on the authenticated DS(CMF)-1 sequence.
    add_executable(hartle_first_order_physics_cmf
        rotation/hartle_first_order_physics_cmf.cpp)
    target_include_directories(hartle_first_order_physics_cmf PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/rotation)
    target_link_libraries(hartle_first_order_physics_cmf PRIVATE CompactStar)
    add_test(NAME hartle_first_order_physics_cmf
             COMMAND hartle_first_order_physics_cmf "${COMPACTSTAR_EOS_DATA_ROOT}")
    set_tests_properties(hartle_first_order_physics_cmf PROPERTIES
        LABELS "rotation;hartle;scientific;external-data;reference"
        TIMEOUT 900)

    # Phase-3D INV-14 / ADR-0004 baryon-number reference on the authenticated
    # DS(CMF)-1 sequence (1.0 / 1.4 / 1.6 / 2.0 Msun). Guards the proper-volume
    # migration under the tolerance predeclared in ADR-0004 SS7.1 before
    # implementation: |dB|/B <= 1.0e-15.
    add_executable(baryon_number_cmf core/baryon_number_cmf.cpp)
    target_link_libraries(baryon_number_cmf PRIVATE CompactStar)
    add_test(NAME baryon_number_cmf
             COMMAND baryon_number_cmf "${COMPACTSTAR_EOS_DATA_ROOT}"
                     --compare
                     "${CMAKE_CURRENT_SOURCE_DIR}/baselines/baryon_number_dscmf1_reference.tsv")
    set_tests_properties(baryon_number_cmf PROPERTIES
        LABELS "core;tov;reference;external-data;scientific"
        TIMEOUT 900)

    # Phase-3E-0 equivalence measurement of the two LIVE TOV paths.
    # Path 1 (TOVSolver::Solve + NStar::Append + FinalizeSurface) vs Path 2
    # (SingleStarSolveToTOVPoints + NStar(points,labels) + BuildFromTOV) at the
    # same frozen central densities. Measurement only: it designates no canonical
    # owner. Criteria are predeclared in the source header.
    add_executable(tov_path_equivalence_cmf core/tov_path_equivalence_cmf.cpp)
    target_link_libraries(tov_path_equivalence_cmf PRIVATE CompactStar)
    add_test(NAME tov_path_equivalence_cmf
             COMMAND tov_path_equivalence_cmf "${COMPACTSTAR_EOS_DATA_ROOT}")
    set_tests_properties(tov_path_equivalence_cmf PROPERTIES
        LABELS "core;tov;scientific;external-data;equivalence"
        TIMEOUT 1800)

    # Phase-3E-I1 workflow/interface contract for TOVSolver::Solve.
    # ADR-0005 Q1/Q2: Solve() stays a public sequence orchestrator and its
    # _Sequence.tsv output is a compatibility contract. All six live callers
    # consume only that file, and one reads it back as program input, so ONE
    # focused test protects the contract -- deliberately not one per caller.
    add_executable(tov_sequence_workflow_cmf core/tov_sequence_workflow_cmf.cpp)
    target_link_libraries(tov_sequence_workflow_cmf PRIVATE CompactStar)
    add_test(NAME tov_sequence_workflow_cmf
             COMMAND tov_sequence_workflow_cmf "${COMPACTSTAR_EOS_DATA_ROOT}")
    set_tests_properties(tov_sequence_workflow_cmf PROPERTIES
        LABELS "core;tov;workflow;interface;external-data"
        TIMEOUT 900)

    # Phase-3E-I4 coverage for TOVSolver::GenTestSequence -- the last remaining
    # caller of the duplicate ordinary-star radial loop. Committed BEFORE the
    # migration it justifies, so its equivalence assertion is evidence about the
    # legacy RadiusLoop implementation rather than a description of the new one.
    add_executable(tov_gen_test_sequence_cmf core/tov_gen_test_sequence_cmf.cpp)
    target_link_libraries(tov_gen_test_sequence_cmf PRIVATE CompactStar)
    add_test(NAME tov_gen_test_sequence_cmf
             COMMAND tov_gen_test_sequence_cmf "${COMPACTSTAR_EOS_DATA_ROOT}")
    set_tests_properties(tov_gen_test_sequence_cmf PROPERTIES
        LABELS "core;tov;scientific;external-data;equivalence"
        TIMEOUT 3600)

    message(STATUS "CompactStar: Tier-B thermal test enabled (EOS data root supplied)")
else()
    message(STATUS "CompactStar: Tier-B thermal test NOT registered "
                   "(set COMPACTSTAR_EOS_DATA_ROOT to enable)")
endif()

# ADR-0009 surface-event contracts. Reports are disposable build outputs only.
add_executable(tov_surface_contract core/tov_surface_contract.cpp)
target_link_libraries(tov_surface_contract PRIVATE CompactStar)
foreach(mode sweep audit derivatives consumers contracts)
    add_test(NAME tov_surface_${mode}_hw COMMAND tov_surface_contract HW ${mode}
        "${CMAKE_CURRENT_BINARY_DIR}/tov_surface_${mode}_hw.tsv")
    set_tests_properties(tov_surface_${mode}_hw PROPERTIES LABELS "core;tov;surface;scientific" TIMEOUT 900)
endforeach()
if(COMPACTSTAR_EOS_DATA_ROOT AND EXISTS "${COMPACTSTAR_EOS_DATA_ROOT}")
    foreach(mode v7a v7b sweep audit derivatives consumers contracts)
        add_test(NAME tov_surface_${mode}_cmf COMMAND tov_surface_contract
            "${COMPACTSTAR_EOS_DATA_ROOT}/DS-CMF-1-with-crust/DS(CMF)-1_with_crust.eos" ${mode}
            "${CMAKE_CURRENT_BINARY_DIR}/tov_surface_${mode}_cmf.tsv")
        set_tests_properties(tov_surface_${mode}_cmf PROPERTIES LABELS "core;tov;surface;scientific;external-data" TIMEOUT 900)
    endforeach()
    add_test(NAME tov_surface_grid_cmf COMMAND grid_convergence_cmf
        "${COMPACTSTAR_EOS_DATA_ROOT}" --surface-dry-run)
    set_tests_properties(tov_surface_grid_cmf PROPERTIES LABELS "core;tov;surface;thermal;scientific;external-data" TIMEOUT 1200)
endif()

# Track-R Structure-1 candidate; all EOS data generated locally.
add_executable(rotochemical_trackr_freegas_structure eos/rotochemical_trackr_freegas_structure.cpp)
target_link_libraries(rotochemical_trackr_freegas_structure PRIVATE CompactStar)
add_executable(rotochemical_trackr_freegas_barotrope eos/rotochemical_trackr_freegas_barotrope.cpp)
target_link_libraries(rotochemical_trackr_freegas_barotrope PRIVATE CompactStar)
add_test(NAME rotochemical_trackr_freegas_barotrope COMMAND rotochemical_trackr_freegas_barotrope)
add_test(NAME rotochemical_trackr_freegas_structure COMMAND rotochemical_trackr_freegas_structure)
set_tests_properties(rotochemical_trackr_freegas_barotrope rotochemical_trackr_freegas_structure
    PROPERTIES LABELS "rotochemical;thermodynamics;track-r;free-gas;structure;self-contained;scientific")
