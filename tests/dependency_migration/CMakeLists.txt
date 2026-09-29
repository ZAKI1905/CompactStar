add_executable(taskmanager_migration taskmanager_driver.cpp)
target_link_libraries(taskmanager_migration PRIVATE CompactStar)
add_executable(taskmanager_selection_replay replay_taskmanager.cpp)
target_link_libraries(taskmanager_selection_replay PRIVATE CompactStar)
if(APPLE)
    target_link_options(taskmanager_migration PRIVATE "-Wl,-map,${CMAKE_CURRENT_BINARY_DIR}/taskmanager.map")
endif()
foreach(mode T1 T2 replay)
    if(mode STREQUAL "replay")
        set(fixture_target taskmanager_selection_replay)
    else()
        set(fixture_target taskmanager_migration)
    endif()
    add_test(NAME taskmanager_${mode}
        COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_fixture.py"
            --exe $<TARGET_FILE:${fixture_target}> --build-mode "${CMAKE_BUILD_TYPE}"
            --mode ${mode}
            --eos "${COMPACTSTAR_EOS_DATA_ROOT}/DS-CMF-1-with-crust/DS(CMF)-1_with_crust.eos"
            --report "${CMAKE_CURRENT_BINARY_DIR}/${mode}-result.json")
    set_tests_properties(taskmanager_${mode} PROPERTIES LABELS "dependency-migration" TIMEOUT 600
        ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1")
endforeach()

set(observed_source "${CMAKE_CURRENT_BINARY_DIR}/TaskManager-observed.cpp")
add_custom_command(OUTPUT "${observed_source}"
    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/instrument_taskmanager.py"
        "${PROJECT_SOURCE_DIR}/CompactStar/Core/src/TaskManager.cpp" "${observed_source}"
    DEPENDS instrument_taskmanager.py "${PROJECT_SOURCE_DIR}/CompactStar/Core/src/TaskManager.cpp")
add_executable(taskmanager_observed taskmanager_driver.cpp "${observed_source}")
target_include_directories(taskmanager_observed PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(taskmanager_observed PRIVATE CompactStar)
add_executable(unavailable_spin_diagnostics unavailable_spin.cpp)
target_link_libraries(unavailable_spin_diagnostics PRIVATE CompactStar)
add_test(NAME unavailable_spin_diagnostics COMMAND unavailable_spin_diagnostics)

add_test(NAME taskmanager_observations
    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_fixture.py"
        --exe $<TARGET_FILE:taskmanager_observed> --build-mode "${CMAKE_BUILD_TYPE}"
        --mode T1 --observe
        --eos "${COMPACTSTAR_EOS_DATA_ROOT}/DS-CMF-1-with-crust/DS(CMF)-1_with_crust.eos"
        --report "${CMAKE_CURRENT_BINARY_DIR}/observations-result.json")
set_tests_properties(taskmanager_observations PROPERTIES TIMEOUT 600
    LABELS "dependency-migration" ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1")
