# HISTORICAL ORACLES — NOT ACTIVE BUILD INPUTS:
# dependencies/include/{Zaki,Confind} and dependencies/lib are retained unchanged.
set(COMPACTSTAR_ZAKI_PREFIX "" CACHE PATH "Authenticated Zaki 2.0.1 package prefix")
set(COMPACTSTAR_CONFIND_PREFIX "" CACHE PATH "Authenticated CONFIND 2.0 package prefix")
if(NOT COMPACTSTAR_ZAKI_PREFIX OR NOT COMPACTSTAR_CONFIND_PREFIX)
    message(FATAL_ERROR "Set both COMPACTSTAR_ZAKI_PREFIX and COMPACTSTAR_CONFIND_PREFIX to qualified packages")
endif()
if(CMAKE_CONFIGURATION_TYPES OR NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Release)$")
    message(FATAL_ERROR "This dependency candidate requires separate Debug or Release builds")
endif()
execute_process(
    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_LIST_DIR}/verify_dependencies.py"
        --mode "${CMAKE_BUILD_TYPE}" --zaki "${COMPACTSTAR_ZAKI_PREFIX}"
        --confind "${COMPACTSTAR_CONFIND_PREFIX}"
        --output "${CMAKE_BINARY_DIR}/dependency-provenance.json"
    RESULT_VARIABLE dependency_status ERROR_VARIABLE dependency_error)
if(NOT dependency_status EQUAL 0)
    message(FATAL_ERROR "Dependency authentication failed: ${dependency_error}")
endif()
# Override stale *_DIR cache entries as well as disabling default search paths.
set(Zaki_DIR "${COMPACTSTAR_ZAKI_PREFIX}/lib/cmake/Zaki" CACHE PATH "Pinned Zaki config" FORCE)
set(CONFIND_DIR "${COMPACTSTAR_CONFIND_PREFIX}/lib/cmake/CONFIND" CACHE PATH "Pinned CONFIND config" FORCE)
find_package(Zaki 2.0.1 EXACT CONFIG REQUIRED PATHS "${Zaki_DIR}" NO_DEFAULT_PATH)
find_package(CONFIND 2.0.0 EXACT CONFIG REQUIRED PATHS "${CONFIND_DIR}" NO_DEFAULT_PATH)
foreach(dependency Zaki CONFIND)
    string(TOUPPER "${CMAKE_BUILD_TYPE}" dependency_mode)
    get_target_property(dependency_archive ${dependency}::${dependency} IMPORTED_LOCATION_${dependency_mode})
    if(NOT dependency_archive)
        message(FATAL_ERROR "${dependency}: requested package build mode unavailable")
    endif()
endforeach()
