cmake_minimum_required(VERSION 3.10)

# ============================================================
# CompactStar (Top-level) — deterministic macOS toolchain + Python
# ============================================================

# -----------------------------
# Toolchain selection (MUST be before project())
# -----------------------------
if(APPLE)
    # Prefer AppleClang on macOS for consistent libc++ runtime behavior.
    # If you want to override, do it on the cmake command line and remove FORCE.
    set(CMAKE_C_COMPILER "/usr/bin/clang" CACHE FILEPATH "C compiler" FORCE)
    set(CMAKE_CXX_COMPILER "/usr/bin/clang++" CACHE FILEPATH "C++ compiler" FORCE)
endif()

# -----------------------------
# Project
# -----------------------------
project(CompactStar VERSION 2.0 DESCRIPTION "Compact Star" LANGUAGES C CXX)

string(TIMESTAMP CompactStar_RELEASE_DATE "%m, %d, %Y")

# -----------------------------
# Language standard
# -----------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Minimal global compile options you were using
add_compile_options(-pthread)

# -----------------------------
# Prefer project venv Python if present (prevents Homebrew python3.14 surprises)
# -----------------------------
# If you want a different Python, override at configure time:
# cmake .. -DPython3_EXECUTABLE=/path/to/python
#
# This block only sets Python3_EXECUTABLE if you did NOT already provide one.
if(NOT DEFINED Python3_EXECUTABLE)
    if(EXISTS "${CMAKE_SOURCE_DIR}/.venv/bin/python")
        set(Python3_EXECUTABLE "${CMAKE_SOURCE_DIR}/.venv/bin/python"
            CACHE FILEPATH "Project virtualenv Python" FORCE)
    elseif(EXISTS "${CMAKE_SOURCE_DIR}/.venv/bin/python3")
        set(Python3_EXECUTABLE "${CMAKE_SOURCE_DIR}/.venv/bin/python3"
            CACHE FILEPATH "Project virtualenv Python3" FORCE)
    endif()
endif()

# -----------------------------
# Dependencies
# -----------------------------
find_package(GSL REQUIRED)

# OpenMP: do NOT inject flags globally; link the imported target instead
find_package(OpenMP REQUIRED)

# Python / NumPy
find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)
message(STATUS "Using Python3 executable: ${Python3_EXECUTABLE}")
message(STATUS "Python3 version: ${Python3_VERSION}")

# -----------------------------
# Local dependency paths
# -----------------------------
set(DEP_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/lib)
set(DEP_INC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/include)

set(ZAKI_LIB ${DEP_LIB_DIR}/Zaki/${CMAKE_SYSTEM_NAME}/${CMAKE_HOST_SYSTEM_PROCESSOR}/libZaki.a)
set(CONFIND_LIB ${DEP_LIB_DIR}/Confind/${CMAKE_SYSTEM_NAME}/${CMAKE_HOST_SYSTEM_PROCESSOR}/libConfind.a)

if(NOT EXISTS "${ZAKI_LIB}")
    message(FATAL_ERROR "Missing Zaki static library: ${ZAKI_LIB}")
endif()

if(NOT EXISTS "${CONFIND_LIB}")
    message(FATAL_ERROR "Missing Confind static library: ${CONFIND_LIB}")
endif()

# -----------------------------
# Configure headers
# -----------------------------
configure_file(
    CompactStarConfig.h.in
    ${CMAKE_CURRENT_SOURCE_DIR}/CompactStar/Core/CompactStarConfig.h
)

# -----------------------------
# Subdirectories / sources
# -----------------------------
add_subdirectory(CompactStar)

# Build the CompactStar static library from sources exported by the subdir
add_library(CompactStar STATIC ${CompactStar_SRC_Files})
target_compile_features(CompactStar PUBLIC cxx_std_17)
set_target_properties(CompactStar PROPERTIES
    VERSION ${PROJECT_VERSION}
    POSITION_INDEPENDENT_CODE ON
)

target_include_directories(CompactStar
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${GSL_INCLUDE_DIRS}
    ${DEP_INC_DIR}
)

target_link_libraries(CompactStar
    PUBLIC
    ${ZAKI_LIB}
    ${CONFIND_LIB}

    GSL::gsl
    GSL::gslcblas

    Python3::Python
    Python3::Module
    Python3::NumPy

    OpenMP::OpenMP_CXX
)

# -----------------------------
# Install rules
# -----------------------------
install(TARGETS CompactStar DESTINATION lib)

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/dependencies/include"
    DESTINATION dependencies
)
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/dependencies/lib"
    DESTINATION dependencies
)

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/License" DESTINATION Docs)

# -----------------------------
# Optional docs
# -----------------------------
option(GEN_DOCS "Generate documents using doxygen." OFF)

if(GEN_DOCS)
    execute_process(
        COMMAND doxygen ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()

# -----------------------------
# Executables
# -----------------------------
add_subdirectory(main)

# -----------------------------
# Packaging
# -----------------------------
include(InstallRequiredSystemLibraries)
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
set(CPACK_PACKAGE_VERSION_MAJOR "${CompactStar_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${CompactStar_VERSION_MINOR}")
include(CPack)