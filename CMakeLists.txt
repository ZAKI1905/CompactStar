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
# Default development build type
# -----------------------------
# Historically no CMAKE_BUILD_TYPE was set, so builds were unoptimized with no
# debug information. Debug is the closest standard configuration to that prior
# behavior while adding symbols. Release/RelWithDebInfo are NOT made the default:
# introducing optimization implicitly, before any numerical baseline exists,
# would change results with nothing to detect the change against.
#
# A user-supplied -DCMAKE_BUILD_TYPE is never overridden. Multi-configuration
# generators (Xcode, Visual Studio) select per-build and are left alone.
#
# Debug and Release are NOT asserted to be numerically equivalent. Any future
# scientific baseline must state the build configuration it was produced under.
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING
        "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        Debug Release RelWithDebInfo MinSizeRel)
    message(STATUS "CompactStar: no CMAKE_BUILD_TYPE given, defaulting to Debug")
endif()

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
# Generated configuration header goes in the BINARY tree, never the source tree.
# The include spelling used by the sources is "CompactStar/Core/CompactStarConfig.h"
# (Prog.hpp:253, Core/src/Banner.cpp:12), so mirror that layout under a generated
# include root and put that root ahead of the source root on the include path.
set(CompactStar_GENERATED_INCLUDE_DIR "${CMAKE_BINARY_DIR}/generated/include")

configure_file(
    CompactStarConfig.h.in
    ${CompactStar_GENERATED_INCLUDE_DIR}/CompactStar/Core/CompactStarConfig.h
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

# CompactStar-owned include roots — ordinary (non-system), so our own headers
# are subject to the warning policy below.
target_include_directories(CompactStar
    PUBLIC
    # Must precede the source root: a stale generated copy is still tracked at
    # CompactStar/Core/CompactStarConfig.h, and this ordering makes the
    # binary-tree header win. See docs/build/MACOS_BUILD.md.
    ${CompactStar_GENERATED_INCLUDE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# External and vendored dependency headers — SYSTEM, so the compiler suppresses
# diagnostics originating inside them (-isystem under AppleClang/GCC). This keeps
# CompactStar's warning inventory a statement about CompactStar's own code.
# It does NOT fix those dependencies: warnings in Zaki/Confind belong in their own
# source repositories. See docs/build/MACOS_BUILD.md.
target_include_directories(CompactStar
    SYSTEM PUBLIC
    ${GSL_INCLUDE_DIRS}
    ${DEP_INC_DIR}
)

# -----------------------------
# Warning policy (CompactStar-owned code only)
# -----------------------------
# Baseline signal floor, not a final quality policy. PRIVATE on purpose: these
# flags must not propagate to consumers of the installed library, nor to the
# manual programs under main/ or the tests under tests/.
#
# -Werror is deliberately NOT enabled. Warnings surfaced by this policy are
# inventoried in docs/build/MACOS_BUILD.md and repaired under their own change
# class -- some touch scientific source and are not engineering-class fixes.
target_compile_options(CompactStar
    PRIVATE
    $<$<COMPILE_LANG_AND_ID:CXX,AppleClang,Clang,GNU>:-Wall;-Wextra>
    $<$<COMPILE_LANG_AND_ID:C,AppleClang,Clang,GNU>:-Wall;-Wextra>
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
# Testing (standard CTest semantics)
# -----------------------------
# include(CTest) defines the BUILD_TESTING option (default ON) and calls
# enable_testing() for us. tests/ is the canonical automated-test root;
# main/Test/ holds manual demo programs and is deliberately NOT registered.
# Configuring with -DBUILD_TESTING=OFF skips tests entirely; the CompactStar
# library still configures and builds, and no test target is ever installed.
include(CTest)

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()

# -----------------------------
# Packaging
# -----------------------------
include(InstallRequiredSystemLibraries)
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
set(CPACK_PACKAGE_VERSION_MAJOR "${CompactStar_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${CompactStar_VERSION_MINOR}")
include(CPack)