# ---------------------------------------------------------------------------
#  Optional historical / private program directories.
#
#  These are gitignored (see .gitignore) and are absent from a clean checkout.
#  Add each only if it is actually present, so a clean checkout configures.
# ---------------------------------------------------------------------------
set(CompactStar_optional_program_dirs
    BNV_2022
    EOS
    SpinDown_2022
    SpinDown_2023
    BNV_2023
    BNV_2024
    LightDM_2024
)

foreach(dir IN LISTS CompactStar_optional_program_dirs)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${dir}/CMakeLists.txt")
        add_subdirectory(${dir})
    else()
        message(STATUS "CompactStar: optional program directory absent, skipping: main/${dir}")
    endif()
endforeach()

# ---------------------------------------------------------------------------
#  Tracked, required. Absence or a misspelling must fail loudly, not be skipped.
# ---------------------------------------------------------------------------
add_subdirectory(Examples)
add_subdirectory(Test)
