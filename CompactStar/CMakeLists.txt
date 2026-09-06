add_subdirectory(Core)
add_subdirectory(Analysis)
add_subdirectory(EOS)
add_subdirectory(Physics)
add_subdirectory(Extensions)
add_subdirectory(Microphysics)

set(CompactStar_SRC_Files

    ${CompactStar_Core_sources}
    ${CompactStar_Analysis_sources}
    ${CompactStar_EOS_sources}
    ${CompactStar_BNV_sources}

    ${CompactStar_Physics_sources}
    ${CompactStar_Microphysics_sources}

    ${CompactStar_Extensions_sources}

    PARENT_SCOPE)
