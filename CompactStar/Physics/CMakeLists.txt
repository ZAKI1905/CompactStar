add_subdirectory(Driver)
add_subdirectory(State)
add_subdirectory(Evolution)
add_subdirectory(Rotochemical)

set(CompactStar_Physics_sources

    ${CompactStar_Physics_Driver_sources}
    ${CompactStar_Physics_State_sources}
    ${CompactStar_Physics_Evolution_sources}
    ${CompactStar_Physics_Rotochemical_sources}

    PARENT_SCOPE)