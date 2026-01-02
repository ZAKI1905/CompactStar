add_subdirectory(Boundary)

set(CompactStar_Physics_Driver_Thermal_headers

    # HeatingFromChem.hpp
    NeutrinoCooling.hpp
    PhotonCooling.hpp
    PhotonCooling_Details.hpp
    NeutrinoCooling_Details.hpp
    NeutrinoCooling_Cache.hpp
)

install(FILES ${CompactStar_Physics_Driver_Thermal_headers} DESTINATION include/CompactStar/Physics/Driver/Thermal)

set(CompactStar_Physics_Driver_Thermal_sources

    CompactStar/Physics/Driver/Thermal/src/PhotonCooling.cpp
    CompactStar/Physics/Driver/Thermal/src/PhotonCooling_Details.cpp
    CompactStar/Physics/Driver/Thermal/src/NeutrinoCooling.cpp
    CompactStar/Physics/Driver/Thermal/src/NeutrinoCooling_Details.cpp

    # PARENT_SCOPE
)

# Append integrator sources
list(APPEND CompactStar_Physics_Driver_Thermal_sources
    ${CompactStar_Physics_Driver_Thermal_Boundary_sources}
)

# NOW export the final list
set(CompactStar_Physics_Driver_Thermal_sources
    ${CompactStar_Physics_Driver_Thermal_sources}
    PARENT_SCOPE
)
