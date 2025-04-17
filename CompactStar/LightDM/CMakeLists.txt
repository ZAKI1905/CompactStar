set(CompactStar_LightDM_headers
    LightDM_Scalar_Density.hpp
)

install(FILES ${CompactStar_LightDM_headers} DESTINATION include/CompactStar/LightDM)

set(CompactStar_LightDM_sources

    CompactStar/LightDM/src/LightDM_Scalar_Density.cpp

    PARENT_SCOPE
)