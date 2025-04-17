set(CompactStar_ChemicalHeating_headers
    SigmaOmegaRho_npemu.hpp
)

install(FILES ${CompactStar_ChemicalHeating_headers} DESTINATION include/CompactStar/ChemicalHeating)

set(CompactStar_ChemicalHeating_sources

    CompactStar/ChemicalHeating/src/SigmaOmegaRho_npemu.cpp

    PARENT_SCOPE
)