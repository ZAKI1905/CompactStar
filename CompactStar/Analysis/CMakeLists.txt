install(FILES ParticleNumberResponse.hpp DESTINATION include/CompactStar/Analysis)
install(FILES ChemicalResponse.hpp DESTINATION include/CompactStar/Analysis)
set(CompactStar_Analysis_sources
    CompactStar/Analysis/src/ParticleNumberResponse.cpp
    CompactStar/Analysis/src/ChemicalResponse.cpp
    PARENT_SCOPE)
