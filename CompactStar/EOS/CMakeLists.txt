set(CompactStar_EOS_headers
    Baryon.hpp
    Common.hpp
    CoulombLattice.hpp
    Fermi_Gas.hpp
    Lepton.hpp
    Model.hpp
    Particle.hpp
    Polytrope.hpp
    SigmaOmega.hpp
    SigmaOmegaPar.hpp
    SigmaOmegaRho.hpp
    CompOSE_EOS.hpp
    SigmaOmegaRho_nstar.hpp
    Fermi_Gas_Many.hpp
)

install(FILES ${CompactStar_EOS_headers} DESTINATION include/CompactStar/EOS)

set(CompactStar_EOS_sources

    CompactStar/EOS/src/Baryon.cpp
    CompactStar/EOS/src/Common.cpp
    CompactStar/EOS/src/CoulombLattice.cpp
    CompactStar/EOS/src/Fermi_Gas.cpp
    CompactStar/EOS/src/Lepton.cpp
    CompactStar/EOS/src/Model.cpp
    CompactStar/EOS/src/Particle.cpp
    CompactStar/EOS/src/Polytrope.cpp
    CompactStar/EOS/src/SigmaOmega.cpp
    CompactStar/EOS/src/SigmaOmegaPar.cpp
    CompactStar/EOS/src/SigmaOmegaRho.cpp
    CompactStar/EOS/src/CompOSE_EOS.cpp
    CompactStar/EOS/src/SigmaOmegaRho_nstar.cpp
    CompactStar/EOS/src/Fermi_Gas_Many.cpp

    PARENT_SCOPE
)