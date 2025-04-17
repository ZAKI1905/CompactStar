# set(CompactStar_SRC_Files      
#     src/Baryon.cpp
#     src/Common.cpp
#     src/CoulombLattice.cpp
#     src/Lepton.cpp
#     src/Model.cpp
#     src/Particle.cpp
#     src/Polytrope.cpp
#     src/Prog.cpp
#     src/SigmaOmega.cpp
#     src/SigmaOmegaPar.cpp
#     src/SigmaOmegaRho.cpp
#     src/TOVSolver.cpp
#     src/RotationSolver.cpp
#     src/NStar.cpp
#     src/Pulsar.cpp
#     src/SigmaOmegaRho_npemu.cpp
#     src/MixedStar.cpp
#     src/Fermi_Gas.cpp
#     src/TOVSolver_Thread.cpp
#     src/Analysis.cpp
#     src/DarkCore_Analysis.cpp
#     src/TaskManager.cpp
#     src/BNV_Analysis.cpp
#     src/Decay_Analysis.cpp

#     PARENT_SCOPE
# )
add_subdirectory(Core)
add_subdirectory(EOS)
add_subdirectory(BNV)
add_subdirectory(MixedStar)
add_subdirectory(ChemicalHeating)
add_subdirectory(LightDM)

set(CompactStar_SRC_Files 

    ${CompactStar_Core_sources}
    ${CompactStar_EOS_sources}
    ${CompactStar_BNV_sources} 
    ${CompactStar_MixedStar_sources} 
    ${CompactStar_ChemicalHeating_sources}
    ${CompactStar_LightDM_sources}

    PARENT_SCOPE)