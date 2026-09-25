install(FILES
    EquilibriumBaryonTangent.hpp
    OrdinaryMatterSource.hpp MovingReferenceSource.hpp ProductFate.hpp
    DirectEnergyLedger.hpp StaticZeroSpinHistory.hpp FrozenBnvValidity.hpp
    FrozenControlledBnvRunContext.hpp ControlledBnvSecularDriver.hpp BnvDiagnostics.hpp
    PassiveCheckpointOutput.hpp
    DESTINATION include/CompactStar/Physics/BNV)
set(CompactStar_Physics_BNV_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/src/MovingReferenceSource.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/DirectEnergyLedger.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/StaticZeroSpinHistory.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FrozenBnvValidity.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FrozenControlledBnvRunContext.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ControlledBnvSecularDriver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/EquilibriumBaryonTangent.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/PassiveCheckpointOutput.cpp
    PARENT_SCOPE)
