install(FILES
    OrdinaryMatterSource.hpp MovingReferenceSource.hpp ProductFate.hpp
    DirectEnergyLedger.hpp StaticZeroSpinHistory.hpp FrozenBnvValidity.hpp
    FrozenControlledBnvRunContext.hpp ControlledBnvSecularDriver.hpp BnvDiagnostics.hpp
    DESTINATION include/CompactStar/Physics/BNV)
install(FILES ../../Analysis/EquilibriumBaryonTangent.hpp DESTINATION include/CompactStar/Analysis)
set(CompactStar_Physics_BNV_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/src/MovingReferenceSource.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/DirectEnergyLedger.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/StaticZeroSpinHistory.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FrozenBnvValidity.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FrozenControlledBnvRunContext.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ControlledBnvSecularDriver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../../Analysis/src/EquilibriumBaryonTangent.cpp
    PARENT_SCOPE)
