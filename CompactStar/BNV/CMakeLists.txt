set(CompactStar_BNV_headers
    BNV_Analysis.hpp
    Decay_Analysis.hpp
    BNV_Chi.hpp
    BNV_B_Chi_Photon.hpp
    BNV_B_Chi_Transition.hpp
    BNV_B_Chi_Combo.hpp
    BNV_Sequence.hpp
    BNV_B_Psi_Pion.hpp
    
)

install(FILES ${CompactStar_BNV_headers} DESTINATION include/CompactStar/BNV)

set(CompactStar_BNV_sources

    CompactStar/BNV/src/BNV_Analysis.cpp
    CompactStar/BNV/src/Decay_Analysis.cpp
    CompactStar/BNV/src/BNV_Chi.cpp
    CompactStar/BNV/src/BNV_B_Chi_Combo.cpp
    CompactStar/BNV/src/BNV_B_Chi_Photon.cpp
    CompactStar/BNV/src/BNV_B_Chi_Transition.cpp
    CompactStar/BNV/src/BNV_Sequence.cpp
    CompactStar/BNV/src/BNV_B_Psi_Pion.cpp
    
    PARENT_SCOPE
)