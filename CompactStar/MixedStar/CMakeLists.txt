set(CompactStar_MixedStar_headers
    DarkCore_Analysis.hpp
)

install(FILES ${CompactStar_MixedStar_headers} DESTINATION include/CompactStar/MixedStar)

set(CompactStar_MixedStar_sources

    CompactStar/MixedStar/src/DarkCore_Analysis.cpp

    PARENT_SCOPE
)