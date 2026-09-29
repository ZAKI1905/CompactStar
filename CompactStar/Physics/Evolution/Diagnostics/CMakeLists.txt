set(CompactStar_Physics_Evolution_Diagnostics_headers
	DiagnosticsJson.hpp
	DiagnosticPacket.hpp
	UnitContract.hpp
	DiagnosticCatalog.hpp
	DiagnosticsCatalogJson
)

install(FILES ${CompactStar_Physics_Evolution_Diagnostics_headers} DESTINATION include/CompactStar/Physics/Evolution/Diagnostics)

set(CompactStar_Physics_Evolution_Diagnostics_sources
	CompactStar/Physics/Evolution/Diagnostics/src/DiagnosticsJson.cpp
	CompactStar/Physics/Evolution/Diagnostics/src/DiagnosticPacket.cpp
	CompactStar/Physics/Evolution/Diagnostics/src/UnitContract.cpp
	CompactStar/Physics/Evolution/Diagnostics/src/DiagnosticCatalog.cpp
	CompactStar/Physics/Evolution/Diagnostics/src/DiagnosticsCatalogJson.cpp
	PARENT_SCOPE
)