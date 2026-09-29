set(CompactStar_Physics_Evolution_Observer_headers
	IObserver.hpp
	DiagnosticsObserver.hpp
	TimeSeriesObserver.hpp
)

install(FILES ${CompactStar_Physics_Evolution_Observer_headers} DESTINATION include/CompactStar/Physics/Evolution/Observer)

set(CompactStar_Physics_Evolution_Observer_sources
	CompactStar/Physics/Evolution/Observers/src/IObserver.cpp
	CompactStar/Physics/Evolution/Observers/src/DiagnosticsObserver.cpp
	CompactStar/Physics/Evolution/Observers/src/TimeSeriesObserver.cpp

	PARENT_SCOPE
)