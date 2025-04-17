set(CompactStar_Core_headers
    Analysis.hpp
    MixedStar.hpp
    NStar.hpp
    Prog.hpp
    Pulsar.hpp
    RotationSolver.hpp
    TaskManager.hpp
    TOVSolver_Thread.hpp
    TOVSolver.hpp
)

install(FILES ${CompactStar_Core_headers} DESTINATION include/CompactStar/Core)

set(CompactStar_Core_sources

    CompactStar/Core/src/Analysis.cpp
    CompactStar/Core/src/MixedStar.cpp
    CompactStar/Core/src/NStar.cpp
    CompactStar/Core/src/Prog.cpp
    CompactStar/Core/src/Pulsar.cpp
    CompactStar/Core/src/RotationSolver.cpp
    CompactStar/Core/src/TaskManager.cpp
    CompactStar/Core/src/TOVSolver_Thread.cpp
    CompactStar/Core/src/TOVSolver.cpp

    PARENT_SCOPE
)