# Bounded early arithmetic checks

While Slurm allocations queued, two tiny login-node probes were checked; no
stellar solve or evolution ran on the login node. The standalone probe is only
a transcription. The second compiles the actual changed LocalThermodynamics.cpp
and TrackRFreeGasThermodynamics.cpp together with tests/linux/arithmetic_probe.cpp
and the unchanged accepted Zaki Debug archive. Compiler: authenticated TC/bin/g++,
-std=c++17 -O0 -ffp-contract=off, includes source root and Debug/Zaki/include.
The initial manual link omitted the Zaki constant archive and failed with undefined
constant symbols; adding that archive repaired the probe invocation. Production
CMake already supplies Zaki transitively. Successful build: 2.61 s, 117212 KiB.

The exact-input Decimal oracle passes the unchanged PE-V6 3e-16 criterion, twelve
inverse-density cases and forty local-EOS quantity checks. See
linux-production-probe.json for every absolute error, conditioning bound,
type width and loaded runtime path. This receipt does not replace the queued
Debug/Release full-library or representative compute-node tests.

Near-tie geometry was checked against the preserved Mac table and already valid
Linux Debug/Release CONFIND tables (dependency algorithm/archive unchanged).
The exact raw count is six; retained curves have four Mac versus six Linux points.
Both retain one simple connected vertical contour with exact endpoints and
independently bounded coordinates. See near-tie-geometry.json. This is scoped
analytic geometry evidence, not a pass of historical text identity or general
stellar/randomized topology qualification.
