# Local Mac dependency migration fixtures

Frozen OLD authority: canonical CompactStar
`812463ac9ed374f64ac9cadd500066ab723d3a6c`, historical arm64 archives,
AppleClang 21, GSL `/opt/local` 2.7.1. Debug and Release are independent
same-mode authorities. All runs use TaskManager(1), DS(CMF)-1, dark mass
0.8 neutron masses, 20x20 logarithmic grid and target mass 2.01.

- `taskmanager_T1` calls actual FindCriticalCurve and FindMtotContour (which
  calls FindBtotContour) using the frozen real stellar sequence. All 13 output
  files must match hashes.
- `taskmanager_T2` generates the grid, runs the contours, Precision_Task and
  FindLimits. All 38 text outputs must match, including all three BNV_tau files.
- `taskmanager_replay` replays 1,510 exact level/intersection/GetIdx/Bisect
  records captured from historical TaskManager. The real contour-4 near-tie
  case distinguishes unpatched Zaki 2.0.0 (22 differing records in each mode)
  from historical OLD and Zaki 2.0.1. Fixture provenance is the exact e263a6e
  qualification record and `tests/fp_preservation/fixtures` in that candidate.
- `taskmanager_observations` uses a generated test-only copy of the actual
  TaskManager TU with read-only observer hooks. It compares 53,484 hex/bit
  records: raw levels, found flags, point counts/order/x/y/z, each argmax
  update/index/mass, selected point, smoothing, mass_curve[0], reimported
  curve, intersections, bounds, GetIdx and ordered Bisect points. Both OLD
  observer runs also reproduce all 13 uninstrumented OLD file hashes.

The production library contains no observer hooks. The ordinary T1/T2
executables link its unmodified TaskManager implementation. Precision_Task
input files are among the exact 13/38 governed outputs. The short temporary
working roots satisfy the historical path limit. Failed runs are retained;
successful temporary outputs are removed after comparison. Golden gzip data
uses deterministic mtime=0. Python runs tests but is never linked into C++.

The candidate's earlier qualification and the integrated replay provide the
non-vacuous unpatched-Zaki control; fixtures must never be regenerated from
NEW to make a mismatch pass.

T1, the observer and replay run from committed fixtures without external EOS
data. T2 is registered only when the DS(CMF)-1 EOS exists under the configured
COMPACTSTAR_EOS_DATA_ROOT; qualified migration runs provide and authenticate
that file and must include T2. An absent EOS is announced during configure.
