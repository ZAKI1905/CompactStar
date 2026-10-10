# CompactStar manual handoff — campus-laptop instructions

**Purpose:** source/input delivery and preparation for a separately authorized Linux
qualification. This package does not authorize a campaign, establish Linux scientific
equivalence, or resolve the numerical gaps in `PHASE6_BNV_CAMPAIGN_READINESS.md:1`.
The Mac task made no EKU connection. No remote address, SSH configuration, key,
credential store, installed Mac library, or disposable build is supplied.

## Verify and unpack

Keep the ZIP and its adjacent `.sha256` receipt together. Compare the receipt to the
SHA-256 provided by the owner; a manifest shipped with a payload is an integrity
check, not an independent signature. On the campus laptop, in their directory:

```sh
sha256sum -c CompactStar-EKU-handoff-20261010.zip.sha256
unzip CompactStar-EKU-handoff-20261010.zip
cd CompactStar-EKU-handoff-20261010
python3 verify_transfer.py --root .
python3 import_sources.py --bundles bundles --output "$PWD/../CompactStar-EKU-sources"
```

On macOS use `shasum -a 256 -c` for the first line. On Windows PowerShell use
`Get-FileHash -Algorithm SHA256` and compare its digest to the receipt; use a
Linux/WSL filesystem for the Git import and subsequent Linux builds. Preserve
case-sensitive paths and executable modes. Python 3, ZIP tools, and Git with
`index-pack --promisor` support are required. The importer was tested with Git
2.52.0 and Python 3.12.10 on this Mac; campus versions remain to be checked.

The importer verifies each bundle SHA-256, original refs, exact commit set,
all supplied checkout file hashes, omitted-file absence, clean sparse status,
and `git fsck --full`. Its only configured remote is a deliberately nonexistent
**local** path named `OMITTED_OBJECTS_NOT_SUPPLIED`. It performs no network fetch.
The source directories must be new; the importer refuses an existing destination.

## Bundle semantics and identities

These are **filtered v3 Git bundles**, with all original reachable commit and tree
objects and source blobs. Commits and refs are not rewritten. They are not complete
binary-history backups. `git clone` rejects filtered bundles on the tested Git;
use the supplied importer. The v3 `filter=blob:none` capability marks the pack as
promisor; additional source blobs allow an offline source checkout. The format and
promisor requirement are specified in <https://git-scm.com/docs/gitformat-bundle>.

`bundles/{CompactStar,Zaki,CONFIND}.json` names every included ref, checkout SHA,
commit-count/hash, supplied file, and omitted blob/path. The CompactStar checkout
is the documentation/transfer branch. Production/CMake/test bytes equal readiness
`5433fc2da71ba380264f490b87926d850ff9f375`; production numerical bytes also equal
accepted breakpoint candidate `0c3f9269c884161c9d89163660a47b141a93ace6`.
Master remains `9a5c1eca40758e4c39f9f80db7324e12ad23e505`.

Zaki is checked out at accepted production
`e263a6e180c5c417198e7778bd21fc9c0a32dc33`. Its master ref is
`7aca56cfdf3875cf2a8c0acea4acdac44f69c9f5`, which adds only the acceptance document.
`refs/handoff/accepted-production` is a transfer-only ref, not a new upstream tag.
CONFIND is checked out at accepted master
`b0cbd510fd3fd0c772fa50499cd749287cb39e7b`.
Sources: `evidence/source-build-tests.json:1` and the bundle manifests.

Full commit logs/ancestry and included source history are inspectable. Commands
that request omitted binary/generated blobs fail; do not disable sparse checkout,
pretend omitted files are deletions, or commit their removal. A broader historical
checkout may need a new sparse selection of already supplied source blobs.
The omitted paths include inactive Mac archives, generated Xcode/Doxygen material,
old profiling/results data, and duplicate review archives. No active source,
CMake or test path was omitted. Original full repositories remain on the Mac.

## Scientific inputs and reference evidence

`scientific/input-map.json` binds relocated paths to original locations, exact
SHA-256 identities, and historical authorities. The six Phase-6 inputs are the
radial certificate, profiles, thermal tables, frozen BNV certificate, coefficients,
and entry manifest. Both CompOSE sets and the three obs32 input files are also
included. File bytes, including embedded original paths, are preserved. Directory
hashes use sorted relative names, 8-byte big-endian name lengths, and raw file
SHA-256 digests, matching `tests/bnv/run_breakpoint_qualification.py:9`.

`verify_transfer.py` checks both the complete transfer manifest and all scientific
file/tree identities. It does not rerun depletion qualification. The historical
BNV certificate is authenticated evidence, not a physical rate/model selection or
a globally accepted trajectory. Original Mac packages were recovered and checked,
but are intentionally excluded from this source transfer.

For qualification, use a separate run manifest that maps the six logical keys to
these relocated paths. Do not edit certificate/coefficients/thermal tables or their
embedded provenance to accommodate a host. The current runner hardcodes Mac paths;
a separate relocation adapter is still required before Linux numerical execution
(`tests/bnv/run_breakpoint_qualification.py:52`). Keep exact input hashes fail-closed.

The readiness report, global accuracy predeclaration, ADR-0018/0019, acceptance
records, baselines and compressed numerical evidence are present in CompactStar's
source bundle. The independent Fable 5.1 report itself was not supplied to this
task; its owner-reported disposition must not be replaced by the older obs32 report.

## Current Linux requirements and blockers

This inventory uses the named production sources, not the older preflight's
superseded vendored-dependency descriptions.

| Requirement | Current basis and disposition |
|---|---|
| C/C++ compiler and standard library | CompactStar C++17, extensions OFF; Zaki C++17, extensions ON; CONFIND C++17. Native Linux compiler/linker/libc/libm identities and behavior are UNQUALIFIED. `CMakeLists.txt:1`, `:62`; Zaki `CMakeLists.txt:8`; CONFIND `CMakeLists.txt:7`. |
| C++17 conformance | `CompactStar/Physics/State/Tags.hpp:72` uses C++20 `using enum`; AppleClang accepts an extension. A strict C++17 build has a real portability defect. Make a separate mechanical candidate repair or explicit language-policy decision; do not silently increase the global standard. |
| CMake/build tools | At least CMake 3.22 for all three projects (CompactStar 3.18, Zaki 3.16, CONFIND 3.22); Make or Ninja, Git, Python interpreter, tar/gzip, a writable isolated prefix. |
| GSL | Exact upstream 2.7.1 source supplied, SHA-256 `dcb0fbd43048832b757ff9942691a8dd70026d5da0ff85601e52687f6deeb34b`. GNU HTTPS download matches original MacPorts distfile and installed Portfile checksum (`evidence/gsl-source.json:1`). Linux library still needs build/authentication. No GPG signature verification claimed. |
| GSL Mac provenance caveat | Portfile is supplied for inspection; its `dynamic_lookup-11.patch` bytes were not recovered in the installed registry. GNU source identity is established, not a reconstruction of every MacPorts build action. Linux GSL ODE/math sources and build policy require review. |
| zlib and threads | zlib headers/library via `ZLIB::ZLIB`; standard Threads. No current C++ Python/NumPy, ROOT or OpenMP requirement in the accepted dependency graph. `Zaki/CMakeLists.txt:12`, `:87`; CONFIND `CMakeLists.txt:4`; CompactStar `CMakeLists.txt:91`. |
| Python | Interpreter for build manifest checks/tests. NumPy and mpmath for one of the 13 tests (`chemical_curved_gc9.py:7`, `chemical_curved_reference.py:7`); SciPy for readiness diagnostics (`readiness_summary.py:5`) and wider qualification. Linux versions/wheels/native dependencies are not supplied or qualified. |
| Packages/lock | Separate Debug and Release Zaki/CONFIND prefixes; exact CONFIG versions 2.0.1/2.0.0. Existing verifier is hardwired to Mac lock (`cmake/verify_dependencies.py:29`). A Linux package cannot pass that lock. See next section. |
| Library layout | Force `CMAKE_INSTALL_LIBDIR=lib` because active target checks expect `prefix/lib`, not `lib64` (`cmake/ExternalDependencies.cmake:21`). |
| Qualification harness | Mac paths, `.dylib` suffixes, `/usr/bin/clang++`, and bit-exact platform metadata exist in producer/comparator paths. Adapt separately without weakening Mac oracles (`tests/rotochemical/fresh_context.py:53`, `tests/bnv/adr0017_production_verify.py:136`; older EKU preflight §12). |
| Runtime/reproducibility | Record compiler flags, FP mode, C/C++ runtimes, GSL/zlib/header/archive/config hashes, CPU/OS, linkage, environment and input hashes. Preserve existing `-ffp-contract=off`, HistoricalMath `-fno-builtin`/`-fno-lto`. No fast-math, LTO or native-CPU tuning by default. Tests with fixed temporary paths run serially. |
| Cluster | Module availability, filesystem, scheduler/partition/account, resource limits and instruction set must be re-inventoried on the campus side. Earlier GCC/CMake/GSL observations are historical, not current facts. No Slurm script is activated by this handoff. |

The old preflight remains evidence, but its Python C-API/OpenMP/vendored-CONFIND
migration concerns were superseded by ADR-0019 and the canonical acceptance.
Its numerical/platform gates also predate breakpoint Amendment 1 and the readiness
predeclaration. Reconcile the qualification plan explicitly; never waive failed
BA12/BA12R contraction/floor gates or promote bounded-window evidence to global PASS.

## Platform lock design (proposal, not implemented)

Keep `cmake/dependency-lock.json` byte-identical at
`e1307a042180034ad028f2ff07dbc7fd80752494ff4b2ff61cc8ca0eec2089c5`.
Do not overwrite it with Linux hashes, forge Mac provenance, disable the verifier,
or reuse a Release prefix for Debug.

A future dependency/build candidate should add explicit Linux-only lock selection,
for example `COMPACTSTAR_LINUX_DEPENDENCY_LOCK`, while Darwin continues to use the
fixed original lock and rejects Linux selection. The Linux lock must record OS,
architecture, compiler/build/toolchain identity, exact source SHAs, mode, clean
source status, archive/header/config digests and authenticated dependency closure
(GSL and zlib included). Authenticate bytes **before** package CMake code executes;
retain exact versions, absolute prefixes, no default package fallback, and imported
archive/include-path checks. Reject absent locks, wrong platform/mode, stale hashes,
unknown dependencies, or provenance mismatches. Isolate prefixes by toolchain and
mode, e.g. `<root>/toolchains/<id>/Debug/{Zaki,CONFIND}`.

Initial Linux packages/lock are candidates. Build success authenticates a build,
not scientific equivalence. Mac-lock preservation and wrong-lock/mode/hash refusal
need explicit tests, followed by dependency tests and the governed platform
qualification. Do not create a Linux PASS merely by hashing a new archive.

## Staged Linux bootstrap (instructions only)

Proceed only in the separately authorized campus/Linux task. First verify the
handoff, re-read authority documents, inventory the host, and predeclare the
qualification matrix. Do not copy Mac binaries or modify master. Required package
manager access, zlib development files and Python wheels are listed in
`unavailable-required-inputs.json`. No automatic download is hidden in this guide.

Set absolute `SRC` to the imported source directory, `HANDOFF` to the unpacked
package, `TC` to a new toolchain root, and `CC`/`CXX` to discovered, recorded Linux
compiler executables. Create separate build/install/log directories. For example,
a **candidate** static GSL bootstrap is:

```sh
mkdir -p "$TC/source" "$TC/build/gsl" "$TC/install"
tar -xzf "$HANDOFF/upstream/gsl-2.7.1.tar.gz" -C "$TC/source"
cd "$TC/build/gsl"
"$TC/source/gsl-2.7.1/configure" --prefix="$TC/install/gsl-2.7.1" \
  --disable-shared --enable-static CC="$CC" CFLAGS='-O2 -g -ffp-contract=off'
make -j2
make check
make install
```

Capture configure/cache/logs, test results, archive and installed-header hashes.
These proposed GSL flags are not a claim of equivalence to MacPorts. Audit them
before numerical qualification. Do not fall back to the historical cluster GSL 2.6.

With authenticated local zlib, set `ZLIB_INCLUDE` and `ZLIB_LIBRARY` explicitly.
For each mode (Debug first; Release is a separate candidate):

```sh
MODE=Debug
GSL="$TC/install/gsl-2.7.1"
ZAKI="$TC/install/$MODE/Zaki"
CONFIND="$TC/install/$MODE/CONFIND"
cmake -S "$SRC/Zaki" -B "$TC/build/Zaki-$MODE" \
  -DCMAKE_BUILD_TYPE="$MODE" -DCMAKE_CXX_COMPILER="$CXX" \
  -DCMAKE_INSTALL_PREFIX="$ZAKI" -DCMAKE_INSTALL_LIBDIR=lib \
  -DGSL_ROOT_DIR="$GSL" -DGSL_INCLUDE_DIR="$GSL/include" \
  -DGSL_LIBRARY="$GSL/lib/libgsl.a" -DGSL_CBLAS_LIBRARY="$GSL/lib/libgslcblas.a" \
  -DZLIB_INCLUDE_DIR="$ZLIB_INCLUDE" -DZLIB_LIBRARY="$ZLIB_LIBRARY" \
  -DZAKI_BUILD_TESTS=ON -DZAKI_BUILD_EXAMPLES=OFF \
  -DCMAKE_FIND_USE_PACKAGE_REGISTRY=OFF -DCMAKE_FIND_USE_SYSTEM_PACKAGE_REGISTRY=OFF
cmake --build "$TC/build/Zaki-$MODE" --parallel 2
ctest --test-dir "$TC/build/Zaki-$MODE" --output-on-failure -j1
cmake --install "$TC/build/Zaki-$MODE"
cmake -S "$SRC/CONFIND" -B "$TC/build/CONFIND-$MODE" \
  -DCMAKE_BUILD_TYPE="$MODE" -DCMAKE_CXX_COMPILER="$CXX" \
  -DCMAKE_INSTALL_PREFIX="$CONFIND" -DCMAKE_INSTALL_LIBDIR=lib \
  -DZaki_DIR="$ZAKI/lib/cmake/Zaki" -DCMAKE_PREFIX_PATH="$GSL" \
  -DGSL_INCLUDE_DIR="$GSL/include" -DGSL_LIBRARY="$GSL/lib/libgsl.a" \
  -DGSL_CBLAS_LIBRARY="$GSL/lib/libgslcblas.a" \
  -DZLIB_INCLUDE_DIR="$ZLIB_INCLUDE" -DZLIB_LIBRARY="$ZLIB_LIBRARY" \
  -DCONFIND_BUILD_TESTS=ON \
  -DCMAKE_FIND_USE_PACKAGE_REGISTRY=OFF -DCMAKE_FIND_USE_SYSTEM_PACKAGE_REGISTRY=OFF
cmake --build "$TC/build/CONFIND-$MODE" --parallel 2
ctest --test-dir "$TC/build/CONFIND-$MODE" --output-on-failure -j1
cmake --install "$TC/build/CONFIND-$MODE"
```

Record every test outcome; platform-bound Mac reference tests may fail and must
remain failures until separately adjudicated. Check Zaki's generated Version.hpp
reports `e263a6e...` and clean source. Authenticate the resulting ELF archives,
headers/configs and complete link closure into a **candidate Linux lock**.

**Stop before CompactStar configuration until the explicit Linux lock selector and
C++17 defect have been addressed on a separately reviewable candidate branch.**
Current CompactStar is expected to refuse Linux artifacts against its Mac lock.
After those repairs, use explicit `COMPACTSTAR_ZAKI_PREFIX`,
`COMPACTSTAR_CONFIND_PREFIX`, GSL libraries/includes, Python executable, and
`COMPACTSTAR_EOS_DATA_ROOT="$HANDOFF/scientific/compose"`; build mode must match the
prefixes. Record source/compiler/build/input/executable identities and run the
13 tests in `evidence/source-build-tests.json` serially, followed by the agreed
broader qualification matrix. Keep Debug/Release and source/control/tolerance
qualifications distinct. No full 500 kyr campaign follows automatically.

## Recommended campus-agent first action

Verify ZIP and input hashes, import the three source repositories, and return an
inventory of the actual Linux host/toolchain/offline dependencies plus a concrete
plan for the Linux lock selector, C++17 repair, and relocated qualification
harness. Obtain the separate platform-qualification authorization required by the
governance records before executing numerical windows or scheduler jobs. The
owner's Mac preparation request alone is not cluster-execution authorization.
