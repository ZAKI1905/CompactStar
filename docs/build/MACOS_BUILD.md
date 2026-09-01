# CompactStar — macOS Build Notes

> **STATUS: PHASE-1 EVIDENCE, not final user documentation.**
> This records what was *observed to work* on one machine during roadmap Phase 1A. The commands
> below are canonical for macOS development today; they are not a supported installation guide and
> the dependency-discovery step is not yet machine-independent (see *Known rough edges*).
>
> **Supported status: initial Mac-first development baseline.** `GOVERNANCE.md` §7 and
> `docs/MODERNIZATION_ROADMAP.md` Phase 1 keep Mac-first as the governed policy.
> **No claim of Linux or Windows support is made or implied.**

## Scope of Phases 1A, 1B and 1C

**Phase 1A** made a clean checkout **configure** out of source without mutating tracked files.
**Phase 1B** repaired a malformed historical merge in `RotationSolver` so the library **builds**.
**Phase 1C** added the CTest mechanism and one build/link smoke test.
None of the three changed a numerical method or added a scientific reference value.

| Goal | Status |
|---|---|
| Clean checkout configures out of source | ✅ Achieved (Phase 1A) |
| Configure leaves tracked source unmodified | ✅ Achieved (Phase 1A) |
| `CompactStar` library builds | ✅ Achieved (Phase 1B) — see *RotationSolver merge repair* |
| CTest plumbing + smoke test | ✅ Achieved (Phase 1C) — see *Automated tests* |
| Warning policy / default build type | ⬜ **Outstanding** — remaining Phase-1 item |

## Canonical commands

Configure, out of source:

```bash
cmake -S . -B build -DPython3_EXECUTABLE="$(command -v python3)"
```

Build everything (library, manual demo programs, and the test executable):

```bash
cmake --build build -j8
```

Build only the library:

```bash
cmake --build build --target CompactStar -j8
```

Run the automated tests:

```bash
ctest --test-dir build --output-on-failure
```

Configure without tests (library only):

```bash
cmake -S . -B build-no-tests -DPython3_EXECUTABLE="$(command -v python3)" -DBUILD_TESTING=OFF
```

`build/` is already gitignored. It is disposable — delete and re-create it freely.

### Why `-DPython3_EXECUTABLE` is currently required

The root `CMakeLists.txt` requires `find_package(Python3 COMPONENTS Interpreter Development NumPy
REQUIRED)` (`CMakeLists.txt:60`). On a machine where the first Python CMake finds lacks NumPy,
configure fails at that line **before reaching anything else**. That is what happened here:
CMake selected Homebrew Python 3.14.2, which has no NumPy, while NumPy 2.3.1 was installed under
miniforge Python 3.12.10.

Two supported ways to resolve it, both already anticipated by the root `CMakeLists.txt`:

1. **Point CMake at a Python that has NumPy** — the `-DPython3_EXECUTABLE=` override documented
   in-file at `CMakeLists.txt:37-38`. Substitute the interpreter that actually has NumPy; on the
   reference machine that was `/Users/keeper/miniforge3/bin/python3`.
2. **Create a project virtualenv** with NumPy at `.venv/` in the repository root. `CMakeLists.txt:41-49`
   prefers `.venv/bin/python` automatically when present, needing no `-D` flag. No `.venv` existed
   on the reference machine, so this path was not exercised.

Making Python selection robust without a machine-specific flag is **not** Phase-1A work and is
recorded as future build work.

## Reference environment — observed 2026-08-31

| Component | Value |
|---|---|
| macOS | 15.6.1 (build 24G90) |
| Architecture | `arm64` (Apple Silicon, T6000) |
| Kernel | Darwin 24.6.0 |
| CMake | 4.2.1 (`/opt/homebrew/bin/cmake`) |
| C / C++ compiler | AppleClang 17.0.0.17000604 (`/usr/bin/clang`, `/usr/bin/clang++`), forced on APPLE by `CMakeLists.txt:10-15` |
| Xcode | `/Applications/Xcode.app/Contents/Developer` |
| C++ standard | C++17 (`CMakeLists.txt:27-29`) |
| GSL | 2.7 — headers resolved at `/opt/local/include` |
| OpenMP | 5.1, via `-Xclang -fopenmp` |
| Python3 | 3.12.10 — `/Users/keeper/miniforge3/bin/python3` |
| NumPy | 2.3.1 |
| Zaki archive | `dependencies/lib/Zaki/Darwin/arm64/libZaki.a` (2,478,624 bytes) |
| Confind archive | `dependencies/lib/Confind/Darwin/arm64/libConfind.a` (703,208 bytes) |

The vendored archive paths are composed from `CMAKE_SYSTEM_NAME` and `CMAKE_HOST_SYSTEM_PROCESSOR`
(`CMakeLists.txt:70-71`) and hard-fail if absent (`:73-79`). Only `Darwin/{arm64,x86_64}` are
shipped, which is why non-macOS configuration is impossible today.

## What Phase 1A changed

### 1. Optional program directories no longer break a clean checkout

`main/CMakeLists.txt` previously called `add_subdirectory()` unconditionally for nine directories.
Seven of them are **explicitly gitignored** — `main/BNV_2022/`, `main/EOS/`, `main/SpinDown_2022/`,
`main/SpinDown_2023/`, `main/BNV_2023/`, `main/BNV_2024/`, `main/LightDM_2024/` (`.gitignore:7-18`)
— and are therefore absent from every clean checkout, producing seven hard CMake errors.

They are now added **only if** the directory carries a `CMakeLists.txt`; otherwise a `STATUS`
message is emitted and configure continues:

```
-- CompactStar: optional program directory absent, skipping: main/BNV_2022
```

`main/Examples` and `main/Test` are tracked, are not gitignored, and remain **unconditional**.
A missing or misspelled required directory therefore still fails loudly — verified by temporarily
misspelling `Examples`, which produced
`CMake Error at main/CMakeLists.txt:28 (add_subdirectory)`.

No directory was fabricated and no placeholder CMake file was added.

### 2. The generated config header no longer lands in the source tree

`configure_file` previously wrote into `${CMAKE_CURRENT_SOURCE_DIR}/CompactStar/Core/CompactStarConfig.h`.
Because `CompactStar_RELEASE_DATE` comes from `string(TIMESTAMP ...)` (`CMakeLists.txt:22`),
**every configure rewrote a tracked file** — even a configure that then failed. Observed directly:

```
 M CompactStar/Core/CompactStarConfig.h
-#define CompactStar_RELEASE_DATE "01, 05, 2026"
+#define CompactStar_RELEASE_DATE "08, 31, 2026"
```

The header is now generated into a binary-tree include root that mirrors the include spelling used
by the sources (`Prog.hpp:253`, `Core/src/Banner.cpp:12`):

```
build/generated/include/CompactStar/Core/CompactStarConfig.h
```

and that root is placed on `CompactStar`'s include path **ahead of** the source root.
No macro name, value semantics, or version definition changed — only the output location.

`CompactStarConfig.h` is not a member of `CompactStar_Core_headers`
(`CompactStar/Core/CMakeLists.txt:1-17`), so **install behavior is unaffected** by the move.

### 3. The `build/` ignore rule is anchored

`.gitignore` carried an unanchored `build/`, which in git matches a directory of that name **at
any depth** — so it silently swallowed `docs/build/` as well, and this very file would not have
been committable. It is now `/build/`, which still ignores the top-level out-of-source build
directory and nothing else.

Blast radius checked before changing it: the only directories named `build` in the tree are the
disposable top-level `./build` and `./docs/build`, and **no tracked file lives under any `build/`
path segment**, so anchoring altered no existing tracked state. (`build_xcode/` is a different
name and was never matched by this rule.)

## Which `CompactStarConfig.h` the compiler uses

A **stale generated copy is still tracked** at `CompactStar/Core/CompactStarConfig.h`, so two
copies coexist. This is deliberate: removing a tracked generated artifact is a separate
artifact/architecture decision (`GOVERNANCE.md` §2, *generated artifact* class) and was out of
scope for Phase 1A.

The two copies are distinguishable — the tracked copy reads `"01, 05, 2026"`, the freshly generated
one `"08, 31, 2026"` — which makes resolution directly testable. **The binary-tree copy wins**,
confirmed two independent ways.

Include order on the real compile line (from `compile_commands.json`, `Core/src/Banner.cpp`):

```
-I.../build/generated/include      <-- generated root, first
-I...                              <-- source root, second
-I.../dependencies/include
```

`clang -H` on that exact command reports the header actually opened:

```
. .../build/generated/include/CompactStar/Core/CompactStarConfig.h
```

Preprocessing a probe TU with the same flags yields `"08, 31, 2026"`; reversing the two `-I` flags
yields `"01, 05, 2026"`, confirming the ordering is what decides it.

**Consequence:** the stale tracked copy is inert for the build but is a trap for anyone reading it
as current. Retiring it is recorded as future work.

## Verification performed

- Fresh configure from a deleted/recreated `build/` succeeds.
- `git status --short` after configure shows **no** modification to any tracked source file.
- Re-configuring an existing `build/` is idempotent and likewise leaves the tree clean.
- A misspelled required directory still fails configure.
- The generated header exists only under `build/generated/include/...`.

## RotationSolver merge repair (Phase 1B)

Phase 1A left the library unbuildable: 11 translation units compiled, then
`CompactStar/Core/src/RotationSolver.cpp` failed with 20 errors. The cause was **not** a coding
mistake but a **malformed historical merge**.

### Origin

Merge `9f70f14` ("Merge branch 'master'", 2026-04-07, three minutes after `3639d71`) combined:

- **`3639d71`** — owner work introducing profile-backed interpolation at the actual GSL RHS radius,
  cached bracket indices, `SafeR0`/finite-radius handling, and the single-pass MixedStar
  master-grid moment-of-inertia integration;
- **`e60e656`** — the PR #1 second-order Hartle / rotochemical candidate lineage
  (`GOVERNANCE.md` §5: UNVERIFIED SCIENTIFIC CANDIDATE).

The merge resolved three hunks the wrong way, in each case taking the candidate side:

1. **Header taken wholesale.** `RotationSolver.hpp` at `9f70f14` is the candidate header verbatim
   except that `fast_p`/`fast_e`/`fast_m` were commented out. All ten owner interpolation members
   (`fast_r_`, `fast_p_`, `fast_e_`, `fast_m_`, `fast_k_`, `fast_r_mix_`, `fast_p_tot_`,
   `fast_e_tot_`, `fast_m_tot_`, `fast_k_mix_`) and four method declarations
   (`SetFastProfilePtrs_`, `SetFastMixedPtrs_`, `EvalFastPEM_`, `EvalFastMixedPEM_`) were dropped —
   while the `.cpp` kept defining and using them.
2. **`FindNMomInertia` setup deleted.** The declarations of `R` and `r`, the `i0`/`kR_EPS_KM`
   start-radius scan, the `P`/`E`/`M` fetches, and the `SetFastProfilePtrs_(*R, *P, *E, *M)` call
   were lost, yet the body still referenced `R` and `r`. The candidate's `HartleResult` population
   was also duplicated verbatim.
3. **A dead function was resurrected.** Seven active lines from `e60e656:1293-1299` were spliced
   over the owner's fully commented-out obsolete implementation (`3639d71:734-740`), making
   `void RotationSolver::FindMixedMomInertia()` and its `{` active again while the matching `}`
   stayed commented. The file therefore ended at **brace depth 1**, illegally nesting the last
   ~370 lines — including the real master-grid `FindMixedMomInertia` — inside a phantom body.

### Repair

Reconstructed the intended **union** of the two parents rather than choosing either:

- Restored the owner's ten members, four method declarations, and `FindNMomInertia` setup from
  `3639d71`; removed the duplicated `HartleResult` block.
- Re-commented the two spliced lines, closing the brace imbalance.
- Restored `fast_p`/`fast_e`/`fast_m` as declarations, because the second-order candidate code
  (`ODE_Hartle2_N_Fast`, `SolveHartle2_N`) reads and writes them. No first-order code touches them.
- Commented out the five genuinely obsolete first-order scalars (`fast_p_v`, `fast_p_d`,
  `fast_e_v`, `fast_e_d`, `fast_m_tot`), matching `3639d71` — verified to have no active use.

**Second-order candidate code was not altered and remains unratified and unreachable.**
`init_omega_bar` was not touched: INV-07 defers its normalization to Phase 4.

### Numerical scope

Engineering class. The first-order O(Ω) path is equivalent to `3639d71`:

- `cpp.3639d71:5-667` is **byte-identical** to the corresponding pre-repair block (663 lines),
  covering `kR_EPS_KM`, `SafeR0`, both `SetFast*Ptrs_`, both `EvalFast*PEM_`, `ODE_N_Fast`,
  `ODE_Mixed_Fast`, every `GetHartle*Coeff*`, and `Solve_Mixed`. Not one arithmetic expression
  differs, before or after the repair.
- The repaired `FindNMomInertia` diffs against `3639d71` with **additions only** — the candidate's
  `stored_omega_bar_`/`stored_domega_bar_` storage and the `HartleResult` copy. These write to new
  members and never feed back into `y[]`, `r`, or the GSL driver.
- The repaired `FindMixedMomInertia` is identical to `3639d71` modulo indentation.

### Build result

```
[100%] Linking CXX static library libCompactStar.a
[100%] Built target CompactStar
```

70 translation units compile (was 11), 0 errors, `libCompactStar.a` ≈ 10 MB. Verified by a clean
rebuild from a deleted `build/`. No warning appears in `RotationSolver`, and there are no
`unused-*` warnings anywhere.

## Automated tests (Phase 1C)

### Layout and mechanism

```
tests/
├── CMakeLists.txt
└── smoke/
    └── compactstar_library_smoke.cpp
```

`tests/` is the **canonical automated-test root**. The top-level `CMakeLists.txt` calls
`include(CTest)` — which defines the standard `BUILD_TESTING` option (default `ON`) and calls
`enable_testing()` — and then adds `tests/` only when `BUILD_TESTING` is on. No custom test switch
was invented, no third-party framework was added, and no test target is installed.

**`main/Test/` is not the automated suite.** It holds manual demo and debugging programs
(`tov_debug_main`, `spin_therm_evol_2_main`, `compose_thermo_main`, …). None of them is registered
with CTest, none was renamed, and none of their output is treated as a baseline.

### What the smoke test proves

`compactstar_library_smoke` constructs `CompactStar::Core::Prog` with a name and calls
`GetName()`. `Prog` is the library's non-scientific base class; `Prog::GetName()` is defined
**out of line** in `CompactStar/Core/src/Prog.cpp`, so the call forces the linker to pull a real
object file out of the archive. An empty `main()` linked against a static library can pull nothing
at all, which would make such a test vacuous — this one is not:

```
$ nm -C build/tests/compactstar_library_smoke | grep 'Prog::GetName'
0000000100003630 T CompactStar::Core::Prog::GetName() const
```

`T` means *defined in the text section*, not an unresolved import. The binary carries 16
`CompactStar::` symbols and is 233 KB against 17 KB for an empty `main()`.

So the test demonstrates that: public headers compile; an executable resolves a genuine
`libCompactStar.a` symbol; the archive and its transitive dependencies (GSL, OpenMP, Python3,
Zaki, Confind) link; and the binary starts and exits zero under CTest.

`Prog` was chosen over `RotationSolver` because it is semantically neutral — nothing about it
suggests a physics result. Its constructor performs only member initialisation, it reads no files,
needs no EOS table or data directory, opens no network connection, and mutates nothing.
Constructing through the naming constructor means `GetName()` takes its non-logging path.
`ShowBannerOnce()` was also considered and rejected: it calls `EnableClearScreen()`, which would
wipe the terminal of anyone running `ctest` interactively.

### What the smoke test explicitly does NOT prove

**It is infrastructure validation, not the CompactStar scientific validation baseline.**

It asserts no stellar mass, no moment of inertia, no cooling luminosity, no EOS value, no Hartle
result, no reaction rate, no thermal result, and no reference number of any kind — and none may be
added to it. It says nothing about whether any physics in the library is correct. Scientific
validation is roadmap **Phase 2A** (independent verification of `C_⋆(T∞)` under `GOVERNANCE.md`
§3.1) and **Phase 2B** (the regression baseline), governed separately.

### Observed result

```
$ ctest --test-dir build -N
  Test #1: compactstar_library_smoke
Total Tests: 1

$ ctest --test-dir build --output-on-failure
1/1 Test #1: compactstar_library_smoke ........   Passed    0.33 sec
100% tests passed, 0 tests failed out of 1
```

Exactly one test is discovered — the manual programs under `main/Test/` do not appear.

With `-DBUILD_TESTING=OFF`: configure succeeds, no `tests/` subdirectory is generated,
`ctest -N` reports `Total Tests: 0`, and `cmake --build build-no-tests --target CompactStar`
still produces `libCompactStar.a`.

Note that `cmake --build build --target CompactStar` does **not** build the test executable, and
CTest does not build it either — a full `cmake --build build` is required first. This was verified
rather than assumed.

`.gitignore` gained `/build-*/` so the documented `build-no-tests/` directory stays untracked. It
does not match the tracked `build_xcode/` (hyphen versus underscore).

## Known rough edges (unchanged, pre-existing)

- **Python discovery is not machine-independent** — see above.
- **A stale generated `CompactStarConfig.h` is still tracked**, plus an editor duplicate
  `CompactStar/Core/CompactStarConfig 2.h`. Both inert for the build.
- **Residual indentation artifact.** The ~370 lines after the former splice point are still
  tab-indented from the historical reformat. This is cosmetic, compiles correctly, and was left
  alone deliberately — reformatting is out of scope.
- **17 pre-existing warnings** now visible because the whole library compiles: variable-length
  arrays (`SigmaOmegaRho.cpp:580`, `SigmaOmegaRho_nstar.cpp:585`), a missing `override`
  (`EnvelopePotekhin2003.hpp:58`), a C++20 `using enum` in a C++17 build (`Tags.hpp:72`), and a
  format-security warning inside the vendored Zaki header. None are in `RotationSolver`. Adopting
  a warning policy and a default build type is **the one outstanding Phase-1 item**; `-Werror` is
  deliberately not enabled.
- **Test coverage is one smoke test.** It proves the mechanism, nothing more; real validation is
  Phase 2A/2B.
