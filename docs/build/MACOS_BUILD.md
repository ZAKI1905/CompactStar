# CompactStar — macOS Build Notes

> **STATUS: PHASE-1 EVIDENCE, not final user documentation.**
> This records what was *observed to work* on one machine during roadmap Phase 1A. The commands
> below are canonical for macOS development today; they are not a supported installation guide and
> the dependency-discovery step is not yet machine-independent (see *Known rough edges*).
>
> **Supported status: initial Mac-first development baseline.** `GOVERNANCE.md` §7 and
> `docs/MODERNIZATION_ROADMAP.md` Phase 1 keep Mac-first as the governed policy.
> **No claim of Linux or Windows support is made or implied.**

## Scope of Phases 1A and 1B

**Phase 1A** made a clean checkout **configure** out of source without mutating tracked files.
**Phase 1B** repaired a malformed historical merge in `RotationSolver` so the library **builds**.
Neither added tests nor changed any numerical method.

| Goal | Status |
|---|---|
| Clean checkout configures out of source | ✅ Achieved (Phase 1A) |
| Configure leaves tracked source unmodified | ✅ Achieved (Phase 1A) |
| `CompactStar` library builds | ✅ Achieved (Phase 1B) — see *RotationSolver merge repair* |
| CTest plumbing | ⬜ Deliberately not yet — next Phase-1 task |

## Canonical commands

Configure, out of source:

```bash
cmake -S . -B build -DPython3_EXECUTABLE="$(command -v python3)"
```

Diagnostic build of the library:

```bash
cmake --build build --target CompactStar -j8
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

## Known rough edges

- **Python discovery is not machine-independent.** The canonical configure command needs
  `-DPython3_EXECUTABLE` (or a project `.venv`) wherever the default Python lacks NumPy.
- **A stale generated `CompactStarConfig.h` is still tracked**, as is an apparent editor duplicate
  `CompactStar/Core/CompactStarConfig 2.h`. Both are inert for the build.
- **No warning policy or default build type yet** — remaining Phase-1 build items.
- **No test plumbing yet** — the next Phase-1 increment.

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
  format-security warning inside the vendored Zaki header. None are in `RotationSolver`; no
  warning policy exists yet (a remaining Phase-1 item).
- **No test plumbing yet** — the next Phase-1 increment.
