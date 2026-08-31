# CompactStar — macOS Build Notes

> **STATUS: PHASE-1 EVIDENCE, not final user documentation.**
> This records what was *observed to work* on one machine during roadmap Phase 1A. The commands
> below are canonical for macOS development today; they are not a supported installation guide and
> the dependency-discovery step is not yet machine-independent (see *Known rough edges*).
>
> **Supported status: initial Mac-first development baseline.** `GOVERNANCE.md` §7 and
> `docs/MODERNIZATION_ROADMAP.md` Phase 1 keep Mac-first as the governed policy.
> **No claim of Linux or Windows support is made or implied.**

## Scope of Phase 1A

Phase 1A made a **clean checkout configure** out of source without mutating tracked files.
It did **not** attempt a full working build, add tests, or touch scientific source.

| Goal | Status |
|---|---|
| Clean checkout configures out of source | ✅ Achieved |
| Configure leaves tracked source unmodified | ✅ Achieved |
| `CompactStar` library builds | ❌ **Blocked** by pre-existing source defect — see *Known next blocker* |
| CTest plumbing | ⬜ Deliberately not in this increment — next Phase-1 task |

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

## Known next blocker — library build fails

`cmake --build build --target CompactStar` was run as a **diagnostic**. It is **not** a Phase-1A
exit criterion and it **failed**. Reported separately from configure success, which passed.

- 11 translation units compiled successfully; the build stops in
  `CompactStar/Core/src/RotationSolver.cpp` with **20 errors** (`-ferror-limit` reached).
- The failure is a **header/implementation mismatch**: the `.cpp` defines a "fast profile pointer"
  fast path whose members and methods are **not declared** in `CompactStar/Core/RotationSolver.hpp`.
  Undeclared in the header but used in the `.cpp`: `fast_r_`, `fast_m_`, `fast_k_`, `fast_r_mix_`,
  `fast_p_tot_`, `fast_e_tot_`, `fast_m_tot_`, `fast_k_mix_`. Three out-of-line definitions match no
  declaration: `SetFastProfilePtrs_` (`:83`), `SetFastMixedPtrs_` (`:96`), `EvalFastPEM_` (`:114`).
  Related type errors at `:102-104` assign `const Zaki::Vector::DataColumn *` to `double`.
- **This is pre-existing and unrelated to the Phase-1A changes.** `RotationSolver.cpp` and
  `RotationSolver.hpp` are unmodified at `1677caa`, and `RotationSolver.cpp` contains **zero**
  references to `CompactStarConfig` — the only thing this increment relocated.
- It plausibly originates in owner commit `3639d71`, which reworked `RotationSolver` and which
  `docs/SCIENTIFIC_INVARIANTS.md` INV-05 already flags ⚠ as requiring re-audit. **Repairing it is
  not authorized by Phase 1A** and must not be done as a drive-by fix: `RotationSolver` is on the
  O(Ω) path (INV-07) and any change to it needs its governing class established first.

One non-fatal warning also appears: `Physics/State/Tags.hpp:72` uses `using enum StateTag;`,
a C++20 construct in a C++17 build (`-Wc++20-extensions`). AppleClang accepts it as an extension.
