# Phase 0 Reconnaissance — CompactStar

| Field | Value |
|---|---|
| **Report date** | 2026-08-31 |
| **Repository path** | `/home/user/CompactStar` |
| **Audited branch** | `claude/hartle-rotochemical-heating-FwAHP` |
| **Audited HEAD** | `d91c31bfcc9919915ebb0e745946f57c7f31f141` |
| **Report status** | **DESCRIPTIVE SNAPSHOT** |
| **Authority** | **NON-NORMATIVE** |
| **Scope limitation** | Conclusions apply **only** to the audited commit `d91c31b`. |

> ### ⚠ Supersession notice
>
> This report was produced against `d91c31b`, which was subsequently found **not** to be the
> latest source state. Owner commit `3639d71` ("updates", ZAKI1905, 2026-04-07) reworked
> `RotationSolver` and `MixedStar` and is **absent** from the audited branch. It is reachable
> from `master` (`9f70f14`).
>
> Later local or branch-specific work may supersede individual findings. A commit-specific
> supersession appendix is provided at the end of this document (§N). **The body below is
> preserved verbatim as written against `d91c31b`** — it has not been softened or strengthened.
> Read every §D claim against §N before acting on it.

---

## A. Repository Authentication

| Field | Value |
|---|---|
| Absolute path | `/home/user/CompactStar` |
| Current branch | `claude/hartle-rotochemical-heating-FwAHP` |
| HEAD | `d91c31b` — "Add ARCHITECTURE.md to .gitignore" |
| Working tree | Clean |
| Remote | `origin` → `https://github.com/ZAKI1905/CompactStar` |
| Default branch | `master` (there is no `main`) |
| Tracked files | 2,634 — of which **234 (8.9%) are C++ source** |

**Commit timeline:**

| Date | Activity |
|---|---|
| Dec 2025 – Jan 6 2026 | Dense human development, 46 commits by ZAKI1905 |
| Jan 6 → Apr 5 2026 | 3-month gap |
| Apr 5 2026 | `675b4a9` — 964-line **AI-authored** commit: Hartle O(Ω²) + rotochemical framework. 6 files, **no CMake edit**, no tests |
| Aug 31 2026 | `d91c31b` — `.gitignore` only |

**Build entry points:** root `CMakeLists.txt` → one STATIC library target `CompactStar`;
`main/CMakeLists.txt` adds 9 subdirectories, **7 of which do not exist** (gitignored).

**Test commands: none exist.** No `enable_testing()`, no `add_test()`, no framework, no CI.

**Governance / AI infrastructure — exhaustively searched, all absent:** `CONTRIBUTING*`,
`AGENTS.md`, `CLAUDE.md`, `GOVERNANCE.md`, `.claude/`, `.codex/`, `.github/`, ADR directories,
provenance records, validation documents, scientific-invariant documents. Present: `README.md`,
`.clang-format` (unenforced), `notes.txt` (aspirational design sketch), `docs/` (Doxyfile +
2,204 committed generated files), untracked `ARCHITECTURE.md`.

### Stop conditions triggered

1. **Code contradicts ARCHITECTURE.md** in scientifically meaningful ways (§B).
2. **Multiple active implementations; authority unestablished** — TOV integrator duplicated by
   acknowledged copy-paste, both live; MixedStar total mass has two definitions.
3. **Units ambiguous** — species columns documented as densities, used as fractions.
4. **Build would mutate source tree** — `configure_file` writes into `CMAKE_CURRENT_SOURCE_DIR`.
   **No build was run.**
5. **External data unavailable** — Zaki/Confind ship only as prebuilt macOS archives.

---

## B. ARCHITECTURE.md Authentication Matrix

| Claim / subsystem | Status | Evidence |
|---|---|---|
| StarProfile column enum + typed access | VERIFIED | `StarProfile.hpp:238-248`, `:453-476` |
| `m_version` + `EditScope` RAII deferral | VERIFIED | `StarProfile.hpp:192-218`, `:142-177` |
| NStar `prof_` canonical, legacy `ds` retained | PARTIALLY VERIFIED | `NStar.hpp:90-91` — `ds` is *commented out*, not retained |
| MixedStar still on legacy DataSet | VERIFIED | `MixedStar.hpp:269-306` |
| Superposition shooting for O(Ω²) | VERIFIED | `RotationSolver.cpp:1591-1678` |
| Shared `*_Details` pattern prevents drift | VERIFIED | `NeutrinoCooling.cpp:89,258` |
| Drivers additive via RHSAccumulator | VERIFIED | `EvolutionSystem.cpp:115,133` |
| `fast_p/e/m` are "GSL spline-backed lookups" | **CONTRADICTED** | `RotationSolver.hpp:227-242`; frozen scalars |
| `HartleResult::Omega` in `[s^-1]` | **CONTRADICTED** | `:105` vs `.cpp:1283` — geometric km⁻¹ |
| O(Ω²) is a live stage | **CONTRADICTED** | Zero call sites; `rot_solver` private |
| Second-order supports NStar **and** MixedStar | **CONTRADICTED** | `:1542-1548`, `:1702-1705` no-op stubs |
| `ProfileVersionedCache` used by 3 caches | **CONTRADICTED** | Exactly one client |
| `OmegaPoint` a live "key type" | STALE | Written only by commented-out code |
| `WV()` "universal" measure | PARTIALLY VERIFIED | Universal in Evolution; Core uses inline equivalents |
| Species columns are density columns | **CONTRADICTED** | Consumers treat as fractions |
| RotochemicalCache computes N,A,B,Z | PARTIALLY VERIFIED | Real bodies; not compiled; `A_i` units defect |
| `Rotochemical` wired into evolution | **CONTRADICTED** | `Driver/Chem/CMakeLists.txt:9-12` empty source list |
| Interpolation is spline-based | **CONTRADICTED** | `DataSet.hpp:572` → `gsl_interp_linear` |

**Net:** ARCHITECTURE.md accurately describes the *intended* design but systematically presents
scaffolding as live. It is a design document mislabeled as an architecture document.

---

## C. Current Canonical Architecture (as of `d91c31b`)

Exercised by exactly one program (`main/Test/spin_therm_evol_2_main.cpp`):

```
CompOSE EOS table
   ↓
TOVSolver::SolveToProfile → SingleStarSolveToTOVPoints   [CGS internally, G and c explicit]
   ↓  (unit conversion at NStar::Append)
NStar::prof_ : StarProfile   [geometric: km, km⁻²]  ← canonical, versioned
   ↓
RotationSolver::FindNMomInertia  → SeqPoint::I   [O(Ω) only]
   ↓
StarContext(prof_)  ──┬─→ GeometryCache (deep copy, NO version gate)
                      ├─→ DirectUrcaMask     (hand-rolled gate)
                      ├─→ HeatCapacity C(T∞) (different hand-rolled gate)
                      └─→ MassDensity g/cm³  (another gate)
   ↓
DriverContext {star, geo, envelope, thermo, cfg}   [5 raw non-owning pointers]
   ↓
StateLayout{Spin:1, Thermal:1}   ← n_eta = 0, Chem block NEVER configured
   ↓
EvolutionSystem  — plain registration-order driver loop, NO dependency graph
   ├─ MagneticDipole   dΩ/dt = −K·sign(Ω)·|Ω|ⁿ
   ├─ NeutrinoCooling  L_ν from placeholder Q₀ normalizations
   └─ PhotonCooling    L_γ / hardcoded C_eff = 1e40
   ↓
GSLIntegrator (MSBDF, rtol 1e-6 / atol 1e-10)
   ↓
TimeSeriesObserver + DiagnosticsObserver
```

Everything else — the Chem layer, second-order Hartle, rotochemical, BNV drivers,
`HeatingFromChem`, `CheckpointObserver` — is scaffolding outside this path.

---

## D. Rotation / Hartle State

> **See §N — this section is substantially superseded by `3639d71`.**

**First-order O(Ω) — PARTIAL, live, load-bearing.**

- *Center regularization:* none in RotationSolver; r=0 avoided only because TOV starts at
  `r_min = 1 cm` (`TOVSolver.hpp:559`). No series expansion.
- *Variables:* `y[0]=ω̄` [km⁻¹], `y[1]=dω̄/dr` [km⁻²].
- *Extraction:* `J = R⁴y[1]/6`, `Ω = y[0] + Ry[1]/3`, `I = J/Ω` (`:1276-1288`).
- **`init_omega_bar = 5e-3` hard-coded** (`:1240`); ω̄ equation is linear-homogeneous, so `Omega`
  and `J` carry an arbitrary normalization. Only `I` is scale-invariant.
- `HartleResult::Omega` documented `[s^-1]` but holds km⁻¹.

**Fast interpolation — mechanism (a) only.** `fast_p/e/m` are mutable member scalars assigned by
index in the outer loop and read by the RHS across all ~12 RK8PD stages. Coefficients are
piecewise-constant, evaluated at the **right endpoint**. An 8th-order integrator with 1e-10
tolerances driving a first-order-accurate source. Not thread-safe, not reentrant.

**MixedStar master-grid totals — DO NOT EXIST.** Only `mass_tot_dc`, present since the initial
import, index-summing across two independent `DataSet`s with no grid-alignment check.
`RotationSolver.cpp:1360-1380` reaches into `ds_vis[3]`, `ds_dar[4]` with raw literals.
Piecewise core/mantle zeroing is live.

**Second-order O(Ω²) — ORPHANED, UNREACHABLE, physically unverified.**

- Zero call sites; `rot_solver` private with no accessor.
- MixedStar: `f[0]=f[1]=0.0` and an empty function body.
- **The j² = e^(−2ν)(1−2m/r) factor is dropped**; `S_m` (`:1514`) uses `1/(1−2m/r)` — the
  reciprocal of what j² requires.
- Ships with `[FIX: confirm exact from textbook]` (`:1444`) and `???` (`:1518`).
- `delta_M = m0[-1]` only — J²/R³ exterior term absent.
- `dε/dp` from finite differences of the **profile**, not the EOS, with a `1.0` fallback
  mislabeled "incompressible limit".
- Central BCs literal `{0,0}` / `{0,1}` — no series expansion.
- All second-order output scaled by ~2.5e-5 relative to physical Ω.
- Superposition shooting **is** correctly implemented (`:1591-1678`).

**Validation: none.** The second-order solver has never been executed.

---

## E. Rotochemical-Heating Readiness

| Component | Status | Evidence |
|---|---|---|
| First-order Hartle ω̄, I | Partial | `RotationSolver.cpp:1234-1288` |
| Second-order Hartle m₀,p₀,ξ₀ | Partial, orphaned | `:1554-1700`; j² wrong |
| HartleResult exposure | Missing | `rot_solver` private, no accessor |
| Enclosed species numbers N_i | Partial | `RotochemicalCache.cpp:25-44` |
| A_i = ∂N_i/∂Ω²\|_εc | Partial | `:93-103` — **never divides by Ω²** |
| B_i = ∂N_i/∂ε_c\|_Ω | Partial | `:160-165` — perturbed grids × unperturbed `WV()` |
| Z_i = A_i − B_i(A_B/B_B) | **Complete in form** | `:167-171` |
| Spin-down forcing 2ΩΩ̇ | Partial | `Rotochemical.cpp:117` |
| η_npe / η_npμ definitions | Missing | Nowhere in code — comments only |
| ChemState DOF count/ordering | Missing | `n_eta = 0`; ordering undefined |
| Redshift convention for η | Missing | Absent |
| Equilibrium MUrca emissivity | Missing | `Q0_MU = 1.0e21`, self-labeled placeholder |
| **Out-of-equilibrium ΔΓ(η,T)** | **Missing** | **No function anywhere takes η** |
| F/H(ξ = η/k_BT) polynomials | Missing | Zero occurrences repo-wide |
| Direct Urca rates | Missing | `Q0_DU = 1.0e27` toy; threshold logic real |
| Superfluid suppression | Missing | `K_PBF = 0.0` |
| WeakRestoration | Missing | **0-byte file** |
| HeatingFromChem | Missing | Header only, no `.cpp`, commented out of CMake |
| Energy bookkeeping | Missing | Two different C_eff in one dT∞/dt |
| Validation vs F&R 2005 | Missing | Zero |

**Where the chain stops — four independent breaks:**

```
Ω → MagneticDipole → Ω̇                                    ✅ works
  ╫ BREAK 0 (BUILD)    Rotochemical.cpp + RotochemicalCache.cpp in NO CMake list
  ╫ BREAK 1 (DATA)     Build() has zero callers → valid=false → driver returns
  ╫ BREAK 2 (STATE)    n_eta=0; ConfigureRHS has no Chem case → AddTo would throw
η → weak rates                                             ❌ TERMINAL
  ╫ BREAK 3 (PHYSICS)  No ΔΓ(η,T) exists. No rate takes η as an argument.
  ├─ restoration       ❌ 0-byte file
  ├─ modified L_ν(η)   ❌ NeutrinoCooling never sees η
  └─ heating           ❌ no .cpp, excluded from build
T∞                                                         ✅ runs as passive cooling
```

The *geometric/structural* half of Fernández & Reisenegger is written, and the Z_i reduction —
the hardest conceptual piece — is correct in form. The *weak-interaction* half has not been
started. `Microphysics/Rates/Urca.hpp` is a mislabeled copy of `RateChannels.hpp` (guard
`..._ChemicalHeating_RateChannels_H`) with three pure-virtual methods, no implementations, no
includer; `# add_subdirectory(Rates)` is commented out.

**README.md:88 claims the formalism is "fully implemented." It is not compiled.**

---

## F. Invariants Requiring Governance (ranked)

1. **Species columns: density or fraction?** Documented fm⁻³ (`StarProfile.hpp:45`); used as
   `Y_i = n_i/n_B` (`StarContext.cpp:544-546,695`). A silent factor-of-n_B error lives here.
2. **Hartle normalization.** `init_omega_bar = 5e-3` makes Ω and J arbitrary.
3. **Fixed-ε_c vs. equilibrium-sequence derivatives.** A_i's missing Ω² division shows the
   convention is not enforced.
4. **Redshift conventions for η and μ.** Entirely undefined.
5. **Units split: CGS inside TOV, geometric outside** (`TOVSolver.cpp:1403`). Undocumented.
6. **Cache invalidation.** Four incompatible rules; GeometryCache has no gate; StarContext binds
   column pointers once and never re-binds.
7. **Interpolation is linear** (`DataSet.hpp:572`), documented as "spline".
8. **Effective heat capacity ownership.** `C_eff = 1e40` vs GR integral, same state slot.
9. **Metric convention ν, λ.** Well documented in three headers.
10. **Proper-volume measure.** `WV = 4πr²e^Λ` canonical in Evolution; Core uses inline copies.
11. **Direct-Urca threshold.** `k_Fn ≤ k_Fp + k_Fe`, **electrons only, no muon channel**.
12. **Radial grid + surface definition.** `r_min=1 cm`, `r_max=70 km`, 5-tier hardcoded step
    ladder; **two competing blanket thresholds**.
13. **Numerical tolerances.** TOV/rotation hardcoded `(1e-1, 1e-10, 1e-10)`.
14. **Baryon number.** Scalar accessor omits the ×1e54 the vectorized path applies.
15. **k_B defined three times at two precisions.**

---

## G. Validation Map and Gaps

| Subsystem | Protection today |
|---|---|
| EOS / CompOSE | **None.** Two historical eyeball checks, both **not in the build** |
| TOV solver | **None** |
| First-order Hartle | **None.** Never verified against anything |
| Second-order Hartle | **None.** Never executed |
| Rotochemical | **None.** Never compiled |
| Thermal evolution | **None.** One demo program, no assertions |
| Whole library | **Zero assertions** |

Software correctness — nothing automated. Numerical convergence — nothing. Physics regression —
nothing (`results/` files are past outputs, not baselines). Analytical validation — two unbuilt
2025-era examples.

---

## H. Modernization Classification

**PRESERVE** — StarProfile versioning/EditScope; `*_Details` shared-physics pattern;
RHSAccumulator additivity; DriverContext explicitness; Z_i reduction; superposition shooting;
`SurfaceGravity`; DUrca threshold rationale; diagnostics catalog/packet design.

**GOVERN** — all 15 invariants in §F.

**CONSOLIDATE** — duplicated TOV integrator; two NStar construction blocks; two MixedStar
total-mass definitions; proper-volume in three places; k_B ×3, km→cm ×2, km³→cm³ ×2, σ_SB ×2;
two baryon-integrand variants differing by 1e54; C_eff dual ownership.

**MODERNIZE** — large commented-out blocks; hardcoded column literals; no warning policy; no
optimization default; `configure_file` into source tree; declared-but-undefined symbols; raw
pointer lifetime discipline; `EditScope` move trap; null-deref ordering at
`NeutrinoCooling_Details.cpp:889` vs `:901`.

**VALIDATE** — TOV against known solutions; first-order I against published values; Hartle
convergence; EOS thermodynamic consistency; cooling-curve regression.

**RETIRE-CANDIDATE** *(after dependency review only)* — `CMakeLists_old.txt`; `build_xcode/`;
`docs/doxygen_output/`; `main/Test/results/` (~400 MiB incl. three ~90 MB profiler traces);
superseded `spin_evol*` chain; `SigmaOmegaRho_npemu.{hpp,cpp}`; `Microphysics/Rates/Urca.hpp`;
7 macOS duplicate files; six 0-byte headers.

**CLARIFY** — 0-byte file intent; MixedStar scope; DUrca muon omission; authoritative blanket
threshold; ARCHITECTURE.md gitignore intent; Zaki/Confind ownership.

---

## I. Recommended Initial Governance Package

Six documents, authority descending: `GOVERNANCE.md` → accepted ADRs →
`SCIENTIFIC_INVARIANTS.md` → `NUMERICAL_POLICY.md` → `VALIDATION_POLICY.md` →
current-architecture documentation.

Seed ADRs: species semantics · Hartle normalization · canonical TOV path · MixedStar scope ·
generated artifacts in VCS · dependency ownership.

---

## J. Recommended Initial AI-Skill Package

Five skills: `repo-authenticate` (read-only) · `scientific-change` · `numerical-audit`
(read-only) · `consolidate-duplication` · `governance-sync`. All share: cite governing
authority; record provenance; never hardcode reaction/campaign/species names into generic
infrastructure; fail closed rather than guess.

---

## K. Ordered Modernization Roadmap

1. Governance foundation
2. Make it buildable and reproducible
3. Validation baseline
4. Architecture consolidation
5. Rotation correctness
6. Rotochemical completion
7. BNV integration

---

## L. Recommended Immediate Next Task

Write `docs/SCIENTIFIC_INVARIANTS.md` and ADR-0001 (species columns: density or fraction?).
It cannot be resolved from repository evidence alone — code and documentation disagree, and only
the owner can say which was intended.

---

## M. Evidence Appendix

**Core:** `StarProfile.hpp` (:142-218, :238-248, :453-476, :746-752, :913) · `NStar.hpp`
(:90-91, :164, :254, :276-310) · `NStar.cpp` (:23, :90-343, :358-632, :102-123, :157, :211-227,
:274-283, :675-738, :823-842, :1053-1068, :1077) · `MixedStar.hpp:269-306` · `MixedStar.cpp`
(:127-152, :472-493, :776-780) · `TOVSolver.hpp` (:456, :468, :559, :564, :629) · `TOVSolver.cpp`
(:1204-1208, :1365-1406, :1708-1712, :1741, :1751, :1839-1858, :1971, :2500-2673, :2553/:2574/:2597)
· `RotationSolver.hpp` (:56-68, :102-124, :105, :130-152, :227-242, :339/:343/:372) ·
`RotationSolver.cpp` (:358-410, :459-494, :1234-1288, :1240, :1259-1272, :1360-1380, :1428,
:1444/:1518, :1449-1538, :1487, :1502-1515, :1542-1548, :1554-1700, :1591-1678, :1696, :1702-1705)
· `SeqPoint.hpp:200` · `Pulsar.cpp:204,212`

**Evolution:** `ProfileCache.hpp:31,109-190` · `StarContext.hpp:215-228` · `StarContext.cpp`
(:55-107, :177-196, :241-281, :412-591, :505-506, :544-546, :636-642, :662, :695, :712-717,
:725-738, :757-817) · `GeometryCache.cpp:118-133,162-237` · `EvolutionSystem.cpp:101-197` ·
`RHSAccumulator.hpp:62,126-131` · `EvolutionConfig.hpp:162-192,295` · `GSLIntegrator.cpp:44-61,461-466`
· `RunBuilder.cpp:39,126-141` · `UnitContract.cpp:32-33` · `RotochemicalCache.{hpp,cpp}`

**Drivers:** `IDriver.hpp:76,89,97,102-103` · `MagneticDipole.cpp:70-154` ·
`NeutrinoCooling.cpp:75-116,129-138` · `NeutrinoCooling_Details.cpp` (:94, :100-102, :121-199,
:214, :217-792, :889/:901, :988-1035) · `PhotonCooling.hpp:229` ·
`PhotonCooling_Details.cpp:50-53,213-216,309-324` · `Rotochemical.cpp:48-124` ·
`HeatingFromChem.hpp:29,74-77` · `TbDefinition.cpp:69-139` · `EnvelopePotekhin1997.cpp:85-120` ·
`EnvelopePotekhin2003.cpp:18-46` · `SurfaceGravity.cpp:38-207` · 0-byte: `AccretionTorque.hpp`,
`BNVSpinTorque.hpp`, `BNVSource.hpp`, `WeakRestoration.hpp`, `Coupling.hpp`, `CheckpointObserver.hpp`

**State/Micro:** `State.hpp:116` · `ChemState.hpp:18,46-57,94-97,133` ·
`ThermalState.hpp:52-66,166` · `SpinState.hpp:175` · `Microphysics/Rates/Urca.hpp:27-121`

**Build/docs:** `CMakeLists.txt` (:10-15, :27-32, :54-60, :70-96, :84-87, :136, :141-147) ·
`CMakeLists_old.txt:18-19,68,76` · `main/CMakeLists.txt:1-9` · `Driver/Chem/CMakeLists.txt:1-12` ·
`Driver/Thermal/CMakeLists.txt:5` · `dependencies/.../DataSet.hpp:572` ·
`dependencies/.../Constants.hpp:111-128,254` · `docs/Doxyfile:502,911,1002,1011` ·
`README.md:31,88-104,171-196,333-394` · `notes.txt`

---

## N. Supersession Appendix (added 2026-08-31, post-reconciliation)

Owner commit `3639d71` ("updates", 2026-04-07) is **not** an ancestor of the audited `d91c31b`.
It reworked `RotationSolver.hpp/.cpp` and `MixedStar.hpp/.cpp` (+755 / −888 lines).
`RotationSolver.cpp` shrank from **1707 → 1265 lines**.

All eighteen probe identifiers are **present in `9f70f14`** and **absent from `d91c31b`**:
`kR_EPS_KM`, `SafeR0`, `SetFastProfilePtrs_`, `EvalFastPEM_`, `SetFastMixedPtrs_`,
`EvalFastMixedPEM_`, `fast_r_`, `fast_r_mix_`, `fast_k_`, `fast_k_mix_`, `r_master_dc`,
`pre_tot_dc`, `eps_tot_dc`, `totals_ready_`, `GetRadius_Master`, `GetPress_Total_Master`,
`GetEps_Total_Master`, `GetMass_Total_Master`.

### Findings SUPERSEDED (must be re-audited against `9f70f14`)

| §D finding (against `d91c31b`) | Superseding evidence in `9f70f14` |
|---|---|
| "No r=0 regularization in RotationSolver" | `kR_EPS_KM = 1e-6` km + `SafeR0()` (`RotationSolver.cpp:31-37`), `r_safe` guards at `:244,260,276,292`, `SafeR0` used at `:375,933` |
| "`fast_p/e/m` frozen at grid right-endpoint; first-order accurate" | `EvalFastPEM_(r, p, e, m)` called from `ODE_N_Fast` at `:318,328` with the **actual GSL RHS radius** |
| "No cached bracket indices" | `fast_k_`, `fast_k_mix_` left-bracket hunt (`:92,105,147-158,197-206`) |
| "MixedStar master-grid totals DO NOT EXIST" | `r_master_dc`, `pre_tot_dc`, `eps_tot_dc`, `totals_ready_` (`MixedStar.hpp:307-339`), built in `SurfaceIsReached()` |
| "No total pressure or energy-density column" | `GetPress_Total_Master()`, `GetEps_Total_Master()`, `GetMass_Total_Master()` (`MixedStar.hpp:471-503`) |
| "Piecewise core/mantle zeroing is live" | `FindMixedMomInertia` (`:891`) integrates over the master grid and guards on `HasTotalMasterProfiles()`; the piecewise branch at `:800` is commented out |
| "~940/1707 lines commented-out legacy" | File reduced to 1265 lines |

### Findings that PERSIST in `9f70f14` (re-verified)

- `init_omega_bar` still present (13 occurrences) — **Hartle normalization invariant unresolved**.
- `SolveHartle2_N` still present — second-order reachability **not** re-verified; treat as open.
- **Rotochemical build gap persists**: `Physics/Driver/Chem/CMakeLists.txt` still declares an
  empty `..._sources` list (`Rotochemical.hpp` appears only in the *headers* install list), and
  `RotochemicalCache.cpp` is still **absent** from `Physics/Evolution/CMakeLists.txt`.
  BREAK 0 of §E is unchanged.
- `ARCHITECTURE.md` is **not tracked** in `9f70f14`.

### Not re-audited

§C, §E (beyond BREAK 0), §F, §G, §H were established against `d91c31b` and have **not** been
re-verified against `9f70f14`. Findings outside `RotationSolver`/`MixedStar` are unlikely to have
changed — `3639d71` touched only those four files — but this has not been proven file-by-file.

**A full re-audit of `RotationSolver` and `MixedStar` against `9f70f14` is required before any
Phase-4 rotation work. The superseded §D claims must not be cited as current.**
