# Passive-cooling regression baseline — Phase 2B-1

> **STATUS: PASSIVE COOLING BASELINE ESTABLISHED.**
> `GOVERNANCE.md` §3.1 **condition 7 is SATISFIED**. The golden baseline is committed at
> `tests/baselines/passive_cooling_cmf_1p6_debug.tsv` and the regression runs under CTest as
> `passive_cooling_regression` whenever `COMPACTSTAR_EOS_DATA_ROOT` is supplied.
>
> **This is a regression baseline, not a validation of the placeholder neutrino emissivity
> normalizations.** `Q0_DU = 1.0e27` and `Q0_MU = 1.0e21` remain self-labelled placeholders; a
> justified change to them is expected to move every value here and must be re-baselined.
>
> **All passive-cooling outputs predating Phase 2A-3 remain superseded as validation references.**

| Field | Value |
|---|---|
| **Purpose** | Freeze `C_⋆(T∞) dT∞/dt = −L_ν,∞ − L_γ,∞` as a regression reference |
| **Governance role** | `GOVERNANCE.md` §3.1 condition 7 — **discharged** |
| **Source commit** | repair `a2c1f62`, baseline this commit (parent `800f245`) |
| **Build / toolchain** | `Debug`, AppleClang 17.0.0.17000604, CMake 4.2.1, GSL 2.7.1, macOS 15.6.1 arm64 |
| **Harness** | `tests/thermal/passive_cooling_regression.cpp` |
| **Runtime** | ~11 s |

## 1. §3.1 condition-7 chain

1. **Authorizing ADR** — ADR-0002.
2. **Correction made (Phase 2A-3)** — `−L_γ/1e40` replaced by `−L_γ/C_⋆(T∞)`.
3. **Independent evidence** — ADR-0002 §V1, **V1 VERIFIED** (`HEAT_CAPACITY_V1.md`).
4. **Baseline** — **established here.** Condition 7 satisfied.

## 2. The two blockers, and how they were repaired

Phase 2B-1 was blocked twice. Both defects were pre-existing and unrelated to ADR-0002.

### 2.1 The default stepper could not run

`Config::stepper` defaulted to `MSBDF`, an implicit BDF method, while
`GSLIntegrator::Integrate` passed `sys.jacobian = nullptr`.

**Official GSL contract** (GSL 2.7.1, the include the build actually uses at
`/opt/local/include`): `gsl_odeiv2.h:48` — *"Some methods require the jacobian function"*; and
`gsl_odeiv2.h:69` defines

```c
#define GSL_ODEIV_JA_EVAL(S,t,y,dfdy,dfdt)  (*((S)->jacobian))(t,y,dfdy,dfdt,(S)->params)
```

which dereferences `(S)->jacobian` with **no null guard**. So GSL calls address `0x0` and the
process dies with `SIGSEGV`. **This is not a macOS quirk** — it is the documented API contract.

**Repaired (commit `a2c1f62`)**: the default is now `RKF45`, chosen from the convergence study
below; `MSBDF` is rejected before any GSL allocation with a message naming MSBDF and the missing
Jacobian, and is **never silently substituted**; `spin_therm_evol_2_main.cpp` now sets the stepper
explicitly. A Jacobian was deliberately **not** implemented — that is a larger numerical feature.

### 2.2 A thermal-only system could not run

`EvolutionSystem::operator()` bound `GetThermal()`/`GetSpin()` in a block whose only consumer was
a commented-out `Z_LOG_INFO`. Because `GetSpin()` throws on a missing tag, a registered `Spin`
block was **mandatory for every evolution**. Nothing read the bindings, so removing them cannot
change any RHS value. **Repaired** in the same commit; guarded by `evolution_stepper_contract`.

## 3. Numerical study

### 3.1 Explicit-stepper comparison (Configuration A, spin+thermal)

Full `100 yr → 1 Myr` trajectory on the authenticated 1.6 M☉ CMF star. No historical result file
was used as expected data.

| Stepper | `T∞(10⁶ yr)` @ rtol 1e-6 | nominal → rtol 1e-7 | runtime |
|---|---|---|---|
| **RKF45** | 4.66966060e5 | **6.35e-8** | 11.08 s |
| RKCK | 4.66966138e5 | 9.48e-8 | 10.99 s |
| RK8PD | 4.66966058e5 | 3.13e-7 | 10.99 s |

Cross-stepper at rtol 1e-6: RK8PD vs RKCK **1.81e-7**, RK8PD vs RKF45 **2.76e-7**.

All three complete, are stable under tightening, emit no GSL errors, and converge to a common
trajectory — four orders of magnitude inside the `1e-3` gate, at indistinguishable runtime.

**Selected: RKF45.** Being numerically indistinguishable, the lower-complexity general-purpose
method is preferred over claiming unsupported superiority for a higher-order scheme; it also had
the smallest nominal-to-tightest difference. **This is a non-stiff choice for the present passive
system only** — see §3.6.

### 3.2 Repeatability

Three runs from clean process starts: **bit-identical** (max rel `ΔT∞` = 0).

### 3.3 Tolerance tightening (RKF45, thermal-only)

| vs nominal rtol 1e-6 | max rel `ΔT∞` | max rel `ΔL` |
|---|---|---|
| rtol 3e-7, atol 1e-11 | 1.19e-7 | 7.15e-7 |
| rtol 1e-7, atol 1e-12 | 6.35e-8 | 3.81e-7 |

Well inside the `1e-3` gate: the nominal configuration is converged.

### 3.4 Cadence sensitivity — and one fragile cadence

| samples/decade | result | `T∞(10⁶ yr)` |
|---|---|---|
| 50 | ok | 4.66966511e5 |
| 75 | ok | 4.66965107e5 |
| 100 | ok | 4.66964486e5 |
| **150 (canonical)** | ok | **4.66966060e5** |
| 200 | ok | 4.66965978e5 |
| 250 | ok | 4.66965893e5 |
| **300** | **FAILS** | state becomes non-finite at t ≈ 302 yr |
| 400 | ok | 4.66966005e5 |

Across the seven working cadences the spread in `T∞(10⁶ yr)` is **4.3e-6** — cadence does *not*
materially change the physical trajectory, which is the property §15 requires.

**The 300-samples/decade failure is recorded as a defect, not explained away.** It is isolated and
non-monotonic — 400 (finer) succeeds — so it is not "finer is worse". Its cause was not diagnosed.
See §3.6 for why the early transient makes this regime fragile.

### 3.5 Continuous vs segmented integration

The Phase-2B-1 harness restarted the integrator at every checkpoint, discarding adaptive step
history. Measured difference against one continuous `Integrate(t0,t1)`:

**max rel `ΔT∞` = 4.03e-3, max rel `ΔL` = 2.39e-2.**

That is ~400× the state tolerance — restarting is *not* the same numerical procedure. **The
continuous trajectory is the baseline authority**, and the canonical harness now performs a single
continuous integration with checkpoints captured by a test-side observer.

### 3.6 The early transient is stiff — a recorded limitation

`T∞` falls from `1.0e9 K` at 100 yr to `1.21e7 K` at 302 yr — a factor of **83 in 200 years** —
driven by `L_ν = 4.45e45 erg/s` at t₀ under the placeholder `Q0_DU`. The system is effectively
**stiff** in that window, which is exactly where an explicit method is least comfortable and where
the 300-samples/decade fragility appears.

This also suggests why `MSBDF` was chosen as the original default: someone wanted a stiff solver.
It was never usable because no Jacobian exists. **RKF45 is adequate for the present system on the
evidence above, but a stiff method with a real Jacobian is likely the right long-term answer**,
particularly once chemical heating adds further timescales. Not decided here.

## 4. Baseline

### Structural fingerprint

| Property | Value |
|---|---|
| Requested mass | `1.6 M☉` |
| **Achieved mass** | **`1.599975834 M☉`** |
| **Radius** | **`13.4683 km`** |
| Central energy density | `7.3125e14` |
| Radial points | `2635` |
| Profile version | `5` |

### Spin decoupling — demonstrated, not assumed

With the thermal-only repair in place, Configuration A (spin + thermal, `MagneticDipole` active)
and Configuration B (thermal only) were run under identical star, EOS, thermo, envelope, drivers,
initial state, interval, tolerances and stepper:

**max rel `ΔT∞` = 0, max rel `ΔL` = 0 — bit-identical.**

Spin is exactly decoupled from passive thermal evolution, so **Configuration B is the canonical
baseline**: it isolates the quantity under test and carries no `Ω`.

### Golden checkpoints

Nine log-spaced epochs, single continuous integration, captured on the LogTime sampling grid.
Full precision is in the committed TSV; abridged here.

| t (yr) | `T∞` (K) | `C_⋆` (erg/K) | `L_ν` (erg/s) | `L_γ` (erg/s) | `L_ν/L_γ` | regime |
|---|---|---|---|---|---|---|
| 1.00e2 | 1.0000e9 | 2.4274e39 | 4.4527e45 | 2.2866e36 | 1.95e9 | ν-dominated |
| 3.02e2 | 1.2087e7 | 2.9340e37 | 1.3882e34 | 5.2284e31 | 2.66e2 | ν-dominated |
| 1.00e3 | 8.3075e6 | 2.0160e37 | 1.4637e33 | 2.1102e31 | 6.94e1 | ν-dominated |
| 3.02e3 | 6.1681e6 | 1.4972e37 | 2.4520e32 | 1.0265e31 | 2.39e1 | ν-dominated |
| 1.00e4 | 4.4966e6 | 1.0919e37 | 3.6805e31 | 4.7771e30 | 7.70 | ν-dominated |
| 3.02e4 | 3.3105e6 | 8.0375e36 | 5.8607e30 | 2.2768e30 | 2.57 | comparable |
| 1.00e5 | 2.2421e6 | 5.4421e36 | 5.6568e29 | 8.8672e29 | 0.638 | γ-dominated |
| 3.02e5 | 1.3290e6 | 3.2269e36 | 2.4537e28 | 2.5012e29 | 0.098 | γ-dominated |
| 1.00e6 | 4.6697e5 | 1.1337e36 | 4.6167e25 | 1.9901e28 | 2.32e-3 | γ-dominated |

The **ν → γ crossover falls between 3×10⁴ and 10⁵ yr** in this model. Because `Q0_DU` and `Q0_MU`
are placeholders, that transition is a property of the current regression model and **not an
astrophysical prediction**.

The TSV also carries `L_nu_DU`, `L_nu_MU`, `Tsurf_K`, `Tb_K` and `dLnTinf_dt_1_s`, so a future
failure can be localised to structure, heat capacity, neutrino, photon, envelope or integrator.

### Energy-accounting identity

At every checkpoint, independently from driver diagnostics:

```
d ln T∞/dt  ==  −(L_ν,∞ + L_γ,∞) / (C_⋆(T∞) · T∞)
```

**Maximum residual: 2.12e-16 — machine precision.** Checked at every checkpoint, and a failure
aborts the run as an accounting defect rather than reporting baseline drift.

The **ADR-0002 Pattern-A invariant** is checked too: the `C_⋆` used by `PhotonCooling` and by
`NeutrinoCooling` must be identical to `1e-12` relative. It holds at every checkpoint.

## 5. Tolerance policy

Derived from measurement, fixed before the golden values were accepted, and applied uniformly —
not hand-tuned per cell.

| Input | Value |
|---|---|
| Repeat-run variation | 0 (bit-identical) |
| Nominal vs tighter integration | 1.19e-7 |
| Cadence 150 → 75 | 2.04e-6 (state), 1.22e-5 (luminosity) |
| Floating-point floor | 1e-14 |
| **max** | **2.04e-6** state, **1.22e-5** luminosity |
| Safety factor (fixed in advance) | ×5 |

```
state       tolerance = 1e-5     (T∞, C_⋆, Tsurf, Tb, d ln T∞/dt)
luminosity  tolerance = 1e-4     (L_ν, L_ν,DU, L_ν,MU, L_γ)
energy identity                  1e-10
Pattern-A C_⋆ agreement          1e-12
```

Both are far tighter than percent level, so the baseline is numerically mature. The
continuous-vs-segmented difference (4.03e-3) and the spin difference (0) are **reported findings,
deliberately not tolerance inputs** — the canonical procedure is fixed, so its alternatives must
not widen the gate.

## 6. Regression-detector proof

A modest, test-side perturbation — `PhotonCooling::Options::global_scale` `1.00 → 1.01` — is
caught decisively:

```
checkpoints exceeding tolerance : 16
max rel ΔT∞                     : 1.04e-2      (1000x the 1e-5 state tolerance)
```

The perturbation lives behind the harness's `--detector` mode and is **never** the canonical
configuration, so it cannot be left committed by accident. Canonical comparison passes; perturbed
comparison fails.

## 7. Historical January-2026 provenance — **HISTORICAL STEPPER UNKNOWN**

`main/Test/results/spin_therm_evol_2/spin_therm_evol_2_main.log` records a run on 2026-01-06
02:07 reaching `OnFinish(ok=true, t_final=3.1556952e13)` — exactly 1 Myr — which the then-default
`MSBDF + null Jacobian` could not have done.

Bounded history search establishes:

| Fact | Evidence |
|---|---|
| The only stepper ever set explicitly in that program was **RKF45** | `0ce9008` (2025-12-18) added `cfg.stepper = StepperType::RKF45;` |
| It was already **commented out** before the run | at `8586741^` the line reads `// cfg.stepper = ...RKF45;` |
| It was **deleted** on 2026-01-05, the day before the run | `8586741` |
| `Config`'s default was **never** anything but `MSBDF` | only `d07ac2a` (2025-11-29) ever touched that line |

So at run time the committed source had **no** stepper override and a default that cannot run.
The run therefore came from a source state not represented by any commit. **Classified
`HISTORICAL STEPPER UNKNOWN`** — RKF45 is strongly suggested but not proven, and the success of
the run is explicitly *not* used to infer it. Reconstructing dirty-worktree history was out of
scope.

This is precisely the metadata hole §24 asks to close: the new baseline records its own stepper,
tolerances, cadence, EOS hashes, source commit and build configuration, so it can never become
unreproducible the same way. The historical outputs remain **superseded as validation references**.

## 8. EOS provenance — dual representation

CompactStar consumes two forms; both are recorded.

**Processed** (consumed by TOV via `NStar::SolveTOV_Profile`):

| SHA256 | File |
|---|---|
| `5747dd73256c0c28bc56be337cbb96d0918a54bc9ed9fc40984c5befd47ae5dd` | `DS-CMF-1-with-crust/DS(CMF)-1_with_crust.eos` |

**Raw CompOSE cold** (the tables the processed form derives from):

| SHA256 | File |
|---|---|
| `416444999ccac569e2c9b34808888949c36d759f30cce25dab0d42c13e900ce3` | `eos.thermo` |
| `d9c8e78c2fcf37fe770fecfc2d3a211d840a28299821a56c77e66f9ff74edef8` | `eos.nb` |
| `1a37b9563c40962b203e7bca1aa3b41e8c8b1427953df68095a51dd2cc17ff96` | `eos.t` |
| `1a37b9563c40962b203e7bca1aa3b41e8c8b1427953df68095a51dd2cc17ff96` | `eos.yq` |
| `4e69b9193e0f075584239d818e1b459791da4d12427531914c86cdd209c898a8` | `README` |
| `8b4472405295655cf530572af7edd7448efd0393d0bf9ad86ead3e4c87228c90` | `eos.init` |

**Raw CompOSE finite-`T`** (consumed directly by `CompOSE_Thermo`):

| SHA256 | File |
|---|---|
| `a456fb8595208ddf3119350a856fbf2b906c0a0e19bb7c716571748d0aa0724b` | `eos.thermo` |
| `3f79dbcc6f8b519696377f89ebc86464bc55cd61d9e2459f6e21e2d9e00f380d` | `eos.nb` |
| `2e4c6ec1feb85b16d0ee7036dce183782a9f681577e79c72315171069aa8513d` | `eos.t` |
| `d98fcd2f7752039c552c2ef2d04ab485b75db47a61f8ae1740875b54bf9824fd` | `eos.yq` |
| `48af68b1b4f6727252ae0051fc35c6445240241d654566973a068d20e1f35222` | `README` |
| `412f6739c769df650b3238a6e4b6f0d0f2d7d4a5df43e7d16c80e913bcaddbfb` | `eos.init` |

Supplied via `COMPACTSTAR_EOS_DATA_ROOT`. **No EOS data is committed** and no absolute path appears
in committed source. The EOS converter itself was not validated — a dedicated task may follow.

## 9. Running it

```bash
cmake -S . -B build -DPython3_EXECUTABLE="$(command -v python3)" \
      -DCOMPACTSTAR_EOS_DATA_ROOT=/path/to/data/compose
cmake --build build -j8
ctest --test-dir build -L regression --output-on-failure
```

Registered only when the data root is supplied; labels `thermal;scientific;external-data;regression`;
~11 s. Without the data root the regression is **not registered at all** — absence of data can
never masquerade as execution of the scientific regression. Data-free suite: 4/4 pass. With data:
6/6 pass.

## 10. Limitations

1. **Not a physics validation.** `Q0_DU = 1e27` and `Q0_MU = 1e21` are placeholders; a justified
   change moves every value here and requires re-baselining.
2. **One cadence fails.** `samples_per_decade = 300` produces a non-finite state at t ≈ 302 yr;
   cause undiagnosed. Seven other cadences from 50 to 400 complete and agree to 4.3e-6.
3. **The early transient is stiff** (§3.6). RKF45 is adequate on the evidence, but a stiff method
   with a real Jacobian is the likely long-term answer.
4. **Single star, single EOS, single platform.** 1.6 M☉ CMF, `Debug`, macOS/arm64, GSL 2.7.1.
5. **No CI.** The regression must currently be run by hand.
6. **Crust thermodynamics remain absent** from `C_⋆` (see `HEAT_CAPACITY_V1.md` §2).
