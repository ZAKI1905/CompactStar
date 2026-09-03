# Phase 4C-I1 — atomic replacement of the O(Ω²) candidate by the governed monopole response

> **FORMAL STATUS: `PHASE-4C-I1 GOVERNED HARTLE MONOPOLE REPLACEMENT COMPLETE — INDEPENDENT VALIDATION REQUIRED`**
>
> This increment **executes** the `GOVERNANCE.md` §3.1 correction authorized by ADR-0007 §9: the
> AI-authored O(Ω²) candidate of commit `675b4a9` is **deleted**, and the governed fixed-`ε_c`,
> `Ω²`-normalized monopole response replaces it **in the same commit**. There is no intermediate
> state in which one exists without the other.
>
> **Nothing here is validated physics.** Conformance to an accepted contract is not verification.
> Every check in this record protects an *implementation contract*; not one of them computes the
> monopole solution by an independent route. That is **Phase 4D**. Until 4D completes, INV-08
> stands at *source conformed, independent physical validation pending*, **no monopole baseline
> exists, and none may be created**.

| Field | Value |
|---|---|
| **Starting HEAD** | `a97f9c529954ed019a130f008d898b2896a89a3d` — upstream equal, tree clean, **7 ahead / 0 behind** `master` = `df859b5a73c4cac0c115f240744d89ce9f830b8d` |
| **Change class** | **scientific-semantic** (equations, perturbation variable, boundary condition, normalization, `δM` definition) **and structural** (public O(Ω²) API replaced) — strictest applies (`GOVERNANCE.md` §2), executed under **§3.1** |
| **Governing authority** | **ADR-0007** (ACCEPTED 2026-09-02) P1–P14 and §9; ADR-0006 (first-order normalization, no implicit spin); ADR-0003 (provenance/`Version()`); ADR-0004 (metric-factor owner); ADR-0005 (`I` and the goldens bitwise); INV-02, INV-05, INV-06, INV-13 |
| **Derivation of record** | `docs/validation/PHASE4C_HARTLE2_DERIVATION.md` — primary-source derivation from Hartle (1967) ApJ 150, 1005, eqs. (97), (100), (105)–(108) |
| **Prerequisite** | `docs/validation/PHASE4C_I0_EOS_DERIVATIVE.md` — the `dε/dp` authority, complete at `a97f9c5` |
| **Baseline** | pre: **25/25** (209.72 s) + **13/13** (15.62 s), seven hashes as `PHASE3_CLOSEOUT.md` §1; post: **27/27** + **14/14**, same seven hashes (§18) |

---

## 1. Starting state

```
worktree  /Users/keeper/Documents/CompactStar/worktrees/CompactStar-rotation
branch    physics/rotation-correctness      HEAD = a97f9c5   upstream equal   clean
master    df859b5                            7 ahead / 0 behind
```

Pre-task suites on that tree: **25/25 PASS** (209.72 s) and **13/13 PASS** (15.62 s), confirmed
from CTest's own `LastTest.log`. Seven durable artifacts unchanged, in particular
`hartle_I_dscmf1_debug.tsv = ddf018579364e9b3bf24f1e3a3e2577e70fb09b8b84eb718237d9da683aa9d15`.
No baseline was regenerated before the mutation.

## 2. `GOVERNANCE.md` §3.1 — the seven conditions, restated and executed

ADR-0007 §9 authorizes this pre-baseline correction. Each condition, and what this increment
actually did about it:

| | Condition | Execution |
|---|---|---|
| **1** | **Invalid / superseded behaviour identified.** `ODE_Hartle2_N_Fast` and `SolveHartle2_N` at `a97f9c5`: dimensionally inconsistent homogeneous terms, an inverted `j²` with `e^{−2ν}` missing, a source term wrong in sign, factor and power of `r`, an equation correct for neither `δp₀` nor `p₀*`, a non-Hartle `δp₀(R) = 0` shooting condition that shifts the central density, a `δM` missing `J²/R³`, profile-FD `dε/dp` with a `1.0` fallback, and every output quadratic in the arbitrary internal seed | **Deleted** (§4, §5) |
| **2** | **Why no baseline of it may be captured.** A golden of `δM/M = −1.6`, `ξ₀(R) = 0` and a shifted `ε_c` would make every later correction register as a regression | **No candidate output was executed, exported or baselined at any point in this increment** |
| **3** | **Minimum correction only.** ADR-0007 P1–P10, ordinary-`NStar` `l = 0`. No `l = 2`, no MixedStar second order, no rotochemical coefficients, no surface-convention change | **Held** — see §20 for the exact non-scope |
| **4** | **Independent verification required.** Deferred to Phase 4D and **mandatory before any scientific monopole baseline** | **Not performed here, and explicitly not claimed** (§13, §21) |
| **5** | **Narrow scope.** The second-order ordinary-`NStar` path and its result/API; first order protected | **Held** — first-order arithmetic byte-identical (§17) |
| **6** | **Historical outputs are not references.** No number ever produced by `SolveHartle2_N` / `GetHartleResult` — including the diagnostics in `PHASE4_ROTATION_ENTRY.md` §12 — is a reference result | **Recorded; nothing was compared against them** |
| **7** | **Baseline policy.** No monopole golden in I1; the first one is permitted only after 4D independent verification | **No `hartle_monopole_dscmf1_debug.tsv` was created** (§18) |

> **§3.1 CORRECTION EXECUTED — INDEPENDENT VERIFICATION PENDING — NO BASELINE YET**

## 3. Candidate authentication before deletion

Extracted from `CompactStar/Core/src/RotationSolver.cpp` at `a97f9c5` and hashed before a single
edit:

| Retired definition | Lines | SHA-256 of the extracted text |
|---|---|---|
| `ODE_Hartle2_N_Fast` | 1140–1226 | `b13328ca0fd25946d84d28284b3b66299036b5e5a771caa0525f776f3a834613` |
| `SolveHartle2_N` | 1242–1386 | `d803c65bf4b68e15a106833c53edced2d677f0b76715ee35a3000d666d985308` |
| `ODE_Hartle2_Mixed_Fast` | 1230–1236 | `4c5779580c1a88a5579f76044ac2f4bf43df69e9b51fafaa9d178f702ff9b3b7` |
| `SolveHartle2_Mixed` | 1390–1393 | `415d1a0941114316b7891aa2a349bec1579b74bbfd2845de1366369228af9e73` |
| `GetHartleResult` | 1061–1064 | `5e89d4ffa78f1a22025ef64c55a4319df37b40608c658463ee47650d99ec53e4` |
| `struct HartleResult` (`RotationSolver.hpp`) | 186–223 | `de5ae2a11b331551df8c9c6380b819020aaa816d4aed66028e9a25df66851a9e` |

These hashes exist so the retirement is auditable, **not** so the behaviour can be restored.

## 4. Reference audit of every associated member — proven before deletion

Nothing was deleted on the assumption that it was candidate-only. Every member was searched
across `CompactStar/`, `tests/` and `main/`:

| Member | Live uses found | Verdict |
|---|---|---|
| `include_m0p0_source_` | candidate only (`:1187`, `:1280`, `:1314`, `:1347`) | **deleted** |
| `fast_dEdP` | candidate only (`:1147`, `:1299`, `:1332`) | **deleted** |
| `fast_dEdP_v`, `fast_dEdP_d` | **declared, never assigned or read** | **deleted** (dead since `675b4a9`) |
| `fast_nu` | **declared, never assigned or read** | **deleted** (dead) |
| `fast_nu_prime` | candidate only (`:1148`, `:1298`, `:1331`) | **deleted** |
| `fast_omega_bar`, `fast_domega_bar` | candidate only (`:1189-1190`, `:1300-1301`, `:1333-1334`) | **deleted** |
| `fast_p`, `fast_e`, `fast_m` | candidate only. Their own declaration comment said so — *"used ONLY by the second-order O(Omega^2) candidate path … Not read by any first-order code, which uses the profile-backed `fast_p_`/`fast_e_`/`fast_m_`"* — and the audit confirmed it: the three other appearances are commented-out legacy lines | **deleted**, with those five dead comment lines removed alongside so nothing names a member that no longer exists |
| `hartle_result_` | candidate, **plus one line in `FindNMomInertia`** (`:807`) that existed solely to hand the candidate a grid pointer | **deleted**, including that line (§17) |
| **`stored_omega_bar_`, `stored_domega_bar_`** | **first-order machinery** — written at `:730-731`, `:753-754` and divided by `Ω_raw` at `:795-797` to build the seed-free response | **KEPT.** These are no longer candidate-only; deleting them would have broken ADR-0006's first-order contract |
| **`fast_r_`, `fast_p_`, `fast_e_`, `fast_m_`, `fast_k_`, `fast_r_mix_`, `fast_p_tot_`, `fast_e_tot_`, `fast_m_tot_`, `fast_k_mix_`** | first-order profile-backed interpolation (`EvalFastPEM_`, `SetFastProfilePtrs_`) | **KEPT** |

## 5. Atomic retirement

Deleted in the same commit that adds the replacement (ADR-0007 Q13 = A): `SolveHartle2_N`,
`ODE_Hartle2_N_Fast`, `SolveHartle2_Mixed`, `ODE_Hartle2_Mixed_Fast`, `GetHartleResult`,
`struct HartleResult`, and every member of §4 marked deleted.

**No compatibility wrapper. No name reused with different semantics. No temporary fail-closed
stub** — the replacement lands with the deletion, so a stub would have described a state that
never exists. Zero repository callers meant no caller migration was required.

The MixedStar second-order stubs (`SolveHartle2_Mixed`, `ODE_Hartle2_Mixed_Fast`) belonged to the
retired candidate API and are deleted with it; **no new MixedStar physics replaces them**.

**Phase-5 replacement debt.** The uncompiled `RotochemicalCache` candidate still names the retired
type textually (`Physics/Evolution/RotochemicalCache.hpp:44,135,143,163`;
`src/RotochemicalCache.cpp:51,109`). It is in **no** CMake source list — verified by searching
every `CMakeLists.txt` — so it is not compiled, has zero callers, and would not compile as written
in any case (INV-09). It was **not** repaired or activated here. Recorded as **Phase-5 replacement
debt governed by `GOVERNANCE.md` §5 and INV-09**: whatever Phase 5 builds must consume the
governed response, not a resurrected `HartleResult`.

## 6. The governed equations, transcribed from ADR-0007 P2

Implemented in `RotationSolver::ODE_HartleMonopole_`, in the normalized variables
`y[0] = m̂₀ = m₀/Ω_geom² [km³]`, `y[1] = p̂₀* = p₀*/Ω_geom² [km²]`, driven by the **verified**
seed-free first-order response `s = ω̄/Ω` and `s' = ω̄'/Ω`:

```
dm̂₀/dr  = 4π r² (ε+p)(dε/dp) p̂₀*                          [term 1]
        + (1/12) r⁴ e^{−2ν}(1 − 2m/r) s'²                  [term 2]
        + (8π/3) r⁴ (ε+p) e^{−2ν} s²                       [term 3]

dp̂₀*/dr = − m̂₀ (1 + 8π r² p)/(r − 2m)²                     [term 4]
         − 4π (ε+p) r² p̂₀*/(r − 2m)                         [term 5]
         + (1/12) r³ e^{−2ν} s'²                            [term 6]
         + (2/3) r e^{−2ν} s (s + r s' − r ν' s)            [term 7]
```

Each term carries its `[term N]` marker in the source, and the block comment above the function
cites ADR-0007 P2 and the derivation record. The code is **deliberately not algebraically
compressed**: every term must stay visibly traceable to the contract it came from. Term 7 is the
analytically expanded form of Hartle's `(1/3) d/dr[r³j²ω̄²/(r−2m)]`, expanded in the derivation
record §7 precisely so that **no numerical differentiation of the source appears in the RHS**.

**Nothing was copied from the candidate.** The candidate's formulas are not an authority even
where a term happens to resemble one.

## 7. RHS inputs at the actual driver radius

The candidate's structural error was to assign per-node scalars before each
`gsl_odeiv2_driver_apply` and let its RHS use them at whatever internal radius GSL chose — a node
value standing in for an inter-node one. That is not repeated.

`RotationSolver::MonopoleBackground_` holds pointers to the nine profile/response columns and
samples **all eight inputs — `p`, `ε`, `m`, `ν`, `ν'`, `dε/dp`, `s`, `s'` — at the true RHS
radius**, through **one shared bracket index** so they always refer to the same radial interval.
The interpolation is **linear**, the order the profile itself carries (INV-13), and the bracket
walk mirrors `EvalFastPEM_`, which the validated first-order solver already uses for exactly this
purpose. Outside the tabulated range it clamps to the end node, the same convention
`EvalFastPEM_` uses.

**Metric factor.** `1 − 2m/r` comes from the canonical owner
`CompactStar::Geometry::MetricDenominator(r, m)` (ADR-0004) — no second definition, **no clamp,
no `1e-15` regularization**; `r − 2m` is formed as `r · D`. That primitive *throws* on
out-of-domain input, and a C++ exception must not cross the GSL C callback boundary, so the RHS
catches and returns `GSL_EBADFUNC`; a non-finite derivative does the same.

**`ν'` units.** `NStar::BuildFromTOV` converts `TOVPoint::nu_der` from cm⁻¹ to km⁻¹ when it fills
`MetricNuPrime` (`NStar.cpp:187`), so the profile column — and therefore the ODE — consumes
`ν' [km⁻¹]`. The original cgs derivative is never touched, and no second `ν'` formula exists.

## 8. Regular-centre start (ADR-0007 P4)

At the first strictly-positive grid radius `r₀` (the first-order convention, INV-05), with
`j² = e^{−2ν(r₀)}(1 − 2m(r₀)/r₀)`:

```
p̂₀*(r₀) = (1/3) j² s(r₀)² r₀²
m̂₀(r₀)  = (4π/15) (ε+p)(r₀) [ (dε/dp)(r₀) + 2 ] j² s(r₀)² r₀⁵
```

**Not** the candidate's literal `{0,0}`. There is no homogeneous admixture, no shooting and no
surface boundary condition of any kind: this is the unique regular solution of the fixed-`ε_c`
family (H67 eq. 108, p. 1009, p. 1022; HT68 §II f). Verified to **relative 0.0 — exact** against
an independently recomputed series in checks `C5a`/`C5b`.

## 9. `dε/dp` consumption

Read **only** through `StarProfile::HasEosDEdP()` / `GetEosDEdP()`, the Phase-4C-I0 authority.
The rotation solver never touches an EOS object, never differentiates the radial profile, and has
no fallback value. A star that carries no authoritative derivative **fails the whole computation**
and publishes nothing (`C1`; detector `D3`). An explicit `0.0` — the correct value for
incompressible matter — remains a valid input and is exercised throughout the contract test.

## 10. Result types and the public API

```
HartleMonopoleResponse            // seed-free, fixed-eps_c, per Omega_geom^2
  m0_over_Omega2(r)               [km^3]
  p0star_over_Omega2(r)           [km^2]
  delta_p0_over_Omega2(r)         [1]      = (eps+p) * p0star_over_Omega2
  xi0_over_Omega2(r)              [km^3]   = p0star_over_Omega2 / nu'
  deltaM_over_Omega2              [km^3]
  surface_shell_mass_over_Omega2  [km^3]
  surface_xi0_over_Omega2         [km^3]
  R_surface [km] , I [km^3] , r_grid , source_profile , source_version , valid
  At(AngularVelocity) -> PhysicalHartleMonopole
```

Every field name encodes the `Ω²` normalization; there is no bare `m0`, `p0` or `deltaM` among
the canonical coefficients. `PhysicalHartleMonopole` carries the materialized `m0`, `p0star`,
`delta_p0`, `xi0`, `delta_M`, `surface_shell_mass`, `surface_xi0` with **one** canonical stored
angular velocity (`Omega_geom`; the physical value is derived by `OmegaRadPerSecond()`, never
stored alongside).

Public entry points on `NStar`:

| Method | Semantics |
|---|---|
| `bool ComputeHartleMonopoleResponse()` | **performs work** — the name says so. One ODE integration; recomputes when absent or stale, reuses otherwise; fails closed |
| `const HartleMonopoleResponse *MonopoleResponse() const` | cheap, read-only, **never integrates**; returns `nullptr` when absent **or stale** |
| `PhysicalHartleMonopole MonopoleAt(AngularVelocity) const` | pure scaling by `Ω_geom²`; throws when no current response exists |

**There is no innocent-looking const getter that secretly runs an expensive ODE.**

## 11. Provenance, staleness and atomic publication

The response records the **identity** of its source `StarProfile` and that profile's
**`Version()`** (ADR-0003). `MatchesSource()` gates every read: a response whose provenance no
longer matches is **stale and is never returned as current**, and `MonopoleAt` refuses to
materialize from it. Recomputation is explicit, through the same public API. No Core → Evolution
dependency was introduced; the mechanism is local to Core.

**Publication is atomic.** Everything is built in local storage; `monopole_response_` is assigned
exactly once, at the end, and only after: the profile is complete, the first-order response is
valid, the EOS derivative is present and correctly sized, the integration completed, every node
value is finite, and every surface quantity is finite. On any failure the cached response is
cleared, so a partial or stale result can never be mistaken for a current one.

**An observation, not a defect introduced here.** `NStar::GetSequence()` has a non-`const`
overload returning `prof_.SeqMutable()`, which calls `Touch()` — so *reading* the sequence
through a mutable star bumps the profile version and correctly invalidates every version-keyed
cache. That is pre-existing ADR-0003 behaviour; the contract test reads through a `const`
reference and documents it inline. **No first-order architecture was broadened** to accommodate
it, and none needed to be.

## 12. Derived fields, surface semantics and `δM`

- `delta_p0_over_Omega2 = (ε+p)·p̂₀*` — dimensionless.
- `xi0_over_Omega2 = p̂₀*/ν'` [km³]. Both vanish at the centre (`p₀* ∼ r²`, `ν' ∼ r`), so the
  ratio is finite there but the division is not always well conditioned on the first node. Where
  it is not, the **regular-centre series** `ξ̂₀ → j_c² s_c² r / [4π(ε_c + 3p_c)]` (derivation
  record §9.2) is used — **not a zero, and not an arbitrary epsilon**, either of which would
  replace a real value with a fabricated one. If a derived field still cannot be produced
  finitely, the response fails as a whole.
- **`R_*` is the last profile node** — the production EOS-floor surface (INV-06), documented as
  such in the type, in the source and in the test output. It is **never** called the exact `p = 0`
  radius, and no extrapolation or surface reconstruction was added.
- `δM` exactly as ADR-0007 P6:

  ```
  shell_hat  = 4π R_*² ε_* ξ̂₀(R_*)
  deltaM_hat = m̂₀(R_*) + shell_hat + I²/R_*³
  ```

  with `I` the **verified first-order moment of inertia** — never a raw `J`, never a seed-dependent
  quantity, never an arbitrary `Ω`. The shell term is **computed, never assumed small**: §13 shows
  it is `~10⁻⁶` of `δM̂` on the CMF stars and **90 %** of it on the constant-density star.

## 13. Validation performed

### 13.1 Self-contained contract — `tests/rotation/hartle_monopole_contract.cpp` (22/22)

Exact Schwarzschild constant-density interior, `dε/dp = 0` supplied explicitly through the
Phase-4C-I0 mechanism. Bounds predeclared from the accepted contract, not chosen here:
seed invariance `1e-10` (ADR-0007 §7-1), quadratic materialization `1e-14` (§7-2), `δM` identity
`1e-14` (arithmetic).

| Check | Result |
|---|---|
| `C1` no authoritative `dε/dp` ⇒ compute fails, no response, `MonopoleAt` throws | pass |
| `C2` explicit `dε/dp = 0` accepted; response computes | pass |
| `C3` one value per radial node in every column (n = 601) | pass |
| `C4` every published value finite | pass |
| `C5a` `p̂₀*(r₀)` equals the governed regular series | **rel 0.0 (exact)** |
| `C5b` `m̂₀(r₀)` equals the governed regular series | **rel 0.0 (exact)** |
| `C6` provenance records profile identity and `Version()` | pass |
| `C7a` a version bump makes the response stale; it is not returned | pass |
| `C7b` materializing from a stale response is refused | pass |
| `C7c` explicit recomputation restores it, with one new solve | pass |
| `C8` zero spin materializes **exact** zeros | pass |
| `C9` `+Ω` and `−Ω` materialize **bit-identical** perturbations | pass |
| `C10` `Q(2Ω) = 4Q(Ω)` | **worst rel 0.0** (bound `1e-14`) |
| `C11a` ordinary construction runs **zero** monopole solves | pass |
| `C11b` ordinary construction still runs the first-order solve | pass |
| `C11c` one explicit request ⇒ exactly **one** integration | pass |
| `C11d` recompute on an unchanged profile ⇒ still one | pass |
| `C11e` 25 further materializations ⇒ **zero** additional integrations | pass |
| `C12` first-order `I` bit-identical | `1.571329385870e+02` |
| `C13a` `δM̂ = m̂₀(R_*) + shell + I²/R_*³` | **rel 0.0** |
| `C13b` the shell term is present and non-zero | `shell = 1.3142e3`, `m̂₀(R_*) = 1.4196e2`, `I²/R_*³ = 1.1238e1` km³ |
| `C14a/b` no response before computing; computing confers no spin | pass |

**Seed invariance** (§25 of the task; ADR-0006 P9): the internal first-order seed was varied over
`{1e-3, 1e-1, 1e1, 1e3}` through the privileged test seam, and every monopole coefficient
compared against the default-seed response:

| seed | worst relative difference |
|---|---|
| `1e-3` | `5.2353e-15` |
| `1e-1` | `3.3863e-15` |
| `1e+1` | `7.1050e-15` |
| `1e+3` | `7.8529e-15` |

**Worst `7.85e-15` against a `1e-10` bound.** The solver consumes only `ω̄/Ω` and `ω̄'/Ω`; the raw
`stored_omega_bar_` arrays are unreachable from it by construction, and detector `D1` shows what
happens if that is violated.

### 13.2 Real-EOS conformance — `tests/rotation/hartle_monopole_cmf.cpp` (37/37)

DS(CMF)-1_with_crust at 1.0/1.4/1.6/2.0 M☉. Each star: derivative present, first-order response
valid, response computes, provenance matches, sizes correct, every value finite, `δM̂` identity
exact, `R_*` equals the last profile node, materialization exactly quadratic.

**DIAGNOSTIC — not a golden, not compared against the retired candidate, not a scientific claim.**

| M [M☉] | R_* [km] | nodes | I [km³] | m̂₀(R_*) | p̂₀*(R_*) | ξ̂₀(R_*) | shell_hat | I²/R_*³ | δM̂ |
|---|---|---|---|---|---|---|---|---|---|
| 1.00 | 13.42632 | 2629 | 8.6996e1 | 9.3099e2 | 3.5338e1 | 3.3647e3 | 9.7842e-4 | 3.1270e0 | 9.3411e2 |
| 1.40 | 13.54532 | 2646 | 1.3562e2 | 8.8553e2 | 3.4252e1 | 2.1118e3 | 6.3334e-4 | 7.4004e0 | 8.9293e2 |
| 1.60 | 13.46832 | 2635 | 1.5959e2 | 8.1555e2 | 3.2446e1 | 1.6171e3 | 6.2260e-4 | 1.0424e1 | 8.2597e2 |
| 2.00 | 12.71232 | 2527 | 1.9372e2 | 5.3788e2 | 2.4053e1 | 7.0452e2 | 2.5921e-4 | 1.8268e1 | 5.5615e2 |

Units: `m̂₀`, `ξ̂₀`, `shell_hat`, `δM̂` in km³; `p̂₀*` in km².

Two **structural** observations, offered to Phase 4D as starting points and **explicitly not as
validation**: the shell term is `~10⁻⁶` of `δM̂` here and `90 %` of it on the constant-density
star, which is the behaviour the derivation record §11.2/§12.2 predicted and the reason ADR-0007
P6 computes it rather than assuming it small; and `δM̂` is **positive** on every star, where the
retired candidate produced `δM/M = −1.6`. Whether the magnitude is right is a 4D question that
this increment does not answer.

## 14. Detectors

Contract/architecture detectors only. ADR-0007's scientific `M1`–`M9` suite belongs to 4D, where
independent oracles exist. Each mutation: one unique site, rebuilt, measured, restored and
verified **byte-identically by SHA-256** (`RotationSolver.cpp` → `db7d213d…`, MATCH), with
`os.utime` after every write.

| | Mutation | Fired — how |
|---|---|---|
| **D1** | feed the raw `stored_omega_bar_`/`stored_domega_bar_` to the solver instead of `s`, `s'` | contract **FAILS**: seed invariance `Sa` at **`4.0e10`** against its `1e-10` bound, plus `C5a`/`C5b` (the centre series no longer matches). The CMF test still passes — a seed leak is invisible to a finiteness check, which is precisely why the seed-invariance contract exists |
| **D2** | scale `m0` linearly in `Ω` instead of `Ω²` on materialization | contract **FAILS** `C9` (± parity) and `C10` at **`5.0e-1`**; CMF **FAILS** on all four stars |
| **D3** | bypass the fail-closed `dε/dp` guard with the candidate's `1.0` fallback | contract **FAILS** `C1` — a derivative-absent star now silently produces a "response" |
| **D4** | drop the provenance check in `MonopoleResponse()` | contract **FAILS** `C7a` and `C7b` — a stale response is served as current |

No detector mutates a governed ODE term: inflating the count by breaking physics that has no
independent oracle yet would prove nothing at this stage.

## 15. Candidate deletion proof

Search of the whole tree for the retired identifiers:

| Symbol | Live definitions | Public declarations | Compiled callers |
|---|---|---|---|
| `SolveHartle2_N` | 0 | 0 | 0 |
| `ODE_Hartle2_N_Fast` | 0 | 0 | 0 |
| `SolveHartle2_Mixed` | 0 | 0 | 0 |
| `ODE_Hartle2_Mixed_Fast` | 0 | 0 | 0 |
| `GetHartleResult` | 0 | 0 | 0 |
| `HartleResult` | 0 | 0 | 0 |
| `include_m0p0_source_` | 0 | 0 | 0 |

Remaining textual occurrences, all classified:

- **Historical documentation** — `AngularVelocity.hpp:26`, `RotationSolver.hpp` (the new type's
  doc comment naming what it replaced), `RotationSolver.cpp` (the section header recording the
  retirement), and the `docs/` records. Deliberate.
- **Uncompiled Phase-5 candidate** — `RotochemicalCache.{hpp,cpp}`, §5. Not compiled, not
  repaired, recorded as Phase-5 debt.

## 16. Automatic-solve and performance contract

| Situation | Monopole integrations |
|---|---|
| ordinary `NStar` construction (TOV → profile → `Find_MomInertia` → `I` + first-order response) | **0** (`C11a`) |
| first explicit `ComputeHartleMonopoleResponse()` | **1** (`C11c`) |
| repeat request, profile unchanged | **0** — cache reused (`C11d`) |
| 25 materializations at different spins | **0** (`C11e`) |
| after a profile version bump, explicit recompute | **1** (`C7c`) |

Existing workflows pay nothing for a second-order solution they never asked for.

## 17. First-order protection

Each first-order function was extracted from `a97f9c5` and from the working tree, stripped of
comment-only and blank lines, and compared:

| Function | Arithmetic identical? |
|---|---|
| `ODE_N_Fast` | **YES** |
| `GetHartleOmegaCoeff_N_Fast` | **YES** |
| `GetHartleDOmegaCoeff_N_Fast` | **YES** |
| `HartleFirstOrderResponse::At` | **YES** |
| `FindNMomInertia` | one code line removed, and only one: `hartle_result_.r_grid = nstar_ptr->Profile().GetRadius();` — a write to the retired struct. **No first-order arithmetic changed.** |

Corroborated by `C12` (`I` bit-identical on the analytic star at `1.571329385870e+02`), by the
unchanged `hartle_I_dscmf1_debug.tsv` hash, and by both Hartle first-order CTests and the two
Phase-4B physics tests passing unchanged. No physical-normalization behaviour changed.

## 18. Existing scientific artifacts

All seven durable artifacts **byte-identical**, `hartle_I_dscmf1_debug.tsv` in particular at
`ddf018579364e9b3bf24f1e3a3e2577e70fb09b8b84eb718237d9da683aa9d15`. No rebaseline. **No
`hartle_monopole_dscmf1_debug.tsv` was created**, and none may be until 4D.

## 19. Test counts

| Configuration | Before | After |
|---|---|---|
| authenticated (external data root) | 25/25 | **27/27** |
| self-contained | 13/13 | **14/14** |

Two new executables: one self-contained contract, one external-data conformance. Not one per
star.

## 20. Exact non-claims and non-scope

**Non-claims.** Nothing in this increment verifies the O(Ω²) physics. The regular-centre start is
verified against a recomputation of the same series, not against an independent solution; the
`δM̂` identity is arithmetic on published fields; the CMF table is a diagnostic. **No number
produced here is validated physics, and no result of it may be cited.**

**Non-scope.** No `l = 2` equation of any kind. No rotochemical coefficients — no `A_i`, `B_i`,
`Z_i`, `dn_i/dp`, particle-number perturbation integrals, chemical imbalance, weak rates or
heating. No public homogeneous/sequence-derivative response (ADR-0007 P11 as modified at
acceptance; 4D may build one test-side). No MixedStar second-order physics. No change to the TOV
integration, `Geometry` mathematics, thermal code, `RotochemicalCache`, or the BNV/Decay/DarkCore
modules. No surface-convention change, no INV-06 change. No monopole baseline.

**Numerical settings are provisional.** `rk8pd`, `h₀ = 1e-1`, absolute and relative tolerance
`1e-10` mirror the first-order driver because that is a reasonable engineering starting point —
**not** because the retired candidate used them, and they are inherited as no kind of scientific
authority. Their adequacy and the convergence of the result under radial refinement are Phase-4D
questions (ADR-0007 §7 item 9). Nothing was tuned to make an output look desirable.

## 21. 4D readiness

In place for Phase 4D: a conforming, seed-free, fixed-`ε_c` implementation with a stable public
API; provenance that makes staleness detectable; an authoritative `dε/dp`; a contract test that
already exercises the fail-closed paths; and a CMF harness that builds four stars and computes the
response in one run.

4D must supply what this increment deliberately does not — independent scientific validation:
regular-centre series against an **independent** solver in different variables (the Eulerian
`(m₀, δp₀)` pair, and/or `(m₀, h₀)`); Newtonian homogeneous-star limits (`δM̂ → R³`,
`3Mξ̂₀/R⁴ → 1`); the exterior identity `m̂₀ + I²/r³` constant across near-vacuum nodes; the
conservation identity `p̂₀* + ĥ₀ − (1/3)r²e^{−2ν}s²` constant; published Hartle & Thorne (1968)
Table 5 models; radial convergence; EOS-derivative sensitivity; the sequence-derivative identity
via a test-side homogeneous solve; and detectors `M1`–`M9`. **Only after that may the first
monopole scientific baseline be created** (§3.1 condition 7).

## 22. Status

> **`PHASE-4C-I1 GOVERNED HARTLE MONOPOLE REPLACEMENT COMPLETE — INDEPENDENT VALIDATION REQUIRED`**

**Exact next task.** Phase 4D — independently validate the new O(Ω²) monopole response against
the primary-source identities, the regular-centre series, an independent solver in different
variables, Newtonian homogeneous-star limits, published Hartle–Thorne results, radial convergence,
EOS-derivative sensitivity, and the planned `M1`–`M9` detectors; only after that may the first
monopole scientific baseline be created.
