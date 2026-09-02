# Phase 3D — Canonical proper-volume ownership: implementation and validation

> **STATUS: `PHASE-3D CANONICAL PROPER-VOLUME OWNERSHIP COMPLETE`**
> with the explicit qualifier
> **`LEGACY PATH-1 / MIXEDSTAR / CANDIDATE MIGRATIONS REMAIN DEFERRED`.**
>
> **INV-04 is NOT fully resolved.** The canonical validated path conforms; four categories of
> legacy site are governed by ADR-0004 but have not migrated, each for a stated reason (§8).

| Field | Value |
|---|---|
| **Starting HEAD** | `5eb6bdd314b34b1acccb2df59a92eb7eb813c145` (8 ahead / 0 behind `master` `d2040d89…`, upstream equal, clean) |
| **Acceptance commit** | `f92df9a` — `docs: accept proper-volume ownership contract` |
| **Governing authority** | **ADR-0004 (ACCEPTED 2026-09-01)**, rank 2; INV-03, INV-04, INV-13, INV-14 |
| **Change class** | **STRUCTURAL / ARCHITECTURE** *and* **SCIENTIFIC-SEMANTIC** (degenerate-input behavior becomes governed) |
| **GOVERNANCE §3.1** | **NOT invoked** |

---

## 1. Owner adjudication, as implemented

| Question | Decision | Where it landed |
|---|---|---|
| **Q1** ownership boundary | **Option B** — three owners: a dependency-neutral mathematical primitive; `GeometryCache` as cached representation; consumers keep their own physics factors and unit conversions | `CompactStar/Geometry.hpp`; `GeometryCache` unchanged in composition; `NStar` keeps `n_B` and `1e54` |
| **Q2** `MixedStar` | **governed now, source migration deferred** | no `MixedStar` source touched (§8) |
| **Q3** degenerate inputs | **hybrid physical-domain contract** — exact regular-centre limit, fail closed otherwise, **no `1e-15` clamp** | `Geometry.hpp` domain rules; `proper_volume_contract` group B |

`Core` was **not** made to depend on `Physics/Evolution`: `NStar.cpp` includes no `GeometryCache`.

## 2. The primitive

`CompactStar/Geometry.hpp`, namespace `CompactStar::Geometry`. Four free scalar functions:

```cpp
double MetricDenominator(double r_km, double m_km);   // f = 1 - 2m/r
double Lambda(double r_km, double m_km);              // -0.5 * log(f)
double ExpLambda(double r_km, double m_km);           // exp(Lambda)
double ProperVolumeWeight(double r_km, double m_km);  // (4*pi) * r^2 * expLambda
```

Includes `<cmath>`, `<stdexcept>`, `<string>` and nothing else — no `Core`, `EOS`, `Physics`,
`Zaki` or Confind. No state, no registry, no class hierarchy, no templates. It therefore adds no
edge to the dependency graph, which is why it sits at the top level beside `Units.hpp`.

`e^Λ` is `exp(-0.5*log(f))`, **not** `f^{-1/2}`: the two agree to ≤ 2 ULP but are not bitwise
equal, and the logarithmic form is what every authenticated profile and cache already stores.
Expression order (`1.0 - 2.0 * m_km / r_km`) is load-bearing — it is what makes the migration
bit-identical on the ordinary domain.

**Domain contract, as accepted:**

| Case | Behavior |
|---|---|
| `r = 0`, `m = 0` | exact regular-centre limit: `f = 1`, `Λ = 0`, `e^Λ = 1`, `w_V = 0` |
| `r < 0` | `std::runtime_error` |
| `r = 0`, `m ≠ 0` | `std::runtime_error` |
| `r > 0`, `f ≤ 0` | `std::runtime_error` (horizon and beyond) |
| non-finite `r` or `m` | `std::runtime_error` |
| `r > 0`, `f > 0` | evaluate; **no clamp, no epsilon** |

No tolerance band exists around `r = 0`. `std::runtime_error` is the project's convention (123
existing uses; no other exception type appears in `CompactStar/`).

**A signed-zero detail, kept rather than papered over.** At the regular centre
`Λ = -0.5 * log(1.0) = -0.0`. That is exactly what the legacy `r <= 0 → denom = 1` branch
produced, so it is bit-identical to pre-3D behavior, and it composes correctly
(`exp(-0.0) = 1.0`, `1/1.0 = 1.0`). The contract test asserts the composition rather than
asserting the sign away.

## 3. Migrated sites — the validated path only

| Site | What changed | Standard | Result |
|---|---|---|---|
| `NStar::BuildFromTOV` — Λ production | clamped inline block → `Geometry::Lambda(r_km, m_km)` | bit-identical | **bit-identical** (all five artifacts byte-identical) |
| `NStar::BuildFromTOV` — baryon integrand | `/= (1 - 2m/r).sqrt()` → `w_V · n_B · 1e54` from the primitive | `\|ΔB\|/B ≤ 1.0e-15` predeclared | **1.368e-16** |
| `GeometryCache::DeriveLambdaFromMR_` | own formula + `eps` clamp → delegates; `eps` parameter removed | bit-identical | **bitwise, 96/96** |

The former baryon integrand divided through `Zaki::Vector::DataColumn::operator/=`, which
**silently substitutes 1 for a zero divisor** and logs — dropping the metric factor entirely at
that node. Nineteen call sites inherited that behavior without asking for it. The canonical path
no longer does.

## 4. Baryon-number migration — measured against a blind prediction

ADR-0004 §7.2 predicted these values from a scratch replication written **before** any
implementation existed. The pre-mutation capture on this machine reproduced all four **bitwise**,
which authenticates the ADR's evidence rather than quoting it.

| Star | `B` pre-3D | `B` post-3D | `\|ΔB\|/B` | `M`, `R`, `ec` |
|---|---|---|---|---|
| 1.0 M☉ | `1.2738887310935454e+57` | `1.2738887310935455e+57` | **1.368e-16** (1 ULP) | bitwise |
| 1.4 M☉ | `1.8321833625787515e+57` | `1.8321833625787515e+57` | **bitwise** | bitwise |
| 1.6 M☉ | `2.124575695479721e+57` | `2.124575695479721e+57` | **bitwise** | bitwise |
| 2.0 M☉ | `2.7457630611447906e+57` | `2.7457630611447906e+57` | **bitwise** | bitwise |

**Worst `|ΔB|/B` = 1.368e-16**, against the `1.0e-15` predeclared in ADR-0004 §7.1 from the
operation counts before implementation. The ADR's blind §7.2 prediction — 1 ULP on 1.0 M☉ only,
bitwise on the other three — was reproduced **exactly**. The bound was not widened.

## 5. Behavior preservation — all five protected artifacts

| Artifact | Required | Result |
|---|---|---|
| `tov_dscmf1_reference.tsv` | byte-identical | **byte-identical** |
| `hartle_I_dscmf1_debug.tsv` | byte-identical | **byte-identical** |
| `passive_cooling_cmf_1p6_debug.tsv` | byte-identical | **byte-identical** |
| `grid_convergence_cmf_1p6_trajectory.tsv` | byte-identical | **byte-identical** |
| `grid_convergence_cmf_1p6_debug.tsv` | byte-identical, or B-only within `1.0e-15` | **byte-identical** |

The conditional B-only exception §20 permitted was **not needed**: the grid-debug artifact
carries a `B` column and it did not move at any of the five resolutions. No golden was
re-baselined in Phase 3D.

## 6. New durable coverage

| Test | Kind | What it pins |
|---|---|---|
| `proper_volume_contract` | self-contained | regular-centre limit; 14 fail-closed cases; `f` bitwise vs the canonical subtraction; `e^Λ` vs `1/√f` (**1 ULP**, bound 2) and vs `exp(-½ln f)` (**bitwise**); `w_V` composition; flat-space limit; a magnitude guard proving the primitive carries no unit conversion |
| `geometry_cache_measure_contract` | self-contained | stored-Λ and derived-Λ routes produce **bitwise** `ExpLambda`/`WV`/`WVExpNu`/`WVExp2Nu`; cached `WV` equals the primitive node by node; `WVExpNu`/`WVExp2Nu` **composed from the one `WV` array**; ADR-0003 provenance undisturbed |
| `baryon_number_cmf` | external-data | `B` on the 1.0/1.4/1.6/2.0 M☉ sequence against a canonical reference, under the predeclared `1.0e-15`; structure under the same bound; a magnitude guard that the `1e54` conversion is applied exactly once |

`tests/baselines/baryon_number_dscmf1_reference.tsv` is a **canonical value table**, not an
old-vs-new diff — the latter would stop being useful once this migration is history. The pre-3D
values are recorded in its comment header, where they document the transition without shaping
the schema. Before Phase 3D, `B` had **no dedicated artifact at all**; it appeared only as one
column inside the grid-convergence debug file, where a change in it was indistinguishable from a
change in anything else.

## 7. Detector proof

Every mutation was reverted and the revert verified **byte-identical by SHA-256**.

| ID | Mutation | Fires in | Evidence |
|---|---|---|---|
| **D1** | metric denominator `1 − 2m/r` → `1 − m/r` | `proper_volume_contract`, `geometry_cache_measure_contract`, `baryon_number_cmf` | all three FAILED |
| **D2** | drop `e^Λ` from the **cached** `w_V` (coordinate volume) | `passive_cooling_regression`, `geometry_cache_measure_contract` | `C_⋆` moved **15–17 %** vs the golden (`rel = 0.150`–`0.175` against a `1e-5` tolerance) |
| **D2b** | drop `e^Λ` from the **primitive's** `w_V` | `proper_volume_contract`, `baryon_number_cmf` | both FAILED |
| **D3** | move the `1e54` baryon conversion **into** the primitive | `proper_volume_contract`, `baryon_number_cmf` | both FAILED |
| **D5** | violate the hybrid domain — reinstate the `1e-15` clamp | `proper_volume_contract` | FAILED (group B) |
| **D6** | compose `WVExpNu` from a **second** re-derived measure | `geometry_cache_measure_contract` | FAILED — `72/96` bitwise, worst rel `2.154e-16` |

**D4 was deliberately not executed.** It substitutes a component mass for the total mass in
`MixedStar`, which is not migrated and has zero coverage — there is nothing for it to fire in.
Manufacturing one would violate ADR-0004 §18.7's own instruction not to build a detector that
cannot legitimately fire.

**Two detectors initially did not fire, and both were diagnosed rather than reported as passes:**

- **D2 against `heat_capacity_v1`** — that fixture pushes `m = 0.0` at **every** node, so
  `e^Λ ≡ 1` identically and the test **structurally cannot** observe the metric factor. It is a
  flat-space synthetic fixture and is not a proper-volume detector.
- **D2 against `heat_capacity_real_star`** — that program contains **no assertions at all**: it
  has no `Report()` calls and fails only on load/solve errors. It is an evidence-emitting
  diagnostic harness, **not** a detector. This is worth recording as a coverage fact.
- **D6 as first written** — the mutation recomposed `WVExpNu` as `((4π·r²)·e^Λ)·e^ν`, which is
  the *same association* as `(w_V)·e^ν` and therefore bit-identical. It was a fake mutation with
  nothing to detect. Replaced with a genuine reassociation, `(area·e^ν)·e^Λ`, which fires.

The material detector for the cached measure is `passive_cooling_regression`, which compares
against a golden artifact.

## 8. Deferred — governed by ADR-0004, NOT yet conformed

These still carry their own inline forms and their own degenerate behavior. **They were not
touched**, and this record does not describe them as complete.

| Site | Why deferred |
|---|---|
| **TOV Path 1** — `NStar::Append` Λ block (`NStar.cpp:726` clamp retained), `NStar::FinalizeSurface` baryon integrand | no equivalence coverage; `TOVSolver::Solve` has three live callers (`spin_therm_evol_main`, `tov_debug_main`, `sig_omega`). Phase 3E-0 must build the harness first. |
| **`NStar::BaryonNumIntegrand(double)`** scalar accessor | carries a **separate INV-14 defect** (missing `1e54`) and has **zero callers**. Deliberately **not repaired** here — repairing it would silently convert its output units inside a structural increment. |
| **`MixedStar`** — six sites, two build paths | COMPILED but UNEXERCISED, zero tests, two-sector mass semantics no detector protects (ADR-0004 §0-Q2, §15). Migration waits on focused coverage; superficial tests must not be manufactured to unlock it. |
| **Candidate scientific code** — `DarkCore_Analysis`, `BNV_*`, `Decay_Analysis` | `GOVERNANCE.md` §5 protection; contract only. |
| `NStar.cpp:845` — `EvaluateNu` boundary condition | not this measure (ADR-0004 §4.4). |

`SurfaceGravity`, the Hartle ODE coefficients and the TOV mass ODE share the token `4πr²` but are
**not** this measure and were not touched.

## 9. Invariant disposition

- **INV-03** — headline unchanged (`VERIFIED CURRENT BEHAVIOR`); the metric convention
  `Λ = −½ ln(1 − 2m/r)` is unchanged. Its former clamp wording is **superseded by ADR-0004**
  (rank 2). The entry now distinguishes three statuses explicitly: the governed contract,
  canonical-path conformance, and legacy nonconformance still deferred.
- **INV-04** — `GOVERNED (ACCEPTED) — CANONICAL VALIDATED PATH CONFORMED; LEGACY MIGRATIONS
  DEFERRED`. **Not marked resolved.** The three ownership roles and the four deferred categories
  are recorded there.
- **INV-14** — untouched. The `1e54` conversion stays with the consumer; the scalar accessor's
  defect is untouched and still separately governed.

## 10. Validation summary

| Check | Result |
|---|---|
| Authenticated suite | **16/16 passed** |
| Self-contained suite | **10/10 passed** |
| Protected artifacts | **5/5 byte-identical** |
| Worst `\|ΔB\|/B` | **1.368e-16** (bound `1.0e-15`, predeclared, unwidened) |
| `GeometryCache` arrays, both Λ routes | **bitwise 96/96** |
| Detectors D1, D2, D2b, D3, D5, D6 | **all fire**; D4 correctly not executed |
| Detector reverts | **byte-identical by SHA-256** |
| Tests deleted / tolerances widened | **none** |
| `Core → Physics/Evolution` dependency | **not created** |
| ADR-0003 provenance contract | **unchanged**; no `Refresh()` added |
| Canonical `1e-15` clamp | **gone** from the migrated path; retained in unmigrated Path 1 |
