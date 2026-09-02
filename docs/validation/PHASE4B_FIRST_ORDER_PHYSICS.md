# Phase 4B — independent physical validation of the normalized first-order Hartle response

> **FORMAL STATUS: `PHYSICAL FIRST-ORDER HARTLE RESPONSE VERIFIED`**
>
> The seed-free response Phase 4A exposed is not merely contract-conformant — its **shape** is
> physically right. Production's `omega_bar(r)/Omega` and `[d omega_bar/dr](r)/Omega` agree with
> an independently derived, independently normalized profile at **every node of every validated
> star**, satisfy the exterior-matching and volume-integral identities, and reproduce two
> **derived** weak-field coefficients.
>
> **No production source was changed.** The diff is tests and documentation only. The O(Ω²)
> candidate is byte-identical and remains an unverified candidate (INV-08); nothing here
> validates, executes for a result, or claims anything about it.

| Field | Value |
|---|---|
| **Starting HEAD** | `36ea3b3f385f0b9233e35162565d3662566b1349` (Phase 4A) |
| **Canonical `master`** | `df859b5a73c4cac0c115f240744d89ce9f830b8d` — 3 ahead / 0 behind at entry, upstream equal, clean tree |
| **Change class** | scientific / numerical **validation** (`GOVERNANCE.md` §2) — no production behavior changes |
| **Governing authority** | ADR-0006 (ACCEPTED); INV-07; INV-05; Hartle (1967) ApJ **150**, 1005 |
| **Builds on** | `PHASE4A_FIRST_ORDER_NORMALIZATION.md` (contract), `HARTLE_MOMENT_INERTIA.md` (2B-4B, `EQUATION MATCH`, `HARTLE-I VERIFIED`) |
| **Baseline at entry** | **21/21** authenticated (206.96 s), **11/11** self-contained (14.34 s); seven artifact hashes unchanged |

---

## 1. Exact scope — what 4B claims that 4A did not

Phase 4A proved the **contract**: the arbitrary seed does not leak, a requested spin comes back,
the units are what they say, `J = IΩ`, zero spin is well defined, the response scales linearly.
**Every one of those checks is satisfiable by a profile of the wrong shape.** They constrain the
surface values, the scaling behavior and the seed independence — not the interior.

Phase 4B asks the different question: **is the shape physically correct?** The object under test is

```
s(r)  = omega_bar(r) / Omega          [dimensionless]
s'(r) = [d omega_bar/dr](r) / Omega   [km^-1]
```

with the frame-dragging fraction `omega(r)/Omega = 1 - s(r)`. These carry the entire physical
content of the first-order solution independent of spin amplitude, and they are exactly what
Phase 5 will consume.

The D1 detector (§10) makes the distinction concrete: a corruption of the interior shape that
leaves `I` and the surface untouched passes **every** Phase-2B and Phase-4A test and is caught
only here.

## 2. The independent solver, and why it is still independent

`tests/rotation/hartle_reference.hpp` integrates the **conservative** Hartle system in different
state variables from production:

```
q = r^4 j omega_bar' ,   j = exp[-(nu + lambda)]
d(omega_bar)/dr = q / (r^4 j)
dq/dr           = 16 pi r^4 (eps + p) exp(lambda - nu) omega_bar
```

Its right-hand side is built from the metric columns `nu`, `lambda` and `(eps + p)` — never from
production's `1/(1 - 2m/r)` coefficient helpers. It does not call `ODE_N_Fast`,
`GetHartleOmegaCoeff_N_Fast`, `GetHartleDOmegaCoeff_N_Fast` or `HartleFirstOrderResponse::At`.

**Phase-4B extension, purely additive.** The solver now also records the derivative profile,
recovered from its own conserved flux as `omega_bar'(r) = q(r)/(r^4 j(r))` — a by-product of the
same integration, not a finite difference of the solution — and normalizes itself:

```
Omega_ref = omega_bar_ref(R) + R omega_bar_ref'(R)/3        (its OWN surface extraction)
s_ref(r)  = omega_bar_ref(r)  / Omega_ref
s'_ref(r) = omega_bar_ref'(r) / Omega_ref
```

**Production's normalization factor is never used.** The two shapes being compared are separately
computed and separately normalized.

The extension changed no integration, no tolerance and no existing output: with it in place the
suite is **21/21 and 11/11** and `hartle_I_dscmf1_debug.tsv` is unchanged, which is a stronger
statement than inspection.

## 3. Predeclared acceptance bounds

Fixed **before** the comparison, from the published Phase-2B record only, by one stated rule:

> **profile bound = 5 × (Phase-2B production/reference `I` agreement for that star class),
> rounded up to the next power of ten.**

The rule is tied to `I` because `I` is a functional of exactly this profile, so the two cannot
disagree in order of magnitude; the factor 5 allows a worst **node** to exceed an **integrated**
quantity's error.

| Quantity | Phase-2B input | Derived | **Bound** |
|---|---|---|---|
| analytic profile `s`, `s'` | `9.455e-9` (2B §10) | `4.7e-8` | **`1e-7`** |
| CMF profile `s`, `s'` | `1.32e-5 … 2.15e-5` (2B §11) | `1.08e-4` | **`1e-4`** |
| volume identity, analytic | `1.24e-7` surf/vol + `9.455e-9` | `6.7e-7` | **`1e-6`** |
| volume identity, CMF | `8.05e-7` surf/vol + `2.15e-5` | `1.1e-4` | **`1e-4`** |
| exterior identity (internal) | algebraic — roundoff only | — | **`1e-12`** |
| weak-field coefficient at `M/R = 0.002` | corrections are `O(M/R)` | few `1e-3` | **`5e-3`** |
| reference admissibility | floor must be ≥3 decades below what it measures | — | **`1e-3` × measured** |

**No bound depends on a Phase-4B measurement, and none was widened.** Scratch measurement did
precede the fixing of these numbers; the rule above is stated so a reader can verify that each
bound follows from prior published evidence alone.

## 4. Reference admissibility — the oracle's own numerical floor

Tightening the reference's tolerance by two decades and moving its seed over six, at the **same
centre convention production uses** (both begin at the first profile node):

| Variation | analytic `s` | analytic `s'` | 1.6 M☉ `s` | 1.6 M☉ `s'` |
|---|---|---|---|---|
| tol `1e-13 / 1e-15` | `4.87e-15` | `1.05e-15` | `0` | `3.44e-17` |
| tol `1e-14 / 1e-16` | `6.17e-15` | `1.26e-15` | `1.98e-16` | `5.74e-17` |
| seed × `1e3` | `3.58e-15` | `9.46e-16` | `3.34e-15` | `3.12e-15` |
| seed × `1e-3` | `2.46e-13` | `1.68e-15` | `2.84e-15` | `2.02e-15` |
| **floor** | **`2.46e-13`** | **`1.68e-15`** | **`3.34e-15`** | **`3.12e-15`** |
| measured production/reference | `2.94e-9` | `9.46e-9` | `9.38e-6` | `2.28e-5` |
| **ratio floor/measured** | **`8.4e-5`** | **`1.8e-7`** | **`3.6e-10`** | **`1.4e-10`** |

The reference's own error is four to ten decades below the difference it is used to measure. It
is admissible as the oracle.

**A separate sensitivity, recorded because it is real.** Starting the reference at `r[1]` or
`r[3]` instead of `r[0]` moves `s'` near the centre by up to `4.9e-4` (analytic) and `1.3e-5`
(CMF). That is **not** a tolerance effect and not a disagreement between production and the
reference — it is the sensitivity of the near-centre derivative to where the integration begins,
and production and the reference deliberately share the same convention (first profile node,
`omega_bar' = 0`). It is reported so the shared modelling choice is visible rather than hidden;
INV-05 already records that neither solver uses a series expansion there.

## 5. Experiment A — analytic star, full normalized profile

Exact Schwarzschild constant-density interior, `M = 2 km`, `R = 13 km` (`M/R = 0.154`), 4001
nodes — the Phase-2B-4B fixture.

| Quantity | max relative | max absolute | RMS | worst radius |
|---|---|---|---|---|
| `s(r)` | **`2.936e-9`** | `2.008e-9` | `1.821e-9` | `3.26e-3 km` |
| `s'(r)` | **`9.455e-9`** (also `9.455e-9` as max abs / peak) | — | `1.331e-10` | `13.0 km` |

Bound `1e-7`; margin ≈ 34× on `s` and ≈ 11× on `s'`. The RMS is within a factor 1.6 of the
maximum, so the agreement is a property of the whole profile and not of a few lucky nodes.

`I_prod = 1.5713287051e2`, `I_ref = 1.5713286903e2`, relative `9.455e-9` — **exactly** the Phase
2B-4B figure, and exactly the worst-node `s'` discrepancy. That coincidence is structural: `I`
is fixed by the surface derivative, so the integrated observable and the worst node of the shape
carry the same error.

## 6. Experiment B — DS(CMF)-1 sequence, full normalized profile

| `M` [M☉] | `R` [km] | nodes | max rel `s` | RMS `s` | worst `r` [km] | max rel `s'` | RMS `s'` |
|---|---|---|---|---|---|---|---|
| 1.0000 | 13.4263 | 2629 | `3.415e-6` | `2.328e-6` | 4.795 | `1.704e-5` | `1.289e-7` |
| 1.4000 | 13.5453 | 2646 | `5.845e-6` | `3.405e-6` | 5.306 | `1.603e-5` | `2.184e-7` |
| 1.6000 | 13.4683 | 2635 | `9.375e-6` | `5.116e-6` | 6.160 | `2.276e-5` | `2.789e-7` |
| 2.0001 | 12.7123 | 2527 | **`1.848e-5`** | `7.412e-6` | 6.419 | **`2.244e-5`** | `5.085e-7` |

Bound `1e-4`; worst margin ≈ 4.4×. Every node of every star, not just the surface.

`I` agreement — `1.613e-5`, `1.322e-5`, `2.146e-5`, `1.968e-5` — reproduces the Phase-2B values
`1.61e-5`, `1.32e-5`, `2.15e-5`, `1.97e-5` to the digits published there. The profile discrepancy
tracks the already-characterized `I` discrepancy, which Phase 2B attributed to the tabulated TOV
background rather than to the Hartle solver and showed to be resolution-independent. **No
unexplained discrepancy remains.**

## 7. Experiment C — exterior-matching identities

Outside the star `omega_bar = Omega - 2J/r^3` with `J = IΩ`, so dividing by `Ω`:

```
s(R) = 1 - 2I/R^3 ,        s'(R) = 6I/R^4
```

**Two forms, deliberately separated.**

| Form | What it tests | Result | Bound |
|---|---|---|---|
| production `s(R)`, `s'(R)` vs **production** `I` | Given production's own definitions this is an **algebraic identity**. It is *internal consistency*, not independent physics — the check that the published profile and the published `I` describe the same solution | `0` … `1.44e-16` (all stars) | `1e-12` |
| production `s(R)`, `s'(R)` vs the **independent** `I_ref` | Genuinely independent: nothing on the two sides shares an origin | analytic `1.58e-9` (`s`), `9.46e-9` (`s'`); CMF worst `2.146e-5` | `1e-7` / `1e-4` |

Recording the first as a tautology is deliberate. It is still worth asserting — the D1 detector
shows a corruption that breaks it — but it must not be presented as independent evidence.

## 8. Experiment D — volume-integral identity from production's own shape

```
I_vol = (8 pi / 3) int_0^R r^4 (eps + p) e^{lambda - nu} [omega_bar(r)/Omega] dr
```

derived by integrating the conservative form (`HARTLE_MOMENT_INERTIA.md` §8). It reads the
**interior** of the shape and uses no `J_raw`, no `Omega_raw` and no seed — and since Phase 4A
made `omega_bar/Omega` the actual production response object, it now tests production's own
published profile rather than a reference copy.

| Star | `I_vol` [km³] | `I` [km³] | relative | bound |
|---|---|---|---|---|
| analytic | `1.5713288819e2` | `1.5713287051e2` | `1.125e-7` | `1e-6` |
| CMF 1.0 M☉ | `8.69940852e1` | `8.69957583e1` | `1.923e-5` | `1e-4` |
| CMF 1.4 M☉ | `1.35613750e2` | `1.35616131e2` | `1.756e-5` | `1e-4` |
| CMF 1.6 M☉ | `1.59582613e2` | `1.59587141e2` | `2.837e-5` | `1e-4` |
| CMF 2.0 M☉ | `1.93717087e2` | `1.93723144e2` | `3.127e-5` | `1e-4` |

## 9. Experiment E — weak-field profile limit, with **derived** coefficients

Phase 2B established the integrated limit `I/(MR²) → 2/5`. Phase 4B adds the profile-level
statement, and derives the coefficients rather than fitting them.

In the Newtonian limit `j → 1` and `j' → −4πrρ`, so the frame-dragging equation becomes
`omega_bar'' + (4/r) omega_bar' = k² omega_bar` with `k² = 16πρ = 12M/R³` for a uniform sphere.
Expanding the regular solution for `kR ≪ 1`:

```
omega_bar(r) = omega_bar_c [ 1 + k^2 r^2/10 + O(k^4 r^4) ]
  =>  omega_bar(R) = omega_bar_c (1 + 1.2 M/R) ,   R omega_bar'(R)/3 = omega_bar_c (0.8 M/R)
  =>  Omega = omega_bar_c (1 + 2 M/R)
  =>  omega(0)/Omega -> 2 (M/R) ,        omega(R)/Omega -> 0.8 (M/R)
```

The surface coefficient is independently the exterior result `2I/R³` with `I → (2/5)MR²`.

| `M/R` | `[1−s(0)]/(M/R)` | deviation from **2** | `[1−s(R)]/(M/R)` | deviation from **0.8** |
|---|---|---|---|---|
| 0.150 | 2.05287472 | `2.644e-2` | 0.92578494 | `1.572e-1` |
| 0.100 | 2.03325180 | `1.663e-2` | 0.87793234 | `9.742e-2` |
| 0.050 | 2.01576005 | `7.880e-3` | 0.83645245 | `4.557e-2` |
| 0.020 | 2.00611710 | `3.059e-3` | 0.81404636 | `1.756e-2` |
| 0.010 | 2.00302892 | `1.514e-3` | 0.80693902 | `8.674e-3` |
| 0.005 | 2.00150719 | `7.536e-4` | 0.80344890 | `4.311e-3` |
| 0.002 | 2.00060115 | **`3.006e-4`** | 0.80137467 | **`1.718e-3`** |

Both deviations fall **monotonically** and essentially linearly in `M/R` — the expected leading
correction, not numerical noise — and both are inside the predeclared `5e-3` at the weakest
field. Frame dragging vanishes as `2(M/R) → 0`: no dragging without relativity.

This is a genuinely independent physical check. It uses no reference solver at all; the
prediction comes from the equation's Newtonian limit in closed form.

## 10. Detector D1 — profile-shape corruption

One mutation, at a unique site: production's normalized response multiplied by `1 + 1e-3` at
**interior nodes only**, leaving the centre node, the surface node and `I` untouched.

| Suite | Expected | Observed |
|---|---|---|
| **4B analytic** (new coverage) | FAIL | **FAIL** — `Aa` max rel `1.000e-3` vs bound `1e-7`; `Ac` RMS `7.383e-4`; `Da` volume identity `9.994e-4` vs `1e-6` |
| **4B CMF** (new coverage) | FAIL | **FAIL** — `Ba` worst `9.987e-4` vs `1e-4`; `Da` `9.824e-4` vs `1e-4` |
| Phase-2B `hartle_moment_inertia_analytic` | PASS | **PASS** |
| Phase-2B `hartle_moment_inertia_cmf` | PASS | **PASS** |
| Phase-4A `hartle_normalization_contract` | PASS | **PASS** |
| Phase-4A `hartle_normalization_cmf` | PASS | **PASS** |

**This is the proof that Phase 4B adds coverage.** A wrong interior shape is invisible to the
validated `I` observable and to every contract check, and is caught only by the full-profile
comparison and the volume identity. Reverted byte-identically (SHA-256 verified); all six suites
green afterwards.

A second detector was not added: it would protect no distinct failure mode that D1, the exterior
identities and the volume identity do not already cover.

## 11. Secondary checks

**Spin reversal (Experiment F).** Materializing at `+100`, `−100` and `+2π·100 rad/s`:
`Ω(−) = −Ω(+)`, `J(−) = −J(+)`, and both profiles negate **bit-exactly** at all 4001 analytic and
2646 real-star nodes; `I` is unchanged; the seed-free response is bit-identical after three
materializations. These are exact consequences of first-order linearity and are recorded as
secondary — they do not substitute for Experiments A–E.

**Frame-dragging fraction (diagnostic).**

| Star | `ω(0)/Ω` | `ω(R)/Ω` |
|---|---|---|
| analytic `M/R = 0.154` | `3.161e-1` | `1.430e-1` |
| CMF 1.0 M☉ | `2.753e-1` | `7.189e-2` |
| CMF 1.4 M☉ | `3.800e-1` | `1.091e-1` |
| CMF 1.6 M☉ | `4.380e-1` | `1.306e-1` |
| CMF 2.0 M☉ | `5.975e-1` | `1.886e-1` |

`0 < ω/Ω < 1` holds on every validated star and dragging decreases outward. **This is asserted as
a measured ordering, not as a derived inequality** — no governing document derives a hard bound,
and this increment does not invent one from intuition.

## 12. Slow-rotation regime — an explicit non-claim

| Star | `Ω_K = (M/R³)^{1/2}` [km⁻¹] | `Ω/Ω_K` at `2π·716 Hz` |
|---|---|---|
| analytic | `3.017e-2` | `0.497` |
| CMF 1.0 M☉ | `2.470e-2` | `0.608` |
| CMF 1.4 M☉ | `2.884e-2` | `0.520` |
| CMF 1.6 M☉ | `3.110e-2` | `0.483` |
| CMF 2.0 M☉ | `3.792e-2` | `0.396` |

> **Normalization correctness is not slow-rotation truncation accuracy.**

The first-order problem is linear, so the normalized response is exact at any `Ω` — that is what
§5–§9 verify. Whether *truncating the stellar structure at Hartle order* is physically adequate
at `Ω/Ω_K ≈ 0.6` is a completely separate question about the size of the neglected O(Ω²) and
higher terms. **Nothing in this record bears on it.** It belongs with the O(Ω²) validation and
regime work (ADR-0006 §10, INV-08).

## 13. Test-seam classification (wording correction, no production change)

The Phase-4A records describe `RotationSolverTestSeam` as reachable only by the validation
harnesses. In C++ that is **too strong**: a friend declaration names a type, and any translation
unit may define that type and gain the same access. The accurate classification is:

> **`RotationSolverTestSeam` = PRIVILEGED TEST BACKDOOR — NOT SUPPORTED SCIENTIFIC API.**

The invariant ADR-0006 Q2 actually requires is unchanged and still holds:

- no supported public seed setter;
- no supported public seed constructor argument;
- no production consumer of the seam.

All three verified at this HEAD. This is a documentation correction only; **no production change
was made for it**, and it is not a Phase-4B scientific finding.

## 14. O(Ω²) candidate — frozen and confirmed

`ODE_Hartle2_N_Fast`, `SolveHartle2_N`, `ODE_Hartle2_Mixed_Fast` and `SolveHartle2_Mixed` are
**byte-identical** to the Phase-4A source. The production diff for this increment is **empty**.
The candidate was not executed for a result, no candidate baseline was created, and **INV-08's
scientific status is unchanged**.

## 15. Artifacts and final suite

No artifact was regenerated or re-baselined. All seven Phase-3 goldens are byte-identical, and
`hartle_I_dscmf1_debug.tsv` re-emitted from the rebuilt binary compares byte-for-byte
(`ddf018579364e9b3bf24f1e3a3e2577e70fb09b8b84eb718237d9da683aa9d15`).

The Hartle golden remains an **`I` reference**, not a normalized-profile golden. No profile
golden was created: the comparison is against a solver that runs in the same process, which is
strictly stronger than freezing a radial table — a frozen table could only detect change, while
the independent solver detects *wrongness*.

| Configuration | Registered | Result | Wall time | Before 4B |
|---|---|---|---|---|
| Full, authenticated data root | **23** | **23/23 PASS** | 214.43 s | 21/21, 206.96 s |
| Self-contained | **12** | **12/12 PASS** | 15.95 s | 11/11, 14.34 s |

Both counts are up by one: `hartle_first_order_physics_analytic` (self-contained) and
`hartle_first_order_physics_cmf` (external data). The `+7.5 s` and `+1.6 s` are those two tests;
nothing pre-existing slowed, and no production code runs differently.

## 16. The exact boundary of the claim

**VERIFIED.**

- The normalized first-order response `s(r)`, `s'(r)` is correct as a *shape*, at every node of
  the analytic star and of the four authenticated CMF stars, against an independently derived
  and independently normalized profile.
- It satisfies the exterior-matching identities against an independent `I`.
- It reproduces `I` through a volume integral that reads only the interior.
- It reproduces two independently derived weak-field coefficients.
- The oracle used is numerically admissible by four to ten decades.

**NOT claimed.**

- Nothing at O(Ω²). INV-08 is untouched.
- No statement about the adequacy of the slow-rotation truncation at any `Ω` (§12).
- No hard inequality for `ω/Ω` — that ordering is diagnostic (§11).
- No claim from empirical universal relations. Breu–Rezzolla and Lattimer–Schutz remain the
  Phase-2B secondary sanity checks and are **not** used as an oracle here.
- No statement about `MixedStar` two-fluid rotation, which remains on its own track.
- The residual production/reference difference on real stars (`≤ 2.3e-5`) is inherited from the
  tabulated TOV background, as Phase 2B established; 4B does not re-open the background's own
  accuracy, which is the Phase-2B-4A grid-convergence question.

## 17. Status

> **`PHYSICAL FIRST-ORDER HARTLE RESPONSE VERIFIED`**

INV-07 becomes **GOVERNED (ADR-0006 ACCEPTED) — FIRST-ORDER PHYSICAL NORMALIZATION CONFORMED AND
PHYSICAL RESPONSE INDEPENDENTLY VERIFIED**. First-order normalization and response work is
complete. **This implies nothing about O(Ω²)**, which remains an unverified candidate.

**Exact next task.** *Phase 4C-G — derive and propose the separate O(Ω²) Hartle monopole ADR from
primary equations, replacing the existing AI candidate as scientific authority rather than
patching it incrementally.*
