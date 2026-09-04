# ADR-0007 — Hartle O(Ω²) monopole structural response at fixed central energy density

| Field | Value |
|---|---|
| **Status** | **ACCEPTED** |
| **Date accepted** | 2026-09-02 — by project-owner adjudication |
| **Date drafted** | 2026-09-02 |
| **Drafted at** | `bb073c8ed0c7ce15f3a8c960e9f76173bde51a39` (branch `physics/rotation-correctness`) |
| **Accepted at** | `bcef5b57e3267fe26d1269efaa472d314818d367` (acceptance commit; implementation follows separately as 4C-I0 / 4C-I1) |
| **Change class** | **scientific-semantic** (equations, perturbation variable, boundary condition, `dε/dp` ownership, `δM` definition, normalization) **and structural** (public O(Ω²) API and result type replaced) — strictest applies (`GOVERNANCE.md` §2). Implementation would proceed under `GOVERNANCE.md` **§3.1** (recorded in §9) |
| **Governing authority** | **Primary:** J. B. Hartle, ApJ 150, 1005 (1967) — eqs. (87)–(88), (97), (99)–(100), (105)–(108), (109), p. 1009, p. 1022. **Secondary:** Hartle & Thorne, ApJ 153, 807 (1968) §II — (7a–c), (13)–(17b), (24b). Repository: ADR-0006 (ACCEPTED; first-order normalization, no implicit spin), ADR-0001 (`n_i = Y_i n_B`), ADR-0004 (`w_V = 4πr²e^{λ}`), ADR-0005 (`I` and the goldens bitwise), INV-02, INV-05, INV-06, INV-07, INV-13, INV-14 |
| **Evidence record** | `docs/validation/PHASE4C_HARTLE2_DERIVATION.md` (Phase 4C-G) — every equation, unit, bound and line citation below is derived or cited there |
| **Affected invariants** | **INV-08** (primary), **INV-09** (the fixed-`ε_c` response Phase 5 consumes), INV-06 (surface semantics for O(Ω²) quantities — no convention change), INV-14 (O(Ω²) extension of the baryon integral) |
| **Blocks** | roadmap 4C-I0 (EOS derivative authority), 4C-I1 (monopole implementation), 4D (validation), 4E (Phase-5 structural fields); transitively Phase 5 |
| **Explicitly does not decide** | the `l = 2` sector (shape, quadrupole) beyond its scope classification (Q11); the chemical-imbalance convention (INV-11); MixedStar rotation; the solar-mass / `G` authority; any change to the TOV surface convention |

---

## 1. Context

The repository carries an AI-authored O(Ω²) monopole candidate (`675b4a9`; `ODE_Hartle2_N_Fast`
`CompactStar/Core/src/RotationSolver.cpp:1140-1226`, `SolveHartle2_N` `:1242-1386`,
`HartleResult` `RotationSolver.hpp:203-217`) that is **publicly callable, has zero repository
callers, and is unvalidated** (INV-08; `PHASE4_ROTATION_ENTRY.md` §13). Against the primary
source it is defective in every part that matters (evidence record §16): three of its seven ODE
terms are dimensionally inconsistent with their own equation; its `j²` factor is inverted and
lacks `e^{−2ν}`; its second frame-dragging source has the wrong sign, a factor 2 and one power
of `r` missing; it integrates a variable (`p0 ≡ δp₀`) whose evolution equation is correct for
neither `δp₀` nor Hartle's `p₀*`; it imposes `δp₀(R) = 0` by shooting — not Hartle's condition —
which forces `ξ₀(R) = 0` and shifts the central density; its `δM` omits `J²/R³`; it takes `dε/dp`
from profile finite differences with a `1.0` fallback; and every output scales as the square of
the arbitrary first-order seed (measured), in breach of ADR-0006 P9. Its own comments cite
H67 "(30)–(33)", which are not the `l = 0` equations, and carry `[FIX: confirm exact from
textbook]` and `???` (`:1135`, `:1207`).

No baseline of its output exists and none may be created: freezing these numbers would make
adjudicated-invalid physics the reference (`GOVERNANCE.md` §3.1 rationale).

Phase 4A/4B settled first order: the canonical seed-free response `I`, `s(r) = ω̄/Ω`,
`s'(r) = ω̄'/Ω` is published by `HartleFirstOrderResponse` (`RotationSolver.hpp:159`) and
independently verified (`PHASE4B_FIRST_ORDER_PHYSICS.md`). Phase 5 (INV-09) needs
`A_i = (∂N_i/∂Ω²)|_{ε_c}` — the O(Ω²) monopole response **at fixed central density, per unit
`Ω²`**. This ADR proposes the contract that replaces the candidate as scientific authority. **It
writes no code and authorizes no source change until accepted.**

## 2. Decision questions

| | Question |
|---|---|
| Q1 | canonical integrated pressure variable |
| Q2 | fixed-`ε_c` centre condition |
| Q3 | canonical `l = 0` equations |
| Q4 | regular-centre numerical start |
| Q5 | EOS derivative authority |
| Q6 | surface / matching convention |
| Q7 | `δM` definition |
| Q8 | `ξ₀` definition |
| Q9 | `Ω²` normalization / coefficient representation |
| Q10 | public result / API semantics |
| Q11 | `l = 0` vs `l = 2` scope |
| Q12 | §3.1 correction / baseline policy |
| Q13 | disposition of the current public candidate API |

## 3. Definitions proposed (normative once accepted)

Conventions: CompactStar metric `g_tt = −e^{2ν}`, `g_rr = e^{2λ} = (1 − 2m/r)^{−1}` (INV-03);
Hartle's `ν_H = 2ν`, `λ_H = 2λ`; geometric km (INV-02); `Ω ≡ Ω_geom = Ω_phys/c` (ADR-0006).

- **`j(r) ≡ e^{−ν−λ} = e^{−ν}(1 − 2m/r)^{1/2}`** (H67 40); `j² = e^{−2ν}(1 − 2m/r)`;
  `dj²/dr = −8πr(ε+p)e^{−2ν}`; `j²/(r − 2m) = e^{−2ν}/r`.
- **`m₀(r)`** [km]: the O(Ω²) `l = 0` perturbation of the enclosed-mass function,
  `g_rr → e^{2λ}[1 + 2m₀/(r − 2m)]` (H67 66–68).
- **`p₀*(r)`** [dimensionless]: Hartle's pressure perturbation factor,
  `p₀* ≡ −ξ₀ (dp/dr)/(ε+p) = ν' ξ₀` (H67 87–88, 99). Equivalently the Eulerian fractional
  enthalpy perturbation.
- **`δp₀(r) ≡ (ε+p) p₀*`** [km⁻²]: the Eulerian monopole pressure perturbation at fixed `r`
  (H67 79–80; HT68 7a); `δε₀ = (ε+p)(dε/dp)p₀*`; `δn_i = (ε+p)(dn_i/dp)p₀*`.
- **`Δp₀ ≡ δp₀ + ξ₀ dp/dr = 0`** identically (Lagrangian change vanishes by construction of the
  constant-density labelling, H67 6, 71). "`δp₀(R) = 0`" is therefore not a boundary condition.
- **`ξ₀(r) ≡ p₀*/ν' = −δp₀/(dp/dr)`** [km]: outward displacement of the isobar at `r`.
- **`R_*`, `ε_*`, `p_*`, `n_*`**: the production surface node and its values (INV-06: the last
  node, at `p ≤ max(1e-15 p_c, eos_tab.pre[0])`). `R_*` is **not** identified with the `p = 0`
  surface; every surface quantity is defined at `R_*` (P7).
- **Normalized coefficient** of any O(Ω²) quantity `Q`: `Q̂ ≡ Q/Ω_geom²`. Units: `m̂₀` km³,
  `p̂₀*` km², `δp̂₀` 1, `ξ̂₀` km³, `δM̂` km³, `δN̂_i` km².
- **Fixed-`ε_c` family**: the solution with `m₀(0) = 0`, `p₀*(0) = 0` (H67 108, p. 1009; HT68 §II f)
  — the rotating star has the same central energy density as the non-rotating one.

## 4. Proposed contract

- **P1 — Canonical variable (Q1).** The integrated pressure variable is Hartle's `p₀*`.
  `δp₀`, `ξ₀`, `δε₀` are derived fields. No other pressure-like variable is integrated.
- **P2 — Governed equations (Q3).** The `l = 0` system is H67 (97), (100) — HT68 (15a–b) —
  written per unit `Ω_geom²` with the verified first-order shape `s`, `s'` and every metric
  factor explicit:

  ```
  dm̂₀/dr  = 4π r² (ε+p) (dε/dp) p̂₀*  +  (1/12) r⁴ e^{−2ν} (1 − 2m/r) s'²  +  (8π/3) r⁴ (ε+p) e^{−2ν} s²
  dp̂₀*/dr = − m̂₀ (1 + 8π r² p)/(r − 2m)²  −  4π (ε+p) r² p̂₀*/(r − 2m)
             + (1/12) r³ e^{−2ν} s'²  +  (2/3) r e^{−2ν} s ( s + r s' − r ν' s )
  ```

  (dimensions km² and km term by term; evidence record §6–§7). Inputs: profile `r`, `m`, `p`,
  `ε`, `ν` (`MetricNu`), `ν'` (`MetricNuPrime`); `s`, `s'` from `HartleFirstOrderResponse`;
  `dε/dp` from the EOS (P5). **Forbidden:** `stored_omega_bar_`, `stored_domega_bar_`, the seed,
  any raw `ω̄`, `J`, `Ω`.
- **P3 — Centre condition (Q2).** Fixed `ε_c`: `m̂₀ = p̂₀* = 0` at the centre, no homogeneous
  admixture, **no surface condition of any kind**. The regular homogeneous solution
  (`p̂₀*(0) = c₀`) is the TOV sequence derivative and may be exposed separately (P11); it is never
  mixed into the fixed-`ε_c` response.
- **P4 — Numerical start (Q4).** At the first grid node `r₀` the integration starts from the
  regular series (evidence record §9):
  `p̂₀*(r₀) = (1/3) j(r₀)² s(r₀)² r₀²`,
  `m̂₀(r₀) = (4π/15)(ε+p)(r₀)[(dε/dp)(p_c) + 2] j(r₀)² s(r₀)² r₀⁵`
  (truncation `O(r₀⁴/R⁴) ≈ 4e-25`). Starting from literal zeros is an admissible documented
  approximation (relative error `≈ 6e-13` on `p̂₀*(R_*)`, `ε_c` perturbed at `≈ 1e-12`), not the
  contract.
- **P5 — EOS derivative authority (Q5).** `dε/dp` is the derivative of the **same** tabulated
  EOS interpolant that built the star (`gsl_interp_steffen`, `TOVSolver.hpp:488`), evaluated by
  the EOS/TOV layer with an explicit unit contract (cgs inside `TOVSolver`, dimensionless
  `dε/dp` outside, INV-02), and carried on the profile as a column; the same rule governs
  `dn_i/dp` for Phase 5. **Profile finite differences are not an authority**; a star without an
  EOS-derivative column (point-constructed `NStar`) must be given one explicitly or fail closed.
  A tabulated `c_s²`, where the EOS source provides one, is a cross-check, not a source.
  **Implementation prerequisite:** no such API exists today (evidence record §10).
- **P6 — Exterior matching and `δM` (Q6, Q7).** With `J = IΩ` and the explicit surface shell,

  ```
  δM̂ = m̂₀(R_*⁻) + 4π R_*² ε_* ξ̂₀(R_*) + I²/R_*³      [km³] ,        δM = δM̂ Ω_geom²  [km]
  ```

  (H67 105–107; HT68 15c–16; shell from H67 (18)/(97) at a density discontinuity). The shell
  term is **computed, never assumed zero** (it is `≲ 1e-5 δM̂` on EOS-floor stars and dominant on a
  constant-density star). `m̂₀ + I²/r³` is constant wherever `(ε+p)r²` is negligible, so the
  matching at `R_*` is exact to `O(ε_*/⟨ε⟩) ≈ 1e-6`.
- **P7 — Surface semantics (Q6).** All O(Ω²) surface quantities are evaluated at `R_*`.
  `ξ₀(R_*)` is the displacement of the `p = p_*` isobar; relative to the `p = 0` surface it
  carries a systematic `O(ΔR/R_*) ≈ 4–7 × 10⁻³` (relative), of the same origin as the
  existing `0.20–0.35 %` radius residual (`TOV_REFERENCE.md` §5). No extrapolation, no surface
  reconstruction, no change to INV-06. **Classification: `SURFACE ADEQUATE AS-IS`.**
- **P8 — `ξ₀` (Q8).** `ξ̂₀ = p̂₀*/ν'` [km³] on the whole grid; `ξ̂₀(R_*)` and `δR/R = ξ₀(R_*)/R_*`
  are reported with the P7 caveat; `ξ₀(R_*)` is generally nonzero and positive.
- **P9 — Normalization and representation (Q9).** The public product is the set of
  **seed-free coefficients per `Ω_geom²` at fixed `ε_c`**: `m̂₀(r)`, `p̂₀*(r)`, `δp̂₀(r)`, `ξ̂₀(r)`,
  `δM̂`, `ξ̂₀(R_*)`, the shell term, on the profile grid, one canonical geometric representation,
  no duplicated physical/geometric state. Physical values are `Q = Q̂ Ω_geom²` from an explicit
  `AngularVelocity` through a materialization method; **no new ODE solve per `Ω`**; zero spin
  materializes zeros; ADR-0006's binding clarification extends: **no `NStar` carries an O(Ω²)
  physical perturbation without an explicit physical angular velocity**.
- **P10 — Public semantics (Q10).** A new result type (`HartleMonopoleResponse`, name
  indicative) reached through a public `NStar` accessor; nothing produced by the candidate sits
  behind it; every field carries its unit in the declaration; `I` and the first-order response
  are untouched.
- **P11 — Homogeneous (sequence-derivative) solution.** The regular homogeneous solution per
  unit `p₀*(0)` equals `∂/∂ε_c` along the non-rotating sequence with
  `δε_c = (ε_c+p_c)(dε/dp)_c p₀*(0)`, and serves Phase 5's `B_i` and validation identity K.
  > **MODIFIED AT ACCEPTANCE (owner, 2026-09-02) — see Decision Q8.** The mathematics stands and
  > may be computed **internally or test-side in 4D** to validate the sequence-derivative
  > identity. **It is NOT exposed as a public production API in Phase 4C.** Public ownership of
  > `B_i` / `∂N_i/∂ε_c` is deferred until Phase 5 adjudicates that architecture.
- **P12 — Scope (Q11).** The `l = 0` sector is sufficient for every scalar count (`A`, `N_i`),
  for `M` and for `I` at O(Ω²) — derived from the angular integration (evidence record §14).
  The `l = 2` sector (shape, quadrupole) is outside this ADR; whether the ratified Phase-4 exit
  criterion "O(Ω²) validated" also demands it is decided by the owner (§11 Q6).
- **P13 — Phase-5 interface.** Phase 4 exposes the fields of P9 (and P11); Phase 5 forms
  `A_i = δN̂_i = ∫₀^{R_*} w_V { (dn_i/dp)δp̂₀ + n_i[m̂₀/(r−2m) + (1/3)r²e^{−2ν}s²] } dr + w_V(R_*) n_i(R_*) ξ̂₀(R_*)`
  with `n_i = Y_i n_B` (ADR-0001), `w_V = 4πr²e^{λ}` (ADR-0004), EOS `dn_i/dp` (P5). Constant
  baryon number is **never** imposed inside the Hartle solve (INV-09 `Z_i` is Phase 5).
- **P14 — Regime.** Results carry `Ω/Ω_K`, `Ω_K = (M/R_*³)^{1/2}` (ADR-0006 P10); correctness of
  the coefficients is not accuracy of the slow-rotation truncation.

## 5. Alternatives

### Q1 — integrated variable

| | Option | Assessment |
|---|---|---|
| **A** | **Hartle's `p₀*`** | primary equations most direct; dimensionless; centre condition a literal zero; `r²` series; `δp₀ = (ε+p)p₀*` a well-conditioned product; `ξ₀ = p₀*/ν'`; conservation identity (H67 90) available. **Recommended** |
| B | Eulerian `δp₀` | adds `(1 + dε/dp)ν'δp₀` to the homogeneous part — the EOS derivative enters the *conditioning* of the system where it is largest (crust); `ξ₀` needs a division by `dp/dr → 0`; centre condition `δp₀(0) = 0` equivalent. Kept as the **independent 4D formulation** (line C) |
| C | `(m₀, h₀)` via (97)+(98) | rigorous, but `h₀` is fixed only by exterior matching of a constant (HT68 17b) — reintroduces a surface step for a metric normalization. Second independent 4D formulation; rejected as canonical |

### Q2 — centre condition

| | Option | Assessment |
|---|---|---|
| **A** | **fixed `ε_c`: `m̂₀(0) = p̂₀*(0) = 0`** | Hartle's family; the derivative Phase 5 needs; unique with regularity. **Recommended** |
| B | fixed total baryon number inside the solve | mixes Phase-5's sequence reduction into Phase 4; obtainable afterwards as particular + `α`·homogeneous. Rejected |
| C | the candidate's `δp₀(R) = 0` | a different family member (central density shifted; `ξ₀(R) = 0`); not Hartle's; **rejected** |

### Q4 — numerical start

| | Option | Assessment |
|---|---|---|
| A | literal zeros at `r₀` | error `≈ 6e-13`; admissible approximation |
| **B** | **regular series at `r₀`** | two one-line expressions; error `≈ 4e-25`; makes "fixed `ε_c`" exact. **Recommended** |

### Q5 — `dε/dp`

| | Option | Assessment |
|---|---|---|
| **A** | **derivative of the star's own Steffen interpolant, EOS/TOV-owned, on the profile** | consistent by construction with INV-13; continuous; no second authority. **Recommended** |
| B | tabulated `c_s²` from the EOS source | not imported; can disagree with A at table resolution; cross-check only |
| C | profile finite differences (candidate) | chain-rule equivalent only where `dp/dr ≠ 0`; noisy (`3.8–4.9e3` range on CMF); silent causal-limit fallback. **Rejected** |
| D | smoother re-interpolation for the derivative only | derivative of a different function than the star used; needs INV-13 re-governed. Rejected |

### Q6 / Q7 — surface and `δM`

| | Option | Assessment |
|---|---|---|
| **A** | **evaluate at `R_*`, explicit shell and boundary terms, recorded `O(ΔR/R)` systematic on `ξ₀(R_*)`** | identities exact to `≈ 1e-6`; nothing silently identified with `p = 0`. **Recommended — `SURFACE ADEQUATE AS-IS`** |
| B | extrapolate the profile to `p = 0` before matching | a numerical-method change moving every existing result (`TOV_REFERENCE.md` §5); no O(Ω²) quantity needs it |
| C | block 4C on a surface-convention ADR | not justified by the bounds (evidence record §11) |

### Q9 / Q10 — representation and API

| | Option | Assessment |
|---|---|---|
| **A** | **coefficients per `Ω_geom²` + materialization at explicit `AngularVelocity`** | ADR-0006 P9 route (b); one solve per star; parallels the verified first-order architecture. **Recommended** |
| B | coefficients only | forces every caller to multiply and convert; loses the no-implicit-spin guard |
| C | solve per requested `Ω` from `ω̄_phys` | re-solves a linear system for a provable multiple; rejected |

### Q11 — scope

| | Option | Assessment |
|---|---|---|
| **A** | **`l = 0` now; `l = 2` as a separate later ADR** | sufficient for Phase 5 and for replacing the (`l = 0`-only) candidate; the exit-criterion reading is put to the owner. **Recommended** |
| B | `l = 0` and `l = 2` together in 4C | doubles the derivation and validation surface for quantities no downstream phase consumes yet |

### Q13 — candidate API

| | Option | Assessment |
|---|---|---|
| **A** | **atomic replacement in 4C-I** (delete `SolveHartle2_N`, `ODE_Hartle2_N_Fast`, the MixedStar stubs, `HartleResult`, `GetHartleResult`, `include_m0p0_source_`, `fast_dEdP*`, `fast_nu` in the commit that adds the governed product) | zero callers ⇒ zero churn; the overloaded name retires. **Recommended** |
| B | fail closed until replacement | a production change; not for this task; available as 4C-I-0 if implementation is delayed |
| C | keep the name, replace internals | continuity of wrong semantics; rejected |

## 6. Consequences (once accepted)

- INV-08 moves from `INTENDED BUT UNVERIFIED` to **GOVERNED (ADR-0007 ACCEPTED)**; conformance
  is 4C-I, verification 4D, tracked separately (ADR-0002/ADR-0006 precedent).
- Forbidden: integrating any variable but `p₀*`; any surface condition on `p₀*`; any `dε/dp`
  from profile differences; any O(Ω²) output not per `Ω_geom²` or not from an explicit
  `AngularVelocity`; any `δM` without the `I²/R_*³` and shell terms; identifying `R_*` with the
  `p = 0` surface in any formula or label; exposing candidate output behind the new accessor.
- Required, split across two governed increments (owner sequencing, 2026-09-02):
  **4C-I0** — the EOS-derivative authority of P5 (the `dε/dp` evaluator on the star's own `ε(p)`
  interpolant and its profile-attached delivery), with its own analytic and real-EOS validation;
  **4C-I1** — the governed monopole solver, the new result type and `NStar` accessor, and the
  atomic deletion of the candidate (Q13 A), plus `CURRENT_ARCHITECTURE.md` rows, INV-08/INV-09
  status, and the tests of §7.
- Not changed: the first-order ODE, `I`, `HartleFirstOrderResponse`, the seven goldens, the TOV
  surface convention.
- Historical outputs: every value ever produced by `SolveHartle2_N` / `GetHartleResult`
  (including the diagnostic numbers in `PHASE4_ROTATION_ENTRY.md` §12) is **not a reference
  result** and must not be compared against (§9 item 6).

## 7. Validation requirements (predeclared; evidence record §21 has the full plan)

| # | Requirement | Bound (basis) |
|---|---|---|
| 1 | seed invariance of every coefficient over `[1e-3, 1e3]` (test seam) | rel `≤ 1e-10` (ADR-0006 §7-1) |
| 2 | materialized quantities quadratic in the requested `Ω`; zero spin ⇒ zeros | `≤ 1e-14` |
| 3 | centre series coefficients `(1/3)j_c²s_c²`, `b₅` recovered | first ten nodes `≤ 1e-8` |
| 4 | node-by-node agreement with an independent solver in different variables | analytic `1e-7`, CMF `1e-4` (to be re-predeclared at 4D entry) |
| 5 | exterior identity `m̂₀ + I²/r³` constant over the near-vacuum nodes; `δM̂` matching-node independent | `≤ 1e-6` |
| 6 | conservation identity `p̂₀* + ĥ₀ − (1/3)r²e^{−2ν}s² = const` | `≤ 1e-9` relative |
| 7 | Newtonian homogeneous-star limits `δM̂ → R³`, `3Mξ̂₀/R⁴ → 1`, monotone in `M/R` | intercept `≤ 5e-3` |
| 8 | published HT68 Table 5 models on the printed HW EOS; Chandrasekhar & Miller (1974) homogeneous star — **sources to be authenticated in 4D** | `≤ 2e-2` (HT68's stated accuracy + interpolation differences) |
| 9 | radial convergence at the profile's interpolation order | Richardson residual reported |
| 10 | EOS-derivative sensitivity (Steffen vs `c_s²` column vs FD) | spread `≤ 1e-3` on `δM̂`, FD noise visible |
| 11 | sequence-derivative identity (homogeneous `δM̂_hom` vs `dM/dε_c` sweep) | `≤ 1e-3` |
| 12 | `I`, both Hartle CTests, the seven goldens | **bitwise** (ADR-0005 §13) |

Then the first monopole baseline (§9 item 7).

## 8. Detectors (each applied, measured, reverted byte-identically — later, in 4D)

`M1` drop `e^{−2ν}` from a source → 4, 3 · `M2` flip the `(8π/3)` source sign → 4, 3 · `M3` omit
`I²/R_*³` → 5, 8 · `M4` omit the shell term → 7, 4 (homogeneous star) · `M5` profile-FD `dε/dp` →
10, 9 · `M6` impose `δp₀(R_*) = 0` → 3, 11, 8 · `M7` literal-zero start (weak at `r₀ = 1e-5 km`)
→ 3 only with an enlarged `r₀` · `M8` seed leak → 1, 2 · `M9` drop `(1 + 8πr²p)` → 4, 8.

## 9. `GOVERNANCE.md` §3.1 record (Q12)

Once accepted, this ADR authorizes the pre-baseline correction:

1. **Invalid / superseded behaviour:** `ODE_Hartle2_N_Fast` and `SolveHartle2_N` at `bb073c8`
   (`RotationSolver.cpp:1140-1226`, `:1242-1386`) — dimensionally inconsistent terms, wrong
   `j²` and source terms, non-Hartle surface shooting, incomplete `δM`, profile-FD `dε/dp`,
   seed-quadratic normalization.
2. **Why capturing it is forbidden:** a golden of `δM/M = −1.6`, `ξ₀(R) = 0` and a shifted
   central density would make every later correction register as a regression.
3. **Minimum correction authorized:** the contract of §4 — P1–P10 — and nothing else: no
   `l = 2`, no rotochemical coefficients, no MixedStar, no surface-convention change.
4. **Independent verification required:** §7 items 1–11 (analytic limits, series, independent
   solver, conservation identities, published results, convergence, EOS sensitivity) — none
   compares against candidate output.
5. **Narrow scope:** the O(Ω²) path and its result type; first order and `I` bitwise (§7-12).
6. **Historical candidate outputs are not references:** see §6.
7. **Baseline policy:** created immediately after 4D validation on DS(CMF)-1 (proposed
   `tests/baselines/hartle_monopole_dscmf1_debug.tsv`), never before.

An agent invoking §3.1 under this ADR must restate these seven items in its report.

## 10. Non-scope

No production change in the proposing task; no `l = 2` derivation; no rotochemical physics;
no INV-11; no MixedStar; no change to the TOV surface convention or to INV-06; no cgs `J`/`I`/`δM`
authority; no acceptance.

## 11. Owner questions

| | Question | Recommendation |
|---|---|---|
| Q1 | Canonical pressure variable? | **`p₀*`** (§5-Q1 A) |
| Q2 | Fixed `ε_c` as the governed family? | **Yes** (§5-Q2 A) |
| Q3 | Coefficients per `Ω_geom²`, materialized at an explicit `AngularVelocity`? | **Yes** (§5-Q9 A) |
| Q4 | EOS owns `dε/dp` (and `dn_i/dp`), via the star's own interpolant, as a 4C-I-0 prerequisite? | **Yes** (§5-Q5 A) |
| Q5 | Does the EOS-floor surface require a Phase-4 correction? | **No** — `SURFACE ADEQUATE AS-IS` with P7's semantics (§5-Q6 A) |
| Q6 | Is `l = 0` sufficient for the governed Phase-5 deliverable, with `l = 2` deferred — and is the ratified exit criterion "O(Ω²) validated" satisfied by the monopole sector? | **`l = 0` sufficient (derived)**; recommend reading the exit criterion as the monopole sector and scheduling `l = 2` separately — **owner to confirm** |
| Q7 | Disposition of the public candidate API? | **A — atomic replacement** (§5-Q13) |
| Q8 | Expose the homogeneous (sequence-derivative) solution as an optional second struct? | **Yes** (P11) |

## 12. Implementation record — Phase 4C-I0 (2026-09-02)

**P5 only.** The `dε/dp` authority of P5 is implemented and independently validated; **no other
clause of §4 is implemented**, no O(Ω²) equation was added, and the candidate is byte-identical
to `master`. Evidence: `docs/validation/PHASE4C_I0_EOS_DERIVATIVE.md`.

| P5 requirement | As implemented |
|---|---|
| the EOS/TOV layer owns the evaluation | `TOVSolver::GetEDensDeriv`, with `HasEDensDeriv`, `EDensDerivPressMin/Max` (`TOVSolver.hpp`) |
| the derivative of the **same** `ε(p)` interpolant | differentiates `visi_eps_p_spline` (`gsl_interp_steffen`, `TOVSolver.hpp:488`) through the same accelerator `GetEDens` uses; no second interpolant exists |
| delivered dimensionless | `× (INV_FM4_2_Dyn_CM2 / INV_FM4_2_G_CM3)`, derived from the same two constants `NStar::BuildFromTOV` uses for `ε` and `p`, and asserted against an independent literal `c = 2.99792458e10 cm/s` to `1.458e-16` |
| never a profile finite difference, no `1.0` fallback | the only finite differences in the tree are a test oracle and a printed diagnostic; detector D1 proves the fail-closed contract catches a substitution |
| carried to the profile so `RotationSolver` needs no EOS object | `StarProfile::HasEosDEdP()` / `GetEosDEdP()`, filled through `TOVPoint::dedp` by both ordinary-`NStar` construction paths. `RotationSolver.{hpp,cpp}` were **not touched** |
| absence fails closed | out-of-domain / non-finite / no-interpolant ⇒ NaN (never a clamped boundary value, never `0.0`, which is the physical value for incompressible matter); a profile publishes the set only if every node has a finite value |
| point-constructed stars may supply their own | `TOVPoint`'s trailing `dedp` (default NaN = "not supplied"), distinct from an explicit `0.0` |

Validated by `tests/eos/eos_derivative_contract.cpp` (15/15, self-contained, bounds predeclared
from the interpolation order) and `tests/eos/eos_derivative_cmf.cpp` (17/17, DS(CMF)-1), with
detectors D1–D3 fired and reverted byte-identically. Suites: **25/25** and **13/13**; the seven
durable artifacts and `I` are byte-identical. `dn_i/dp` remains **GOVERNED PRINCIPLE
ESTABLISHED — IMPLEMENTATION DEFERRED TO PHASE 5**. `GOVERNANCE.md` §3.1 is **AUTHORIZED, NOT
YET EXERCISED**: 4C-I0 changed no scientific output and created no monopole baseline.

## 13. Implementation record — Phase 4C-I1 (2026-09-03)

**The §3.1 correction is EXECUTED.** The candidate of `675b4a9` — `SolveHartle2_N`,
`ODE_Hartle2_N_Fast`, the two MixedStar stubs, `HartleResult`, `GetHartleResult` and every
proven candidate-only member — is **deleted**, and the governed response replaces it **in the
same commit**. Evidence: `docs/validation/PHASE4C_I1_MONOPOLE_IMPLEMENTATION.md`.

| Contract clause | As implemented |
|---|---|
| P1 canonical variable | `p₀*` integrated; `δp₀`, `ξ₀` derived |
| P2 governed equations | `RotationSolver::ODE_HartleMonopole_`, each term carrying its `[term N]` marker and citing this ADR; term 7 analytically expanded so no numerical differentiation of the source appears in the RHS; **all eight background inputs interpolated at the ACTUAL driver radius through one shared bracket** (INV-13), unlike the candidate's per-node scalars; `1 − 2m/r` from `Geometry::MetricDenominator` (ADR-0004), no clamp |
| P3 fixed `ε_c` | no homogeneous admixture, no shooting, **no surface condition** |
| P4 regular-centre start | series at `r₀`; verified exact (rel `0.0`) against an independent recomputation |
| P5 EOS derivative | consumed **only** through `StarProfile::HasEosDEdP()`/`GetEosDEdP()`; absence fails the whole computation |
| P6/P7 `δM`, surface | `m̂₀(R_*) + 4πR_*²ε_*ξ̂₀(R_*) + I²/R_*³` with `I` from the verified first-order response; `R_*` documented as the EOS-floor node, never as `p = 0` |
| P9/P10 representation | `HartleMonopoleResponse` (all fields `_over_Omega2`) + `PhysicalHartleMonopole` via `At(AngularVelocity)`; one canonical geometric `Ω`; **explicit** `NStar::ComputeHartleMonopoleResponse()`, never automatic on construction |
| P11 (as modified) | the homogeneous response is **not** exposed |
| provenance | source profile identity + `Version()`; a stale response is never returned as current |

Measured: seed invariance **`7.85e-15`** (bound `1e-10`); quadratic materialization and the `δM̂`
identity exact (`0.0`); zero solves on ordinary construction, one per explicit request, none per
materialization; first-order arithmetic byte-identical and `I` bit-identical; the seven durable
artifacts unchanged. Suites **27/27** and **14/14**. Detectors D1–D4 fired and were reverted
byte-identically.

> **§3.1 status: CORRECTION EXECUTED — INDEPENDENT VERIFICATION PENDING — NO BASELINE YET.**
> No number produced by this implementation is validated physics; that requires Phase 4D
> (§7 items 3–11, detectors M1–M9), after which the first monopole baseline may be created.

---

## 14. Validation record — Phase 4D (2026-09-03)

Evidence: `docs/validation/PHASE4D_MONOPOLE_VALIDATION.md`. Production diff **NONE**; `I` and the seven
durable artifacts byte-identical; detectors M1–M9 applied, measured and reverted byte-identically.

| §7 item | Outcome |
|---|---|
| 1, 2 | seed invariance and materialization contract re-run (`hartle_monopole_contract`); `±Ω` bitwise, `Q(2Ω) = 4Q(Ω)` exact |
| 3 | centre series `≤ 5.0e-10` over the first ten nodes of a fine-centre fixture (bound `1e-8`) — **met** |
| 4 | independent `(m₀, h₀)` solver, (97)+(98)+(90): analytic `9.7e-9` (bound `1e-7`), exact `h²` scaling against the continuum solver; DS(CMF)-1 four stars `≤ 3.8e-7` isolated, `≤ 3.7e-5` fully independent (bound `1e-4`) — **met** |
| 5 | matching arithmetic exact; near-vacuum identity on the outermost nodes `≤ 1.7e-8` after the matter-source correction (raw spread `1e-6`–`4e-6` over 3–10 nodes); no profile node lies beyond `R_*` — **met as testable** |
| 6 | continuum `6.1e-15` (bound `1e-9`) — **met**; tabulated form is an exact-`h²` input floor (`1.9e-9` at N = 16001, `4.8e-10` at 32001), recorded |
| 7 | linear-in-`M/R` intercepts `1.1e-6` (`δM̂/R³`) and `4.9e-7` (`3Mξ̂/R⁴`), monotone over `M/R = 0.15…0.001` (bound `5e-3`); omitting the shell gives `4e-4` instead of 1 — **met** |
| 8 | Chandrasekhar & Miller 1974 Table I: 19/19 configurations, `I/MR²`, `ϖ₁`, `ξ₀` (incl. its sign change) and `δM/M` in C&M's shell-excluded convention to `≤ 7.3e-4`; Hartle & Thorne 1968 Tables 3/5 on the printed HW EOS: 8/8, `R`, `M`, `R_g/R`, `ω_c/Ω`, `ω_s/Ω`, `δR/R ≤ 4.9e-3`, `δM/M ≤ 1.1e-2` (bound `2e-2`) — **met** |
| 9 | radial convergence 5000/10000/20000 at fixed `ε_c`: `δM̂` order `2.35`, `ξ̂(R_*)` order `4.9`, Richardson residual `7.7e-4` on `δM̂` at 20000 — **reported** |
| 10 | no independent derivative source (`c_s²`: CONDITIONAL CHECK UNAVAILABLE); the retired profile FD moves `δM̂` by `5.0e-2` — **not independently testable** |
| 11 | analytic star: homogeneous `δM̂` vs exact `dM/dp_c` `3.0e-9` — **met**; DS(CMF)-1: `1.17e-3`/`1.02e-3`/`1.04e-3` at res 10000/20000/40000 vs `1e-3`, in both the `p_c` and `ε_c` forms — **NOT MET**, resolution-independent; the nodal `dε/dp` column integrated over the crust misses 17 % of the crust's own `Δε` at every resolution, i.e. density steps of the crust EOS that no sampled derivative represents (Hartle's `dE/dP` delta functions at internal discontinuities, which this ADR did not adjudicate). Summing the same source against the profile's own `ε` steps reproduces the TOV-sequence derivative to `7e-5`, and applied to the SOURCED solution it quantifies the omitted shells at **≈ `4.6 %` of `δM̂`** on DS(CMF)-1 — a substantive physical discrepancy of this contract, reported and not repaired |
| 12 | `I`, both Hartle CTests, seven goldens — **bitwise** |

Two 4D-entry re-scopings are recorded verbatim in the harness and the record: item 7 was first
asserted at the weakest field (`4.5e-3`, `5.3e-3` at `M/R = 0.002`) before being taken as the
*intercept* this table specifies; item 6's tabulated form was first asserted at N = 16001 (`1.9e-9`)
before being taken at the continuum level. Neither bound was widened.

> **Status: `HARTLE MONOPOLE VALIDATION FAILED`** — precisely: the implementation of this contract
> is independently verified, but the contract itself (P2 term 1 with the smooth P5 derivative)
> omits Hartle's internal delta-function shells on a tabulated crust with density steps, worth
> ≈ `4.6 %` of `δM̂` on DS(CMF)-1. **§3.1: CORRECTION EXECUTED — IMPLEMENTATION INDEPENDENTLY
> VERIFIED — PHYSICAL VALIDATION FAILED ON STEPPED CRUSTS — NO BASELINE YET.** Next: amend this ADR
> for internal density discontinuities, implement, re-run 4D Experiment J, then the first baseline.

**Post-validation note (2026-09-03, Phase 4D-RG).** The amendment is drafted as `docs/adr/ADR-0008-measure-complete-eos-energy-density-source.md`
(**PROPOSED**): the EOS energy-density source of (97) is taken as the measure `−4πr²ξ̂₀ dε` of Hartle's
eq. (93), with the nodal `dε/dp` column demoted from source representation to diagnostics and the
surface shell the terminal atom of the same measure. It amends/supersedes **P2, P5 and P6 only where EOS
energy-density variation is represented by a measure rather than a pointwise radial derivative**; P1,
P3, P4, P7–P14 and the Decision below are unaffected. Evidence: `docs/validation/PHASE4D_R_EOS_MEASURE_DERIVATION.md`.

---

## Decision

**Adjudicated by the project owner on 2026-09-02.** The contract of §4 (**P1–P14**) is
**ACCEPTED**, together with the `GOVERNANCE.md` §3.1 record of §9, with **one modification**
(P11, recorded below). The owner's eight adjudications map onto this ADR's thirteen decision
questions (§2) as follows.

| Owner decision (as adjudicated) | ADR question(s) | Outcome |
|---|---|---|
| **Canonical integrated pressure variable = Hartle `p₀*`** | Q1 | §5-Q1 **A**. `δp₀`, `ξ₀`, `δε₀` are derived fields (P1) |
| **Governed family = fixed central energy density `ε_c`** | Q2 | §5-Q2 **A**. `m̂₀(0) = p̂₀*(0) = 0`, no homogeneous admixture, **no surface condition of any kind** (P3) |
| **Second-order scientific products = seed-free coefficients per `Ω_geom²`, materialized only at an explicit `AngularVelocity`** | Q9, Q10 | §5-Q9 **A** (P9, P10). The ADR-0006 no-implicit-spin clarification extends to O(Ω²): no `NStar` carries a physical second-order perturbation without an explicit spin |
| **`dε/dp` authority = derivative of the SAME EOS `ε(p)` interpolant that constructed the star, owned by the EOS/TOV layer** | Q5 | §5-Q5 **A** (P5). Profile finite differences are not an authority; no second interpolant; no `1.0` fallback. Implemented and validated as increment **4C-I0** |
| **EOS-floor surface = `SURFACE ADEQUATE AS-IS`** for this monopole / Phase-5 scope, with explicit `R_*` semantics, the shell and boundary terms, and the documented surface-displacement systematic. **`R_*` must never be labelled as the exact `p = 0` surface** | Q6, Q7, Q8 | §5-Q6 **A** (P6, P7, P8). No change to INV-06 or to the TOV surface convention |
| **Validated `l = 0` O(Ω²) structural response is the Phase-4 second-order deliverable that unlocks Phase 5.** `l = 2` shape / quadrupole physics is a **separate future rotation extension** and does **not** block Phase-4 completion, Phase 5, or the BNV thermal program | Q11 | §5-Q11 **A** (P12). **No claim is made that `l = 2` physics is validated** — it is out of scope, not verified |
| **Current public second-order candidate API = atomic replacement in 4C-I1** | Q13 | §5-Q13 **A**. It is **not** deleted in 4C-I0; until 4C-I1 it remains publicly callable, zero-caller, and marked an unverified candidate |
| **Homogeneous sequence-derivative response = NOT a public production API in Phase 4C** | — (modifies P11) | The mathematics remains valid and **may be computed internally or test-side in 4D** to validate the sequence-derivative identity (§7 item 11). Public ownership of `B_i` / `∂N_i/∂ε_c` is **deferred until Phase 5** adjudicates that architecture |
| **`GOVERNANCE.md` §3.1 record explicitly accepted** | Q12 | §9 items 1–7 are in force. §3.1 is **AUTHORIZED**; the correction itself is **not yet executed** — 4C-I0 adds no O(Ω²) physics and creates no monopole baseline |

Questions **Q3** (canonical `l = 0` equations, P2) and **Q4** (regular-centre numerical start, P4)
are accepted as proposed in §4; they were not varied.

**Implementation sequencing (owner, binding).** 4C-I0 establishes the `dε/dp` authority alone and
must leave every existing scientific output byte-identical; 4C-I1 then performs the atomic
replacement. The second-order solver must not execute without authoritative `dε/dp`, but ordinary
first-order/TOV construction must **not** begin to require it.

**Scope of acceptance.** This settles the O(Ω²) *monopole* contract. It ratifies no number: no
O(Ω²) quantity in this repository is validated physics until 4C-I1 conformance and 4D independent
validation are complete. It decides nothing about the `l = 2` sector, the chemical-imbalance
convention (INV-11), MixedStar rotation, or the solar-mass / `G` authority.

## Provenance

Drafted 2026-09-02 by an AI agent (Claude, Phase 4C-G) from a page-by-page reading of the
journal scans of Hartle (1967) and Hartle & Thorne (1968) — not from the candidate, its
comments, or memory — with every derivation, dimensional audit, surface bound and repository
citation recorded in `docs/validation/PHASE4C_HARTLE2_DERIVATION.md`. The agent presents
alternatives and recommendations; **the project owner decides.** Two misprints in the primary
source (H67 (115) sign, (117) factor) are recorded in the evidence record §2.3; neither enters
this contract.

**Accepted 2026-09-02 by the project owner** (adjudication recorded verbatim in the Decision
section above), at `bcef5b57e3267fe26d1269efaa472d314818d367`. The owner varied one item of the
drafted contract — P11's optional public exposure of the homogeneous response — and accepted the
remainder, including the §3.1 record, as drafted.

**Post-revalidation pointer (2026-09-03, Phase 4D-RV).** The corrected independent revalidation of this
contract as amended by ADR-0008, on the ADR-0009 surface event, is recorded in
`docs/validation/PHASE4D_CORRECTED_MONOPOLE_REVALIDATION.md`: status `HARTLE O(OMEGA^2) MONOPOLE RESPONSE CHARACTERIZED — INDEPENDENT VALIDATION INCOMPLETE` — every physics line of §7 passes on the migrated backgrounds
(analytic `9.7e-9`, DS(CMF)-1 `≤ 5.3e-5` against an independent EOS-knot measure route, C&M 1974 and
HT68 1968 reproduced, M1–M10 fire), while ADR-0008 Validation D's monotonicity clause is not met at the
node-placement floor. §9 item 7's baseline remains deferred pending the owner's adjudication of that clause.


## Post-closeout pointer — Phase 4E (2026-09-04)

The accepted Decision is unchanged. Correction → independent revalidation → owner
convergence adjudication → VERIFIED → first trusted baseline (`2e2f016`) → **Phase-4
closeout / Phase-5 structural-interface ratification**. GOVERNANCE §3.1 condition 7
was discharged by the separate baseline commit; this closeout adds no scientific layer.

The existing `HartleFirstOrderResponse` and `HartleMonopoleResponse` are ratified as
normalized structural inputs, with their existing units and provenance/lifetime rule.
INV-08 is CLOSED / VERIFIED / REGRESSION-PROTECTED for ordinary-NStar fixed-ε_c l=0
O(Ω²) on ADR-0009 backgrounds. This is not a fixed-baryon-number response; INV-09
remains unresolved. l=2 is unimplemented/unvalidated and nonblocking; MixedStar and
high-spin/O(Ω⁴) accuracy are not covered. **PHASE 5 NOT YET BEGUN.**
Current contract, unchanged eight-artifact authority, final suites and integration gate:
`docs/validation/PHASE4_CLOSEOUT.md:33`, `docs/validation/PHASE4_CLOSEOUT.md:64`, `docs/validation/PHASE4_CLOSEOUT.md:158`,
`docs/validation/PHASE4_CLOSEOUT.md:183`, `docs/validation/PHASE4_CLOSEOUT.md:201`, `docs/validation/PHASE4_CLOSEOUT.md:219`.
Historical failed/CHARACTERIZED/no-baseline statements above retain their dated meaning.
