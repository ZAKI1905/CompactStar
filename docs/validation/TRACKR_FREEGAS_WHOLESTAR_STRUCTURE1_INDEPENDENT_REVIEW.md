# Track-R Structure-1 — independent scientific review

**Reviewed SHA:** `3b0499e3cb5fb39db63d77fb136b4098b8084b4f`
**Disposition:** **B — TRACK-R STRUCTURE-1 ACCEPTABLE WITH QUALIFIED CLAIM.**
Static free-gas structure and the common-central-state Table-1 numbers are
independently verified; the source-qualified `M_max` meaning is **not**
reproduced, and under the supremum reading it is independently falsified as a
description of the printed row. No blocking physics, numerical, surface,
source-comparison or oracle failure was found. Review only; no production or
test file was modified.

This review was performed by one agent against the primary PDF, the source
code, and independently coded numerics. Structure-1's own documentation was
read but was **not** used as source authority.

---

## 1. Reviewed SHA

| Item | Value |
|---|---|
| Reviewed HEAD | `3b0499e3cb5fb39db63d77fb136b4098b8084b4f` |
| Branch | `physics/trackr-freegas-wholestar-structure1` |
| Worktree | `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-trackr-freegas-structure1` |
| Commit subject | `feat: reproduce Track-R free-gas static structure` |

## 2. Authentication

All checks performed before any scientific conclusion.

| Check | Result |
|---|---|
| Working tree clean | **PASS** — `git status --porcelain` empty, before and after review |
| HEAD exact | **PASS** — `3b0499e3cb5fb39db63d77fb136b4098b8084b4f` |
| local == upstream | **PASS** — `@{u}` = same SHA |
| live remote == local | **PASS** — `git ls-remote origin` returns `3b0499e3…` for the branch |
| canonical master unchanged | **PASS** — local, `origin/master` and live remote all `f4ae22d971e25bdd74530aec184f3fe0c3440b95` |
| master is exact merge-base | **PASS** — `git merge-base HEAD master` = `f4ae22d9…` |
| exactly one implementation commit ahead | **PASS** — `git rev-list --left-right --count master...HEAD` = `0  1` |

Governing authority read and applied: `AGENTS.md`, `GOVERNANCE.md` (§1 hierarchy,
§2 classes, §3 fail-closed, §5 AI-candidate status), ADR-0005, ADR-0009,
ADR-0010, `PHASE5A5_TRACKR_LOCAL_RATIFICATION.md`,
`TRACKR_FREEGAS_WHOLESTAR_INTERFACE_AUDIT.md`, `docs/SCIENTIFIC_INVARIANTS.md`,
`docs/MODERNIZATION_ROADMAP.md`, `docs/architecture/CURRENT_ARCHITECTURE.md`.

## 3. Direct source ledger

Primary PDF read directly:
`/Users/keeper/Documents/CompactStar/literature/rotochemical/2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf`,
SHA-256 `f184d7d1d7030b61a021eb5c7ac14b1f1b30c7ea69e9d53473d153cfb069ea88`
(hash independently recomputed; matches the Structure-1 record).

**Table 1, p. 30 — free-gas row, transcribed from the PDF:**

| EOS | M_max (M☉) | ρ_c (10¹⁵ g cm⁻³) | R (km) | R_∞ (km) | P_K (ms) |
|---|---|---|---|---|---|
| Fermi gas | 0.62 ᶜ | 1.10 | 12.77 | 13.80 | 0.98 ᵈ |

Footnote c: *"Corresponds to maximum mass before appearance of Σ− hyperons"*.
Footnote d: *"Calculated with empirical formula (see text), except the last
value, which was adopted from Haensel et al. (1995)"* — the last value **is**
the Fermi-gas 0.98 ms. §3.1 p. 9 confirms it is *"the exact value calculated for
a pure neutron gas"*. **P_K is therefore correctly excluded** from the static
reproduction: it is neither this EOS nor a static quantity.

**§3.1, pp. 8–9 — the load-bearing sentences:**

- *"a non-interacting Fermi gas (see, e.g., Shapiro & Teukolsky 1983) for the
  **whole star**"* — no crust splice. Confirms the model applied here.
- *"In Table 1 we list, for each EOS, the mass M, **central energy density** ρ_c,
  coordinate radius R, and effective radius as seen from infinity R_∞ of the
  maximum-mass non-rotating configuration."*
- *"The values for the Fermi gas EOS listed in Table 1 correspond to a central
  density **slightly below** the appearance of Σ− hyperons."*

ρ_c is explicitly **central energy density**, i.e. ε/c² in g cm⁻³, not a baryon
rest-mass density. Structure-1's interpretation is correct.

**Metric / radius definitions, pp. 3 and 6:**
`g_tt = −e^{2Φ}`; `R` is the stellar coordinate radius; `Φ_s = Φ(R)`;
**`R_∞ = R e^{−Φ_s}`**; `g_rr = e^{2Λ} = (1 − 2M/r)^{−1}`. Eq. (6),
`L_γ^∞ = 4πσR²T_s⁴e^{2Φ_s} = 4πσR_∞²(T_s^∞)⁴`, is consistent. With an exterior
Schwarzschild match this is `R_∞ = R/√(1 − 2GM/Rc²)` — the formula used.

**p. 16, §4.3 (hyperons):** hyperons were **not** included in the calculations;
*"The first particles to appear after muons are probably the Σ− and Λ0
hyperons"*; eq. (71) `η_nnΣp = 2μ_n − μ_Σ − μ_p`. I verified independently that
in this free gas Σ⁻ does appear before Λ⁰ (Σ⁻ needs μ_n ≳ 1073.4 MeV given
μ_e ≈ 124 MeV; Λ⁰ needs μ_n ≥ m_Λ = 1115.683 MeV), so the code's choice of the
Σ⁻ ceiling is right. The code's threshold residual `m_Σ + μ_p − 2μ_n` is
algebraically identical to `μ_Σ = μ_n + μ_e` under npeμ beta equilibrium.

Supporting sources consulted: R1995 and R2006. Neither introduces a different
static free-gas EOS; R2006's electrostatic correction governs response scope,
not this static barotrope.

**No source-authority defect was found in the Structure-1 source ledger.**

## 4. Implementation-diff scope

`git diff master...HEAD` touches 16 files, +3126/−3.

| Group | Files |
|---|---|
| Production (additive only) | `CompactStar/EOS/LocalThermodynamics.{hpp,cpp}` (+38), `CompactStar/EOS/TrackRFreeGasThermodynamics.{hpp,cpp}` (+96) |
| Tests / tools | 2 executables, 3 test-only headers, `collect.py`, `tests/CMakeLists.txt` (+10) |
| Documentation | Structure-1 record, results JSON, roadmap, invariants, architecture |

**Scope integrity — independently verified:**

- `TOVSolver.cpp`, `TOVSolver.hpp`, `NStar.cpp`, `NStar.hpp`, `StarProfile.cpp`,
  `StarProfile.hpp`, `RotationSolver.cpp`, `ThermalSolver.cpp`:
  `git diff --stat master...HEAD` is **empty** — byte-identical.
- No file under `Physics/Evolution`, rotation, thermal, or BNV is touched.
- `tests/baselines/` unchanged; **8** baselines, no ninth.
- `tests/CMakeLists.txt` diff contains **no deletions** — no test removed or
  reclassified.
- The production diff contains **no `try` and no `catch`** — no exception can be
  caught and replaced by a neighbouring branch or state.

## 5. Independent barotrope review (value-only path)

`TrackRFreeGasThermodynamicProvider::BarotropeAt` was read line by line and its
algebra re-derived from scratch.

**Hessian independence.** `BarotropeAt` never calls `EquilibriumAt`,
`EvaluatePe`, `EvaluateNpe` or `SolveNpeEquilibrium`. The npe and npe-µ branches
use only rest masses, `ħc`, `NumberDensityForChemicalPotentialFm3` and the new
`Values`. The `else` branch (pe, exact neutron onset, exact muon onset) calls
`EquilibriumStateAt`, which for those three domains constructs the state
directly or via `EvaluateNeutronThreshold`/`EvaluateMuonThreshold`; no
`ChemicalHessian` is built and no N-1 guard is on that path. *(Minor
technicality: those two threshold helpers call `Evaluate`, which always computes
`dμ/dn` internally; nothing consumes it and nothing can fail on it.)*

**N-1 / H semantics preserved.** The `2^30`-ULP guard in `EvaluateNpe`
(`TrackRFreeGasThermodynamics.cpp:136-142`) is untouched by the diff, and all
response methods are unchanged.

**Fail-closed above Σ⁻.** `BarotropeAt` calls `EquilibriumDomainAt` first, which
throws for `n_B ≥ σ_onset`, `n_B < 0`, or non-finite. Verified by execution:
`BarotropeAt(σ_onset)` and `BarotropeAt(0.7)` both throw.

**npe δ-parameterisation — re-derived independently.**
With `h = μ_n`, `n_p = n_e` ⟹ `k_p = k_e`, so `μ_e = (h² − m_p² + m_e²)/(2h)`.
Writing `h = m_n + δ` and `f(h) = h/2 + (m_e²−m_p²)/(2h)`:

```
Δμ_e = f(m_n+δ) − f(m_n) = (δ/2)[1 + (m_p²−m_e²)/(m_n(m_n+δ))]
```

which is exactly the code's `dt`. Then `k_e² = k_0² + dt(2t_0+dt)`, so
`n_e/n_e0 = (1 + dt(2t_0+dt)/k_0²)^{3/2}` and the increment is
`n_e0·expm1(1.5·log1p(…))` — exactly the code. `n_e0` is `neutron_onset_n_B_fm3_`,
correct because `n_n = 0` at onset makes `n_B = n_p = n_e`. The neutron branch is
`n_n = [δ(2m_n+δ)]^{3/2}/(3π²ħ³c³)`, exact. The bisection target is
`n_n + Δn_e = n_B − n_onset`, which removes the catastrophic `n_B − n_e`
cancellation that defeats the response path. **The algebra is correct**, both
summands are monotone in δ, and the bracket `[0, √(m_n²+k²) − m_n]` is valid
(written cancellation-free as `k²/(hypot(m_n,k)+m_n)`).

**npe-µ common-potential inversion — re-derived.** With `t = μ_e = μ_µ`,
`n_p = n_e(t)+n_µ(t)`, `μ_n = μ_p(n_p)+t`, `n_n = n(μ_n)`. `n_n+n_p` is monotone
increasing in `t`; the bracket `[m_µ, μ_n(n_B) − m_p]` is valid because the true
`μ_n` uses `n_n < n_B`. Correct. The final `n_µ ≤ 0` check fails closed rather
than manufacturing a positive muon density.

**Numerical verification against 60-digit mpmath** (my own EOS, no CompactStar
code), over 29 densities spanning `10⁻¹⁵ → 0.6173 fm⁻³` including both
thresholds and their `nextafter` neighbours:

| Quantity | Worst relative error |
|---|---|
| ε | **3.96e-13** |
| P | **4.06e-13** |
| n_n (bulk, n_B ≥ 1e-6) | ≤ 2.6e-16 |
| n_e, n_µ (bulk) | ≤ 3.4e-15 |

Vacuum returns exactly zero energy, pressure and densities; pe gives
`n_p = n_e = n_B`, `n_n = n_µ = 0`; exact thresholds retain their boundary
identity; **no species floor and no branch substitution anywhere.**

**Independent threshold reproduction** (mpmath, from the dumped binary
constants, solving the physics myself):

| Threshold | Independent value | Stored value | Rel. diff |
|---|---|---|---|
| neutron onset | 7.356728903732990e-9 | 7.3567289037328318e-9 | **2.15e-14** |
| muon onset | 0.45698480541241896 | 0.45698480541241987 | **1.98e-15** |
| Σ⁻ ceiling | 0.61735520796652978 | 0.61735520796653498 | **8.43e-15** |

**Verdict: the value-only barotrope is correct.**

## 6. Stable pressure implementation review

Derived independently from the T=0 spin-1/2 (g=2) phase space:

```
n = k³/(3π²b³),   b = ħc
ε = (π²b³)⁻¹ ∫₀ᵏ q²√(m²+q²) dq
P = (3π²b³)⁻¹ ∫₀ᵏ q⁴/√(m²+q²) dq
```

matching the task statement and the record. The closed form follows from
`∫₀ˣ t⁴(1+t²)^{−1/2}dt = ⅛[x√(1+x²)(2x²−3) + 3 asinh x]`, which I verified by
differentiation (`d/dx = x⁴/√(1+x²)`), giving

```
P = m⁴[x√(1+x²)(2x²−3) + 3 asinh x] / (24π²b³)
```

— **exactly the production expression**, with correct degeneracy, dimensions and
overall factor.

**Small-x series.** `(1+t²)^{−1/2} = Σ C_j t^{2j}`, `C_j = Π_{i≤j}[−(2i−1)/(2i)]`,
so `∫ = x⁵[1/5 + Σ_{j≥1} C_j x^{2j}/(5+2j)]` and
`P = m⁴x⁵·(…)/(3π²b³) = n k²/m · (…)`. The loop's `coefficient`, `power` and
`sum` reproduce `C_j`, `x^{2j}` and the sum term-for-term; `pressure = n·k²/m·sum`
is algebraically exact; the leading term is `n k²/(5m)` as claimed. Truncation at
`j = 12` leaves `|C₁₃|x²⁶/31 ≈ 5e-29` at `x = 0.1` — **utterly negligible**;
vacuum behaviour is correct and no rest-energy subtraction occurs.

**Independent high-precision quadrature check** (mpmath, 42 species/x samples,
`x` from 1e-4 to 400 across the switch):

| Region | Worst relative P error |
|---|---|
| `x < 0.1` (series) | ≤ 2.1e-16 |
| `x ≳ 0.3` (closed form) | ≤ 3.8e-16 |
| **`x ∈ [0.1, 0.11]` (closed form, at the switch)** | **1.96e-12** |

This is expected cancellation: at `x = 0.1` the bracket is `≈1.6e-5` while its
two terms are each `≈0.3`, so ~4 digits are lost. The pre-existing energy
bracket shows the same at its own `x = 1e-2` switch (2.6e-14).

**Materiality:** the neutron crosses `x = 0.1` at `n_n ≈ 3.65e-3 fm⁻³`, where
`P_n ≈ 6.9e-3 MeV fm⁻³` against a central `P_c = 30.68 MeV fm⁻³` — a ≲2e-16
relative perturbation of the star's pressure scale. Muon and electron crossings
are smaller still. **No material local-physics defect is introduced.** The S1
tolerance (1e-8) has four orders of margin.

## 7. P-7 adjudication

Executed directly against the built library, comparing `BarotropeAt` (value) with
`EquilibriumAt` (response) at relative offsets above neutron onset:

| rel. offset | `BarotropeAt` | n_n returned | `EquilibriumAt` |
|---|---|---|---|
| 1e-6 | VALUES OK | 6.03e-15 | OK (above the 2^30-ULP guard) |
| **1e-7** | **VALUES OK** | 4.87e-16 | **`EquilibriumResolutionError`** |
| **1e-9** | **VALUES OK** | 1.68e-18 | **`EquilibriumResolutionError`** |
| 1e-10 … 1e-15 | VALUES OK | 6.8e-20 … 2.96e-27 | `EquilibriumResolutionError` |
| `nextafter(onset,∞)` | VALUES OK, domain = Npe | 9.37e-29 | (refusal region) |

Pushed to the limit of double representability: **value resolution never fails**
anywhere in the representable npe interval. ε and P at
`nextafter(onset,0)`, `onset` and `nextafter(onset,∞)` are bit-identical —
correct, since the newborn neutron contributes ~8.8e-26 MeV fm⁻³ against an ulp
of ~8.7e-22. Domains are correctly distinct (1 → 2 → 3) and no exception is
swallowed.

**Classification:**

- **VALUE availability: complete** over the whole source domain, including
  arbitrarily close to onset.
- **CHEMICAL RESPONSE availability: unchanged and still refused** below the
  `2^30`-ULP criterion.

**P-7 is genuinely resolved for structure, with the chemical-H guards intact.**

Accuracy caveat: near onset the returned `n_n` inherits the stored threshold's
own 2.15e-14 relative error, amplified as `≈1.5·ε_thr/offset`. Measured against
fully exact physics: 2.29e-8 (offset 1e-6), 2.42e-7 (1e-7), 2.89e-5 (1e-9),
3.08e-4 (1e-10). These are **larger than §17's reported 1.40e-9 / 1.67e-7 /
1.78e-6**, which evidently use a reference anchored on the stored double
threshold rather than exact physics; §17 does not state the convention. The
mechanism §17 names is correct. Impact on ε and P is ≤3e-14 relative —
immaterial (Finding 7).

## 8. Generated-table review

Regenerated the tables by running the executable into a fresh directory, then
audited the 10375-row accepted table independently.

| Check | Result |
|---|---|
| Exact header | `eps(g/cm^3) ⇥ pre(dyn/cm^2) ⇥ rho(1/fm^3) ⇥ 10 ⇥ 11 ⇥ 0 ⇥ 1` — **as documented** |
| Species ids | `10`=n, `11`=p, `0`=e, `1`=µ — matches `StarContext.cpp:494-497` |
| Columns are Y_i, not n_i | **confirmed** (`n_i/n_B` written) |
| Strictly increasing P | **PASS** |
| Strictly increasing ρ | **PASS** |
| Strictly increasing n_B | **PASS** |
| P > 0 on every row (no vacuum row) | **PASS**, min P = 1.0850284294719746e15 dyn cm⁻² |
| Threshold knots present exactly | **PASS** (both) |
| Strictly sub-Σ | **PASS**, margin 1.3552e-3 fm⁻³ |
| Baryon identity `|Y_n+Y_p−1|` | ≤ **2.0e-15** |
| Charge identity `|Y_p−Y_e−Y_µ|` | ≤ **2.4e-18** |
| Causality `max dP/dρc²` | **0.07583** ≤ 1/3 |
| Round-trip precision | 17 significant digits |
| Upper ρ / clamp | 1.1211542516959692e15 / 1.1200330974442732e15; **clamp margin above 1.105e15 = 1.5033e13** |

Units re-derived independently: `ρ = ε·1.602176634e33/c²` gives
`ρ/ε = 1.7826619216278975e12`; `P_cgs = P·1.602176634e33`. Correct.

**Not self-referential:** I regenerated every table from the executable, and
separately recomputed the EOS at 29 off-grid points in 60-digit mpmath (§5) and
the full star from an analytic EOS with no table at all (§10). The committed
numbers reproduce exactly.

## 9. Derivative / interpolation review (S2)

Reproduced the nested off-grid ladder exactly:

| grid | max rel ε(P) err | max rel n_B(P) err | max abs fraction err | **max rel interpolant slope err** |
|---:|---:|---:|---:|---:|
| 1024 | 4.1981e-5 | 4.1981e-5 | 6.5629e-6 | 0.12808 |
| 2048 | 1.0256e-5 | 1.0256e-5 | 6.2762e-7 | 0.041902 |
| 4096 | 2.5336e-6 | 2.5336e-6 | 2.3705e-7 | 0.013329 |
| 8192 | 6.2958e-7 | 6.2958e-7 | 2.8784e-8 | **0.0040951 (0.410 %)** |

**Production trace (decisive).** I read the canonical RHS directly:

- `TOVSolver::ODE` (`TOVSolver.cpp`) computes
  `f[0] = 4πr²ε`, `f[1] = −(G/r²)(ε+P/c²)(m+4πr³P/c²)/metric`
  — **`GetEDens(p)` only; `dε/dP` appears nowhere.**
- `PressureCoordinateODE` (the surface locator) **reuses `ODE`** — no `dε/dP`.
- `GetNuDer` (the lapse derivative) uses **`GetEDens` only**.
- `GetEDensDeriv` feeds only the published `dEdP` profile column, whose sole
  consumer is `RotationSolver` (ADR-0007 P5).

**Contribution of the 0.410 % to M, R, the surface, and profile construction:
exactly zero** — the quantity is not in the static computation path.

I also identified *why* the error is large and localised: `dε/dP` has a
square-root cusp at neutron onset. Measured one-sided approach at relative
offsets 1e-2 → 1e-10: **402035 → 88504 → 20162 → 5668 → 2975.6** against the
common limit **2581.06**. This is structural physics that no polynomial
interpolant resolves at finite table resolution — not a coding defect.

**Disposition: harmless / nonblocking for Structure-1 static reproduction, and a
must-quantify item before any consumer of `dEdP`.** I did not demand a tighter
derivative.

## 10. Independent TOV oracle: derivation and reproduction (S12)

**Independence audit of `tests/eos/structure1/tov_oracle.hpp`:**

| Shares with production? | Verdict |
|---|---|
| Production TOV RHS | **No** — enthalpy form in `(u=r², v=m/r³)` vs radial `(m,P)` |
| Production ε(P) Steffen spline | **No** — direct 96-node phase-space quadrature |
| `SingleStarSolveToTOVPoints` | **No** — never called |
| Production surface locator | **No** — integrates to matched H or to vacuum H=0 |
| Production center initialisation | **No** — regular-center series in δ vs a 1-cm uniform start |
| Integrator | **No** — GSL RKF45 vs canonical RK8PD |
| Physical constants | **Yes** — acceptable per the review mandate |

**Equations re-derived from scratch.** With `H = ln(h/h₀)`, `h₀ = m_p+m_e`
(correct: `h = (ε+P)/n_B = μ_n`, or `μ_p+μ_e` in pe, and `h → m_p+m_e` as P→0),
`x = H_c − H`, `u = r²`, `v = m/r³`, and `A = (1−2uv)/(v+4πP)`:

`uv = m/r`, so `A = r²(r−2m)/(m+4πr³P)`. From `dr/dH = −r(r−2m)/(m+4πr³P)`:

```
du/dx = 2r²(r−2m)/(m+4πr³P) = 2A                    ✔ matches code
dv/dx = A[4πε − 3v]/u                               ✔ matches code
```

**Center expansion re-derived.** `P ≈ P_c − (2π/3)(ε_c+P_c)(ε_c+3P_c)r²` and
`P_c − P = (ε_c+P_c)δ` give `u = 3δ/(2π(ε_c+3P_c))` ✔; integrating
`m = ∫4πr²ε` to O(r⁵) gives `v = 4πε_c/3 − (4π/5)(ε_c+P_c)(dε/dP)_c δ` ✔. Both
exactly as coded.

**Units and mass convention.** `to_geo = 1.602176634e33·G/c⁴·1e10` converts
MeV fm⁻³ → km⁻² (dimensions check: `[erg cm⁻³]·[s² g⁻¹ cm⁻¹] = cm⁻²`, ×1e10 →
km⁻²) ✔. `solar_length = G·M☉/c²/1e5` km ✔; `M = v r³/solar_length` ✔.

**Critically, I verified the mass convention is consistent**, because a mismatch
would be invisible in a passing test: `SingleStarSolveToTOVPoints` emits
`y[0]/GSL_CONST_CGSM_SOLAR_MASS` (`TOVSolver.cpp:2591`) and integrates with
GSL's `G`. The oracle divides its geometric mass by `G·M☉_GSL/c²`. Both are
therefore `m_grams/M☉_GSL`. Had the oracle used `SUN_M_KM` instead, the S12 mass
comparison would have carried a **3.85e-5 M☉ offset** — 380× the gate. It does
not. The `1e-7 M☉` gate is convention-consistent.

**Independent reproduction — a third implementation.** I wrote my own TOV
integrator: Python + SciPy `DOP853`, enthalpy variable but state `(r, m)`
(different from the oracle's `(u,v)`), my own analytic EOS, my own center
expansion, integrating **to exactly H = 0** so the P=0 radius is obtained
directly with no cutoff or tail model.

| ladder (δ/H_c, rtol) | M [M☉] | R₀ [km] |
|---|---|---|
| 1e-9, 1e-11 | 0.623635569243 | 12.7681549015 |
| 1e-10, 1e-12 | 0.623635569159 | 12.7681549010 |
| 1e-11, 1e-13 | 0.623635569158 | 12.7681549010 |
| 1e-12, 1e-13 | 0.623635569157 | 12.7681549010 |

**Three-way agreement at the midpoint star:**

| Source | M [M☉] | R₀ [km] |
|---|---|---|
| Production canonical (8192 table, 40000 partition) | 0.623635569063 | 12.7681549034 (upper bound) |
| Repo test oracle (finest ladder) | 0.623635569137 | 12.7681549010 |
| **This review, independent** | **0.623635569157** | **12.7681549010** |

Spread: **≤9.4e-11 M☉** and **≤2.4e-9 km**.

**Challenge to the claimed δR < 6e-10 km / δM < 3e-11 M☉ agreement:** it is
genuine, not a hidden shared assumption. My fresh run gives, at the finest
partition, `|ΔM| = 1.18e-11 M☉` and `|ΔR_cut| = 7.54e-11 km`; the recorded
`2.19e-11 / 5.42e-10` are the **max over the five partitions** (the field name
`finest_midpoint_oracle_disagreement` is imprecise, but the recorded values are
the more conservative ones — nothing is overstated). Across all 159 emitted rows
the maxima are `3.41e-8 M☉` and `1.60e-6 km`, both from the coarse 1024/2048
tables, which are convergence evidence rather than the accepted result.

**The S12 oracle is genuinely independent and load-bearing.**

## 11. Surface / zero-pressure tail review (S7)

**Canonical cutoff unchanged.** `TOVSolver::PressureCutoff()` is
`max(1e-15·GetInitPress(), eos_tab.pre[0])` — and `TOVSolver.cpp` is
byte-identical to master. ADR-0009 coefficients untouched.

**Tail bounds re-derived.** From `−dH/dr = (m+4πr³P)/[r(r−2m)]`, positive mass
growth and monotone P, ε in the tail, the four inequalities `r_upper`,
`m_upper`, `g_upper`, `r_lower` follow with correct signs and geometric-unit
conversions; both denominators are checked positive in the code.

**Plateau authenticated** (regenerated run):

| effective cut/P_c | R_cut [km] | [R₀ lower, R₀ upper] | width [km] |
|---:|---:|---|---:|
| 1e-11 | 12.7006956998 | [12.7677684382, 12.7681549014] | 3.865e-4 |
| 1e-12 | 12.7388701287 | [12.7680820716, 12.7681549005] | 7.283e-5 |
| 1e-13 | 12.7559892544 | [12.7681423329, 12.7681549015] | 1.257e-5 |
| 1e-14 | 12.7632197893 | [12.7681528315, 12.7681548998] | 2.068e-6 |
| **1e-15** | **12.7661747597** | **[12.7681545694, 12.7681549024]** | **3.330e-7** |
| floor 1e-16, cut 1e-15 | 12.7661747587 | [12.7681545684, 12.7681549014] | 3.330e-7 |
| floor 1e-18, cut 1e-15 | 12.7661747574 | [12.7681545671, 12.7681549001] | 3.330e-7 |

The plateau at cut/P_c = 1e-15 is real: once the table floor drops below the
ADR-0009 coefficient, the coefficient binds. `R_cut ≈ 12.7661747607`, tail
`≈ 1.980 m`, bracket width `3.330e-7 km` — **all confirmed.** `R_cut` is never
relabelled as R₀. M_cut stays at 0.623635569 ± 1e-10.

**The independent vacuum star lies inside the bracket — but at its very top.**
My direct H→0 integration gives **12.7681549010**, versus
`R₀_upper = 12.7681549024` (2.4e-9 km below it) and
`R₀_lower = 12.7681545694`. Across the whole ladder `R₀_upper` is essentially
constant and converges to the true P=0 radius: **`r_upper` is near-tight;
`r_lower` is loose.**

**Is the bracket midpoint legitimate?** It is *sound* — the true value is
provably inside and the stated half-width (1.665e-7 km) covers the residual —
but it is a **biased** estimator: the reported `R₀ = 12.7681547369` sits
1.65e-7 km below the truth, i.e. essentially the full half-width. That bias is
60× below the 1e-5 km comparison margin and ~3e4× below the bin half-width, so
it is **immaterial to every claim made**. It is nonetheless a hardening item
(Finding 3): `r_upper`, or the direct-vacuum value, is the better estimate.

**Propagated uncertainties.** At the midpoint `dR_∞/dR ≈ 1.081`, so the 1.65e-7
km radius bias maps to 1.78e-7 km in R_∞ — I confirmed this directly (my R_∞ vs
the record's differ by 1.63e-7 km in both conventions). Mass is insensitive to
the tail at this precision (M_cut varies by ≤1.5e-10 M☉ across the whole ladder).

## 12. Midpoint-star reproduction

**Central state, independently inverted in 60-digit mpmath** from the dumped
binary constants, with no CompactStar code and no table:

| Quantity | This review (independent) | Structure-1 | Rel. diff |
|---|---|---|---|
| n_B,c [fm⁻³] | 0.6049252413288510711336 | 0.6049252413288510735586 (90-digit) | **4.01e-18** |
| n_B,c vs the double used by the test | — | 0.6049252413288515 | 7.09e-16 (≈3 ulp) |
| ε_c [MeV fm⁻³] | 617.05474641848975168 | 617.05474641848974704 | **7.51e-18** |
| P_c [MeV fm⁻³] | 30.681092089959359537 | 30.681092089959357232 | **7.51e-17** |
| achieved ρ | 1.100000000000000e15 | — | — |
| composition Y_n, Y_p, Y_e, Y_µ | 0.984339, 0.015661, 0.013736, 0.0019252 | — | — |
| (dP/dε)_c | 0.075233 | — | ≤ 1/3 ✔ |

**§17's 90-digit numbers are independently confirmed.** No table interpolation or
pressure-inversion bias was found; the reported achieved-vs-requested analytic
difference of **9.2625e5 g cm⁻³** (sign: achieved slightly *below* requested,
`analytic_rho = 1.099999999073754e15`) reproduces exactly and is 1000× below the
declared 1e9 g cm⁻³ envelope. Its structural effect is `9.0e-11 M☉` and
`2.5e-9 km`.

**No target-mass inversion is used anywhere.** `SolveToProfile` (the inverse-mass
solver) is never called; the tests use `SingleStarSolveToTOVPoints` at fixed
central density only, and the midpoint was **predeclared** at the parent commit
(`TRACKR_FREEGAS_WHOLESTAR_INTERFACE_AUDIT.md:601`).

**Midpoint results (all four printed quantities):**

| Quantity | Independent | Structure-1 | Printed bin | Verdict |
|---|---|---|---|---|
| M [M☉] | 0.623635569157 | 0.623635569063 | [0.615, 0.625) | **IN** |
| R₀ [km] | 12.7681549010 | 12.7681547369 | [12.765, 12.775) | **IN** |
| R_∞,₀ profile [km] | 13.8023680743 | 13.8023679117 | [13.795, 13.805) | **IN** |
| R_∞,₀ raw-GSL [km] | 13.8024398765 | 13.8024397139 | [13.795, 13.805) | **IN** |

## 13. Entire printed ρ interval (S8/S10)

Independently integrated 9 states across `[1.095e15, 1.105e15)`:

| Image | This review | Structure-1 | Agreement |
|---|---|---|---|
| M [M☉] | [0.623147392100, 0.624120126217] | [0.623147391953, 0.624120126175] | ≤1.5e-10 |
| R₀ [km] | [12.7547529615, 12.7816288631] | [12.7547527962, 12.7816286996] | 1.65e-7 (the tail-midpoint bias) |
| R_∞,₀ [km] | [13.7900059205, 13.8147957579] | [13.7900057567, 13.8147955958] | 1.64e-7 |

M rises and R₀ falls monotonically across the interval, confirmed independently.

**Common matching subset — independently determined.** Bisecting on the *exact*
bins with no numerical margins, the continuous window in which **one same central
state** satisfies all three printed bins is

```
rho_c ∈ [1.098939e15, 1.101175e15] g cm^-3
```

— width **2.236e12 g cm⁻³, only 22.4 % of the printed ρ bin**, bounded below by
`R_∞ < 13.805` and above by `R₀ ≥ 12.765`. Structure-1's tested discrete subset
`[1.0990625e15, 1.10109375e15]` is a conservative subset of this (it carries the
extra 1e-7 M☉ / 1e-5 km margins). **The printed ρ_c = 1.10e15 sits at 47.5 % of
the window — essentially dead centre.**

This is the single strongest piece of evidence in the whole exercise: the
simultaneous three-observable match is available on only ~22 % of the printed
density bin, and the source's printed density lands in the middle of it. The
midpoint passes robustly after all uncertainties and both metric conventions.

**No separate density is used for different observables.** The `bins()` predicate
is evaluated on one `Result` from one solve, and neither the bins nor any
constant is widened — the margins *narrow* the bins.

## 14. High-density pre-Σ mass analysis

Recomputed the sub-Σ branch independently (analytic EOS, direct-vacuum stars):

| n_B [fm⁻³] | ρ_c [g cm⁻³] | M [M☉] | R₀ [km] |
|---:|---:|---:|---:|
| 0.590 | 1.071537e15 | 0.620807150 | 12.845841 |
| 0.600 | 1.090601e15 | 0.622714920 | 12.793542 |
| 0.610 | 1.109690e15 | 0.624571368 | 12.742247 |
| **0.615** | 1.119243e15 | **0.625480892** | 12.716965 |
| 0.617 | 1.123066e15 | 0.625841278 | 12.706919 |
| 0.6173552 | 1.123745e15 | 0.625905080 | 12.705139 |
| **Σ⁻ ceiling 0.617355208** | **1.123745e15** | **0.625905081** | **12.705138** |

`M(0.615) = 0.625480892` independently confirms Structure-1's `0.625480949269`
(difference 5.7e-8 M☉, fully explained by the recorded ≤6.03e8 g cm⁻³
central-state recovery offset on the high branch). **The mass is still rising at
every sampled state, including the last one below the ceiling.**

**One-sided sub-Σ mass limit (not a source-authenticated M_max):**
**M → 0.6259051 M☉ as ρ_c → ρ_Σ⁻ = 1.123745e15 g cm⁻³.** There is no interior
maximum below Σ⁻; the supremum is attained only in the limit at the domain
boundary.

**The hypothetical "true supremum" row in our model:**

| | M | ρ_c | R | R_∞ |
|---|---|---|---|---|
| At the Σ⁻ ceiling | 0.625905 → **0.63** | 1.123745e15 → **1.12** | 12.7051 → **12.71** | 13.7442 → **13.74** |
| **Source prints** | **0.62** | **1.10** | **12.77** | **13.80** |
| At the printed ρ_c = 1.10e15 | 0.623636 → **0.62** ✔ | 1.10e15 ✔ | 12.7682 → **12.77** ✔ | 13.8024 → **13.80** ✔ |

**Not one of the four printed values matches the Σ-threshold configuration; all
four match at the printed ρ_c.**

## 15. R_∞ convention analysis

**Formula verified.** `R_∞ = R e^{−Φ_s} = R/√(1 − 2GM/Rc²)` is exactly FR2005
p. 3. The production profile route (`NStar::BuildFromTOV` surface lapse) agrees
with the explicit Schwarzschild identity to 1e-14 (S9), which I confirm.

**Origin of the 7.18 cm offset — identified exactly.** It is *entirely* the
solar-mass constant:

- `SingleStarSolveToTOVPoints` emits `M = y[0]/GSL_CONST_CGSM_SOLAR_MASS`
  (`TOVSolver.cpp:2591`), and integrates the ODE with `G = 6.673e-8` (GSL).
- `NStar` then maps that number to km with `Zaki::Physics::SUN_M_KM`
  (`NStar.cpp:182, 767, 782`).
- `G·M☉_GSL/c² = 1.4767161818921164 km` vs `SUN_M_KM = 1.4766250380501249 km`
  — **6.1725e-5 relative**.
- Propagating: `d ln R_∞ = ½·(2GM/Rc²)/(1−2GM/Rc²)·d ln M = 0.0843 × 6.1725e-5
  = 5.20e-6`, i.e. `13.80 × 5.20e-6 = 7.18e-5 km`. **Exactly the reported
  7.18 cm**, and my two independent values differ by 7.180e-5 km. Confirmed.

**Adjudication.**

1. *Is separating a deterministic convention offset from numerical convergence
   scientifically legitimate?* **Yes.** It is a reproducible, sign-definite bias
   with an identified algebraic cause, not scatter. Folding it into a
   convergence envelope would misrepresent both. Structure-1 reports it
   separately and does not absorb it — correct handling.

2. *Which convention should govern the FR2005 comparison?* **The raw-GSL one.**
   `2·solar_length_GSL·M/R` is algebraically identical to `2Gm/(Rc²)` with the
   *same* `G` used to build the star, so it reproduces the computed star's true
   compactness. The profile route mixes a GSL-normalised mass number with a
   different `GM_⊙/c²`, which is equivalent to evaluating the metric with an
   effective `G = 6.67259e-8` — not the `G = 6.673e-8` that produced the star.
   Structure-1 designates the *profile* value as principal; on the physics, that
   is the less defensible of the two. **Both pass**, so this changes no verdict.

3. *Does the 7.18 cm block the reproduction claim?* **No.** The bin margins are
   ~2.4e-3 km, 33× the offset. Both values round to 13.80.

4. *Technical debt?* **Yes, and pre-existing** — recorded in the parent
   commit's interface audit and living in `NStar.cpp`, which this change does not
   touch. It is correctly outside this task's scope.

## 16. S1–S12 classification

| ID | Classification | Assessment |
|---|---|---|
| **S1** | **Genuinely independent / load-bearing** | GL-96 phase-space quadrature vs closed form + series, plus the enthalpy identity restricted to well-conditioned subtraction. I verified the GL-96 rule is itself machine-accurate (5.6e-16 – 1.8e-15 relative for x from 0.1 to 440 — my initial concern about near-axis branch points was unfounded). Reproduced: energy 2.0910e-11, pressure 2.8908e-12, identity 3.9168e-12. *Gap:* the 15-point sample set contains no species' `x ≈ 0.1` switch, so the 1.96e-12 branch-switch loss is not exercised (I measured it separately). |
| **S2** | **Genuinely independent, but not load-bearing for this claim** | Both halves are real (differentiated quadrature 4.7393e-6; interpolant slope 0.12808 → 0.0040951, strictly decreasing). But `dε/dP` is provably absent from the static path (§9), so S2 constrains a future consumer, not M/R/R_∞. |
| **S3** | **Genuinely independent / load-bearing** | Vacuum, exact and `nextafter` domains, 201-point threshold bridges vs the oracle (neutron 3.0e-13 / 5.164e-10 / 3.657e-10; muon ≤1e-14), one-sided slope approach, and the P-7 rows with H refusal. Reproduced exactly. The one-sided-slope check only requires monotone decrease (no rate) — appropriate given the √δ cusp. |
| **S4** | **Load-bearing, and stronger than its samples** | I independently **re-proved** the analytic `0 < dP/dε ≤ 1/3` bound: with `s_i = 3n_iμ_i/k_i² ≥ 3n_i/μ_i`, `1/s_p + 1/L ≤ μ_n/(3n_p)` gives `hC ≥ 3n_B`, hence `dP/dε = n_B/(hC) ≤ 1/3` on every branch including the one-sided threshold limits. The proof, not the sampled 0.0759, supplies the global bound — as the record says. Identities reproduced (1.1435e-13). |
| **S5** | **Path-equivalence only** — correctly self-classified in the parent audit as "not independent physics validation". Bit-identical columns; useful anti-drift regression. |
| **S6** | **Regression / convergence** | Five partitions × four tables; observed M spread 1.05e-11, R_cut spread 6.18e-10 km, lapse 2.52e-12 — far inside the unwidened ADR-0009 contract. Genuine but not independent physics. |
| **S7** | **Genuinely independent / load-bearing** | Seven floor levels, analytic enthalpy bounds derived from the exact EOS, plateau at the unchanged coefficient. All reproduced. *Weakness:* see S12 gap below. |
| **S8** | **Mixed — one assertion is tautological** | `minR`/`maxR` are attained at the *fixed* interval endpoints, identical in all three grids, so `|minR − previous_min| < 1e-7` is bit-exactly `0` and **can never fail**. The load-bearing S8 evidence is `curve_error`, the midpoint-insertion residual (1.4183e-7 → 3.5586e-8 km), which is genuine, plus the monotonicity requirement. Reclassify the endpoint check as regression, not convergence evidence. |
| **S9** | **Useful, partially shared** | Profile lapse vs the explicit identity at 1e-14 is a genuine internal check; it cannot detect the constant-convention issue it sits next to, because both sides use `SUN_M_KM`. The raw-GSL value is reported but not asserted. |
| **S10** | **Load-bearing** | One predeclared state, three bins, no widening. I confirmed the bins are the correct ordinary-rounding intervals and that the margins *narrow* them. Mutation testing (§17) shows real discriminating power. |
| **S11** | **Genuinely load-bearing** | Σ-ceiling and generator guards verified by execution; clamp margin 1.5033e13 g cm⁻³ above the printed interval; no super-Σ padding. |
| **S12** | **Genuinely independent / load-bearing — the strongest falsifier** | Independence audited line by line (§10) and confirmed by a *third* implementation. Mass convention verified consistent, which the passing test alone would not show. *Gap:* the canonical/oracle gate is applied only to grid-8192 radial rows and rho-bin rows — **not** to the floor sweep or the 17 high-branch rows. |

**Missing material falsifier (one, identified).** §16 asserts narratively that
"the independent full enthalpy star lies inside the bracket", but **no
executable assertion tests it**: line 97 only self-refines the vacuum oracle, and
lines 122/188 compare canonical to oracle at the *finite cutoff*. I verified the
claim manually and it holds (12.7681549010 ∈ [12.7681545694, 12.7681549024]) —
and this is precisely the check that would have surfaced the midpoint bias of
Finding 3.

**No "PASS" label was accepted without reading its assertion**; every numeric
claim in §8–§20 of the record was regenerated by execution.

## 17. Numerical-budget audit

I measured the **total** production error directly, using `R₀_upper` (shown in
§11 to be near-tight) against my independent analytic-EOS vacuum star:

| grid | R₀_upper − R₀(independent) | M − M(independent) |
|---:|---:|---:|
| 1024 | −7.051e-5 km | +2.528e-6 M☉ |
| 2048 | −2.110e-7 km | +5.382e-9 M☉ |
| 4096 | −1.907e-8 km | +5.583e-10 M☉ |
| **8192 (accepted)** | **+2.418e-9 km** | **−9.414e-11 M☉** |

The residual at the accepted table is at the noise floor of the comparison. So:

| Claimed envelope | Directly measured | Margin |
|---|---|---|
| 9e-6 km total radius | **2.4e-9 km** | ~3700× conservative |
| 1e-7 M☉ | **9.4e-11 M☉** | ~1000× conservative |
| 1e9 g cm⁻³ central density | **9.26e5 g cm⁻³** (fine grid) | ~1000× conservative |

**Component check — nothing omitted or double-counted.** The seven allocations
cover equilibrium values/pressure (measured ≤4e-13 relative → ≲1.4e-12 km),
EOS interpolation (6.3e-7 relative off-grid, oscillating in sign, aggregate
effect bounded by the table ladder above), central inversion (9.26e5 g cm⁻³ →
2.5e-9 km; high branch 6.03e8 → 1.6e-6 km, inside its 3e-6 allocation), radial
integration (measured spread 6.2e-10 km), central sampling (3.56e-8 km),
surface tail (half-width 1.665e-7 km, and the actual bias is 1.65e-7 km — the
5e-7 allocation covers it), and metric arithmetic. The R_∞ propagation
(`dR_∞/dR ≈ 1.081`, mass contribution <2e-7 km) is handled.

**Verdict: these are evidence-based envelopes, not merely selected tolerances.**
They are conservative by three orders of magnitude, and the 7.18 cm convention
offset is correctly excluded from them rather than hidden inside.

**Discriminating power of the bins (mutation testing).** I perturbed my
independent model and re-ran the whole star:

| Mutation | M | R₀ | All bins |
|---|---|---|---|
| baseline | 0.623636 | 12.76815 | PASS |
| **omit muons (npe only)** | 0.623962 | 12.76216 | **FAIL** |
| **pressure × 2 (degeneracy defect)** | 0.588054 | 12.88016 | **FAIL** |
| **ε × (1 + 1e-3)** | 0.623255 | 12.76461 | **FAIL** |
| ε × (1 + 1e-4) | 0.623598 | 12.76780 | PASS |
| P × (1 + 1e-3) | 0.623597 | 12.76827 | PASS |

The comparison genuinely discriminates against omitted species, degeneracy and
factor errors. Its resolution is ~1e-4 relative in ε — the intrinsic limit of a
two-significant-figure published comparison. **The binding evidence for
correctness is therefore the independent oracles (1e-11–1e-15), not the bins**,
which is the right structure and is what the record claims.

## 18. Regression / baseline evidence

| Check | Result |
|---|---|
| `TOVSolver.cpp` / `.hpp` | **byte-identical to master** |
| `NStar.cpp` / `.hpp` | **byte-identical to master** |
| `StarProfile.cpp` / `.hpp` | **byte-identical to master** |
| `RotationSolver.cpp`, `ThermalSolver.cpp` | **byte-identical to master** |
| Evolution / BNV code | **untouched** |
| Governed baselines | **8**, all SHA-256 unchanged, `git diff` empty |
| Ninth baseline | **NO** |
| Tests removed or reclassified | **NONE** |
| Full test inventory | **51** tests (authenticated by `ctest -N`) |
| Data-free test inventory | **28** tests |
| Results-JSON hashes (11 source + 4 protected + 8 baseline = 23) | **all match the working tree** |

**Focused rerun** (`ctest -R "rotochemical_trackr|rotochemical_local|eos_derivative_contract|tov_reference_analytic"`):
**8/8 passed**, 5.47 s — including `rotochemical_trackr_freegas_barotrope`
(0.39 s), `rotochemical_trackr_freegas_structure` (4.10 s),
`rotochemical_trackr_freegas_local`, `rotochemical_trackr_npe`,
`rotochemical_trackr_pe`, `rotochemical_local_thermodynamics`,
`eos_derivative_contract`, `tov_reference_analytic`.

**Complete data-free suite rerun:** configured and built fresh into
`build-selfcontained/`; **28/28 passed in 95.57 s** (record: 28/28 in 94.87 s).

**Full 51-test suite: not rerun.** Justified under the review mandate: every
focused result is identical, the inventory is authenticated at 51, the protected
files are byte-identical, the baselines are unchanged, and the data-free suite
reproduced. There is no reason to distrust the recorded 51/51 in 670.29 s.

The working tree remained clean throughout; all my run artifacts were written to
a session scratchpad, none into the repository.

## 19. Ranked remaining findings

**BLOCKER BEFORE STRUCTURE-1 ACCEPTANCE — none.**

**SOURCE AMBIGUITY**

1. **Source `M_max` semantics.** Adjudicated in §21 below. The label cannot be
   reconciled with a continuous supremum in this model; the printed row is a
   sub-threshold configuration at ρ_c = 1.10e15. Irreducible from the paper.
2. **Unknown source retained mesh / exact author input.** FR2005 publishes no
   scan grid, full-precision input, constants inventory, surface algorithm or
   error budget. Not resolvable. Correctly disclaimed by Structure-1.

**MUST FIX BEFORE THE NEXT GLOBAL ROTOCHEMICAL STAGE**

3. **Threshold-localised derivative interpolation error (0.410 %).** Harmless
   here (§9), but `dEdP` is published into `StarProfile` and consumed by
   `RotationSolver`; the underlying √δ cusp at neutron onset means table
   refinement alone converges slowly. Any future consumer of `dEdP` — Hartle
   O(Ω²) response, and the eventual B̃/Z/W integrals — needs its own accuracy
   analysis in this region before it can rely on the published column.
4. **INV-09 unresolved** (fixed-ε_c vs equilibrium-sequence derivatives) and
   **INV-11 unresolved** (chemical-imbalance convention). Correctly recorded as
   unresolved; neither is touched or silently closed by Structure-1. Both block
   the particle-number spin-down response, not this static result.

**NONBLOCKING HARDENING**

5. **Tail-bracket midpoint is a biased estimator.** `r_upper` is near-tight
   (2.4e-9 km from truth), `r_lower` is loose, so the midpoint under-reports R₀
   by 1.65e-7 km — essentially the full half-width. Immaterial (60× below the
   comparison margin). Prefer `r_upper` or the direct-vacuum value.
6. **Missing executable falsifier:** no assertion places the direct-vacuum
   oracle inside `[R₀_lower, R₀_upper]`. Verified true by hand; should be a test.
7. **S8 endpoint-image convergence assertion is tautological** (bit-exactly zero
   by construction). Reclassify; the genuine evidence is `curve_error`.
8. **S12 gate coverage gap:** canonical/oracle agreement is not asserted on the
   floor-sweep or high-branch rows.
9. **Pressure closed form loses ~4 digits for `x ∈ [0.1, ~0.15]`** (measured
   1.96e-12 relative). Immaterial to structure (≲2e-16 of the star's pressure
   scale), not exercised by S1's sample set. Consider extending the series range
   or a cancellation-free closed form.
10. **§17's near-onset neutron-density errors are ~170× optimistic** relative to
    a fully-exact-physics reference (measured 2.42e-7 / 2.89e-5 / 3.08e-4 at
    offsets 1e-7 / 1e-9 / 1e-10). The mechanism §17 names is right; the
    comparison convention is unstated. Impact ≤3e-14 relative in ε and P.
11. **`finest_midpoint_oracle_disagreement` is a max over partitions**, not the
    finest partition (which gives 1.18e-11 / 7.54e-11). Conservative, so nothing
    is overstated; the field name is imprecise.
12. **Adjacent-double threshold ambiguity — confirmed and contained.** At the
    muon onset, P at `nextafter(onset,∞)` is 2.1e-14 MeV fm⁻³ *below* P at the
    onset. Disclosed in §9. Contained: the generator's minimum relative distance
    is 1e-11 (neutron) / 1e-8 (muon), it asserts strict monotonicity, and the
    emitted table is strictly increasing in P, ρ and n_B.

**NO ISSUE**

13. **Raw-GSL / Zaki constant-convention debt (7.18 cm).** Pre-existing, in
    protected `NStar.cpp`, recorded in the parent audit, correctly out of scope
    here — with the caveat in §15 that the *raw-GSL* value is the more defensible
    one for the source comparison. Both pass.
14. **`Values()` computes `k` in `double` while `Evaluate()` uses `long double`**,
    so the two production entry points can differ by a few ulps in ε for the same
    n (~1e-16 relative). Cosmetic.
15. **Residual P-7 / value-resolution limit: none found.** Value resolution never
    failed anywhere in the representable domain.

## 20. Level-1 / Level-2 / Level-3 claim adjudication

**LEVEL 1 — ESTABLISHED.** The static implementation is numerically correct and
independently validated. Evidence: analytic re-derivation of every formula
(pressure closed form and series, npe δ-parameterisation, npe-µ inversion,
enthalpy TOV RHS, center expansion, tail bounds, the ≤1/3 causality proof);
60-digit agreement on thresholds (≤2.2e-14), central state (4e-18) and the EOS
(≤4.1e-13); three-way structural agreement to 2.4e-9 km and 9.4e-11 M☉;
convergent table/partition/cutoff ladders; a measured total error 3 orders below
the claimed budget.

**LEVEL 2 — ESTABLISHED.** A source-compatible central state inside the printed
ρ interval — specifically the *printed* value ρ_c = 1.10e15 g cm⁻³, predeclared
before any structure result existed — reproduces the printed M, R and R_∞
simultaneously, in **both** metric conventions, with no bin widened and no
constant tuned. Strengthened by the finding that the simultaneous match holds on
only ~22 % of the printed ρ bin, with the printed ρ_c essentially at its centre.

**LEVEL 3 — NOT ESTABLISHED, and independently falsified under the supremum
reading.** In this model the sub-Σ⁻ mass supremum is 0.6259051 M☉ at
ρ_c = 1.123745e15, R = 12.7051 km, R_∞ = 13.744 km — printing as
0.63 / 1.12 / 12.71 / 13.74, matching **none** of the four printed values. The
source's `M_max` label therefore cannot mean the continuous supremum. Whether it
means "the largest mass in the authors' retained sequence" (their mesh is
unpublished) or "a deliberately chosen representative sub-threshold state" is not
resolvable from the paper, so Level 3 is not reproducible in principle from the
published information.

Level 2 does **not** imply Level 3, and this review does not treat them as
equivalent.

## 21. The §2 adjudication — what has actually been reproduced

**Q1. What does the source's "maximum mass before appearance of Σ− hyperons"
mean?**

**Answer: C** — a representative sub-threshold configuration whose central
density is reported as 1.10e15 — **with B as the plausible but unstated
mechanism** by which that particular state was selected. **A is excluded**, and
D applies only to the selection mechanism, not to the semantics.

This is settled by the source's own numbers, not by inference from the existence
of a numerical mesh:

- FR2005 §3.1 p. 9 states directly that the tabulated values *"correspond to a
  central density slightly below the appearance of Σ− hyperons"* — the values are
  **at** a specific sub-threshold density, not at a limit.
- Independently, ρ_Σ⁻ = **1.123745e15 g cm⁻³** in this model. A supremum
  configuration would have been printed as ρ_c = **1.12**, M = **0.63**,
  R = **12.71**, R_∞ = **13.74**. The paper prints 1.10 / 0.62 / 12.77 / 13.80.
  Reading A is inconsistent with all four printed entries.
- Footnote c is then a **domain-truncation label**: within the model's validity
  region (no hyperons), the mass increases monotonically to the boundary, so
  "the maximum mass before Σ⁻" is a statement about *which branch was retained*,
  not about a stationary point. There is no `dM/dρ_c = 0` anywhere below Σ⁻ —
  I verified the mass rises monotonically all the way to the ceiling.

**Q2. Does our sub-Σ sequence reaching M > 0.625 conflict with the printed
0.62?**

**It conflicts with reading A only; it does not conflict with the printed row.**
The printed row is a (ρ_c, M, R, R_∞) tuple. At the printed ρ_c our M is
0.6236 → **0.62**, matching. Our M(0.615 fm⁻³) = 0.6255 → 0.63 belongs to a
*different* central state that the source does not report. The gap between the
printed state and the ceiling is only **0.00227 M☉ (0.36 %)** — the free-gas
M(ρ_c) curve is very flat here, which is exactly why a modest mesh truncation or
a deliberate "slightly below" choice shifts the printed mass across the 0.625
rounding boundary.

**Q3. Can convention differences account for the discrepancy? Quantified — no.**

| Candidate | Effect on M | Enough to explain 0.0023–0.0059 M☉? |
|---|---|---|
| Solar-mass convention (6.17e-5 rel.) | 3.9e-5 M☉ | No — 60× too small |
| `G` (GSL 6.673e-8 vs CODATA 6.67430e-8; M ∝ G^{−3/2}) | ~2.4e-4 M☉ | No — 10× too small |
| Particle masses m_n, m_p, m_e, m_µ | ≪1e-5 M☉ | No |
| Surface convention (R_cut vs R₀) | ~1e-9 M☉ | No |
| Source numerical precision | unknown, but the printed row is self-consistent to 4/4 | No |
| **m_Σ⁻ (to move the ceiling to ρ_c = 1.10e15)** | would need **−3.031 MeV (0.253 %)**; PDG uncertainty is ~0.03 MeV | **No — 100× beyond any plausible value** |

No combination of plausible conventions turns the printed row into the
Σ-threshold configuration, and none moves the Σ⁻ threshold down to the printed
density.

**Q4. What is the defensible claim?**

> **The printed FR2005 Table-1 free-gas row is reproduced at a common
> source-compatible central state — the source's own printed ρ_c = 1.10e15
> g cm⁻³ — with M, R and R_∞ simultaneously inside their printed rounding bins in
> both available metric conventions. The source-qualified pre-Σ⁻ `M_max`
> supremum is *not* reproduced, and is not reproducible from the published
> information.**

The stronger Structure-1 headline — *"FR2005 FREE-GAS STATIC TABLE-1 STRUCTURE
REPRODUCED"* — is **acceptable only as qualified**, and the record does in fact
qualify it in place (§2, §18, §22, §25 all state the boundary explicitly, and
the document header disclaims identifying the authors' input or mesh). It is not
scientifically invalid; it is a Level-2 claim whose Level-3 gap is disclosed. It
should not be promoted, restated, or cited without that qualification.

No source bin was broadened to rescue anything; on the contrary, the executable
bins are *narrower* than the printed rounding intervals.

## 22. Final disposition

> **B — TRACK-R STRUCTURE-1 ACCEPTABLE WITH QUALIFIED CLAIM —
> STATIC FREE-GAS STRUCTURE AND COMMON-STATE TABLE-1 NUMBERS VERIFIED;
> READY FOR HUMAN RATIFICATION AND CANONICAL INTEGRATION**

Level 1 and Level 2 are established by independent derivation, independent
high-precision numerics and an independent third TOV implementation. Level 3 is
not established, because the source-qualified `M_max` meaning refers to a
retained sub-threshold configuration whose selection rule the paper does not
publish — and, under the supremum reading, is positively inconsistent with all
four printed entries.

Disposition A is withheld solely because the `M_max` qualification is material
and cannot be discharged from the source. Disposition C is not warranted: no
blocking physics, numerical, surface, source-comparison or independent-oracle
failure was found, the protected surfaces and baselines are untouched, and no
production or test change is scientifically required.

This remains an AI-authored scientific candidate under `GOVERNANCE.md` §5.
Human ratification is not implied by this review. No merge is performed.
