# Phase 5A-5 — Independent review of final Track-R free-gas local thermodynamic coverage

> **PHASE 5A-5 ACCEPTABLE WITH NONBLOCKING FOLLOW-UP —
> LOCAL COVERAGE COMPLETE; READY FOR HUMAN RATIFICATION AND CANONICAL INTEGRATION**
>
> The narrow claim — that the local source model supplies the correct physical
> state and an honest active response at every equilibrium density from the
> `P=0` surface to the `Sigma-minus` ceiling — is **justified**. The active-set
> ladder, the p-e thermodynamics, the vacuum boundary, the neutron threshold,
> the p-e/npe limiting identity and the N-1 fail-closed policy were each
> re-derived from first principles and from the local source PDFs, and checked
> against independent 90-digit arithmetic. No blocking local scientific or
> numerical defect was found.
>
> This review is **evidence, not ratification**. Per `GOVERNANCE.md` §5 the
> reviewed code remains an **UNVERIFIED SCIENTIFIC CANDIDATE** until a human
> with domain authority ratifies it. No production or test file was modified by
> this review.

## 1. Reviewed SHA

| Item | Value |
|---|---|
| Reviewed HEAD | `933494d86daf2cf8965079ece49fabd66d9390e5` |
| Commit subject | `feat: complete Track-R low-density thermodynamics` |
| Branch | `physics/phase5a-trackr-pe-branch` |
| Review worktree | `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-rotochemical-trackr-pe` |
| Canonical master | `41ab66bd6e6f351691173a1e2a033c646ffd3772` |
| Merge base | `41ab66bd6e6f351691173a1e2a033c646ffd3772` |
| Change class reviewed | scientific-semantic + numerical policy + structural/API + test/build registration + documentation |
| Change class of this review | documentation |

## 2. Authentication

Performed before any file was read for scientific content.

| Check | Result |
|---|---|
| `git status --porcelain` | empty (clean, including untracked) |
| Local HEAD | `933494d86daf2cf8965079ece49fabd66d9390e5` |
| Local upstream `@{u}` | `933494d86daf2cf8965079ece49fabd66d9390e5` — **equal** |
| Live remote `refs/heads/physics/phase5a-trackr-pe-branch` | `933494d86daf2cf8965079ece49fabd66d9390e5` — **equal** |
| Local `master` | `41ab66bd6e6f351691173a1e2a033c646ffd3772` |
| Live remote `refs/heads/master` | `41ab66bd6e6f351691173a1e2a033c646ffd3772` — **unchanged** |
| `git merge-base HEAD master` | `41ab66bd6e6f351691173a1e2a033c646ffd3772` |
| `git rev-list --left-right --count master...HEAD` | `0 1` — **exactly 1 ahead / 0 behind** |
| Diff scope | 10 files: 3 production, 3 test/build, 4 documentation |

All four expectations stated in the review task hold exactly.

## 3. Source ledger

The three governing PDFs were read directly from the shared read-only library
and their SHA-256s recomputed live; all three match the manifest quoted in the
implementation record (`PHASE5A5_TRACKR_PE_BRANCH.md:55-57`).

| File | SHA-256 | What this review took from it, verbatim where it matters |
|---|---|---|
| `2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf` | `f184d7d1…069ea88` | §3.1 EOS list, item 3: **"a non-interacting Fermi gas (see, e.g., Shapiro & Teukolsky 1983) for the whole star."** Items 1–2 (APR, BPAL) are explicitly *"supplemented with"* separate inner/outer crust EOSs; item 3 is not. §3.4 (eq. 50): "The sum is performed over all particle species and the integral is done over the region where each particle species exists. We take only free particles into account, neglecting the heat capacity due to the lattice of ions in the crust." Table 1's free-gas entry is at a central density "slightly below the appearance of Σ− hyperons" |
| `2006-Reisenegger-…-Electrostatic-Potential-Perturbations.pdf` | `a286f15e…5c4d23ae8` | Eq. (2): the total redshifted potential contains `q_i φ`. Eq. (4): local charge neutrality `Σ n_i q_i = 0` replaces the Poisson equation. §3: the imbalances `η_npl = μ_n − μ_p − μ_l` are intrinsic combinations from which `φ` cancels. Eq. (6): the reaction volume integral runs over the core |
| `1995-Reisenegger-Rotochemical-Heating.pdf` | `9af85e37…91d0ff` | §2: `n_p = n_e` "in order to preserve charge neutrality"; the equilibrium composition "is found by minimizing the total energy per baryon `E(n,x,s)` (including the contribution from the electrons) with respect to `x`" — the variational principle used in §4 |

ADR-0010 remains the accepted contract. Nothing in the reviewed diff adds an
electrostatic variable, a second projection, a GR term, or a stellar integral,
so the 2006 ownership split is untouched.

**Constant authority probed independently**, not read from any record: a scratch
program linked against `libZaki.a` printed `m_n = 939.56542051999997511`,
`m_p = 938.27208815999995295`, `m_e = 0.51099894999999995182`,
`m_μ = 105.65837550000000533`, `m_Σ− = 1197.4490000000000691`,
`hbar c = 1/MEV_2_INV_FM = 197.32698045930249009`, and
`__LDBL_MANT_DIG__ = 53` — re-confirming Phase-5A-3 finding **R8** (arm64 macOS
`long double` carries no extra precision). All independent arithmetic below uses
Python `decimal` at 80–90 digits seeded with the **exact binary values** of those
doubles; no repository routine is used as a numerical oracle anywhere in §§4–11.

## 4. Active-set ladder

Re-established from the 1995 §2 variational principle, treating the prior
reviews' conclusion as a hypothesis to be challenged rather than an input.

Minimize `ε = Σ_i ε_i(n_i)` at fixed `n_B > 0` subject to `n_n + n_p = n_B`,
`n_p = n_e + n_μ`, `n_i ≥ 0`. Each `ε_i` is convex, the feasible set is a convex
polytope, so the minimizer is unique. KKT gives `λ_B = μ_n`, `−λ_Q = μ_e = μ_μ`
where those species are active, and `μ_p = λ_B + λ_Q`, i.e. `μ_n = μ_p + μ_e`.

- **Protons are active at every `n_B>0`.** `n_p = 0` forces `n_e = n_μ = 0` and
  requires `m_p + m_e ≥ μ_n(n_B) ≥ m_n`, false because
  `m_n − (m_p+m_e) = +0.78233341 MeV`.
- **Electrons are active whenever protons are.** `n_e = 0` with `n_μ > 0` needs
  `m_e ≥ μ_μ ≥ m_μ`, false.
- **Neutrons are inactive iff `μ_p(n_B) + μ_e(n_B) < m_n`,** with `n_p = n_e = n_B`.
  The left side is strictly increasing and starts below `m_n`, so it crosses once.
- **Muons are inactive iff `μ_e < m_μ`.**

**Independent variational falsifier, assuming no KKT at all.** I scanned
`E(x) = ε_n(x) + ε_p(n_B−x) + ε_e(n_B−x)` over `x = n_n > 0` at fixed `n_B`,
asking directly whether any positive neutron density lowers the energy:

| `n_B/n₀` | any `x>0` with `E(x) < E(0)`? | `dE/dx|₀₊ = m_n − μ_p − μ_e` |
|---|---|---|
| `1e-6`, `1e-2`, `0.5`, `0.9`, `0.99999` | **no** | `+7.82e-1 … +3.64e-6` |
| `1.0` | **no** | `−4e-78` (zero to working precision) |
| `1.00001`, `1.1`, `2.0`, `10.0` | **yes** (`x = 7.36e-14 … 4.14e-8`) | `−3.64e-6 … −1.32e+0` |

The minimizer moves off the boundary **exactly** at `n₀`. A companion test that
replaces `δ` electrons by muons at fixed `n_p` gives `ΔE = +7.7e-16 … +7.7e-10
MeV fm^-3 > 0` everywhere on the p-e branch, where `μ_e ≤ 1.292581 MeV ≪ m_μ`.

**Independently confirmed ladder** for the FR2005 free-particle inventory
`{n,p,e,μ}` with `Σ−` as the declared ceiling:

```text
vacuum  ->  p-e  ->  n-p-e  ->  n-p-e-mu
```

with, as the task asks to confirm:

| Region | Confirmed content |
|---|---|
| `0 < n_B < n₀` | `n_n = 0`, `n_μ = 0`, `n_p = n_e = n_B` — **forced by the two constraints, not by optimization** (zero remaining freedom) |
| `n_B = n₀` | `m_n = μ_p + μ_e` |
| `n₀ < n_B < n_μ` | `η_npe = 0` |
| `n_B = n_μ` | `μ_e = m_μ` |
| `n_μ < n_B < n_Σ` | `η_npe = η_npmu = 0` |

**No additional active-set transition exists between vacuum and neutron
appearance.** Protons and electrons are active at every positive density,
neutrons and muons are the only other species in the model, and each appears
exactly once. Atomic binding, nuclei, a Coulomb lattice or a crust prescription
would define a different source model; FR2005 §3.4 explicitly counts "only free
particles" and assigns no crust supplement to this EOS.

## 5. Independent p-e derivation

On `0 < n_B < n₀` the constraints leave one independent direction, so

```text
z_pe = (n_B)                              [fm^-3]
epsilon_pe(n_B) = epsilon_p(n_B) + epsilon_e(n_B)   [MeV fm^-3]
h_pe = d epsilon_pe / d n_B = mu_p + mu_e           [MeV]
H_pe = d h_pe / d n_B = D_p + D_e                   [MeV fm^3]
```

**All four required properties confirmed.**

- **Dimension is genuinely 1.** It is not a projection or a padded 2×2: the
  physical manifold at `n_n = n_μ = 0` *is* one-dimensional, because neutrality
  and baryon accounting each remove one of the three canonical coordinates.
- **Units.** `h_pe` is a chemical potential (MeV); `H_pe` is `MeV fm^3`, the same
  units as every other entry of the ADR-0010 primitive.
- **`H_pe > 0` throughout.** `D_i = p_{Fi}²/(3 μ_i n_i) > 0` for every positive
  density, so `H_pe = D_p + D_e > 0` **analytically**, not merely on a grid. It
  is a genuine 1×1 positive-definite response.
- **No inactive derivative belongs here.** `D_n` and `D_μ` are undefined at zero
  density (both diverge as `n^{-1/3}`); neither appears in `PeConjugates` or
  `PeChemicalHessian`. The inactive-neutron condition is carried only as the
  separately named value diagnostic `eta_npe_threshold_diagnostic_MeV = m_n − h_pe`.

**Independent numerical comparison** — production `EvaluatePe` versus separately
coded 90-digit ideal-gas formulae (`p_F = ħc(3π²n)^{1/3}`, `μ = √(m²+p_F²)`,
`ε = m⁴[x√(1+x²)(1+2x²) − asinh x]/(8π²(ħc)³)`, `D = p_F²/(3μn)`), at 13
densities spanning **nine decades** from `1e-15 n₀` to `0.999999 n₀`:

| Quantity | Worst relative error over all 13 densities |
|---|---|
| `epsilon_pe` | **6.32e-16** |
| `h_pe` | **1.02e-16** |
| `H_pe` | **3.52e-16** |
| composition (`n_p = n_e = n_B`, `n_n = n_μ = 0`) | **exact** (bit-equal at every density) |

The inactive-neutron diagnostic is the one quantity whose *relative* accuracy
degrades — to `2.5e-7` at `0.999999 n₀` — because it is a cancellation of two
`~939 MeV` numbers producing `3.64e-7 MeV`. Its **absolute** accuracy stays at
the `~1e-13 MeV` roundoff floor, and PE-V5 correctly bounds it absolutely
(`3e-13 MeV`). Recorded as finding **P-6**, not a defect.

## 6. Vacuum review

At `n_B = 0` the physical facts are unambiguous: every species density is zero
and `ε = 0`. The one-sided limits along the p-e branch, computed independently:

| `n_B/n₀` | `epsilon_pe` | `h_pe − (m_p+m_e)` | `H_pe` |
|---|---|---|---|
| `1e-6` | `6.9064e-12` | `1.3799e-04` | `1.25027e+10` |
| `1e-12` | `6.9064e-18` | `1.3801e-08` | `1.25061e+12` |
| `1e-18` | `6.9064e-24` | `1.3801e-12` | `1.25061e+14` |
| `1e-24` | `6.9064e-30` | `1.3801e-16` | `1.25061e+16` |
| `1e-30` | `6.9064e-36` | `1.3801e-20` | `1.25061e+18` |

So `ε → 0`, **`h_pe → m_p + m_e` finitely**, and **`H_pe → +∞` as `n_B^{-1/3}`**
(each factor-`10^6` density drop multiplies `H_pe` by exactly `100`). The
divergence is analytic: `D_i = p_F²/(3μ_i n) ∝ n^{2/3}/n = n^{-1/3}`.

**Therefore returning value information, `response_dimension = 0`, no Hessian,
and an explicit `VacuumBoundary` status is the correct representation.** A finite
1×1 vacuum Hessian would be a fabrication; omitting the values would discard a
finite, physically meaningful one-sided limit. Distinguishing `VacuumBoundary`
from `SpeciesThreshold` is also right: this boundary is the edge of the density
domain, not a species-appearance point.

**No code path fabricates a finite vacuum derivative** — verified by inspection
(`VacuumBoundaryEvaluation` has no `hessian` member;
`static_assert(!HasHessian<VacuumBoundaryEvaluation>::value)`) and by live probe:

```text
VAC eps=0  h=938.78308710999999676  diag=0.78233340999997835752  dim=0
    n_B=0 n_n=0 n_p=0 n_e=0 n_mu=0
```

`h` equals the double `m_p + m_e` and `diag` equals `m_n − (m_p+m_e)`.

**Negative-density and invalid-input rejection**, live-probed:
`EquilibriumAt(-1e-30)`, `(-1.0)`, `(NaN)`, `(+inf)` are all rejected;
`EvaluateActive({0.0, 1e-300, 0.0})` — a non-vacuum composition at `n_B = 0` — is
rejected rather than silently snapped to the vacuum result; `EvaluatePe(0.0)` and
`EvaluatePe(-1e-30)` are rejected. `EquilibriumStateAt(0.0)` fails closed with a
message directing callers to the typed boundary result, which preserves the
pre-existing strictly-positive `ChargeNeutralCompositionState` invariant instead
of weakening it — the reason `VacuumCompositionState` is a separate type.

## 7. Neutron-threshold derivation

At neutron appearance `n_p = n_e` and both species have the same degeneracy, so
their Fermi momenta coincide at `p_*`, and `√(m_p²+p_*²) + √(m_e²+p_*²) = m_n`
closes in one step:

```text
mu_e,* = (m_n^2 - m_p^2 + m_e^2)/(2 m_n),  p_* = sqrt(mu_e,*^2 - m_e^2),
n_0    = p_*^3 / (3 pi^2 (hbar c)^3).
```

Evaluated independently at 90 digits from the exact binary constants:

```text
mu_e,*                    = 1.2925811676744569648806698E+0  MeV
n_B,n-onset               = 7.3567289037328326656351713E-9  fm^-3
mu_p(n_0)+mu_e(n_0) - m_n = -1e-87   (zero to working precision)
```

| | Value | Relative agreement |
|---|---|---|
| **Independent (90 digits, this review)** | `7.3567289037328326656351713e-9` | — |
| Production `NeutronOnsetBaryonDensityFm3()` | `7.356728903732831754e-9` | **1.24e-16** (≈ 0.9 ULP) |
| Task reference value | `7.3567289037328326656352e-9` | `< 1e-22` |

**The threshold is independently verified.** My 90-digit value reproduces the
task's reference digit for digit; it was recomputed here from the source
condition rather than carried over, so the chain is not circular.

## 8. Independent onset pressure

`P_i = μ_i n_i − ε_i` at `T = 0` (equivalently `P = n dε/dn − ε`), with the
explicit relativistic ideal-gas energy integral independently coded at 90 digits.
At the onset `n_n = 0` contributes nothing:

```text
P(n_B,n-onset) = 1.8964875026317865961252521E-9  MeV fm^-3   > 0
```

matching the task reference `1.8964875026317866e-9` exactly.

**`P` is positive for every `n_B > 0` on the p-e branch, strictly increasing, and
reaches zero only as `n_B → 0`:**

| `n_B/n₀` | `P` (MeV fm^-3) | `P/n_B` (MeV) |
|---|---|---|
| `1e-15` | `4.0611e-34` | `5.5202e-11` |
| `1e-12` | `4.0611e-29` | `5.5202e-9` |
| `1e-9` | `4.0611e-24` | `5.5202e-7` |
| `1e-6` | `4.0603e-19` | `5.5192e-5` |
| `1e-3` | `3.9852e-14` | `5.4171e-3` |
| `0.1` | `6.5177e-11` | `8.8596e-2` |
| `0.9` | `1.6340e-9` | `2.4678e-1` |
| `1.0` | `1.8965e-9` | `2.5779e-1` |

`P/n_B → 0` as well, confirming `P` vanishes faster than linearly (`P ∝ n_B^{5/3}`
for the degenerate gas). **The p-e branch therefore genuinely reaches the
zero-pressure stellar surface, and the free-gas star's surface lies on it, five
orders of magnitude in density below neutron appearance.** This is the fact that
makes the p-e branch mandatory rather than optional for a whole-star EOS.

## 9. Threshold semantics

Live probe of `EquilibriumAt(n₀)`:

```text
NTHRESH n_B=7.356728903732831754e-09  n_n=0  n_p=n_e=7.356728903732831754e-09
        n_mu=0  eps=6.9102315985847432509e-06
        h=939.56542051999997511 (= double m_n, bit-exact)
        diag=0  dim=0
```

Every required property holds exactly: `n_n = 0`, `n_μ = 0`, `n_p = n_e = n_B`,
and `m_n − μ_p − μ_e = 0`.

**Returning values with `response_dimension = 0` and no smooth 1×1 or 2×2
Hessian is correct.** The state is a genuine point of the manifold with finite
energy and a finite one-sided p-e conjugate, so the values are real physics. But
the active chart itself changes dimension there (1-D below, 2-D above), and the
neutron derivative `D_n` diverges from above: no single smooth response exists at
that point. A 1×1 result would assert the p-e chart still applies; a 2×2 would
require a finite `D_n`. Both would be fabrications. `NeutronThresholdEvaluation`
has no `hessian` member at all
(`static_assert(!HasHessian<NeutronThresholdEvaluation>::value)`), so the
mistake is not expressible.

**Immediate floating-point neighbours**, live-probed:

| Input | Classification | Result |
|---|---|---|
| `nextafter(n₀, 0)` | `ProtonElectron` | `PeThermodynamicEvaluation` (smooth 1-D) |
| `n₀` | `NeutronOnset` | `NeutronThresholdEvaluation` (value-only) |
| `nextafter(n₀, +∞)` | `Npe` | `EquilibriumResolutionError` (explicit) |

`EvaluatePe(n₀)` and `EvaluateNpe({n₀, n₀})` both throw. **Classification depends
on the physical threshold, not a tolerance**: `EquilibriumDomainAt` uses ordered
`<`, `==`, `>` against the stored double with no comparison window, and the
`nextafter` neighbours land on distinct branches one ULP apart.

**No neutron-density floor is used**: on the p-e branch `n_n` is *exactly* `0.0`
by construction (`MakeChargeNeutralCompositionState({n_B, n_B, 0.0})`), never a
small positive substitute; above the threshold the composition solve returns
`n_n` down to `6.79e-20` unmodified (see §11).

## 10. p-e ↔ npe response-limit derivation

Derived algebraically here before reading the test, because the task asks not to
accept it by inspection.

The p-e branch sits inside the npe chart `z = (n_B, n_e)` along `n_e = n_B`, so
its tangent is `t = (1,1)^T` (moving along the branch means `dn_e = dn_B`, which
keeps `n_n = n_B − n_e` at zero). With `A ≡ D_p + D_e`,

```text
H_npe = [[ D_n,   -D_n ],
         [-D_n, D_n + A ]].
```

**Tangential restriction.**

```text
t^T H_npe t = D_n - D_n - D_n + (D_n + A) = A = H_pe.
```

Production reproduces this to a maximum relative error of **4.44e-16** over the
four densities PE-V9 tests, and my independent 90-digit construction gives
`t^T H t / A = 1.0000000000` at every `r` from `1e-3` to `1e-20`.

**Inverse.** `det H_npe = D_n(D_n + A) − D_n² = D_n A`, hence exactly

```text
H_npe^{-1} = (1/(D_n A)) [[D_n + A,  D_n],
                          [D_n,      D_n]]
           = [[1/A + 1/D_n,  1/A],
              [1/A,          1/A]].
```

As `n_n → 0⁺`, `D_n ∝ n_n^{-1/3} → +∞`, so

```text
H_npe^{-1}  ->  (1/A) [[1,1],[1,1]]  =  (t t^T) / H_pe.
```

**There is no factor of `1/2`, and it is not a normalization accident.** The
off-diagonal and `(1,1)` entries are *exactly* `1/A` at every finite `D_n`, so no
factor can appear in the limit. The reason is that `t` maps the scalar
displacement `δn_B` to `δz = t δn_B`, and `t^T H t` is the second derivative with
respect to that same scalar `n_B`. The result is also normalization-independent:
using the unit tangent `û = t/√2` gives
`û û^T/(û^T H û) = (t t^T/2)/(A/2) = t t^T/A`, the same matrix. A spurious `1/2`
could only arise from mixing a normalized tangent with an unnormalized stiffness.

**Independent numerical confirmation** (90 digits, deviation reported as
`max|H^{-1} − t t^T/A| · A`, which equals `A/D_n` exactly):

| `r = n_B/n₀ − 1` | `n_n` | `D_n` | `t^T H t / A` | deviation |
|---|---|---|---|---|
| `1e-3` | `7.21e-12` | `6.84e+05` | `1.0000000000` | `7.229e+01` |
| `1e-6` | `6.03e-15` | `7.26e+06` | `1.0000000000` | `6.812e+00` |
| `1e-9` | `1.68e-18` | `1.11e+08` | `1.0000000000` | `4.453e-01` |
| `1e-12` | `7.74e-23` | `3.10e+09` | `1.0000000000` | `1.595e-02` |
| `1e-16` | `7.86e-29` | `3.09e+11` | `1.0000000000` | `1.603e-04` |
| `1e-20` | `7.86e-35` | `3.09e+13` | `1.0000000000` | `1.603e-06` |

**This is the correct physical collapse** from a 2-D active response to the 1-D
p-e tangent: the transverse (neutron-creating) direction becomes infinitely
stiff, so the susceptibility retains support only along `t`, exactly as
ADR-0010 V10's absent-species embedding requires.

**One honest limitation, recorded as P-3.** The convergence constant is large:
`A/D_n ≈ 6.81` at `r = 1e-6` and `3.60` at the N-1 policy boundary
`r ≈ 1.71e-7`, which is the closest a *production* response can be obtained.
The numerical collapse is therefore **not demonstrable from production values
anywhere in the numerically available npe domain**. It is established instead by
(a) the exact algebraic inverse identities, which PE-V9 does check tightly at
production values, and (b) an independent construction carried to `r = 1e-12`
(and here to `1e-20`). That is scientifically sufficient — but PE-V9's printed
`production final normalized error = 6.8117` is a number that should not be read
as a convergence measure. I reproduce that exact value independently, so it is
correct, merely far from its limit.

## 11. N-1 numerical-policy assessment

The rule, in `EvaluateNpe`: compute the **downward** local spacing
`ulp_B = n_B − nextafter(n_B, 0)` and return a smooth 2×2 response only when
`n_n ≥ 2^30 · ulp_B`; otherwise throw `EquilibriumResolutionError`.

**Is it mathematically well defined?** Yes, and the choice of the *downward*
spacing is the principled one. `n_e` is a double lying just below `n_B`, so
representable `n_e` values sit on a grid of spacing exactly `ulp_B`; therefore
`n_n = n_B − n_e` is an exact integer multiple of `ulp_B` (the subtraction itself
is exact by Sterbenz's lemma, since `n_e ∈ [n_B/2, 2n_B]`). The rule is thus a
**representational-resolution guarantee**: it requires the reconstructed `n_n` to
span at least `2^30` grid steps, i.e. to carry at least 30 bits of resolution.
That is precisely the quantity the previous review found to be silently
collapsing.

**Is the boundary exactly where it is documented?** Yes. A live bisection on the
availability flag locates the transition at `r* = 1.71068909…e-7`, where the
first available response has `n_n / ulp_B = 1.07374e9 = 2^30` to six figures. A
direct off-equilibrium sweep confirms a hard, tolerance-free cut:

| `n_n / ulp_B` | `1e12` | `2e9` | `1.1e9` | `1e9` | `1e6` | `1` |
|---|---|---|---|---|---|---|
| `n_B = 1e-6` | response | response | response | **error** | error | error |
| `n_B = 1e-3` | response | response | response | **error** | error | error |

The guard is applied to **all** npe responses, not only equilibrium ones.

**Reproduced response error as a function of distance above onset** — production
`H00` against my independent 90-digit `D_n` (which obtains `n_n` directly from
`n_n(μ_n)` and never subtracts):

| `r` above onset | `n_n/ulp_B` | rel. err `n_n` | **rel. err `H00`** | status |
|---|---|---|---|---|
| `1.71068909374e-7` (**first allowed**) | `1.07374e9` | `1.51e-6` | **`5.03e-7`** | response |
| `1.72779598e-7` | `1.09e9` | `1.24e-6` | `4.13e-7` | response |
| `1e-6` | `7.29e9` | `2.64e-7` | `8.80e-8` | response |
| `1e-5` | `8.09e10` | `1.44e-8` | `4.81e-9` | response |
| `1e-4` | `8.51e11` | `1.01e-9` | `3.37e-10` | response |
| `1e-3` | `8.71e12` | `5.72e-12` | `1.91e-12` | response |
| `1e-2` … `1e-1` | `≥8.8e13` | `≤7.8e-13` | `≤2.6e-13` | response |
| `1.71068909365e-7` (**first refused**) | `<2^30` | — | — | `EquilibriumResolutionError` |
| `1e-7`, `1e-8`, `1e-9`, `1e-10` | `<2^30` | — | — | `EquilibriumResolutionError` |
| `1e-12` | — | — | — | composition itself unrepresentable (older guard) |

This exactly reproduces the recorded behaviour: **responses available through
`1e-6` relative above onset, fail closed from `1e-7` through `1e-10`.**

**Assessment against the five criteria in the task.**

| Criterion | Verdict |
|---|---|
| Mathematically well defined | **Yes** — an exact ULP-count predicate on an exactly representable quantity, using the correct (downward) spacing |
| Conservative enough | **Yes.** The dominant error is the solver's residual noise in `n_e` (`δn_e ≈ δF/\|F'\|` with `δF ≈ 1e-13 MeV` from the `~939 MeV` cancellation), which is roughly constant near onset at `~10^3 ulp_B`. Pinning `n_n ≥ 2^30 ulp_B` therefore caps `δn_n/n_n` at `~10^{-6}` and `δD_n/D_n ≈ (1/3)δn_n/n_n` at `~5e-7`. Measured worst case over the whole allowed region — **at the boundary itself, which is where it must be** — is **`5.03e-7`**. Without the rule the error was unbounded (`6.0e-5` at `r=1e-9`, diverging as `r→0`) |
| Not excessively arbitrary | **Acceptable.** `2^30` is a round power of two, but it is exposed as the public, testable constant `MinimumResolvedNpeNeutronUlps()` rather than a buried literal, PE-V11 asserts the returned states satisfy it, and the resulting `~5e-7` cap is a sensible engineering margin (30 bits of representation minus ~10 bits of solver noise ≈ 20 bits ≈ `1e-6`) |
| Tied to actual floating-point resolution | **Yes** — `ulp_B` is recomputed per call from `nextafter`, so the rule scales automatically with `n_B` and never bites where it should not (at `n_B = 0.4`, `2^30·ulp_B = 5.96e-8 ≪ n_n = 0.396`) |
| Honest about what it guarantees | **Mostly** — see finding **P-1**. The code comment and record describe the *rule* accurately, but the record's phrase "the maximum allowed `H_00` relative error is `3.600691649285892e-7`" is the maximum over PE-V11's four probes, not over the allowed region. The true worst allowed value, measured here at the policy boundary, is `5.03e-7`. Both remain below PE-V11's asserted `1e-6`, so nothing in the code or the tests is wrong |

**No classification shift, no floor, no substitution** — verified live at every
refused density: `EquilibriumDomainAt` still returns `Npe`;
`EquilibriumStateAt` still returns a valid composition with `n_n > 0`
(`4.87e-16`, `3.26e-17`, `1.68e-18`, `6.79e-20` at `r = 1e-7 … 1e-10`),
`n_p == n_e` and `n_μ == 0`; and no p-e result is ever returned above onset.

**Disposition of N-1: CLOSED.** The unbounded near-onset `H00` error identified
by the Phase-5A-4 review is now bounded at `~5e-7` by an explicit, inspectable,
fail-closed rule that changes no physical classification.

## 12. Muon / high-density regression

| Check | Result |
|---|---|
| Muon onset | `0.45698480541241986996` — **bit-identical** to the value probed in the Phase-5A-4 review |
| `Sigma-minus` ceiling | `0.6173552079665349801` — **bit-identical** |
| Exact muon threshold | still `MuonThresholdEvaluation`, value-only; probe returns `n_μ=0`, `ε=460.15924880094979699`, `μ_n=1049.860653941171222`, `η_npe=1.42e-14`, `η_npμ=2.84e-14` — **all bit-identical to the previous review's probe** |
| Finite `D_μ` at `n_μ = 0` | `Muon().Evaluate(0).dchemical_potential_dn_MeV_fm3` is still **absent** |
| npe-μ branch equations | unchanged; `Evaluate` body untouched in the diff. `EquilibriumAt(0.5)` gives `n_μ=1.8516e-4`, `H00=148.80`, `H11=6991.32`, `H22=20455.51`, `ε=505.457` |
| Near-muon fail-closed window | unchanged (`SolveActiveEquilibriumUnchecked` guard untouched); the 15-ULP declared-vs-true band remains fully fail-closed (finding N-6, carried) |
| `Sigma-minus` source ceiling | `EquilibriumAt(n_Σ)` still rejected |

**Byte-identical regression evidence.** I built the `41ab66b` (master) tree and
diffed its ladder stdout against the reviewed branch's:

```text
diff pre_rfg.out post_rfg.out   ->  BYTE-IDENTICAL   (RFG1-RFG11)
diff pre_v.out   post_v.out     ->  BYTE-IDENTICAL   (V1-V10)
```

Every previously reported RFG metric, both onsets, and the whole Phase-5A-2
ladder are numerically untouched. `R1-V1` through `R1-V8` also print identically;
`R1-V9`/`R1-V10` differ only in the text and assertions that *had* to change
because p-e, the neutron threshold and vacuum are now implemented — and those
changes are **strengthenings** (see §14).

## 13. Active-response API assessment

The variant now carries six alternatives:

| Alternative | Active particles | `response_dimension` | `domain_status` | Hessian member |
|---|---|---:|---|---|
| `LocalThermodynamicEvaluation` | n,p,e,μ | 3 | `SmoothInterior` | 3×3 |
| `NpeThermodynamicEvaluation` | n,p,e | 2 | `SmoothInterior` | 2×2 |
| `PeThermodynamicEvaluation` | p,e | 1 | `SmoothInterior` | 1×1 |
| `MuonThresholdEvaluation` | n,p,e | 0 | `SpeciesThreshold` | **none** |
| `NeutronThresholdEvaluation` | p,e | 0 | `SpeciesThreshold` | **none** |
| `VacuumBoundaryEvaluation` | none | 0 | `VacuumBoundary` | **none** |

- **Type safety: sound.** Static assertions across the two test files cover all
  six: `tuple_size` 1/2/3 for the three Hessian ranks, no conversion
  `Pe→Npe` or `Npe→full`, and `!HasHessian<>` for all three boundary types. A
  generic `std::visit` that touches `.hessian` **fails to compile** rather than
  misbehaving; `std::get<>` throws and `std::get_if<>` returns `nullptr` on the
  wrong alternative. The conjugate types (`PeConjugates`, `NpeConjugates`,
  `NeutralConjugates`) hold arrays of different rank, so none converts to another.
- **Can a caller accidentally read a boundary result as Hessian-bearing? No.**
  The boundary types have no such member, so the mistake is inexpressible; and
  the members they do carry are named `limiting_pe_conjugates` /
  `limiting_npe_conjugates`, which reads as a one-sided limit rather than a value.
- **Active-particle metadata is meaningful, not decorative.** The header comment
  "A zero-density species at its appearance threshold is not yet active" fixes the
  convention, and it is applied consistently: `NeutronThresholdEvaluation` reports
  `{n=false,p=true,e=true,μ=false}` even though it is the neutron's appearance
  point, and `MuonThresholdEvaluation` likewise excludes the muon.
- **`response_dimension = 0` remains a slightly weak label** for "no smooth
  response returned" (it is documented in-header and made unambiguous by
  `domain_status`), but with no Hessian member on any 0-dimension type it cannot
  be misused. Carried, unchanged, from the Phase-5A-4 review.
- **`VacuumCompositionState` is a separate type** rather than a relaxation of
  `ChargeNeutralCompositionState`'s strictly-positive `n_B` invariant. That is the
  right trade: the invariant that protects every other code path is preserved, and
  the vacuum state still exposes the same duck-typed accessors
  (`BaryonDensityFm3()`, …, `Coordinates()`), so a generic visitor reads
  composition uniformly across all six alternatives.
- **Default generic-provider behaviour is unchanged** and still carries the
  latent hazard recorded as **N-3**: `ILocalThermodynamicProvider::EvaluateActive`
  defaults to `Evaluate(x)`, which statically asserts full npe-μ activity. True
  for both current implementors; a hazard only for a future third provider.
- **ADR-0010 consistency: preserved, no reopening required.** The canonical
  three-coordinate chart, its conjugates `g = (μ_n, −η_npe, −η_npμ)` and the full
  `H_x` are exactly what `Evaluate` still returns. ADR-0010 requires that a
  species boundary be "a domain boundary to report, not a reason to *silently*
  change the meaning or dimension of returned derivatives" (`:120-124`) and that
  at "absent-species boundaries … critical points, and phase transitions, report
  the domain restriction rather than fabricate a full-rank Hessian" (`:204`).
  Six distinct types with declared species, declared dimension and declared
  status is the opposite of silent. The vacuum point — which ADR-0010 does not
  contemplate, its chart C being explicitly degenerate at `n_B = 0` — is handled
  as a declared boundary rather than by extending any accepted chart to zero.
  **No accepted owner decision (Q1–Q6) is reopened.**

**Sufficiency for the next phase: yes**, and one property worth recording that
neither ADR-0010 nor the implementation record states. A TOV integration needs
`ε(n_B)` and `P(n_B)`; the provider returns `ε` but no pressure. The equilibrium
pressure is nonetheless **exactly reconstructible from the returned values on
every branch**, because at equilibrium

```text
sum_i mu_i n_i = mu_n(n_n + n_e + n_mu) = mu_n n_B   (using eta_npe = eta_npmu = 0),
so   P = h n_B - epsilon,
```

with `h = h_pe` on the p-e branch and `h = μ_n` on npe and npe-μ — both returned,
including at the two value-only thresholds. Verified independently at 80 digits
against `P = Σ μ_i n_i − ε` on all three branches: agreement `≤ 3e-71` relative.
At vacuum, `ε = 0` and `n_B = 0` give `P = 0`.

## 14. PE-V1 – PE-V13 assessment

Oracle independence was checked by reading the test's own primitives: it codes
`p_F`, `μ`, `dμ/dn` (as the reciprocal of `dn/dμ`), the energy by 16384-panel
Simpson quadrature with Kahan compensation, the pressure identity, and a
common-momentum equilibrium inversion — none of which calls a production fermion
or threshold routine.

| Test | Classification | Note |
|---|---|---|
| **PE-V1** | **Load-bearing independent falsifier** | Exact (bit-equal) p-e composition, inactive n/μ, 1-D active-set metadata at seven logarithmic densities |
| **PE-V2** | **Load-bearing independent falsifier** | Independent Simpson energy and `μ_p+μ_e`; I reproduce both at 90 digits to `≤6.3e-16` |
| **PE-V3** | **Load-bearing independent falsifier** | Independent `D_p+D_e` via the reciprocal-derivative route |
| **PE-V4** | **Load-bearing independent falsifier** | Centered finite perturbations converge at order `1.9993`, pinning `H_pe = dh_pe/dn_B` rather than re-checking the same formula |
| **PE-V5** | **Load-bearing independent falsifier** | Inactive-neutron diagnostic value, positivity and strict monotonicity; correctly uses an **absolute** tolerance for a cancellation-limited quantity |
| **PE-V6** | **Useful but partially shared oracle** — see **P-2** | Its closed form is genuinely independent of production, but `long double` is 53-bit here (R8), so it cannot exceed double precision, and it anchors to the previous review's 80-digit literal. Its own comment says so. I re-derived the value at 90 digits in §7, so the anchor is now independently re-established |
| **PE-V7** | **Load-bearing, with the same anchor caveat (P-2)** | The pressure identity, positivity below onset and `P/n_B → 0` are genuinely independent; the exact onset constant `1.8964875026317866e-9` is carried from the prior review. Re-derived here in §8 |
| **PE-V8** | **Load-bearing independent falsifier** | Value-only threshold, both `nextafter` sides, and refusal of both smooth evaluators at the threshold |
| **PE-V9** | **Load-bearing for the identity; weak for the limit** — see **P-3** | `t^T H t = H_pe` to `4.4e-16` and the *exact* inverse entries `1/A+1/D_n, 1/A, 1/A, 1/A` are the real falsifiers, and they alone settle the "no 1/2 factor" question. The `inverse_error` monotonicity check is only directional and never becomes small at production values |
| **PE-V10** | **Load-bearing, with two tautological assertions** — see **P-4** | The one-sided `h → m_p+m_e` / `H_pe` growth sweep is a genuine falsifier. The two `Near(…, 0.0)` checks re-evaluate production's own double expressions bit-for-bit, so they verify the expression but cannot detect a shared rounding-order effect (I measured the vacuum diagnostic to be `6.9e-14` relative from `m_p+m_e` being rounded before subtraction — physically immaterial for a value diagnostic) |
| **PE-V11** | **Load-bearing independent falsifier** | Independent `D_n` oracle, the ULP rule re-asserted on returned states, and — importantly — a check that composition remains available where the response is refused. Its `max_h00_relative` is measured over four probes, not the region (P-1) |
| **PE-V12** | **Load-bearing dispatch falsifier** | The complete eight-way ladder, cross-checked against the domain enum so the two cannot drift apart |
| **PE-V13** | **Useful regression check** | Muon/`Σ` values are pinned only to `2e-15`/`2e-14` absolute, looser than a bit-exact check — but RFG10's independent `Σ` solve and R1-V2's independent muon onset are tighter, and I confirmed both onsets bit-identical to master in §12 |

**Missing material falsifiers** (none blocking):

1. No test evaluates the N-1 policy **at its boundary**, where the worst allowed
   `H00` error actually occurs (`5.03e-7`, measured here). PE-V11's grid stops at
   `1e-6`, comfortably inside.
2. No falsifier for the **muon-side** absent-species susceptibility embedding —
   PE-V9 supplies the neutron-side analogue only. This is carried finding **N-4**,
   and the record says so explicitly.
3. **R2 remains uncovered**: `IntrinsicChemicalPotentialsMeV` still has zero
   callers and zero assertions repository-wide, so ADR-0010's `S_x^T μ = g`
   obligation is still untested.

**Do the pre-existing ladders remain trustworthy? Yes.** RFG1–RFG11 and V1–V10
stdout are byte-identical to master (§12). The only changes to
`rotochemical_trackr_npe.cpp` are five edits in R1-V9/R1-V10, and all five are
**strengthenings or necessary updates**, none a weakening:

- two `Throws(...)` expectations replaced by *positive* type assertions
  (`PeThermodynamicEvaluation`, `NeutronThresholdEvaluation`) now that those
  branches exist;
- a new positive `PeThermodynamicEvaluation` assertion at `nextafter(n₀, 0)`;
- `EquilibriumAt(0.0)` moved from the "must throw" list to a positive
  `VacuumBoundaryEvaluation` assertion;
- **a new `Require(pressure > 0)` at the exact neutron onset**, which is precisely
  the one-line gap the Phase-5A-4 review recorded as **N-5**.

The `Throws<EquilibriumResolutionError>` guard at `nextafter(n₀, +∞)` is retained.

## 15. Focused reruns

Both configurations were reconfigured and rebuilt at the reviewed SHA from the
canonical commands (`Debug`, C++17, AppleClang,
`Python3_EXECUTABLE=/Users/keeper/miniforge3/bin/python3`, and
`COMPACTSTAR_EOS_DATA_ROOT=/Users/keeper/Documents/CompactStar/data/compose` for
the full config).

| Run | Result |
|---|---|
| `cmake -S . -B build …` + `cmake --build build -j8` | exit 0 |
| `cmake --build build-selfcontained -j8` | exit 0 |
| `rotochemical_trackr_pe` | **PE-V1 – PE-V13 PASS**, exit 0 |
| `rotochemical_trackr_npe` | **R1-V1 – R1-V10 PASS**, exit 0 |
| `rotochemical_trackr_freegas_local` | **RFG1 – RFG11 PASS**, exit 0 |
| `rotochemical_local_thermodynamics` | **V1 – V10 PASS**, exit 0 |

Every metric in the implementation record's PE table was reproduced exactly:
PE-V2 `6.6613381477509392e-16` / `1.1368683772161603e-13 MeV`;
PE-V3 `6.6613381477509392e-16`;
PE-V4 min order `1.9993158682997674`, final `2.2060625748387963e-07`;
PE-V5 final `3.6402309433469782e-07 MeV`;
PE-V6 independent `7.356728903732833408335e-09` vs production
`7.356728903732831753973e-09`;
PE-V7 `P_onset = 1.8964875026324584e-09 MeV fm^-3`, final `P/n_B` ratio
`5.5191587679462291e-05`;
PE-V9 `4.4408920985006262e-16` and `6.8117014736069059 / 0.015520519611682909`;
PE-V11 `3.600691649285892e-07`.

**Registered inventories, authenticated live from the built trees:**

| Configuration | `ctest -N` | Record | Match |
|---|---:|---:|---|
| Self-contained | **26** | 26 | ✓ |
| Full (data root set) | **49** | 49 | ✓ |

| Complete suite | Result |
|---|---|
| Self-contained, `-j1` | **26/26 passed, 0 failed**, 89.59 s (record: 26/26, 86.46 s) |
| Full | **not rerun** (policy). The record's 49/49 in 669.86 s stands; nothing found here gives a reason to distrust it — the four focused ladders reproduce exactly, RFG and V stdout are byte-identical to master, no stellar/evolution/data path is touched, and all eight baselines are unchanged |

## 16. Baseline hashes

All eight governed artifacts under `tests/baselines/` re-hashed live;
`git status --porcelain tests/baselines/` is empty and
`git diff --stat 41ab66b..HEAD -- tests/baselines/` is empty. **All eight
unchanged**, matching the implementation record and both prior reviews:

```text
7b036942f2ae599ace3b2fc8b9a7d91f6d11b5899ab7e8a88d2bb4ee6686493b  baryon_number_dscmf1_reference.tsv
2c68e2f7e871192e00322bcc12b6d8b13eb7fdbda4a6d8d7050f26f0271ce5eb  grid_convergence_cmf_1p6_debug.tsv
7c84557742ec0ec747756d118b2837fb54095a94c10194d4ccc2d0a778f5b04f  grid_convergence_cmf_1p6_trajectory.tsv
a21d4c3f6d89322cb4cca7a073e0e142f180e529332d704b7b16a145eae741c9  hartle_I_dscmf1_debug.tsv
bd49e5a091ebcc59f7c4899422200181d4e71ecf552284840454d01aac4b8d52  hartle_monopole_dscmf1_debug.tsv
afcad5d078fe458bdb441278ce56caa2d025becc6b489d00e9baad91a101c0de  passive_cooling_cmf_1p6_debug.tsv
3d9af9129a6a4ffde9e0f8c5507a160f968a861c0cf9f3b089cceecab86b701a  tov_dscmf1_reference.tsv
5c0f4b3bdb70921f8f2a869af10edc4d8f5ae3963a9d150e11ef859d21e1c678  tov_path_equivalence_dscmf1.tsv
```

No scientific baseline is created by the reviewed commit or by this review.

## 17. Remaining findings

### New in this review — all nonblocking

| # | Finding | Class | Smallest required response |
|---|---|---|---|
| **P-1** | The record's "maximum allowed `H_00` relative error is `3.600691649285892e-7`" is the maximum over PE-V11's four probes, not over the region the N-1 rule allows. Measured at the policy boundary (`r = 1.71068909374e-7`, the first allowed density) the true worst allowed value is **`5.03e-7`**. Still inside PE-V11's asserted `1e-6`, so no code or test is wrong | **DOCUMENTATION PRECISION** | State the boundary value, or add one PE-V11 probe at `r*` so the asserted bound is measured where it binds |
| **P-2** | PE-V6 and PE-V7 hardcode the reference constants `7.3567289037328326656352e-9` and `1.8964875026317866e-9` taken from the Phase-5A-4 independent review, and their own oracles use `long double`, which is 53-bit on this platform (R8). They are therefore regression checks against a previously reviewed number plus a same-precision closed form, not fresh higher-precision falsifiers. The test comment states this honestly, and **this review re-derived both at 90 digits**, so the evidence chain is not circular | **FUTURE HARDENING** — test-provenance note | Keep the literals, but cite this review (§7, §8) as their independent origin so a future platform change cannot quietly invalidate them |
| **P-3** | The rank-one collapse `H_npe^{-1} → t t^T/H_pe` is **not** numerically demonstrable from production values: the convergence constant is `A/D_n ≈ 6.81` at `r = 1e-6` and `3.60` at the N-1 boundary, the closest production can reach. The result is fully established by the exact algebraic inverse identities (PE-V9 checks these tightly) and by an independent construction, so the science is sound — but PE-V9's printed `production final normalized error = 6.81` invites misreading as a convergence measure | **DOCUMENTATION / TEST WORDING** | Relabel that diagnostic, or print `A/D_n` with its `√r` scaling so its size is self-explanatory |
| **P-4** | PE-V10's two `Near(…, 0.0)` vacuum checks re-evaluate production's own double expressions bit-for-bit (`mp+me`, `mn-(mp+me)`), so they are tautological. Immaterial: the vacuum diagnostic carries `6.9e-14` relative error purely from rounding `m_p+m_e` before subtraction, which is irrelevant for a value diagnostic | **NO ISSUE** — test nit | Optional: compare against a differently ordered expression |
| **P-5** | **N-2 extended to the p-e chart.** `EvaluatePe` requires `EquilibriumDomainAt(n_B) == ProtonElectron`, so a physically valid *off-equilibrium* p-e-shaped state (`n_n = 0`, `n_p = n_e = n_B`) at `n_B > n₀` is refused — verified live at `n_B = 2n₀` and `0.1`. Fail-closed and declared; irrelevant to TOV structure, which needs only the equilibrium branch, but it will matter for Layer-B off-equilibrium perturbation work | **FUTURE HARDENING** (carried, widened) | When Layer B opens, replace the equilibrium-classification gate with the response domain's own conditions |
| **P-6** | The inactive-neutron diagnostic `m_n − μ_p − μ_e` is a `~939 MeV` cancellation; its **relative** accuracy degrades to `2.5e-7` at `0.999999 n₀` while its absolute accuracy stays at the `~1e-13 MeV` floor. PE-V5 correctly bounds it absolutely | **NO ISSUE** — consumer note | Record that this diagnostic is an absolute-accuracy quantity |
| **P-7** | Inside the N-1 refused window the **energy density is also unavailable**, because `ε` is only returned alongside the response. Verified live at `r = 1e-7 … 1e-10`: composition available, `ε`/response not. The window is `n_B ∈ (7.3567289037e-9, 7.3567301622e-9]`, a relative width of `1.71e-7` at `ε ≈ 6.9e-6 MeV fm^-3` (`≈1.2e7 g/cm^3`), where the p-e value of `ε` at the same `n_B` differs by at most `4.6e-13 MeV fm^-3` (`6.7e-8` relative). Negligible for structure, but the next task must not assume `EquilibriumAt` returns values at every sampled density | **FUTURE HARDENING** — flag for the whole-star task | Either have the whole-star integrator skip/step past the window, or later add a value-only npe alternative for it |

### Carried findings — dispositions authenticated

| Finding | Claimed | Independent verdict |
|---|---|---|
| **N-1** | CLOSED | **CONFIRMED CLOSED.** The `2^30`-ULP rule caps the `H00` relative error at `5.03e-7` (measured at the boundary) where it was previously unbounded; classification, composition and the absence of any floor are all preserved (§11) |
| **N-5** | CLOSED | **CONFIRMED CLOSED.** `Require(pressure > 0)` is now asserted at the exact neutron onset in *both* PE-V7 and R1-V9 — the precise one-line gap recorded by the Phase-5A-4 review |
| **N-7** | OBSOLETE | **CONFIRMED.** `TrackRFreeGasThermodynamicProvider::EquilibriumAt` is no longer a duplicate of the base body; it performs the six-way typed dispatch, which the base cannot (the vacuum result has no `ChargeNeutralCompositionState` to round-trip through) |
| **R2** | REMAINS | **CONFIRMED OPEN.** `grep` finds `IntrinsicChemicalPotentialsMeV` only in its own declaration and definition — still zero callers, zero assertions |
| **R3** | REMAINS | **CONFIRMED OPEN** as written; the RFG ladder is untouched. The PE oracle independently counts `g=2` phase space, narrowing the underlying defect further |
| **R5** | REMAINS | **CONFIRMED OPEN.** RFG7 untouched |
| **R9** | REMAINS | **CONFIRMED OPEN.** `CompactStar/EOS/src/LocalThermodynamics.cpp` is not in the diff |
| **N-2** | REMAINS | **CONFIRMED OPEN**, and now widened to the p-e chart — recorded above as **P-5** |
| **N-3** | REMAINS | **CONFIRMED OPEN.** The default `EvaluateActive` still returns `Evaluate(x)` and still statically asserts full npe-μ activity |
| **N-4** (muon side) | REMAINS | **CONFIRMED OPEN.** PE-V9 adds the neutron-side embedding falsifier only; the muon-side `χ_x → embed(χ_npe)` obligation is untested, as the record states |
| **N-6** | REMAINS | **CONFIRMED OPEN and fully mitigated.** The `~15`-ULP muon-onset construction rounding is unchanged, and the fail-closed guard still covers the whole misclassification window |
| **INV-09** | INTENDED BUT UNVERIFIED | **CONFIRMED.** No global sequence response exists |
| **INV-11** | UNRESOLVED | **CONFIRMED.** No evolved/redshifted-state convention is introduced; the p-e branch adds only a local conjugate |

**Does any remaining finding block (A) acceptance of local free-gas
thermodynamic coverage? No.** Every open item is either test hardening (R2, R3,
R5, P-2, P-4), a latent design hazard with no current trigger (N-3, R9), a
declared conservative restriction that fails closed (N-2/P-5, N-6), or a
documentation-precision note (P-1, P-3, P-6). None identifies a wrong number,
sign, unit, domain rule or convention.

**Does any block (B) beginning whole-star TOV structure reproduction? No** — with
one item to carry into that task's plan rather than fix first: **P-7**, the
`1.71e-7`-wide window above neutron onset where `ε` is unavailable. Its physical
effect on a structure integration is bounded by `6.7e-8` relative in `ε` at
`~1.2e7 g/cm^3`. N-2/P-5 does not bite either, because TOV needs only the
equilibrium branch; it becomes relevant at Layer B.

## 18. Scope integrity

`git diff --name-status 41ab66b..HEAD` — exactly 10 files:

- **Production (3):** `CompactStar/EOS/LocalThermodynamics.hpp` (purely additive
  — new `Pe*`, `NeutronThreshold`, `Vacuum*` types, a third
  `LocalResponseDomainStatus` enumerator, and three new variant alternatives; **no
  existing field, name, unit, dimension or status meaning changed**),
  `CompactStar/EOS/TrackRFreeGasThermodynamics.hpp`,
  `CompactStar/EOS/src/TrackRFreeGasThermodynamics.cpp`.
- **Test/build (3):** `tests/eos/rotochemical_trackr_pe.cpp` (new),
  `tests/eos/rotochemical_trackr_npe.cpp` (five edits, all strengthenings — §14),
  `tests/CMakeLists.txt` (one new target).
- **Documentation (4):** roadmap, invariants, current architecture, and the new
  Phase-5A-5 record. Each carries the narrow-claim caveat explicitly.

**Nothing outside `CompactStar/EOS/`, `tests/` and `docs/` is touched at all** —
verified by filtering the changed-file list. Per the task's checklist:

| Component | Files in tree | Changed |
|---|---:|---:|
| `TOVSolver` | 4 | **0** |
| `StarProfile` | 2 | **0** |
| `NStar` | 4 | **0** |
| `RotationSolver` | 2 | **0** |
| `ThermalSolver` | 0 (no such file) | **0** |
| `RotochemicalCache` | 2 | **0** |
| `ChemState` | 2 | **0** |
| BNV sources, EOS data, baselines, `src/LocalThermodynamics.cpp` | — | **0** |

Grepping the **added** lines of the production and test diff for `APR`, `BPAL`,
`CMF`, `Btilde`, paper `B`/`Z`/`W`, susceptibility, integration, spin-down,
redshift, neutrino, heat capacity, evolution, `TOVSolver`, `StarProfile`,
`RotationSolver`, `NStar`, `ThermalSolver`, `RotochemicalCache`, `ChemState`,
`BNV` returns **nothing**.

Confirmed **NOT implemented**: whole-star TOV reproduction; stellar
susceptibility integration; paper `B`/`Z`/`W`; APR/BPAL; DS(CMF) off-equilibrium
physics; rotochemical evolution; BNV. There is no second electrostatic
projection and no tolerance-based active-set reclassification.

**SCOPE INTEGRITY: CLEAN.**

## 19. Whole-star-local-readiness decision

The central question is whether the local FR2005 free-gas model is now complete
and trustworthy enough to serve as the EOS/thermodynamic source for a **separate**
whole-star structure reproduction task.

| Requirement | Verdict | Evidence |
|---|---|---|
| Vacuum handled correctly | **YES** | §6 — `ε=0`, finite `h` limit, divergent `H_pe`, value-only, negative/NaN/∞ rejected |
| p-e correct | **YES** | §5 — `h_pe`, `H_pe`, `ε_pe` match independent 90-digit formulae to `≤6.3e-16` over nine decades |
| Neutron appearance correct | **YES** | §7 — `1.24e-16` relative, condition re-derived and variationally falsified in §4 |
| npe correct | **YES** | §11–§12 — unchanged from the Phase-5A-4 review, now with a bounded near-onset accuracy policy |
| Muon appearance correct | **YES** | §12 — bit-identical, value-only, no finite `D_μ` |
| npe-μ correct | **YES** | §12 — `Evaluate` body untouched; RFG stdout byte-identical |
| `Sigma` ceiling correct | **YES** | §12 — `0.6173552079665349801` bit-identical; requests at/above it rejected |
| Every active response dimension honest | **YES** | §13 — 3/2/1/0 with no Hessian member on any 0-dimension type; enforced at compile time |
| No missing density interval | **YES** | Every `n_B ∈ [0, n_Σ)` is classified and returns either a smooth response of the correct dimension or an explicit value-only/fail-closed result |
| Numerical fail-closed regions explicit | **YES** | §11 — the N-1 rule is a public constant with a sharp, measured boundary; the near-muon 15-ULP band and the 1-ULP threshold neighbourhoods remain explicit |
| No blocking local-physics defect | **YES** | §17 — every open finding is hardening, a declared conservative restriction, or documentation precision |

**DECISION: YES.** The local Track-R free-gas thermodynamic coverage is complete
and trustworthy enough to be the source for a separate whole-star structure
reproduction task. Two supporting facts established by this review and worth
carrying into that task:

1. **Equilibrium pressure is exactly reconstructible** from the returned values on
   every branch as `P = h·n_B − ε` (§13), verified to `≤3e-71` relative — so the
   absence of a `P` accessor is not a gap for TOV.
2. **P-7**: `ε` is unavailable in a `1.71e-7`-relative window above neutron
   appearance; the structure task must step past it rather than assume totality.

This is **not** approval of whole-star reproduction results, which do not exist.

## 20. Final disposition

**PHASE 5A-5 ACCEPTABLE WITH NONBLOCKING FOLLOW-UP —
LOCAL COVERAGE COMPLETE; READY FOR HUMAN RATIFICATION AND CANONICAL INTEGRATION**
(disposition **B**).

The narrow claim in the reviewed record — that the local source model supplies
the correct physical state and an honest active response at every equilibrium
density from the `P=0` surface to the `Sigma-minus` ceiling — **is justified**.
The broader readings the record explicitly disclaims (FR2005 whole-star benchmark
reproduced, stellar coefficients, `Z`/`W`, thermal benchmark) are correctly
disclaimed and none of them is implemented.

Disposition **B** rather than **A** because of the seven nonblocking findings in
§17 and the eleven carried ones, chief among the new items being **P-1** (the
recorded worst-case N-1 accuracy is a four-probe maximum, not a regional bound;
the true boundary value is `5.03e-7`) and **P-7** (`ε` unavailable in the N-1
window). No finding blocks acceptance of local coverage, and none blocks
beginning a separate whole-star structure reproduction task.

**Review evidence versus ratification.** Everything above is *review evidence*
produced by an AI reviewer. Per `GOVERNANCE.md` §5 and `AGENTS.md` §11, the
reviewed code remains an **unverified scientific candidate**; this document does
not ratify it, and the model does not ratify its own scientific result. Human
domain ratification by the project owner is the separate, required step.

No production source or test was modified by this review. The only permanent
change is this file.

**Exactly one recommended next action:** present Phase 5A-5 and this review to
the project owner for human ratification, and on ratification fast-forward
`933494d86daf2cf8965079ece49fabd66d9390e5` plus this review to canonical
`master`. Do not begin whole-star TOV structure reproduction automatically; open
it as the next governed task after ratification, carrying **P-7** into its plan.
