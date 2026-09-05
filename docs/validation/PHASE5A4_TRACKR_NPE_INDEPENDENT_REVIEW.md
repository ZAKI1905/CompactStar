# Phase 5A-4 — Independent review of the Track-R npe branch and active-species semantics

> **PHASE 5A-4 ACCEPTABLE WITH NONBLOCKING FOLLOW-UP —
> READY FOR CANONICAL INTEGRATION; PE BRANCH IS NEXT GOVERNED GATE**
>
> R1 is genuinely resolved for the npe branch. The neutron-appearance threshold,
> the active 2x2 conjugates and Hessian, and the value-only muon-threshold
> semantics were re-derived from the local source PDFs and from first principles
> and independently reproduced numerically. No fabricated full-rank response
> exists at zero muon density. No production or test file was modified by this
> review. Per `GOVERNANCE.md` §5 the reviewed code remains an **unverified
> scientific candidate**; this review is evidence, not ratification.

## 1. Reviewed SHA

| Item | Value |
|---|---|
| Reviewed HEAD | `b1b736a4f268de95a54630e3a5c09759ac5ac812` |
| Commit subject | `feat: support muon-free Track-R thermodynamics` |
| Branch | `physics/phase5a-trackr-npe-threshold` |
| Review worktree | `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-rotochemical-trackr-npe` |
| Canonical master | `d3cbc005c53d03194407ed1080e8181568bbf1bf` |
| Merge base | `d3cbc005c53d03194407ed1080e8181568bbf1bf` |
| Change class of the reviewed commit | scientific-semantic + numerical-method + structural/API + test/build registration + documentation |
| Change class of this review | documentation |

## 2. Authentication

Performed before reading any file for scientific content.

| Check | Result |
|---|---|
| `git status --porcelain` in the review worktree | empty (clean, including untracked) |
| Local HEAD | `b1b736a4f268de95a54630e3a5c09759ac5ac812` |
| Local upstream `@{u}` | `b1b736a4f268de95a54630e3a5c09759ac5ac812` — **equal** |
| Live remote `git ls-remote origin refs/heads/physics/phase5a-trackr-npe-threshold` | `b1b736a4f268de95a54630e3a5c09759ac5ac812` — **equal** |
| Local `master` | `d3cbc005c53d03194407ed1080e8181568bbf1bf` |
| Live remote `refs/heads/master` | `d3cbc005c53d03194407ed1080e8181568bbf1bf` — **unchanged** |
| `git merge-base HEAD master` | `d3cbc005c53d03194407ed1080e8181568bbf1bf` |
| `git rev-list --left-right --count master...HEAD` | `0 1` — **exactly 1 ahead / 0 behind** |
| Diff scope vs master | 11 files: 3 production, 3 test/build, 5 documentation |

All four expectations stated in the review task hold exactly.

## 3. Direct-source ledger

The three governing PDFs were read directly from the shared read-only library and
their SHA-256s recomputed live; all three match the manifest quoted in the
implementation record (`PHASE5A4_TRACKR_NPE_BRANCH.md:66-69`).

| File under `/Users/keeper/Documents/CompactStar/literature/rotochemical/` | SHA-256 | What this review took from it |
|---|---|---|
| `2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf` | `f184d7…9ea88` | §3.1 item 3, verbatim: **"a non-interacting Fermi gas (see, e.g., Shapiro & Teukolsky 1983) for the whole star."** Items 1–2 (APR, BPAL) are explicitly *"supplemented with"* separate inner/outer crust EOSs; item 3 is not. §3.4 eq. (50): "The sum is performed over all particle species and the integral is done over the region where each particle species exists. We take only free particles into account, neglecting the heat capacity due to the lattice of ions in the crust." §3.5: the rotochemical integrals of eq. (12) are taken "only over the core". Table 1 free-gas entry is "slightly below the appearance of Σ− hyperons" |
| `2006-Reisenegger-…-Electrostatic-Potential-Perturbations.pdf` | `a286f1…3ae56` | Eq. (2): total redshifted potential contains `q_i φ`. Eq. (4): local charge neutrality `Σ n_i q_i = 0`, adopted as an excellent approximation replacing the Poisson equation. §3: imbalances `η_npl = μ_n − μ_p − μ_l`, `l ∈ {e, μ}` — intrinsic combinations from which `φ` cancels. Eq. (6): the reaction volume integral runs over the core |
| `1995-Reisenegger-Rotochemical-Heating.pdf` | `9af85e…1d0ff` | §2: `n_p = n_e` "in order to preserve charge neutrality"; the equilibrium composition "is found by minimizing the total energy per baryon `E(n,x,s)` (including the contribution from the electrons) with respect to `x`". This is the exact variational principle used in §4 below |

ADR-0010 is treated as accepted authority for the local cold charge-neutral
contract, active-species/domain honesty, no fabricated response at thresholds,
and no second electrostatic projection
(`docs/adr/ADR-0010-…-thermodynamic-contract.md:120-267`, `:364-374`).
No scientific conclusion in this review is taken from the Phase-5A-4
implementation record.

**Constant authority independently probed**, not read from the record: a scratch
program linked against `libZaki.a` printed
`m_n = 939.56542051999997511`, `m_p = 938.27208815999995295`,
`m_e = 0.51099894999999995182`, `m_μ = 105.65837550000000533`,
`m_Σ− = 1197.4490000000000691`, `hbar c = 1/MEV_2_INV_FM = 197.32698045930249009`
MeV (fm). The same probe printed `__LDBL_MANT_DIG__ = 53`, independently
re-confirming Phase-5A-3 finding **R8** (no extra `long double` precision on
arm64 macOS).

## 4. Independent active-set derivation

Derived from the source variational principle (1995 §2), **not** from the
implementation.

Minimize `ε(n_n,n_p,n_e,n_μ) = Σ_i ε_i(n_i)` at fixed `n_B > 0` subject to
`n_n + n_p = n_B` (multiplier `λ_B`), `n_p − n_e − n_μ = 0` (multiplier `λ_Q`),
and `n_i ≥ 0`. Each `ε_i` is convex (`μ_i` strictly increasing), the feasible set
is a convex polytope, so the minimum is unique. KKT gives

```text
n_n>0 : mu_n = lambda_B        n_n=0 : m_n     >= lambda_B
n_p>0 : mu_p = lambda_B+lambda_Q   n_p=0 : m_p >= lambda_B+lambda_Q
n_e>0 : mu_e = -lambda_Q       n_e=0 : m_e     >= -lambda_Q
n_mu>0: mu_mu = -lambda_Q      n_mu=0: m_mu    >= -lambda_Q
```

**Protons are present at every `n_B>0`.** If `n_p = 0`, neutrality forces
`n_e = n_μ = 0`, and the KKT conditions require `m_p ≥ λ_B + λ_Q` with
`λ_Q ≥ −m_e` and `λ_B = μ_n(n_B)`. The most permissive case gives
`m_p + m_e ≥ μ_n(n_B) ≥ m_n`, which is false because
`m_n − (m_p + m_e) = +0.78233341 MeV > 0`. Pure neutron matter is therefore
never the minimum, at any density.

**Electrons are present whenever protons are.** If `n_e = 0` and `n_μ > 0`, then
`−λ_Q = μ_μ ≥ m_μ`, while `n_e = 0` requires `m_e ≥ −λ_Q ≥ m_μ`. Contradiction,
since `m_e < m_μ`.

**Neutrons need not be present.** On the p-e boundary `n_n = n_μ = 0`,
`n_p = n_e = n_B`, and `λ_B = μ_p + μ_e`. The KKT inactivity condition is
`m_n ≥ μ_p(n_B) + μ_e(n_B)`. As `n_B → 0⁺` the right side tends to
`m_p + m_e < m_n`, so the condition holds; `μ_p + μ_e` is strictly increasing in
`n_B`, so it crosses `m_n` exactly once. Below the crossing neutrons are absent,
above it they are present — sign verified explicitly on both sides in §5.

**Muons are impossible on the p-e branch.** Their onset needs `μ_e ≥ m_μ`, but
`μ_e` on that branch is bounded above by its value at neutron appearance,
`1.29258 MeV ≪ 105.65838 MeV`.

**Independently derived equilibrium active-species sequence** for the source
free-gas species inventory `{n,p,e,μ}` (with `Σ−` as the declared ceiling):

```text
vacuum  ->  p-e  ->  n-p-e  ->  n-p-e-mu
```

This **matches** the sequence asserted by the implementation. There is **no
additional particle threshold between vacuum and neutron appearance**: protons
and electrons are active at every positive density, neutrons and muons are the
only other species in the model, and each appears exactly once.

Transition conditions, each with its inequality side:

| Boundary | Equality at the boundary | Inequality below | Inequality above |
|---|---|---|---|
| vacuum → p-e | `n_B = 0` | — | `n_p = n_e = n_B > 0` |
| p-e → npe | `m_n = μ_p(n_B) + μ_e(n_B)` at `n_p = n_e = n_B` | `μ_p+μ_e < m_n` ⇒ `n_n = 0` | `μ_p+μ_e > m_n` ⇒ `n_n > 0` |
| npe → npe-μ | `μ_e = m_μ` with `η_npe = 0` | `μ_e < m_μ` ⇒ `n_μ = 0` | `μ_e > m_μ` ⇒ `n_μ > 0` |

## 5. Neutron-onset derivation and independent numerical result

At neutron appearance `n_p = n_e` (neutrality) and both are spin-1/2 with the
same degeneracy, so their Fermi momenta are equal, `p_*`. The condition is

```text
sqrt(m_p^2+p_*^2) + sqrt(m_e^2+p_*^2) = m_n.
```

Because `A+B = m_n` and `A²−B² = m_p²−m_e²`, this closes in one step:

```text
mu_e,* = (m_n^2 - m_p^2 + m_e^2)/(2 m_n),
mu_p,* = (m_n^2 + m_p^2 - m_e^2)/(2 m_n),
p_*    = sqrt(mu_e,*^2 - m_e^2),
n_B,n-onset = p_*^3 / (3 pi^2 (hbar c)^3).
```

**The threshold condition `m_n = μ_p(n_B) + μ_e(n_B)` checked in the task is
correct.** It is the `n_n → 0⁺` limit of `η_npe = μ_n − μ_p − μ_e = 0`, since a
cold ideal neutron gas has `μ_n(0) = m_n`.

Computed at 80-decimal-digit precision in an independent Python `decimal`
script, using the exact binary values of the Zaki double constants and no
repository routine (`p_F`, `μ`, `n(μ)`, `ε`, `dμ/dn`, `P` all re-implemented):

```text
mu_e,*                        = 1.2925811676744569648807 MeV
mu_p,*                        = 938.27283935232551814814 MeV
mu_p,* + mu_e,* - m_n         = -2e-77  (exact to working precision)
p_*  (from the electron leg)  = 1.1872852008365808735371 MeV
p_*  (from the proton leg)    = 1.1872852008365808735371 MeV  (agree to 1.1e-74)
n_B,n-onset                   = 7.3567289037328326656352e-9 fm^-3
```

Direct re-evaluation gives `μ_p(n_onset) + μ_e(n_onset) − m_n = 0` to working
precision. Stability signs on both sides:

| `n_B / n_onset` | `μ_p+μ_e−m_n` (MeV) | neutron |
|---|---|---|
| 0.5 | `−2.2088e−1` | absent |
| 0.9 | `−3.7583e−2` | absent |
| 0.999 | `−3.6413e−4` | absent |
| 1.0 | `0` | onset |
| 1.001 | `+3.6391e−4` | present |
| 1.1 | `+3.5339e−2` | present |
| 2.0 | `+2.8862e−1` | present |

**Agreement with the reviewed implementation.**

| Quantity | Independent (80 digits) | Provider (probed live) | Relative agreement |
|---|---|---|---|
| `n_B,n-onset` | `7.3567289037328326656e-9` | `7.356728903732831754e-9` | **1.24e-16** (≈0.9 ULP) |
| `μ_e,*` | `1.2925811676744569649` | `1.292581167674457` (test print) | `< 1e-16` |
| `n_B,μ-onset` | `4.5698480541241904498e-1` | `4.5698480541241986996e-1` | **1.81e-15** (15 ULP) |
| `n_B,Σ−-onset` | — | `6.1735520796653498e-1` | test oracle `6.1735520796652943e-1`, 9e-15 |

The record's `7.3567289037328318e-9` is the 17-significant-digit print of the
provider's double `7.356728903732831754e-9`; my independent value agrees with it
to sub-ULP. **The neutron threshold is independently verified.**

The 15-ULP construction error in the declared *muon* onset is examined in §9.

## 6. Independent pressure check at neutron onset

The cold ideal single-species pressure is `P_i = μ_i n_i − ε_i`, with
`ε_i = m_i^4[x√(1+x²)(1+2x²) − asinh x]/(8π²(ħc)³)`, `x = p_F/m_i`, both
re-implemented independently in the 80-digit script. At neutron appearance
`n_n = 0` contributes nothing, so the total is the p-e pressure at
`n_p = n_e = n_B,n-onset`:

```text
P(n_B,n-onset) = 1.8964875026317865961253e-9 MeV fm^-3   (> 0)
  proton part  = 2.2105269382616265861e-12
  electron part= 1.8942769756935249695e-9   (dominant, as expected: p_* ~ 2.3 m_e)
```

The record's `1.89648750263e-9` and the test's Simpson value
`1.8964875026317504e-9` both agree (1.9e-14 relative, i.e. quadrature error).

Pressure immediately below onset and along the p-e branch:

| `n_B/n_onset` | `P` (MeV fm^-3) |
|---|---|
| `1e-6` | `4.0603e-19` |
| `1e-3` | `3.9852e-14` |
| `0.1` | `6.5177e-11` |
| `0.5` | `7.0584e-10` |
| `0.9` | `1.6340e-9` |
| `0.999999` | `1.8964848246e-9` |
| `1.0` (onset) | `1.8964875026e-9` |

`P` is a sum of strictly positive integrals for every positive density and is
strictly increasing; on the p-e branch it behaves as `n_B^{5/3}` (nonrelativistic
protons plus a mildly relativistic electron gas) and reaches zero **only** in the
limit `n_B → 0`.

**Conclusion: pressure at neutron appearance is positive and ~5 orders of
magnitude above nothing on this branch's own scale; the zero-pressure surface of
the source free-gas model lies strictly below neutron appearance, on the p-e
branch.** The record's claim is independently confirmed.

### Is the p-e branch genuinely part of the FR2005 source model?

Answering the task's A/B/C decomposition from the PDF text quoted in §3:

- **A — what FR2005 explicitly says.** EOS set 3 is "a non-interacting Fermi gas
  … **for the whole star**". Unlike sets 1 and 2, it is given **no** crust
  supplement. §3.4 further instructs that each species be integrated "over the
  region where each particle species exists" and that "only free particles" be
  counted, "neglecting the heat capacity due to the lattice of ions in the crust".
- **B — what follows mathematically.** Solving TOV for that EOS requires
  `P(n_B)` down to `P = 0`. §5–§6 above show `P = 0` is reached only as
  `n_B → 0`, and that for `0 < n_B < 7.3567e-9 fm^-3` the unique energy-minimizing
  charge-neutral state is p-e matter. The p-e branch is therefore not an optional
  refinement: it is the outermost part of the source model's own EOS.
- **C — what would be added realism absent from the source.** Atomic hydrogen,
  Coulomb binding, molecular formation, nuclei, a lattice, and any crust
  prescription. FR2005 explicitly declines the last of these for this EOS, and
  ADR-0010 Q1 forbids silently substituting another EOS.

**Answer: yes.** Faithful reproduction of the stated ideal free-gas EOS requires
the p-e branch down to the zero-pressure surface, and there is **no source-backed
alternative**. The only qualification worth recording for the later stellar layer
is that FR2005 §3.5 / R2006 eq. (6) integrate the *rotochemical reaction* terms
over the **core** only; the structural (TOV) integration and the eq.-(50) heat
capacity are not so restricted. That distinction narrows where the p-e response
is needed, but not whether the p-e EOS is needed.

## 7. npe solver review

The production equation is
`F(n_e; n_B) = μ_n(n_B−n_e) − μ_p(n_e) − μ_e(n_e) = 0`, derived here
independently: substituting `n_p = n_e`, `n_n = n_B − n_e` into `η_npe = 0` gives
exactly this, and it is the `∂ε/∂n_e = 0` condition of the 1995 §2 minimization.

| Property | Independent finding |
|---|---|
| Bracket | `[0, n_B]` is the exact physical range of `n_e` at fixed `n_B` |
| Left endpoint | `F(0) = μ_n(n_B) − m_p − m_e ≥ m_n − m_p − m_e = +0.78233341 MeV > 0` — **strictly positive for every `n_B>0`**, no near-cancellation |
| Right endpoint | `F(n_B) = m_n − μ_p(n_B) − μ_e(n_B) < 0` **iff** `n_B > n_B,n-onset` — the sign condition *is* the neutron-onset condition, so the bracket is valid exactly on the declared branch and nowhere else |
| Monotonicity | `F' = −D_n − D_p − D_e < 0` strictly, with `D_i = p_{Fi}²/(3μ_i n_i) > 0` |
| Uniqueness | Existence + strict monotonicity ⇒ exactly one root |
| Domain restriction | Guarded to the open interval `(n_B,n-onset, n_B,μ-onset)`; `EvaluateNpe` additionally requires `n_n>0`, `n_e>0`, `μ_e(n_e) < m_μ` |
| Density floor / active-set substitution | **None.** Bisection runs until the endpoints are adjacent representable doubles, then picks the smaller-residual endpoint; no absolute density tolerance, no clipping, no fallback to another active set. Confirmed by reading `TrackRFreeGasThermodynamics.cpp:245-296` and by R1-V7's unfloored `n_μ = 1e-18, 1e-24, 1e-30` checks |
| Near neutron onset | Fail-closed via `EquilibriumResolutionError` when `F(n_B) ≥ −roundoff_bound`; `nextafter(n_onset,+)` throws (verified live) |
| Near muon onset | Fail-closed when `F(n_e(μ_e=m_μ)) ≥ −roundoff_bound`; see §9 — this guard covers the entire declared-vs-true threshold window |

**Independent reproduction of equilibrium compositions.** I solved the same
equilibrium by a *different* parameterization — invert the strictly increasing
map `n_B(q) = n_n(μ_p(q)+μ_e(q)) + n_e(q)` in the common p/e Fermi momentum `q`,
which obtains `n_n` directly rather than by subtraction — at 80 digits, and
compared against a live probe of the provider at 15 densities:

| `n_B` (fm^-3) | rel. err `n_e` | rel. err `n_n` | rel. err `ε` | rel. err `μ_n` | rel. err `H00` | rel. err `H11` | independent `\|η_npe\|` at the returned state |
|---|---|---|---|---|---|---|---|
| `7.3567e-9` (`1+1e-7`) | 3.8e-14 | **5.8e-7** | 8.8e-17 | 5.9e-17 | **1.9e-7** | **4.9e-8** | 1.9e-14 |
| `7.4303e-9` (`1.01×`) | 7.1e-15 | 7.1e-13 | 1.7e-16 | 2.3e-17 | 2.4e-13 | 3.0e-15 | 2.6e-15 |
| `1e-8 … 4.5e-1` (13 states) | ≤ 1.8e-13 | ≤ 1.6e-14 | ≤ 6.2e-14 | ≤ 1.3e-16 | ≤ 2.2e-16 | ≤ 1.1e-13 | ≤ 1.4e-13 |

**The solver is valid.** Over the whole tested and declared branch the returned
composition, energy, conjugates and Hessian agree with an independent solve to
`≲1e-13` relative. The single degraded row is analysed as finding **N-1** in §15.

## 8. Active 2x2 thermodynamic derivation

Derived independently, before looking at the production matrix initializer.

At exact `n_μ = 0`, local neutrality gives `n_p = n_e` and `n_n = n_B − n_e`, so
the physical manifold is genuinely two-dimensional with chart `z = (n_B, n_e)`
and

```text
epsilon_npe(z) = epsilon_n(n_B - n_e) + epsilon_p(n_e) + epsilon_e(n_e).
```

Differentiating at fixed remaining coordinate, without imposing beta equilibrium:

```text
d epsilon_npe = mu_n d n_B + (-mu_n + mu_p + mu_e) d n_e
              = mu_n d n_B - eta_npe d n_e.
```

**Therefore `h = (μ_n, −η_npe)` is CONFIRMED.** It is exactly the ADR-0010
conjugate vector `g = (μ_n, −η_npe, −η_npμ)` restricted to the active
coordinates, with no re-definition of any surviving component.

```text
dh/dz : d mu_n / d n_B = D_n            d mu_n / d n_e = -D_n
        d(-eta)/d n_B  = -D_n           d(-eta)/d n_e  = D_n + D_p + D_e

H_npe = [[ D_n,        -D_n      ],
         [-D_n,  D_n+D_p+D_e ]]      MeV fm^3.
```

**CONFIRMED**, entry for entry, against
`TrackRFreeGasThermodynamics.cpp:345`.

| Property | Independent verification |
|---|---|
| Units | `D_i = ∂μ_i/∂n_i` in MeV fm³; `h` in MeV, `z` in fm^-3 ⇒ `H` in MeV fm³ ✓ |
| Signs | Off-diagonal negative because `∂n_n/∂n_e = −1`; diagonal positive ✓ |
| Symmetry | `H_npe = S_z^T diag(D_n,D_p,D_e) S_z` with `S_z = [[1,−1],[0,1],[0,1]]`, hence symmetric by construction. Live probe: `H01 − H10 = 0` **exactly** at every state tested; R1-V5 reports max asymmetry `0` over 121 states |
| Positive definiteness | `H00 = D_n > 0` and `det = D_n(D_p + D_e) > 0` **analytically everywhere** on the chart (all `D_i > 0`), not merely on a grid. Verified numerically at all probed states |
| Finite perturbations | R1-V4 remainder halves quadratically (min observed order `1.9996021829197526`), and an independent Simpson energy integral pins `h = ∂ε/∂z` to `3.4023e-8 MeV` — so a mutually consistent but wrong `(h,H)` pair is excluded |
| Relation to the full response | See below |

**Relation to the full npe-μ response as `n_μ → 0⁺`.** The full chart's Hessian,
which I re-derived independently and found to match
`TrackRFreeGasThermodynamics.cpp:414-418`, is
`H_x = S_x^T diag(D_n,D_p,D_e,D_μ) S_x`. By Cauchy–Binet
`det H_x = D_n(D_p D_e + D_p D_μ + D_e D_μ)`. Since
`D_μ ∝ n_μ^{−1/3} → ∞`, the correct statements are:

- the **upper-left 2×2 block of `H_x` is identically `H_npe`**, so the shared
  block converges trivially and continuously — live probe at `r = ±1e-10`:
  `H00 = 154.20030099352…` (above) vs `154.20030100574…` (below) vs the exact
  limit `154.2003009996…`, agreeing to `~8e-11`;
- `H_x` itself has **no finite limit** (`H22 → ∞`; probed `H22` grows as
  `r^{−1/2}`: `6.34e4, 6.28e5, 6.27e6, 6.26e7, 6.26e8` at `r = 1e-2 … 1e-10`),
  and `det H_x → +∞`. `H_x` is *unbounded*, not singular;
- the **inverse response does converge**: computing `H_x^{-1}` at 60 digits, the
  maximum deviation from the npe response `H_npe^{-1}` embedded with an
  identically zero muon mode falls as `2.41e-4, 2.41e-5, 2.41e-6, 2.41e-7,
  2.41e-9, 2.41e-12` for `r = 1e-4 … 1e-20`.

This is exactly the behaviour ADR-0010 V10 anticipates ("embed absent-species
responses explicitly"), and it is why returning the 2×2 rather than a padded 3×3
is the scientifically correct choice. The implementation's own wording — `D_mu`
"diverges as the positive muon Fermi sea empties. At zero, the underlying
derivative remains unavailable, not merely large" — is accurate; nowhere does the
record incorrectly claim the full matrix becomes singular.

## 9. Muon-threshold review

**Independent onset.** Constructing the npe equilibrium at `μ_e = m_μ`
(`n_e = n(μ=m_μ)`, `n_p = n_e`, `μ_n = μ_p(n_p) + m_μ`, `n_n = n(μ_n)`,
`n_B = n_n + n_e`) at 80 digits gives

```text
n_e            = 5.1846105465770679024e-3 fm^-3
mu_p           = 944.20227844117111311 MeV
mu_n           = 1049.8606539411711184 MeV
n_n            = 4.5180019486584197708e-1 fm^-3
n_B,mu-onset   = 4.5698480541241904498e-1 fm^-3
eta_npe        = -3e-77  (zero to working precision)
```

against the provider's `0.45698480541241986996` — **1.81e-15 relative
(15 ULP)**. The record's own supplementary 75-digit value,
`0.45698480541241916514`, sits between the two, consistent with its use of
decimal literals rather than binary doubles.

**Exact-threshold semantics — all four claims independently confirmed by live
probe of `EquilibriumAt(muon_onset)`:**

| Claim | Observed |
|---|---|
| `n_μ = 0` | `0` exactly |
| `μ_e = m_μ` | R1-V2 `Near(μ_e, m_μ, 4e-14)` passes |
| composition well-defined | `n_B=0.45698480541241986996`, `n_e=n_p=0.0051846105465770684`, `n_n=0.45180019486584277` |
| `η_npe = 0` | `1.42e-14 MeV` |
| `η_npμ = 0` as a value diagnostic | `2.84e-14 MeV`, carried in `eta_npmu_threshold_diagnostic_MeV`, outside `NpeConjugates` |
| `D_μ` is not a finite active derivative | `MuonThresholdEvaluation` has **no `hessian` member at all** (`static_assert(!HasHessian<MuonThresholdEvaluation>)`); `response_dimension = 0`; `Evaluate` throws; `Muon().Evaluate(0).dchemical_potential_dn_MeV_fm3` is `std::nullopt` |

**Is returning composition/values but no smooth Hessian scientifically correct
here? Yes.** The state is a perfectly well-defined point of the charge-neutral
manifold with a finite energy, finite `μ_n`, and finite one-sided imbalances — so
suppressing values would discard true physics. But no first derivative of the
conjugates exists there in the muon direction: `∂μ_μ/∂n_μ` diverges from above
and is undefined at `n_μ = 0`, while the two-sided `n_B`-derivative of the
equilibrium composition is itself one-sided-different (npe below, npe-μ above).
Returning a 2×2 there would also be wrong in a subtler way: it would assert a
smooth interior, when the correct classification is a boundary. `response_dimension = 0`
plus `LocalResponseDomainStatus::SpeciesThreshold` is the honest encoding, and it
is exactly what ADR-0010 requires ("report the domain restriction rather than
fabricate a full-rank Hessian", `:204`).

**One-sided limits.** Live probe at relative separations `1e-2 … 1e-10` on both
sides: below, `η_npμ_diagnostic → 0` as `−6.2e-1, −6.2e-3, −6.2e-5, −6.2e-7,
−6.2e-9` (linear in `δn_B`); above, `n_μ → 0` as `δn_B^{3/2}`; the shared 2×2
block converges from both sides to the same independently-computed limit
(§8). R1-V6 measures this two-sidedly and passes with final scaled physical error
`1.9599853578025517e-9` and shared-block relative error `1.1027947444119945e-9`.

**Robustness of the exact classification, including the 15-ULP construction
error.** Because the provider's declared constant lies **above** the true onset,
every `n_B` in the 15-ULP window `[true onset, declared onset)` is *classified*
`Npe` while physically carrying muons. I scanned that entire window one ULP at a
time:

```text
declared onset = 0.45698480541241986996
exact  onset   = 0.45698480541241904498      ULP gap = 15
window scan (16 representable densities): accepted = 0,
                                          EquilibriumResolutionError = 16,
                                          other errors = 0
```

**Not a single misclassified density returns a state.** The guard
`residual(n_e(μ_e=m_μ)) ≥ −roundoff_bound` fires throughout, because above the
true onset that residual is positive. The same argument and the R1-V9 evidence
cover the sub-ULP window at neutron appearance. The threshold classification is
therefore **fail-closed by construction**, not merely by luck: the declared
constant's accuracy bounds *which* densities are refused, never whether a wrong
active set is returned. This is the single most important robustness property in
the increment, and it holds.

## 10. Active-response API assessment

Reviewed against the six questions in the task.

| # | Question | Assessment |
|---|---|---|
| 1 | Does `std::variant` preserve type safety? | **Yes.** Three types with no implicit conversion (`static_assert(!std::is_convertible_v<NpeChemicalHessian, ChargeNeutralChemicalHessian>)`, likewise for the evaluations); `NpeConjugates` and `NeutralConjugates` hold arrays of different rank so neither converts; `MuonThresholdEvaluation` has no `hessian` member, so a generic `std::visit` that touches `.hessian` **fails to compile** rather than misbehaving; `std::get<>` throws and `std::get_if<>` returns `nullptr` on the wrong alternative. R1-V7 exercises all of these |
| 2 | Is `response_dimension = 0` semantically clear? | **Acceptable with a caveat.** It is documented in-header ("Zero means no smooth response returned, not a zero-dimensional matrix") and is *redundant with* `domain_status == SpeciesThreshold`, which already carries the meaning unambiguously. Since the type has no Hessian member, the `0` cannot be misused; it is a slightly weak label rather than a hazard |
| 3 | Can the default `ILocalThermodynamicProvider` active methods create hidden assumptions? | **Yes, mildly — recorded as N-3.** The default `EvaluateActive` returns `Evaluate(x)`, i.e. a `LocalThermodynamicEvaluation`, which now *statically* asserts `active_particles = {n,p,e,μ}`, `response_dimension = 3`, `SmoothInterior`. Any provider implementing only `Evaluate` inherits that claim silently. Today only two implementors exist (the Track-R provider, which overrides; and the Phase-5A-2 analytic toy fixture, for which the claim is true), so nothing is currently wrong — but the default is a claim, not a neutral fallback |
| 4 | Is exact threshold classification robust? | **Yes** — see §9. Exact `<`/`==`/`>` against declared doubles, no comparison tolerance, no snapping, with a fail-closed roundoff guard proven to cover the entire declared-vs-true window |
| 5 | Is `EquilibriumResolutionError` scientifically/numerically sensible? | **Yes.** It cleanly separates "this physical branch exists and `EquilibriumDomainAt` still says so" from "double precision cannot certify the bracket sign here", and it returns no substitute state or density. Its scale, `64·ε_double·(sum of positive chemical-potential scales) ≈ 2.7e-11 MeV`, is a roundoff margin for subtracting rest-mass-dominated values, correctly described in the record as neither an interval enclosure nor a physical cutoff. That it is a distinct type (not a bare `runtime_error`) is what lets tests and future callers distinguish precision failure from domain violation |
| 6 | Can the p-e branch fit as a 1-D response without changing existing meanings? | **Yes.** On that branch `n_n = n_μ = 0` and `n_p = n_e = n_B`, so the manifold is one-dimensional: `z_pe = (n_B)`, `h_pe = (μ_p + μ_e)`, `H_pe = (D_p + D_e)`, `active_particles = {false,true,true,false}`, `response_dimension = 1`, plus the two inactive-channel value diagnostics `m_n − μ_p − μ_e` and `μ_e − m_μ`. Adding a `PeThermodynamicEvaluation` alternative changes **no** field, meaning, or dimension of `LocalThermodynamicEvaluation`, `NpeThermodynamicEvaluation`, or `MuonThresholdEvaluation`, and generic visitors that read only the common static members keep compiling. A `NeutronOnsetEvaluation` value-only alternative would mirror `MuonThresholdEvaluation` exactly |

**Minimality and generality.** The addition is small (≈90 header lines, all
additive; no existing field redefined — confirmed by reading the header diff) and
scientifically explicit: every alternative names its species, dimension and
domain status. It is **not overfitted to Track R** — nothing in the generic
header mentions FR2005, free gas, or the provider; the same types would serve a
DS(CMF) npe-μ provider. It is, however, **specific to npe-μ physics**:
`ActiveParticleContent` is a fixed four-bool struct, and the type names encode
`npe`/`muon`. Under AGENTS.md §6 that would ordinarily be a hardcoding concern,
but ADR-0010 Q5 explicitly restricts thermodynamic-response v1 to a verified
npe-μ domain and requires separate governance for additional responding species,
so the specificity is governed rather than unauthorized. Adding `Σ−` later will
require a header change; that is the accepted consequence of Q5, not a defect
here.

**Does the p-e branch require reopening an accepted ADR-0010 owner decision?
No.** ADR-0010 Q2 requires an explicit three-coordinate chart with exact
four-species reconstruction and forbids *silently* changing derivative meaning or
dimension at a boundary; it simultaneously mandates that "a muon threshold or
phase boundary is a domain boundary to report" (`:120-124`) and that at
"absent-species boundaries … report the domain restriction rather than fabricate
a full-rank Hessian" (`:204`). Distinct types with declared active species,
declared dimension and declared domain status are the opposite of silent. The
canonical `x = (n_B,n_e,n_μ)` chart and the full `H_x` retain their accepted
meanings and remain what `Evaluate` returns. A 1-D p-e alternative fits under the
same accepted active-species/domain semantics. **No owner decision needs to be
reopened.** (The p-e *implementation* still needs its own governed task and
validation ladder — that is a task gate, not an ADR amendment.)

## 11. R1-V1 – R1-V10 assessment

Oracle independence was checked by reading the test's own primitives: it counts
phase space as `g·(4πp³/3)/(2πħc)³` with `g = 2` explicit, obtains `dμ/dn` as the
reciprocal of `dn/dμ`, integrates the energy by 4096-panel Simpson quadrature
with Kahan compensation, and solves equilibrium in a *different* independent
variable (common Fermi momentum, or common lepton chemical potential) from
production's `F(n_e)`. No production fermion, threshold helper or equilibrium
method appears in any oracle. I verified the Simpson normalization algebraically
(`n·Σ/panels = 3n∫₀¹u²√(m²+p²u²)du = ε`).

| Test | Classification | Note |
|---|---|---|
| **R1-V1** | **Load-bearing independent falsifier** | Independent momentum-parameterized equilibrium + independent Simpson energy + independent `μ` at 9 densities; also asserts `n_μ == 0` and `n_p == n_e` exactly |
| **R1-V2** | **Load-bearing independent falsifier** | Muon onset rebuilt from `p = √(m_μ²−m_e²)`; strict inequalities on both sides plus rest-mass equality at onset |
| **R1-V3** | **Load-bearing independent falsifier** | `S_z^T diag(D) S_z` contraction rather than a copy of the production initializer, at three deliberately off-equilibrium fixtures (`\|η\|>1e-4` enforced, so it cannot degenerate into an equilibrium-only check) |
| **R1-V4** | **Load-bearing independent falsifier** | Quadratic-remainder convergence in two axes and a mixed direction, *plus* an independent energy-gradient pin — the combination is what excludes a mutually consistent wrong `(h,H)` |
| **R1-V5** | **Useful; partially shared oracle** | Adds the 121-state grid, Sylvester minors and an independent determinant, but reuses V3's `Derivative` oracle. Positive definiteness is analytic anyway (§8), so this is breadth rather than new falsification |
| **R1-V6** | **Load-bearing independent falsifier** | The strongest test in the ladder: two-sided approach over seven decades, independent above-onset composition, shared-block convergence, monotone `D_μ` divergence, and extraction of `D_μ = H22 − H21` compared against the independent derivative |
| **R1-V7** | **Load-bearing independent falsifier** | Compile-time type separation + generic-interface dispatch + unfloored tiny `n_μ`. This is the test that actually falsifies "fabricated 3×3 at `n_μ=0`" |
| **R1-V8** | **Load-bearing independent falsifier** | One-ULP classification, typed resolution errors, no threshold snapping of neighbouring compositions, smoothness restored at `1e-10` relative |
| **R1-V9** | **Load-bearing independent falsifier — the key new result** | Neutron threshold from an independent momentum bisection on `√(m_p²+p²)+√(m_e²+p²)=m_n`; independent p-e pressure positivity; distinct lower-boundary classification and explicit unimplemented gate |
| **R1-V10** | **Load-bearing regression falsifier** | Bitwise equality of the full-chart path through old and new entry points at three densities, plus invalid-input and source-ceiling rejection |

None is redundant in the strict sense; V5 is the only one whose unique
contribution is modest.

**Missing material falsifiers** (none blocking):

1. **No falsifier of the near-neutron-onset accuracy floor.** V1's lowest density
   is `1.01×n_onset` and V5's grid starts there too, so the degradation quantified
   in §7/N-1 is entirely untested and unbounded by any assertion.
2. **No falsifier of the susceptibility limit** `H_x^{-1} → embed(H_npe^{-1})`
   with a zero muon mode. ADR-0010 V10 names this obligation; it is formally a
   later-layer duty, but the local provider is now the only place it can be
   cheaply pinned (I verified it holds — §8).
3. **V9 computes the pressure at exactly the onset and prints it, but does not
   assert it is positive** (line 491 has no `Require`); positivity is asserted
   only at `0.01, 0.1, 0.99 ×` onset. Trivial to close.
4. **R2 remains uncovered:** `IntrinsicChemicalPotentialsMeV` still has zero
   callers and zero assertions anywhere in the repository (verified by grep); the
   ADR-0010 `S_x^T μ = g` obligation is still untested.

## 12. Focused rerun results

Rebuilt both configurations at the reviewed SHA and ran the required focused
targets. Every metric the implementation record reports was reproduced
**exactly**, digit for digit.

| Run | Result |
|---|---|
| `cmake --build build -j8` (full config, data root set) | exit 0 |
| `cmake --build build-selfcontained -j8` | exit 0 |
| `rotochemical_trackr_npe` | **R1-V1 – R1-V10 PASS**, exit 0 |
| `rotochemical_trackr_freegas_local` | **RFG1 – RFG11 PASS**, exit 0 |
| `rotochemical_local_thermodynamics` | **V1 – V10 PASS**, exit 0 |

Reproduced R1 metrics: `max |η| = 1.4210854715202004e-14 MeV`;
independent density error `1.6479873021779667e-17 fm^-3`;
energy relative `1.9428902930940239e-14`;
independent muon onset `0.4569848054124197 fm^-3`;
V3 `6.6613381477509392e-16`;
V4 min order `1.9996021829197526`, max final `1.1783489977112984e-07 MeV`,
energy-gradient `3.4023093675727978e-08 MeV`, and the mixed-direction sequence
`7.5385158785588047e-6, 1.8850459990287227e-6, 4.7131362399555071e-7,
1.1783489977112984e-7`;
V5 `121` states, min `H00 = 154.12025913936742`, min det `1166449.1837603361`,
asymmetry `0`;
V6 final `1.9599853578025517e-09`, block `1.1027947444119945e-09`,
`D_μ = 198115051.42562252 MeV fm^3`;
V9 neutron onset `7.3567289037328318e-09 fm^-3`, `μ_e = 1.292581167674457 MeV`,
`P = 1.8964875026317504e-09 MeV fm^-3`.

**Independent before/after regression check.** The Phase-5A-3 record's claim that
both original ladders print byte-identical output was verified directly rather
than accepted: the master worktree at `d3cbc00` was built and its
`rotochemical_trackr_freegas_local` and `rotochemical_local_thermodynamics`
binaries run, then `diff`ed against the reviewed branch's:

```text
diff master_rfg.out branch_rfg.out  -> IDENTICAL
diff master_v.out   branch_v.out    -> IDENTICAL
```

This independently confirms that both onsets are unchanged — in particular
`Sigma-minus onset = 0.61735520796653498 fm^-3` survives the R7 bracket-endpoint
change — and that every previously reported RFG and V metric is preserved.

**Complete suites.** Per the review validation policy the ~11-minute full suite
was **not** rerun; the record's `48/48 in 677.45 s` stands, and I found no reason
to distrust it (the focused ladders reproduce bitwise, the diff touches no
stellar, evolution or data path, and all eight baselines are unchanged). The
cheap self-contained suite **was** rerun as extra evidence:

| Validation | Registered | Result |
|---|---|---|
| Self-contained inventory (`ctest -N`) | **25** | matches the record's 25 |
| Full inventory (`ctest -N`, `COMPACTSTAR_EOS_DATA_ROOT=/Users/keeper/Documents/CompactStar/data/compose`) | **48** | matches the record's 48 |
| Complete self-contained suite (`-j1`) | **25/25 passed, 0 failed**, 90.01 s | reproduces the record's 25/25 (91.49 s) |
| Complete full suite | not rerun (policy) | record's 48/48 in 677.45 s stands |

## 13. Baseline hashes

All eight governed artifacts under `tests/baselines/` re-hashed live and
`git status --porcelain tests/baselines/` is empty. **All eight are unchanged**
and match both the Phase-5A-4 record and the Phase-5A-3 review table:

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

No new scientific baseline is created by the reviewed commit or by this review.

## 14. Scope integrity

`git diff --name-status d3cbc00..HEAD` — exactly 11 files:

- **Production (3):** `CompactStar/EOS/LocalThermodynamics.hpp` (purely additive:
  new types, a `<variant>` include, and two virtual methods with defaults — no
  existing field, name, unit or meaning changed),
  `CompactStar/EOS/TrackRFreeGasThermodynamics.hpp`,
  `CompactStar/EOS/src/TrackRFreeGasThermodynamics.cpp`.
- **Test/build (3):** `tests/CMakeLists.txt` (one new target),
  `tests/eos/rotochemical_trackr_npe.cpp` (new),
  `tests/eos/rotochemical_trackr_freegas_local.cpp` (two superseded RFG10
  expectations replaced — and *strengthened*, from "must throw" to a positive
  comparison against an independent composition solve; every other RFG falsifier
  retained).
- **Documentation (5):** roadmap, invariants, current architecture, the
  Phase-5A-3 provider record's status header, and the new Phase-5A-4 record.

**No behaviour change** to `TOVSolver`, `StarProfile`, `NStar`, `RotationSolver`,
thermal evolution, neutrino cooling, legacy `Rotochemical`/`ChemState`, or BNV —
none of those files appears in the diff, and the byte-identical master-vs-branch
ladder outputs in §12 confirm no numerical drift in the shared EOS primitives.
`CompactStar/EOS/src/LocalThermodynamics.cpp` is untouched (so R9's wrapper
hazard is genuinely unchanged, as claimed).

Grepping the added (`+`) lines of the production and test diff for
`APR`, `BPAL`, `DS(CMF)`, `Btilde`, paper `B`, `Z`/`W`, susceptibility,
spin-down, redshift, evolution, reaction rate, neutrino, heat capacity, BNV,
`TOVSolver`, `StarProfile`, `RotationSolver`, `NStar`, `ChemState` returns
**nothing**. The full smooth `Evaluate` body, the optional
`IntrinsicChemicalPotentialsMeV` body, the muon-onset construction and the whole
of the generic fermion implementation are unchanged.

**SCOPE INTEGRITY: CLEAN.**

## 15. Remaining material findings

Ranked. **None is blocking.**

| # | Finding | Class | Smallest required response |
|---|---|---|---|
| **N-1** | **Near-neutron-onset relative accuracy of `n_n` and `H00` is unstated and untested.** The npe root is bisected in `n_e` over `[0,n_B]` and accepted on `\|η_npe\| ≤ 5e-11 MeV` alone, so the derived `n_n = n_B − n_e` inherits cancellation error `~ε·n_B/n_n`, and `D_n = H00 ∝ n_n^{−1/3}` inherits about a third of it. Measured against an independent 80-digit solve: at `(n_B/n_onset − 1) = 1e-3, 1e-6, 1e-7, 1e-9` the `H00` relative error is `1.9e-12, 8.8e-8, 1.9e-7, 6.0e-5`, while `n_e` stays accurate to `~1e-13` throughout. The returned state is always a valid charge-neutral state and no record claim is falsified — V1/V5 begin at `1.01×` onset, where the error is `2.4e-13` | **FUTURE HARDENING** — accuracy limitation, currently outside every declared and tested region | Either state the `n_n`-relative accuracy floor in the provider's declared domain, or solve in the common p/e Fermi momentum (as the test oracle already does) and return `n_n` from `n_n(μ_n)` directly instead of by subtraction |
| **N-2** | **`EvaluateNpe`'s off-equilibrium response domain is tied to the *equilibrium* `n_B` classification.** It requires `EquilibriumDomainAt(n_B) == Npe`, so a physically valid off-equilibrium npe state (`n_n>0`, `n_e>0`, `μ_e<m_μ`) is refused whenever `n_B` lies below neutron onset or at/above muon onset — even though the 2×2 response is well defined there. The record declares this as deliberately conservative, and it fails *closed*, never returning a wrong response | **FUTURE HARDENING** (robustness) — will matter when a stellar layer perturbs composition at fixed `n_B` near either boundary | Replace the equilibrium-classification gate with the response domain's own conditions (`n_n>0`, `n_e>0`, `μ_e(n_e)<m_μ`, `n_B <` Σ− ceiling) when the p-e gate is opened |
| **N-3** | **The default `ILocalThermodynamicProvider::EvaluateActive` silently asserts full npe-μ activity.** It returns `Evaluate(x)`, whose type now statically claims `{n,p,e,μ}` active, `response_dimension = 3`, `SmoothInterior`. A future provider that implements only `Evaluate` and has an absent-species boundary would mislabel its results without any diagnostic. True today for both existing implementors | **FUTURE HARDENING** — latent design hazard | When a second real provider appears, consider making `EvaluateActive` pure virtual (or defaulting it to throw), so each provider declares its own active-set semantics |
| **N-4** | **No falsifier for the absent-species susceptibility embedding.** ADR-0010 V10 requires that absent-species responses be embedded explicitly and verified before inversion. This review verified it holds (`H_x^{-1} → embed(H_npe^{-1})` with a zero muon mode, converging as `√δn_B`), but no repository test pins it | **FUTURE HARDENING** — formally a Layer-B obligation, cheapest to pin here | One local assertion in the R1 ladder, or explicit deferral in the Layer-B task |
| **N-5** | R1-V9 computes and prints the pressure at exactly the neutron onset but asserts positivity only at `0.01/0.1/0.99 ×` onset (`rotochemical_trackr_npe.cpp:491`) | **FUTURE HARDENING** — one-line gap | Add the `Require(pressure > 0)` at the onset itself |
| **N-6** | The declared muon-onset constant carries `1.81e-15` (15 ULP) of construction error from double round-trips through `μ ↔ n`, and is printed to 17 digits in the provider metadata, where it could be over-read as exact. **Mitigated completely**: every one of the 16 representable densities in the misclassification window fails closed (§9) | **NO ISSUE for correctness** — documentation nit | Note in the metadata or record that the declared thresholds are `O(10 ULP)` constructions and that classification near them is fail-closed |
| **N-7** | `TrackRFreeGasThermodynamicProvider::EquilibriumAt` overrides the base method with an identical body | **NO ISSUE** — cosmetic | Optional: drop the redundant override |

**No finding is a blocking scientific or active-species-contract defect.** No
production formula, sign, unit, index order, domain rule, conjugate, Hessian
entry, threshold condition or convention was found to be wrong.

### Phase-5A-3 R2–R10 dispositions authenticated

Each stated disposition was checked against the code and the Phase-5A-3 review
text, not accepted on assertion.

| Finding | Record's disposition | Independent verdict |
|---|---|---|
| **R1** | Resolved for npe; p-e gate revealed | **CONFIRMED.** R1's three required responses — (a) expose the npe equilibrium branch and its thermodynamics, (b) supply per-result active-species/threshold status, (c) without fabricating a full-rank 3×3 at `n_μ=0` — are all delivered and independently verified, and the status wording was qualified in the same change (roadmap, INV-09, INV-11, architecture, and the 5A-3 record header). R1's own framing was whole-star readiness, and that is **not** restored: §4–§6 show the p-e branch is required to reach the surface |
| **R2** | Unchanged; still future hardening | **CONFIRMED OPEN.** `grep` finds `IntrinsicChemicalPotentialsMeV` only in its own declaration and definition — still zero callers, zero assertions, `S_x^T μ = g` still untested |
| **R3** | Existing finding retained; new R1 oracle counts phase space independently, RFG ladder not broadened | **CONFIRMED — and conservatively stated.** The RFG ladder is indeed untouched. The new test's `DensityFromMomentum` does count `g=2` phase space from first principles, and via `IndependentFull` it now also pins npe-μ compositions above onset — so the underlying defect is materially narrowed even though R3 as written is not closed. The record understates the coverage gained rather than overstating it |
| **R4** | No-issue/documentation note retained | **CONFIRMED OPEN** as a note; RFG11 remains the structural falsifier |
| **R5** | Unchanged; conditioning bound remains a grid statement | **CONFIRMED OPEN.** RFG7 untouched |
| **R6** | Unchanged; not represented as a new runtime check | **CONFIRMED OPEN.** RFG9 is still an unconditional print, and the Phase-5A-4 record explicitly says so |
| **R7** | **Partially** addressed: near-muon part fixed, high-density bracketing observation remains | **CONFIRMED, and the fix was necessary rather than optional.** The Σ− outer bracket's lower endpoint moved from `nextafter(muon_onset,+∞)` to the exact onset with a special-cased value-only state; without it, the new resolution guard would have made the *constructor* throw deterministically. The second half of R7 — the doubling loop still evaluating the model up to `~0.914 fm^-3`, past the source ceiling — is untouched. The Σ− onset value is unchanged, independently confirmed by the byte-identical master-vs-branch RFG output in §12 |
| **R8** | Note retained; guard and tests assume no extra precision | **CONFIRMED.** Independently re-measured `__LDBL_MANT_DIG__ = 53` |
| **R9** | Unchanged; wrapper not touched | **CONFIRMED OPEN.** `LocalThermodynamics.cpp` is not in the diff |
| **R10** | Addressed for the directly affected API description | **CONFIRMED DONE.** Three rows added to `CURRENT_ARCHITECTURE.md` §2 (`ActiveLocalThermodynamicEvaluation` + npe value types, `ColdRelativisticIdealFermion`, `TrackRFreeGasThermodynamicProvider`), each carrying an explicit `SCIENTIFIC CANDIDATE` marking per `GOVERNANCE.md` §5 |

R2/R3/R5/R9 (and R4/R6 as notes) **remain genuinely open**, exactly as claimed.
The increment did not become a cleanup sweep.

## 16. Exact whole-star-local-coverage status

Coverage of the source model's equilibrium density axis, as implemented at
`b1b736a` (all densities in fm^-3):

| Interval | Active set | Response available | Status |
|---|---|---|---|
| `n_B = 0` | vacuum | — | No API; `n_B > 0` required |
| `0 < n_B < 7.3567289037328318e-9` | p, e | **none** | **NOT IMPLEMENTED.** `EquilibriumDomainAt` classifies `ProtonElectron`; `EquilibriumStateAt`/`EquilibriumAt` throw with an explicit "separate unimplemented gate" message |
| `n_B = 7.3567289037328318e-9` | p, e (n at appearance) | **none** | **NOT IMPLEMENTED.** Classified `NeutronOnset`; throws |
| `nextafter(n-onset, +)` (≈1 ULP) | n, p, e | none | Fail-closed `EquilibriumResolutionError` |
| `7.3567e-9 < n_B < 0.45698480541241987` | n, p, e | **2-D**: `z=(n_B,n_e)`, `h=(μ_n,−η_npe)`, `H_npe` 2×2, plus `η_npμ` value diagnostic | **IMPLEMENTED AND VALIDATED** |
| 15-ULP band below the declared muon onset | n, p, e, μ (physically) | none | Fail-closed `EquilibriumResolutionError` — no wrong active set returned (§9) |
| `n_B = 0.45698480541241987` | n, p, e (μ at appearance) | **0-D**: state, energy, limiting `h`, both value diagnostics; **no Hessian member** | **IMPLEMENTED AND VALIDATED** as a value-only boundary |
| `nextafter(μ-onset, +)` (≈1 ULP) | n, p, e, μ | none | Fail-closed `EquilibriumResolutionError` |
| `0.45698480541241987 < n_B < 0.61735520796653498` | n, p, e, μ | **3-D**: canonical `x=(n_B,n_e,n_μ)`, `g`, `H_x` 3×3 | **UNCHANGED from Phase 5A-3**, bitwise |
| `n_B ≥ 0.61735520796653498` | Σ− would appear | — | Refused; declared source ceiling |

**WHOLE-STAR LOCAL COVERAGE IS INCOMPLETE.** The gap is the p-e branch and the
neutron-appearance boundary, i.e. everything from the stellar surface inward to
`n_B = 7.357e-9 fm^-3`. §6 establishes independently that this gap is real: the
pressure at neutron appearance is `1.896e-9 MeV fm^-3 > 0`, and `P → 0` only as
`n_B → 0`, so a free-gas star cannot be integrated to its surface without it.
The record's withdrawal of the Phase-5A-3 whole-star readiness claim is correct
and necessary.

Also still outside this increment, as claimed: no embedding rule for the active
subspace in a global response, no radial susceptibility integral, no TOV
reproduction, no spin-down particle-number response, no paper `Btilde/Z/W`, no
reaction rate, no rotochemical evolution, no APR/BPAL, no DS(CMF) off-equilibrium
physics, no superfluidity, no BNV. INV-09 remains **INTENDED BUT UNVERIFIED**;
INV-11 remains **UNRESOLVED**.

## 17. Final disposition

**PHASE 5A-4 ACCEPTED WITH NONBLOCKING FOLLOW-UP —
READY FOR CANONICAL INTEGRATION; PE BRANCH IS NEXT GOVERNED GATE**
(disposition **B**).

The A/B preconditions are each independently satisfied:

- **R1 genuinely resolved** for the npe branch — §7, §8, §11, §15;
- **neutron threshold independently verified** — `7.3567289037328326656e-9` vs
  the provider's `7.356728903732831754e-9`, `1.24e-16` relative (§5);
- **active 2×2 thermodynamics correct** — `h` and `H_npe` re-derived from the
  restricted energy and confirmed entry for entry, exactly symmetric,
  analytically positive definite, `det = D_n(D_p+D_e)` (§8);
- **no fake 3×3 threshold response** — no Hessian member exists on the threshold
  type, `Evaluate` throws at `n_μ=0`, `D_μ` is genuinely unavailable and provably
  divergent, and the *inverse* response embeds with a zero muon mode (§8, §9);
- **active-response API scientifically defensible** — type-safe, minimal,
  additive, generic across providers, and compatible with a future 1-D p-e
  alternative without reopening ADR-0010 (§10);
- **no blocking regression** — both original ladders byte-identical to master,
  25/25 self-contained, 25/48 inventories, all eight baselines unchanged (§12–§14).

Disposition B rather than A solely because of the seven nonblocking findings in
§15, of which **N-1** (an unstated `n_n`/`H00` accuracy floor within `~1e-3`
relative of neutron appearance) is the only one carrying a measurable numerical
consequence, and it lies outside every currently declared and tested region.

No production source or test was modified by this review. The only permanent
change is this file. Per `GOVERNANCE.md` §5 the reviewed code remains an
**unverified scientific candidate** until a human with domain authority ratifies
it; this independent review is evidence, not ratification.

**Exactly one recommended next action:** fast-forward the Phase-5A-4 npe branch
and this independent review to canonical `master`, then open the governed local
**p-e branch and neutron-appearance-boundary** implementation/validation gate as
the next task. Do not begin that gate automatically, and do not begin whole-star
reproduction before it closes.
