# Phase 5A-4 — R1 muon-free npe branch and active-species thresholds

> **PHASE 5A-4 R1 MUON-FREE NPE BRANCH IMPLEMENTED AND VERIFIED —
> ADDITIONAL LOW-DENSITY ACTIVE-SPECIES GATE IDENTIFIED**
>
> Local analytic validation only. The source model also requires proton-electron
> matter below neutron appearance; that branch is identified but not implemented.
> Whole-star local coverage is incomplete. This code remains an **UNVERIFIED
> SCIENTIFIC CANDIDATE** pending human domain ratification (`GOVERNANCE.md:152-170`).

## 1. Authentication, scope and authority

| Item | Value |
|---|---|
| Starting SHA | `d3cbc005c53d03194407ed1080e8181568bbf1bf` |
| Canonical checkout | `/Users/keeper/Documents/CompactStar/repo/CompactStar` |
| Branch | `physics/phase5a-trackr-npe-threshold` |
| Fresh worktree | `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-rotochemical-trackr-npe` |
| Entry authentication | Canonical `master`, `origin/master`, HEAD and live remote master equal the starting SHA; canonical and fresh trees clean including untracked files |
| Branch ancestry | Every existing local branch is an ancestor of canonical master; no remote branch has a competing commit touching task files; other worktrees have no uncommitted task-file changes |
| Execution | One agent; serial validation; no suite overlap; validation completed 2026-09-05 |
| Change class | Scientific-semantic + numerical-method + structural/API + test/build registration + documentation |
| Governing contract | Accepted ADR-0010 Q2-Q4, required API semantics and threshold restrictions (`docs/adr/ADR-0010-rotochemical-off-equilibrium-thermodynamic-contract.md:120-267`, `:364-374`) |
| R1 authority | Phase-5A-3 independent review, R1 (`docs/validation/PHASE5A3_TRACKR_FREEGAS_LOCAL_INDEPENDENT_REVIEW.md:581`) |

This adds a previously unavailable analytic branch. It does not replace any
scientific baseline or invoke the pre-baseline exception in GOVERNANCE section
3.1. Existing local analytic tests passed before editing, and all eight governed
baselines were authenticated before editing and after validation. The accepted
canonical full chart and its owner decisions are unchanged: the added npe chart
is explicitly a different, lower-dimensional active response. No accepted owner decision is reopened; no existing field is redefined.

### R1 independently reproduced before editing

A standalone probe compiled the two original production translation units at the
starting SHA and obtained:

- Muon onset `0.45698480541241987 fm^-3`; Sigma-minus onset
  `0.61735520796653498 fm^-3`.
- `EquilibriumStateAt` rejected `1e-8`, `0.01`, `0.99*muon_onset`, and exact
  muon onset with `outside the strict active npe-mu source domain`.
- `Evaluate({0.3,0.003,0})` rejected zero muon density.
- `Muon().Evaluate(0).dchemical_potential_dn_MeV_fm3` was absent; the scalar
  derivative accessor threw. No finite `D_mu` was available.
- Existing V1-V10 and RFG1-RFG11 passed before changes.

The original behavior matches R1. The full 3x3 Hessian cannot be evaluated at
zero muon density by inserting zero, a small density, a large finite number,
NaN, or a guessed derivative. For a cold ideal species,
`D_mu ~ n_mu^(-1/3)` diverges there. The retained full evaluator still rejects
this boundary (`CompactStar/EOS/src/TrackRFreeGasThermodynamics.cpp:391-401`).

## 2. Local source ledger

The local library README and catalog were read before deriving the branch
inventory. The PDFs were read directly; extracted equations with ambiguous
fonts were checked against rendered pages. No memory or external EOS was used
as scientific authority.

| Exact file under `/Users/keeper/Documents/CompactStar/literature/rotochemical/` | Source locations and use |
|---|---|
| `2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf` | PDF p. 2: n,p,e,mu composition. PDF pp. 8-9 section 3.1: noninteracting Fermi gas for the whole star; no APR/BPAL crust substitution; endpoint slightly below Sigma-minus appearance. PDF p. 12 section 3.4/equation (50): each particle species contributes only in the region where it exists. |
| `2006-Reisenegger-Rotochemical-Heating-of-Neutron-Stars-Rigorous-Formalism-with-Electrostatic-Potential-Perturbations.pdf` | PDF p. 2/equations (2),(4): local charge neutrality. PDF p. 3/equations (9)-(13): intrinsic beta combinations, electrostatic cancellation and singular charged response. The neutral local Hessian already incorporates the constraint; no second projection. |
| `1995-Reisenegger-Rotochemical-Heating.pdf` | PDF pp. 14-15/section 2, equation (3): n_p=n_e and minimization of energy including electrons with respect to composition; `mu_p+mu_e-mu_n=0` for active npe matter. |

The three PDFs also match the shared library's SHA256SUMS manifest:
1995 `9af85e37c7a52fd5b704c0ba07cc0ad89741d23b049df31cb6867d501d91d0ff`,
2005 `f184d7d1d7030b61a021eb5c7ac14b1f1b30c7ea69e9d53473d153cfb069ea88`,
2006 `a286f15e083e52becd95b3000cbb5ec3ed97148681cf10a43f1a1cc5c4d23ae8`.

The 2006 source governs its correction scope. The lower threshold below is a
new derivation from those source-model assumptions, **not a number quoted from
FR2005**. The source supplies no local numerical free-gas table. The repository
constant authority remains Zaki: `m_n=939.56542052`, `m_p=938.27208816`,
`m_e=0.51099895`, `m_mu=105.6583755` MeV and
`hbar*c=197.32698045930249 MeV fm`; total energy includes rest masses
(`CompactStar/EOS/src/LocalThermodynamics.cpp:68-187`).

## 3. Complete sub-muon equilibrium active-set derivation

Minimize the sum of cold ideal-species energies at fixed positive baryon density
and zero net charge, with nonnegative species densities. Each active species
has `mu_i=sqrt(m_i^2+p_Fi^2)` and `D_i>0`. Convexity makes the constrained
minimum unique; the endpoint conditions are inequalities, not extra active
beta equalities.

**Protons are present for every n_B>0.** If protons were absent, neutrality
would also remove electrons and muons. Converting a neutron into p+e would
change energy by `m_p+m_e-mu_n(n_B)<0`, since already
`m_n>m_p+m_e`. Pure neutron matter is therefore not the minimum.

**Electrons are present whenever protons are present.** If electrons were
absent, muons would have to carry the negative charge. Replacing a muon by an
electron would change energy by `m_e-mu_mu<0`, since `mu_mu>=m_mu>m_e`.
Thus p and e remain active at all positive densities in the sub-onset inventory.

**Neutrons need not be present.** On the p-e boundary `n_n=n_mu=0` and
`n_p=n_e=n_B`. The allowed creation of a neutron by removing p+e has energy
cost `m_n-mu_p(n_B)-mu_e(n_B)`. Neutrons remain absent when that cost is
nonnegative. The sum of the two active chemical potentials increases strictly
with n_B, so it crosses m_n once.

At neutron appearance, the proton and electron Fermi momenta are equal, `p_*`,
by neutrality and identical spin degeneracy. Hence

```text
sqrt(m_p^2+p_*^2) + sqrt(m_e^2+p_*^2) = m_n,
mu_e,* = [m_n^2-m_p^2+m_e^2]/(2 m_n),
p_* = sqrt(mu_e,*^2-m_e^2),
n_B,n-onset = p_*^3/[3 pi^2 (hbar c)^3].
```

The implementation factors `m_n^2-m_p^2` as `(m_n-m_p)(m_n+m_p)`.
Numerically:

```text
mu_e,*             = 1.292581167674457 MeV,
n_B,n-onset        = 7.3567289037328318e-9 fm^-3,
P(n_B,n-onset)     = 1.89648750263e-9 MeV fm^-3 > 0.
```

A supplementary 75-digit decimal calculation, using the **displayed decimal
constants** rather than their compiled binary representations, gives
`p_*=1.1872852008365567573 MeV`,
`n_B,n-onset=7.3567289037323843851e-9 fm^-3`, and
`P_*=1.8964875026316233684e-9 MeV fm^-3`. The same independent construction
gives `n_B,mu-onset=0.45698480541241916514 fm^-3`. Differences from the
provider's declared double thresholds are rounding/constant-representation
scale; they are not new physical branches or replacement endpoint constants.

These are independently checked by a separate test-side momentum root, with
no production threshold helper (`tests/eos/rotochemical_trackr_npe.cpp:465`).
Above this boundary, a strictly interior npe equilibrium exists. Its
`d n_e/d n_B = D_n/(D_n+D_p+D_e)>0` and
`d n_n/d n_B=(D_p+D_e)/(D_n+D_p+D_e)>0`; neither active population disappears
again before muon onset.

**Muon appearance.** With active npe equilibrium, adding a muon through
`n -> p+mu` costs `mu_p+m_mu-mu_n=m_mu-mu_e`. The muon is absent for
`mu_e<m_mu`, with equality at appearance. On the p-e branch the same absence
condition follows by exchanging electron charge for muon charge; its electron
chemical potential never exceeds `1.29259 MeV`, far below `m_mu`.

| Density | Active particles | Condition/status | Implemented here? |
|---|---|---|---|
| `n_B=0` | none | Vacuum; all densities and pressure vanish | No vacuum API |
| `0<n_B<n_B,n-onset` | p,e | `mu_p+mu_e<m_n`, `mu_e<m_mu`; n and mu absent | Classification only; p-e branch is next gate |
| `n_B=n_B,n-onset` | p,e; n at appearance | `mu_p+mu_e=m_n`, `n_n=0`; no regular npe Hessian | Classification only; value/response treatment deferred with p-e gate |
| `n_B,n-onset<n_B<n_B,mu-onset` | n,p,e | `eta_npe=0`, `mu_e<m_mu`, exact `n_mu=0` | Yes: 2-coordinate active response, subject to explicit numerical-resolution failures |
| `n_B=n_B,mu-onset` | n,p,e; mu at appearance | `mu_e=m_mu`, exact `n_mu=0` | State/values only; no smooth Hessian |
| `n_B,mu-onset<n_B<n_B,Sigma-onset` | n,p,e,mu | `mu_n-mu_p=mu_e=mu_mu>m_mu` | Existing 3-coordinate response retained |

This exhausts the sub-muon active sets of the stated free-particle model.
Ideal-gas pressure is a sum of positive integrals
`P_i=[3 pi^2(hbar c)^3]^-1 integral_0^p_F p^4/sqrt(m_i^2+p^2) dp` for positive
densities. Thus neutron appearance occurs **before** the zero-pressure
surface, which is reached only as `n_B -> 0`. There is no further positive-density
p/e disappearance. Atomic binding, nuclei or a separate crust prescription
would change the source model and are not silently introduced.

## 4. Npe equations, active chart and Hessian

At exact zero muon density, local neutrality gives

```text
z=(n_B,n_e),   n_p=n_e,   n_n=n_B-n_e,   n_mu=0,
epsilon_npe(z)=epsilon_n(n_B-n_e)+epsilon_p(n_e)+epsilon_e(n_e).
```

Differentiate at fixed remaining active coordinate, before imposing beta
equilibrium:

```text
d epsilon_npe = mu_n d n_B + (-mu_n+mu_p+mu_e) d n_e,
h = (mu_n,-eta_npe),              eta_npe=mu_n-mu_p-mu_e,

H_npe = d h / d z = [[ D_n,             -D_n ],
                     [-D_n, D_n+D_p+D_e ]],       units MeV fm^3.
```

Equivalently, `S_z=[[1,-1],[0,1],[0,1]]` in species order `(n,p,e)`, and
`H_npe=S_z^T diag(D_n,D_p,D_e) S_z`. The independent test contracts this
Jacobian rather than copying the production matrix initializer. Symmetry is
exact, and `det(H_npe)=D_n(D_p+D_e)>0` with `H_00=D_n>0`. Positive definiteness
holds analytically everywhere on the declared smooth chart, beyond the tested
grid (`tests/eos/rotochemical_trackr_npe.cpp:226`, `:285`).

The optional inactive-channel diagnostic is
`eta_npmu_threshold_diagnostic_MeV=mu_n-mu_p-m_mu`. It lies outside
`NpeConjugates` and outside the 2x2 Hessian. No muon derivative is taken or
assigned. The canonical state still carries `(n_B,n_e,0)` in fm^-3; the
canonical full conjugates and full matrix fields retain their old meanings.

## 5. API representation and domains

The generic local header adds a small variant, without any stellar policy
(`CompactStar/EOS/LocalThermodynamics.hpp:109-196`, `:229-239`):

| Alternative | Active particles | `response_dimension` | Domain status | Response members |
|---|---|---:|---|---|
| `LocalThermodynamicEvaluation` | n,p,e,mu | 3 | `SmoothInterior` | Existing `NeutralConjugates`, full 3x3 `ChargeNeutralChemicalHessian` |
| `NpeThermodynamicEvaluation` | n,p,e | 2 | `SmoothInterior` | Distinct `NpeConjugates`, 2x2 `NpeChemicalHessian`, separately named inactive diagnostic |
| `MuonThresholdEvaluation` | n,p,e | 0 (no response returned) | `SpeciesThreshold` | State, energy, limiting npe conjugates and diagnostic; **no Hessian member** |

`0` at onset is an availability marker, not a zero matrix or a claim that the
physical state has no tangent directions. The response dimensions and particle
content are fixed by the selected type. There is no conversion from a 2x2 to a
3x3 response. Static assertions and generic-interface runtime tests verify that
neither the npe nor threshold alternative yields a full evaluation
(`tests/eos/rotochemical_trackr_npe.cpp:388-424`).

- Existing `Evaluate(x)` remains the full smooth npe-mu evaluator and rejects
  `n_mu=0`, even below equilibrium onset. Its off-equilibrium domain remains
  positive n,p,e,mu with n_B below Sigma-minus onset.
- `EvaluateNpe(z)` evaluates independent off-equilibrium composition; it never
  calls the equilibrium solver. Its deliberately declared domain is
  `n_B,n-onset<n_B<n_B,mu-onset`, `0<n_e<n_B`, and `mu_e(n_e)<m_mu`.
  This is a conservative response domain, not a claim that every such state is
  in beta equilibrium or satisfies the inactive weak-channel inequality away
  from equilibrium (`TrackRFreeGasThermodynamics.cpp:329-348`).
- `EvaluateActive(x)` validates the physical coordinates first, then returns the
  full alternative for strictly positive muon density, the value-only threshold
  alternative for the exact published-by-provider threshold composition, or the
  explicit npe alternative. No input density is changed (`:351-361`).
- `EquilibriumAt(n_B)` returns the active variant through the generic provider
  interface. The legacy `EquilibriumStateAt(n_B)` remains composition-only and
  is no longer an implicit promise that a full Hessian exists (`:299-314`, `:364`).
- `EquilibriumDomainAt(n_B)` reports the source-model active-set classification,
  including the unimplemented p-e and neutron-appearance gate (`:228-242`).

This is explicit active-domain treatment under ADR-0010, not an amendment to
its smooth full-chart coordinate or derivative authority. No response inverse,
electrostatic projection, redshift term or stellar embedding is added.

## 6. Solver, bracket, uniqueness and numerical failure semantics

Production solves the scalar equation

```text
F(n_e)=mu_n(n_B-n_e)-mu_p(n_e)-mu_e(n_e)=0
```

on the closed physical bracket `[0,n_B]`. At its left endpoint,
`F(0)=mu_n(n_B)-m_p-m_e>0`. At its right endpoint,
`F(n_B)=m_n-mu_p(n_B)-mu_e(n_B)<0` precisely above neutron appearance.
Inside, `F'=-D_n-D_p-D_e<0`; existence and uniqueness follow. The endpoint
values require chemical potentials only, not divergent endpoint derivatives.
Below muon onset, the root has `n_e<n_e(mu=m_mu)`; continuity and monotonicity
establish approach from below.

Bisection continues until its density endpoints are adjacent representable
numbers, with a 256-iteration limit. There is **no absolute density stopping
tolerance**. It selects the endpoint with smaller residual and requires
positive active densities, strictly sub-threshold electron potential and
`|eta_npe|<=5e-11 MeV`; no endpoint/root is clipped
(`TrackRFreeGasThermodynamics.cpp:245-296`).

Thresholds are stored as doubles computed from the same mass authority. Exact
classification uses `<`, `==`, `>` against those declared values, with no
comparison tolerance. Rounding of the constructed threshold itself is not
interpreted as an independently measurable sub-ULP species transition.

A separate conservative numerical guard checks the sign separation from the
neutron endpoint and, when applicable, from the muon endpoint before solving.
Its scale is `64*epsilon_double*(sum of positive chemical-potential scales)`,
approximately a few `1e-11 MeV`. This is a roundoff margin for subtracting
rest-mass-dominated values, not a scientific density cutoff, a rigorous
interval-arithmetic enclosure, or a replacement value for any density or
Hessian. If the side cannot be certified by that margin, the typed
`EquilibriumResolutionError` reports numerical unavailability while
`EquilibriumDomainAt` still reports the requested side. No nearby state is
returned. The full active solver receives the same narrow sign-separation
check (`:150-161`). Its actual common-potential equation, bisection and tolerances
remain unchanged elsewhere.

At **exact muon onset**, the analytical construction is

```text
n_e=n_e(mu=m_mu), n_p=n_e,
n_n=n_n(mu=mu_p(n_e)+m_mu), n_mu=0,
n_B,mu-onset=0.45698480541241987 fm^-3.
```

The state and value limits are available, but neither smooth Hessian is
returned. The public npe and full smooth methods both reject that state.
Immediately below and above, `nextafter` inputs receive distinct Npe/NpeMuon
classification and explicit resolution errors. Inputs at relative separation
`1e-10` on either side return the appropriate smooth response. Neighboring
composition coordinates at exactly the threshold baryon density are not
silently snapped to the threshold state (`tests/eos/rotochemical_trackr_npe.cpp:425`).
This is an explicit limitation of finite-precision equilibrium recovery, not
whole-star local coverage or arbitrary closeness to the threshold.

The Sigma-minus construction previously called the inner solver at one ULP
above muon onset. To avoid making provider construction depend on the now
explicitly unresolvable bracket, its outer lower endpoint uses the exact
value-only onset state. This directly addresses the near-onset part of R7;
no hyperon thermodynamics is added (`TrackRFreeGasThermodynamics.cpp:191-206`).
The Sigma-minus onset remains `0.61735520796653498 fm^-3` on this platform.

## 7. Threshold continuity and independent R1 falsifiers

The npe oracle inverts the increasing equilibrium baryon density as a function
of common proton/electron Fermi momentum. The full oracle independently
inverts total equilibrium baryon density as a function of common lepton
chemical potential. Neither calls a production equilibrium method. Direct
phase-space formulas, Simpson energy quadrature, a Jacobian contraction,
finite perturbations, and Sylvester minors supply distinct checks
(`tests/eos/rotochemical_trackr_npe.cpp:45-170`).

| Test | Defect falsified | Final result |
|---|---|---|
| R1-V1 | Wrong sub-onset equilibrium, charge imbalance, artificial muons or wrong energy | PASS at 9 densities from `1.01*n_B,n-onset` through `0.99*n_B,mu-onset`; max `|eta_npe|=1.4210854715202004e-14 MeV`, independent n_e error `1.6479873021779667e-17 fm^-3`, energy relative error `1.9428902930940239e-14`; exact n_p=n_e and n_mu=0 |
| R1-V2 | Wrong source threshold inequality or onset formula | PASS; independent muon onset `0.4569848054124197 fm^-3`, below/at/above electron-potential checks |
| R1-V3 | Wrong 2x2 entries, signs or active-coordinate reduction | PASS at three off-equilibrium states; independent `S_z^T K S_z` max relative error `6.6613381477509392e-16` |
| R1-V4 | Incorrect finite response or mutually inconsistent energy/conjugates/Hessian | PASS in both axes and a mixed direction; minimum remainder order `1.9996021829197526`; max final error `1.1783489977112984e-7 MeV`; independent energy-gradient error `3.4023093675727978e-8 MeV` |
| R1-V5 | Asymmetry or loss of positive definiteness | PASS on 121 admissible logarithmic-density/composition states from `1.01*n_B,n-onset` to `0.999*n_B,mu-onset`; max asymmetry `0`, minimum H00 `154.12025913936742`, minimum determinant `1166449.1837603361` |
| R1-V6 | Discontinuous physical values/shared active block or clipped muon derivative | PASS from both sides at relative density separations `1e-2` through `1e-9`; final scaled physical error `1.9599853578025517e-9`, shared-block relative error `1.1027947444119945e-9`; D_mu grows monotonically to `198115051.42562252 MeV fm^3` |
| R1-V7 | Invertible full response fabricated at absent muons | PASS: compile-time type separation/no threshold Hessian member; generic full Evaluate throws at n_mu=0; active variant never contains full response there; positive densities `1e-18,1e-24,1e-30` remain exact and unfloored |
| R1-V8 | One-ULP branch ambiguity or threshold snapping | PASS: exact onset is value-only; nextafter sides classify distinctly and raise typed numerical-resolution errors; relative `1e-10` sides are smooth; neighboring threshold compositions are rejected |
| R1-V9 | Npe silently extended below neutron appearance or wrong surface/domain inventory | PASS: independently derived neutron threshold, rest-mass equality and p-e inequalities; positive pressure below/on neutron appearance; distinct lower-boundary classification, explicit unimplemented p-e gate |
| R1-V10 | Regression of full chart, generic dispatch or invalid-state behavior | PASS: generic/full active-path values agree bitwise on three fixed densities; invalid coordinates and source-ceiling requests rejected; V1-V10 and RFG1-RFG11 separately pass |

R1-V4 mixed-direction remainder sequence is
`7.5385158785588047e-6, 1.8850459990287227e-6, 4.7131362399555071e-7,
1.1783489977112984e-7 MeV` under four step halvings.

R1-V6 checks n_n, n_p, n_e, total energy, mu_n, eta_npe, the inactive-channel
value diagnostic, both active beta imbalances above onset, common lepton
potential, and the shared 2x2 block against independently constructed onset
limits. It does **not** demand bounded full 3x3 conditioning: `D_mu` diverges as
the positive muon Fermi sea empties. At zero, the underlying derivative remains
unavailable, not merely large.

### Development failures retained as evidence

The first build request for the new target preceded CMake regeneration and
failed with `No rule to make target 'rotochemical_trackr_npe'`. Reconfiguring
registered the target; compilation then succeeded.

The first R1-V1 run failed its independent energy comparison at relative
`3.4989011687969196e-11` (bound `3e-12`), while beta residuals were already
`1.4210854715202004e-14 MeV`. The test-side inverse mu_n-to-n_n construction
lost relative precision near neutron appearance and changed the total baryon
density being compared. The energy oracle was corrected to use the prescribed
n_B and its independently solved n_e (`n_n=n_B-n_e`); the bound was unchanged.
The final error is `1.9428902930940239e-14`. This was a test-oracle correction,
not a relaxed tolerance or a production-physics adjustment.

## 8. Regression and build evidence

The required order is focused R1, then `rotochemical_trackr_freegas_local`, then
`rotochemical_local_thermodynamics`, then complete self-contained, then complete
authenticated full suite. Complete suites run with `-j1` and do not overlap.
The registered inventory is queried with `ctest -N`, not copied from a prior
count. Commands use CMake Debug, C++17, AppleClang and
`Python3_EXECUTABLE=/Users/keeper/miniforge3/bin/python3` as documented in
`docs/build/MACOS_BUILD.md:23-81`.

| Validation | Registered/passed | Runtime | Result |
|---|---:|---:|---|
| Focused R1 | R1-V1 through R1-V10 | direct run | PASS |
| Existing Track-R | RFG1-RFG11 | direct run | PASS |
| Existing generic local | V1-V10 | direct run | PASS |
| Complete self-contained | 25/25 | 91.49 s | PASS |
| Complete authenticated full | 48/48 | 677.45 s | PASS |

The full configuration uses
`COMPACTSTAR_EOS_DATA_ROOT=/Users/keeper/Documents/CompactStar/data/compose`.
All nine source files recorded by `docs/validation/HEAT_CAPACITY_V1.md:333-341`
match their SHA-256s; the four cold reference files in
`docs/validation/TOV_REFERENCE.md:132-135` also match (one overlapping file).
No EOS file is edited or reinterpreted for the new branch.

Before/after stdout of both original grouped ladders is byte-identical,
including every previously reported RFG metric and both onsets. RFG10's two
obsolete expectations that composition recovery reject sub-onset/onset requests
are replaced by the new typed expectations and an independent composition
comparison. Its zero-density full-Hessian rejection and all other existing
falsifiers remain. RFG9 remains a historical mutation-evidence print; no new
runtime falsifier is claimed for that line.

## 9. Eight governed baselines and scope audit

All eight governed artifacts under `tests/baselines/` remain byte-identical
to the starting Git objects. No analytic free-gas scientific baseline is created.

| Artifact | SHA-256 before and after |
|---|---|
| `baryon_number_dscmf1_reference.tsv` | `7b036942f2ae599ace3b2fc8b9a7d91f6d11b5899ab7e8a88d2bb4ee6686493b` |
| `grid_convergence_cmf_1p6_debug.tsv` | `2c68e2f7e871192e00322bcc12b6d8b13eb7fdbda4a6d8d7050f26f0271ce5eb` |
| `grid_convergence_cmf_1p6_trajectory.tsv` | `7c84557742ec0ec747756d118b2837fb54095a94c10194d4ccc2d0a778f5b04f` |
| `hartle_I_dscmf1_debug.tsv` | `a21d4c3f6d89322cb4cca7a073e0e142f180e529332d704b7b16a145eae741c9` |
| `hartle_monopole_dscmf1_debug.tsv` | `bd49e5a091ebcc59f7c4899422200181d4e71ecf552284840454d01aac4b8d52` |
| `passive_cooling_cmf_1p6_debug.tsv` | `afcad5d078fe458bdb441278ce56caa2d025becc6b489d00e9baad91a101c0de` |
| `tov_dscmf1_reference.tsv` | `3d9af9129a6a4ffde9e0f8c5507a160f968a861c0cf9f3b089cceecab86b701a` |
| `tov_path_equivalence_dscmf1.tsv` | `5c0f4b3bdb70921f8f2a869af10edc4d8f5ae3963a9d150e11ef859d21e1c678` |

A source-level before/after comparison also confirms that the entire full
`Evaluate` body, optional intrinsic-potential body, muon-onset construction and
generic fermion implementation are byte-identical. The new branch and the
explicit near-threshold numerical failure path are the scientific impact;
no existing smooth full-response equation or result is replaced. The final
scope audit finds exactly 11 changed files (3 production, 3 test/build,
5 documentation), and `git diff --check` passes.

Production scope is limited to `CompactStar/EOS/LocalThermodynamics.hpp`,
`CompactStar/EOS/TrackRFreeGasThermodynamics.hpp` and
`CompactStar/EOS/src/TrackRFreeGasThermodynamics.cpp`.
`LocalThermodynamics.cpp`, all stellar/evolution sources and all EOS data remain
unchanged. The new test is `tests/eos/rotochemical_trackr_npe.cpp`; test
registration and the two superseded RFG10 expectations are the only existing
test changes. Current-state documentation changes are the Phase 5A-3 provider
record, roadmap, invariants and architecture, plus this record.

## 10. R2-R10 disposition

The independent Phase-5A-3 review is preserved without edits. R1 has not become
a cleanup sweep; the following are deliberately not all marked resolved
(`docs/validation/PHASE5A3_TRACKR_FREEGAS_LOCAL_INDEPENDENT_REVIEW.md:582-590`).

| Finding | Disposition |
|---|---|
| R2: optional intrinsic-potential API lacks RFG tests | Unchanged; still future hardening |
| R3: independent spin-degeneracy normalization missing in RFG ladder | Existing finding retained. The new R1 oracle counts phase space independently, but the original RFG ladder was not broadened to close R3 |
| R4: RFG5 structure duplicates production | Existing no-issue/documentation note retained; RFG11 supplies structural independence |
| R5: RFG7 condition-number grid is narrower than full off-equilibrium domain | Unchanged; its finite conditioning bound remains a grid statement |
| R6: RFG9 prints historical mutation evidence | Unchanged; not represented here as a new runtime check |
| R7: one-ULP inner bracket at Sigma onset initialization | Near-muon part addressed because R1's explicit resolution guard directly requires it: use exact threshold values as the outer lower endpoint. The separate internal high-density bracketing observation remains; no general R7 cleanup or hyperon implementation |
| R8: long double has no extra precision on arm64 macOS | Note retained; numerical guard and tests do not assume extra precision |
| R9: free-lepton wrapper latent stored-member mismatch | Unchanged; wrapper implementation not touched |
| R10: missing architecture rows/candidate marking | Addressed while updating the directly affected current API description: primitive/provider rows and candidate markings added |

## 11. Local-coverage status and non-goals

The original Phase-5A-3 whole-star readiness claim is withdrawn, while its
historical implementation and measurements remain identifiable. Phase 5A-4
supplies the smooth npe branch, retains npe-mu and exposes the muon threshold
honestly. **Whole-star local coverage is still incomplete** because p-e matter
and neutron-appearance value/response semantics remain a separate local gate.
The roundoff-unresolved equilibrium neighborhood is also explicit; no claim is
made that every representable input arbitrarily close to a transition returns
a smooth numerical response.

The active subspace must later be embedded explicitly when constructing the
global response. No embedding rule, radial susceptibility integral, TOV
reproduction, spin-down particle-number response, paper Btilde/Z/W, reaction
rate or rotochemical evolution is implemented. INV-09 remains **INTENDED BUT
UNVERIFIED** and INV-11 remains **UNRESOLVED**. APR/BPAL, DS(CMF)
off-equilibrium physics, superfluidity and BNV remain outside this increment.

NO STELLAR SUSCEPTIBILITY INTEGRATION, PAPER Z/W, ROTOCHEMICAL EVOLUTION,
APR, DS(CMF) OFF-EQUILIBRIUM PHYSICS, SUPERFLUIDITY, OR BNV IS IMPLEMENTED
IN PHASE 5A-4.

**Exactly one recommended next action:** open the governed local p-e branch and
neutron-appearance-boundary implementation/validation gate. Do not begin it
as part of this task.
