# Phase 5A-3 — Track-R cold free-gas local provider

> **Status (2026-09-04): PHASE 5A-3 TRACK-R FREE-GAS LOCAL PROVIDER IMPLEMENTED
> AND VERIFIED — WHOLE-STAR FREE-GAS REPRODUCTION READY.**
>
> **TRACK-R FREE-GAS LOCAL THERMODYNAMICS IMPLEMENTED;**
> **WHOLE-STAR FREE-GAS REPRODUCTION NOT YET PERFORMED.**

This record covers only the local cold noninteracting `npe-mu` thermodynamic
model used beneath the Fernandez--Reisenegger Track-R free-gas benchmark. It
does not claim a TOV star, Table 1 reproduction, rotation sequence, stellar
susceptibility integral, particle-number response, paper `Z/W`, rate, thermal
evolution, superfluidity, or BNV result.

## 1. Authentication, authority, and change class

| Item | Authenticated value |
|---|---|
| Starting SHA | `8658fee45a77c483c2fbab6b250b9be6d5b4acdf` |
| Branch | `physics/phase5a-trackr-freegas-local` |
| Worktree | `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-rotochemical-trackr-freegas` |
| Starting canonical state | clean `master`; local, `origin/master`, and live remote equal to starting SHA |
| Governing contract | accepted ADR-0010 (`docs/adr/ADR-0010-rotochemical-off-equilibrium-thermodynamic-contract.md:178-251`) |
| Change class | scientific-semantic + numerical-method + structural/API + dependency/build + engineering/test + documentation; strictest evidence applied |

This is new source-model functionality, not a correction of a pre-existing
production result, so `GOVERNANCE.md` section 3.1 is not invoked. No historical
free-gas local output exists to preserve, and the task explicitly requires an
analytic/source-model validation rather than a scientific baseline. The new
code remains an unverified scientific candidate pending human domain review as
required by `GOVERNANCE.md:152-170`.

## 2. Source ledger read before implementation

The shared source catalog identifies the 2005 paper as the primary
non-superfluid reproduction benchmark, the 2006 paper as the correcting
thermodynamic/charge-neutrality authority, and the 1995 paper as supporting
conceptual authority (`/Users/keeper/Documents/CompactStar/literature/catalog.tsv`).
The following PDF locations were read directly:

| Authority | Exact location used | Local conclusion |
|---|---|---|
| Fernandez & Reisenegger 2005, arXiv preprint | PDF p. 2, introduction | The chosen simplest core composition is neutrons, protons, electrons, and muons (`npe-mu`). |
| Same | PDF pp. 8-9, section 3.1 | The third benchmark EOS is a noninteracting Fermi gas for the whole star. |
| Same | PDF p. 9, section 3.1 discussion | The free-gas Table 1 configuration stops slightly below `Sigma-minus` appearance; the adopted Kepler period is whole-star input, not local thermodynamics. |
| Same | PDF p. 12, section 3.4, equation (50) discussion | For the noninteracting Fermi gas, every particle species uses `m_i*=mu_i/c^2`, and the species sum covers the region where each exists. |
| Same | PDF p. 16, section 4.3, equations (71)-(74) | The first particles expected after muons are `Sigma-minus` and `Lambda0`; `eta_nnSigma p=2mu_n-mu_Sigma-mu_p`, with `n+n <-> Sigma-minus+p`. Hyperons were not included. |
| Same | PDF p. 26, Figure 1 | The free-gas curve shown is a whole-star spin-down-compression result; it supplies no additional local EOS convention. |
| Same | PDF p. 30, Table 1 and footnotes | The listed free-gas values are whole-star maximum-configuration quantities: `M_max=0.62 M_sun`, central mass-energy density `1.10e15 g cm^-3`, radii, and period; footnote c says the mass is before `Sigma-minus` appearance. |
| Same | PDF p. 33, Figure 8 | The free-gas curve shown is a whole-star quasi-equilibrium surface-temperature result; it belongs to later reproduction/evolution work. |
| Reisenegger et al. 2006 | PDF p. 2 / journal p. 569, equations (2), (4), section 3 | Local charge neutrality is imposed on intrinsic densities; `eta_npl=mu_n-mu_p-mu_l`; the electrostatic potential is not part of the imbalance. |
| Same | PDF p. 3 / journal p. 570, equations (10)-(18) | Corrected constrained local response eliminates charge, yields a reduced three-variable problem, and does not license a second projection inside the local provider. |
| Reisenegger 1995 | PDF p. 15 / preprint p. 5, equation (3) | At fixed baryon density, beta equilibrium minimizes total energy including electrons and gives `mu_p+mu_e-mu_n=0` for npe matter. |

### Source assumptions established

- Particle content: neutron, proton, electron, muon.
- All four components in the free-gas benchmark are cold noninteracting
  degenerate spin-1/2 fermions. This follows jointly from the 2005 `npe-mu`
  declaration, its whole-star noninteracting-Fermi-gas choice, and section
  3.4's all-species free-particle prescription.
- Local charge neutrality: `n_p=n_e+n_mu`.
- Beta equilibrium: `mu_n-mu_p=mu_e` on the muon-free branch and
  `mu_n-mu_p=mu_e=mu_mu` on the active branch.
- The source gives no alternate numerical particle masses or alternate
  rest-energy subtraction. Per the task authority, the implementation uses
  Zaki's `NEUTRON_M_MEV`, `PROTON_M_MEV`, `ELECTRON_M_MEV`, `MUON_M_MEV`,
  `SIGMA_MINUS_M_MEV`, and `MEV_2_INV_FM` constants. Total energy includes
  every species rest mass.
- The `npe-mu` source model ends before `Sigma-minus` appearance. No hyperon
  energy, susceptibility, rate, or composition is implemented.
- The source gives a whole-star central mass-energy density just below onset,
  not a local baryon-density ceiling. This implementation therefore derives
  the local equilibrium onset from equation (71) using the same authoritative
  ideal-particle masses; it does not silently reinterpret Table 1's
  `g cm^-3` value as `n_B`.
- Stellar structure, rotation, heat capacity, emissivity, and Table 1 data are
  not inputs needed by this local provider.

No source-critical local convention remains unresolved under the task's stated
constant authority.

## 3. Provider architecture and provenance

`CompactStar/EOS/TrackRFreeGasThermodynamics.hpp:34-92` defines
`TrackRFreeGasThermodynamicProvider`, a concrete
`ILocalThermodynamicProvider`. It is implemented in
`CompactStar/EOS/src/TrackRFreeGasThermodynamics.cpp` and registered only as an
EOS library source. There is no connection to `TOVSolver`, `NStar`,
`StarProfile`, `RotationSolver`, `RotochemicalCache`, or `ChemState`.

`CompactStar/EOS/LocalThermodynamics.hpp` adds the source-independent
`ColdRelativisticIdealFermion` primitive with neutron, proton, electron, and
muon factories. `ColdRelativisticFreeLepton` retains its public interface and
delegates to that exact implementation
(`CompactStar/EOS/src/LocalThermodynamics.cpp`). RFG1 verifies bit-identical
electron/muon results, including the zero-density optional derivative.

Concrete metadata records:

- Track-R / Fernandez--Reisenegger 2005 purpose and R2006 corrected response
  interpretation;
- local implementation revision `local-r1`;
- exact `n,p,e,mu` ideal particle content;
- canonical coordinates and neutrality identities;
- cold-only scope;
- total-energy/rest-mass convention;
- provider ownership of analytic leptons exactly once;
- numeric smooth-domain and strict active-equilibrium interval, with no
  threshold smoothing.

The typed methods `MuonOnsetBaryonDensityFm3()` and
`SigmaMinusOnsetBaryonDensityFm3()` carry the two numerical boundaries outside
free-form metadata. This addresses the prior review's provenance/domain warning
for this concrete provider without redesigning generic ADR-0010 metadata.

## 4. Ideal-species formulas and rest-mass convention

For each species `i` and `n_i>0`, with `M_i=m_i c^2` in MeV and
`hbar_c=1/MEV_2_INV_FM` in MeV fm,

```text
p_Fi = hbar_c (3 pi^2 n_i)^(1/3),
mu_i = sqrt(M_i^2+p_Fi^2),
D_i  = dmu_i/dn_i = p_Fi^2/(3 mu_i n_i),

epsilon_i = 1/[pi^2 hbar_c^3]
            integral_0^p_Fi k^2 sqrt(M_i^2+k^2) dk.
```

The authoritative masses used are `M_n=939.56542052 MeV`,
`M_p=938.27208816 MeV`, `M_e=0.51099895 MeV`, and
`M_mu=105.6583755 MeV`. `epsilon_i` includes `M_i n_i`. The shared long-double implementation retains
the previously validated small-`x=p_F/M` series to avoid cancellation but does
not clip or floor density. At zero density, `p_F=0`, `mu=M`, `epsilon=0`, and
`D_i` is unavailable. The neutron and proton implementations use exactly this
same authoritative primitive; historical `Fermi_Gas` classes are not used.

## 5. Neutral potential, conjugates, and Hessian

For canonical `x=(n_B,n_e,n_mu)`,

```text
n_p = n_e+n_mu,
n_n = n_B-n_e-n_mu,

epsilon_CN(x) = epsilon_n(n_n)+epsilon_p(n_p)
              + epsilon_e(n_e)+epsilon_mu(n_mu).
```

Direct differentiation gives

```text
g = d epsilon_CN/dx
  = (mu_n,
     -mu_n+mu_p+mu_e,
     -mu_n+mu_p+mu_mu)
  = (mu_n,-eta_npe,-eta_npmu).
```

Differentiating once more gives

```text
H = [ D_n,             -D_n,                     -D_n ]
    [ -D_n, D_n+D_p+D_e,               D_n+D_p ]
    [ -D_n,       D_n+D_p,        D_n+D_p+D_mu ].
```

This is derived from the potential, not introduced as an independent matrix
rule. Equivalently, with species order `(n,p,e,mu)`, diagonal intrinsic
`K=diag(D_n,D_p,D_e,D_mu)`, and the density Jacobian

```text
S_x = [[1,-1,-1], [0,1,1], [0,1,0], [0,0,1]],
H_x = S_x^T K S_x.
```

The displayed `H` is exactly symmetric. For positive active-species `D_i`,
`v^T H v=(S_x v)^T K(S_x v)>0` for every nonzero `v`, because `S_x` has full
column rank. Thus `H` is positive definite throughout the declared smooth
domain. No electrostatic projection and no inverse are present in production.

## 6. Equilibrium recovery and uniqueness

The active branch is reduced to one scalar common lepton chemical potential
`lambda=mu_e=mu_mu`. Test-independent inverse free-gas relations give

```text
n_e(lambda)  = [(lambda^2-M_e^2)^(1/2)/hbar_c]^3/(3 pi^2),
n_mu(lambda) = [(lambda^2-M_mu^2)^(1/2)/hbar_c]^3/(3 pi^2),
n_p(lambda)  = n_e+n_mu,
n_n(lambda)  = n_B-n_p,
F(lambda)    = mu_n(n_n)-mu_p(n_p)-lambda.
```

The bracket is

```text
lambda_low  = M_mu,
lambda_high = mu_n(n_B)-M_p.
```

For `n_B` strictly above muon onset and below the source ceiling, the lower
residual is positive and the upper residual is nonpositive. Increasing
`lambda` increases `n_p`, decreases `n_n`, decreases `mu_n`, and increases both
`mu_p` and the explicit subtracted `lambda`; therefore `F` is strictly
decreasing and the root is unique. A bracketed bisection stops at residual
`2e-12 MeV` or a width of 64 double epsilons times the chemical-potential
scale, with a hard 256-iteration limit. Nonfinite/unbracketed/nonconverged
states throw `std::runtime_error`.

`Evaluate(x)` never calls this solve and never enforces beta equilibrium.

## 7. Muon threshold and source ceiling

Three equilibrium regions are distinguished:

1. **Muon-free npe branch:** `n_mu=0`,
   `mu_n-mu_p=mu_e<M_mu`. It exists below onset but is not a regular 3-by-3
   active-species Hessian state.
2. **Threshold:** `mu_e=mu_n-mu_p=M_mu`, `n_mu=0`. The free-muon `D_mu`
   diverges, so no full Hessian is returned.
3. **Active npe-mu branch:** `n_mu>0` and
   `mu_n-mu_p=mu_e=mu_mu`; the full Hessian is regular.

At threshold the npe solution has the closed construction

```text
n_e = n_e(mu=M_mu),
n_p = n_e,
n_n = n_n(mu=mu_p(n_e)+M_mu),
n_B,onset = n_n+n_e.
```

Using the Zaki constants gives
`n_B,onset=0.45698480541241987 fm^-3`. `EquilibriumStateAt` supports only the
strict interval above this value. It fails closed below and on the threshold;
it does not return a singular fixed-dimension Hessian or insert a muon floor.

Off-equilibrium `Evaluate` has a distinct and explicit domain: it accepts any
charge-neutral state with `n_n,n_p,n_e,n_mu>0` below the source ceiling, even
when its `n_B` lies below the equilibrium muon onset. This preserves genuine
off-equilibrium semantics while refusing the `n_mu=0` boundary.

For the upper endpoint, FR2005 equation (71) gives first `Sigma-minus`
appearance when

```text
2 mu_n-mu_p-mu_Sigma = 0,
mu_Sigma(n_Sigma=0)=M_Sigma.
```

A second bracketed bisection follows the active beta-equilibrium solution and
finds `n_B,Sigma=0.61735520796653498 fm^-3`; an independent test solve finds
`0.61735520796652943 fm^-3`. The local source model requires `n_B` strictly
below this endpoint. Neither `Sigma-minus` nor `Lambda0` physics is added.

## 8. Dedicated RFG1-RFG11 validation

The self-contained executable is `rotochemical_trackr_freegas_local`, labelled
`rotochemical;thermodynamics;track-r;free-gas;self-contained;scientific;contract`
in `tests/CMakeLists.txt`. It uses no external EOS data and creates no
scientific baseline.

| Line | Result and measured evidence |
|---|---|
| RFG1 | PASS. Zaki masses `(n,p,e,mu)=(939.56542052,938.27208816,0.51099895,105.6583755) MeV`. Independent 32,768-panel long-double Simpson integrals plus separate `p_F`, `mu`, and reciprocal response formulas check n,p,e,mu at `x={1e-3,0.1,1,10}`. Maximum energy relative error `3.0231150144907791e-15`; other maximum `0`. Legacy electron/muon public results are bit-identical at four densities including zero. |
| RFG2 | PASS. Three arbitrary non-equilibrium states agree with an independent four-integral energy sum to maximum relative error `3.8594513925747193e-16`. |
| RFG3 | PASS. At the declared non-equilibrium fixture, independent/actual `eta_npe=-92.634662115730947 MeV` and `eta_npmu=-72.631533467752774 MeV`; public accessors are called directly. |
| RFG4 | PASS. Independent total-energy quadrature differentiated along every canonical axis. Four centered step halvings converge; final absolute errors are `3.1164108804659918e-7`, `4.7153037030511769e-7`, and `1.7942276997473527e-7 MeV`. |
| RFG5 | PASS. Every Hessian entry matches a separately evaluated `D_n,D_p,D_e,D_mu` oracle; maximum absolute error `0 MeV fm^3`. |
| RFG6 | PASS. Every held-fixed coordinate column uses four centered step halvings. Final maximum errors by column are `8.0916561273625121e-7`, `2.5147096948785475e-5`, and `2.7530119496077532e-5 MeV fm^3`. |
| RFG7 | PASS. On a 36-state active-domain grid, maximum asymmetry `0`; minimum first/second leading minors and determinant are `140.14373125628921`, `237340.63311412412`, and `555320814.43648303`; maximum infinity-norm condition number `51.4353890897873`. |
| RFG8 | PASS. Four active equilibrium densities have exact constructed charge neutrality, maximum `|eta|=1.1510792319313623e-12 MeV`, and positive symmetric energy rises in three independent composition directions; minimum rise `1.3918548802394071e-7 MeV fm^-3`. |
| RFG9 | PASS, load-bearing detector. Temporarily changing the production root relation from `mu_n-mu_p-lambda` to `mu_n-mu_p-1.01 lambda` made the focused executable fail with `Track-R free-gas equilibrium solve does not have a finite decreasing-sign bracket.` The source was restored byte-identically to SHA-256 `41891e1821fd9f3de6c326a1ef667e0ec2756a1e3e2f86aba24387e7d69ca459`. The mutation is not committed. |
| RFG10 | PASS. Independent npe recovery at `0.99 n_B,onset` gives `mu_e=105.03583860356072 MeV<M_mu`; independent onset values agree; threshold/below-threshold equilibrium requests fail; a state just above succeeds; `n_mu=0` has no derivative; `n_mu=1e-30 fm^-3` remains unfloored; off-equilibrium active evaluation below onset succeeds; the source ceiling fails closed. Metadata assertions pass. |
| RFG11 | PASS. At three physical-density states, independently diagonal `K` and explicit `S_x` give `S_x^T K S_x` with maximum error `0 MeV fm^3`. No production projector exists. |

The RFG4/RFG6 tolerances are below their observed converged errors by factors of
about 4 and 7, respectively; they are numerical differentiation tolerances,
not broad scientific tolerances. Analytic and direct formula comparisons use
floating-point-scale bounds.

## 9. Direct source numerical comparison boundary

**NO DIRECT LOCAL NUMERICAL SOURCE TABLE EXISTS IN THE CURRENT AUTHORITY SET.**

FR2005 Table 1 contains whole-star maximum mass, central mass-energy density,
radii, and Kepler period. Comparing this local provider directly to those
numbers would require a TOV/barotropic construction and rotation treatment,
which are expressly outside Phase 5A-3. No paper-number comparison is
fabricated. The provider is instead validated analytically as the exact local
source model.

## 10. Regression, artifacts, and scope audit

The required validation order was observed: clean focused Track-R RFG1-RFG11,
existing Phase 5A-2 V1-V10, complete self-contained suite, then complete full
suite. The complete suites used `-j1` and did not run concurrently.

| Validation | Authenticated result | CTest real time |
|---|---:|---:|
| Focused `rotochemical_trackr_freegas_local` | RFG1-RFG11 PASS | direct run |
| Existing `rotochemical_local_thermodynamics` | V1-V10 PASS | direct run |
| Self-contained inventory/suite | 24 registered; 24/24 PASS | 92.23 s |
| Full inventory/suite | 47 registered; 47/47 PASS | 684.10 s |

The full build cache records
`COMPACTSTAR_EOS_DATA_ROOT=/Users/keeper/Documents/CompactStar/data/compose`.
No test was skipped. The one-test increase in each inventory is exactly the new
self-contained Track-R executable.

All eight governed scientific artifacts are byte-identical to the starting
Git object:

| Artifact under `tests/baselines/` | SHA-256 |
|---|---|
| `baryon_number_dscmf1_reference.tsv` | `7b036942f2ae599ace3b2fc8b9a7d91f6d11b5899ab7e8a88d2bb4ee6686493b` |
| `grid_convergence_cmf_1p6_debug.tsv` | `2c68e2f7e871192e00322bcc12b6d8b13eb7fdbda4a6d8d7050f26f0271ce5eb` |
| `grid_convergence_cmf_1p6_trajectory.tsv` | `7c84557742ec0ec747756d118b2837fb54095a94c10194d4ccc2d0a778f5b04f` |
| `hartle_I_dscmf1_debug.tsv` | `a21d4c3f6d89322cb4cca7a073e0e142f180e529332d704b7b16a145eae741c9` |
| `hartle_monopole_dscmf1_debug.tsv` | `bd49e5a091ebcc59f7c4899422200181d4e71ecf552284840454d01aac4b8d52` |
| `passive_cooling_cmf_1p6_debug.tsv` | `afcad5d078fe458bdb441278ce56caa2d025becc6b489d00e9baad91a101c0de` |
| `tov_dscmf1_reference.tsv` | `3d9af9129a6a4ffde9e0f8c5507a160f968a861c0cf9f3b089cceecab86b701a` |
| `tov_path_equivalence_dscmf1.tsv` | `5c0f4b3bdb70921f8f2a869af10edc4d8f5ae3963a9d150e11ef859d21e1c678` |

No baseline was created or regenerated by this increment.

Scope inspection must remain clean of any change to `TOVSolver`, `NStar`,
`StarProfile`, `RotationSolver`, thermal/evolution code, `RotochemicalCache`,
`ChemState`, DS(CMF), APR/BPAL, or BNV. INV-09 remains `INTENDED BUT UNVERIFIED`
for global sequence/particle-number response; INV-11 remains `UNRESOLVED` for
the evolved/redshifted chemical state. This local implementation resolves
neither invariant.

TRACK-R FREE-GAS LOCAL THERMODYNAMICS IMPLEMENTED;
WHOLE-STAR FREE-GAS REPRODUCTION NOT YET PERFORMED.
