# Phase 5A-2 — Local thermodynamic contract implementation

> **Status (2026-09-04): PHASE 5A-2 LOCAL THERMODYNAMIC CONTRACT IMPLEMENTED AND VERIFIED —
> TRACK-R EOS REPRODUCTION PREREQUISITES READY.**
>
> This increment is local and self-contained. It supplies no realistic hadronic EOS and performs
> no stellar or secular-evolution calculation.

## 1. Starting point and authority

| Item | Authenticated value |
|---|---|
| Starting SHA | `1e6dcd8c5cdb28c9d6d959e42e9e9745d54aeac8` |
| Branch | `physics/phase5a-local-thermodynamics` |
| Worktree | `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-rotochemical-local-thermo` |
| Governing decision | `docs/adr/ADR-0010-rotochemical-off-equilibrium-thermodynamic-contract.md` — **ACCEPTED** |
| Change class | Scientific-semantic + structural/API + engineering/test + documentation; strictest governance applied |

The canonical master and `origin/master` were equal to the starting SHA and clean before the fresh
branch/worktree was created. ADR-0010 is the accepted authority; this task did not reopen its owner
decisions.

## 2. Production boundary introduced

`CompactStar/EOS/LocalThermodynamics.hpp` and
`CompactStar/EOS/src/LocalThermodynamics.cpp` introduce the smallest reusable local boundary:

- `ChargeNeutralCoordinates` — canonical `x=(n_B,n_e,n_mu)` coordinates;
- `ChargeNeutralCompositionState` — validated state with reconstructed species densities;
- `NeutralConjugates` — `g=(mu_n,-eta_npe,-eta_npmu)`;
- `ChargeNeutralChemicalHessian` — `H_ab=partial g_a/partial x_b`;
- `LocalThermodynamicEvaluation` — validated state, energy density, `g`, and `H`;
- `LocalThermodynamicProviderMetadata` — model identity, revision, composition, coordinates,
  temperature scope, rest-mass/lepton ownership, and smooth domain;
- `ILocalThermodynamicProvider` — evaluation plus equilibrium-state recovery;
- `ColdRelativisticFreeLepton` — analytic cold electron or muon component.

The interface is independently usable. It has no dependency on `TOVSolver`, `NStar`,
`StarProfile`, `RotationSolver`, `ThermalSolver`, `RotochemicalCache`, or `ChemState`.
It contains no matrix inversion or electrostatic projector and exposes no full paper-`B` object.
Individual charged intrinsic potentials are intentionally not generic outputs.

## 3. State semantics and units

The canonical coordinates are

`x=(n_B,n_e,n_mu)` in fm^-3,

with exact construction identities

`n_p=n_e+n_mu`, and `n_n=n_B-n_e-n_mu`.

Construction requires finite values, `n_B>0`, `n_e>=0`, `n_mu>=0`, and reconstructed
`n_p>=0`, `n_n>=0`. Inputs are never clipped. Concrete providers can further restrict their
smooth active-species domain.

The returned energy density is in MeV fm^-3; `g` is in MeV; and
`H=partial g/partial x` is in MeV fm^3, with each derivative holding the other two canonical
coordinates fixed. `NeutralConjugates::EtaNpeMeV()` and `EtaNpmuMeV()` apply the accepted signs:
`eta_npe=-g_1`, `eta_npmu=-g_2`.

## 4. Analytic free electrons and muons

The implementation uses the repository's authoritative Zaki constants
`ELECTRON_M_MEV`, `MUON_M_MEV`, and `MEV_2_INV_FM`; hence
`hbar c=1/MEV_2_INV_FM` in MeV fm. For `n_l>0`, it evaluates

`p_F=hbar c (3 pi^2 n_l)^(1/3)`,

`mu_l=sqrt((m_l c^2)^2+p_F^2)`,

`d mu_l/d n_l=p_F^2/(3 mu_l n_l)`,

and the closed-form zero-temperature relativistic Fermi integral for total energy density,
including rest mass. A small-`p_F/m` analytic series prevents cancellation in that integral; it
does not clip the density.

At `n_l=0`, the meaningful value limit is explicit: `p_F=0`, `mu_l=m_l c^2`, and energy density
is zero. The singular/non-smooth derivative is `std::nullopt`; requesting it as a scalar throws.
Negative, NaN, and infinite density fail closed. No threshold smoothing is present.

## 5. Test-only analytic charge-neutral toy EOS

The grouped test owns the toy. **Its coefficients are test physics, not a neutron-star EOS and not
a scientific baseline.** With `d=x-x_0`,

- `x_0=(0.20,0.04,0.02)` fm^-3;
- `epsilon_0=200` MeV fm^-3;
- `g_0=(950,0,0)` MeV;
- `lambda=10` MeV fm^6;
- `H_0=[[2,-2,-2],[-2,10,5],[-2,5,12]]` MeV fm^3;
- `epsilon(x)=epsilon_0+g_0 dot d+(1/2)d^T H_0 d+(lambda/6)d_B^3`.

Therefore `g=g_0+H_0 d+(lambda/2)d_B^2 e_B` and
`H=H_0+lambda d_B e_B e_B^T`. The nonzero cross derivatives detect diagonal-only errors; the
cubic term gives a measurable second-order linearization remainder. The declared test domain is
`0.15<=n_B<=0.25` fm^-3 with `n_n,n_e,n_mu>0`. The equilibrium composition is analytic:
`n_e=0.04+(14/95)(n_B-0.20)`, `n_mu=0.02+(10/95)(n_B-0.20)`.

The toy additionally offers a concrete optional intrinsic-potential split solely to verify
`S_x^T mu=g`; that split is not part of `ILocalThermodynamicProvider`.

## 6. Coordinate transformation checks

For `y=(n_n,n_e,n_mu)` the test uses the accepted linear map

`y=T x`, `T=[[1,-1,-1],[0,1,0],[0,0,1]]`, and `U=T^-1`.

It verifies the matrices, not only observables:

`H_x=T^T H_y T`, and `C_y=T H_x^-1 T^T=H_y^-1`.

For nonlinear fraction coordinates `z=(n_B,Y_e,Y_mu)`, the test constructs the complete Hessian
transformation. Away from equilibrium it proves that the mixed entries contain the nonzero
chain-rule terms `g_e` and `g_mu`; `J^T H J` alone is therefore rejected.

## 7. Independent corrected-2006 projection fixture

The V9 oracle is separate from the provider API. In species order `(n,p,e,mu)`, it defines the
positive intrinsic Hessian `K=diag(2,3,5,7)` MeV fm^3, analytic
`chi=K^-1`, and charge vector `q=(0,1,-1,-1)^T`. Test-side code independently constructs

`C_CN=chi-chi q (q^T chi q)^-1 q^T chi`.

The charge-neutral density map `S_x` produces the provider toy's `H_0=S_x^T K S_x`. A transparent
test-only 3-by-3 cofactor inverse then checks

`C_CN=S_x H_0^-1 S_x^T`,

`C_CN q=0`, and row(`p`)=row(`e`)+row(`mu`). This detects missing projection, double projection,
basis reduction errors, and charge-sign errors. No form of this projection is present in the
generic provider or production source.

## 8. V1–V10 results

The grouped executable is `rotochemical_local_thermodynamics`; its labels are
`rotochemical;thermodynamics;self-contained;scientific;contract`.

| Line | Result | Defect falsified |
|---|---|---|
| V1 | PASS | Unit declaration or `eta`/neutral-conjugate sign reversal |
| V2 | PASS; charge/baryon identities exact to arithmetic precision | Non-neutral or clipped reconstruction |
| V3 | PASS; maximum equilibrium composition-conjugate residual `2.78e-17` MeV | Failure to recover beta equilibrium |
| V4 | PASS; independent free-lepton maximum scaled error `8.92e-16` | Wrong constants, rest mass, relativistic formula, or derivative |
| V5 | PASS; toy energy, gradient, and Hessian exact | Provider/energy inconsistency or missing cross derivatives |
| V6 | PASS; maximum Hessian asymmetry `0` | Maxwell/mixed-partial violation |
| V7 | PASS; errors `5.0000e-4, 1.2500e-4, 3.1250e-5, 7.8125e-6` MeV; minimum observed order `1.999999994` | Incorrect local linear response or absent convergence |
| V8 | PASS; `H` transform error `0`, inverse-response error `5.55e-17`; nonlinear chain terms nonzero | Wrong basis map or omission of nonlinear coordinate-Hessian term |
| V9 | PASS; projection equivalence, charge-null, proton-row residuals each `1.39e-17` or smaller | Missing/double electrostatic elimination, wrong reduction, wrong charge signs |
| V10 | PASS over the interval endpoints/reference point; minimum Sylvester values `1.5,11,94.5`, maximum `kappa_inf=23.925926`; invalid/singular/boundary cases rejected | Unstable/rank-deficient fixture, silent invalid-state acceptance, or muon-onset regularization |

Analytic comparisons use floating-point-scale bounds where cancellation is not intrinsic. V7 uses
convergence rather than one permissive tolerance. V10 reports the actual infinity-norm condition
metric. The only inverse is transparent and test-side; it fails closed on the singular fixture.

## 9. Complete regression evidence

The required order was observed: new local test, then complete self-contained suite, then complete
authenticated full suite. The suites were serial (`-j1`) and were not run concurrently.

| Configuration | Registered/passed | CTest real time | Result |
|---|---:|---:|---|
| Final local executable before complete suites | V1-V10 | Direct run | PASS |
| Complete self-contained | 23/23 | 91.14 s | PASS |
| Complete full, authenticated `COMPACTSTAR_EOS_DATA_ROOT=/Users/keeper/Documents/CompactStar/data/compose` | 46/46 | 686.34 s | PASS |

The full count is authenticated from `ctest -N`, not assumed. No test skipped and no scientific
baseline was created or regenerated.

## 10. Existing governed artifacts

All eight existing files under `tests/baselines/` remain byte-identical to their starting hashes:

| Artifact | SHA256 |
|---|---|
| `baryon_number_dscmf1_reference.tsv` | `7b036942f2ae599ace3b2fc8b9a7d91f6d11b5899ab7e8a88d2bb4ee6686493b` |
| `grid_convergence_cmf_1p6_debug.tsv` | `2c68e2f7e871192e00322bcc12b6d8b13eb7fdbda4a6d8d7050f26f0271ce5eb` |
| `grid_convergence_cmf_1p6_trajectory.tsv` | `7c84557742ec0ec747756d118b2837fb54095a94c10194d4ccc2d0a778f5b04f` |
| `hartle_I_dscmf1_debug.tsv` | `a21d4c3f6d89322cb4cca7a073e0e142f180e529332d704b7b16a145eae741c9` |
| `hartle_monopole_dscmf1_debug.tsv` | `bd49e5a091ebcc59f7c4899422200181d4e71ecf552284840454d01aac4b8d52` |
| `passive_cooling_cmf_1p6_debug.tsv` | `afcad5d078fe458bdb441278ce56caa2d025becc6b489d00e9baad91a101c0de` |
| `tov_dscmf1_reference.tsv` | `3d9af9129a6a4ffde9e0f8c5507a160f968a861c0cf9f3b089cceecab86b701a` |
| `tov_path_equivalence_dscmf1.tsv` | `5c0f4b3bdb70921f8f2a869af10edc4d8f5ae3963a9d150e11ef859d21e1c678` |

## 11. Explicit limits and next gate

This increment does not supply a Track-R realistic literature EOS or a Track-P source-matched
DS(CMF) provider. It does not integrate a local response over a star, calculate particle-number
response, reduce a baryon-conserving sequence, construct paper `Z/W`, define an evolved chemical
state, or add weak reactions/heating/evolution. INV-09 and INV-11 remain unresolved.

**NO APR, DS(CMF), STELLAR SUSCEPTIBILITY, Z/W, ROTOCHEMICAL EVOLUTION,
SUPERFLUIDITY, OR BNV IMPLEMENTED IN PHASE 5A-2.**

## 12. Post-review validation correction

**Current status (2026-09-04): PHASE 5A-2 VALIDATION GAPS CORRECTED — READY FOR
CANONICAL INTEGRATION.** The independent review's material F1-F3 validation
findings were corrected without a permanent production diff. The load-bearing
detectors, independent y/z and lepton-energy oracles, direct public Eta accessor
checks, optional all-axis V5-V7 strengthening, artifact hashes, and final serial
23/23 self-contained and 46/46 full suites are recorded in
`docs/validation/PHASE5A2_LOCAL_THERMODYNAMIC_VALIDATION_CORRECTION.md`.

The historical V9 wording above is narrowed by that correction: V9 validates
the correct single local charge-neutral reduction/projection, response
amplitude, charge null, proton-row identity, and charge-sign convention. Source
inspection separately establishes that production contains no second
projection; an identical idempotent projector cannot be universally detected
by this fixture.
