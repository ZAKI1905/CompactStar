# Phase 5A-5 — final Track-R low-density local thermodynamic coverage

> **TRACK-R FREE-GAS LOCAL THERMODYNAMIC COVERAGE COMPLETE — WHOLE-STAR STRUCTURE REPRODUCTION READY**
>
> This is a local source-model result only. The Fernandez--Reisenegger 2005
> whole-star benchmark has **not** been reproduced. The implementation remains
> an **UNVERIFIED SCIENTIFIC CANDIDATE** pending human domain ratification under
> `GOVERNANCE.md` section 5.

## 1. Authentication and starting state

| Item | Authenticated value |
|---|---|
| Starting SHA | `41ab66bd6e6f351691173a1e2a033c646ffd3772` |
| Canonical checkout | `/Users/keeper/Documents/CompactStar/repo/CompactStar` |
| Canonical state | clean `master`; HEAD = `master` = `origin/master` = starting SHA after fetch |
| Branch | `physics/phase5a-trackr-pe-branch` |
| Fresh worktree | `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-rotochemical-trackr-pe` |
| Execution | one agent; complete suites serial and nonoverlapping |
| Change class | scientific-semantic, numerical policy, structural/API, tests/build registration, documentation |
| Governing contract | accepted ADR-0010; Phase-5A-4 independent review is the low-density gate authority |

The branch and worktree did not exist before creation. The fresh worktree was
clean at the exact canonical SHA. No divergent branch commit touching the task
files was found. This increment adds analytic local physics and no governed
scientific baseline; `GOVERNANCE.md` section 3.1 is not invoked.

Before editing, the existing executables authenticated:

```text
n_B,n-onset     = 7.3567289037328318e-9 fm^-3
n_B,mu-onset    = 0.45698480541241987 fm^-3
n_B,Sigma-onset = 0.61735520796653498 fm^-3
```

The pre-change classifier returned `ProtonElectron`, `NeutronOnset`, `Npe`,
`MuonOnset`, and `NpeMuon` on the ordered positive-density intervals and exact
thresholds. The first two alternatives deliberately threw as the missing gate.
All three threshold values are preserved.

## 2. Direct source authority

The three local PDFs were read directly. Extracted equations were checked
against rendered pages where layout or signs mattered.

| Source | Direct use |
|---|---|
| Fernandez & Reisenegger 2005 | PDF pp. 8-9, section 3.1 states that the noninteracting Fermi gas is used for the whole star. Unlike APR/BPAL, no separate crust supplement is assigned to this free-gas case. PDF p. 12 requires each species contribution only over the region where that species exists. |
| Reisenegger et al. 2006 | PDF pp. 2-3, equations (2), (4), (9)-(13): total chemical potentials include electrostatic terms, local charge neutrality removes the charged singular direction, and the intrinsic beta combinations cancel the electrostatic potential. The existing neutral local Hessian therefore receives no second projection. |
| Reisenegger 1995 | PDF pp. 14-15, section 2 and equation (3): `n_p=n_e`; minimizing the electron-inclusive energy gives `mu_p+mu_e-mu_n=0` when neutrons are active. |

Authenticated PDF SHA-256 values are:

```text
1995  9af85e37c7a52fd5b704c0ba07cc0ad89741d23b049df31cb6867d501d91d0ff
2005  f184d7d1d7030b61a021eb5c7ac14b1f1b30c7ea69e9d53473d153cfb069ea88
2006  a286f15e083e52becd95b3000cbb5ec3ed97148681cf10a43f1a1cc5c4d23ae8
```

Atomic binding, nuclei, Coulomb lattices, and a realistic crust are not added.
They would define a different source model.

## 3. Complete active-set ladder

The implemented source-specific equilibrium dispatch is now:

| Density domain | Physical state | Returned response |
|---|---|---|
| `n_B<0` or nonfinite | invalid | fail closed |
| `n_B=0` | vacuum | `VacuumBoundaryEvaluation`, values only |
| `0<n_B<n_B,n` | p,e | `PeThermodynamicEvaluation`, 1D smooth response |
| `n_B=n_B,n` | p,e; neutron at appearance | `NeutronThresholdEvaluation`, values only |
| `n_B,n<n_B<n_B,mu` | n,p,e | `NpeThermodynamicEvaluation`, 2D smooth response when numerically resolved |
| `n_B=n_B,mu` | n,p,e; muon at appearance | `MuonThresholdEvaluation`, values only |
| `n_B,mu<n_B<n_B,Sigma` | n,p,e,mu | existing 3D `LocalThermodynamicEvaluation` |
| `n_B>=n_B,Sigma` | outside FR2005 free-gas source model | fail closed |

Classification uses ordered `<`, `==`, and `>` comparisons to the stored
physical thresholds. No tolerance widens or relabels an active set. Numerical
unavailability is reported separately as `EquilibriumResolutionError`.

## 4. p-e derivation and active coordinate

For `0<n_B<n_B,n`, positivity, baryon accounting, and neutrality give

```text
n_p=n_e=n_B,  n_n=n_mu=0.
```

There is one independent thermodynamic direction:

```text
z_pe=(n_B)  [fm^-3],
epsilon_pe(n_B)=epsilon_p(n_B)+epsilon_e(n_B)  [MeV fm^-3].
```

Using `d epsilon_i/d n_i=mu_i`,

```text
h_pe = d epsilon_pe/d n_B = mu_p+mu_e  [MeV],
H_pe = d h_pe/d n_B = D_p+D_e         [MeV fm^3],
D_i  = d mu_i/d n_i.
```

`PeCoordinates`, `PeConjugates`, `PeChemicalHessian`, and
`PeThermodynamicEvaluation` encode exactly this 1x1 response. Active particles
are `{n=false,p=true,e=true,mu=false}`, response dimension is 1, and status is
`SmoothInterior`. No inactive row, zero padding, fake 2x2/3x3 matrix, neutron
floor, muon floor, or off-equilibrium p-e composition coordinate is introduced.

## 5. Inactive-neutron diagnostic

Below appearance, neutron beta equilibrium is not an active equation. The p-e
result carries only the value diagnostic

```text
eta_npe,threshold = m_n-mu_p-mu_e.
```

It is positive below onset, strictly decreases, and reaches zero at onset. It
is not part of `PeConjugates` or `PeChemicalHessian`. No finite value is assigned
to `d mu_n/d n_n` at zero neutron density.

## 6. Neutron-appearance semantics

At `n_B=n_B,n`,

```text
n_n=n_mu=0,  n_p=n_e=n_B,n,  mu_p+mu_e=m_n.
```

`NeutronThresholdEvaluation` carries the physical state, energy density,
limiting p-e conjugate, zero neutron diagnostic, active-particle status, and
`response_dimension=0`. It has no Hessian member. The zero dimension means
"no smooth response returned" while the active chart changes from 1D to 2D;
it is not a zero matrix. Exact lower and upper `nextafter` values remain p-e and
npe classified respectively. The upper value fails explicitly when unresolved.

The independent test-side threshold expression gives
`7.356728903732833408335e-9 fm^-3` at platform precision, consistent with the
Phase-5A-4 review authority's independent 80-digit
`7.3567289037328326656352e-9 fm^-3` and the unchanged provider double
`7.356728903732831753973e-9 fm^-3` at order-one-ULP scale.

## 7. Vacuum semantics

At `n_B=0`, `VacuumBoundaryEvaluation` contains an all-zero
`VacuumCompositionState`, `epsilon=0`, and the one-sided limiting conjugate

```text
h_pe -> m_p+m_e.
```

It has `response_dimension=0`, no Hessian member, no active particles, and
status `VacuumBoundary`, not `SpeciesThreshold`. `H_pe=D_p+D_e` grows without
bound as `n_B->0+`; no finite boundary derivative is fabricated. The historical
strictly-positive `ChargeNeutralCompositionState` invariant is unchanged, so
the composition-only `EquilibriumStateAt(0)` directs consumers to the typed
`EquilibriumAt(0)` boundary result.

## 8. p-e to npe limiting response

For `z_npe=(n_B,n_e)`, the p-e tangent is `t=(1,1)^T`. Write
`A=D_p+D_e`. The existing npe Hessian is

```text
H_npe = [[ D_n,       -D_n ],
         [ -D_n, D_n + A ]].
```

Therefore, without a normalization factor,

```text
t^T H_npe t = A = H_pe.
```

Its determinant is `D_n A`, and direct 2x2 inversion gives

```text
H_npe^-1 = [[1/A+1/D_n, 1/A],
            [1/A,       1/A]].
```

As `n_n->0+`, `D_n->infinity`, hence

```text
H_npe^-1 -> (t t^T)/H_pe.
```

The factor is exactly one because `t` maps the scalar displacement
`delta n_B` to `delta z=t delta n_B`; using a unit tangent would instead move
the normalization into that coordinate definition. PE-V9 tests the production
tangent contraction and exact inverse entries, then continues an independent
stable momentum construction closer to onset to verify rank-one collapse. It
adds the neutron-side absent-species embedding falsifier. It does not close the
distinct muon-side Layer-B embedding obligation recorded as N-4.

## 9. N-1 numerical treatment

N-1 is closed by a quantified fail-closed response rule, not a density floor.
The canonical state reconstructs `n_n=n_B-n_e`; near neutron appearance,
relative information in that subtraction is eventually insufficient for
`D_n=H_00 proportional to n_n^(-1/3)`.

For any requested npe response, production computes the downward local spacing

```text
ulp_B = n_B-nextafter(n_B,0).
```

It returns a smooth 2x2 response only when

```text
n_n >= 2^30 ulp_B.
```

Otherwise the density remains physically `Npe` classified and an
`EquilibriumResolutionError` states that the Hessian is unavailable. No nearby
density, p-e response, threshold result, or positive neutron floor is returned.
The composition-only API remains available where the equilibrium root itself is
representable, including the geometric probes rejected only for Hessian
accuracy.

Against an independent long-double common-momentum inversion, responses at
`10^-3`, `10^-4`, `10^-5`, and `10^-6` relative above onset are available; the
maximum allowed `H_00` relative error is `3.600691649285892e-7`. Points from
`10^-7` through `10^-10` fail explicitly under the exact rule. The existing
stricter bracket/residual guard continues to reject even closer unrepresentable
roots. Classification is never tolerance-based.

## 10. PE-V1 through PE-V13

Every test states the defect it falsifies and uses test-side phase-space
formulae, the independently supplied 80-digit threshold reference plus a
separate test-side closed-form expression, the pressure identity, finite
perturbations, analytic 1x1/2x2 algebra, or stable common-momentum inversion.
No production onset routine is used as the PE-V6 oracle.

| Test | Result and independent evidence |
|---|---|
| PE-V1 | PASS — seven logarithmic densities; exact p/e composition, inactive n/mu, 1D status |
| PE-V2 | PASS — independent quadrature energy maximum relative error `6.67e-16`; maximum conjugate difference `1.14e-13 MeV` |
| PE-V3 | PASS — independent `D_p+D_e`; maximum relative error `6.67e-16` |
| PE-V4 | PASS — centered perturbation minimum order `1.9993158683`; final relative derivative error `2.2060625748e-7` |
| PE-V5 | PASS — neutron diagnostic positive and strictly decreasing; `3.6402309433e-7 MeV` at `0.999999` onset |
| PE-V6 | PASS — independent expression agrees with the review's 80-digit `7.3567289037328326656352e-9 fm^-3` reference; production unchanged |
| PE-V7 | PASS — independent pressure identity gives `1.8964875026324584e-9 MeV fm^-3`, consistent with the review's `1.8964875026317866e-9`; positive for every tested density and `P/n_B->0` |
| PE-V8 | PASS — exact threshold value-only; deliberate `nextafter` sides; no 1x1 or 2x2 threshold Hessian |
| PE-V9 | PASS — `t^T H_npe t=H_pe` to `4.45e-16`; exact inverse normalization and independent rank-one limit |
| PE-V10 | PASS — all-zero vacuum, `epsilon=0`, finite conjugate limit, divergent one-sided Hessian, no boundary Hessian |
| PE-V11 | PASS — four resolved/four explicit failed geometric probes; allowed `H_00` error below `1e-6`; composition remains available where reliable |
| PE-V12 | PASS — invalid/vacuum/p-e/neutron/npe/muon/npe-mu/Sigma dispatch ladder |
| PE-V13 | PASS — muon/Sigma values and established 2D/3D paths preserved |

PE-V7 closes N-5 by asserting positive pressure at exact neutron onset. N-7 is
also obsolete because the source-specific `EquilibriumAt` now performs the
six-way typed dispatch rather than duplicating the generic base body.

## 11. Existing regression and complete-suite evidence

The required execution order was focused PE, then R1, RFG, V, then the complete
self-contained suite, then the complete authenticated full suite. Complete
suites used `-j1` and did not overlap.

| Validation | Result | Runtime |
|---|---:|---:|
| PE-V1-PE-V13 | PASS | direct focused run |
| R1-V1-R1-V10 | PASS | direct run |
| RFG1-RFG11 | PASS | direct run |
| V1-V10 | PASS | direct run |
| Complete self-contained | 26/26 PASS | 86.46 s CTest (86.48 s wall) |
| Complete authenticated full | 49/49 PASS | 669.86 s CTest (669.87 s wall) |

The full configuration used
`COMPACTSTAR_EOS_DATA_ROOT=/Users/keeper/Documents/CompactStar/data/compose`.
The nine inputs governed by `HEAT_CAPACITY_V1.md` and the four cold files
governed by `TOV_REFERENCE.md` (one overlaps) matched their recorded SHA-256
values before the full run. CTest inventories were authenticated from the live
builds. The added analytic test is self-contained; no scientific baseline was
added.

## 12. Governed baselines

All eight governed artifacts are byte-identical before and after validation:

| Artifact | SHA-256 |
|---|---|
| `baryon_number_dscmf1_reference.tsv` | `7b036942f2ae599ace3b2fc8b9a7d91f6d11b5899ab7e8a88d2bb4ee6686493b` |
| `grid_convergence_cmf_1p6_debug.tsv` | `2c68e2f7e871192e00322bcc12b6d8b13eb7fdbda4a6d8d7050f26f0271ce5eb` |
| `grid_convergence_cmf_1p6_trajectory.tsv` | `7c84557742ec0ec747756d118b2837fb54095a94c10194d4ccc2d0a778f5b04f` |
| `hartle_I_dscmf1_debug.tsv` | `a21d4c3f6d89322cb4cca7a073e0e142f180e529332d704b7b16a145eae741c9` |
| `hartle_monopole_dscmf1_debug.tsv` | `bd49e5a091ebcc59f7c4899422200181d4e71ecf552284840454d01aac4b8d52` |
| `passive_cooling_cmf_1p6_debug.tsv` | `afcad5d078fe458bdb441278ce56caa2d025becc6b489d00e9baad91a101c0de` |
| `tov_dscmf1_reference.tsv` | `3d9af9129a6a4ffde9e0f8c5507a160f968a861c0cf9f3b089cceecab86b701a` |
| `tov_path_equivalence_dscmf1.tsv` | `5c0f4b3bdb70921f8f2a869af10edc4d8f5ae3963a9d150e11ef859d21e1c678` |

## 13. Scope audit and changed files

Production changes are limited to:

- `CompactStar/EOS/LocalThermodynamics.hpp`
- `CompactStar/EOS/TrackRFreeGasThermodynamics.hpp`
- `CompactStar/EOS/src/TrackRFreeGasThermodynamics.cpp`

Test/build changes are limited to:

- `tests/eos/rotochemical_trackr_pe.cpp`
- `tests/eos/rotochemical_trackr_npe.cpp`
- `tests/CMakeLists.txt`

Documentation changes are this record plus the current-state pointers in
`docs/MODERNIZATION_ROADMAP.md`, `docs/SCIENTIFIC_INVARIANTS.md`, and
`docs/architecture/CURRENT_ARCHITECTURE.md`.

No `TOVSolver`, `StarProfile`, `NStar`, `RotationSolver`, `ThermalSolver`,
`RotochemicalCache`, `ChemState`, BNV source, EOS data, or baseline changed.
There is no second electrostatic projection, no tolerance-based active-set
reclassification, and no premature claim that the whole-star benchmark was
reproduced.

## 14. Remaining review findings

- R2, R3, R5, and R9 remain separate future hardening.
- N-2 remains the conservative off-equilibrium npe response-domain restriction.
- N-3 remains a latent default-interface hazard for a future second provider.
- N-4 remains the distinct muon-side absent-species susceptibility embedding
  obligation; PE-V9 adds the neutron-side analogue but does not replace it.
- N-6 remains the documented approximately 15-ULP muon-onset construction
  rounding, fully mitigated by existing fail-closed behavior.
- N-1 is closed by the explicit `2^30`-ULP response rule and PE-V11.
- N-5 is closed by exact-onset pressure assertion in PE-V7 and R1-V9.
- N-7 is closed because the override is no longer redundant.

INV-09 remains **INTENDED BUT UNVERIFIED**. INV-11 remains **UNRESOLVED**.
No global particle-number response, redshift convention, or evolution
convention is supplied by local coverage.

## 15. Exact coverage claim and non-goals

The local Track-R free-gas model now supplies the source-faithful active-set
ladder from the zero-pressure vacuum boundary through the pre-Sigma-minus
ceiling. It can support a subsequent whole-star structure reproduction task.
It does **not** mean that a TOV star, FR2005 global benchmark, stellar
susceptibility, or rotochemical coefficient has been constructed or verified.

NO WHOLE-STAR TOV REPRODUCTION, STELLAR SUSCEPTIBILITY INTEGRATION,
PAPER B/Z/W, ROTOCHEMICAL EVOLUTION, APR/BPAL, DS(CMF)
OFF-EQUILIBRIUM PHYSICS, SUPERFLUIDITY, OR BNV IS IMPLEMENTED IN
PHASE 5A-5.
