# Phase 5A-2C — Local thermodynamic validation correction

> **Status (2026-09-04): PHASE 5A-2 VALIDATION GAPS CORRECTED —
> READY FOR CANONICAL INTEGRATION.**
>
> This is a test/evidence correction. It changes no production thermodynamic
> API or formula and adds no realistic EOS, stellar integration, coefficient,
> evolution, or BNV implementation.

## 1. Authentication, authority, and scope

| Item | Authenticated value |
|---|---|
| Starting SHA | `46c425457b8ccb9e29d5fc35f15edb6f95796891` |
| Starting branch/upstream/live remote | `physics/phase5a-local-thermodynamics`; all equal after live fetch |
| Starting master/origin master | `1e6dcd8c5cdb28c9d6d959e42e9e9745d54aeac8`; both equal |
| Starting relationship | 2 ahead / 0 behind master |
| Starting tree | clean, including untracked files |
| Governing authority | accepted ADR-0010 local thermodynamic contract (`docs/adr/ADR-0010-rotochemical-off-equilibrium-thermodynamic-contract.md:5`, `:178-199`, `:236-273`) |
| Correction authority | Astra findings F1-F3 (`docs/validation/PHASE5A2_LOCAL_THERMODYNAMIC_INDEPENDENT_REVIEW.md:247-313`, `:469-486`) |
| Change class | engineering/test validation plus documentation; no scientific-semantic production change |

Astra found no incorrect production thermodynamic physics. The integration
blocker was inadequate load-bearing validation: self-confirming V8 coordinate
checks, no V4 coverage of the production small-x lepton-energy branch, and no
direct V1 calls to the public imbalance accessors
(`docs/validation/PHASE5A2_LOCAL_THERMODYNAMIC_INDEPENDENT_REVIEW.md:3-10`,
`:247-313`). ADR-0010 and its accepted local contract were not redesigned.

## 2. F1 — independent coordinate oracles

### Linear source basis

The corrected V8 independently encodes the Hessian obtained by differentiating
the declared toy energy in `y=(n_n,n_e,n_mu)`:

```text
H_y = [[2,0,0],[0,8,3],[0,3,10]] + 10 (n_B-0.20) ones(3,3).
```

No x-basis Hessian or inverse coordinate map enters this oracle
(`tests/eos/rotochemical_local_thermodynamics.cpp:640-655`). At three
non-equilibrium states, V8 compares production-toy `H_x` with
`T^T H_y T`, and compares `T H_x^-1 T^T` with the independently formed
`H_y^-1` (`tests/eos/rotochemical_local_thermodynamics.cpp:673-711`). The final
focused run measured maximum Hessian error `2.2204460492503131e-16` and maximum
inverse-response error `1.1102230246251565e-16`.

### Nonlinear fraction basis

For `z=(n_B,Y_e,Y_mu)`, V8 independently encodes the complete closed-form
Hessian obtained after substituting `n_e=n_B Y_e` and `n_mu=n_B Y_mu` directly
into the toy energy (`tests/eos/rotochemical_local_thermodynamics.cpp:657-671`).
It compares that oracle with the complete chain rule
`J^T H_x J + sum_a g_a d2x_a/dz2`, then separately compares the naive
`J^T H_x J` (`tests/eos/rotochemical_local_thermodynamics.cpp:713-741`). The
complete transform agreed to `8.3266726846886741e-17`; the naive transform
disagreed by `0.084000000000000019` at the chosen non-equilibrium state.

### V8 load-bearing detectors

| Detector | Temporary change | Required observed result | Reversion |
|---|---|---|---|
| D8-1 | Replaced accepted `T` with identity in V8 | focused executable failed: `V8 x Hessian disagrees with independently differentiated y oracle` | restored accepted `T`; corrected test file returned to its pre-detector SHA-256 |
| D8-2 | Removed all four nonlinear `g_e/g_mu` second-coordinate additions | focused executable failed: `V8 complete nonlinear transform disagrees with independent z oracle` | restored all four additions; corrected test file returned to its pre-detector SHA-256 |

The permanent identity-map negative fixture also differs from the independent
oracle by `2.3599999999999999`; it is not used to construct either expected
matrix (`tests/eos/rotochemical_local_thermodynamics.cpp:675-712`). F1 is
corrected.

## 3. F2 — free-lepton small-x branch

V4 now generates density from each independently specified target
`x=p_F/m` by

```text
n = [m x/(hbar c)]^3/(3 pi^2),
```

using test-side arithmetic (`tests/eos/rotochemical_local_thermodynamics.cpp:359-367`).
Both electrons and muons use the same explicit target grid:

| Region | Electron target x | Muon target x |
|---|---|---|
| well below branch | `1e-5`, `1e-3` | `1e-5`, `1e-3` |
| just below | `0.009999` | `0.009999` |
| branch point | `0.01` | `0.01` |
| just above | `0.010001` | `0.010001` |
| relativistic | `1` | `1` |

The energy oracle numerically evaluates the normalized momentum integral with
32,768-panel composite Simpson quadrature and compensated summation. It copies
neither production's small-x series nor its cancellation-prone antiderivative
(`tests/eos/rotochemical_local_thermodynamics.cpp:315-344`). The measured
maximum true relative energy error was
`6.945933048592142e-13`; the committed guard is `8e-13`, narrowly above the
independently measured platform maximum. Maximum relative error among momentum,
chemical potential, and reciprocal-`dn/dmu` derivative checks was
`9.4154034448246174e-16` (`tests/eos/rotochemical_local_thermodynamics.cpp:450-528`).

Independent centered finite differences of the quadrature energy verify
`mu=d epsilon/dn` with step-halving convergence for both species at `x=1e-3`
and `x=1`. The final relative errors are below `1e-8`; the actual absolute-error
sequences are printed by the focused executable
(`tests/eos/rotochemical_local_thermodynamics.cpp:493-521`).

For the load-bearing detector, production's small-x `x^5` coefficient was
temporarily changed from `4/5` to `8/5`. V4 failed with
`V4 analytic lepton energy-density mismatch`. The production file was restored
to starting SHA-256
`0e57ed70a3b11cf588a4dc672a1607d9aadff4bf9e2c8e1715f79fbe86292ca3`
and has no final diff (`CompactStar/EOS/src/LocalThermodynamics.cpp:21-31`). F2
is corrected.

## 4. F3 — public imbalance accessors

At `x=(0.21,0.043,0.021) fm^-3`, direct differentiation of the declared toy
polynomial independently predicts

```text
eta_npe  = -0.015 MeV,
eta_npmu = -0.007000000000000027 MeV.
```

V1 calls `NeutralConjugates::EtaNpeMeV()` and `EtaNpmuMeV()` directly and
compares them with those nonzero analytic values; actual values matched exactly
in the final focused run (`tests/eos/rotochemical_local_thermodynamics.cpp:380-414`).
The optional charged-potential consistency identity remains a separate check.

For the load-bearing detector, `EtaNpeMeV()` was temporarily changed to return
`+value_MeV[1]`. V1 failed with actual `+0.015000` versus expected `-0.015000`.
The header was restored to starting SHA-256
`3ed370b1d18b9a2c63adbc320e7635e75a17de3cb0524147011b9bf9e5d3311f`
and has no final diff (`CompactStar/EOS/LocalThermodynamics.hpp:89-94`). F3 is
corrected.

## 5. Optional V5-V7 strengthening

The small strengthening was straightforward and included. V5 now differentiates
the evaluated energy independently along all three canonical coordinate axes;
maximum gradient error was `2.3489974410040304e-09`
(`tests/eos/rotochemical_local_thermodynamics.cpp:531-572`). V7 separately holds
the other two coordinates fixed and checks every Hessian column; maximum error
was `5.0495148062879025e-09`, while the existing mixed-direction convergence
retained minimum observed order `1.9999999944156268`
(`tests/eos/rotochemical_local_thermodynamics.cpp:585-638`).

## 6. V9 claim correction

V9's valid algebra is unchanged. The supported statement is: V9 validates the
correct single local charge-neutral reduction/projection, response amplitude,
charge null, proton-row identity, and charge-sign convention. Source inspection
separately proves production contains no second projection. An identical
idempotent projector cannot be universally detected by this fixture. The
corrected test comment states this precise boundary
(`tests/eos/rotochemical_local_thermodynamics.cpp:782-789`). The earlier broad
claim in the implementation record remains historical and is superseded by its
post-review correction pointer; the Astra review is unchanged.

## 7. Final validation

The required order was followed: permanent focused test; D8-1; D8-2; small-x
series detector; Eta-accessor detector; clean focused rerun; complete
self-contained suite; complete full suite. The complete suites were serial and
were not run concurrently.

| Validation | Authenticated result | Runtime |
|---|---:|---:|
| Clean focused `rotochemical_local_thermodynamics` | V1-V10 PASS | direct run |
| Complete self-contained inventory/suite | 23 registered; 23/23 PASS | 91.38 s |
| Complete full inventory/suite | 46 registered; 46/46 PASS | 682.88 s |

The full build cache retained the authenticated
`COMPACTSTAR_EOS_DATA_ROOT=/Users/keeper/Documents/CompactStar/data/compose`.
No test skipped and no new scientific baseline was created.

## 8. Governed artifacts and production diff

All eight existing scientific baselines remain byte-identical:

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

The final permanent production diff relative to
`46c425457b8ccb9e29d5fc35f15edb6f95796891` is empty. Permanent changes are
limited to tests and documentation. Detector mutations are fully reverted.

## 9. Explicit non-implementation statement

APR implemented: **NO**. DS(CMF) off-equilibrium physics implemented: **NO**.
Stellar susceptibility integration: **NO**. Paper Z/W: **NO**. Evolution:
**NO**. BNV: **NO**. No free-gas star was implemented.

**PHASE 5A-2 VALIDATION GAPS CORRECTED —
READY FOR CANONICAL INTEGRATION**

**Exactly one recommended next action:** fast-forward the complete Phase 5A-2
implementation, independent review, and validation correction to canonical
master; do not begin Track-R free-gas work automatically.
