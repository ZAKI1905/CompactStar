# Phase-6A-1 controlled abstract BNV implementation campaign

Status: **CONTROLLED TRAJECTORY / NUMERICAL VALIDATION FAILED — CANDIDATE NOT
ACCEPTABLE**

Classification: **CONTROLLED MATHEMATICAL / ARCHITECTURE BNV CAMPAIGN; NOT A
PHYSICAL BNV MODEL; NOT A GOVERNED NUMERICAL BASELINE**.

## 1. Entry and immutable scope

`PHASE6A1I_ENTRY_SHA = 961dfa0de6f76df71df4cb98edc8e1b35a5c21b1`.

At entry, `HEAD`, local `master`, `origin/master`, and live
`refs/heads/master` were equal to that SHA. It is the canonical integration of the
owner-accepted Phase-6A-1 implementation plan. The implementation branch is
`physics/phase6a1-controlled-bnv-implementation` in
`/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase6a1-controlled-bnv-implementation`.

The complete entry hashes are recorded in
`docs/validation/phase6a1_controlled_bnv_entry_hashes.json`. All eleven governed
baselines, all 33 Phase-5D protected paths, the required Phase-5B/C/D authority and
candidate artifacts, the Structure-1 EOS/data inputs, and all 22 authenticated
literature-manifest entries passed their entry checks. These bytes are consumed,
not rewritten.

No physical BNV lifetime, rate, mixing, matrix element, particle mass, decay
channel, or microscopic model is selected. The values below are dimensionless
mathematical test inputs only. No A18, superfluidity, Regime-II, MixedStar,
sliding-background, variable-Z, or baseline-promotion work is authorized.

## 2. Exact implementation map

This table is the scratch/validation checklist required before production mutation.
Every new production owner is fail-closed and retains the accepted identities and
currentness dependencies.

| Proposed path | Change | Semantic owner / dependencies | State or RHS / units | BA gates and mutation detectors |
|---|---|---|---|---|
| `CompactStar/Analysis/EquilibriumBaryonTangent.hpp/.cpp` | new | sole typed `t(B0)` view of current Phase-5B sequence derivative | no ODE state; dimensionless count/count | BA2; M2/M17/M18 |
| `CompactStar/Physics/BNV/OrdinaryMatterSource.hpp` | new | atomic prescribed `B,Bdot,S_y`, event measure, fate/partition and identities | count, count/s; no RHS | BA1; M16/M19 |
| `CompactStar/Physics/BNV/MovingReferenceSource.hpp/.cpp` | new | sole `Sigma=S-tBdot`, `sigma=P Sigma` owner | count/s; no RHS | BA3-BA5; M1-M5/M17-M19 |
| `CompactStar/Physics/BNV/ProductFate.hpp` | new | terminal product fate and once-only event identity | dimensionless weights; MeV remains in partition | BA1/BA7/BA9; M15/M16 |
| `CompactStar/Physics/BNV/DirectEnergyLedger.hpp/.cpp` | new | actual-potential, generic event ledger, GR sum and direct MeV-to-erg boundary | typed direct power erg/s; no chemical RHS | BA7-BA9; M8-M10/M13-M16/M21 |
| `CompactStar/Physics/BNV/StaticZeroSpinHistory.hpp/.cpp` | new | current identity-bearing zero-spin history | `Omega=OmegaDot=0` | BA10b and identity/currentness faults |
| `CompactStar/Physics/BNV/FrozenBnvValidity.hpp/.cpp` | new | certificate and trial/checkpoint validity stop | diagnostics/refusal only | BA13; M17/M18/M20 |
| `CompactStar/Physics/BNV/FrozenControlledBnvRunContext.hpp/.cpp` | new | wrapper of shared const Phase-5D context; transferred AnalyticControl qualifications | same three-state layout | BA1-BA14/BA16-BA17 |
| `CompactStar/Physics/BNV/ControlledBnvSecularDriver.hpp/.cpp` | new | sole augmented RHS owner | `eta_dot=-Z(R+sigma)` and direct thermal term | BA4-BA12; RHS mutants |
| `CompactStar/Physics/BNV/BnvDiagnostics.hpp` | new | immutable owner-assembled diagnostics | accepted typed schema | BA14; no serializer recomputation |
| `CompactStar/Physics/BNV/CMakeLists.txt`, `CompactStar/Physics/CMakeLists.txt` | new/modified | build registration only | no legacy BNV dependency | build/dependency audit |
| `tests/bnv/controlled_neutron_sink_fixture.hpp` | new | controlled neutron history | exact `(Bdot,0,0)` | O1-O6; M1-M7/M17-M20 |
| `tests/bnv/energy_partition_fixtures.hpp` | new | P0/P1/P2, R18 and R-a/b/c independent fixtures | MeV/event | O7-O15; M8-M16/M21 |
| `tests/bnv/source_projection.cpp`, `direct_energy.cpp`, `thermal_ledger.cpp` | new | BA1-BA9 executable oracles | no artifact | M1-M21 as mapped |
| `tests/bnv/frozen_sensitivity.cpp` | new | 21-star static certificate | no trajectory | BA13/M20 |
| `tests/bnv/matched_control.cpp`, `coupled_trajectory.cpp`, `ode_refinement.cpp` | new | BA10-BA12 controls/trajectories | exact common grids | M11/M12/M20/M21 |
| `tests/bnv/qss_beta_bounds.cpp` | new | O16/O17 | diagnostics only | BA16/BA17 |
| `tests/bnv/produce_candidate.py`, `compare_candidate.py` | new | deterministic candidate producer/comparator | JSON candidate only | post-gate only |
| `tests/CMakeLists.txt` | modified | named BA registration | no scientific state | BA15 |
| this record and candidate JSON | new/updated | durable evidence; candidate only after gates | no production RHS | governance/audit |

The legacy `Physics/BNV.hpp`, `BNVState`, old BNV drivers, Microphysics BNV
sources, MixedStar, all governed Z/W/Q mathematics, and every protected Phase-5
path remain untouched. A need to modify any of them is an immediate plan-mismatch
stop.

## 3. Fixed fixture and frozen policy

- Structure-1 central density: `rho_c=1.1000000000000000e15 g cm^-3`.
- Radial resolution: `80000`; EOS resolution: `8192`.
- Initial `B0 = N_n+N_p = 7.6169065187188905e56 count`.
- Initial state: `Tinf=1.0000000000000000e8 K`, `eta_npe=0 MeV`,
  `eta_npmu=0 MeV`.
- State layout remains `(ln(Tinf/1e8 K),eta_npe_inf,eta_npmu_inf)`.
- Spin is OFF through `StaticZeroSpinHistory`; `RunPurpose::AnalyticControl`.
- Frozen quantities are `t`, `Z`, `Cstar`, enabled `Ltilde`, metric/background,
  surface/envelope data, `mu_B^infinity(B0)`, equilibrium direct-energy inputs,
  species support and domain. The actual potential still uses current `eta`.
- `B(t)` and `Bdot` are prescribed atomically; `B` is not an ODE state.
- Hard runtime ceiling: `abs(DeltaB)/B0 <= 1.0e-6`. The first failed
  certificate/currentness/depletion check invalidates that trial state and all
  later evolution; no such point may be serialized.
- Enabled processes are exactly `Me,Mmu`; `De,Dmu` are disabled.
- Required metric identity is
  `qualified Structure-1 radial80000 canonical nu/lambda`.
- Required normalization identities are
  `predeclared mathematical benchmark SMe=1e-51 erg cm^-3 s^-1 K^-8` and
  `predeclared mathematical benchmark SMmu=2e-51 erg cm^-3 s^-1 K^-8`.

## 4. Tangent and pre-result arithmetic

The only production tangent is
`t=(partial N_y^eq/partial B)_Omega` from the Phase-5B sequence derivative, using
the canonical denominator `B_B=B_n+B_p`. The reduced `B_n+B_e+B_mu` expression is
only a closure check.

| component | raw fixture value | propagated absolute numerical error |
|---|---:|---:|
| `t_n` | `0.9657700849496014` | `1.4428870413124226e-7` |
| `t_e` | `0.030852171225661786` | `4.3530476723815904e-9` |
| `t_mu` | `0.0033777438247248118` | `2.1785265203189970e-9` |

`B_B=1.6831408136820063e59 count/(geometric central-energy unit)` and its
conservative error is `1.2652133015974696e52`. The predeclared raw tangent closure
budget is `tau_t=1.508202854293702e-7`; the closed accessor derives
`t_n=1-t_e-t_mu` and must close within two binary64 ulps.

The governed linearized modified-Urca response at `Tinf=1e8 K` gives relaxation
eigen-times `3.526142343534e4 yr` and `2.877568989741e4 yr`. These are feasibility
estimates from the already-governed `Z`, `Ltilde`, and analytic `H_M'(0)`, not a
trajectory result or physical BNV inference. The accepted reaction-free estimate
is `(Z t_lepton)B0 approximately (108.9,276.3) MeV`.

## 5. Immutable mathematical run cards

One Julian year is exactly `31557600 s`. Every history has constant fractional
drive, uniform proper neutron-sink support proportional to `n_n`, exact source
`S_y=(Bdot,0,0)`, no external input, and the same initial state above. The
normalization `gamma=abs(Bdot)/integral_D e^Phi n_n dV` is mathematical only and
must never be labelled a physical rate.

| ID | reactions / partition | `(Bdot/B0)` | exact `Bdot` [count/s] | duration | predicted final `abs(DeltaB)/B0` | role |
|---|---|---:|---:|---:|---:|---|
| `RF-P0-TRANSIENT-v1` | reactions OFF; P0; thermal state held at its analytic-control initial value | `-1.0e-13 yr^-1` | `-2.4136520263641377e36` | `1.0e6 yr = 3.15576e13 s` | `1.0e-7` | reaction-free two-channel transient |
| `CPL-P0-TRANSIENT-v1` | Me/Mmu ON; P0; full thermal RHS | `-1.0e-13 yr^-1` | `-2.4136520263641377e36` | `1.0e5 yr = 3.15576e12 s` | `1.0e-8` | coupled small-imbalance P0 transient |
| `CPL-P1-TRANSIENT-v1` | Me/Mmu ON; P1; full thermal RHS | `-1.0e-13 yr^-1` | `-2.4136520263641377e36` | `1.0e5 yr = 3.15576e12 s` | `1.0e-8` | relativistic uniform-sea P1 trajectory |
| `CPL-P2-LINEAR-QSS-v1` | Me/Mmu ON; P2; full thermal RHS | `-1.0e-12 yr^-1` | `-2.4136520263641375e37` | `5.0e5 yr = 1.57788e13 s` | `5.0e-7` | full-retention trajectory and reachable linear-QSS fixture |

Each coupled card has a matched no-BNV control with identical initial state,
zero-spin owner, process selection, background, solver and checkpoints; only the
BNV source/direct bundle is zeroed.

The reaction-free card predicts
`abs(eta_e)=1.089e-5 MeV`, `abs(eta_mu)=2.763e-5 MeV`, far more than 100 times the
component ODE resolution and only one tenth of the hard depletion ceiling. The P2
QSS card lasts at least 14.2 and 17.4 initial-temperature linear relaxation times,
respectively, while using only half of the depletion ceiling. Its linear steady
estimate is approximately `eta=(-3.78e-6,-7.62e-6) MeV`, or
`xi=(-4.39e-4,-8.84e-4)` at `1e8 K`, safely below the required
`max abs(xi)<=0.5`. P2 direct retention is expected to prevent passive cooling from
making QSS less reachable; QSS is nevertheless accepted only from measured O16
metrics on the terminal interval `[4.0e5,5.0e5] yr`.

The old `xi approximately 4.9097` root-crossing target is expressly absent. Even
the maximum allowed reaction-free depletion gives only
`abs(xi_e) approximately 0.0126` at `1e8 K`. No card may be added, retuned, or
extended in response to a result. An impossible card causes STOP.

## 6. Energy partitions and finite-temperature omission

- P0: `Eesc_fluid=mu_n,actual` in the same local-frame, rest-mass-inclusive energy
  zero; cold direct residual is exactly zero.
- P1: local occupied relativistic Fermi-sea weight `w(p)=3p^2/pF^3` and the exact
  ADR-0015 R10 average. The `2/5 E_F,kin` expression is used only in a separate
  `pF/m_n << 1` nonrelativistic-limit fixture.
- P2: `Eesc_fluid=Eesc_star=E_X=0`; direct residual is exactly
  `mu_n,actual`.

All have `finite_T_direct_terms_included=false`. P0 reports the one-sided event
floor `(pi^2/6)(kBT)^2/E_F,kin`; P1 reports
`(pi^2/3)(kBT)^2/E_F,kin`. P2 uses the latter, larger smooth-weight Sommerfeld
scale as a conservative omission floor. Each is integrated through the same event
measure and reported in erg/s. A sign claim requires the complete matched-control
signal to exceed ten times the sum of numerical, quadrature, endpoint,
finite-temperature, partition, and frozen-background error bounds; otherwise its
classification is `SIGN_UNRESOLVED`.

## 7. Solver, grids, certificate, R18 and R20

The baseline RKF45 tolerances are `rtol=1e-7` and
`atol=(1e-12,1e-18,1e-18)`. Refined tolerances are `rtol=1e-9` and
`atol=(1e-14,1e-20,1e-20)`. Component convergence uses
`abs(y_base-y_refined)/(atol_base+rtol_base*max(abs(y_base),abs(y_refined)))<=1`.

`RF-P0-TRANSIENT-v1` uses 1025 uniformly spaced checkpoints. Each coupled
transient uses 2049 uniformly spaced checkpoints. `CPL-P2-LINEAR-QSS-v1` uses
8193 uniformly spaced checkpoints. Both tolerance levels use the identical grids.
R20 luminosity integrals use the 8193-point composite trapezoid and an independent
4097-point even-index recomputation. The quadrature estimate is their absolute
difference and must be `<=5e-5 N_R20`; endpoint-state propagation from the
baseline/refined difference must also be `<=5e-5 N_R20`. No missing or nonfinite
term is permitted. The raw acceptance remains
`abs(R20)/N_R20<=2e-4`, where

```text
N_R20=max(1 erg,
          abs(Delta U_th),
          abs(MeVToErg Delta Echem),
          integral (abs(Lnu_full)+abs(Lgamma)+abs(Lother)) dt).
```

Neither `abs(Delta Eeq)` nor `integral abs(L_out,fluid)dt` enters the normalizer.
The finite-T omission floor is reported separately and does not redefine the cold
mathematical fixture. `Eeq_dot=MeVToErg*mu_B^infinity*Bdot`.

R18 is tested at nonzero two-channel `eta` for P2 and a non-neutron
charge-balanced stoichiometric fixture:
`Pdir(actual)=Pdir(eq)+MeVToErg*eta^T sigma`. Its tolerance is the propagated
direct/source/tangent budget plus 64 ulps. R-a/R-b/R-c use the same event, rate,
energy zero and terminal fate and must agree within
`max(64 ulps,10*independent_quadrature_error)`. Adding R-b terms on top of R-c is
a required rejection.

The frozen certificate uses 21 independent stars at
`DeltaB/B0=0,-5e-8,...,-1e-6`. Both target residual and final bracket width obey
`tau_B,target=5e-11 B0`. Achieved `B_solved` values are regression abscissae.
At every star, `Cstar(Tinf)` and `Tsurface_inf(Tinf)` are evaluated on the exact
13-knot grid `log10(Tinf/K)=6.00,6.25,...,9.00`. The same knots are used for all
stars. The free-gas thermal source is linear in temperature, so this grid also
checks the temperature-independent `Cstar/Tinf` coefficient; the envelope result
is retained separately at every knot. A later controlled checkpoint outside
`1e6 K <= Tinf <= 1e9 K` invalidates the frozen certificate and causes STOP rather
than extrapolation. This grid completion was recorded before any BNV trajectory;
it changes no run-card drive or acceptance tolerance.
Uncertainty-weighted fits use the outer half `[-1e-6,-5e-7]`; every fit residual
must obey `abs(residual)<=3u_res+0.10T_X`, with full monotone drift, numerical and
nonlinearity envelopes still below the accepted threshold. Runtime validity is
checked within every RHS call, including trial states, and at every checkpoint.

## 8. BA1-BA17 immutable acceptance budgets

| Gate | Predeclared numerical requirement |
|---|---|
| BA1 | `abs(Bdot-b^T S)<=tau_S`, `tau_S=32 eps max(1,abs(Bdot),sum abs(S_i))`; constant-history closure exact and nonconstant closure `<=max(64 ulps,10qerr)` |
| BA2 | component oracles inside propagated errors above; raw closure `<=tau_t`; closed sum within two ulps; all identity/currentness faults refuse |
| BA3 | at least 100 finite cases; `abs(b^T Sigma)<=tau_b`; lift infinity residual `<=max_i tau_Sigma,i` with the accepted formulas |
| BA4 | both sigma lower bounds positive; both eta-dot upper bounds negative; normalized drive ratios agree with `(-1.4303e-55,-3.6281e-55) MeV/count` within propagated Phase-5B/C errors |
| BA5 | sliding source zero within projection/ODE budgets at both signs and three scales; raw-G and k mutants nonzero/detected |
| BA6 | RHS within 64 ulps plus input errors; integrated component scaled norm `<=1` |
| BA7 | P0/P2 exact; P1 `<=max(1e-12 relative,10qerr)`; R18 within propagated budgets plus 64 ulps |
| BA8 | named identities exact; R-a/b/c `<=max(64 ulps,10qerr)`; beta roots within `5e-10` absolute xi; Echem finite difference `<=max(1e-10 relative,10 step/roundoff estimate)` |
| BA9 | every M8-M16 and M21 changes/refuses beyond its BA7/BA8 budget on a nonzero fixture |
| BA10a | expected RHS identities exact; repeated bytes identical; governed-trajectory component scaled difference `<=1` |
| BA10b | every relaxed qualification refuses; zero-BNV wrapper RHS equals untouched Phase-5D RHS under the same zero-spin history |
| BA11 | every sample/currentness/ledger check passes; max frozen utilization `<=1`; R20 conditions in section 7 |
| BA12 | component scaled difference `<=1`; refined R20 `<=2e-4`; no final-decade step-collapse trend |
| BA13 | target/bracket `<=5e-11 B0`; fit rule above; all accepted quantity thresholds from the plan; just-over-limit trial/checkpoint refuses before serialization |
| BA14 | schema complete, typed and unit-tagged; serialized owner values bit-identical; independent ledger reconstruction meets BA8/BA11 budgets |
| BA15 | 33 protected paths and 11 baselines byte-identical; focused Phase-5B/C/D plus complete data-free/authenticated suites return zero; clean `git diff --check` |
| BA16 | on `[4.0e5,5.0e5] yr`, each component `abs(R+sigma)/max(abs(sigma),R_resolution)<=0.05` and `tau_relax/elapsed<=0.10` |
| BA17 | each enabled MU process obeys `-DeltaPbeta<=0.467659 Lnu_eq,M` with allowance `<=1e-10 max(Lnu_eq,M,1 erg/s)` |

The frozen thresholds are exactly those in the accepted preflight: componentwise
`abs(delta t_i)<=max(5u_ti,1e-4 abs(t_i))`; row-scaled Z infinity drift
`<=max(5E_Z,row,1e-4)` using `u_Z=E_Z`; `Cstar`, each enabled `Ltilde`, reference
potential/event averages, `N_i`, metric/structure, radius/surface/envelope and
`Tsurface_inf` each use their declared `1e-4`-plus-governed-error budgets; species
support and process/channel ordering are exact. Muon quantities are monitored
separately; `d ln N_mu/d ln B approximately 54` is only a fixture sensitivity
warning.

## 9. Mutation matrix and production exclusions

M1-M21 are required individually with the accepted primary detectors: M1 raw-G
(BA5); M2 k-for-t (BA2/BA5); M3 wrong sign (BA3/BA4); M4 omitted moving term
(BA3/BA5); M5 omitted channel (BA3/BA4/BA6); M6 Z transpose/cross (BA4/BA6); M7
channel swap (BA4/BA6); M8 double hole (BA7-BA9); M9 Echem heat (BA8/BA9); M10
PdV/gravity heat (construction/BA9); M11 omitted DeltaLnu (BA8/BA10b); M12
double equilibrium neutrinos (BA8/BA10a/b); M13/M14 missing/double conversion
(O12/BA7/BA8); M15 fluid/star confusion (BA7); M16 fate double-book (BA1/BA7);
M17 stale t (BA2/BA13); M18 domain mismatch (BA3); M19 Bdot mismatch (BA1/BA3);
M20 ignored validity (BA13); M21 equilibrium neutron potential at nonzero eta
(R18/O13 primary). R20 is not claimed as a primary M21 detector.

Production may contain neither a raw `G_y S_y` route nor a diagnostic-k route.
There is one source-to-sigma construction. There is no independent Fermi-hole,
Echem, PdV, or gravitational heat. Beta power retains its existing governed
MeV-to-erg boundary; direct event power has one separate typed conversion after
global MeV/s integration.

## 10. Candidate schema and trajectory gate

The candidate schema ID is
`compactstar.phase6a1.controlled-bnv-candidate.v1`. Rows contain all accepted typed
epoch/background, source, tangent/projection, chemical state/rates, chemical-energy
reservoir components, actual potentials, direct/fate and R18/R-a/b/c values,
separate beta/thermal luminosities, temperature/control differences, QSS metrics,
R20 residual/normalizer, frozen-validity metrics and provenance identities. Run
metadata includes entry/predeclaration/implementation SHAs, all source/input hashes,
exact run cards, solver/grid configuration, baseline/refined results, individual
BA and M results, and physical-rate/model flags `false`.

No BNV trajectory may be generated until this declaration is committed and
BA1-BA10b, all applicable pretrajectory mutations, and the static BA13 certificate
pass with protected bytes unchanged. A durable PRETRAJECTORY PASS/STOP entry will
be appended before candidate production. Drive values, tolerances, grids, error
budgets and classifications in this declaration are immutable after this commit.

Any retained artifact is only:
`docs/validation/phase6a1_controlled_bnv_candidate.json`, visibly classified
`CONTROLLED MATHEMATICAL / ARCHITECTURE BNV CANDIDATE; NOT PHYSICAL BNV MODEL; NOT
GOVERNED BASELINE; NOT OWNER-RATIFIED NUMERICAL RESULT`.

## 11. Post-run disposition

The initial form of this record was committed at
`PHASE6A1I_PREDECLARATION_SHA=6a3c3e8653f407ba9734968bb7a1802cb78a5549`
before any BNV trajectory. The immutable cards, grids and tolerances above were not
changed afterward. The durable pretrajectory gate is
`docs/validation/PHASE6A1_CONTROLLED_BNV_PRETRAJECTORY.md` at commit
`460ca713a6bd2633b4afbf3926e959b751cea0c6`.

All four declared cards and their matched controls completed at both tolerance
levels with the exact declared row counts. Independent BA12 comparison then found
a maximum component-scaled difference of `1.7255917120989046` in `x_state` for
`CPL-P2-LINEAR-QSS-v1`, at `t=23113476562.5 s = 732.421875 yr`. The baseline and
refined values were `0.051132258160285306` and `0.05113226698535251`. The immutable
BA12 limit is one. An independent compiled comparator reproduced the refusal.

The trajectory campaign stopped with disposition:

**CONTROLLED TRAJECTORY / NUMERICAL VALIDATION FAILED — CANDIDATE NOT
ACCEPTABLE**.

No tolerance was loosened, no drive was retuned, and no success-labelled candidate
or governed baseline was created. BA16 and trajectory BA17 were not run after the
mandatory BA12 stop. Complete failure evidence is in
`docs/validation/PHASE6A1_CONTROLLED_BNV_NUMERICAL_FAILURE.md`.

Final governed regression then exposed a second, material architectural conflict:
the accepted exact `CompactStar/Analysis/EquilibriumBaryonTangent.hpp/.cpp` paths
necessarily enter Phase-5D's recursive `CompactStar/Analysis` scientific-source
provenance, so the fresh Phase-5D artifact cannot equal its governed baseline.
Changing that path, producer, comparator or baseline is outside this authority.
The single final disposition is therefore:

**IMPLEMENTATION EXPOSED A MATERIAL SCIENTIFIC / ARCHITECTURE CONFLICT — RETURN
TO OWNER**.
