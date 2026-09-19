# Phase-6A-1 controlled BNV recovery preflight

Status: **PROPOSED — OWNER ACCEPTANCE REQUIRED**

Disposition: **PHASE-6A-1 RECOVERY PREFLIGHT COMPLETE — BA15 ARCHITECTURE
FIX AND BA12R NUMERICAL REQUALIFICATION PLAN READY FOR OWNER REVIEW**.

Classification: **DOCUMENTATION-ONLY RECOVERY PLAN; NOT AN IMPLEMENTATION;
NOT A PHYSICAL BNV MODEL; NOT A TRAJECTORY; NOT A GOVERNED BASELINE**.

## 1. Authenticated entry and non-authority

The canonical authority remains
`961dfa0de6f76df71df4cb98edc8e1b35a5c21b1`.  The failed implementation evidence
is retained at
`e56e6e50040dbcd9dcbee1acecf58843f3dddf1c` on
`physics/phase6a1-controlled-bnv-implementation`; its implementation record identifies
the same canonical entry and branch/worktree (`docs/validation/PHASE6A1_CONTROLLED_BNV_IMPLEMENTATION.md:9-24`),
and its failure record preserves the implementation/evidence lineage
(`docs/validation/PHASE6A1_CONTROLLED_BNV_NUMERICAL_FAILURE.md:16-27`).

This proposal is prepared on
`analysis/phase6a1-bnv-recovery-preflight` in
`/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase6a1-bnv-recovery-preflight`,
created from the exact failed implementation SHA.  It changes no equation, source,
test, build file, baseline, EOS/data byte, literature byte, run card, or numerical
result.  Governance classifies this permanent diff as documentation and requires
line-cited factual claims (`GOVERNANCE.md:43-57`).  No trajectory was run.

The authenticated entry state was:

| checkout/ref | authenticated value | state |
|---|---|---|
| canonical local `master`, `origin/master`, live `refs/heads/master` | `961dfa0de6f76df71df4cb98edc8e1b35a5c21b1` | clean |
| failed local branch, upstream, live branch ref | `e56e6e50040dbcd9dcbee1acecf58843f3dddf1c` | clean |
| recovery branch entry | `e56e6e50040dbcd9dcbee1acecf58843f3dddf1c` | clean |

The successful scientific implementation is not rejected.  ADR-0015 continues to
govern the moving-reference seam and explicitly preserves Phase-5B/C/D coefficient
mathematics and baselines (`docs/adr/ADR-0015-bnv-open-system-thermal-ledger.md:14-21`,
`docs/adr/ADR-0015-bnv-open-system-thermal-ledger.md:162-179`).  This plan addresses
only the two recorded blockers.

## 2. Immutable failure evidence

### 2.1 BA12 remains FAIL

The accepted two-level BA12 test was

```text
abs(y_baseline-y_refined)
----------------------------------------------- <= 1,
atol_baseline + rtol_baseline max(abs(y_baseline),abs(y_refined))
```

with baseline `rtol=1e-7`, `atol=(1e-12,1e-18,1e-18)` and refined
`rtol=1e-9`, `atol=(1e-14,1e-20,1e-20)`
(`docs/validation/PHASE6A1_CONTROLLED_BNV_IMPLEMENTATION.md:166-180`).  The
compiled comparator implements exactly that componentwise comparison
(`tests/bnv/ode_refinement.cpp:22-28`).

For `CPL-P2-LINEAR-QSS-v1` at `t=23113476562.5 s = 732.421875 yr`, the
immutable values are:

| quantity | value |
|---|---:|
| baseline `x_state` | `0.051132258160285306` |
| refined `x_state` | `0.05113226698535251` |
| absolute difference | `8.8250672047873735e-9` |
| relative difference in `x_state` | `1.7259291881808068e-7` |
| accepted scaled difference | `1.7255917120989046` |
| accepted limit | `1` |

The stored values, failed norm, independent-comparator refusal, and the fact that all
other cards and controls passed are historical evidence
(`docs/validation/PHASE6A1_CONTROLLED_BNV_NUMERICAL_FAILURE.md:53-79`).  This
proposal does not relabel them.  In particular, changing the old limit from one to two
would select a tolerance from the observed `1.7255917` result and is prohibited.

The state definition is `x_state=ln(Tinf/1e8 K)`
(`docs/validation/PHASE6A1_CONTROLLED_BNV_IMPLEMENTATION.md:67-76`).  Direct
binary64 evaluation therefore gives

```text
Tinf_baseline = 105246208.07772338 K
Tinf_refined  = 105246209.00652823 K
Delta Tinf    = 0.92880484461784363 K
relative Delta Tinf = 8.8250669870705891e-9.
```

This is a small physical-temperature separation, not evidence of a large physical
discrepancy.  It is also not a numerical PASS: the predeclared BA12 metric failed.

### 2.2 BA15 remains FAIL

The Phase-5D governed artifact is
`tests/baselines/phase5d1_controlled_evolution.json`, SHA-256
`2606916915b2da5c051b1a06637c0a371a74751c6e77b2775bc629a63bc9f6dd`.
The governed producer, comparator and fresh regression are respectively
`fresh_context.py`, `compare_artifacts.py`, and
`phase5d1_controlled_evolution_regression.py`
(`docs/validation/PHASE5D1_CONTROLLED_ROTOCHEMICAL_EVOLUTION_INTEGRATION.md:35-59`).
The fresh producer records production files under `CompactStar/Physics` and
`CompactStar/Analysis` as equality-bearing scientific provenance
(`tests/rotochemical/validate_trajectory.py:102-105`,
`tests/rotochemical/artifact_schema.py:108-120`).  The comparator intentionally has
no generic ignore-path escape and requires scientific source-provenance equality
(`docs/validation/PHASE5D1_GOVERNED_ARTIFACT_PREPARATION_PLAN.md:150-176`).

The implementation added exactly these two entries to the recursively authenticated
`CompactStar/Analysis` tree:

| added path | SHA-256 |
|---|---|
| `CompactStar/Analysis/EquilibriumBaryonTangent.hpp` | `7218b1b337037ed11697dca7c321f595be201679fa8683a21c03bcd5c2931771` |
| `CompactStar/Analysis/src/EquilibriumBaryonTangent.cpp` | `a26b6e45f6c3e42a49e42d2d84293b9cf3bbc5657ab114ee274940772aa2d265` |

The fresh artifact consequently contained 93 production hashes rather than the
governed 91, while no existing governed entry changed
(`docs/validation/PHASE6A1_CONTROLLED_BNV_NUMERICAL_FAILURE.md:117-140`).  The
protected baselines, Phase-5D paths, Phase-5 authorities, EOS/data and literature
remained byte-identical (`docs/validation/PHASE6A1_CONTROLLED_BNV_NUMERICAL_FAILURE.md:140-142`).
BA15 is therefore an **architectural ownership/provenance conflict**, not a changed
Phase-5 scientific result.

## 3. BA15 owner-decision candidate

### 3.1 Proposed owner and paths

Move the typed Phase-6 adapter to the existing Phase-6 BNV production module and
rename its namespace owner:

```text
CompactStar/Physics/BNV/EquilibriumBaryonTangent.hpp
CompactStar/Physics/BNV/src/EquilibriumBaryonTangent.cpp
CompactStar::Physics::BNV::EquilibriumBaryonTangent
```

The `src/` placement follows the module's existing header/source layout: installed
BNV interfaces are in `CompactStar/Physics/BNV`, while implementation units are in
`CompactStar/Physics/BNV/src` (`CompactStar/Physics/BNV/CMakeLists.txt:1-14`).  It is
cleaner than putting the `.cpp` beside the public header and makes the object visibly
a Phase-6 consumer of Phase-5B authority rather than a new Phase-5 `Analysis` owner.

The future recovery must remove the two new `Analysis`-path files, update only the
Phase-6 BNV includes/type names/build registration and corresponding Phase-6 tests,
and add no provenance exclusion.  It must not change the Phase-5D producer,
comparator, governed baseline, or protected Phase-5 source.  The Phase-5D baseline
contract permits no source-provenance change
(`docs/validation/PHASE5D1_CONTROLLED_ROTOCHEMICAL_EVOLUTION_INTEGRATION.md:73-94`).

### 3.2 Mathematics that must remain identical

The relocated class remains a typed adapter over the current governed
`Analysis::EquilibriumSequenceNumberDerivative`; it does not solve a sequence and
does not consume `G_y`.  Its mathematics remains

```text
t = (partial N_y^eq / partial B)_Omega,
B_B = B_n + B_p,
t_n = B_n/B_B,
t_e = B_e/B_B,
t_mu = B_mu/B_B.
```

ADR-0011 owns the whole-star baryon count `N_B=N_n+N_p`, domain-qualified sequence
objects, and complete-star `B_i` derivative
(`docs/adr/ADR-0011-particle-number-structural-response.md:44-84`), while ADR-0015
requires `t_i=B_i/B_B` with the same domain/surface policy, errors and currency
(`docs/adr/ADR-0015-bnv-open-system-thermal-ledger.md:181-203`).  The reduced
`B_n+B_e+B_mu` remains only the charge-closure check.  `k` and raw `G_y S_y` remain
forbidden production routes (`docs/adr/ADR-0015-bnv-open-system-thermal-ledger.md:189-203`).

The future recovery must demonstrate, preferably bit-identically and otherwise
within the already accepted component budgets:

- `t=(0.9657700849496014, 0.030852171225661786,
  0.0033777438247248118)` with propagated errors
  `(1.4428870413124226e-7, 4.3530476723815904e-9,
  2.1785265203189970e-9)`;
- the same canonical `B_B`, closure budget, closed representation, sequence/star/
  domain/currentness refusals;
- bit-identical source-to-`sigma` results where the compiler permits it;
- unchanged BA2, BA3, BA4 and BA5 results; and
- unchanged M2, M17 and M18 detection.

These fixture values and budgets were committed before trajectories
(`docs/validation/PHASE6A1_CONTROLLED_BNV_IMPLEMENTATION.md:88-110`).  No physical
precision beyond those budgets is implied.

### 3.3 Phase-5D proof obligation

After relocation, a fresh empty-scratch Phase-5D governed regression must produce

```text
fresh scientific_production_source_hashes
    == governed scientific_production_source_hashes
fresh governed artifact SHA-256
    == 2606916915b2da5c051b1a06637c0a371a74751c6e77b2775bc629a63bc9f6dd
```

with no baseline edit.  The regression already refuses self-comparison and compares
only after a fresh producer run (`tests/rotochemical/phase5d1_controlled_evolution_regression.py:32-42`,
`tests/rotochemical/phase5d1_controlled_evolution_regression.py:166-226`).  Failure
of either equality stops recovery before any BNV trajectory.

### 3.4 ADR disposition

ADR-0015 already fixes the physics, equations, `t` authority, seam and separation of
future owners (`docs/adr/ADR-0015-bnv-open-system-thermal-ledger.md:137-179`,
`docs/adr/ADR-0015-bnv-open-system-thermal-ledger.md:397-409`).  It does not fix a C++
namespace or filesystem path.  Nevertheless the proposed recovery changes module
ownership from `Analysis` to `Physics::BNV`; governance classifies moving ownership
or boundaries as structural/architecture and requires an ADR plus a same-change
`CURRENT_ARCHITECTURE` update (`GOVERNANCE.md:43-57`).

Therefore **ADR-0016 is required before R0**, unless the owner instead chooses a
formally accepted ADR-0015 amendment that explicitly authorizes this exact ownership
relocation.  A new, narrow ADR-0016 is preferred for audit clarity.  It must state
that no ADR-0015 physics, Phase-5B derivative authority, Phase-5D producer/comparator,
or governed baseline changes.  This preflight proposes that decision but does not
make it.

## 4. BA12 numerical analysis

`ScaledRKF45` allocates GSL's embedded RKF45 stepper and scaled controller; its exact
per-component requested level is `D_i=atol_i+rtol*abs(y_i)`
(`CompactStar/Physics/Rotochemical/ScaledRKF45.hpp:12-18`,
`CompactStar/Physics/Rotochemical/ScaledRKF45.hpp:27-50`).  The controller's estimate
is a **local per-step** estimate used to accept/reject and resize a step.  It does not
mathematically imply

```text
abs(y_base-y_refined) /
(atol_base+rtol_base max(abs(y_base),abs(y_refined))) <= 1
```

for accumulated global trajectory error.  Global error also depends on the number
and distribution of steps, propagation/amplification of local defects through the
nonlinear RHS, the embedded estimator's constants, different adaptive step
sequences, output constraints, and floating-point error.  The old BA12 comparison is
a useful conservative empirical cross-solution heuristic, but it is not a theorem
about GSL RKF45.  The official
[GSL ODE documentation](https://www.gnu.org/software/gsl/doc/html/ode-initval.html)
likewise describes RKF45 as embedded (4,5) and the controller as applying the
requested level to the stepper's local error estimate.  Its failure remains
dispositive under the accepted campaign.

## 5. Proposed BA12R hierarchy

BA12R is a new gate requiring explicit owner acceptance before execution.  It does
not replace or retroactively alter BA12.

| tier | `rtol` | `atol` for `(x_state,eta_e,eta_mu)` | role |
|---|---:|---|---|
| BASELINE | `1e-7` | `(1e-12,1e-18,1e-18)` | preserved historical run |
| REFINED | `1e-9` | `(1e-14,1e-20,1e-20)` | proposed nominal candidate tier |
| ULTRA | `1e-11` | `(1e-16,1e-22,1e-22)` | convergence witness only |

The baseline/refined tiers are immutable campaign inputs
(`docs/validation/PHASE6A1_CONTROLLED_BNV_IMPLEMENTATION.md:166-175`).  The exact
ULTRA tier is numerically sensible for the current binary64 state scales and is
therefore adopted without alteration.  At the failed checkpoint the refined and
ULTRA `x_state` requested scales are `5.1142266985352513e-11` and
`5.1142266985352506e-13`; 64 ulps of the state magnitude are
`4.440892098500626e-16`.  Across the stored P2 trajectory the lepton states reach
approximately `4.86e-7 MeV` and `1.11e-6 MeV`, for which the proposed ULTRA requested
scales remain hundreds of 64-ulp floors.  The tier is not a claim that roundoff or
stiffness is absent: the BA12R resolution and step-trend refusals below test that.

No ULTRA run is authorized by this document.  If the exact tier later produces a
roundoff plateau, step collapse, or solver refusal, BA12R fails; the tolerances are
not changed after observing it.

## 6. Exact BA12R acceptance logic

All comparisons use the existing identical 8193-point P2 output grid
(`docs/validation/PHASE6A1_CONTROLLED_BNV_IMPLEMENTATION.md:173-179`).  For every
state component `i` and common checkpoint `j`, define

```text
d_BR(i,j) = abs(y_BASELINE(i,j) - y_REFINED(i,j))
d_RU(i,j) = abs(y_REFINED(i,j)  - y_ULTRA(i,j))
M(i,j)    = max(abs(y_BASELINE), abs(y_REFINED), abs(y_ULTRA))

D_R(i,j) = atol_REFINED,i + rtol_REFINED M(i,j)
D_U(i,j) = atol_ULTRA,i   + rtol_ULTRA   M(i,j)

Q(i,j) = sum of the half-ulp uncertainties introduced by parsing the serialized
         BASELINE, REFINED and ULTRA binary64 values; Q=0 if values are compared
         before lossy formatting.

F(i,j) = max(D_U(i,j), 64 ulp(M(i,j)), Q(i,j)).
```

For `M=0`, `ulp(M)` is the binary64 spacing at zero.  BA12R passes the state gate
only if **all** of the following hold:

1. **Refined-to-ULTRA stability**

   ```text
   d_RU(i,j) / D_R(i,j) <= 1
   ```

   at every component/checkpoint.  This is a predeclared empirical refined-accuracy
   budget, not a claim that local tolerance guarantees global error.

2. **Resolvable contraction**: where `d_BR(i,j) > 10 F(i,j)`,

   ```text
   d_RU(i,j) / d_BR(i,j) <= 0.10.
   ```

3. **Floor-limited behavior**: where `d_BR(i,j) <= 10 F(i,j)`, no contraction ratio
   is formed; instead require

   ```text
   d_RU(i,j) <= 10 F(i,j).
   ```

The tolerance ratio between adjacent tiers is `0.01`.  For an adaptive embedded
RKF 4(5) method in an asymptotic truncation regime, the standard local-to-global
scaling estimate is approximately `tol^(4/5)`, predicting an adjacent-tier factor
`0.01^(4/5)=0.0251188643150958`.  The proposed `0.10` permits a factor of about four
for non-asymptotic step placement, nonlinear amplification and output constraints,
while still requiring decisive contraction.  It is chosen from the method/tolerance
hierarchy, not from a future ULTRA result.  GSL identifies the configured stepper as
embedded Runge-Kutta-Fehlberg (4,5), and the repository selects that exact GSL type
(`CompactStar/Physics/Rotochemical/ScaledRKF45.hpp:44-48`).

### 6.1 Independently converged ledger observables

State convergence alone is insufficient.  For each checkpoint define the gross
power scale

```text
G_P = max(1 erg/s,
          abs(P_dir_actual)+abs(L_H)+abs(DeltaLnu)+abs(Lnu_eq)
          +abs(Lgamma)+abs(Lother)+abs(L_out_fluid)).
```

For each of `P_dir_eq`, `P_dir_actual`, `L_H`, `DeltaLnu`, `DeltaPbeta`, `Lnu_eq`,
`Lnu_full`, `L_out_fluid`, `Lgamma`, `Lother`, and `Pnet`, require

```text
abs(O_REFINED-O_ULTRA) / G_P <= 1e-9.
```

Also apply the same `d_RU/d_BR <=0.10` rule when the baseline/refined difference is
larger than
`10 max(1e-11 G_P,64 ulp(max(abs(O_BASELINE),abs(O_REFINED),abs(O_ULTRA))))`;
otherwise require the refined/ULTRA difference to be at or below that floor.  The
power list is the accepted separated ledger rather than an ambiguous BNV-heating
scalar (`docs/validation/PHASE6A1_CONTROLLED_BNV_IMPLEMENTATION_PREFLIGHT.md:1035-1050`).

For endpoint energies use `max(N_R20_REFINED,N_R20_ULTRA)` as the scale and require
the refined/ULTRA difference in each of `Delta Eeq`, converted `Delta Echem`, and
`Delta Uth` to be `<=1e-9` of that scale.  Independently require, for both REFINED
and ULTRA:

- `abs(R20)/N_R20 <= 2e-4`;
- luminosity quadrature error `/N_R20 <=5e-5`;
- thermal-energy quadrature error `/N_R20 <=5e-5`; and
- refined/ULTRA endpoint-state propagation error `/max(N_R20) <=5e-5`.

These retain the accepted R20 and subsidiary budgets
(`docs/validation/PHASE6A1_CONTROLLED_BNV_IMPLEMENTATION.md:176-190`).  R18 and
R-a/R-b/R-c must retain their existing nonzero-eta/roundoff budgets, and every
REFINED/ULTRA trial and checkpoint must remain frozen-valid; the accepted gate maps
those identities to BA7/BA8/BA11/BA13
(`docs/validation/PHASE6A1_CONTROLLED_BNV_IMPLEMENTATION_PREFLIGHT.md:967-975`).

Matched source-minus-control values of `Tinf`, `Tsurface_inf`, `Lgamma`, `Uth`, and
time-integrated total power use the same REFINED/ULTRA stability and resolvable-
contraction logic, with `G_P` or `N_R20` as dimensionally applicable.  No heating or
cooling label follows unless the already predeclared sign-resolution rule is met
(`docs/validation/PHASE6A1_CONTROLLED_BNV_IMPLEMENTATION.md:157-164`).

### 6.2 Solver-trend refusal

Define the final temporal decade as `[t_final/10,t_final]`.  From the `.steps`
tables, count accepted steps in each output interval.  ULTRA fails BA12R if any of
the following occurs:

- any GSL failure, nonfinite state, currentness/frozen refusal, or the existing
  100000-step per-checkpoint guard fires
  (`CompactStar/Physics/Rotochemical/ScaledRKF45.hpp:40-54`);
- rejection fraction in the final temporal decade exceeds `0.10`;
- the median accepted-step count per output interval in the last half of that
  decade exceeds twice the median in its first half;
- the 95th percentile in the last half exceeds four times the 95th percentile in
  the first half; or
- the minimum accepted step in the final temporal decade is below `1e-6` times the
  fixed output spacing.

These are diagnostic refusal thresholds, not altered integrator physics.  No
thread-level solver or scientific change is proposed.

## 7. Nominal tier and minimum ULTRA scope

If BA12R passes, **REFINED becomes the nominal candidate trajectory and ULTRA is
only its convergence witness**.  This choice is fixed before the future run.  It is
the least expensive tier whose error is directly bounded by a still-finer solution;
the witness and all comparison metrics remain retained evidence.

Only these new trajectories are required:

1. `CPL-P2-LINEAR-QSS-v1` source at ULTRA; and
2. its exact matched no-BNV control at ULTRA.

The control is included so source-minus-control ledger observables have a same-tier
witness even though the existing P2 control was bit-identical across BASELINE and
REFINED (`docs/validation/PHASE6A1_CONTROLLED_BNV_NUMERICAL_FAILURE.md:63-69`).
No ULTRA P0, P1 or reaction-free run is required because only the P2 source failed
the old comparison and every other source/control passed.  All existing baseline/
refined files and hashes remain immutable.  No card, drive, duration, grid,
partition, initial state, frozen boundary, or baseline/refined tolerance changes;
the four exact mathematical cards remain those at
`docs/validation/PHASE6A1_CONTROLLED_BNV_IMPLEMENTATION.md:112-145`.

## 8. Evidence reuse and focused reruns

### 8.1 Reused unchanged as historical evidence

The following completed results remain valid evidence of the failed implementation
and must not be rerun merely because a recovery branch exists:

- BA1-BA10b results and the complete M1-M21 detection record;
- P0/P1/P2 local/event oracles, R18 and R-a/R-b/R-c;
- the accepted static BA13 21-star sample/certificate evidence;
- every existing BASELINE/REFINED source and control trajectory and its hash;
- matched-control results;
- protected-path/baseline/EOS/data/literature entry hashes; and
- the exact BA12 and BA15 failure evidence.

Their achieved dispositions and the later unfinished gates are recorded at
`docs/validation/PHASE6A1_CONTROLLED_BNV_NUMERICAL_FAILURE.md:81-101`.  Reuse does
not convert a historical result into proof of a modified binary.

### 8.2 Reruns required by tangent relocation

After R0, the future implementation must run only the focused dependency-linked
checks before expensive trajectory work:

1. clean configure/build and the tangent/projection executable covering BA2-BA5;
2. explicit M2 (`k` substituted), M17 (stale tangent), and M18 (domain mismatch);
3. the nonzero-eta direct-energy/R18 fixture because actual potentials consume `t`;
4. BA10a/BA10b focused matched-control construction/RHS checks because the wrapper's
   tangent type and include path changed;
5. certificate parsing/refit/currentness against the retained 21-star raw samples,
   confirming identical tangent values and frozen budgets; rerun the expensive star
   solves only if the recovered certificate cannot authenticate the retained samples;
6. fresh Phase-5D governed producer/comparator regression; and
7. focused Phase-6 source/projection/direct/wrapper tests after Phase-5D identity is
   restored.

The relocation does not require rerunning unrelated mutation families, P0/P1 event
quadrature, or the other BASELINE/REFINED trajectories.  If any purportedly neutral
value is not bit-identical, the future task must stop and return to the owner rather
than widening a budget.

## 9. Ordered recovery execution

Owner acceptance must freeze this plan before any code or ULTRA result.

| order | future action | stop/pass condition |
|---|---|---|
| R-1 | accept a narrow ADR-0016 (or explicit ADR-0015 amendment) for Phase-6 BNV ownership of the adapter | required before structural mutation |
| R0 | relocate/rename the tangent adapter; remove only the two new `Analysis` files; update Phase-6 includes/build/tests | no mathematical or Phase-5 byte change |
| R1 | clean build plus focused BA2-BA5, M2/M17/M18, R18 and wrapper/control checks | exact/accepted semantic neutrality |
| R2 | fresh empty-scratch Phase-5D governed regression | exact governed artifact/provenance equality; no baseline mutation |
| R3 | focused Phase-6 regression and retained-certificate revalidation | affected tests pass; retained sample hashes/currentness valid |
| R4 | only now run predeclared P2 ULTRA source and matched control | exact ULTRA tolerances/grid/card; no retuning |
| R5 | evaluate all BA12R state, ledger, R18/R20, frozen and solver-trend criteria | all exact thresholds pass |
| R6 | only after R5, complete unfinished BA11/R20 acceptance, BA16, trajectory BA17, BA14 finalization and full BA15 suites | every original applicable gate passes |
| R7 | produce a candidate artifact only after all gates | remains nonphysical, nongoverned and owner-unratified |

This staging prevents another full campaign before the two blockers are cleared.
The original campaign's complete data-free suite took `3033.11 s`, its coupled
Phase-5D oracle took about `1712 s`, and the failed fresh Phase-5D regression took
`2445.96 s`; these timings are retained execution evidence, not acceptance criteria
(`/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase6a1-controlled-bnv-implementation/build/phase6a1-final-debug/phase6a1-final-datafree.log:106-155`,
`/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase6a1-controlled-bnv-implementation/build/phase6a1-final-auth/Testing/Temporary/LastTest.log:14439-14470`).

## 10. Future process-level parallelism

Process-level concurrency is permitted only with distinct output/scratch roots and
unchanged scientific inputs:

- after one clean build, independent focused executables for BA2-BA5, direct R18,
  and BA10 controls may run concurrently if their fixtures write separate roots;
- P2 ULTRA source and its matched control may run concurrently after R0-R3 pass;
- if static sequence points ever genuinely require regeneration, independent target-
  `B` stars may be distributed by sample, followed by one serialized certificate
  assembly/refit;
- independent Phase-5B and Phase-5C focused regressions may run concurrently with
  Phase-6 data-free focused tests when they do not share generated roots.

The fresh Phase-5D governed producer/comparator is serialized because it owns one
scratch artifact and comparison sequence
(`tests/rotochemical/phase5d1_controlled_evolution_regression.py:166-226`).
Certificate assembly, BA12R comparison, candidate serialization, hash manifests and
any full suite using shared CTest/output roots are also serialized.  No thread-level
change to a star solve, RHS, GSL, coefficient calculation, or ledger is allowed.

## 11. Stop conditions

Recovery stops and returns to the owner if any of the following occurs:

- owner acceptance/ADR authority for the structural relocation is absent;
- the relocated adapter cannot consume the unchanged Phase-5B derivative with the
  same `B_B`, values, errors, closure and currentness;
- any existing governed Phase-5 source, producer, comparator or baseline would need
  modification;
- fresh Phase-5D provenance or artifact bytes do not equal the governed baseline;
- BA2-BA5, M2, M17, M18, R18, BA10a or BA10b changes unexpectedly;
- retained BA13 samples/certificate cannot be authenticated without a separately
  approved rerun;
- the ULTRA tier refuses, reaches a roundoff/floor plateau, or violates a step-trend
  threshold;
- any BA12R state, ledger, R18/R20, matched-control or frozen-validity criterion fails;
- a run card, drive, duration, grid, tolerance or boundary would need post-result
  adjustment;
- the P2 run would exceed `abs(DeltaB)/B0<=1e-6`; or
- interpretation would require a physical BNV rate/model, A18, superfluidity,
  Regime-II, MixedStar, sliding background or variable `Z`.

No failure may be repaired by changing a Phase-5D provenance exclusion, comparator,
baseline, or by substituting `k`/raw `G_y` for `t`.

## 12. Explicit review answers

1. **Is BA15 scientific or architectural?** Architectural ownership/provenance; no
   existing governed scientific source byte changed.
2. **Is relocation outside `Analysis` sufficient?** Yes, contingent on exact fresh
   Phase-5D provenance equality and semantic-neutrality gates.
3. **Does relocation preserve governed `t`?** Yes by construction; the Phase-5B
   derivative owner and ADR-0015 formula remain unchanged.
4. **Must the Phase-5D producer/comparator change?** **NO.**
5. **Must the Phase-5D baseline change?** **NO.**
6. **Is BA12 a large physical discrepancy?** No: `Delta Tinf` is about `0.929 K` at
   `Tinf about 1.052462e8 K`, but physical smallness does not reverse numerical FAIL.
7. **Is old BA12 `<=1` guaranteed by RKF45?** **NO.** It compares two global
   numerical trajectories, while the controller regulates estimated local step error.
8. **Should it be relaxed to two?** **NO.** That would be post-result tuning.
9. **Exact ULTRA tier?** `rtol=1e-11`,
   `atol=(1e-16,1e-22,1e-22)`.
10. **Exact BA12R criteria?** Sections 6-6.2: refined-normalized state difference
    `<=1`; resolvable contraction `<=0.10`; floor-limited difference `<=10F`;
    ledger convergence `<=1e-9` of its declared scale; retained R18/R20/frozen
    budgets; exact step-trend refusals.
11. **Nominal candidate?** REFINED, with ULTRA only as retained witness.
12. **Cards needing ULTRA?** P2 source and its matched control only.
13. **Reusable old evidence?** Section 8.1; all raw evidence remains immutable.
14. **What reruns after relocation?** Section 8.2 focused dependency-linked gates,
    fresh Phase-5D regression, P2 ULTRA source/control, then unfinished gates.
15. **Can the full 77-test campaign wait?** Yes; it belongs only after R0-R5 pass.
16. **What can run concurrently?** The distinct-root jobs in section 10; shared
    producer/comparator/artifact operations remain serialized.
17. **Is ADR-0016 required?** **YES**, under the structural/architecture rule, unless
    the owner formally amends ADR-0015 to authorize the exact relocation instead.
18. **Does ADR-0015 physics change?** **NO.** Its equations, seam, energy ledger,
    scope and exclusions remain unchanged.

## 13. Owner decision requested

The requested owner decision is whether to accept this exact bounded recovery:

1. authorize the narrow Phase-6 BNV ownership relocation through ADR-0016 (preferred)
   or an explicit ADR-0015 amendment;
2. preserve the old BA12 result as FAIL;
3. accept the exact three-level BA12R hierarchy, equations and thresholds above;
4. fix REFINED as the nominal candidate and ULTRA as witness;
5. authorize only the P2 ULTRA source/control pair after Phase-5D provenance is
   restored; and
6. retain all physical-model and baseline exclusions.

No implementation follows automatically from acceptance.  The exact next action is
to return this recovery plan to the human owner for explicit acceptance.
