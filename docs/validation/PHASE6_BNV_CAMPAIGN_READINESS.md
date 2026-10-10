# Phase-6 BNV campaign readiness — 2026-10-10

**NOT READY: authenticated dependency/input recovery and nonlinear/global numerical
qualification remain material prerequisites. No full campaign was launched.**

Work is on `physics/phase6-bnv-campaign-readiness`, created from accepted candidate
`0c3f9269c884161c9d89163660a47b141a93ace6`. Master remains
`9a5c1eca40758e4c39f9f80db7324e12ad23e505`. Ratification is separate commit
`1b19d54`; accuracy and experimental predeclaration is `d99ebb4`, preceding all
new numerical probes. Evidence: `evidence/phase6-campaign-readiness/environment.json:1`.

## Ratification, authority and exact scope

The owner ratified ADR-0017 Amendment 1 for the bounded breakpoint-aware
architecture, not global accuracy, a trajectory/baseline, or the full campaign.
The exact acceptance and owner-reported Fable B/PASS disposition are in
`PHASE6_BREAKPOINT_RATIFICATION.md:1`. No independent Fable report bytes were
available here; the older archived obs32 review is not substituted for that report.

Changes in this task are documentation/governance, diagnostic test code, and
generated evidence. Numerical-method policy is proposed only. The accepted
production source, Cstar behavior, protected Phase-5D physics, historical evidence,
dependency lock and baseline bytes are preserved. The local audit authenticates
all 140 archived evidence files, 11 baseline files, and the accepted bytes of all
33 historical protected paths. Thirty match their original Phase-5D pins; three
already differ in the accepted candidate: the two approved full-byte file-reader
changes and `tests/analysis/chemical_production_contract.cpp`. The audit does not
misreport these as 33 original-pin matches. It records both hashes for each path.
Governing authority: `GOVERNANCE.md:36`, ADR-0015/0016/0017 and the explicit owner
readiness instruction; predeclaration: `PHASE6_CAMPAIGN_ACCURACY_PREDECLARATION.md:1`.

Canonical is an ancestor, but the candidate includes 24 intervening commits.
The prior `770d768` certificate-temperature enforcement change in particular
predates the amendment work. Those commits require their own authority review
before any integration proposal. No merge/rebase or master update was performed.

## Reproducibility and recovered inputs

macOS lists only the internal disk; Scriptorium is not available to mount. The
qualified Debug/Release package prefix and original BNV build inputs point through
broken links to `/Volumes/Scriptorium/Research/CompactStar/space_recovery/20261002T001142Z/`.
The drive was requested through the task UI; no replacement packages were built
or installed. The fresh qualified Debug CMake configure fails during dependency
authentication, before tests (`evidence/phase6-campaign-readiness/configure.log:1`).

A no-symlink-following local search covered Documents, Downloads, .local and
Library/Caches. Its 128 matching files contain only nonmatching older Zaki
archives, no qualified CONFIND archive/config package, and identical copies of
the radial certificate. Further tree hashing recovered matching profile and
thermal inputs. Four of six inputs were copied, then rehashed, under
`/Users/keeper/Documents/CompactStar/data/phase6-campaign-readiness/recovered-inputs`.
The original links and source copies were not modified.

| Input | SHA-256 | Fresh status |
|---|---|---|
| radial certificate | `7fc892b50bfbcdd2c5d963e0ff866777f3628f40fc282e1ce5a0b0364200a453` | PASS, recovered |
| profile tree | `233114862a2ab6826151519114e4bcb74a01f6a8e739bac0502ee62884688d72` | PASS, recovered |
| thermal tree | `1cdb98507958b50fbd1e1eff3eb3bfa4f48c790fef86a8740855884fc43b9d94` | PASS, recovered |
| entry manifest | `5cfbf4b1d623b58fa2a94101631dce20c602950eec8eb514dc1fcd3c7a04620d` | PASS, copied |
| frozen BNV certificate.tsv | `9217e278eb0a4b06b193a2a02f4d59f285d74e4033427b87ccd51181c3f7beab` | UNQUALIFIED, unavailable |
| coefficients.tsv | `3efa060d99a7a92266679aadb31885abbed9c940404247f1d3009cafa94832d6` | UNQUALIFIED, unavailable |

Package requirements remain source-built Zaki 2.0.1 at `e263a6e180c5c417198e7778bd21fc9c0a32dc33`
and CONFIND 2.0 at `b0cbd510fd3fd0c772fa50499cd749287cb39e7b`; exact Debug/Release
archive/header/config hashes are retained in `cmake/dependency-lock.json:1` and
copied into the environment evidence, with actual hashes null where unavailable.
Current compiler is Apple clang 21.0.0 (`clang-2100.3.34.2`), arm64; standalone
probes use C++17, `-g`, no optimization flag, GSL from `/opt/local`, and no
Zaki/CONFIND link. Compiler, CMake, OS, GSL version, dylib/header/executable hashes,
linkage, command and Python identities are recorded, not inferred from filenames.

The required 13-test qualified-stack rerun is **UNQUALIFIED / NOT EXECUTED**:
chemical_production_contract, chemical_exact_oracles, chemical_curved_gc9,
particle_number_analytic, passive_checkpoint_output_contract, phase6a1_direct_energy,
phase6a1_thermal_ledger, heat_capacity_v1, cache_contract, cache_thermal_contract,
phase5d_component_tolerances, phase6_breakpoint_contract, phase6_currentness_bytes.
The accepted historical 13/13 remains preserved. A fresh standalone build of the
unchanged `phase6_breakpoint_contract` passes, but is explicitly a different build
identity and does not replace that required suite.

## Tail stability: bounded diagnostic evidence, not nonlinear qualification

From the exact Fehlberg tableau, independent rational arithmetic gives
R5(z)=1+z+z²/2+z³/6+z⁴/24+z⁵/120+z⁶/2080 and
R4(z)=1+z+z²/2+z³/6+z⁴/24+z⁵/104. The first negative-real stability radii are
3.6777066213 and 3.0200175440. These coefficients match the GSL implementation
([GSL rkf45 source](https://raw.githubusercontent.com/ampl/gsl/master/ode-initval2/rkf45.c));
polynomial matrix evaluation reproduces the retained direct-GSL matrices within
5.90e-13 absolute. Derivation and exact input authentication are executable in
`tests/bnv/readiness_stability.py:1`.

The chemical spectral decay scale is about 3.861e-11/s. Applying the predeclared
half-radius margin, the smaller embedded radius, twice the measured Jacobian
refinement change and a 1% screening margin nominates
**h_max=3.8723807757e10 s (~1227 years)**. This is more conservative than the review's
4.8e10 s starting scale. It is a candidate fixed ceiling for bounded experiments;
endpoint samples do not bound the Jacobian along a nonlinear campaign. Production
`Rkf45Configuration::maximum_step_s` already exists; its default is unchanged.

A separate diagnostic uses the frozen archived full 3x3 Jacobian, perturbations
of ±0.25 of each tier's state budget along each coordinate, the retained tail state,
and a 5e11 s interval. There are 72 cases: 3 tiers × 3 axes × 2 signs × 4 policies.
The reference is a 70-digit matrix exponential, checked at 100 digits, of the
exact binary64 matrix/initial state. This quantifies exponential arithmetic only;
it does not bound error in the Jacobian or replace the nonlinear RHS.

| Policy | Worst error / state budget | Worst rejection | Accepted steps/case |
|---|---:|---:|---:|
| unchanged adaptive RKF45 | 0.422287 | 25% | 5–6 |
| candidate ceiling | 0.00867250 | 0% | 13 |
| half ceiling | 0.000158918 | 0% | 26 |
| quarter ceiling | 0.0000529156 | 0% | 52 |

Reference arithmetic uncertainty is <9.3e-48 of a state budget. Capped probes
show progress and exact final times. Their microsecond runtimes are **not** BNV
runtime estimates; no thermal/chemical ledger is defined for this linear surrogate.
No stiffness replacement is justified by these results. Nonlinear source/control
comparisons, refreshed Jacobians and ledger scoring remain UNQUALIFIED.

The full scaled Jacobian is nonnormal. At the candidate ceiling its amplification
spectral radius is 0.8314, but the 2-norm is 10.26 for scales (1,1e-7,1e-7). The
exact frozen flow has norm 10.14, so most of this transient amplification is real
coupled response, not numerical instability. A spectral ceiling alone cannot
transport componentwise error. See `evidence/phase6-campaign-readiness/retained-summary.json:1`.

## Global accuracy and retained numerical margins

The committed predeclaration requires original-initial-state independent references,
explicit maximum-step refinement, quantified reference/roundoff error, separate
source/control budgets, and conservative source-minus-control uncertainty. It
retains R18 bookkeeping and R20's 2e-4 finite-interval closure threshold, now
requiring explicit main/reconstruction/event/quadrature uncertainty. No cancellation
or identical-step agreement supplies missing accuracy evidence.

It also specifies an alternative verified error-transport recurrence through every
prefix, including full variational flow, nonlinear remainder, roundoff and event
error. Local embedded errors, sampled Jacobians and warm restarts cannot establish
that recurrence. No new global accuracy PASS is claimed. The reference/transport
machinery for a full campaign has not been implemented or qualified.

The following are authenticated **retained** bounded results, not new runs or
scores against the stronger global predeclaration (`retained-summary.json:1`):

| Window | worst state ratio | ledger ratio R/U combined | max abs(R20/N) |
|---|---:|---:|---:|
| warm source | 0.929636 | 0.00347480 | 1.28730e-5 |
| warm control | 0.276599 | 0.936051 | 1.24901e-5 |
| source from zero, bounded | 0.911509 | 0.00398206 | 5.25076e-7 |
| warm source tail | 0.345048 | 2.11035e-13 | 8.74729e-10 |

Source eta_mu ratios at BASELINE/REFINED/ULTRA are 0.809145/0.900860/0.929636;
eta_e ratios are 0.476907/0.594452/0.649666. Thus ULTRA mu chemistry has only 7.04%
remaining bounded margin. The control ledger has only 6.39% margin. The archived
ledger maximum combines REFINED/ULTRA; it cannot be decomposed by tier from this
summary, and independent BASELINE ledger accuracy was not scored by that harness.
Maximum retained R18/G_P is 1.319e-16. Maximum retained reconstruction/N is
2.164e-10 versus 5e-6. The worst local reconstruction/event allocation utilization
is 0.2686 of its allowed budget. These margins omit unknown full-prefix error.

Historical BA12/BA12R remain **FAIL**. New bounded old-form control gates also
fail (state c2/c3: 23/12; ledger c2/c3: 98/2); tail state c2 has two failures and
rejection fractions are 14.29%, 22.22%, 15.38%. Their separately accepted bounded
successor accuracy PASS remains a different result. Neither is relabeled under
the proposed global criterion. Original source/control ledger/R20 differences
were not qualified as a matched pair from a common full prefix.

## Breakpoint edge cases and schedule coverage

The accepted implementation is retained. New GSL-only probes exercise both
methods, both directions, exact knot starts and residual-strip starts at t=1e12 s.
Exact starts generate no events; a residual starting just across the knot creates
one finite micro-event, preserves state continuity and the time bracket, and
reaches the endpoint without a loop. Such residual events are not evidence of
new physical crossings. Suppressing them by snapping state or transferring a
shared root would require new error accounting, so no such change was made.
Specific stationary-knot and reversing-trial refusals are verified. All 14 edge
cases pass (`evidence/phase6-campaign-readiness/edge-probe.tsv:1`).

The archived harness checks passive-schedule invariance for **BASELINE** in each
window (`tests/bnv/breakpoint_qualification.cpp:56`); it does not establish a fresh
REFINED/ULTRA schedule-invariance result. Stored endpoints and adjacent boundary
reconstruction were tested at every tier in retained windows. The new ceiling
requires fresh all-tier nonlinear schedule, endpoint and uncertainty checks.

## Readiness matrix

PASS is limited to the explicitly named scope. FAIL denotes observed failed
criteria; UNQUALIFIED denotes missing evidence; DEFERRED denotes intentionally
out-of-scope action. “Retained PASS” is not a fresh rerun.

| Requirement | Status | Evidence or remaining work |
|---|---|---|
| Exact candidate, owner ratification, separate docs commit | PASS | exact SHA; no production delta |
| 140 archived identities / 11 baseline identities | PASS | fresh hash audit |
| 33 protected paths versus accepted candidate | PASS | 30 original pins + 3 accepted deltas |
| Recover qualified packages and all six inputs | UNQUALIFIED | four inputs recovered; two and packages missing |
| Required 13-test qualified-stack suite | UNQUALIFIED | configure blocked before execution |
| Standalone accepted breakpoint contract | PASS | fresh GSL-only build |
| Fresh exact-start/residual/reversal edge probes | PASS | 14 cases; fail-closed behavior retained |
| Existing nonlinear equilibrium-tail rejection/stability | FAIL | 14–22% rejection; amplification radius 8.704 |
| Capped frozen-linear perturbation diagnostics | PASS | diagnostic scope only, 72 cases |
| Capped nonlinear tail source and control | UNQUALIFIED | no qualified stack/input access |
| Main accumulated error / global reference uncertainty | UNQUALIFIED | criterion predeclared, no global proof/run |
| Matched source-minus-control accuracy and signal resolution | UNQUALIFIED | no matched common-prefix qualification |
| R18/R20 under new policy and global error budget | UNQUALIFIED | retained local evidence insufficient |
| Frozen BNV 21-star certificate currentness | UNQUALIFIED | exact table unavailable |
| Existing fixture depletion arithmetic | PASS | 5e-7 / 1e-6 = 0.5; not certificate validation |
| Historical BA12/BA12R | FAIL | preserved without reinterpretation |
| Global successor acceptance | UNQUALIFIED | no new global experiment |
| Physical BNV model/rate, Phase-5D changes, cluster | DEFERRED | excluded from this task |
| Full campaign and canonical integration | DEFERRED | require explicit owner authority |

| Role / tier | Retained local state accuracy | Retained schedule invariance | New policy nonlinear qualification | Global accuracy |
|---|---|---|---|---|
| source BASELINE | PASS, source/bounded/tail windows | PASS | UNQUALIFIED | UNQUALIFIED |
| source REFINED | PASS, source/bounded/tail windows | UNQUALIFIED | UNQUALIFIED | UNQUALIFIED |
| source ULTRA | PASS, source/bounded/tail windows | UNQUALIFIED | UNQUALIFIED | UNQUALIFIED |
| control BASELINE | PASS, warm control window | PASS | UNQUALIFIED | UNQUALIFIED |
| control REFINED | PASS, warm control window | UNQUALIFIED | UNQUALIFIED | UNQUALIFIED |
| control ULTRA | PASS, warm control window | UNQUALIFIED | UNQUALIFIED | UNQUALIFIED |

## Cost and exact next action

At this candidate ceiling a 500 kyr path would need at least 408 main steps;
accuracy and events can require more. The retained main solves are subsecond per
bounded window, but diagnostic/reconstruction/reference work dominates the full
qualification. The 1389.58 s final07 run scored 1212 observations. Scaling its
all-in work proxy (1.1465 s/observation) to six ×8193 observations gives **15.7 h**.
This is a planning scenario, not a measured campaign time, upper bound or launch
recommendation. Fresh reference trajectories and nested R20 grids can cost much
more; the full criterion has no qualified runtime bound. Recalibrate after stack
recovery with separately timed assembly, main, reconstruction, diagnostics and
references (`retained-summary.json:1`; run-card source `tests/bnv/campaign_fixture.hpp:90`).

**Next action:** reconnect the existing Scriptorium volume and authenticate its
identity and the exact locked packages, frozen certificate and coefficient table.
Then rerun the 13 affected regressions and the predeclared bounded nonlinear
source/control/tail experiments using 3.8723807757e10 s and its refinements,
including all-tier passive schedules, endpoints and uncertainty propagation.
Resolve any failure without relaxing the declared budgets. Obtain independent
review/owner acceptance of any qualified new policy and global assessment before
requesting separate full-campaign authorization. Review intervening ancestry
separately before proposing canonical integration. Readiness is not established;
no owner campaign authorization is recommended now.
