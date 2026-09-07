# Phase-5C-0RAT — corrected chemical-coefficient contract ratification

> **PHASE-5C CORRECTED CHEMICAL-COEFFICIENT CONTRACT — HUMAN-RATIFIED;
> PRODUCTION IMPLEMENTATION NOT YET AUTHORIZED**

**Date:** 2026-09-06
**Change class:** scientific-semantic and structural/architecture contract revision;
documentation/governance only.
**Canonical master:** `49ab2b8c2881b6ef7b9309307d18cea51d557f72`
**PHASE5C0_SHA:** `54ec7abac38fa0a32c5fb3a82e424b496361966a`
**Branch:** `analysis/phase5c-corrected-chemical-coefficients-preflight`
**Worktree:**
`/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase5c-chemical-preflight`
**Accepted decision:** `docs/adr/ADR-0013-corrected-rotochemical-chemical-coefficients.md`

This record ratifies the corrected coefficient contract after independent review. It contains
no production implementation, test, baseline, EOS/data, literature, rate, evolution, thermal,
or BNV change. The eventual commit SHA is reported by Git after this document is committed; the
document does not make a self-referential commit claim.

## 1. Authenticated authority and review identity

The ratification started from the exact history

```text
49ab2b8c2881b6ef7b9309307d18cea51d557f72
  ->
54ec7abac38fa0a32c5fb3a82e424b496361966a
```

At entry, local `master`, cached `origin/master`, and live remote `master` all equaled the
canonical SHA. Candidate HEAD, cached upstream, and live remote candidate all equaled
`PHASE5C0_SHA`. Candidate parent and merge-base equaled canonical master; the candidate was
one ahead, zero behind, clean, and changed only the proposed ADR and preflight record. The
detached read-only review worktree was
`/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase5c-preflight-review`. No review
commit exists.

Independent-review disposition:

> **PHASE-5C PREFLIGHT INDEPENDENT REVIEW PASS WITH NONBLOCKING FINDINGS —
> ADR-0013 READY FOR HUMAN RATIFICATION WITH EXPLICIT CAVEATS.**

Blocking findings: **NONE**. Material findings: four, all incorporated. Nonblocking
clarifications/plan-strengthening findings: nine, all incorporated.

Human-owner disposition:

> **RATIFIED WITH EXPLICIT REVIEW REVISIONS**

## 2. Owner-ratified scientific facts

The owner accepts that R2006 governs where it supersedes F2005; the electrostatic correction is
required in the intrinsic derivation and is already incorporated by the neutral-Hessian route;
the full corrected four-species response is singular in the charge/gauge direction; and a full
paper-`B` inverse or silent pseudoinverse is forbidden. Local charge neutrality is local, while
baryon conservation is applied globally only after stellar integration. The unreduced source
basis `G_y` is canonical, `Z` is derived after global baryon reduction, and `W` combines that
corrected `Z` with the governed Phase-5B physical structural driver.

Active species use explicit lower-dimensional physical branches, never padded singular
Hessians. Coefficient-level redshift semantics are fixed without closing INV-11. Track R closes
on A18 + delta-v + UIX*, not CMF. Production implementation remains separately gated.

## 3. Exact Q1–Q8 owner decisions

### Q1 — canonical global chemical number-response authority

`GlobalChemicalNumberResponse` owns the canonical reduced source-basis
`y=(N_n,N_e,N_mu)` matrix `G_y`, with units `count / MeV` and named row/column axes. It remains
the independent source-aligned authority because later non-fixed-baryon sources need its
unreduced physical information. An x-basis transformed view may exist only as a derived
diagnostic/accessor. Do not expose a full corrected 4x4 inverse, pseudoinverse, or bare class
named `B`.

### Q2 — global baryon-reduction ownership and supported-mode policy

The required order is local charge-neutral susceptibility, then GR integration into unreduced
`G_y`, then global baryon reduction, then chemical imbalance response. Global reduction is not
optional and no baryon constraint is imposed pointwise. Retain `G_y` as the canonical authority.
Derive `Q=D-h h^T/a`; do not store `Q` as a duplicate independent authority unless future
implementation evidence requires a separately governed change. Represent supported channel rank
explicitly. Retain a channel only when its global response is supported on positive measure and
its smallest eigenvalue exceeds its uncertainty enclosure; otherwise return a separately named
lower-rank channel result or refuse. This order was already fixed by ADR-0010 and R2006 footnote
4, not selected among competing scientific alternatives here.

### Q3 — canonical corrected chemical Z owner

`ChemicalImbalanceResponse` contains one immutable symmetric 2x2 matrix in named channel order
`BetaChannel::Npe`, `BetaChannel::NpMu`, with
`eta^infinity=-Z delta N_l` and units `MeV / count`. `Z_npe`, `Z_np`, and `Z_npmu` are read-only
views onto that single matrix and are never independently stored. Output-channel rows and
input-lepton columns should use distinct semantic types even though symmetry would hide a
transpose/index error numerically.

### Q4 — W ownership and structural coupling

`RotochemicalSpinDrive` is a separate immutable result depending on
`ChemicalImbalanceResponse` and the complete governed semantic `FixedBaryonNumberResponse`.
It owns `W=Z(I_phys,e,I_phys,mu)^T` in `MeV s^2` and retains complete chemical and structural
provenance. A copied `I` vector is insufficient. The initial whole-star implementation consumes
`WholeStarIPhysical()` only after currency validation. Runtime baseline-JSON coupling is
forbidden. Spin down is one specific source adapter; the universal chemical source must not be
hard-coded as spin down.

### Q5 — domains, active branches, thresholds, and interfaces

The physical embeddings are npemu 3D, npe 2D, pe 1D, vacuum 0D/value-only, plus explicit
value-only `MuonThresholdEvaluation` and `NeutronThresholdEvaluation`. For an active branch,
`C_y=E H_active^-1 E^T` with authenticated branch conjugates. No padded Hessian, fabricated
absent-species response, density floor, hidden threshold extrapolation, or Hessian inversion on
threshold objects is permitted. A continuous onset has a one-sided susceptibility tending to
zero and no density-jump atom. Every finite provider-refusal window crossing support is bounded
or handled through a validated limit adapter; otherwise refuse. A genuine first-order phase
transition requires authenticated interface/phase metadata and a chemical interface-motion law.
The Phase-5B structural jump formula is not an automatic chemical-interface response.

### Q6 — concrete correction-sensitive Track-R benchmark instrument

A18 + delta-v + UIX* remains the mandatory realistic Track-R closing model, as already fixed by
ADR-0010 Q1. CMF, BPAL, and free gas cannot substitute. The preferred coefficient-layer
benchmark is R2006 Figure 1 using authenticated matching A18 authority; authenticated author
arrays are preferred if available, and governed figure extraction is allowed if not. A later
transient comparison may supplement but not replace the coefficient benchmark. Quasi-steady
temperature alone is insufficient. Realistic closure remains blocked until matching A18
authority exists; the owner ratifies the benchmark instrument here, not a new EOS selection.

### Q7 — coefficient-level redshift and source semantics

At coefficient level, `eta^infinity=e^nu eta_local`, with
`eta_npe=mu_n-mu_p-mu_e` and `eta_npmu=mu_n-mu_p-mu_mu`. Global `G_y` contains exactly one
`e^-nu` factor multiplying proper volume. `Z` acts on redshifted imbalances. `W` has units
`MeV s^2`, and for frozen coefficients the source sign is
`dot eta^infinity=-Z R+2 W Omega Omega_dot`. This records source semantics only, not a secular
evolution implementation. ADR-0013 partially resolves INV-11(a) only for coefficient objects;
INV-11 remains UNRESOLVED.

### Q8 — provenance, lifetime, and refusal

Future chemical results are immutable, dependency-complete, fail-closed scientific values.
Global provenance retains profile identity/version, metric/redshift normalization, provider
identity/revision and exact bytes where applicable, component constants, equilibrium anchor,
active branch map, basis/orientation, domain/reservoir, interface/tail policy, realized radial
partition, quadrature rule, onset splits, achieved node count, method/version, structural-zero
register, eigenspectrum/conditioning, error budgets, and R2006 convention. `W` additionally
retains all Phase-5B structural provenance. A changed dependency refuses before scientific
access; no lazy stale science or label-only reconstruction. Raw pointers without guaranteed
lifetime or a validated lifetime token are insufficient. No universal condition-number threshold
is ratified: a smallest supported eigenvalue must exceed its error enclosure before inversion.

## 4. Exact ratified mathematics

```text
x = (n_B,n_e,n_mu)^T
g_x = (mu_n,-eta_npe,-eta_npmu)^T
H_x = partial g_x / partial x

y = (n_n,n_e,n_mu)^T
T = [[1,-1,-1],
     [0, 1, 0],
     [0, 0, 1]]
y = T x
H_x = T^T H_y T
C_x = H_x^-1
C_y = T C_x T^T = H_y^-1
```

Independent validation may use only

```text
C_projected = chi - chi q (q^T chi q)^-1 q^T chi.
```

The public neutral path is already corrected; a second projection is forbidden.

```text
G_y = 10^54 integral_D
      4 pi r^2 e^-nu C_y(r) / sqrt(1-2m/r) dr
                                            [count / MeV]

b = (1,1,1)^T
L = [[-1,-1],
     [ 1, 0],
     [ 0, 1]]
delta N_y = L delta N_l
eta^infinity = -L^T delta g_y^infinity
Z = L^T G_y^-1 L = Q^-1             [MeV / count]
Q = D - h h^T/a                      [count / MeV]
eta^infinity = -Z delta N_l

W = Z (I_phys,e,I_phys,mu)^T         [MeV s^2]
```

The `Q` equality is a global Schur reduction after integration. It is not a pointwise
constraint, and `Q` is not a second stored authority.

## 5. Domain and embedding contract

| Physical branch/result | Dimension | Embedding into `y=(n_n,n_e,n_mu)` |
|---|---:|---|
| npemu | 3D | `T` |
| npe | 2D | `[[1,-1],[0,1],[0,0]]` |
| pe | 1D | `[[0],[1],[0]]` |
| vacuum | 0D/value-only | none; no Hessian inversion |
| `MuonThresholdEvaluation` | value-only | none; no Hessian inversion |
| `NeutronThresholdEvaluation` | value-only | none; no Hessian inversion |

First-order chemical interfaces are never inferred from the structural interface response.
Absent support lowers the returned physical rank or causes refusal; it is not repaired by
padding, extrapolation, a density floor, or a pseudoinverse.

## 6. M1–M4 material findings

### M1 — ADR-0010 V1–V12 crosswalk

`GC1`–`GC14` extend and operationalize the accepted ADR-0010 validation ladder; they do not
supersede it.

| ADR-0010 gate | Phase-5C carry-forward/discharge |
|---|---|
| V1 units/rest-mass/index/sign | GC1 |
| V2 exact neutral reconstruction | GC2 |
| V3 beta equilibrium/threshold conditions | inherited Phase-5A provider gates + GC8 when global support is consumed |
| V4 analytic lepton checks | inherited Phase-5A Track-R validation, not re-credited as new Phase-5C validation |
| V5 analytic toy Hessian/susceptibility | GC3–GC6 as applicable |
| V6 Hessian symmetry/integrability | GC3 |
| V7 finite perturbation vs linear response | GC3 / inherited local-provider validation |
| V8 x/y equivalence | GC5 |
| V9 intrinsic projection vs neutral route; null/rank/proton identity | GC6 + GC7 |
| V10 rank/support/stability/active species | GC8 + conditioning portions of GC10 |
| V11 corrected global Z/source response | GC9 + GC10 + GC11 + GC13 |
| V12 published non-superfluid thermal benchmark | **not discharged by coefficient implementation; later evolution-layer gate** |

### M2 — Q2 and Q6 reframing

Q2 now ratifies ownership and supported-mode policy around an already-fixed global ordering; it
does not present local-before-global baryon reduction as an alternative. Q6 now ratifies the
concrete correction-sensitive benchmark instrument around the already-fixed A18 + delta-v + UIX*
closing identity; it does not reselect the EOS.

### M3 — sign-flipped lapse and proper-volume controls

`GC9` includes **M19 — sign-flipped coefficient lapse**:

```text
G_wrong = integral e^(+nu) C dV
```

instead of `integral e^(-nu) C dV`. An independent curved-GR exact-star oracle is the detector.
The review's approximately `-35%` to `-38%` relative error is detector-separation evidence,
never a tolerance.

`GC9` also records the **M20 companion proper-volume inversion control**: using
`sqrt(1-2m/r)` instead of `1/sqrt(1-2m/r)`. It may remain a GC9 subcase rather than a separately
scored required mutation, but must be exercised.

### M4 — free-gas correction-separation rung

The free-gas fixture must demonstrate numerical separation of old uncorrected F2005 response
from corrected R2006 response by more than numerical uncertainty. Review diagnostics were
approximately `+1.55%` (`Z_npe`), `+82.8%` (`Z_np`), `+3.75%` (`Z_npmu`), and `3.6%` matrix
Frobenius norm. These are diagnostics only, not literature targets or tolerances. Free gas can
detect that the correction machinery is active; it cannot replace the realistic A18 gate.

## 7. N1–N9 nonblocking findings

### N1 — quadrature diagnosis

The current approximately `1.9e-6` scratch difference is localized mainly to trapezoidal
quadrature on the profile partition. Two independent non-trapezoid routes agreed at approximately
`1.8e-8`, while each differed from the trapezoid diagnostic by approximately `1.9e-6`.
These are diagnostics only. Quadrature is the dominant characterized scratch error, and the
production implementation must select and validate an onset-aware quadrature policy before
acceptance budgets are set.

### N2 — tolerance methodology

Quadrature truncation is a first-class uncertainty. Record the realized radial partition, onset
split locations, quadrature rule, achieved node count, and structural-zero register. Exact-zero
entries use absolute budgets; free-gas `G_ne` and `G_nmu` must never be assigned relative error
against zero.

### N3 — provider refusal window

The reviewed fixture bounded `Delta G_nn/G_nn` by approximately `1.7e-14` and induced relative
Z effect by approximately `9.8e-18`. They are not universal tolerances. The ratified policy is:
every finite response-refusal window crossing physical support needs an explicit response-measure
bound or validated limit adapter; otherwise refuse. The review shows the policy is executable.

### N4 — chemical tail

For the reviewed Track-R fixture, actual `Delta G_ee=2.685e45 count/MeV` was enclosed by the
`3.682e45 count/MeV` static chemical-tail bound. `R_upper` agreed with the review's exact `P=0`
radius at approximately `2e-10` relative. The companion shell-mass correction is below
independently demonstrated background reproducibility and is not dominant. This tail result is
not the full coefficient error budget.

### N5 — Figure-1 mass selection

Future extraction should predeclare masses approximately in `1.0`–`1.2 M_sun`, subject to
authenticated extraction evidence, where curve identity and old/corrected separation are clearer.
The review found crossings/overlaps around `1.4`–`1.6 M_sun`; Figure 1 spans approximately
`1.0`–`2.0 M_sun`. No numerical curve values are frozen and visual estimates are not promoted
to source data.

### N6 — governed extraction protocol

Future digitization requires: axis-calibration closure with residuals included in uncertainty;
a curve-identity association rule and rejection of masses inside a predeclared ambiguity zone;
rasterization independence with separate resolutions and preferably separate renderers; and a
source-text round-trip check that extracted trends remain compatible with R2006's qualitative
correction scale. Independent extraction cannot mean two people using the same raster only. No
digitization is performed here.

### N7 — INV-11 partial resolution

Coefficient-level semantics are now governed: `eta^infinity=e^nu eta_local`, named beta
channels, one `e^-nu` in `G_y`, Z action on redshifted imbalance, and W units/sign/source.
Secular-evolution ordering, representation, storage units/conversion boundary, stoichiometry,
net-rate sign ownership, changing coefficients/background, thermal/neutrino partition, and solver
coupling remain unresolved. INV-11 remains UNRESOLVED.

### N8 — value-only thresholds and lifetime safety

`MuonThresholdEvaluation` and `NeutronThresholdEvaluation` are value-only; no Hessian inversion
occurs on them. Future result objects must retain safe lifetime semantics for every dependency
needed by stale-input validation. Raw pointers without guaranteed lifetime or validated tokens
are insufficient. Existing Phase-5B provenance is not redesigned here.

### N9 — convergence/quadrature gate

`GC9a` covers the independent exact integral, lapse, proper volume, conversion, center, and tail.
Required `GC9b` separately covers radial refinement, quadrature-rule comparison, onset-aware
partition refinement, provider/table resolution where applicable, and surface/tail remainder.
One correct manufactured integral cannot make GC9 pass.

## 8. Revised validation ladder and benchmark blocker

The required coefficient ladder is:

1. exact analytic/toy corrected projection;
2. Track-R free-gas whole-star mechanics;
3. explicit free-gas old-vs-corrected separation; and
4. authenticated A18 + delta-v + UIX* Figure-1 or author-array comparison.

The realistic gate is blocked on authenticated matching A18 equilibrium/composition and exact
fit lineage, arbitrary-composition response and lepton conventions, crust/joins/core/phase and
chemical-interface treatment, matching source configuration, and corrected coefficient arrays or
governed Figure-1 extraction with quantified uncertainty. Free gas, CMF, BPAL, or a generic
APR-labelled barotrope cannot close that gap.

## 9. R2006 source dispositions

R2006 printed eq. (11) appears to omit the `e^(-Phi)` required by substitution into eq. (10) and
consistent with eq. (13). Classification is exactly **INFERRED PRINTED/SOURCE OMISSION**, not a
published erratum. No source quotation is rewritten and no Fable adjudication is needed.

R2006 Figure 2's right panel says `2.14 M_sun`; its caption says `2.13 M_sun`. Independent
review verified both. F2005 gives `2.14 M_sun` specific physical significance as its highest
described causal model, favoring the panel label, but this remains inference rather than an
authenticated correction. It is an unresolved source note, not an ADR blocker. No quantitative
Figure-2 benchmark may select either value without later adjudication/source authority. Figure 1
can carry coefficient-layer validation.

## 10. INV-11 boundary and downstream prohibitions

ADR-0013 partially resolves INV-11(a) only for coefficient-object semantics. INV-11 remains
UNRESOLVED for evolution, and no state/storage/rate/thermal/solver choice is made.

The following remain prohibited in this task and unauthorised by this ratification:

- production `G_y`, `Z`, or `W`;
- eta evolution or chemical-state storage;
- weak rates;
- heating/cooling or neutrino partition;
- time-dependent coefficients/backgrounds;
- superfluidity;
- BNV;
- tests, baselines, build changes, EOS/data, or literature changes; and
- canonical merge.

## 11. Final ratification status

| Item | Status after owner ratification |
|---|---|
| ADR-0013 | **ACCEPTED** |
| Phase-5C corrected chemical-coefficient contract | **HUMAN-RATIFIED** |
| Production `G_y/Z/W` | **NOT IMPLEMENTED** |
| Track-R realistic A18 closure | **BLOCKED ON AUTHORITY** |
| INV-09 | **VERIFIED / RESOLVED** in its governed structural scope |
| INV-11 | **UNRESOLVED** |
| BNV | **NOT BEGUN** |

No test suite is run because this is documentation-only and tests are outside the task's
authorized change and execution scope. Documentation integrity, exact allowlist, ancestry, local/
upstream/live equality, and `git diff --check` are the required closure checks. The branch is
pushed non-force and is not merged by this task.
