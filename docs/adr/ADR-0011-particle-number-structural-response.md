# ADR-0011: particle-number structural response and equilibrium-sequence ownership

**Status: ACCEPTED — 2026-09-05**

> **PHASE-5B STRUCTURAL-RESPONSE CONTRACT ACCEPTED —**
> **PRODUCTION IMPLEMENTATION BLOCKED ON UNIT-BOUNDARY RECONCILIATION.**

Acceptance of this structural contract is not implementation authorization. No Phase-5B
production implementation may begin until the separate GSL/Zaki unit-boundary prerequisite
in Decision Q4 has been accepted and completed. INV-09 remains **INTENDED BUT UNVERIFIED**;
INV-11 remains **UNRESOLVED**.

## Context and authority

The Phase-5B-0 preflight derives the complete species-number measure, sequence derivative,
whole-baryon constraint and domain-qualified FR2005 spin-down mapping. It also records the
source inventory, scratch controls, proposed numerical budgets and the demonstrated GSL/Zaki
background-convention mismatch
(`docs/validation/PHASE5B0_INV09_GLOBAL_RESPONSE_PREFLIGHT.md:1`).

ADR-0007 P11 deferred the equilibrium-sequence owner; ADR-0010 separated structural inputs
from the later chemical susceptibility layer. This owner decision settles the Phase-5B
structural contract and validation plan while retaining the governance fail-closed boundary
for ambiguous units (`GOVERNANCE.md:65`;
`docs/adr/ADR-0007-hartle-second-order-monopole-response.md:418`;
`docs/adr/ADR-0010-rotochemical-off-equilibrium-thermodynamic-contract.md:273`).

## Decision

### Q1 — Phase-5B validation star

The primary free-gas Phase-5B implementation and validation fixture is the human-ratified
Structure-1 midpoint, `rho_c = 1.10e15 g cm^-3`, with the already ratified corresponding
static Track-R star. It is an already validated, source-compatible static configuration and
avoids the unresolved source-qualified maximum-mass selection semantics. No free-gas numerical
`I_Omega` target is published. The midpoint is therefore not a published `I_Omega` benchmark
(`docs/validation/PHASE5B0_INV09_GLOBAL_RESPONSE_PREFLIGHT.md:132`).

The distinct FR2005 Figure-1 configuration, `P_c = 4.5e34 dyn cm^-2`, is secondary
source-comparison material only. It neither defines the implementation fixture nor becomes a
benchmark until its domain/core-radius and `Omega_K` conventions are explicitly authenticated
(`docs/validation/PHASE5B0_INV09_GLOBAL_RESPONSE_PREFLIGHT.md:132`).

### Q2 — count and driver domain

The fixed-total-baryon constraint uses the **whole-star** baryon number,
`N_B = sum_baryons b_i N_i`. For the Track-R npe-mu free-gas fixture,
`N_B = N_n + N_p` (`docs/validation/PHASE5B0_INV09_GLOBAL_RESPONSE_PREFLIGHT.md:191`).

The generic structural objects `ParticleNumbers`, `FixedCentralEnergyNumberResponse`,
`EquilibriumSequenceNumberDerivative` and `FixedBaryonNumberResponse` are domain-qualified.
They may represent the whole star or an explicitly declared fixed-isobar subdomain. INV-09
implementation and validation use whole-star free-gas counts and responses.

Whole-star `K_i` is not automatically FR2005 core `I_Omega,i`. For an explicit core or
fixed-isobar domain, production must implement the source-backed PN7/PN8 mapping including
the boundary-flux term. A core result requires an explicit domain boundary. No nuclear-density
cutoff, crust-core pressure, reaction threshold or unpublished FR2005 core radius may be
invented (`docs/validation/PHASE5B0_INV09_GLOBAL_RESPONSE_PREFLIGHT.md:392`).

Because no reproducible free-gas core cutoff is published, Phase 5B may close INV-09 using
whole-star structural `K_i`; PB14 may validate the exact source mapping, formula and domain
machinery. No numerical claim of reproducing FR2005 free-gas core `I_Omega` is authorized
unless a source-authenticated domain becomes available. Conditional free-gas equivalence
outside all composition-changing layers may be tested, but it cannot supply the missing
source boundary.

### Q3 — equilibrium-sequence derivative owner

Production `B_i` is obtained by **complete canonical-star differentiation**:

```
B_i = (partial N_i / partial epsilon_c)_(Omega=0).
```

Complete independently solved neighboring equilibrium stars use the same EOS and the same
domain/surface policy. Prefer `x = ln(rho_c/rho_ref)` as the numerical differentiation
coordinate for conditioning; use symmetric steps at multiple scales and require a higher-order
estimate such as a five-point centered derivative or an equivalently justified local polynomial.
Record achieved central states. Each star uses its own grid, metric, composition and surface;
no geometry from the central star is reused on neighboring stars. No inverse-mass fit is
allowed. Branch, clamp or source-domain crossings fail closed. The production result records
the differentiation variable and materializes the governed derivative with respect to
`epsilon_c` in named units (`docs/validation/PHASE5B0_INV09_GLOBAL_RESPONSE_PREFLIGHT.md:330`).

The regular homogeneous linearized solution is a mandatory **independent oracle** under PB7,
not the production owner. ADR-0011 creates no public homogeneous Phase-4 API. PB7 remains
blocked until the unit-convention prerequisite in Q4 is accepted and completed.

### Q4 — GSL/Zaki unit-convention mismatch

A separate governed unit-boundary reconciliation is required before Phase-5B production
implementation. Canonical TOV uses the existing GSL gravitational and solar-mass convention;
NStar/Profile use Zaki mass, energy and pressure conversions for stored geometry; carried
`nu'` originates from the GSL solve. The resulting stored background does not exactly satisfy
one conventional geometric TOV system. The raw homogeneous/complete-star `B_i` comparison
disagrees by approximately `9e-4` for the leading species, while a convention-matched scratch
diagnostic reduces the discrepancy to approximately `1e-6`
(`docs/validation/PHASE5B0_INV09_GLOBAL_RESPONSE_PREFLIGHT.md:457`).

The convention-matched scratch result is diagnostic only. It neither waives the mismatch nor
constitutes the accepted PB7 independent oracle. The mismatch must not be hidden inside a
looser `K_i` tolerance. ADR-0011 authorizes no modification to TOVSolver, NStar or
RotationSolver.

Before Phase-5B production implementation, a separate unit-boundary task must:

1. identify one authoritative internal convention across the
   TOV -> NStar/Profile -> RotationSolver boundary;
2. determine the minimum correction;
3. preserve explicit public and unit semantics;
4. independently revalidate affected static and Hartle quantities;
5. recheck the qualified Structure-1 result; and
6. demonstrate the homogeneous/complete-star sequence comparison under the accepted single
   convention.

That task may propose its own ADR if needed. This ADR does not resolve the mismatch, select
the correction or invoke the pre-baseline correctness exception (`GOVERNANCE.md:88`).

## Accepted formula contract

For a declared domain `D`, the nonrotating count is

```
N_i^D = C integral_D w n_i dr,
n_i = Y_i n_B,
w = 4 pi r^2 / sqrt(1-2m/r),
C = 1e54
```

where `C` converts the numerical product `fm^-3 * km^3` to a count. There is no lapse factor
(`docs/validation/PHASE5B0_INV09_GLOBAL_RESPONSE_PREFLIGHT.md:191`).

With `q = Omega_geom^2`, the fixed-central-energy response is
`A_i^D = (partial N_i^D / partial q)_(epsilon_c)`. For a single moving outer isobar,

```
A_i = C { - integral_[0,R) w xihat dn_i
          + integral_0^R w n_i
              [ mhat/(r-2m) + r^2 exp(-2nu) s^2/3 ] dr
          + w(R) n_i(R-) xihat(R) }.
```

A shell domain subtracts the corresponding inner-boundary term. The density measure is
`dn_i`, not `dY_i` and not a nodal `dn_i/dp` column. No steepness detector is allowed.
Continuous neutron and muon onsets create no atom; true declared jumps use exact Stieltjes
atoms; the terminal term is counted exactly once
(`docs/validation/PHASE5B0_INV09_GLOBAL_RESPONSE_PREFLIGHT.md:223`).

The sequence derivative is

```
B_i = (partial N_i / partial epsilon_c)_(q=0),
```

with the Q3 production owner. Let `A_B = sum_baryons b_i A_i` and
`B_B = sum_baryons b_i B_i`. The fixed-baryon reduction is

```
(d epsilon_c/dq)_(N_B) = -A_B/B_B,
K_i = A_i - B_i A_B/B_B,
sum_baryons b_i K_i = 0,
sum_i charge_i K_i = 0.
```

The baryon and charge identities must hold within the propagated numerical budget. A
near-zero `B_B` fallback is forbidden (`docs/validation/PHASE5B0_INV09_GLOBAL_RESPONSE_PREFLIGHT.md:374`).

`I_Omega` is domain-qualified. At the whole-star `P=0` boundary,
`I_geom,i^whole = K_i^whole` only under the derived boundary assumptions. For a fixed-isobar
core,

```
I_geom,i^core = K_i^(core) - Y_i,c Q_B(P_core),
```

with the full PN7/PN8 domain and boundary-flux semantics
(`docs/validation/PHASE5B0_INV09_GLOBAL_RESPONSE_PREFLIGHT.md:392`).

## Omega normalization and units

The explicit normalization is

```
q = Omega_geom^2,
Omega_geom = Omega_phys / c.
```

The existing AngularVelocity conversion remains the owner. `A_i`, `K_i` and `I_geom,i` have
units count km^2. The physical source coefficient is

```
I_phys,i = I_geom,i / c^2                         [count s^2],
dot N_i^eq = 2 I_phys,i Omega_phys dot(Omega_phys).
```

No implicit `c^2` conversion is allowed
(`docs/validation/PHASE5B0_INV09_GLOBAL_RESPONSE_PREFLIGHT.md:429`).

## d epsilon/dP disposition

The Structure-1 approximately `0.410%` threshold-localized imported `d epsilon/dP` error is
not a blocker to the accepted measure plus complete-star path at the midpoint fixture.
First-order Hartle does not consume the column; the accepted monopole interior EOS source is
measure-based; fixed-`epsilon_c` `A_i` is measure-based in `dn_i`; production `B_i` uses complete
independently solved stars; and the threshold-region derivative mutation left the relevant
governed fields and scratch `A/K` unchanged
(`docs/validation/PHASE5B0_INV09_GLOBAL_RESPONSE_PREFLIGHT.md:457`).

Central pointwise derivative requirements, EOS value/table convergence and PB12 dependency
testing remain mandatory. ADR-0011 authorizes no `d epsilon/dP` correction.

## Validation status and implementation boundary

PB1–PB14 in the preflight are the **accepted validation plan**, not achieved validation.
Their numerical budgets remain proposed implementation acceptance targets except where Q4
prevents execution. PB7 is explicitly **BLOCKED** until the unit-boundary reconciliation is
accepted and completed (`docs/validation/PHASE5B0_INV09_GLOBAL_RESPONSE_PREFLIGHT.md:673`).

INV-09 cannot close without implemented PN1–PN8, all PB1–PB14, independent review and human
ratification. It remains **INTENDED BUT UNVERIFIED**. INV-11 remains **UNRESOLVED**.

Structural `B_i` and `K_i` are not corrected R2006 chemical `B/Btilde`, chemical `Z` or `W`.
No legacy RotochemicalCache is activated. No Btilde, paper Z/W, evolved chemical state, weak
rate, evolution, realistic-EOS off-equilibrium extension or BNV is begun or authorized.

## Consequences

The Phase-5B particle-number structural-response contract is defined, including its fixture,
domain semantics, complete-star derivative owner, independent homogeneous oracle, number
measure, fixed-baryon reduction, source mapping and explicit angular-frequency conversion.
Production implementation remains blocked on the separate unit-boundary reconciliation.

**Exactly one recommended next action:** open a fresh unit-boundary reconciliation audit from
canonical master to identify and adjudicate one consistent GSL/Zaki geometric convention across
TOVSolver -> NStar/StarProfile -> RotationSolver, quantify the minimal repair, and define the
required static/Phase-4/Structure-1 revalidation before any Phase-5B production implementation.
Do not begin it automatically.


## Implementation-status addendum — Phase-5B-R candidate, 2026-09-06

The accepted decision and original prerequisite/validation status above are preserved as
history. After canonical ADR-0012 A1 integration at
`a43d02227bf53c3242d3212f81dd71963804f3aa`, the owner explicitly authorized this implementation
and resumption of its authenticated PB6 stop. PN1–PN8 are now implemented in the candidate;
PB1–PB14 passed, including fresh PB7 and PB11–PB14. All 21 required mutations fired and all
eight endpoint micro-falsifiers passed. The representation repair uses one stored endpoint at
each continuous boundary; exact continuity, declared jumps and terminal atoms remain unchanged.

Focused tests passed 8/8, complete data-free tests 38/38, and the complete external-data suite
61/61 serially. Two independent generations of the structural candidate artifact are byte-identical.
The original stop evidence and all governed historical baselines remain unchanged. This addendum
changes no accepted formula, tolerance, domain boundary or scientific scope.

**PHASE-5B GLOBAL PARTICLE-NUMBER RESPONSE IMPLEMENTED AND CANDIDATE VALIDATED —
INDEPENDENT REVIEW REQUIRED.** Independent review and human ratification remain **PENDING**.
INV-09 remains **INTENDED BUT UNVERIFIED**; INV-11 remains **UNRESOLVED**. No Btilde, chemical
Z/W, eta evolution, weak-rate evolution, BNV or inferred free-gas core boundary is implemented.
No governed structural baseline, human ratification or canonical integration is conferred here.

Evidence: `docs/validation/PHASE5B_INV09_GLOBAL_RESPONSE_IMPLEMENTATION.md`,
`docs/validation/phase5b_resume_evidence.json`, and the explicitly nongoverned
`docs/validation/phase5b_structural_response_candidate.json`.

## Human-ratification status addendum — Phase-5B-RAT, 2026-09-06

**PHASE-5B IMPLEMENTATION AND INDEPENDENT REVIEW HUMAN-RATIFIED — 2026-09-06**

| Item | Disposition |
|---|---|
| `PHASE5B_SHA` | `fe08c94ed5ae525a9cb78331c5dd69d9c617d591` |
| Independent review | **PASS WITH EXPLICIT CLAIM-NARROWING CAVEATS** |
| Human owner | **RATIFIED** |
| PN1-PN8 | **RATIFIED** |
| PB1-PB14 candidate package | **RATIFIED WITH EVIDENTIARY CAVEATS** |
| Candidate structural artifact | **RATIFIED FOR GOVERNED INTEGRATION** |
| Governed structural baseline | **NOT YET INSTALLED** |
| Canonical integration | **NOT YET PERFORMED** |
| INV-09 | **INTENDED BUT UNVERIFIED** |
| INV-11 | **UNRESOLVED** |

The accepted decision and the historical candidate status above remain unchanged as records of
their respective decisions. Independent Opus scientific/numerical review of
`fe08c94ed5ae525a9cb78331c5dd69d9c617d591` is now **COMPLETE: PASS WITH EXPLICIT
CLAIM-NARROWING CAVEATS**, and the human owner has ratified that exact implementation.

**PHASE-5B STRUCTURAL RESPONSE HUMAN-RATIFIED — GOVERNED INTEGRATION REQUIRED BEFORE
INV-09 CLOSURE.** PN1-PN8 are human-ratified for the reviewed ordinary-`NStar` / Track-R scope.
PB1-PB14 are human-ratified as the candidate validation package, subject to the evidentiary
caveats in the ratification record. The deterministic candidate artifact is ratified as the
authorized content for a later governed installation. This addendum installs no governed
baseline and performs no canonical integration.

Accordingly, governed integration is **NEXT**; INV-09 remains **INTENDED BUT UNVERIFIED pending
integration**; INV-11 remains **UNRESOLVED**. Source-qualified `M_max` selection semantics and
the separately recorded `MixedStar` unit-boundary debt remain unresolved. No corrected R2006
`Btilde`, chemical `Z`/`W`, evolved `eta`, weak-rate or rotochemical evolution, or BNV work is
authorized by this ratification.

The controlling disposition, reviewed numerical evidence, PB6 shared-endpoint repair,
measure-complete `dn_i` representation, caveats, domain limits, and exact next action are in
`docs/validation/PHASE5B_INV09_GLOBAL_RESPONSE_RATIFICATION.md:1`.

## Integration-status addendum — Phase-5B-I, 2026-09-06

**PHASE-5B STRUCTURAL RESPONSE — GOVERNED INTEGRATION COMPLETE.**

| Item | Integrated disposition |
|---|---|
| Implementation | `fe08c94ed5ae525a9cb78331c5dd69d9c617d591` |
| Human ratification | `d67ac9be27677aa038a9ac7cdf5f262a32a008ef` |
| Governed regression artifact | `tests/baselines/phase5b_structural_response.json` |
| Artifact SHA-256 | `7588f0e9cd62f5b6be48bb725e4d0ba6e47b64f50d843117e70c17a78beeb5fa` |
| Post-installation validation | **PASS — focused 9/9; data-free 39/39; complete 62/62** |
| INV-09 | **VERIFIED / RESOLVED** |
| INV-11 | **UNRESOLVED** |

INV-09 closure is for the structural particle-number response only: ordinary `NStar`, the
slow-rotation domain-qualified PN1-PN8 contract, whole-star fixed-baryon validation, the Track-R
free-gas fixture, and explicit fixed-isobar mapping machinery. It does not close the later
chemical-state/evolution invariant.

All human-ratified qualifications remain controlling. In particular, PB9's baryon identity is
tautological; mutation totals mix production-discriminating and local algebraic/null checks;
PB9/PB10 budgets are conservative; `B_B conditioning = 1` is a cancellation ratio; PB8 can be
floor-limited; generic jumps require adapter metadata while ordinary Track-R onsets are
continuous; PB13 is an enclosure; EOS/tail authority is adapter-owned; and PB6 is a
partition/refinement validation rather than an independent physics oracle. No authenticated
free-gas core numerical `I_Omega` benchmark or source-qualified `M_max` is claimed.

The installation, fresh-generation authentication, exact hashes, suite results, baseline
immutability, and retained caveats are recorded in
`docs/validation/PHASE5B_INV09_GLOBAL_RESPONSE_INTEGRATION.md`.

## Owner-ratified uncertainty-semantics clarification — Phase-5C-1R-UQ, 2026-09-07

**ADR-0011 remains ACCEPTED. This is a semantic clarification, not a reopening of the
decision or INV-09.** The Phase-5B implementation record already distinguishes the values
reported by `Errors()` from the separately measured EOS/profile/refinement/oracle evidence:
the former include the declared computation's quadrature/roundoff, tail, stencil, and propagated
ingredient estimates, while PB1/PB6/PB7/PB12 characterize effects outside that universal claim
(`docs/validation/PHASE5B_INV09_GLOBAL_RESPONSE_IMPLEMENTATION.md:114`).

For the declared Phase-5B discrete representation and computation:

- `NumberResult::Errors()`, including the governed `K_error`, is the propagated
  **`numerical_error`**. It is not a claim that the exact continuum `K_i` lies within
  `K_error_i`, and it is not a universal **`certified_bound`** on background or reconstruction
  effects.
- **`certified_bound`** is reserved for a mathematically demonstrated enclosure under stated
  hypotheses, such as an individually established analytic surface/tail or refusal-window
  majorant. The term is not applied to `K_i`, `I_phys,i`, or downstream `W` merely because an
  error field exists.
- A downstream consumer that requires stability against reconstruction, refinement, or
  independent-route variation must carry that separately measured evidence as a
  **`validation_envelope`**: a conservative predeclared empirical envelope, not a probability
  distribution, confidence interval, formal truncation remainder, or certified continuum bound.

Downstream Phase-5C consumes the whole-star physical response
`I_phys,i=K_i/c^2` under the existing unit owner (`docs/validation/PHASE5B_INV09_GLOBAL_RESPONSE_IMPLEMENTATION.md:53`).
The accepted downstream uncertainty requirement is therefore placed on the consumed `K/I`
direction. It does not require a complete componentwise deterministic interval for `A_i` and
`B_i` separately; such a requirement could discard correlated cancellation in
`K_i=A_i-B_i A_B/B_B`. `A_i` and `B_i` remain scientifically meaningful parts of the
Phase-5B construction and retain all of their validation evidence.

This clarification changes no Phase-5B central value, formula, implementation, test, tolerance,
baseline byte, or domain. Phase-5B INV-09 closure remains **VERIFIED / RESOLVED**, with all nine
review qualifications retained (`docs/validation/PHASE5B_INV09_GLOBAL_RESPONSE_INTEGRATION.md:99`,
`docs/validation/PHASE5B_INV09_GLOBAL_RESPONSE_INTEGRATION.md:127`). INV-11 remains
**UNRESOLVED**.

## Governed-regression compiler-portability addendum — 2026-09-10

**ADR-0011 remains ACCEPTED and integrated; INV-09 is not reopened.** The human owner ratifies
the exact field `provenance.build.compiler` as mandatory execution-environment provenance rather
than a cross-compiler scientific-identity equality field. Every fresh Phase-5B artifact must record
its actual, nonempty compiler string. That value may not be normalized, fabricated, removed, or
copied from the historical baseline.

When compiler strings are identical, the governed regression continues to require complete
raw-byte identity. When they differ, the comparator may exempt only
`provenance.build.compiler`, must retain and report both strings, and must require exact parsed
equality everywhere else. No other `provenance.build.*` field is exempt. Any numerical,
scientific, source, EOS, architecture, configuration, or other provenance difference remains a
regression failure; no tolerance is widened.

The historical Phase-5B governed baseline remains byte-identical. The authenticated
Apple-LLVM-17 baseline and Apple-LLVM-21 regeneration agree in every equality-bearing field and
differ only at the compiler path. This is positive portability evidence, not proof of compiler
independence, and it changes no Phase-5B central value, error, formula, domain, physics, or INV-09
status. The detailed decision and comparator controls are recorded in
`docs/validation/PHASE5B_GOVERNED_REGRESSION_PORTABILITY_RATIFICATION.md`.
