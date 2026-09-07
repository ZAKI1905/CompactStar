# Phase-5C-1R-UQ — structural uncertainty semantics ratification

> **END-TO-END K/I VALIDATION ENVELOPE IS SCIENTIFICALLY SUFFICIENT FOR
> PHASE-5C GENERIC/FREE-GAS CANDIDATE IMPLEMENTATION — REVISE N3 SEMANTICS
> BEFORE PRODUCTION OUTPUT.**

**Date:** 2026-09-07

**Disposition:** **INDEPENDENT UNCERTAINTY ADJUDICATION ACCEPTED BY THE HUMAN OWNER**

**Change class:** scientific-semantic governance clarification; documentation-only permanent
change.

**Implementation boundary:** no production code, test, baseline, EOS/data, literature, build,
or numerical-result generation; no merge.

This record freezes the inherited Phase-5B structural uncertainty semantics before any production
Phase-5C `G_y`, `Z`, or `W` result exists. It is review authority for the subsequent generic/free-
gas candidate task; it is not a retrospective claim that Phase-5B performed formal interval
numerical analysis.

## 1. Authenticated decision object

| Item | Authenticated value / disposition |
|---|---|
| Canonical master | `49ab2b8c2881b6ef7b9309307d18cea51d557f72` |
| Accepted Phase-5C preflight, `PHASE5C0_SHA` | `54ec7abac38fa0a32c5fb3a82e424b496361966a` |
| Accepted ADR-0013 ratification, `PHASE5C0_RATIFICATION_SHA` | `4780121f21010374da2eb50898e90a067795b6e5` |
| Reviewed numerical plan, `PHASE5C1P_SHA` | `09d1b3c935919ec85f3c629607797ae87c42d8e5` |
| Exact starting history | canonical -> preflight -> ADR ratification -> numerical plan |
| Planning topology at entry | 3 ahead / 0 behind canonical master; local/upstream/live equal |
| Documentation branch | `docs/phase5c-uncertainty-semantics-ratification` |
| Documentation worktree | `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase5c-uq-ratification` |
| Starting SHA | `09d1b3c935919ec85f3c629607797ae87c42d8e5` |
| Production branch preserved | `physics/phase5c-corrected-chemical-coefficients` |
| Production worktree preserved | `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase5c-chemical-coefficients` |

The canonical and preflight identities are the accepted ADR-0013 history
(`docs/adr/ADR-0013-corrected-rotochemical-chemical-coefficients.md:6`,
`docs/adr/ADR-0013-corrected-rotochemical-chemical-coefficients.md:7`). The reviewed numerical
plan records the accepted ratification/start, linear parents, and deliberately unmerged state
(`docs/validation/PHASE5C1_NUMERICAL_BUDGET_AND_IMPLEMENTATION_PLAN.md:13`). Git authentication for
this task established the exact linear parent chain and clean production worktree before any
documentation edit.

## 2. Independent and owner dispositions

The completed independent Phase-5B to Phase-5C uncertainty adjudication found the end-to-end K/I
validation envelope scientifically sufficient for generic/free-gas candidate implementation,
provided N3 semantics are revised before any production output. The human owner **ACCEPTS** that
adjudication without broadening Phase-5C-owned numerical-error policy.

ADR-0013 Q1-Q8, `G_y/Z/W` mathematics, ownership, active branches, source boundaries, and N1-N9
remain accepted unchanged (`docs/adr/ADR-0013-corrected-rotochemical-chemical-coefficients.md:208`,
`docs/adr/ADR-0013-corrected-rotochemical-chemical-coefficients.md:347`). This decision
supersedes only the earlier reading that every inherited structural-I contribution must be a
complete deterministic certified bound.

## 3. Why the two prior implementation starts stopped

Both earlier Phase-5C implementation starts stopped before production changes or output because
the plan's single `E_structural_I` row and the sentence “inability to certify a needed term is
refusal” could be read as requiring a complete deterministic interval for the inherited
structural input (`docs/validation/PHASE5C1_NUMERICAL_BUDGET_AND_IMPLEMENTATION_PLAN.md:557`,
`docs/validation/PHASE5C1_NUMERICAL_BUDGET_AND_IMPLEMENTATION_PLAN.md:611`). Phase-5B's governed
`Errors()` claim does not cover every EOS/profile reconstruction effect; those effects were
measured separately (`docs/validation/PHASE5B_INV09_GLOBAL_RESPONSE_IMPLEMENTATION.md:114`).

Stopping was therefore the correct fail-closed outcome: proceeding would either have silently
promoted `K_error` into an unsupported continuum certification or assigned an unmeasured inherited
term zero. Governance requires a stop on ambiguous scientific semantics rather than a guessed
interpretation (`GOVERNANCE.md:63`). The stops did not disprove Phase-5B or reopen INV-09; they
identified a downstream uncertainty-language mismatch that required owner adjudication.

## 4. Actual Phase-5B uncertainty claim and INV-09 boundary

Phase-5B INV-09 closure established a validated central structural response using PN1-PN8,
PB1-PB14, independent review, nonlinear closure, refinement evidence, and governed regression.
Its integration record retains nine claim-narrowing qualifications and closes only the ordinary-
`NStar` structural layer (`docs/validation/PHASE5B_INV09_GLOBAL_RESPONSE_INTEGRATION.md:99`,
`docs/validation/PHASE5B_INV09_GLOBAL_RESPONSE_INTEGRATION.md:127`).

It did **not** assert that exact continuum `K_i` lies within `K_error_i`. Existing
`NumberResult::Errors()` / Phase-5B `K_error` is the propagated numerical uncertainty of the
declared discrete computation. The governed baseline itself describes radial/EOS and independent-
oracle convergence as separately recorded evidence
(`tests/baselines/phase5b_structural_response.json:82`).

**INV-09 remains VERIFIED / RESOLVED, with all nine qualifications retained.** No central value,
formula, implementation, or baseline changes. INV-11 remains **UNRESOLVED**.

## 5. Ratified terminology

### 5.1 `numerical_error`

A propagated numerical uncertainty associated with the declared discrete representation and
computation. It may include arithmetic, quadrature on that representation, finite-difference or
stencil estimates, local solve residual, roundoff, and a mathematically proved remainder
explicitly included by the computation. Phase-5B `K_error` has this classification.

### 5.2 `certified_bound`

A mathematically demonstrated enclosure under explicitly stated hypotheses. The term is used only
where a proof exists, such as a specific analytic surface/tail or refusal-window majorant after its
hypotheses are established. Neither `I_phys` nor `W` receives a certified-bound claim here.

### 5.3 `validation_envelope`

A conservative, predeclared envelope assembled from measured discrepancies against independent
numerical/analytic routes and controlled refinement or representation variants. It is an empirical
validation quantity used to test candidate stability. It is not a probability distribution,
confidence interval, error bar, formal truncation remainder, or mathematically certified continuum
bound. These three terms are never interchangeable.

## 6. Downstream structural authority and consumed species

The Phase-5C structural input is

```text
I_phys,i = K_i/c^2.
```

The initial whole-star mapping and units are governed by ADR-0011
(`docs/validation/PHASE5B_INV09_GLOBAL_RESPONSE_IMPLEMENTATION.md:53`).
`RotochemicalSpinDrive` directly consumes only `I_phys,e` and `I_phys,mu`:

```text
W = Z (I_phys,e,I_phys,mu)^T.
```

(`docs/adr/ADR-0013-corrected-rotochemical-chemical-coefficients.md:182`). Neutron and proton
coefficients remain necessary to construct and validate the fixed-baryon response, but they are not
direct entries of W. The large neutron cancellation amplification is not itself W's uncertainty
amplification; this statement does not remove neutron/proton evidence.

`A_i` and `B_i` remain scientifically meaningful Phase-5B decomposition quantities. Because the
consumed downstream direction is `K/I`, Phase-5C does not require a complete separately certified
componentwise deterministic interval for A and B. Such separate propagation can discard correlated
cancellation in `K_i=A_i-B_i A_B/B_B`.

## 7. Two-track inherited structural uncertainty

For consumed species `i in {e,mu}`:

```text
E_I_numerical,i = K_error_i/c^2.
```

The existing governed unit owner is applied exactly once. `E_I_numerical` is inherited
`numerical_error`.

`V_I_validation` is the separately named, componentwise K/I `validation_envelope` assembled from
the predeclared evidence below. It is never exposed or described as an error bar, confidence
interval, certified/rigorous bound, or formal remainder.

## 8. Required validation-envelope evidence and classification

| Evidence class | Classification | Required disposition |
|---|---|---|
| Existing stored K error | `NUMERICAL_ERROR` | Convert once into `E_I_numerical`; do not duplicate in V |
| PB11 finite-q `q -> 0` direct K discrepancy | `VALIDATION_ENVELOPE_INGREDIENT` | Retain raw per-species discrepancy |
| PB11 fixed-baryon central shift / `Delta N_B ~ q^2` | `FALSIFIER` | Keep distinct from Richardson/direct K comparison |
| PB7 independent-background B discrepancy | `VALIDATION_ENVELOPE_INGREDIENT` | Transfer into consumed K with documented sensitivity/amplification |
| PB12 K-level EOS/table-resolution variation | `VALIDATION_ENVELOPE_INGREDIENT` | Retain direct K variation |
| Direct K-level radial ladder | `VALIDATION_ENVELOPE_INGREDIENT` | M1; new predeclared measurement |
| PB6 partition/knot variation | `VALIDATION_ENVELOPE_INGREDIENT` | Transfer remaining A-level evidence to K with documented sensitivity |
| PB10 raw per-species direct K discrepancies | `VALIDATION_ENVELOPE_INGREDIENT` | M2; new durable recording |
| Separately reviewed independent K comparison | `FALSIFIER` | Corroborate without duplicate pathway counting |
| W algebra alone | `CONSISTENCY_ONLY` | Contract evidence, not candidate stability validation |

PB6, PB7, PB10, PB11, and PB12 have distinct existing evidentiary roles
(`docs/validation/PHASE5B_INV09_GLOBAL_RESPONSE_RATIFICATION.md:99`,
`docs/validation/PHASE5B_INV09_GLOBAL_RESPONSE_RATIFICATION.md:225`).

## 9. M1 — K-level radial ladder

No direct K-level radial-resolution ladder presently exists. PB6 varies radial resolution
primarily at A level, while PB12 varies EOS/table resolution at K level. Before candidate
acceptance, the production-validation task must measure K directly across a predeclared radial
ladder on the same governed fixture with unchanged physical semantics. This is implementation-time
validation, not new Phase-5B physics and not a prerequisite to begin candidate code. It is
**PREDECLARED / NOT YET RUN**.

## 10. M2 — PB10 raw per-species discrepancies

PB10 already validates independent reconstruction of every K and is stronger evidence than its
conservative budget inequality (`docs/validation/PHASE5B_INV09_GLOBAL_RESPONSE_RATIFICATION.md:103`).
It does not durably record every raw per-species discrepancy needed for the new envelope. The
implementation task must add that durable test-side reporting. It requires no Phase-5B production
source change and is **PREDECLARED / NOT YET RECORDED**.

## 11. M3 — PB6 knot-shift interpretation

Large ratios formed by dividing neutron/muon PB6 knot shifts by stored `A_error` do not show that a
certified bound was violated, because stored `A_error` does not claim that reconstruction class.
PB6 knot/profile variations are retained as raw `VALIDATION_ENVELOPE_INGREDIENT` representation
sensitivity. For p/e, the existing tail contribution can dominate `A_error`, explaining why the
same ratio may be below unity. No variation is erased or hidden.

## 12. Nonconvergent and floor-limited policy

1. Contracting ladder without demonstrated asymptotic order: use the finest measured discrepancy;
   do not Richardson extrapolate.
2. Nonmonotone or sign-alternating ladder consistent with a numerical floor: record the floor as an
   envelope ingredient; do not invent an order.
3. Truly noncontracting or unbounded behavior without a stable envelope: refuse candidate
   acceptance.
4. Retain every raw variant, including inconvenient variants, in durable evidence.
5. Fit no arbitrary safety factor after seeing candidate W.

The existing plan already forbids inferred Richardson order from nonmonotone differences
(`docs/validation/PHASE5C1_NUMERICAL_BUDGET_AND_IMPLEMENTATION_PLAN.md:93`).

## 13. Double-counting policy

- Do not add Phase-5B K-error/tail components already included in K error again as envelope terms.
- PB11 direct-K/Richardson and fixed-baryon closure constrain different aspects; retain both with
  dependence stated.
- PB7 B discrepancy and central-shift/B_B effects may be correlated. Without proved covariance or
  cancellation, combine conservatively and state the possible overlap.
- Use a conservative maximum/envelope for PB6 radial/knot variants in one representation class;
  do not blindly sum duplicates.
- Independent-review discrepancies are falsifiers/corroboration; one numerical path is not counted
  twice under different labels.

## 14. W numerical error and structural validation envelope

The numerical track remains

```text
E_W_numerical
  <= |Z| E_I_numerical
   + E_Z |I|
   + E_Z E_I_numerical
   + E_W_arithmetic.
```

This is `numerical_error`; validation-envelope quantities are excluded. The separate inherited
structural track is

```text
V_W_validation <= |Z| V_I_validation.
```

This is `validation_envelope`, not numerical error, certified bound, or confidence interval. A
future Z validation envelope requires separate governance before it contributes.

## 15. Dual predeclared acceptance and GC12

Every required W component must satisfy both

```text
E_W_numerical <= G_W_numerical
V_W_validation <= G_W_validation.
```

Both goals are componentwise, absolute, derived from pre-production authority, predeclared before
the first production G/Z/W output, and immutable afterward. Either failure causes refusal without
silent relaxation. Structural failure uses a distinct class such as
`StructuralValidationEnvelopeUnmet`, separate from `AccuracyGoalUnmet`. **No numerical or
validation-envelope W goal is selected in this task.**

GC12 now validates both W numerical accuracy and inherited structural validation stability. Its
future inputs include `E_I_numerical`, `V_I_validation`, the new K-level radial ladder, and raw
PB10 per-species discrepancies. W=ZI algebra alone remains `CONTRACT / CONSISTENCY_ONLY`; the
reviewed plan already identifies the algebraic fixture and future real-owner requirement
(`docs/validation/PHASE5C1_NUMERICAL_BUDGET_AND_IMPLEMENTATION_PLAN.md:291`).

## 16. Phase-5C-owned numerics remain fail-closed

The two-track clarification applies only to inherited Phase-5B structural input. Every required
Phase-5C-owned numerical uncertainty term for G, Q, and Z remains measured/bounded under the
accepted fail-closed policy. An unmeasurable required owned numerical term causes refusal. The
structural envelope is not a general escape from numerical error accounting.

## 17. Governed-baseline semantics

A later governed Phase-5C regression artifact may be installed after candidate acceptance under
this two-track model. It certifies that accepted deterministic candidate bytes are reproducibly
regenerated under the ratified contract. It does not certify a formal interval-enclosed continuum
solution. Phase-5B likewise distinguishes byte-identical production from independent scientific
validation (`docs/validation/PHASE5B_INV09_GLOBAL_RESPONSE_INTEGRATION.md:44`).

## 18. Realistic A18 and publication scope

The method separating numerical error from validation envelope transfers to later realistic Track
R. The numerical free-gas envelope does not transfer to A18. A18 requires its own adapter
authority, EOS/table-resolution evidence, phase/interface evidence, radial-resolution evidence,
correction-sensitive source benchmark, and validation envelope. GC13 remains
**SOURCE-LIMITED / BLOCKED ON A18 AUTHORITY**
(`docs/validation/PHASE5C1_NUMERICAL_BUDGET_AND_IMPLEMENTATION_PLAN.md:615`).

Publication-level claims may report a clearly labelled validation envelope and numerical-error
ledger. They must not call a validation envelope a certified bound, confidence interval, or formal
truncation error without separately establishing those semantics.

## 19. Final scientific and scope status

| Item | Status after ratification |
|---|---|
| ADR-0011 | **ACCEPTED; uncertainty semantics clarified** |
| ADR-0013 | **ACCEPTED; two-track inherited structural uncertainty ratified** |
| Phase-5C numerical plan | **REVIEWED / REVISED UNDER INDEPENDENT UQ ADJUDICATION** |
| Production G/Z/W | **NOT IMPLEMENTED** |
| K-level radial ladder | **PREDECLARED / NOT YET RUN** |
| PB10 per-species envelope values | **PREDECLARED / NOT YET RECORDED** |
| INV-09 | **VERIFIED / RESOLVED** |
| INV-11 | **UNRESOLVED** |
| GC13/A18 | **BLOCKED ON SOURCE AUTHORITY** |
| BNV | **NOT BEGUN** |

## 20. Pre-production remainder and explicit output statement

Before candidate acceptance, the implementation task must measure the direct K-level radial ladder
and raw PB10 per-species K discrepancies; record the other predeclared K/I evidence classes;
construct the componentwise `V_I_validation`; propagate `E_W_numerical` and `V_W_validation`; and
freeze both componentwise absolute W goals before the first production output. All retained N1-N9
requirements and Phase-5C-owned G/Q/Z numerical terms remain binding.

**NO PRODUCTION G/Z/W RESULT EXISTS AT RATIFICATION.**

No production source, test, test registration, baseline, EOS/data, literature, or build file is
modified by this ratification. The clean production branch remains at the reviewed plan commit and
is not fast-forwarded or merged here.
