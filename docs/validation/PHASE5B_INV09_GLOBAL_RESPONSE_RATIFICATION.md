# Phase-5B INV-09 global particle-number structural-response ratification

> **PHASE-5B STRUCTURAL RESPONSE HUMAN-RATIFIED —**
> **GOVERNED INTEGRATION REQUIRED BEFORE INV-09 CLOSURE.**

Date: **2026-09-06**.

This is the durable human-owner ratification record for the Phase-5B structural candidate.
It is a documentation/governance decision. It changes no production code, test, tolerance,
build file, governed baseline, or candidate-artifact byte, and it performs no canonical merge.

## 1. Authenticated decision object

| Item | Authenticated value / disposition |
|---|---|
| Canonical master | `a43d02227bf53c3242d3212f81dd71963804f3aa` |
| Candidate branch | `physics/phase5b-global-particle-response` |
| Candidate worktree | `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase5b-global-response` |
| Ratified implementation, `PHASE5B_SHA` | `fe08c94ed5ae525a9cb78331c5dd69d9c617d591` |
| Candidate parent / merge base | canonical master above |
| Pre-ratification topology | candidate exactly 1 ahead / 0 behind master |
| Independent review | **PASS WITH EXPLICIT CLAIM-NARROWING CAVEATS** |
| Human owner | **RATIFIED** |
| PN1-PN8 | **HUMAN-RATIFIED for the reviewed scope** |
| PB1-PB14 candidate package | **HUMAN-RATIFIED WITH EVIDENTIARY CAVEATS** |
| Candidate structural artifact | **RATIFIED FOR GOVERNED INTEGRATION** |
| Governed structural baseline | **NOT YET INSTALLED** |
| Canonical integration | **NOT YET PERFORMED** |
| INV-09 | **INTENDED BUT UNVERIFIED pending governed integration** |
| INV-11 | **UNRESOLVED** |

The authenticated implementation, validation, provenance, suite, and preserved-history record is
`docs/validation/PHASE5B_INV09_GLOBAL_RESPONSE_IMPLEMENTATION.md:1`; the complete machine-readable
candidate evidence is `docs/validation/phase5b_resume_evidence.json:1`. The owner ratifies the
exact commit and no broader state.

## 2. Human-owner disposition and exact structural claim

The human owner ratifies `PHASE5B_SHA` as the accepted implementation of the ADR-0011 structural
particle-number response, subject to every qualification in this record.

CompactStar implements a **domain-qualified structural particle-number response for slowly
rotating ordinary neutron stars under ADR-0011**. For declared species and a declared domain it
provides

```text
N_i = nonrotating particle number,
A_i = (partial N_i / partial Omega_geom^2)_(epsilon_c),
B_i = (partial N_i / partial epsilon_c)_(Omega=0),
K_i = A_i - B_i A_B/B_B
```

where production `B_i` is obtained from complete, independently solved equilibrium stars and
`K_i` is the response along a fixed-total-baryon sequence. The result also provides the
domain-qualified `I_Omega` mapping and an explicit geometric/physical angular-frequency
conversion. The accepted formula and ownership contract is
`docs/adr/ADR-0011-particle-number-structural-response.md:119`.

The `A_i` displaced-composition contribution uses the measure-complete numerical representation

```text
- integral w xihat dn_i
```

including smooth variation, sharp continuous variation, declared true jumps, moving boundaries,
and terminal semantics. This claim is ratified for the implemented ordinary-`NStar` / Track-R
validation scope only.

The owner specifically accepts the PN1-PN8 implementation, PB1-PB14 candidate validation, PB6
shared-endpoint representation repair, measure-complete `dn_i` formulation, complete-star `B_i`
ownership, PB7 homogeneous/sensitivity oracle, fixed-baryon reduction, PB11 nonlinear finite-spin
closure, PB13 treatment for the validated Track-R tail adapter, domain-qualified PN7/PN8 mapping,
provenance/stale-input refusals, and deterministic candidate artifact. No broader claim is
accepted.

## 3. PN1-PN8 ratification summary

| Contract | Ratified content and boundary |
|---|---|
| PN1 | `N_i = 1e54 integral 4 pi r^2 n_i/sqrt(1-2m/r) dr`, with `n_i=Y_i n_B`, count units, and no lapse |
| PN2 | Fixed-`epsilon_c` `A_i` includes the signed `dn_i` measure, metric, velocity, and moving-boundary terms per `Omega_geom^2` |
| PN3 | `B_i` is the complete canonical-star sequence derivative with each neighbor owning its grid, metric, composition, and surface |
| PN4-PN6 | `d epsilon_c/dq=-A_B/B_B` and `K_i=A_i-B_i A_B/B_B`; near-zero `B_B` refuses without fallback |
| PN7-PN8 | Whole-star and explicit fixed-isobar core/shell mappings retain the required domain and boundary-flux terms |
| Units | `q=Omega_geom^2`, `Omega_geom=Omega_phys/c`, and `I_phys=I_geom/c^2` |

The detailed source traceability and implemented API mapping are recorded at
`docs/validation/PHASE5B_INV09_GLOBAL_RESPONSE_IMPLEMENTATION.md:45`.

## 4. PB1-PB14 ratification summary

| Gate | Human-ratified candidate result; evidentiary interpretation |
|---|---|
| PB1 | PASS: analytic, same-representation, and independent count checks; independent counts agree with the candidate at approximately `1e-15` |
| PB2 | PASS: lapse cancellation, angular average, `l=2` orthogonality, and integration-by-parts checks; some are local algebraic/null checks rather than production mutations |
| PB3 | PASS: signed jumps, continuous-ramp limit, terminal accounting, and exact endpoint micro-falsifiers |
| PB4 | PASS: independent term-by-term `A_i` recomputation, including density-measure, metric, velocity, and surface terms, agrees at approximately `1e-15` |
| PB5 | PASS: nonlinear current limit and explicit `q`/physical-unit normalization |
| PB6 | PASS: partition, radial-resolution, and threshold-refinement sensitivity after the ratified shared-endpoint repair; **not an independent physics oracle** |
| PB7 | PASS: the independent homogeneous/sensitivity oracle agrees by approximately `[9.94e-7, 8.43e-7, 8.03e-7, 1.21e-6]` for `[n,p,e,mu]`, below the predeclared `2e-4` criterion |
| PB8 | PASS as an error-controlled derivative study; no claim of clean uniform asymptotic convergence for every species |
| PB9 | PASS for the conditioned reduction and negative controls; the raw baryon identity is algebraically guaranteed and is **not independent physical validation** |
| PB10 | PASS: strongest evidence is independent reconstruction of all four `K_i`; the charge residual normalized to the charged assembly scale is approximately `2.0e-13` |
| PB11 | PASS: independent nonlinear finite-spin closure gives `Delta N_B proportional to q^2` and species quotient errors linear in `q` |
| PB12 | PASS: dependency audit, EOS-resolution ladder, and threshold-center probes; no hidden `d epsilon/dP` dependency in PN1-PN5 |
| PB13 | PASS for the validated Track-R adapter; the response bound is a conservative enclosure, not a precision estimate |
| PB14 | PASS: formula, sign, units, whole/core/shell domain behavior, and boundary-flux mapping; no core boundary is invented |

The candidate values and original gate records are at
`docs/validation/PHASE5B_INV09_GLOBAL_RESPONSE_IMPLEMENTATION.md:75` and
`docs/validation/phase5b_resume_evidence.json:1`. The qualifications below govern how those
passes may be cited.

## 5. Independently reviewed numerical evidence

The independent review verified the following without requiring a code correction, tolerance
change, or additional computation:

- `N_i`: independent counts agreed with the candidate at approximately `1e-15`.
- `A_i`: independent term-by-term recomputation agreed at approximately `1e-15`.
- `B_i`: an independent complete-star calculation using a different step ladder agreed at
  approximately `1e-7` to `1e-9` scale.
- PB7 is genuinely independent in differential variable, EOS evaluation, integrator, center
  expansion, and derivative owner. The species discrepancies were approximately
  `[9.94e-7, 8.43e-7, 8.03e-7, 1.21e-6]` against `2e-4`.
- PB11 is genuinely nonlinear and independent. Local powers for `Delta N_B` were approximately
  `1.995`, `2.003`, and `1.997`; species quotient errors were linear in `q`. Independent
  `q -> 0` extrapolation reproduced `K_i` at approximately
  `[5.3e-6, 5.2e-8, 9.5e-8, 2.9e-7]` relative for `[n,p,e,mu]`.
- The deterministic artifact was regenerated twice with SHA-256
  `7588f0e9cd62f5b6be48bb725e4d0ba6e47b64f50d843117e70c17a78beeb5fa` both times.
- The reviewed suite results were focused `8/8`, complete data-free `38/38`, and complete
  external-data `61/61`, all PASS.

## 6. Material claim-narrowing findings

### M1 — PB9 baryon identity is algebraic

The raw identity `sum_baryons b_i K_i=0` follows identically from the construction:

```text
sum b_i K_i
= A_B - B_B A_B/B_B
= 0
```

Therefore PB9's small raw baryon residual is principally floating-point cancellation roundoff,
not independent validation of fixed-baryon physics. The identity itself must not be cited as
evidence that the fixed-baryon reduction is physically correct. Independent support instead
comes from PB11 nonlinear finite-spin closure, PB10 independent `K_i` reconstruction, the
fixed-baryon-sign and wrong-baryon-species mutations, and fail-closed `B_B` conditioning. This
narrows PB9's evidentiary interpretation; it does not invalidate the formula or implementation.

### M2 — mutation/null-check campaign is heterogeneous

The reviewed campaign contains **34 fired mutation/null-check labels plus 8 endpoint
micro-falsifiers**. A subset directly discriminates production behavior; a subset is local
algebraic or toy-arithmetic null checking. Examples in the latter category include some checks
of the angular `1/3`, `l=2` algebra, lapse cancellation, count-unit factor, invented onset atom,
seed-squared contamination, `q` normalization/omitted division, and measure sign.

The raw count is not an evidence-strength metric and must not be described as though every label
independently validates production. The relevant physical properties are independently supported
by PB1, PB3, PB4, PB5, PB7, PB10, and PB11 as applicable.

### M3 — PB9/PB10 budgets are conservative

The propagated PB9/PB10 budgets are defensible in scale but substantially larger than the
achieved residuals. The PB9 residual is below its budget by a ratio of approximately `3e-10`,
and the PB10 charge residual by approximately `1e-6`; those inequalities alone therefore have
limited independent falsification power. For charge, the stronger achieved statement is a
normalized residual of approximately `2.0e-13` relative to the charged `K` assembly scale.
PB10's strongest evidence is its independent reconstruction of all four `K_i`, not merely the
budget inequality. No code or tolerance change is required or authorized.

## 7. Nonblocking technical and scope findings

### `B_B` conditioning metric

The reported `conditioning=1` is the cancellation ratio
`|sum b_i B_i| / sum |b_i B_i|`, not a general numerical condition number. It is exactly one
for this fixture because `B_n` and `B_p` are both positive. More informative values are the
approximately `7.5e-8` relative uncertainty in `B_B` and approximately `33` cancellation
amplification in `K_n`. The production near-zero-`B_B` refusal remains accepted and fail-closed.

### PB8 convergence

Neutron shows the expected high-/low-order step convergence and electron is broadly consistent.
Proton and muon enter finite-difference cancellation/noise floors. Production error estimates
correctly inflate to include adjacent-step differences, and PB7 independently checks the values.
PB8 therefore passes as an error-controlled derivative study, not because every species displays
a clean asymptotic ratio.

### Declared jumps and the ordinary profile adapter

The low-level measure machinery supports exact declared Stieltjes jump atoms. The current
ordinary profile adapter does not automatically discover and materialize first-order phase
transitions; such inputs fail closed. The validated Track-R neutron and muon onsets are
continuous, so this limitation does not affect this Phase-5B result. A future EOS with a genuine
first-order density/composition jump must supply authenticated jump metadata through an
appropriate adapter. Automatic jump discovery is not claimed.

### Exact continuity guard after PB6 repair

Canonical shared endpoints are now constructed once and shared exactly. The exact continuity
guard remains strict and was not weakened. For internally generated ordinary continuous segments
the corrected construction makes that guard structurally unreachable; it remains useful as a
structural invariant and an externally supplied partition guard. It is not described as the
primary production discontinuity detector.

### PB13 response enclosure

PB13 is accepted for the validated Track-R tail adapter. The count bound is relatively tight.
The response bound is intentionally conservative: approximately `4900` times the realized
correction in the reviewed fixture, while still only approximately `7e-8` of the relevant `A_p`
response scale. It is a conservative enclosure, not a precision estimate.

### Adapter-owned EOS and tail completion

Tail/EOS completion is adapter-owned. The generic Phase-5B API derives no universal tail model.
The validated Track-R free-gas fixture has explicit validated tail authority. A future realistic
EOS or a different table/domain requires its own validated tail/EOS adapter before the generic
result's error bounds may be applied. This is a scope boundary, not an implementation defect.

### PB6 classification

PB6 is partition/radial-resolution/threshold-refinement sensitivity validation. It detects lost
sharp continuous composition, node-placement sensitivity, and partition pathology. It is not an
independent physics oracle; independent physics support comes especially from PB4, PB5, PB7,
PB10, and PB11.

## 8. PB6 shared-endpoint repair ratification

The owner explicitly ratifies the PB6 representation repair. The original refusal was
`STOP unrepresented gap/jump in number measure` for neutron, segment 2517, canonical node 2516.
The radius differed by 0 ULP and the reconstructed density by 5 ULP, approximately `7.5e-16`
relative. This was floating-point endpoint reconstruction, not a physical jump
(`docs/validation/PHASE5B_INV09_GLOBAL_RESPONSE_IMPLEMENTATION.md:17`).

The accepted repair constructs each canonical endpoint once and shares it between adjacent
segments. Exact continuity, declared-jump, continuous-onset, terminal-atom, and domain-boundary
semantics are unchanged. No epsilon continuity tolerance was introduced.

## 9. Measure-complete `dn_i` representation

The owner explicitly ratifies `- integral w xihat dn_i` as the production representation of the
displaced-composition contribution. Here `n_i=Y_i n_B`; no lapse appears; continuous onsets have
no atom; true declared jumps use exact signed atoms; and the terminal contribution appears
exactly once. The independent review found this compatible with FR2005 Eq. (24) and Appendix-B
distribution/jump treatment. This is an accepted rigorous numerical representation of
source-backed physics, not a publication-novelty claim. The governing derivation is
`docs/validation/PHASE5B0_INV09_GLOBAL_RESPONSE_PREFLIGHT.md:223`.

## 10. Ratified Track-R primary-fixture values

Species order is `[n,p,e,mu]`.

```text
N [count] =
[7.5683115419394785e56,
 4.8594976779411839e54,
 4.8119204549793175e54,
 4.7577222961866916e52]

A [count km^2] =
[3.9647718619373394e59,
 1.8404826326625264e57,
 1.8292534582866583e57,
 1.1229174375868009e55]

B [count km^2] =
[1.6255270466118125e59,
 5.7613767070193853e57,
 5.1928548580616960e57,
 5.6852184895566917e56]

A_B = 3.9831766882639647e59 count km^2
B_B = 1.6831408136820063e59
(d epsilon_c/dq)_(N_B) = -2.3665142309456835

K [count km^2] =
[+1.1793897334337809e58,
 -1.1793897334337830e58,
 -1.0459711462551772e58,
 -1.3341858717812756e57]

conditional whole-star I_phys [count s^2] =
[+1.3122480530141585e47,
 -1.3122480530141607e47,
 -1.1637998545112904e47,
 -1.4844819850233820e46]
```

These are ratified implementation results for the declared fixture and domain, not published
FR2005 numerical targets. The exact machine-readable values and units remain in
`docs/validation/phase5b_structural_response_candidate.json:1`.

## 11. Source, domain, and Structure-1 qualifications

The whole-star mapping `I_geom,i^whole=K_i` is ratified only under the documented `P=0`
boundary assumptions. The explicit fixed-isobar core/shell mapping is ratified with its
boundary-flux terms. No source-qualified free-gas core `I_Omega` benchmark, source-authenticated
free-gas core boundary, or reproduced FR2005 Figure-1 core coefficient is claimed. No core
cutoff is invented. The governing mapping is
`docs/validation/PHASE5B0_INV09_GLOBAL_RESPONSE_PREFLIGHT.md:392`.

Phase-5B preserves the qualified Structure-1 common-state claim at
`rho_c=1.10e15 g/cm^3`, with the printed source-compatible result
`0.62 / 12.77 / 13.80`. Source-qualified `M_max` remains **UNRESOLVED**; the candidate uses no
maximum-mass tuning (`docs/validation/TRACKR_FREEGAS_WHOLESTAR_STRUCTURE1_RATIFICATION.md:1`).

## 12. Candidate artifact, suites, and historical immutability

The owner ratifies
`docs/validation/phase5b_structural_response_candidate.json` as the authorized content for a
later first governed Phase-5B structural regression artifact. Its SHA-256 is
`7588f0e9cd62f5b6be48bb725e4d0ba6e47b64f50d843117e70c17a78beeb5fa`.
Independent review regenerated it twice byte-identically. This task does not install, move, copy,
or modify it; ratification authorizes only a later governed installation step.

The reviewed results are focused `8/8 PASS`, data-free `38/38 PASS`, and complete external-data
`61/61 PASS`. Every pre-existing governed baseline remains byte-identical, and the historical
PB6 stop record remains immutable (`docs/validation/phase5b_resume_evidence.json:1`).

## 13. Exact remaining boundary

After this decision:

- ADR-0011 Phase-5B structural implementation is **HUMAN-RATIFIED**.
- PN1-PN8 are **HUMAN-RATIFIED** for the reviewed scope.
- PB1-PB14 are **HUMAN-RATIFIED AS THE REVIEWED CANDIDATE VALIDATION PACKAGE**, with the
  evidentiary caveats above.
- The deterministic candidate artifact is **RATIFIED FOR GOVERNED INTEGRATION**.
- INV-09 remains **INTENDED BUT UNVERIFIED** until governed artifact installation,
  post-installation validation, and canonical integration complete.
- INV-11 remains **UNRESOLVED**.

This ratification does not authorize corrected R2006 `Btilde`, chemical `Z` or `W`, evolved
`eta`, Urca/weak-rate evolution, rotochemical-heating evolution, BNV sources, BNV energy
deposition, BNV-induced chemical disequilibrium, or endothermic/exothermic BNV accounting.
Those activities remain blocked until INV-09 is canonically closed and their own prerequisites
are satisfied.

## 14. Exactly one recommended next action

Perform a separate governed Phase-5B integration task from the ratified branch: install the
ratified deterministic structural-response artifact as the first governed Phase-5B regression
artifact through its canonical producer, rerun the complete data-free and external-data suites,
independently authenticate the artifact and all historical baselines, fast-forward canonical
master only if all post-installation validation passes, and only then mark INV-09
VERIFIED/RESOLVED.
