# ADR-0013 — corrected rotochemical chemical coefficients

**Status:** ACCEPTED
**Decision:** HUMAN-RATIFIED WITH EXPLICIT INDEPENDENT-REVIEW REVISIONS
**Date:** 2026-09-06
**Starting canonical SHA:** `49ab2b8c2881b6ef7b9309307d18cea51d557f72`
**Preflight SHA:** `54ec7abac38fa0a32c5fb3a82e424b496361966a`
**Change class:** scientific-semantic and structural/architecture decision; documentation-only
ratification.
**Evidence companions:**
`docs/validation/PHASE5C0_CORRECTED_CHEMICAL_COEFFICIENT_PREFLIGHT.md` and
`docs/validation/PHASE5C0_CORRECTED_CHEMICAL_COEFFICIENT_RATIFICATION.md`.
**Implementation state:** no production `G_y`, `Z`, `W`, chemical evolution, weak rate,
heating/cooling, superfluid, BNV, EOS, test, baseline, or source-data implementation is
authorized or added by this decision.

> **PHASE-5C CORRECTED CHEMICAL-COEFFICIENT CONTRACT — HUMAN-RATIFIED;
> PRODUCTION IMPLEMENTATION NOT YET AUTHORIZED.**

## 1. Context, authority, and review disposition

ADR-0010 is ACCEPTED and governs the cold local charge-neutral Hessian, active-species domains,
and the corrected R2006 interpretation. ADR-0011 is ACCEPTED and its ordinary-`NStar`
structural particle-number response is canonically integrated; INV-09 is VERIFIED / RESOLVED
within that structural scope. INV-11 remains UNRESOLVED. This ADR governs only the chemical
coefficient layer between those capabilities and a future secular-evolution layer
(`docs/adr/ADR-0010-rotochemical-off-equilibrium-thermodynamic-contract.md:5`,
`docs/adr/ADR-0010-rotochemical-off-equilibrium-thermodynamic-contract.md:212`,
`docs/validation/PHASE5B_INV09_GLOBAL_RESPONSE_INTEGRATION.md:127`,
`docs/SCIENTIFIC_INVARIANTS.md:958`).

The governing corrected primary source is
`/Users/keeper/Documents/CompactStar/literature/rotochemical/2006-Reisenegger-Rotochemical-Heating-of-Neutron-Stars-Rigorous-Formalism-with-Electrostatic-Potential-Perturbations.pdf`,
SHA-256 `a286f15e083e52becd95b3000cbb5ec3ed97148681cf10a43f1a1cc5c4d23ae8`,
journal pp. 569–571 / PDF pp. 2–4, especially eqs. (10)–(19). F2005 supplies the retained
non-superfluid framework and closing benchmark where R2006 does not supersede it. Exact source
hashes, pages, supersession boundaries, and scratch evidence are recorded in the preflight
sections 2–15.

The independent Opus review examined the preflight commit in the read-only detached worktree
`/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase5c-preflight-review`; no review
commit exists. Its disposition was:

> **PHASE-5C PREFLIGHT INDEPENDENT REVIEW PASS WITH NONBLOCKING FINDINGS —
> ADR-0013 READY FOR HUMAN RATIFICATION WITH EXPLICIT CAVEATS.**

The review reported no blocking finding, four material contract/documentation findings, and
nine nonblocking clarifications or plan-strengthening findings. The owner ratifies this ADR only
with the revisions recorded below. Proposal history is preserved in the preflight and in the
revision ledger; acceptance is not implementation validation.

## 2. Ratified scientific scope

The accepted contract covers:

- cold bulk npe-mu chemical-coefficient mathematics;
- the corrected R2006 electrostatic response;
- local neutral Hessian to susceptibility;
- GR global integration and global baryon reduction;
- chemical `Z` and spin-drive `W` semantics;
- active-species branch embedding;
- source, domain, redshift, provenance, conditioning, and refusal semantics; and
- the future validation architecture.

It does not cover production implementation, realistic A18 closure, reaction rates,
eta-state storage/evolution, thermal or neutrino accounting, changing coefficients/backgrounds,
superfluidity, or BNV. No production `Btilde/G/Z/W`, evolution equation, rate, baseline, EOS/data,
or literature source is created or modified by this decision.

## 3. Ratified mathematical contract

### 3.1 Local corrected response

Use canonical local coordinates and conjugates

```text
x = (n_B,n_e,n_mu)^T,
g_x = (mu_n,-eta_npe,-eta_npmu)^T,
H_x = partial g_x / partial x.
```

Use the source basis

```text
y = (n_n,n_e,n_mu)^T,
T = [[1,-1,-1],
     [0, 1, 0],
     [0, 0, 1]],
y = T x,
H_x = T^T H_y T,
C_x = H_x^-1,
C_y = T C_x T^T = H_y^-1.
```

The inverses above are qualified solves on the declared active physical branch, not a padded
full-species inverse. The accepted neutral-Hessian route already incorporates electrostatic
elimination. No second projection is allowed.

Only an independent validation fixture may construct a full intrinsic susceptibility `chi` and
compare against

```text
q = (0,1,-1,-1)^T,
C_projected = chi - chi q (q^T chi q)^-1 q^T chi.
```

The full corrected four-species response is singular in the charge/gauge direction. A full
corrected paper-`B` inverse, a silent pseudoinverse, and a public charged-direction response are
forbidden.

### 3.2 Global number response

For a declared connected diffusive chemical domain `D`, with `nu=Phi`, `r,m` in km and
`C_y` in `fm^-3 MeV^-1`, the canonical global object is

```text
G_y = 10^54 integral_D
      4 pi r^2 e^-nu C_y(r) / sqrt(1-2m/r) dr.
```

`G_y` is represented by `GlobalChemicalNumberResponse` in the reduced source basis
`y=(N_n,N_e,N_mu)`, with named row and column axes and units `count / MeV`. It is the canonical
unreduced authority. An x-basis transform may exist only as a derived diagnostic or accessor.
A bare class named `B` is forbidden.

Exactly one inverse lapse occurs in `G_y`; the proper-volume factor is
`1/sqrt(1-2m/r)`. The positive lapse used for reaction-rate conversion is a different future
operation and does not alter the coefficient integral.

### 3.3 Global baryon reduction and chemical response

Local charge neutrality is imposed locally. Global baryon conservation is imposed only after
stellar integration:

```text
local charge-neutral susceptibility
  -> integrate G_y
  -> global baryon reduction
  -> chemical imbalance response.
```

Define

```text
b = (1,1,1)^T,
L = [[-1,-1],
     [ 1, 0],
     [ 0, 1]],
delta N_y = L delta N_l,
eta^infinity = -L^T delta g_y^infinity,
Z = L^T G_y^-1 L,
eta^infinity = -Z delta N_l.
```

Equivalently, after integration transform to
`G_x=[[a,h^T],[h,D]]`, derive the conditioning/reduction object

```text
Q = D - h h^T/a,
Z = Q^-1.
```

`Q` is derived, not a second independently stored authority. No pointwise baryon constraint or
local Schur reduction is allowed. A retained beta channel must be supported on positive measure,
and its smallest eigenvalue must exceed its uncertainty enclosure before inversion. Otherwise a
separately named lower-rank channel result is returned or the operation refuses.

`ChemicalImbalanceResponse` owns one immutable symmetric 2x2 matrix in named order
`BetaChannel::Npe`, `BetaChannel::NpMu`:

```text
Z = [[Z_npe, Z_np],
     [Z_np,  Z_npmu]]       [MeV / count].
```

The paper scalars are read-only accessors onto that matrix and are never independently stored.
Typed row/output-channel and column/input-lepton orientation remains required even though the
matrix is symmetric.

### 3.4 Structural drive and coefficient redshift semantics

`RotochemicalSpinDrive` is a separate immutable result that depends on the complete governed
`FixedBaryonNumberResponse`, not on a copied `I` vector:

```text
W = Z (I_phys,e,I_phys,mu)^T       [MeV s^2].
```

The initial whole-star implementation will consume `WholeStarIPhysical()` only after currency
validation. Runtime coupling to a baseline JSON is forbidden. A core or other reservoir requires
its own governed structural mapping and boundary/source semantics.

At coefficient level,

```text
eta^infinity = e^nu eta_local,
eta_npe = mu_n - mu_p - mu_e,
eta_npmu = mu_n - mu_p - mu_mu,
dot eta^infinity = -Z R + 2 W Omega Omega_dot
```

for frozen coefficients and the declared source sign. The last line fixes source sign semantics
only; it is not a secular-evolution implementation contract. `W` is a spin-drive adapter, not a
universal chemical source. Future externally supplied particle-number sources, including any
baryon-changing source, require their own conservation/source mapping and must not be forced into
the fixed-baryon lift `L`.

## 4. Owner-ratified decisions Q1–Q8

| ID | Ratified answer |
|---|---|
| **Q1 — canonical global chemical number response** | `GlobalChemicalNumberResponse` owns named-axis `G_y` in the reduced source basis `(N_n,N_e,N_mu)`, units `count / MeV`. `G_y` is retained as the source-aligned unreduced physical authority for later non-fixed-baryon sources. An x-basis view is derived only. No full corrected 4x4 inverse, pseudoinverse, or bare `B` class. |
| **Q2 — global baryon-reduction ownership and supported-mode policy** | Integration produces unreduced `G_y` first; global baryon reduction occurs afterward. No pointwise baryon constraint. Retain `G_y`; derive `Q` rather than store a second authority. Supported channel rank is explicit; positive-measure support and an uncertainty-resolved smallest eigenvalue are required before inversion, otherwise return a named lower-rank result or refuse. This ordering was already fixed by ADR-0010 and R2006 footnote 4 and is not a free scientific choice. |
| **Q3 — canonical corrected chemical Z** | `ChemicalImbalanceResponse` owns one immutable symmetric 2x2 named-channel matrix with canonical relation `eta^infinity=-Z delta N_l`, units `MeV / count`. `Z_npe`, `Z_np`, and `Z_npmu` are views only. Typed output-channel/input-lepton orientation is required. |
| **Q4 — W ownership and structural coupling** | `RotochemicalSpinDrive` is separate and immutable, depends on `ChemicalImbalanceResponse` and the complete governed semantic `FixedBaryonNumberResponse`, and uses `WholeStarIPhysical()` for the initial whole-star path after currency validation. No copied `I` vector or runtime baseline JSON is sufficient provenance. Spin down is one source adapter, not the universal chemical source. |
| **Q5 — domains, branches, thresholds, and interfaces** | Use explicit 3D npemu, 2D npe, 1D pe, and 0D/value-only vacuum branches, plus value-only `MuonThresholdEvaluation` and `NeutronThresholdEvaluation`. Form `C_y=E H_active^-1 E^T` using authenticated branch conjugates. No padded H, absent-species response, density floor, hidden extrapolation, or threshold inversion. Continuous onsets have no density-jump atom. A finite refusal window is bounded or handled by a validated limit adapter, otherwise refuse. A genuine first-order transition needs authenticated interface/phase metadata and a chemical interface-motion law; the Phase-5B structural jump formula is not automatically reusable. |
| **Q6 — correction-sensitive Track-R benchmark instrument** | A18 + delta-v + UIX* remains the mandatory realistic Track-R closing model already fixed by ADR-0010 Q1; CMF, BPAL, and free gas do not substitute. The preferred coefficient benchmark is R2006 Figure 1 with authenticated matching A18 authority; authenticated author arrays are preferred, with governed figure extraction permitted if unavailable. A later transient may supplement, not replace, the coefficient benchmark; quasi-steady temperature alone is insufficient. Realistic closure remains blocked on missing authenticated A18 authority. |
| **Q7 — coefficient redshift and source semantics** | Fix `eta^infinity=e^nu eta_local`, named beta signs, exactly one `e^-nu` in `G_y`, Z acting on redshifted imbalance, and W units/sign/source semantics now. This partially resolves INV-11(a) only for coefficient-object semantics. It does not choose an evolved layout, storage, rates, changing-coefficient treatment, thermal/neutrino accounting, or solver coupling. INV-11 remains UNRESOLVED. |
| **Q8 — provenance, lifetime, and refusal** | Future results are immutable, dependency-complete, fail-closed values. Every scientific access validates currency. Raw pointers without guaranteed lifetime or a validated lifetime token are insufficient. No universal condition-number cutoff is accepted; stability/refusal is uncertainty-aware, and the smallest supported eigenvalue must exceed its error enclosure before inversion. |

## 5. Active-branch and interface contract

For each active chart `z`, form `C_y=E H_active^-1 E^T`:

| Branch/result | Active dimension | Embedding into `(n_n,n_e,n_mu)` | Boundary rule |
|---|---:|---|---|
| npemu | 3D | full `T` | qualified smooth branch |
| npe | 2D | `[[1,-1],[0,1],[0,0]]` | no fabricated muon response |
| pe | 1D | `[[0],[1],[0]]` | active conjugate is `mu_p+mu_e` |
| vacuum | 0D/value-only | none | no Hessian inversion |
| `MuonThresholdEvaluation` | value-only | none | no Hessian inversion exactly at threshold |
| `NeutronThresholdEvaluation` | value-only | none | no Hessian inversion exactly at threshold |

Continuous onset uses the one-sided susceptibility limit and has no density-jump atom. Every
positive-width response-refusal window crossing physical support requires an explicit
response-measure bound or validated limit adapter; it may not be silently skipped. A genuine
first-order phase boundary requires authenticated phase/interface metadata and a chemical
interface-motion law. Structural interface terms from Phase 5B do not automatically define the
chemical interface response.

## 6. Provenance and lifetime contract

Every future global chemical-response result retains:

- `StarProfile` identity and version;
- metric/redshift identity and normalization;
- provider identity and revision;
- exact provider/data bytes where applicable and component constants;
- equilibrium anchor and active branch map;
- source/result basis and named row/column orientation;
- chemical domain/reservoir and interface/tail policy;
- realized radial partition, quadrature-rule identity, onset split locations, and achieved
  node count;
- numerical method and version;
- a structural-zero entry register;
- matrix eigenspectrum, conditioning, and error budgets; and
- the R2006 source convention.

`RotochemicalSpinDrive` additionally retains all governed Phase-5B structural provenance. A
changed dependency refuses before scientific access. No lazy stale science or reconstruction
from labels is allowed. Lifetime safety is normative: a provenance object containing raw
pointers without a guaranteed lifetime or validated lifetime token is insufficient. This ADR
does not redesign the existing Phase-5B provenance implementation.

## 7. Accepted validation architecture

### 7.1 Validation ladder

The future coefficient validation order is:

1. exact analytic/toy corrected projection;
2. Track-R free-gas whole-star coefficient mechanics;
3. a free-gas old-F2005-versus-corrected-R2006 separation gate; and
4. authenticated A18 + delta-v + UIX* R2006 Figure-1 or author-array coefficient comparison.

Free gas can establish that the correction machinery is active; it cannot validate realistic
nuclear interactions or replace the A18 benchmark. A later transient comparison belongs
primarily to the evolution ladder and may supplement coefficient validation. Quasi-steady
temperature alone does not validate the correction.

The preflight's `GC1`–`GC14` ladder is accepted with these binding revisions:

- `GC9a` checks the analytic/manufactured and independent curved-GR integral, including lapse,
  proper volume, count conversion, center, and tail terms.
- `GC9b` is a separate required convergence/numerical-budget subgate: radial refinement,
  quadrature-rule comparison, onset-aware partition refinement, table/provider resolution where
  applicable, and surface/tail remainder. One correct manufactured integral cannot by itself
  make `GC9` pass.
- `GC9` negative controls include `M19` and the companion proper-volume inversion control below.
- `GC13` contains both the intermediate free-gas correction-separation requirement and the
  separately mandatory realistic A18 source benchmark. Passing the former never discharges the
  latter.
- `GC14` includes lifetime-token/dangling-dependency refusal and the full numerical provenance
  set in section 6.

### 7.2 ADR-0010 V1–V12 to Phase-5C GC crosswalk

`GC1`–`GC14` **extend and operationalize** the accepted ADR-0010 ladder; they do not supersede it.

| ADR-0010 gate | Phase-5C discharge/carry-forward |
|---|---|
| **V1** units/rest-mass/index/sign | `GC1` |
| **V2** exact neutral reconstruction | `GC2` |
| **V3** beta equilibrium and threshold conditions | Existing validated Phase-5A provider gates, plus `GC8` when global branch support is consumed |
| **V4** analytic lepton checks | Existing Phase-5A Track-R validation; inherited dependency, not re-credited as new Phase-5C validation |
| **V5** analytic toy reduced Hessian/susceptibility | `GC3`–`GC6`, as applicable |
| **V6** Hessian symmetry/integrability | `GC3` |
| **V7** finite perturbation versus linear response | `GC3` and inherited local-provider validation |
| **V8** x/y response equivalence | `GC5` |
| **V9** full intrinsic electrostatic projection versus neutral route; charge null/rank/proton identity | `GC6` + `GC7` |
| **V10** rank/support/stability/active-species handling | `GC8` plus conditioning portions of `GC10` |
| **V11** corrected global Z/source response | `GC9` + `GC10` + `GC11` + `GC13` |
| **V12** end-to-end published non-superfluid thermal benchmark | **Not discharged by Phase-5C coefficient implementation; remains a later evolution-layer gate** |

### 7.3 Required negative controls added by review

**M19 — sign-flipped coefficient lapse.** The wrong route is

```text
G_wrong = integral e^(+nu) C dV
```

instead of `G=integral e^(-nu) C dV`. The detector is an independent curved-GR exact-star
oracle. The independent review found approximately `-35%` to `-38%` relative error for this
wrong route. That separation is diagnostic evidence, not a production tolerance.

**M20 — inverted proper-volume factor.** `GC9` must also detect use of
`sqrt(1-2m/r)` instead of `1/sqrt(1-2m/r)`. This companion control may remain a `GC9`
subcase rather than a separately scored top-level mutation, but it must be recorded and run.

## 8. Review findings and accepted dispositions

### 8.1 Material findings M1–M4

| Finding | Accepted disposition |
|---|---|
| **M1 — V/GC crosswalk** | Section 7.2 explicitly preserves and maps ADR-0010 V1–V12. The GC ladder extends/operationalizes rather than replaces it; V12 remains future evolution validation. |
| **M2 — Q2/Q6 not free scientific choices** | Q2 is reframed around ownership/support policy while retaining the already-fixed integrate-first/global-reduction order. Q6 ratifies the concrete correction-sensitive benchmark instrument while retaining the already-fixed A18 + delta-v + UIX* closing identity. |
| **M3 — sign-flipped lapse route** | `M19` and the companion inverted proper-volume control are binding `GC9` negative controls, with the independent curved-GR exact-star oracle as detector. |
| **M4 — free-gas correction sensitivity** | The explicit intermediate free-gas old-vs-corrected separation gate is required, while `GC13` still requires authenticated A18 realistic-source comparison. Free gas does not replace A18. |

The review's free-gas diagnostics found approximate old-versus-corrected changes of `+1.55%`
for `Z_npe`, `+82.8%` for `Z_np`, `+3.75%` for `Z_npmu`, and `3.6%` in matrix Frobenius norm.
These values show detector separation only; they are not literature targets, golden data, or
production tolerances.

### 8.2 Nonblocking findings N1–N9

| Finding | Accepted disposition |
|---|---|
| **N1 — quadrature diagnosis** | Independent review localized the approximately `1.9e-6` scratch difference mainly to trapezoidal quadrature on the realized profile partition. Two independent non-trapezoid routes agreed at approximately `1.8e-8` while each differed from the trapezoid diagnostic by approximately `1.9e-6`. This is diagnostic evidence only. Quadrature is the dominant characterized scratch error; production must choose and validate an onset-aware policy before acceptance budgets are set. |
| **N2 — tolerance methodology** | Quadrature truncation is a first-class uncertainty. Record realized radial partition, onset splits, quadrature rule, achieved nodes, and structural-zero register. Use absolute budgets for analytically exact zeros such as free-gas `G_ne` and `G_nmu`; never relative error against zero. |
| **N3 — provider refusal window** | The reviewed fixture gave `Delta G_nn/G_nn <=` approximately `1.7e-14` and induced relative Z effect `<=` approximately `9.8e-18`. These are not universal tolerances. The ratified rule is bound every finite response-refusal window crossing support or use a validated limit adapter; otherwise refuse. |
| **N4 — chemical tail** | The reviewed Track-R fixture gave actual `Delta G_ee=2.685e45 count/MeV` within the proposed `3.682e45 count/MeV` bound, and `R_upper` agreed with the review's exact `P=0` radius to approximately `2e-10` relative. The shell-mass correction lies below independently demonstrated background reproducibility and is not dominant. The tail bound is not the full coefficient error budget. |
| **N5 — Figure-1 mass selection** | Future governed extraction should predeclare masses approximately in the `1.0`–`1.2 M_sun` region, subject to authenticated evidence, because old/corrected identity and separation are clearer there; the review found overlaps/crossings around `1.4`–`1.6 M_sun` in a figure spanning roughly `1.0`–`2.0 M_sun`. No numerical curve value is frozen and visual estimates are not source data. |
| **N6 — extraction protocol** | Require axis-calibration closure with residuals in uncertainty; a curve-identity rule and predeclared ambiguity-zone rejection near crossings; rasterization independence using separate resolutions and preferably renderers; and a source-text round-trip trend check before promotion to benchmark data. Independent extractors must not merely inspect the same raster. No digitization occurs here. |
| **N7 — INV-11 boundary** | ADR-0013 partially resolves INV-11(a) only for coefficient-object semantics: redshifted eta, named beta channels, one inverse lapse in G, Z action, and W units/sign/source. Evolved ordering, representation, units/conversion boundary, stoichiometry, rate-sign ownership, changing coefficients/background, thermal/neutrino partition, and solver coupling remain unresolved. INV-11 is not marked resolved. |
| **N8 — thresholds and lifetime** | `MuonThresholdEvaluation` and `NeutronThresholdEvaluation` are explicitly value-only; no threshold-object Hessian inversion occurs. Lifetime safety is normative for all dependency-complete result objects. Existing Phase-5B provenance is not redesigned here. |
| **N9 — convergence/quadrature gate** | `GC9b` is an explicit required subgate. A manufactured-integral pass alone cannot establish production `GC9` success. |

## 9. Source notes and realistic benchmark blocker

### 9.1 R2006 equation (11)

The owner ratifies ADR-0010's existing classification: printed R2006 eq. (11) appears to omit
`e^(-Phi)` required by substitution into eq. (10) and consistent with eq. (13). Classification:
**INFERRED PRINTED/SOURCE OMISSION**, not a published erratum. No source quotation is rewritten,
and no Fable adjudication is required.

### 9.2 R2006 Figure-2 mass label

R2006 Figure 2's right panel says `2.14 M_sun`, while its caption says `2.13 M_sun`.
Independent review verified both directly. F2005 gives `2.14 M_sun` specific physical
significance as the highest causal model described there, which favors the panel label, but that
is inference rather than an authenticated correction. This remains an unresolved **SOURCE NOTE**,
not an ADR blocker. Neither value may be selected for quantitative Figure-2 benchmarking without
later adjudication/source authority. Figure 1 can carry the coefficient-layer benchmark.

### 9.3 Exact A18 blocker

Track-R realistic closure is blocked until all matching A18 + delta-v + UIX* authority exists:

- authenticated equilibrium EOS/composition and exact model/fit lineage;
- arbitrary-composition nuclear energy/response with consistent lepton conventions;
- authenticated crust, joins, core boundary, and phase/interface treatment;
- matching mass/configuration and corrected Figure-1 author arrays or governed extraction; and
- source-compatible uncertainty sufficient for the declared comparison.

CMF, BPAL, a generic APR-labelled barotrope, and free gas do not substitute for that authority.

## 10. INV-11 partial-resolution boundary

Coefficient-object semantics now governed are:

- `eta^infinity=e^nu eta_local`;
- named `npe` and `np-mu` beta channels and signs;
- exactly one `e^-nu` factor in `G_y`;
- `Z` acting on redshifted imbalance;
- `W` units `MeV s^2`; and
- the frozen-coefficient source sign.

Still unresolved for secular evolution are evolved-state ordering, representation, storage
units and conversion boundary, reaction stoichiometry, net-rate sign ownership, changing
coefficients/background, thermal/neutrino partition, and solver coupling. Therefore INV-11
remains **UNRESOLVED** and continues to block chemical-state/evolution implementation.

## 11. Decision and revision ledger

**ACCEPTED — 2026-09-06.** The project owner ratifies Q1–Q8 and the mathematical, domain,
provenance, validation, and source-boundary contract above. Acceptance governs semantics and
ownership only. It does not certify scratch numbers, claim any `GC` pass, authorize production
implementation, close realistic Track R, or resolve secular evolution.

| Stage | Status and effect |
|---|---|
| Phase-5C-0 proposal, commit `54ec7abac38fa0a32c5fb3a82e424b496361966a` | ADR-0013 proposed; no owner answer, review disposition, or implementation authorization. |
| Independent Opus review | PASS WITH NONBLOCKING FINDINGS; no blockers, M1–M4 and N1–N9 recorded; detached read-only review, no review commit. |
| Phase-5C-0RAT owner revision and acceptance | Q1–Q8 ratified with Q2/Q6 reframed, V/GC crosswalk, `M19`/proper-volume control, free-gas separation rung, numerical/provenance strengthening, source notes, and INV-11 partial boundary. Documentation only. |

Downstream prohibitions remain explicit: no production `G_y/Z/W`; no eta evolution; no weak
rates; no heating/cooling; no superfluidity; no BNV; no tests, baselines, EOS/data, or literature
changes; and no canonical merge under this decision.

## 12. Owner-ratified two-track inherited structural uncertainty addendum — 2026-09-07

**ADR-0013 remains ACCEPTED. Q1-Q8 are accepted unchanged.** The mathematics and ownership of
`G_y`, `Q`, `Z`, and `W` are unchanged. In particular, `RotochemicalSpinDrive` still consumes
only the electron and muon whole-star structural inputs,

```text
I_phys,i = K_i/c^2,
W = Z (I_phys,e,I_phys,mu)^T,
```

through the existing governed unit owner and complete structural provenance
(`docs/adr/ADR-0013-corrected-rotochemical-chemical-coefficients.md:182`,
`docs/adr/ADR-0013-corrected-rotochemical-chemical-coefficients.md:215`). Neutron and proton
responses remain necessary to construct and validate the fixed-baryon response, but they are not
direct entries of the `W` vector. The large neutron cancellation amplification is therefore not
itself the uncertainty amplification of `W`; it does not remove any neutron/proton evidence
required to validate the construction.

### 12.1 Ratified terminology

The following terms are distinct and must not be interchanged:

- **`numerical_error`** is propagated uncertainty associated with the declared discrete
  representation/computation: arithmetic, quadrature on that representation, finite-difference
  or stencil estimates, local solve residual, roundoff, and a proven remainder explicitly
  included by the computation.
- **`certified_bound`** is a mathematically demonstrated enclosure under explicitly stated
  hypotheses. It is used only where such a proof exists. Neither `I_phys` nor `W` is assigned a
  certified bound by this addendum.
- **`validation_envelope`** is a conservative, predeclared envelope assembled from measured
  discrepancies against independent numerical/analytic routes and controlled refinement or
  representation variants. It is empirical validation evidence, not a probability distribution,
  confidence interval, formal truncation remainder, or mathematically certified continuum bound.

Phase-5B INV-09 closure remains **VERIFIED / RESOLVED** and retains all nine qualifications
(`docs/validation/PHASE5B_INV09_GLOBAL_RESPONSE_INTEGRATION.md:99`,
`docs/validation/PHASE5B_INV09_GLOBAL_RESPONSE_INTEGRATION.md:127`). Existing Phase-5B
`K_error` is `numerical_error`, not `certified_bound`; the baseline itself already separates
reported numerical errors from radial/EOS and independent-oracle evidence
(`tests/baselines/phase5b_structural_response.json:82`). This narrows no accepted structural
claim and changes no Phase-5B value, formula, or baseline.

### 12.2 Two inherited structural quantities

For each consumed species `i in {e,mu}`:

```text
E_I_numerical,i = K_error_i/c^2.
```

`E_I_numerical` is the Phase-5B `numerical_error` converted exactly once by the existing unit
owner. Separately, `V_I_validation` is the componentwise `validation_envelope` assembled from
predeclared end-to-end `K/I` evidence. It is publicly described as a validation envelope, never
as an error bar, confidence interval, certified bound, or formal remainder. A complete separate
componentwise `A/B` certified interval is not required because Phase-5C consumes `K/I`, although
`A/B` evidence remains part of validating that construction.

Before candidate acceptance, the `V_I_validation` evidence must include the existing stored K
numerical error as a separately classified input, PB11 direct `q -> 0` K discrepancy, PB11
fixed-baryon central-shift / `Delta N_B ~ q^2` evidence kept distinct from the Richardson K
comparison, PB7 independent-background B discrepancy transferred into the consumed K direction
with documented sensitivity/amplification, PB12 K-level EOS/table-resolution variation, a new
direct K-level radial-resolution ladder, PB6 partition/knot variation transferred to K with a
documented sensitivity factor where still A-level, raw per-species PB10 direct-K discrepancies,
and any separately reviewed independent K comparison retained without duplicate counting. The
Phase-5B records establish the distinct existing roles of PB6, PB7, PB10, PB11, and PB12
(`docs/validation/PHASE5B_INV09_GLOBAL_RESPONSE_RATIFICATION.md:99`,
`docs/validation/PHASE5B_INV09_GLOBAL_RESPONSE_RATIFICATION.md:225`).

No direct K-level radial-resolution ladder currently exists; it is required as an
implementation-time measurement before candidate acceptance, not before production code may be
written. PB10 must durably record raw per-species direct-K discrepancies, which may be test-side
and requires no Phase-5B production-source change. Neither measurement is performed by this
documentation addendum.

PB6 knot/profile shifts are raw representation sensitivity and
`VALIDATION_ENVELOPE_INGREDIENT` evidence, not violations of a certified bound. Dividing the
neutron/muon shifts by stored `A_error` does not create such a violation because `A_error` does
not claim that reconstruction class. For p/e, the existing tail contribution may dominate the
stored `A_error`, explaining why the same ratio may be below one. The variations must remain
visible in durable evidence.

### 12.3 Downstream propagation and dual acceptance

The Phase-5C numerical propagation remains

```text
E_W_numerical
  <= |Z| E_I_numerical
   + E_Z |I|
   + E_Z E_I_numerical
   + E_W_arithmetic.
```

`E_W_numerical` is `numerical_error`; no validation-envelope quantity is included in it. The
separate inherited structural stability envelope is

```text
V_W_validation <= |Z| V_I_validation.
```

`V_W_validation` is a `validation_envelope`. A future validation envelope for Z may contribute
only after separate governance.

For every required W component, candidate acceptance requires both

```text
E_W_numerical <= G_W_numerical
V_W_validation <= G_W_validation.
```

Both goals are componentwise, absolute, derived from pre-production authority, and predeclared
and immutable before the first production `G_y/Z/W` output. Failure of either causes refusal;
the structural-envelope failure is separately classified as
`StructuralValidationEnvelopeUnmet` (or an exactly equivalent repository-consistent name), not
as `AccuracyGoalUnmet`. This addendum selects neither numerical goal.

The clarification applies only to inherited Phase-5B structural representation stability. All
required Phase-5C-owned numerical errors in `G_y`, `Q`, and `Z` remain measured/bounded and
fail-closed under the accepted plan. An unmeasurable required owned numerical uncertainty causes
refusal; the validation-envelope distinction is not an escape from owned numerical accounting.

Ladder handling is also frozen before candidate results: use the finest measured discrepancy for
a contracting ladder with no demonstrated asymptotic order; record a nonmonotone/sign-alternating
numerical floor without inventing an order; refuse truly noncontracting or unbounded behavior;
retain every raw variant; and fit no after-the-result safety factor. Already included Phase-5B
`K_error`/tail terms are not added twice. Distinct PB11 K and fixed-baryon constraints may both be
retained with dependence stated; potentially correlated PB7 B and central-shift effects are
combined conservatively absent proved cancellation; PB6 radial/knot variants use a conservative
max/envelope within one representation class; and one numerical route is never relabelled and
counted twice.

### 12.4 Regression, realistic-source, and status boundary

A later governed Phase-5C regression artifact may be installed after acceptance under this
two-track contract. It certifies reproducible regeneration of accepted deterministic candidate
bytes, not a formally interval-certified continuum solution; this is the same reproducibility
distinction recorded for Phase-5B (`docs/validation/PHASE5B_INV09_GLOBAL_RESPONSE_INTEGRATION.md:44`).

The method separating `numerical_error` from `validation_envelope` transfers to realistic Track
R, but the numerical free-gas envelope does not transfer to A18. A18 requires its own adapter,
EOS/table and phase/interface evidence, radial evidence, correction-sensitive source benchmark,
and validation envelope. GC13 remains **SOURCE-LIMITED / BLOCKED** on authenticated A18 authority
(`docs/adr/ADR-0013-corrected-rotochemical-chemical-coefficients.md:379`). Publication claims may
report a clearly labelled numerical error ledger and validation envelope, but may not call the
latter a certified bound, confidence interval, or formal truncation error without separately
establishing those semantics.

At this ratification point **NO PRODUCTION `G_y`, `Z`, OR `W` RESULT EXISTS**. Production code,
tests, numerical-result generation, and baseline changes are outside this addendum. INV-11 remains
**UNRESOLVED**; eta evolution, weak rates, heating/cooling, realistic A18 closure, and BNV remain
outside the authorized scope.

## 13. Phase-5C-2RAT implementation-ratification addendum — 2026-09-07

**ADR-0013 remains ACCEPTED.** The production implementation candidate is commit
`4d78bf4000848ddecc2127daa2f2840872f266f5`. Independent Opus review returned disposition B:
**PASS WITH NONBLOCKING FINDINGS — candidate ready for human ratification with explicit
caveats**, with 0 blocking and 0 material findings. The human owner ratifies the implementation
with the explicit independent-review caveats recorded in
`docs/validation/PHASE5C2_CORRECTED_CHEMICAL_COEFFICIENT_RATIFICATION.md`.

For v1, the use of one `ImbalanceChannel` type for both `PaperZ` axes is accepted as a narrow
typing concession because the current symmetric two-channel matrix has a one-to-one
`Npe <-> Electron`, `NpMu <-> Muon` mapping and independent review found no numerical ambiguity.
This does not weaken the v1 numerical result and is not a general API precedent. Before output
and input spaces differ, the matrix can become nonsymmetric, species/channels are added, or the
ordering ceases to be one-to-one, the API must use distinct semantic `BetaChannel` output and
`LeptonInput` input types, or an equivalently strong typed representation. No API change is
authorized by this addendum.

`ChargeNeutralNumberSusceptibility::NumericalError()` is the local congruence,
factorization, and solve arithmetic uncertainty for a supplied local thermodynamic response; it
is not the total EOS/provider/background physical-model uncertainty. Provider and background
characterization remains separate at global level, especially in `E_background`. Future prose
must preserve that decomposition, and a realistic EOS provider must supply or separately govern
nonzero provider uncertainty where available.

The two-track UQ contract in sections 12.1-12.3 remains controlling: `numerical_error` and
`validation_envelope` are distinct, both immutable componentwise W gates passed, and
`V_W_validation` is not a certified bound, confidence interval, formal truncation error, error
bar, or achieved-accuracy estimate. The review caveats concerning PB11 dominance, the immaterial
PB7 transfer micro-difference, equivalent-but-not-literal M20 construction, current cellwise-linear
refusal extrema, and the thin positive old-M tail margin are retained in the ratification record.

GC13 remains **SOURCE-LIMITED / BLOCKED**; the generic/free-gas fixture is not realistic A18
closure. INV-11 remains **UNRESOLVED** for evolved chemical-state ownership, storage, coefficient
evolution, rates, and evolution coupling. This addendum authorizes no eta evolution, weak rates,
neutrino/heating evolution, realistic A18 closure, superfluidity, or BNV. The human-ratified
candidate still requires a separate governed canonical-integration task before Phase-5C closes.
