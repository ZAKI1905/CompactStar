# Immutable Phase-5C-2 production acceptance predeclaration

Starting production SHA: `d3b102d2225c44bea581832b27cdc33eb5874e83`.
Authority: ADR-0013, reviewed Phase-5C-1 plan, and Phase-5C-1R ratification.
Scope: Structure-1 midpoint rho_c=1.10e15 g/cm^3, whole-star generic/free-gas
candidate only. No A18, eta evolution, weak rates, heating/cooling, or BNV.

**NO PRODUCTION G/Z/W RESULT HAS BEEN COMPUTED YET.**

This document and its accompanying evidence are to be committed with message
`docs: predeclare chemical coefficient acceptance` before the first production
coefficient output. The commit then defines PHASE5C2_PREDECLARATION_SHA; its
hash will be recorded in the later implementation record. Goals are immutable
after that first output: no amendment, rewriting or force-push over the commit.
Failure of a goal is refusal, never permission to fit a replacement goal.

## Exact absolute goals

Entrywise inequalities apply; no universal relative tolerance. G uses axes
(Neutron, Electron, Muon), units count/MeV. Q uses (Electron, Muon), units
count/MeV. Z uses named rows (Npe, NpMu), lepton-number columns (Electron,
Muon), units MeV/count. W goals use (Npe, NpMu), units MeV s^2.

G numerical_error absolute goal:

```text
[2.0000000000000002e+50, 1e+40, 1e+40]
[1e+40, 4e+46, 6.9999999999999998e+44]
[1e+40, 6.9999999999999998e+44, 6.0000000000000002e+45]
```

Q numerical_error absolute goal:

```text
[5.0000000000000002e+46, 1.9999999999999999e+45]
[1.9999999999999999e+45, 7.9999999999999994e+45]
```

Z numerical_error absolute goal:

```text
[1.9999999999999999e-60, 1.9999999999999999e-60]
[1.9999999999999999e-60, 9.0000000000000002e-59]
```

G's source-qualified structural zeros are n-e/e-n/n-mu/mu-n. Each has
absolute numerical goal 1e40 count/MeV; zero is never evaluated with a relative
goal. Exact lower-dimensional absent support remains zero by embedding, not
by padding or inversion of an inactive Hessian. Synthetic fixtures use their
own independently analytic absolute roundoff/reference budgets, not these
large whole-star dimensional scales.

E_I_numerical (e, mu), count s^2: `[2.082871228648889e+40, 9.660006381851622e+39]`.

V_I_validation (e, mu), count s^2: `[4.880818755395441e+42, 5.482892414134075e+41]`.

G_W_numerical = `[1e-12, 4e-12]` MeV s^2.

G_W_validation = `[1e-10, 3e-10]` MeV s^2.

Both independently required:

```text
E_W_numerical <= G_W_numerical
V_W_validation <= G_W_validation
E_W = |Z| E_I + E_Z |I| + E_Z E_I + E_W_arithmetic
V_W = |Z| V_I
```

Numerical and validation tracks must never be merged. Runtime results consume
the complete current structural result and owned dependencies, not baseline JSON.
Frozen structural validation inputs are explicitly versioned fixture authority.

## Pre-production basis for the goals

The accepted planning artifact has authenticated SHA-256
`722ac4162bad21799ed9fc813698665d3c9bd86cc26b75c4dc646a2c019b8e88`.
Its reported E_G is approximately [[1.02424e50,1.83494e37,1.12886e36],
[1.37970e37,1.78521e46,2.98931e44],[5.49369e35,2.98931e44,2.80209e45]].
E_Q is approximately [[2.38626e46,5.55704e44],[5.55704e44,2.81301e45]],
and E_Z approximately [[4.92117e-61,4.38990e-61],[4.38990e-61,2.96826e-59]].
The goals above are independently chosen rounded absolute ceilings above this
pre-production characterization, allowing explicit global-solve/Schur/Z
arithmetic terms and production endpoint/partition implementation differences.
They are not a promise that the candidate will pass.

The supporting plan includes GL8/16/32 and partition bisection characterization,
radial/table and independent-background comparisons, GC9 analytic accuracy,
provider characterization, refusal windows, center and tail estimates, and
conditioning/rho methodology. Production must remeasure its own achieved
errors and the full convergence ladder. Historical precursor errors are not
silently assigned as achieved production uncertainties.

The known old/corrected named Z separations are approximately
(6.98847e-56,2.34281e-55,3.71025e-54), far larger than these numerical goals.
The independent old-route uncertainty must also enter the candidate separation
gate; its tolerance is not fitted to historical percentages.

Accepted W diagnostics were approximately (-5.40618e-7,-1.58457e-6) MeV s^2.
Combining accepted Z/E_Z with freshly measured E_I gives pre-production
E_W estimates (1.64167e-13,1.49446e-12), before the separately required
arithmetic term. Combining accepted Z with V_I gives empirical envelopes
(2.26344e-11,5.88269e-11). The rounded W ceilings preserve substantial
separation from the percent-scale known free-gas correction while allowing
conservative structural validation; they were selected without production W.
The evidence JSON retains exact inputs and these diagnostic estimates.

## N1–N9 binding requirements

| Requirement | Candidate criterion |
|---|---|
| N1 | Actual production integration kernel accepts injected metric/background and local C evaluators; GC9 runs that kernel with an independent analytic reference |
| N2 | A positive-width segment continuous-onset at both endpoints splits into two halves, each mapped from its own onset; measure/support/convergence tested; no threshold H; first-order discontinuities refuse |
| N3 | Ratified separate numerical_error and validation_envelope; every required owned error term present; dual W gates |
| N4 | Explicit E_global_solve, E_Schur_arithmetic, E_Z_arithmetic, E_W_arithmetic and cross terms retained |
| N5 | Conservative tail radius with consistent mass-upper enclosure where required, analytic inequality, direct same-cut and independent high-precision checks |
| N6 | Every refusal interval has authenticated containing cells, full containment, branches, availability authority and whole-interval geometry extrema; negative tests |
| N7 | Outward-rounded/ULP-aware neutron source enclosure covers provider boundary with source-derived representational margin; otherwise refuse |
| N8 | Resolvable highly conditioned SPD may pass absolute goal; uncertain smallest mode, indefinite, nonsymmetric beyond budget and unresolved global Q refuse |
| N9 | Both finite-refusal edges are explicit partition boundaries; excluded interval contributes no quadrature values, only bounded additive error; provenance retained |

## GC1–GC14 binding candidate criteria

Independent expected mathematics in the preflight and plan is preserved; actual
candidate operations are redirected through production. Algebraic identities
retain CONTRACT status and never substitute for an independent physics oracle.

| Gate | Required candidate evidence |
|---|---|
| GC1 | Unit, index and signed eta/number conventions through typed production access |
| GC2 | Independent neutral species reconstruction; wrong map separated |
| GC3 | Coupled analytic energy Hessian, finite-perturbation ladder and stability |
| GC4 | Production solve action versus independent reference; N8 positive and refusal fixtures |
| GC5 | Independent x/y energy-response equivalence |
| GC6 | Independent full-intrinsic corrected oracle versus production neutral path |
| GC7 | Full singular inverse and pseudoinverse forbidden; second projection impossible/refused |
| GC8 | npemu 3D, npe 2D, pe 1D; vacuum and both thresholds value-only; no padded/threshold H; onset/refusal controls |
| GC9 | Correct injected production curved-star kernel passes independent reference; omit lapse, extra inverse lapse, exp(+nu), nu/2Phi, omitted and inverted curvature mutants separated by analytic accuracy goals |
| GC9b | GL8/16/32, sealed partition refinement, table/profile refinement, independent background, center, refusal, tail and roundoff; no single-grid/order pass |
| GC10 | Exact two-zone integrate-before-Schur oracle; local-Schur mutant fails; Q derived only, never independent authority |
| GC11 | Source traceability of named R2006 Z from independently supplied synthetic G/M; semantic channels explicit; not independent physics validation |
| GC12 | Both W gates pass; signed Omega_dot and missing/double c^-2, omit e/mu and swap-channel controls; W=ZI alone remains CONTRACT |
| GC13 | Realistic A18 SOURCE-LIMITED / BLOCKED; no arbitrary APR or figure substitution |
| GC14 | All material dependency mutations refuse before values: profile identity/version, provider revision/bytes and same-label changed bytes, domain, partition, onsets/refusal/tail, accuracy goal, G/Z, structural/central/sequence sources, EOS bytes, lifetime owner/token; alive/mutated/released/foreign cases safe |

The separate free-gas correction gate requires each named Z_npe, Z_np,
Z_npmu and matrix separation to exceed combined candidate/old-route numerical
uncertainties, with reported margin. Free gas never substitutes for A18.

## Refusal and release gates

Missing/stale/foreign source, unavailable lifetime ownership, invalid basis or
metric, illegal active set, genuine discontinuity without interface authority,
unbounded refusal/tail/center term, indefinite or unresolved supported mode,
uncertainty-inclusive denominator failure, nonsymmetry beyond budget,
unresolved quadrature or representation evidence, missing numerical terms,
AccuracyGoalUnmet, and StructuralValidationEnvelopeUnmet all refuse.
No universal kappa cutoff, epsilon denominator, or fabricated rank is permitted.

Candidate release additionally requires all focused, complete data-free and
complete authenticated external-data suites with raw rc=0; no raw failure is
absorbable. Reproduce the candidate JSON independently twice, byte-identical,
with all required false scope flags and no timestamp/absolute path. Preserve
all nine historical baselines byte-for-byte; no candidate baseline installation.
No ratification, integration, merge or downstream evolution is authorized.
Independent review is required after candidate implementation.
