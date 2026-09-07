# Phase-5C-2I pre-production structural evidence

Status: M1/M2 measured successfully (raw rc=0); component envelopes fixed
before any production chemical computation. Candidate empirical evidence only.

Starting SHA: `d3b102d2225c44bea581832b27cdc33eb5874e83`, reached by ff-only
from `09d1b3c935919ec85f3c629607797ae87c42d8e5` in the existing production worktree.
Authority: ADR-0011, ADR-0013 and
`PHASE5C1R_STRUCTURAL_UNCERTAINTY_SEMANTICS_RATIFICATION.md` sections 9–14.
Change class: additive test-side measurement and validation documentation.

## Protocol fixed before measurement

M1 uses the governed Structure-1 whole-star midpoint, rho_c = 1.10e15 g/cm^3,
canonical 8192-interval table, and radial resolutions **20000, 40000, 80000**.
At each rung recompute the central star's A, all 15 sequence stars' counts and
B with log steps (0.001, 0.0005, 0.00025), then production fixed-baryon K.
The EOS source identity, achieved-abscissa derivative, Hartle equations,
whole-star reduction and PB13-comparison-v1 surface policy are unchanged.
Record all four raw species, A/B/reduction/K, inherited numerical error, and
the governed physical I conversion. No Richardson extrapolation is authorized.

For each consumed component, a contracting ladder without an independently
demonstrated order contributes the finest adjacent absolute K discrepancy.
A nonmonotone ladder requires independently supported floor characterization;
its stable floor may be retained under ratification section 12. A genuinely
unstable envelope stops chemical production. No retrospective safety factor
will be selected from production W.

M2 reports the existing PB10 reconstruction at the finest rung, using its
unchanged finite-current and independent count quadratures, q=(1e-7,5e-8),
and its unchanged acceptance inequalities. Delta_K = K_PB10 - K_prod;
relative discrepancy = abs(Delta_K)/abs(K_prod). Species order is n,p,e,mu.
The additive test translation unit imports the unchanged Phase-5B fixture;
the historical source and baseline producer remain byte-identical.

K and its errors/envelopes have units count km^2; physical I and its errors/
envelopes have units count s^2. K_error is NUMERICAL_ERROR. Validation
ingredients are empirical stability evidence, not certified continuum bounds,
confidence intervals, probabilities, or certified truncation errors.

No production G/Z/W output existed when this measurement protocol was fixed.
The following sections record the subsequent measurements and fixed envelopes;
chemical goals are in the separate production acceptance predeclaration.

## Raw measurements and classification

All raw A, B, A_B, B_B, central shift, K, K_error, I and E_I for all
12 species/rung combinations, plus the four raw PB10 comparisons, are retained
in `phase5c2_preproduction_evidence.json`. The freshly generated EOS has SHA-256
`7cd44c92e1e7206e0e68e3fed7e3f0ca68e79ab4517d02b96ff78b9be23d3f1a`.
The finest K and numerical errors reproduce the governed fixture values.

| Radial resolution | K_n | K_p | K_e | K_mu |
|---|---:|---:|---:|---:|
| 20000 | 1.17939603087605025e+58 | -1.17939603087605554e+58 | -1.04597773050993669e+58 | -1.33418300365856961e+57 |
| 40000 | 1.17938971789129435e+58 | -1.17938971789129755e+58 | -1.04597108532747557e+58 | -1.33418632563559813e+57 |
| 80000 | 1.17938973343378087e+58 | -1.17938973343378296e+58 | -1.04597114625517719e+58 | -1.33418587178127557e+57 |

All signed successive gaps reverse sign, but their magnitudes contract:
electron ratio 0.0091687025878886205, muon ratio 0.13662175224743806 (full precision in JSON).
This does not demonstrate a radial asymptotic order. The selected empirical
radial ingredient is the finest measured gap, with no Richardson extrapolation
or certified remainder claim. The observed contraction, rather than a fitted
floor or an assumed order, supports this finite-ladder stability envelope.

| Species | K_prod | K_PB10 | Delta_K | abs(Delta_K/K_prod) |
|---|---:|---:|---:|---:|
| n | 1.17938973343378087e+58 | 1.17938985147312800e+58 | 1.18039347127289654e+51 | 1.00085106543720137e-07 |
| p | -1.17938973343378296e+58 | -1.17938973327099101e+58 | 1.62791955457836279e+48 | 1.38030670305963174e-10 |
| e | -1.04597114625517719e+58 | -1.04597114606537625e+58 | 1.89800945486406679e+48 | 1.81459064302049565e-10 |
| mu | -1.33418587178127557e+57 | -1.33418587178503534e+57 | -3.75976626081477225e+45 | 2.81802284099673253e-12 |

PB10 acceptance is unchanged: each absolute discrepancy <= its production
K_error, independent charged reconstruction <= charge_budget, and the wrong
charge-map mutation remains separated. All passed; raw rc=0.

| Ingredient | Classification | Combination / limitation |
|---|---|---|
| K_error | NUMERICAL_ERROR | Convert once to E_I; excluded from V |
| PB11 q->0 direct K | VALIDATION_ENVELOPE_INGREDIENT | Last two independently reconstructed quotients; linear-q remainder demonstrated by PB11; abs(2 Q(q/2)-Q(q)-K) |
| PB11 fixed-baryon central shift / Delta N_B ~ q^2 | VALIDATION_ENVELOPE_INGREDIENT | Current task expressly requests inclusion beyond ratification's falsifier classification; retain conservative finite-q residual as described below |
| PB7 independent background B | VALIDATION_ENVELOPE_INGREDIENT | Componentwise ratio sensitivity, including B_B denominator |
| PB12 table ladder | VALIDATION_ENVELOPE_INGREDIENT | Finest adjacent K gap, no radial duplication |
| New radial ladder | VALIDATION_ENVELOPE_INGREDIENT | Finest adjacent K gap |
| PB6 knots/partition | VALIDATION_ENVELOPE_INGREDIENT | Transfer A and A_B sensitivity, then maximum with direct radial ingredient |
| PB10 raw K | VALIDATION_ENVELOPE_INGREDIENT | Absolute measured discrepancy |
| Prior review K | FALSIFIER / CORROBORATION | Not added as another sample of the same path |
| Charge/baryon identities and W=ZI | CONSISTENCY_ONLY | No envelope credit |
| PB13 tail | NUMERICAL_ERROR ingredient | Already inside K_error; not added to V |

The legacy raw evidence is `phase5b_resume_evidence.json` (PB6-repaired,
PB7-fresh, PB9-11, PB12); the controlling ratification establishes the
independence and qualifications. Historical provisional tail entries in that
legacy log are not imported: current K_error is freshly measured.

For PB7 let d_i=abs(B_prod,i-B_oracle,i), d_B=d_n+d_p. Transfer
abs(A_B) [d_i/(abs(B_B)-d_B) + abs(B_i)d_B/(abs(B_B)(abs(B_B)-d_B))].
For PB6 let a_i=abs(A_i)*recorded_knot_relative_i. Transfer
 a_i + abs(B_i)*(a_n+a_p)/abs(B_B).
For the PB11 central constraint retain the finite-q baryon quotient residual
R_B=abs(Delta N_B(q))/q at q=1.25e-7 km^-2,
Delta N_B=1.6812991431370884e48 count, and transfer abs(B_i/B_B)*R_B.
This deliberately retains the finite-spin nonlinear residual, rather than
claiming it is a q->0 linear-coefficient error. It is empirical validation
sensitivity only, and is the dominant conservative ingredient.

PB7 and PB11 central reduction can be dependent; no covariance cancellation
or RMS is claimed. The two PB11 observations share a current calculation but
constrain species reconstruction and the fixed-baryon condition respectively;
both are conservatively retained with this dependence explicit. The distinct
prior-review Richardson comparison is corroboration of the same PB11 path,
not another summand. PB6 radial and M1 are one representation class; use the
maximum of transferred knot sensitivity and direct K radial gap. Thus:

V_K = max(M1_radial, PB6_knots_transferred) + PB7_transferred
      + PB11_direct + PB11_constraint + PB12_table + abs(PB10_delta).

No K_error or tail term is added again; no statistical combination is used.

| Retained ingredient (count km^2) | e | mu |
|---|---:|---:|
| PB10_raw_absolute | 1.89800945486406679e+48 | 3.75976626081477225e+45 |
| PB11_Richardson_discrepancy | 9.98817276167329117e+50 | 3.88283318975645844e+50 |
| PB11_finite_q_baryon_constraint_transferred | 4.14973832364991122e+53 | 4.54319823859713402e+52 |
| PB12_table | 6.40332224250557062e+49 | 4.10577721877086919e+49 |
| PB6_knots_transferred | 1.17562980582457551e+50 | 1.31600281097907683e+49 |
| PB7_transferred | 2.20182553994228560e+52 | 2.96259795713440236e+51 |
| direct_radial | 6.09277016282893338e+50 | 4.53854322562759455e+50 |
| representation_max | 6.09277016282893338e+50 | 4.53854322562759455e+50 |

| Fixed quantity | e | mu |
|---|---:|---:|
| V_K_validation | 4.38666113288744123e+53 | 4.92777795165981210e+52 |
| V_I_validation | 4.88081875539544120e+42 | 5.48289241413407494e+41 |
| E_I_numerical | 2.08287122864888904e+40 | 9.66000638185162170e+39 |

Physical conversion uses AngularVelocity::FromRadPerSecond(1).GeomKmInverse()
squared, the same governed owner as WholeStarIPhysical(), exactly once.
The ledger's arithmetic value of c=299792.458 km/s is authenticated against
that owner; no new production unit constant or owner is introduced.

These are validation envelopes, not certified bounds. They apply to this
whole-star mathematical fixture, not a realistic A18 source target, confidence
interval, probability statement, or rigorous continuum/truncation error.

**NO PRODUCTION G/Z/W OUTPUT EXISTED WHEN THESE VALUES WERE FIXED.**
