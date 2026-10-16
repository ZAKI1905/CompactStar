# BA13 onset pilot: independent review of the failed predictions

| Item | Value |
|---|---|
| Candidate | `candidate/linux-ba13-onset-pilot-20261013` at `dad05b17056f3ccb1e3871e6e4bb24ebe5db5df4` |
| Evidence | `/mnt/sdd/zaki/CompactStar/linux-transition-20261010/ba13-onset-pilot` (read only) |
| Previous review | `ba13-uncertainty-independent-review-20261012/` (read only; still verifies) |
| Mode | Read-only scientific review. No source edit, commit, push, merge, baseline change, certificate field, production solve or BNV trajectory. |
| Calculations | Scripts in this directory; outputs `*.json`; hashes in `SHA256SUMS` (§11). |

Every calculation in this review is **post hoc**: it was run on retained data after the pilot's outcomes were known. These calculations explain mechanisms and can close an error model retrospectively. They do not validate an estimator, and nothing here may be presented as a prospective PASS.

---

## 0. Disposition

**I concur with disposition B: further bounded evidence is required. Nothing qualifies, and there is no 21-star expansion.** The pilot's own disposition is correct in outcome. Its explanation of both failure classes is incomplete, and one of the failed bands came from my own previous review.

1. **The P7 phase failure is quantitatively explained by an a-priori term.**
   - The linear-in-ν inverse interpolation carries an O(h²) error, −(ν″/2ν′)θ(1−θ)h², which is always outward. The coefficient depends only on M and r_e.
   - That term predicts the observed ranges: edge 1.30×10⁻⁹ (observed 1.24–1.30×10⁻⁹), L 3.16×10⁻⁹ (observed 3.00–3.13×10⁻⁹), and paired L 1.881/2.933×10⁻¹⁰ (observed 1.898/2.925×10⁻¹⁰).
   - The 10⁻⁹ / 10⁻¹⁰ bands were therefore set below a term that could be computed before the run. My previous review proposed those bands; the error is mine.
   - **The linear-ν estimator class is ruled out as phase-free.** A coefficient-free Hermite inverse interpolation that uses the exact TOV ν′ at the nodes removes the term by construction.
2. **The residual phase noise left after the correction is a common-B shooting effect.**
   - The production baryon-number quadrature has a grid-phase-dependent error. Common-B shooting carries it into ρ_c (about 10⁻¹⁰ variation across jitter), and from there into every quantity.
   - The same B functional also explains the "clean h²" radial convergence of R, M, N and ρ_c at fixed B, to 2–3% at every level from 20k to 640k.
   - My previous review attributed that convergence to TOV truncation; that attribution was wrong.
3. **The P5 EOS forecast was mathematically unsupported at the ~10⁻¹² scale.**
   - It extrapolated a single observed ratio two levels ahead, with no model for components that do not contract under EOS refinement.
   - The target-20 excess is localized: an outward pressure offset (4.6×10⁻¹⁰ in ln P) created where the integrator crosses the neutron-onset EOS kink. It is absent at targets 0 and 10, and absent in ρ_c, M and N_μ.
   - The common-B shooting residual (≤1.1×10⁻¹⁴) and binary64 representation are both ruled out.
4. **Frozen-relative drifts are robust and independently confirmed.**
   - The target-20 Ltilde_Me drift is +1.809209×10⁻⁶ from the Hermite-ν estimator (±4.5×10⁻¹² across six grid phases) and from the independent oracle (both tolerances).
   - The production value, −1.2266×10⁻⁶, is wrong because the edge is pinned to the grid.
   - None of the newly found terms comes within 10⁻⁴ of any frozen threshold. The largest is the known production pinning error, about 3% of the threshold.
5. **Absolute accuracy.**
   - The production error model, closed post hoc, agrees with the independent oracle to ≤3.4×10⁻¹² for r_e, R and M, and to ≤1.2×10⁻¹¹ for L.
   - C is the exception at 5×10⁻¹¹. That is a sharp, falsifiable target for the oracle's own convergence.
   - The oracle is an independent formulation but is **not** a certified reference: it has two tolerances only and an unregularized endpoint singularity at the onset.

**Recommendations**

| Question | Recommendation |
|---|---|
| ADR-0021 | Option A remains appropriate. Interpretation (2) with the right side frozen remains the recommended reading, pending an owner ruling. |
| Production amendments | Warranted as separately governed changes, in priority order: thermal cache law, then neutron-support edge, then envelope base, then onset-aware Cv quadrature. They are not prerequisites for frozen-relative qualification of historical BA13. |
| 0.5 CPU-hour continuation | **NO-GO** as proposed. |
| Next experiment | One zero-production-solve oracle-anchor experiment, about 2 CPU-minutes, with frozen numerical predictions (§9). |

---

## 1. Authenticated evidence summary

Verified by [`authenticate_pilot.py`](authenticate_pilot.py) and recorded in `authentication.json`:

- **Candidate state.** HEAD is `dad05b1` and the worktree is clean.
- **Frozen files, unchanged against their commits:**

  | File | Commit | SHA-256 (prefix) |
  |---|---|---|
  | PREDECLARATION.md | `8bb72e2` (14:58:00) | `63ccc49d…` |
  | predictions.json | `e657efb` (15:00:09) | `1803bc8d…` |
  | eos-predictions.json | `650871e` (15:04:41) | `6b9dbc08…` |

- **Gate timing.** The first controller output (`table-2048-11.log`) is 2 s after the final gate commit `650871e`. No estimator, runner or controller file changed after `650871e`. The later changes were limited to:
  - `score.py`: added metrics and one stricter extra envelope test; no frozen band changed;
  - `publish.py`;
  - `diagnose.py` and `geometry_diagnostic.py`, which are post-hoc.
- **Receipts.** All 31 run receipts carry one executable (`e574d959…`, equal to the retained binary) and the frozen predeclaration and prediction hashes. All are `COMPLETE_NOT_QUALIFIED`.
  - *Governance note:* the receipts record four different HEADs (3 / 3 / 22 / 3 runs), because commits continued in the same worktree during the fixed runs. Executed code was unaffected, but future pilots should run from a detached, pinned worktree and hash each executed script.
- **Manifest.** `PILOT_COMPLETE.json` (SHA-256 `3d508ffd…`) lists 650 files, 1,130,004,091 bytes, with zero mismatches. Source and protected-file hashes match the worktree. Compact repository copies are byte-identical. The bundle hash `a8884a62…` matches the record. The previous review directory still verifies.

**Classes of evidence**

| Class | Members |
|---|---|
| Prospective | The 31 new stars and 6 oracle solves, scored by P1, P2, P3 (640k and jitter), P5, P6 and P7 |
| Discovery | My previous ν ladders (20k–320k at targets 0, 10, 20); the old e2048, e4096 and nominal profiles; the nominal P3 controls |
| Secondary | P4 (targets 5 and 15 on retained profiles) |
| Post hoc | `diagnose.py`, `geometry_diagnostic.py`, and everything in this review |

**Failures retained exactly as recorded** (observed / band):

| Test | Result |
|---|---|
| P5, 32768c holdout, target 20 | L 1.037, C 1.138: **FAIL** (r_e 0.974 passes) |
| P5, 65536c consistency, target 20 | L 6.47, C 6.95, r_e 6.08: **FAIL** |
| P7, absolute | L 3.00 / 3.13 / 3.09; edge 1.24 / 1.30 / 1.28: **FAIL** |
| P7, paired | L 1.90 / 2.92; edge 1.21 (target 20): **FAIL** |
| P7, C | passes (absolute 0.04–0.05, paired 0.07–0.10) |

**An internally inconsistent predeclared rule, applied in degenerate form without disclosure.**
- P5 specified a three-level EOS discovery basis (2048, 4096, 8192) but invoked the P1 eligibility rule, which needs three increments and two ratios. Three levels give only two increments and **one** ratio (q ≈ 0.0605, `eos-predictions.json`). The code silently evaluated the rule on what was available.
- The very broad bands are labelled honestly, but their observed fractions show how little they constrain: P2 stability (10⁻⁷) was met at ≤3×10⁻⁴ of the band, P4 at ≤4×10⁻³, and P6 at ≤10⁻³. P3, a mechanism-falsification test, passed at ≤0.14% of its 2% band.

---

## 2. Resolution-phase failure (P7)

### 2.1 Derivation

Let ν be smooth on the cell [r_a, r_b] = [r_a, r_a + h]. The true onset r* satisfies ν(r*) = ν_on, and θ = (r* − r_a)/h. Linear interpolation gives

ν(r) − ν_lin(r) = ½ ν″(ξ)(r − r_a)(r − r_b).

Solving ν_lin(r_lin) = ν_on and Taylor-expanding gives

**r_lin − r\* = −(ν″/2ν′) θ(1−θ) h² + θ(1−θ)(1−2θ) h³ [(ν″/ν′)²/4 − ν‴/(6ν′)] + O(h⁴).**

At the neutron onset, 4πr³P/m ≈ 6×10⁻¹¹ and the soft-layer ε contributes only about 7×10⁻⁷ of ν″. To that accuracy, ν′ = M/D with D = r(r − 2M), so that:
- ν″/(2ν′) = −(r − M)/D;
- the bracketed h³ coefficient is −M²/(3D²) exactly.

**Applicability and sign.** ν is C¹. The soft layer changes ν″ only at the ~10⁻⁶ relative level. Because ν″ < 0 and ν′ > 0, the linear edge is always **outward** of the true edge, so L_lin ≥ L_true, with the maximum bias at θ = ½.

### 2.2 Quantitative agreement (`phase_model.json`)

All three targets give ν″/(2ν′) = −0.085837 km⁻¹ and ν′ = 6.7435×10⁻³ km⁻¹, matching the pilot's own TOV-based evaluation. The predicted maxima are:

- edge: (r − M)h²/(4Dr) = 1.2995×10⁻⁹;
- L: s·(r − M)h²/(4D) = 3.159×10⁻⁹, with s = ∂ln L/∂r_e = 0.1923 km⁻¹.

| Quantity | A-priori model | Observed (six grid phases) |
|---|---|---|
| Edge range, max over θ | 1.30×10⁻⁹ | 1.24 / 1.30 / 1.28×10⁻⁹ |
| L range, max over θ | 3.16×10⁻⁹ | 3.00 / 3.13 / 3.09×10⁻⁹ |
| Paired L range, target 10 | 1.881×10⁻¹⁰ (per-case θ) | 1.898×10⁻¹⁰ |
| Paired L range, target 20 | 2.933×10⁻¹⁰ | 2.925×10⁻¹⁰ |
| Linear-minus-Hermite edge difference vs formula | identity | ≤4.1×10⁻⁵ relative |

The 640k and r_e withheld ratios (0.273–0.292), which the pilot noted depart from ¼, are this phase term. With a phase-free edge, the 20k–640k ladders become **0.2485, 0.2509, 0.2497, 0.2501** at all three targets (`hermite_ladder.json`).

### 2.3 The residual is common-B shooting noise, not a curvature remainder

After the leading term is removed, the edge still varies by 2.5–2.8×10⁻¹¹ and L by 3.7–4.3×10⁻¹¹ across jitter. The h³ term is at most 10⁻¹⁵ km, so it cannot account for this. The cause is the following:

- **ρ_c moves with the grid.** Across the six jitter grids, ρ_c varies by up to 1.1×10⁻¹⁰ and N_μ by 1.2×10⁻⁹. That is about 20 times the smooth h² change from a 6×10⁻⁵ step change.
- **The B functional drives it.** The production B quadrature has a phase-dependent error relative to a phase-free ν-coordinate B. That error predicts the large ρ_c excursions to about 10%, through dlnB/dlnρ_c = 0.1805 (`b_functional.json`). For example:

  | Case | Predicted ρ_c shift | Observed |
  |---|---|---|
  | r80003, target 0 | 7.75×10⁻¹¹ | 7.20×10⁻¹¹ |
  | r80004, target 0 | 8.62×10⁻¹¹ | 7.99×10⁻¹¹ |
  | r80003, target 20 | 9.46×10⁻¹¹ | 8.65×10⁻¹¹ |

  Smaller excursions agree only to about ±4×10⁻¹¹, so at least one further component remains.
- **Removing it closes the residual.** Subtracting dlnL/dlnρ_c · δρ_c (−0.3265, a sequence derivative, not fitted) leaves 0.3–9.3×10⁻¹² in L. In paired drifts the ρ_c noise is nearly common, and the Hermite paired range is only 6.1 / 8.9×10⁻¹².

### 2.4 Validity and limits of the curvature correction

1. **It is not validation.** It was derived after the failure.
2. **Its coefficient is independent of the approximate profile.** It needs only M and r_e at the 10⁻⁹ level, and both are available from the independent oracle. The phase θ needs the true edge, but iterating with the corrected edge converges, and an error δθ ~ 10⁻⁵ changes the correction by about 10⁻¹⁴.
3. **The omitted terms are bounded.** By Taylor's theorem with explicit TOV derivative bounds, the h³ term is at most 1×10⁻¹⁵ km and the soft-layer ν″ jump is at most 7×10⁻⁷ of the leading term. Both are far below every other component.
4. **It leaves the B-functional ρ_c noise untouched** (~3×10⁻¹¹ absolute, ~10⁻¹¹ paired). Claiming "phase-free to 10⁻¹¹ absolute" would need a phase-free B functional as well.
5. **A better method exists.** Hermite inverse interpolation with the exact nodal ν′ = (m + 4πr³P)/(r(r − 2m)) has no fitted coefficient and is accurate to O(h⁴). It reproduces the corrected values exactly (`hermite_edge` columns in `phase_model.json`). It is the estimator to adopt; another jitter sweep would add nothing.

---

## 3. EOS-refinement failure (P5)

### 3.1 Separation of components (`eos_decomposition.json`)

| Component | Evidence | Size at 32768c / 65536c |
|---|---|---|
| EOS tabulation and Steffen interpolation | M, N_μ and ρ_c increments are identical across the three targets to 2–3 digits. L, C, r_e and R contract by about 1/12–1/16 per doubling at targets 0 and 10. N_μ and ρ_c contract irregularly (0.05–0.16, with ρ_c alternating in sign), but identically at every target. M reaches about 2×10⁻¹⁴ | L, C, r_e: ~3×10⁻¹² then ~1.7×10⁻¹³ (targets 0, 10) |
| Binary64 representation | Capped tables keep ≥56 ulp gaps in P and ≥72 in e at the audited near-onset pairs. Depth 12 vs 11 changes **no output bit**. Historical refused pairs have true gaps of 5.1–7.2 ulp (P) and 6.4–9.2 ulp (e), consistent with my previous review | not implicated |
| TOV solution error at the onset crossing | δln P(32768c→65536c) at fixed nodes is ~10⁻¹¹ in the interior and ~5×10⁻¹⁰ just inside the onset. It leaves an **outward offset** of +2.4×10⁻¹¹, −1.6×10⁻¹¹ and **+4.57×10⁻¹⁰** (targets 0, 10, 20), giving R changes of 1.2×10⁻¹³, −2.8×10⁻¹⁴ and **1.36×10⁻¹²** | target-20 excess ≈1–3×10⁻¹² in R, ν_s, r_e, L, C |
| Central-density / common-B shooting | Maximum \|residual\|/B is 1.08×10⁻¹⁴ over all new stars. With sensitivity ~1.8 that is ≤2×10⁻¹⁴. ρ_c increments are target-independent | ruled out |
| Neutron-onset position (estimator) | At a fixed radial grid, θ changes by ~10⁻⁸ between EOS levels, so the curvature term cancels | negligible |
| Paired-difference cancellation | Not used by P5, which scores absolute values. Paired EOS effects are ≤6×10⁻¹³ except the target-20 excess (~3×10⁻¹²) | irrelevant at 10⁻⁴ thresholds |

**Most probable mechanism for the target-20 excess.**
- The production TOV integrates adaptively (rk8pd, relative tolerance 10⁻¹⁰) between output nodes, with no event location at the onset pressure P_on.
- At P_on the derivative dε/dP jumps by about 100×. The crossing error is therefore tolerance-limited and erratic. It changes whenever ε(P) near P_on is perturbed, and table refinement is such a perturbation.
- This is supported by three observations: the perturbation is localized to the onset, it differs by target, and it leaves ρ_c, M and N_μ unchanged. It is not proven. A decisive test would need event-located integration, which is a production-solver change and should not be made merely to pass P5.

### 3.2 Was the forecast supported at the 10⁻¹² scale?

**No.** Three reasons:
- It rested on one ratio from three levels, below the declared eligibility rule.
- It assumed pure geometric contraction for two more levels, into a band of 4.1×10⁻¹² at 32768c and about 5×10⁻¹³ at 65536c, with no model for the non-contracting component.
- The only floor checked was 128ε ≈ 2.8×10⁻¹⁴, two orders below the onset-crossing component the same pilot realized.

The passes at targets 0 and 10 for 65536c reflect small realizations of that component, not support for the model. I propose no numerical floor. The defensible responses are either to remove the component by method (event-located crossing, or an enthalpy-coordinate structure solver) or to measure it independently before forecasting.

### 3.3 Is the capped family a convergent approximation to the same EOS?

**Yes, to binary64 resolution.**
- Nodes are exact provider values.
- The base grid and the cluster density refine together.
- The only interval that is never refined is |n/t − 1| < 10⁻¹¹. It lies in the stiff regime (Γ ≈ 1.41), where interpolation error is about (10⁻¹¹)² ≈ 10⁻²² relative.
- The depth-12 experiment confirms this: it changed no output bit at any target.

The historical depth-13 refusal stays a refusal; it is a representational artefact of that grid family.

---

## 4. Independent solver

**Formulation.** The oracle (`onset_structure_oracle.cpp`, using `tests/eos/structure1/tov_oracle.hpp`) works as follows:
- the independent variable is the enthalpy H, with state variables (r², m/r³);
- the EOS comes from momentum integrals, with no table and no Steffen spline;
- integration is RKF45, stopping exactly at H_μ, H_b, H_n and H_s;
- it runs its own common-B secant, aligned to 5×10⁻¹³ at tolerance 2×10⁻¹².

It shares with production only the constants, units, surface pressure and B targets. **It is an independent mathematical formulation, not a certified reference.**

**Tolerance sensitivity** (relative change from tolerance 2×10⁻¹⁰ to 2×10⁻¹²):

| L | C | M | R | r_e | μ∞ |
|---|---|---|---|---|---|
| −7.05×10⁻¹⁰ | −2.59×10⁻⁹ | −3.32×10⁻¹⁰ | +9.99×10⁻¹¹ | +9.50×10⁻¹¹ | +3.64×10⁻¹¹ |

These changes are identical across targets to four digits, so they are common-mode, and the oracle's paired drifts agree to seven digits at both tolerances. C's larger sensitivity is consistent with the unregularized √(H − H_n) endpoint singularity in its integrand. Two tolerances cannot establish a convergence order.

**Sign of the target-20 Ltilde_Me response.** Three determinations agree:

| Source | Value |
|---|---|
| Hermite-ν estimator, six grid phases | +1.809209×10⁻⁶ ± 4.5×10⁻¹² |
| Independent oracle, both tolerances | +1.809209×10⁻⁶ |
| Previous review's ν value | +1.809229×10⁻⁶ |
| Production (pinned edge) | −1.2266×10⁻⁶ |

The positive sign is physical: the outward edge motion (+1.59 cm, worth +3.06×10⁻⁶) outweighs the change in the e^(−6ν) integrand. The P3 pinning identity explains the production error to 2.8×10⁻⁵ relative.

**Absolute closure, post hoc** (`anchor.json`). The production limit X_inf is built from the Hermite edge, Richardson extrapolation with p = 2 fixed, and the measured EOS-table offset. The existing oracle at 2×10⁻¹² differs from it by:

| r_e | R | M | L | C |
|---|---|---|---|---|
| −0.7 to −2.7×10⁻¹² | +0.6 to +2.1×10⁻¹² | +3.3×10⁻¹² | +0.7 to +1.2×10⁻¹¹ | +4.7 to +5.3×10⁻¹¹ |

**Required for a defensible absolute bound:**
1. at least four oracle tolerances with monotone contraction;
2. regularized onset segments (s = √(H − H_on));
3. compensated accumulation;
4. prospective agreement with X_inf inside model-derived bands.

That is the next experiment (§9).

---

## 5. Thermal and envelope consequences

| Quantity | Absolute error (nominal) | Frozen-relative (paired, target 20) | Cancels in equal-knot ratios? |
|---|---|---|---|
| Cstar cache law F(T) | +3.30 to +6.56×10⁻⁴; log-slope 0.965–1.037 | 0 exactly | yes |
| Cstar spatial (nodal trapezoid vs ν) | +1.71 to 1.82×10⁻⁷ | −1.10×10⁻⁸ | no |
| Ltilde_Me, production edge | +9.29×10⁻⁵ (target 0); +8.98×10⁻⁵ (target 20) | −3.04×10⁻⁶ (wrong sign) | no |
| Envelope base T_b (node snapping) | +4.12 to 4.23×10⁻⁶ | +1.07×10⁻⁷ (57% of the drift; invisible to refinement) | no |
| T_s, T_s∞ (iron, slope 0.605) | +2.49 to 2.56×10⁻⁶ | +6.48×10⁻⁸ | no |
| Photon luminosity ∝ T_s⁴ | ≈ +1.0×10⁻⁵ | ≈ +2.6×10⁻⁷ | no |
| B functional, via common-B shooting | B +7.5×10⁻⁹, ρ_c −4.2×10⁻⁸; propagates to all fixed-B quantities | ≤10⁻¹⁰ (near-common) | no |

**Can historical BA13 remain a bounded approximation?** Yes, in substance, for frozen-relative use. Every paired discrepancy is mechanistically explained and independently checked. The largest is the Ltilde_Me pinning error, 3.04×10⁻⁶ (3% of 10⁻⁴); it raises Ltilde_Me utilization from 0.012 to about 0.043 and changes no decision. What is still missing is an enclosure *construction* (§6), not new resolution data.

**Should production accuracy be improved independently?** Yes, separately, for different reasons:
- For physical thermal predictions, every absolute bias above is far below envelope, EOS and superfluid modelling uncertainty. The production rules are therefore *physically* adequate.
- They are **numerically inadequate** in three respects: as accuracy claims (the governed errors are 0 or 64ε), for derivative and sensitivity quantities (dL/dB has the wrong sign; the T_b drift is 57% wrong), and as grid-independent definitions.

These two answers must not be combined: improving production is neither required for, nor a substitute for, the frozen-relative enclosure.

---

## 6. ADR-0021 and governance

**Option A remains appropriate.**
- The historical fitter is not the error source.
- Every newly found term enters through reconstruction rules or common-B shooting, and all of them can be enclosed by independent analysis.

**Interpretation (2), frozen right side.** The left side is a joint, correlated, paired enclosure plus the nonlinear and inter-sample remainder. The right side is the historical threshold, evaluated from the frozen governed u. This remains the recommended reading. The pilot sharpens it in three ways:
- the enclosure must be built from mechanism-level models (pinning identity, a-priori curvature, B-functional shooting, onset-crossing allowance), anchored by an independent reference;
- bands inferred from refinement differences are unreliable wherever phase, pinning or kink terms exist;
- the frozen right-hand allowances (≤1×10⁻¹³) stay frozen.

**Unresolved mathematical questions**
- The inter-sample and inner-half remainder.
- Paired Z enclosure including the resolvent Jacobian term: Ĝ⁻¹ − G⁻¹ = −Ĝ⁻¹EG⁻¹. The pilot's correction of my previous §8.5 is right. The term is about 10⁻⁸ here, but it must be in the construction.
- 21-target coverage of the phase and pinning terms. These are deterministic and computable from the existing nominal profiles.
- Tangent-stencil, background and tail paired budgets.
- An absolute reference certification, needed only for absolute claims.

**Unresolved governance questions**
- An owner ruling on (2)-frozen.
- Relabelling the governed errors (0, 64ε, |order16 − order32|) as "not error bounds" in any successor record.
- Whether a successor artifact or runtime monitor should carry enclosed bounds. This is not needed if the enclosed utilization stays ≤1 over the whole window.
- Nomenclature. The pilot's "C" means a separate *production* proposal. It is not ADR-0021 Option C, which is the fitter amendment, and nothing here triggers that.

No ratification is implied by this review.

---

## 7. Decision on production-method amendments

These are separately governed proposals. **Do not implement them under BA13, and do not make them prerequisites for frozen-relative qualification.**

| Rule | Decision | Reason | Method to propose |
|---|---|---|---|
| Thermal cache law | **Warranted; highest priority by magnitude** | +6.6×10⁻⁴ one-signed bias, ±3.6% log-slope error, 160 slope kinks in the right-hand side | Interpolate ln C against ln T (exact for C ∝ T), or cache C/T |
| Neutron-support edge | **Warranted** | Wrong-sign drift; O(h) bias of 9.3×10⁻⁵; pinned to the grid | Equilibrium-coordinate edge with Hermite inverse interpolation using the TOV ν′, or support integrals in enthalpy |
| Envelope base | **Warranted; low priority** | O(h) snapping (4×10⁻⁶ in T_b); a paired error that refinement cannot see | Interpolate ν at ε_b, or use the identity T_b = T∞ μ_b/μ∞ |
| Onset-aware Cv quadrature | **Warranted only inside a thermal-owner amendment** | 1.8×10⁻⁷, well below the cache bias | Onset-split quadrature in the equilibrium coordinate |
| Baryon-number functional | **Awareness item** | Sets the fixed-B targeting bias (ρ_c −4.2×10⁻⁸) and ρ_c phase noise (~10⁻¹⁰) | Integrate N alongside the TOV ODE |

---

## 8. Assessment of the proposed 0.5 CPU-hour continuation: **NO-GO**

| Element | Scientific value | Reason |
|---|---|---|
| Jitter 80006–80011 (18 stars) | Very low | The phase dependence is already explained a priori (curvature) and by the B functional. The linear-ν estimator is already ruled out, and the Hermite method needs no coefficient to validate. More phases would resample known deterministic functions. |
| Three 1,280,000-point profiles | Low | The phase-free ladders already show 0.2500 ± 0.0005 through 640k. Absolute accuracy is better addressed by the oracle. |
| Whole-table 32768c high-precision audit | Low for this question | Binary64 representation is implicated in neither failure. Postpone it to an absolute-EOS enclosure, using interval arithmetic if that is ever needed. |
| Independent solver checks | **High, if redesigned** | This is the only missing component (an absolute reference). As proposed it lacks regularization, a tolerance ladder and frozen predictions. |

The existing data rule out the proposed estimator class (linear-ν edge, refinement-derived bands). The response is a different method (Hermite edge plus a certified oracle anchor), not another resolution sweep. Three-target success would not justify 21-star expansion.

---

## 9. Smallest next numerical experiment: X1 "oracle anchor", zero production solves

**Purpose.** Certify the independent enthalpy oracle as an absolute reference and test, prospectively, the production absolute error model:

X_inf = Hermite-edge value + p = 2 Richardson tail + EOS-table offset.

**Runs.** All three targets (0, 10, 20) with the retained B targets and surface pressure:
- the unchanged oracle at new tolerances 2×10⁻¹¹ and 2×10⁻¹³ (6 solves);
- a validation-only regularized oracle at 2×10⁻¹⁰, 2×10⁻¹¹, 2×10⁻¹² and 2×10⁻¹³ (12 solves). It differs only in integrating the neutron- and muon-onset segments in s = √(H − H_on) and accumulating B, L and C with compensated summation. It must carry a distinct identity.

**Freeze before any run.** Commit `anchor.json` from this directory, or regenerate it bit-identically, as the prediction record.

### Error models (a priori)

| Model | Statement |
|---|---|
| Regularized oracle | RKF45 global error ∝ tol^α with α ∈ [0.7, 1]. Prediction: each decade of tolerance contracts \|ΔX\| by at least 3.3× for L, C, R, M and r_e. |
| Production limit | Band B_X = q-term + \|last EOS increment\| + ρ_c phase at 640k + onset-crossing allowance (10×rtol × H_P = 4.5×10⁻¹¹ km). No component uses the oracle-versus-prediction differences. **The onset-crossing allowance is the weakest-founded component.** It is a heuristic safety factor of 10× over the integrator tolerance, sized with knowledge of the single observed realization (4.57×10⁻¹⁰ in ln P at target 20). It is not derived from a crossing-error theorem. |

Frozen bands B_X (relative):

| L | C | r_e | R | M |
|---|---|---|---|---|
| 1.32 / 1.34 / 1.64×10⁻¹¹ | 3.7 / 4.0 / 7.3×10⁻¹² | 6.7 / 6.8 / 8.1×10⁻¹² | 6.9 / 6.9 / 8.1×10⁻¹² | 2.4 / 2.4 / 2.5×10⁻¹² |

### Frozen acceptance criteria

- **A1, oracle certification.** For every X and target, the regularized ladder contracts by at least 3.3× per decade over the last two decades. U_oracle(X) = \|ΔX(2e-12 → 2e-13)\|·ρ/(1 − ρ), where ρ is the observed last ratio (≤0.3).
  - If A1 fails, the outcome is **INCONCLUSIVE** (oracle not certified) and A2 is not scored.
- **A2, absolute anchor.** A2 is scored for a quantity X **only if U_oracle(X) ≤ B_X**. Otherwise the outcome for X is **INCONCLUSIVE**, never PASS, because a large oracle uncertainty must not be able to absorb the test. When scored, \|X_oracle(2e-13, regularized)/X_inf − 1\| ≤ B_X + U_oracle gives **PASS**; anything larger is **FAIL** (a missing common systematic).
  - Sharp prediction: oracle C must fall by about 4–5×10⁻¹¹ from its 2×10⁻¹² value, the direction its own tolerance trend already shows. If it stays near its 2×10⁻¹² value, a C component is missing.
  - Likewise, oracle M must move toward X_inf by about 1×10⁻¹², or else carry U_oracle(M) of that size. Oracle L must end within B_L + U_oracle of X_inf.
  - *Expected decisiveness:* at the minimum allowed contraction (3.3× per decade), the retained 2.59×10⁻⁹ two-decade C change implies a last increment of about 1.8×10⁻¹⁰ and U_oracle(C) ≈ 8×10⁻¹¹ > B_C. C would then be INCONCLUSIVE. Making C decisive needs a 2×10⁻¹⁴ level or a fixed-step, known-order cross-method (for example Gauss–Legendre collocation in s with step halving). A C FAIL would be informative but not alarming, since frozen validity is unaffected.
- **A2s, sign prediction (frozen separately, independent of magnitude).** From the 2×10⁻¹² value to the regularized 2×10⁻¹³ value, oracle C, L and M must **decrease** and oracle r_e must **increase**, toward X_inf. For these four quantities the existing difference and the oracle's own tolerance trend point the same way. R is excluded from A2s: its difference is already within B_R and its tolerance trend points away from X_inf, so it is scored by A2 only. A move in the wrong direction is **FAIL** for that quantity, whatever the magnitude.
- **A3, paired invariance.** Every oracle variant's paired drifts (L, C, r_e, R, M, μ∞ at targets 10 and 20) agree with the retained 2×10⁻¹² oracle to ≤1×10⁻¹². They agree with the Hermite-ν production paired drift (nominal, all six phases) to ≤1×10⁻¹¹. That band is model-derived: ρ_c paired ≤5×10⁻¹² plus the onset-crossing allowance.

### Conditions and limits

**Independence checks (all required)**
- The regularized and unregularized oracles at 2×10⁻¹³ agree within 10·U_oracle.
- B alignment ≤5×10⁻¹³ for every solve.
- No production table, profile or spline is read by the oracle (verified by code inspection and hashes).

**INCONCLUSIVE if:**
- A1 fails;
- any solve fails or cannot be aligned;
- regularized and unregularized results disagree beyond 10·U_oracle;
- the resource ceiling is reached.

**Early stop:**
- stop at the first solve exceeding 60 CPU-s or failing alignment;
- stop scoring after an INCONCLUSIVE A1;
- add no new tolerances or targets after any outcome.

**Resources.** About 0.1 CPU-s per measured oracle solve and at most ~10× that for the regularized variant, so 18 solves need ≤20 CPU-s. Scoring needs ≤60 CPU-s. **Ceiling: 120 CPU-s** numerical, Slurm, one core, ≤1 GiB. No production solve, EOS sweep or 21-target run.

**What a result would mean**

| Result | Consequence |
|---|---|
| PASS | A certified absolute reference exists at about 10⁻¹¹ for three targets, and the production absolute error model is prospectively confirmed. |
| FAIL | A common component of at least 10⁻¹¹ is missing. This is relevant to absolute claims only; frozen validity is unaffected. |

**Before any 21-target step:** the owner rules on (2)-frozen, X1 passes, and a predeclared enclosure construction exists. The step itself would be analysis plus oracle runs on existing nominal profiles, still with no production solves.

---

## 10. Corrections to my previous review (2026-10-12)

1. The P7 bands I proposed (10⁻⁹ absolute, 10⁻¹⁰ paired) ignored the ν estimator's own a-priori O(h²) curvature term. They were mathematically unsupported, and P7's failure belongs to them.
2. I described the ν-coordinate edge as phase-free. It is phase-free only with curvature-aware (Hermite) inverse interpolation.
3. I attributed the clean radial convergence of R, M, N and ρ_c to TOV truncation. It is predominantly the B-functional quadrature error propagated through common-B shooting.
4. §8.5 ("apply the majorant to a paired E_G") was incomplete. The paired resolvent term must be included, as the pilot's EQUATIONS_AND_LIMITS states.

My frozen-relative conclusions are unchanged. The four errors above concern error attribution and band design, not any utilization outcome.

---

## 11. Artifacts and reproduction

Run from this directory: `/mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain/python/bin/python -I -B <script>`. The scripts read the evidence trees read-only and write only here. The measured cost is below 70 CPU-s in total, mostly B quadratures.

| Script | Output | Content |
|---|---|---|
| [`authenticate_pilot.py`](authenticate_pilot.py) | `authentication.json` | Hashes, frozen files, receipts, timing |
| [`phase_model.py`](phase_model.py) | `phase_model.json` | Curvature model, Hermite edge, ρ_c attribution, paired comparison |
| [`hermite_ladder.py`](hermite_ladder.py) | `hermite_ladder.json` | 20k–640k ladders, linear vs Hermite |
| [`b_functional.py`](b_functional.py) | `b_functional.json` | B quadrature error vs ρ_c ladder and jitter |
| [`eos_decomposition.py`](eos_decomposition.py) | `eos_decomposition.json` | EOS increments, shooting residuals, onset localization |
| [`anchor.py`](anchor.py), then `anchor.py bands` | `anchor.json` | Absolute closure and frozen-candidate X1 predictions and bands |
| [`common.py`](common.py) | — | Shared read-only helpers |
