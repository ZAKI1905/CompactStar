# CompactStar — Independent Scientific Review of the BA13 Final Decision Package

| | |
|---|---|
| Candidate | `bb0f9c5af8ebaf86b5a59aece5551bb9dafbd264`, branch `candidate/linux-ba13-final-decision-20261030`, worktree `worktrees/ba13-final-decision`. Clean before and after this review, with 0 ignored files. |
| Entry / predeclaration | Entry `84e1062` (2026-10-26 evidence package). Predeclaration `d5b40b3`. Assessment frozen at `5957154`; controls-v2 at `328ccf4`. |
| Primary documents | `docs/validation/ba13-final-decision/REPORT.md`, `OWNER_DECISION.md` |
| Also read | `METHODS.md`, `MATRIX.md`, `PREDECLARATION.md`, the 462-row matrix, ADR-0021, preflight §11 (`PHASE6A1_CONTROLLED_BNV_IMPLEMENTATION_PREFLIGHT.md:655`–`762`) with §§8–9, the source and evidence manifests, the conditional-interval and tangent models, the 2026-10-26 evidence (holdouts, diagnostics, continuation), the X2 report and ablations, and the 2026-10-20 independent closure review |
| Review directory | `/mnt/sdd/zaki/CompactStar/linux-transition-20261010/ba13-final-decision-independent-review-20261030/`. It is a sibling of the evidence trees and outside every Git repository. Contents are listed in `SHA256SUMS`. |
| Mode | Read-only. No source edits, commits, pushes or merges. No solver runs, stellar evaluations or trajectories. No certificate promotion and no ratification of any owner decision. Only static Python on retained files (toolchain Python 3.11.9 with `-I -B`) and a re-run of the existing unit tests. |
| Date | 2026-10-30 |

This is the first independent review of the 2026-10-26 evidence package and of the 2026-10-30 decision package. The previous independent review (2026-10-20) covered `093d121`.

---

## 0. Verdict

**Authentication is clean.**
- 120 of 120 independent checks pass (§2).
- All 55 unit tests pass on re-run.
- All 1,993 inherited manifest entries and all 1,386 historical drift/threshold/utilization fields verify value-for-value against the immutable certificate.
- The original certificate, fitter, production code and dependency locks are unchanged.
- Every earlier FAIL, INCONCLUSIVE and UNQUALIFIED record is preserved.

**The package's disposition E is correct under its own predeclared rules.** I do not rescore it. My assessment in the mission's vocabulary differs for one reason:

> **The six remainder types that remain null in every approximate row are missing *declared error models*, not missing *evidence*.** Existing retained evidence supports each of them as either an independently validated numerical estimate (class B) or a conditional estimate with explicit quantitative content (class C). None needs a new stellar calculation.

The package treated "no declared model" as "no estimate". In two cases it applied a standard stricter than the one it proposes ("uniform stage enclosure"). That is correct procedure for an agent forbidden to choose an error model on the owner's behalf. It is not evidence that the estimates are unsupported.

**Answer to the central question — yes.** An adequately justified frozen-background approximation for the controlled BNV experiment can be established without a mathematical proof of every numerical operation. Three facts support this:

1. **The criterion has drift structure.** §11 bounds the change of each frozen quantity across a 10⁻⁶ window. A smooth error shared by both implementations enters only through its derivative with respect to B. For it to matter in the muon number, two independent codes would have to share an error of ≥ 85 % in the B-sensitivity d ln N_μ/dx itself. That is a model-specification failure, not a numerical-accuracy question (§3.1).
2. **Implementation independence is real.** The X2 reference includes only the authenticated physical-constants header (verified in source, §3.2). Three implementations (production, X1 and X2) agree far below the required accuracies.
3. **The inter-sample model has been tested prospectively, quantitatively and successfully.**
   - The three-target curvature predicts both withheld holdout discrepancies within a factor 0.64–1.61 in every group, and within 1.0–6.4 % in the muon sector.
   - What remains after curvature is one common abscissa offset, |δx| = 1.5–1.7×10⁻¹⁴. It is the same in 14 of 15 groups whose slopes span a factor of 290, and it equals the recorded X2 B-alignment residual, 1.6–1.9×10⁻¹⁴ (§3.7).

**Disposition: B — ready for a narrower, explicitly conditional acceptance.** It is narrower than full BA13 in three explicit ways:
- The wording is "qualified under named validated numerical estimates and itemized conditionals", not "rigorously enclosed".
- The claim is frozen-relative only, not absolute thermal-evolution accuracy.
- `metric_structure` is covered through its certified consumer functionals, not by a nodewise reference. This needs a consumer-inventory confirmation (§6).

Full BA13 qualification (A) additionally needs one bounded calculation: the passive nodewise metric observer (§6). If the owner instead accepts the consumer route for the metric row, the claim stays at B as an explicit narrowing of the §11 table, however well justified.

**Smallest justified next action: owner adjudication.**
1. The owner decides Option A, interpretation (2) and A–E, with the three amendments in §4.
2. The owner accepts or rejects each of the seven quantitative conditions stated in Appendix A.
3. A single **mechanical, static sidecar completion** then follows. It writes the six term models into the rows, commits the registry separately from the assessor, and re-runs the assessor and tests: ≤ 60 CPU s and no stellar evaluation.
4. A short acceptance check follows.

No further solver campaign, holdout or refinement is warranted.

---

## 1. Findings at a glance

### 1.1 Confirmed

| Item | Result | Where |
|---|---|---|
| Identity, history, manifests, bundle, archive, receipts, certificate, fitter, locks, protected trees | 120/120 checks pass | §2, `outputs/authentication.json` |
| Tests | 55/55 pass on re-run (3 + 14 + 13 + 13 + 12) | `outputs/rerun-test_ba13_*.out` |
| Historical thresholds | All 1,386 historical fields equal the certificate. All thresholds are unchanged. | §2 |
| Required-accuracy algebra | Every row's cap and joint right-hand side reproduces from `model.feasible_cap`. Every finite ledger equals the sum of its finite terms. | `outputs/quantitative-review.json` A |
| Fit operator norm | The package's f = 5.171–5.185 matches my unweighted ‖I−H‖∞ = 5.182 over all 21 rows. Over fitted rows only it would be 1.909, so the package is about 2.7× conservative. | A |
| Tangents from raw outputs | All nine N512/1024/2048 tangents recomputed from the 18 raw JSON files at 60 digits match the package to < 10⁻⁴⁰. Orders lie in 3.948–4.003. | C |
| Holdouts | Curvature prediction is confirmed at both withheld targets for all 15 groups. The residual equals the recorded X2 B-alignment offset. | §3.7 |
| Chemical-inverse replay, thermal common factor, μ_n identity | Consistent with prior review and preflight | §3.8 |
| Disposition E under predeclared rules; refusal behavior; mutation controls | Correctly implemented | §2 |

### 1.2 Defects, gaps and mis-specifications (ranked by consequence)

1. **Proposed owner decision C does not name the surface boundary law.**
   - `METHODS.md` defines Q "with the declared finite-pressure boundary law", but the ratification text C never declares it.
   - Without it, `moving_cut_input_remainder` cannot close: there is no object to compare to.
   - The tangent object (fixed-cut H-derivative versus derivative along the cut law) is also unstated.
   - **Amendment:** write the historical cut law into C (§4).
2. **The six null terms carry no proposed quantitative statements.**
   - `OWNER_DECISION.md` lists C-REF, C-TAN, C-CONT and C-INPUT as "pending" with no quantitative content.
   - The owner is therefore asked to accept conditions that have not been written.
   - Appendix A supplies proposed statements for owner consideration.
3. **Two refusals rest on a standard stricter than the proposed standard B.**
   - `reference_stage_remainder` is refused as "not quantity-specific uniform stage enclosures".
   - `reference_arithmetic_remainder` is refused as "do not bound every passive functional".
   - Standard B asks for validated estimates, not uniform enclosures. Quantity-specific evidence exists: clean fourth-order Richardson sequences for each referenced quantity at three targets, plus precision ablations (§3.2).
4. **`reference_parameter_remainders` partly duplicates charged terms.**
   - EOS and centre sensitivities are already charged as validated estimates.
   - A separate "remainder beyond the sensitivity" is a remainder of a remainder. It needs only the convergence-model statement (exponential Gauss convergence; known centre order), not new data (§3.4).
5. **The caps ignore drift structure.** This is conservative, not wrong.
   - The required pointwise caps charge a paired error as 2·sup|e|. That is valid for arbitrary errors.
   - For the smooth component of a common error, the drift contamination is ≤ sup|e′|·10⁻⁶. That makes the real requirement on smooth shared errors O(1) relative to the physical B-sensitivity (§3.1).
   - The caps are therefore far more demanding than necessary for the error type that is genuinely unobservable.
6. **The original missing-error control was detected only by a TypeError crash.**
   - The original mutation (`5957154`) removed the `None` guard, so the test errored inside `math.isfinite(None)` rather than failing an assertion.
   - Controls-v2 (`328ccf4`) substitutes zero and is caught semantically.
   - The package discloses this as "strengthened". It should say the original control was incidental. This is a transparency point with no scientific consequence.
7. **The historical tool `tools/linux/ba13_closure/assess.py` was edited in `2cac8b0`** (2026-10-26 package, before this entry).
   - The edit means a nonzero unclassified value now becomes `unresolved`, never `analytic_bound`.
   - It is conservative, and the closure evidence files are unchanged.
   - It is still an edit to a historical assessor and should be named in the sidecar provenance.

### 1.3 Statements that understate or overstate the evidence

| Statement in package | Assessment |
|---|---|
| Holdouts "remain successful … largest criterion ratios 0.117275 and 0.116597" | **Understates.** The ratios cluster at about 0.1 because K = 10 × measured curvature. They are quantitative curvature confirmations, not marginal passes (§3.7). |
| "Four common missing models account for the approximate rows" | Correct as bookkeeping. Scientifically they are missing *model statements*. The evidence for each exists. |
| C-REF: "no quantitative bound follows from that comparison alone" | True of the comparison alone. Implementation independence (verified) plus V&V solution verification supports a class-B/C estimate (§3.2). |
| "Option A remains numerically plausible" | **Understates.** Option A is viable and evidentially supported. Nothing measured threatens any threshold. |
| "A/C clarify the old scientific contract; B is evidence policy, D record format, E validation method" | Agreed (§4). |
| "Investigation category shares do not become new historical physical thresholds" | Agreed, and important (§4). |
| "Full metric_structure also lacks nodewise reference data" | Correct. This is the only true evidence gap (§6). |

---

## 2. Authentication

Script: `scripts/authenticate.py` → `outputs/authentication.json` (120 checks, 0 failed). It is independent of the package's `verify_delivery.py` and `authenticate_entry.py`.

**Git identity and history**
- HEAD = `bb0f9c5`. The working tree is clean, including untracked and ignored files. The named branch points at HEAD.
- Linear chain entry `84e1062` → `d5b40b3` → `5957154` → `328ccf4` → `bb0f9c5`: exactly four commits, no merges.
- `093d121` (previously reviewed) and the certificate-producing recovery commit `8e6ac58` are ancestors.
- Every local branch tip (14) is an ancestor of the final commit. Master is unchanged at `9a5c1ec`.

**Change boundary**
- Entry..final adds 84 files, all under the new package, tool and test paths. It modifies exactly one file: ADR-0021, append-only since both the entry and `093d121`, and still `PROPOSED`.
- The predeclaration is unchanged after `d5b40b3`. Its SHA-256 `4beffcf3…` is bound in all three receipts.
- These protected paths are unchanged since the entry: production `CompactStar/`, `GOVERNANCE.md`, `AGENTS.md`, `SCIENTIFIC_INVARIANTS.md`, the preflight, both dependency locks, `ExternalDependencies.cmake`, the fitter, X2 source, and every earlier validation tree, tool and test.
- Since `093d121`, the only protected-path change is the conservative `ba13_closure/assess.py` hardening (finding 7).
- Production code and dependency locks are unchanged since the certificate-producing commit `8e6ac58`. Master differs only by the pre-existing Linux-transition lineage (25 files), which predates this work.

**Certificate and fitter**
- Certificate SHA-256 `4212cad0…41fb`.
- Fitter blob SHA-256 `9b3ef1ce…2370`, identical at `8e6ac58`, `093d121`, entry and final.

**Evidence and manifests**
- `evidence/SHA256SUMS`: 84 of 84 entries match both the working tree and the committed blobs.
- External manifest: 79 of 79 entries match on hash and size, with no unlisted files. The manifest cannot list itself; its copy matches the committed one.
- The committed `matrix.json.gz`, `matrix.tsv.gz` and `conditional-intervals.json.gz` decompress to the external raw files. Ten other compact files are byte-identical to their external copies.
- `source-manifest.json`: all inputs authentic.
- Inherited manifests named in `authentication.json`: 1,993 entries re-verified. These are five external manifests, three `inputs.json` groups, the prior review, the 2026-10-26 evidence and the X1 SHA256SUMS. The total matches the package figure.

**Receipts**
- Three PASS receipts: Slurm jobs 19138–19140, zero stellar runs, CPU 2.657819 + 5.131496 + 0.574553 = **8.363868 s** (ceiling 600 s).
- Every retained stdout/stderr matches its recorded hash.
- The benchmark and validate receipts bind `5957154` sources. `assess.py` and `model.py` there equal the final blobs, so the assessment output is valid at HEAD.
- The test file and `run.py` changed in `328ccf4`. The controls-v2 receipt binds the final versions.
- After `328ccf4` only packaging tools (`publish.py`, `verify_delivery.py`), docs and evidence were added.

**Mutation controls**
- All 14 mutation runs (7 original, 7 controls-v2) exit nonzero.
- Controls-v2 detects every mutation by a semantic assertion. The original missing-error control was a crash (finding 6).

**Handoff**
- SHA256SUMS: 3 of 3 entries match.
- HANDOFF.json: head, base and branch correct; `owner_ratified=false`; no GitHub publication; disposition E.
- The incremental bundle verifies against this repository, and its head is the final commit.
- All 80 archive members equal the manifest entries plus the manifest itself.

**Matrix**
- 462 unique rows (21 × 22).
- All 1,386 historical `drift_bound`, `threshold` and `utilization` values equal `certificate.json` `samples[i]` exactly.
- Status counts are 21 PASS (exactly the P2 rows), 0 FAIL and 441 UNQUALIFIED.
- No non-P2 row carries a proposed total.
- The 2026-10-26 status is preserved in every row.

**Registry, owner record and prior review**
- `acceptance-registry.json` is empty and `PROPOSED_NOT_ACCEPTED`. The owner record reads "PROPOSED — NOT ACCEPTED" with owner and date UNSET.
- The prior independent review's 14 SHA256SUMS entries still verify.

**Preserved failures.** These remain byte-unchanged: the original X1 failures, X2 reference INCONCLUSIVE, the binary64/Newton N512 failure (retained in `ablations.json` as `"status": "FAIL"`), BA12/BA12R, and every earlier BA13 outcome.

**Scope limit.** This is content and ancestry authentication. It does not establish signed authorship.

---

## 3. Remaining scientific uncertainty — quantitative review

Script: `scripts/quantitative_review.py` → `outputs/quantitative-review.json`, with a printed summary in `outputs/quantitative-review.txt`. The outputs keep four data types apart: **measured** differences, **estimates** under a stated model, the package's **requirements** (caps), and **multipliers** (requirement ÷ evidence scale). A multiplier measures how wrong an estimate would have to be to matter. **It is not an error bound and was not used as one.**

### 3.1 What the requirement actually constrains

At achieved production abscissa x, the package writes the error as

    P − Q = (P − I_X) + (I_X − I_Q) + (I_Q − Q).

- (P − I_X) is measured at all 21 samples.
- (I_X − I_Q) is the reference error at targets 0, 10 and 20, interpolated.
- (I_Q − Q) is the inter-sample physical remainder.

The frozen-background criterion is a **drift** criterion. What enters the left side is the variation of the error across the window, e(x) − e(0), not e itself. Split e into two parts:

- **Non-smooth part.** This covers adaptive-mesh jitter, support-edge pinning and envelope-node snapping. It lives in production and is **measured** at every sample in (P − I_X). Under decision C it is never needed between samples.
- **Smooth part.** This is any error of either code that varies smoothly with B. It contributes at most sup|e′|·|x| ≤ sup|e′|·10⁻⁶. To exhaust the joint right-hand side, its relative B-derivative must reach `rhs/D` times the physical drift rate. The values are listed below.

| Group | `rhs/D`: smooth shared-error derivative needed, in multiples of the physical B-sensitivity |
|---|---:|
| N_μ | 0.849 |
| t_μ | 1.65 |
| Z_row_npμ | 1.85 |
| Ltilde_Mμ | 2.89 |
| t_e | 14.9 |
| N_e | 19.5 |
| every other X2-referenced group | 26–530 |

Example: for N_μ to fail, a shared error would have to change d ln N_μ/dx = 54.07 by ≥ 46 per unit x, i.e. by 85 %. Production and X2 compute t_μ = dN_μ/dB independently and agree to 2.2×10⁻¹⁰ absolute, or 6.6×10⁻⁸ relative.

The package's caps instead bound an *arbitrary* uniform error e by 2·sup|e| + f·sup|e|. That is valid but far more demanding than needed for the smooth component. The two treatments agree that nothing measured or estimated approaches the requirement.

**Implication for the evidence standard.** The only errors that are both unobservable and consequential are those shared by independent implementations that alter the B-*derivative* by O(1). Such an error is a specification or physical-model error. It is outside BA13's numerical scope and outside what any numerical proof could address. Numerical V&V suffices for everything else.

### 3.2 Reference integration and common numerical bias (C-REF)

| Component | Evidence (retained) | Scale | Class |
|---|---|---|---|
| GL4 discretization (all X2-referenced scalars) | Ladders N = 32…4096 at 3 targets; p = 4.000 ± 0.0003; signed, consistent increments. Manufactured order 4.00025. | Relative tails ≤ 2.5×10⁻¹⁵. N_μ charged 1.9×10⁻¹⁵. | **B** |
| Implicit stages | Newton residual monitored. Quantity-specific orders stay clean down to tails ~10⁻¹⁵ (N_μ order 3.99997 at target 20). Non-asymptotic stage or roundoff noise larger than the finest increment would break the order. | ≤ finest Richardson increment (~10⁻¹⁴ relative) | **B** (needs the stated model) |
| Floating point and summation | Plain vs compensated sums ≤ 5.9×10⁻¹⁸. Binary64 vs long double ≤ 1.4×10⁻¹⁵ on primary outputs, implying ~7×10⁻¹⁹ for long double if error scales with ε. Finest tangent increments agree across 3 targets to half-spread ≤ 1.3×10⁻¹⁶ (t_μ). | ≤ 10⁻¹⁵ relative | **B** for primary outputs and tangents. Passive functionals are covered by the order-cleanliness argument above. |
| Unseen common bias | X2 shares only `Zaki/Physics/Constants.hpp` (the sole non-standard include in `tests/bnv/x2_fixed_step.cpp`, verified). X1–X2 ≤ 9.1×10⁻¹² on primary outputs. Production–X2 ≤ 1.5×10⁻¹⁰ (N_μ) and ≤ 2.2×10⁻¹⁰ absolute (t_μ). | Inter-implementation agreement 10⁻¹²–10⁻¹⁰. Exception: Ltilde_Me, where production's measured pinning error is 3.0×10⁻⁶ but the X1/X2 agreement on L is 9.1×10⁻¹². | **C, irreducibly** |

The unseen common bias is the one genuinely irreducible conditional. Its quantitative content is: *no shared error source between X2 and production other than the physical constants and the model specification; any residual shared numerical error is ≤ 10× the observed inter-implementation agreement.*

Robustness:
- In the arbitrary-error treatment, the requirement exceeds that allowance by a factor of at least about 10³. For example, the N_μ cap 6.39×10⁻⁶ is 4×10³ × (10 × 1.5×10⁻¹⁰).
- In the smooth-error treatment, the robustness is O(1) of the physical sensitivity (§3.1).

**Adjudication with the 2026-10-20 review.** That review said the ablations resolve stage and arithmetic. The package says they are primary-output only. The package is right about the ablations' scope. It is wrong to conclude that no quantity-specific evidence exists: each referenced quantity's own Richardson sequence is that evidence. The gap is one sentence of model statement.

### 3.3 Tangent differentiation and correlated errors (C-TAN)

Recomputed from the 18 raw outputs at 60 digits (`outputs/quantitative-review.json` C):

| | t_n | t_e | t_μ |
|---|---:|---:|---:|
| GL4 orders (targets 0/10/20) | 3.953 / 3.948 / 3.959 | 3.987 / 4.003 / 3.997 | 3.993 / 4.002 / 3.997 |
| Finest tail \|Δ₂\|/15 | 1.57×10⁻¹⁴ | 1.26×10⁻¹⁴ | 2.68×10⁻¹⁵ |
| Cross-target half-spread of finest increment (empirical noise ceiling) | 9.3×10⁻¹⁶ | 1.0×10⁻¹⁵ | 1.3×10⁻¹⁶ |
| Stencil-width tail (3 widths, ratio 4.0000 ± 10⁻⁵) | 5.7×10⁻¹² | 1.4×10⁻¹¹ | 8.0×10⁻¹² |
| Moving-cut vs fixed-cut (cross-geometry residual, max over 3 pairs) | 2.7×10⁻¹³ | 1.5×10⁻¹³ | 3.4×10⁻¹⁵ |
| EOS / centre sensitivities propagated into t (t_μ row) | — | — | 4.6×10⁻¹⁷ / 9.5×10⁻¹⁷ |
| Requirement (single-axis cap), absolute dN/dB | 1.34×10⁻⁵ | 4.03×10⁻⁷ | 2.93×10⁻⁸ |

**The package's concern.** It requires a joint residual J ≥ |a − t·d|. It worries that correlation of one truncation term does not cancel the derivatives of EOS, stage or common-mode errors.

**Why the concern is already answered.** Every component above was estimated **on the ratio t itself**, not on numerator and denominator separately:
- the Richardson sequence is a sequence of secant ratios;
- the stencil widths are ratios;
- the EOS and centre ablations were propagated into t.

J can therefore be stated as the absolute sum of tangent-level estimates, with no cancellation assumed. For t_μ that is ≈ 1.6×10⁻¹¹ + 5.4×10⁻¹⁵ + 3.4×10⁻¹⁵ + 1.4×10⁻¹⁶ ≈ 1.6×10⁻¹¹. That is 1.8×10³ below the cap.

The denominator D = B₊ − B₋ = 1.26×10⁵² carries a relative error ~10⁻¹⁵. Separation uD < |D| is trivial.

**Class.** B, given the asymptotic h⁴ (grid) and h² (width) models. Both are tested: by order at three targets and by width ratio at three widths. The moving-cut entry is *definitional* once decision C names which tangent is certified (§4). Until then it is a measured difference of ≤ 3×10⁻¹³.

### 3.4 EOS and centre remainders

| Remainder | Evidence | Class |
|---|---|---|
| EOS quadrature | 96 → 192 points changes primary outputs by ≤ 1.0×10⁻¹⁸. Gauss quadrature of the analytic T = 0 momentum integrands converges exponentially, so the 96-point error ≈ \|Q₉₆ − Q₁₉₂\|. Thermodynamic, charge and β residuals < 4.7×10⁻¹⁷. | **B** (A-type theorem, conditional on integrand analyticity) |
| Equilibrium inversion | Bisection to long-double resolution, with residuals as above | **B** |
| Centre series | Analytic omitted order O(x₀^(5/2)), with x₀ ∝ N⁻². x₀ factors of 10⁻² and 10⁻⁶ change outputs by ≤ 7.9×10⁻¹⁸. Under the known order, the nominal remainder is at or below rounding. The centre error falls faster (N⁻⁵) than GL4 (N⁻⁴), so it is also inside the Richardson tail. | **B** (analytic order + measured coefficient) |
| "reference_parameter_remainders" | These are the above, restated as remainders beyond the charged sensitivities | Duplicate for primary outputs and tangents. It closes with the convergence-model statement. |

### 3.5 Common-B alignment

- The computed-B alignment residual is ≤ 5.1×10⁻¹⁷ relative. The B-functional GL4 tail is ≈ 1.5×10⁴³ counts, about 2×10⁻¹⁴ relative. Together they give an H error of ≈ 1.5×10⁻¹⁴.
- The local secant slope B_H ≈ 1.918×10⁵⁷ is positive at all three targets. It varies by 5.8×10⁻⁶ relative across the window, and by ≤ 3.3×10⁻⁹ between resolutions at one target.
- The term (residual + u_B)/m is an established bound (A), given a uniform m > 0. The hypothesis m = 0.999 × the sampled secant is a conditional (C). Its 0.1 % slack is 170× the observed variation across the window, and m is far from zero.
- It is already charged: N_μ 4.2×10⁻¹², t_μ 9.8×10⁻¹⁵.

### 3.6 Nonlinear fit-input propagation (C-INPUT)

- **Fixed design and weights.** The response is exactly (I − H)e (A). The package charges it with f = 5.18. That is the all-21-row norm; only the 11 fitted rows enter the §11 residual test, where f = 1.91.
- **Abscissa perturbation.** The per-gap coordinate estimate is the sum of the two endpoint design errors, from production/X2 B differences in `ba13-closure` m3. It is ≤ 5.1×10⁻¹² in x. Its first-order effect is slope × δx. For N_μ that is ≤ 54 × 5.1×10⁻¹² ≈ 2.7×10⁻¹⁰. Relative to the design, δx/Δx ≤ 1×10⁻⁴, so the second-order hat-matrix remainder is ~10⁻⁸ × |residual| — negligible. This is A (perturbation theory with an explicit second-order term), given the δx range (B).
- **Weights.** Under interpretation (2) the historical weights belong to the frozen fitter, not to the physical error model. The only value-dependent weights, the 64ε floor, change by ~10⁻¹⁰ relative under the error sizes in question — a second-order effect.
- **B₀ normalization.** Analytic propagation (A) through the exact coordinate identity, with δB₀ from the B residual (B). Already charged: N_μ 1.8×10⁻¹¹.
- **What closes `fit_nonlinear_input_remainder`.** One derivation paragraph plus the two input ranges. No calculation.

### 3.7 Continuous inter-sample curvature (C-CONT)

The package correctly states two things:
- The interpolation bound K h²/8 is an established theorem given |Q_xx| ≤ K.
- No finite sample set can prove the supremum. The C² bump counterexample makes the point.

The question is how well-supported the hypothesis is. The evidence is unusually strong.

**1. Holdouts are a quantitative prospective confirmation, not a marginal pass.** The criterion ratios cluster at 0.10–0.117 because K = 10 × measured curvature in 13 of 15 groups. Dividing by the measured curvature gives the prediction accuracy of the curvature model at an unseen point:

| Group | Observed ÷ predicted (target 5 / target 15) | K from |
|---|---|---|
| N_μ | 1.014 / 1.014 | 10 × target curvature |
| Z_row_npμ | 0.990 / 0.990 | wider stencil |
| Ltilde_Mμ | 1.064 / 1.062 | 10 × target |
| N_e | 1.084 / 1.081 | 10 × target |
| surface_gravity, P0, P1, μ_B, envelope, Tsurface_inf, Cstar, radius | 1.08–1.22 | 10 × target |
| Ltilde_Me | 1.613 / 1.595 | 10 × target |
| N_n | 0.701 / 0.642 | 10 × stencil |
| Z_row_npe | 0.663 / 0.673 | wider stencil |

For **every** group and both holdouts, |observed − predicted| is ≤ 17 % of the separately frozen non-interpolation part of the criterion (the X2 bands plus root alignment). The Ltilde_Me excess of 3×10⁻¹⁴ sits inside its 1.8×10⁻¹³ root band.

A pure quadratic with the curvature measured at targets 0, 10 and 20 therefore explains both withheld stars to within their numerical noise. The local stencil curvature at three widths agrees with the window-scale target curvature to 5×10⁻⁴ for N_μ.

The residual is not noise. Dividing each group's post-curvature excess by its own slope gives one common abscissa offset:

- |δx| = 1.60–1.66×10⁻¹⁴ in 12 groups; 1.49–1.69×10⁻¹⁴ including P0, P1 and μ_B. The sign is consistent with each slope.
- The slopes span 0.19–54, a factor of 290.
- The offset equals the recorded X2 B-alignment residual at the reference targets: computed-B residual/B₀ = 1.56–1.87×10⁻¹⁴; B-functional GL4 tail/B₀ = 1.99×10⁻¹⁴.
- N_n is the only exception (−3×10⁻¹⁵). Its curvature prediction comes from the stencil, and its excess is at the noise of that prediction.

The holdouts are therefore explained, to about three digits, by two recorded quantities: the measured curvature and the already-charged X2 alignment term. Nothing is left over that a hidden inter-sample feature or an unmodeled reference error could account for (`outputs/quantitative-review.json` B2).

**2. The analytic structure is known.**
- Outputs are analytic in H_c away from species thresholds (analytic EOS between onsets, plus analytic dependence of ODE solutions on parameters).
- The nearest threshold is the muon onset, at H_c − H_μ = 0.01972. The H_c window is 3.97×10⁻⁷, so the onset is **5×10⁴ window widths** away in H, the directly established statement. A local-linear extrapolation of B(H_c), whose slope is measured only inside the window, puts it at x ≈ 0.0497.
- A single power-law onset model N_μ ∝ (x − x_on)^p with the fitted slope 54.075 gives p = 2.685. It predicts a relative curvature of 1,835, against **2,039 measured** (ratio 0.90).
- A Cauchy estimate in the same model gives 9.6×10³, below the K = 2.04×10⁴ used.
- This is an illustration under a model (C), not a bound. It shows that the factor-10 K is consistent with the physical analytic radius.

**3. Failure would need a gross anomaly.** For the N_μ joint right-hand side to be consumed by the inter-sample term alone, the curvature would have to be 1.5×10¹¹. That is **7×10⁷ × measured**, and would require a non-analytic feature within ~10⁻⁵ in x of the window. The exactly known thresholds exclude one.

**Class.** C: a named, prospectively validated hypothesis with explicit quantitative content (K per group as in `conditional-intervals.json`) and a failure multiplier of ≥ 10⁷ for the muon number. It cannot become A without a computer-assisted analyticity bound, which is disproportionate here. Owner acceptance is the right instrument.

**Endpoint sliver.**
- The achieved x₂₀ = −9.999999990056×10⁻⁷ falls short of −10⁻⁶ by 9.94×10⁻¹⁶.
- That is 2×10⁻⁵ of the preflight's own targeting tolerance τ_B = 5×10⁻¹¹ B₀, which the accepted contract already admits.
- Its effect is ≤ 8.4×10⁻⁹ of any cap.
- It is closed by the same C-CONT derivative statement and needs no root campaign. The package agrees.

### 3.8 Moving stellar surface, support boundaries and remaining groups

- **Cut law.**
  - ln P_s is linear in x: slope 8.3839 per unit x, maximum linear-fit residual 2.07×10⁻¹², total change over the window 8.38×10⁻⁶.
  - The surface moves 1.64×10⁻⁵ km over the window, with a 6× margin to a terminal-cell crossing.
  - Both holdouts used their own retained P_s, so the across-target curvature already includes the real cut motion.
  - If decision C writes this law into Q's definition, the moving-cut "remainder" is only the deviation of the historical cuts from the declared law (≤ 2×10⁻¹² in ln P_s) times the surface sensitivity. That is effectively A.
- **Species support (no crossing).**
  - H_c − H_μ = 0.0197 and H_c − H_n = 0.1307, against an H_c window of 3.97×10⁻⁷.
  - Given monotone B(H_c) on the stable branch (B_H > 0, as sampled), no threshold can be crossed. This is A conditional on monotonicity (C; sampled B_H varies by 5.8×10⁻⁶ across the window).
  - Neutron-edge and envelope-node selections are constant at all 21 samples. They belong to production representation, not to Q.
- **μ_n_profile.** The equilibrium identity μ_n^∞ = μ_B^∞ at η = 0 is an exact contract identity within the equilibrium/reference definition of §11 (A). It does not apply to the separate η-dependent actual potential, which §8 evaluates at run time.
- **Z rows.** The direct X2 comparison is the primary route. The 50-digit independent-endpoint majorant (npe 1.96×10⁻⁷, npμ 5.60×10⁻⁷) is an alternative, conditional on its inherited background set, and is correctly not double-counted.
- **Cstar and thermal knots.** The cache factor F(T) cancels exactly in same-knot relative drift (A), so the frozen-relative row is unaffected. The absolute bias, up to 6.555×10⁻⁴, is an absolute-accuracy limitation only (§5).

### 3.9 Classification of every remainder

Key:
- **A** — established mathematical bound
- **B** — independently validated numerical estimate
- **C** — conditional estimate with explicit assumptions
- **D** — missing

| Remainder | Class achievable from existing evidence | Evidence scale | Requirement (tightest) | Multiplier | What converts it to a row value |
|---|---|---|---|---:|---|
| X2 GL4 discretization | B | ≤ 2.5×10⁻¹⁵ relative | 6.4×10⁻⁶ (N_μ) | > 10⁹ | already charged |
| X2 implicit stages | B | ≤ finest increment, ~10⁻¹⁴ | 6.4×10⁻⁶ | ~6×10⁸ | model statement |
| X2 arithmetic and summation | B | ≤ 10⁻¹⁵ | 6.4×10⁻⁶ | > 10⁹ | model statement |
| EOS evaluation and inversion | B (A-type for quadrature) | ≤ 10⁻¹⁸ | 6.4×10⁻⁶ | > 10¹² | already charged; convergence statement |
| Centre series | B (analytic order) | ≤ 8×10⁻¹⁸ | 6.4×10⁻⁶ | > 10¹¹ | already charged; order statement |
| Common-B alignment | A given m > 0 (m: C) | 4.2×10⁻¹² (N_μ) | 6.4×10⁻⁶ | 1.5×10⁶ | already charged |
| Unseen common reference bias | **C (irreducible)** | inter-implementation agreement ≤ 1.5×10⁻¹⁰ | 6.4×10⁻⁶ | 4×10⁴ (or O(1) of the physical sensitivity, §3.1) | named condition C-REF |
| Tangent GL4 on ratio | B | 2.7×10⁻¹⁵ (t_μ) | 2.9×10⁻⁸ | 10⁷ | already charged |
| Tangent stencil width | B | 8.0×10⁻¹² (t_μ) | 2.9×10⁻⁸ | 3.7×10³ | already charged (conditional) |
| Tangent joint J | B | ≈ 1.6×10⁻¹¹ (t_μ) | 2.9×10⁻⁸ | 1.8×10³ | model statement |
| Moving vs fixed cut (tangent) | definitional, or B measured | ≤ 3.4×10⁻¹⁵ (t_μ) | 2.9×10⁻⁸ | 10⁷ | decision C amendment |
| Moving-cut input | A once the law is declared | 2×10⁻¹² in ln P_s | — | — | decision C amendment |
| Continuous curvature | **C (validated)** | K h²/8: N_μ 6.4×10⁻¹² | rhs 4.6×10⁻⁵ | 7×10⁶ (curvature × 7×10⁷) | named condition C-CONT |
| Species support | A given B_H > 0 (C) | ΔH margin 5×10⁴ windows | categorical | — | named condition C-SUPPORT |
| Endpoint sliver | A given C-CONT | ≤ 8×10⁻⁹ of cap | — | 10⁸ | covered by C-CONT |
| Fit fixed projection | A | charged (N_μ 3.8×10⁻¹¹) | — | — | already charged |
| Fit design and weight nonlinearity | A given input ranges (B) | ≲ 10⁻¹³ second-order | 6.4×10⁻⁶ | > 10⁷ | derivation |
| B₀ normalization | A given δB₀ (B) | 1.8×10⁻¹¹ (N_μ) | 6.4×10⁻⁶ | 3.5×10⁵ | already charged |
| μ_n_profile identity | A (contract identity) | paired ≤ 5.5×10⁻¹¹ | — | — | owner acceptance |
| Cstar common factor F(T) | A (exact cancellation) | 0 in relative drift | — | — | none |
| metric_structure nodewise reference | **D** | only P − R = 1.27×10⁻⁶ (not independent) | ≈ 1.3×10⁻⁵ | unknown | observer, or the consumer route (§6) |

Only `metric_structure` is D. Every other term is A, B or C on existing evidence. The C items are the four named conditions plus the inherited m > 0 slope hypothesis.

---

## 4. Option A, interpretation (2) and the five owner decisions

None of these is ratified here. For each, the table gives consistency with the accepted preflight, its nature, the claim acceptance authorizes, what still needs itemized acceptance, whether the evidence is adequate, and my recommended amendment.

| | Consistent with preflight §11 | Clarifies or changes the contract | Authorizes | Still needs itemized acceptance | Evidence adequate | Recommended amendment |
|---|---|---|---|---|---|---|
| **Option A** (historically faithful static candidate; independent error completion; fixed grids, targets, knots, thresholds, fit window) | Yes. No production, fitter or threshold change; verified value-for-value. | Selects among ADR-0021 alternatives. Not a contract change. | Qualification may be sought on the unchanged certificate. | — | Yes | None |
| **Interpretation (2)** (frozen error-aware right side; full numerical contribution on the left) | Yes. It is the most literal reading of "the full drift plus numerical and nonlinearity terms must still remain below T_X". | Clarification. It removes the historical operational reading in which the 5u allowance is the entire numerical accommodation. | Left side = historical envelope + all numerical and inter-sample terms. Right side = historical T. | — | Yes. For most rows u = 0 or 64ε, so T = 10⁻⁴ or 10⁻⁴\|t\| with no hidden allowance. | None |
| **A. Full interval** | Yes. The runtime monitor uses continuous \|ΔB\|/B₀ ≤ 10⁻⁶, so continuous coverage is implied. | Clarification (strengthening). It makes inter-sample coverage explicit. Thresholds, norms and support requirements are unchanged (verified). | The certificate statement holds at every x in [−10⁻⁶, 0], not only at 21 points. | C-CONT (curvature, B_H > 0, topology); endpoint sliver | Yes, under C-CONT (§3.7) | State that the endpoint sliver (9.94×10⁻¹⁶ = 2×10⁻⁵ τ_B) is closed by the C-CONT derivative statement. The preflight's τ_B already admits larger target mismatch. |
| **B. Evidence standard** | Yes. §11 already propagates "governed numerical error", which historically were estimates, not enclosures. | **Substantive policy selection**: §11 defines no evidence classes. It formalizes existing practice. | Admissibility of analytic bounds, validated estimates and itemized conditionals. **By itself it fills no row.** | None by itself | Not applicable | Make "independently validated numerical estimate" operational (Appendix A, item S). Require every conditional to state its failure multiplier. Fix the qualification wording "qualified under named validated numerical estimates and itemized conditionals; not a rigorous enclosure". |
| **C. Certified object** | Yes. §11 bounds the drift of governed equilibrium quantities and evaluates the η-dependent actual potential separately. The μ_n restriction matches the §11 table row exactly. | Clarification | Q is the governed model's equilibrium response; production's per-sample numerical error is included; production's discrete mesh, pinning and envelope-node rules are not certified between samples. | None, once amended | Yes | **Name the finite-pressure surface law** (ln P_s linear in x through the 21 retained cuts; slope 8.3839, residual ≤ 2.1×10⁻¹²). **Name the certified tangent**: the derivative along that law, or the fixed-cut H-derivative. They differ by ≤ 3×10⁻¹³ absolute. Without the law, `moving_cut_input_remainder` cannot close. |
| **D. Record** | Yes. The certificate is unchanged, and FAIL/INCONCLUSIVE/UNQUALIFIED cannot be overwritten. | New record format. Not a physics change. | An authenticated sidecar as the §11 discharge record. Operational attachment needs separate authorization. | — | Yes | State that the runtime monitor keeps consuming the historical U_X, while the sidecar reports qualified utilization (t_μ 0.376276 → 0.377301; N_μ 0.540753 → 0.540802). Name the `2cac8b0` assessor edit in provenance. Commit the registry separately from the assessor (the model already refuses self-authority). |
| **E. Target transfer** | Yes. It is a validation-layer method; §11's per-sample governed recomputation is unchanged. | Validation-method selection | Direct comparison of all 21 samples with X2 interpolated at achieved B, with matched P_s. | C-CONT; C-REF for the X2 band at targets | Yes: two prospective holdouts confirm the curvature model to within bands (§3.7) | None, beyond keeping "two holdouts do not prove a supremum" (correct as written) |

**Investigation criteria are not requirements.** The package is right that the 12 category shares, the factor-ten K, the 10 %-of-share stop and any observer order gate are investigation criteria, not §11 acceptance thresholds. Exceeding a share is not disposition C. An unconditional independent tangent solver would be an added requirement.

---

## 5. Certificate scope

| Possible scope | Supported? | Reason |
|---|---|---|
| Full BA13 frozen-background qualification (all 22 groups, contract-literal) | **Not yet** | `metric_structure` lacks a nodewise reference (D). The six null terms need written models and owner acceptance. |
| **Narrower, explicitly conditional qualification** | **Yes, after owner adjudication and mechanical sidecar completion** | Every other group's remainders are A/B/C on existing evidence, with multipliers ≥ 10³ (§3.9). The holdouts validate the one substantive continuity hypothesis. |
| Only continued mathematical benchmarking | Not required | Further benchmarking has no decision value. Every estimated term is orders below its requirement. |
| No defensible qualification | Not the case | Nothing measured threatens a threshold. No production correction is indicated. |

**The narrower claim, stated exactly:**

> For the authenticated Structure-1 star on the stable branch, at every depletion |ΔB|/B₀ ≤ 10⁻⁶, each frozen quantity in §11 other than the nodewise metric profile has an equilibrium drift that stays below its unchanged historical threshold. This includes production's numerical error at the 21 samples, the reference error and the inter-sample remainder. The statement holds under the named validated numerical estimates and the itemized conditions C-REF, C-TAN, C-CONT, C-INPUT, C-SUPPORT, C-CUT and C-MUN, and under the owner's decision on how the metric is covered.

**Frozen-relative versus absolute accuracy.**
- BA13 certifies only the **frozen-relative** statement: that freezing at B₀ changes each quantity by less than its threshold across the window.
- It does **not** certify absolute accuracy of the B₀ background. Examples outside its scope:
  - Cstar carries an exact cache-law factor up to 6.555×10⁻⁴ in absolute terms. It cancels in same-knot relative drift.
  - Production differs from X2 by up to 3.0×10⁻⁶ in Ltilde_Me through pinning. This is charged in the relative row.
  - EOS-model and physical-input uncertainty are entirely out of scope.
- Absolute thermal-evolution accuracy claims need their own budget. The matched-control design (BA10b) cancels the common absolute background error only to first order in target-minus-control comparisons, and that is not a BA13 result.

**Muon sector — measured headroom, not bounds:**

| Group | Historical U | U including finite ledger | Remaining joint rhs | Finite ledger (largest term) | Single-axis cap | `rhs/D` |
|---|---:|---:|---:|---|---:|---:|
| N_μ | 0.540753 | 0.540802 | 4.59×10⁻⁵ | 4.83×10⁻⁹ (fit reference inputs 3.3×10⁻⁹) | 6.39×10⁻⁶ | 0.849 |
| t_μ | 0.376276 | 0.377301 | 2.10×10⁻⁷ abs. | 3.46×10⁻¹⁰ (P − I_X 2.23×10⁻¹⁰) | 2.93×10⁻⁸ abs. | 1.65 |
| Z_row_npμ | 0.350922 | 0.350964 | 6.49×10⁻⁵ | 4.22×10⁻⁹ (fit reference inputs 3.0×10⁻⁹) | 9.04×10⁻⁶ | 1.85 |
| Ltilde_Mμ | 0.256811 | 0.256824 | 7.43×10⁻⁵ | 1.26×10⁻⁹ (P − I_X 7.2×10⁻¹⁰) | 1.03×10⁻⁵ | 2.89 |

- The muon sector is **physically** stressed: d ln N_μ/dx = 54 consumes 54 % of the N_μ threshold at the ceiling. It is **not numerically fragile**: finite numerical terms use ≤ 0.16 % of any muon headroom.
- The linear envelope would reach U = 1 only at |ΔB|/B₀ ≈ 1.85×10⁻⁶, so the hard ceiling binds first.
- These figures are measurements and estimates under their stated models. They are not converted into certified bounds.
- Experiment interpretation: any muon-sensitive conclusion at depletion near 10⁻⁶ carries a frozen-background error of up to ≈ 5.4×10⁻⁵ relative in N_μ. That is within the contract, but not negligible against 10⁻⁴.

---

## 6. The metric-structure profile and the passive nodewise observer

**What exists.**
- P − R = 1.27×10⁻⁶. Here R is a reconstruction from production's own profiles, so this is representation error only, not P − Q.
- Historical utilization 0.033 (drift 3.30×10⁻⁶ scaled).
- Joint requirement by the package's algebra ≈ (10⁻⁴ − 3.30×10⁻⁶ − 1.27×10⁻⁶)/7.18 ≈ 1.3×10⁻⁵ scaled per node.

**Is it needed for the scientific claim?** Probably not, subject to one check.

- In the controlled experiment (§§8–9), the frozen metric enters through:
  1. the immutable source event measure, with infinity-frame event energies and rates (P0/P1/P2 averages, μ_n event weights);
  2. global integrals: N_i, t, Z, Ltilde, Cstar;
  3. surface quantities: radius, surface gravity, envelope, Tsurface_inf.
- Every one of these is a separately certified group.
- My partial inspection of `DirectEnergyLedger.cpp`, `FrozenControlledBnvRunContext.cpp` and `FrozenBnvValidity.cpp` found two metric uses: an identity-string check (`MetricIdentity()`), and a local/infinity consistency check on per-event `exp_phi` against a 64-ulp tolerance. Neither is a nodewise r, m, ν, λ consumer in the right-hand side.
- If a complete consumer inventory confirms this, the nodewise metric row is **redundant for the experiment's scientific validity**. It is required only by the literal §11 table.
- Treating it as "covered through certified consumers" is a narrowing that needs an explicit owner decision under D. It must not be called full BA13 qualification. The package says the same.

**If the owner wants full qualification, the minimal observer is:**

| Item | Specification |
|---|---|
| Quantities | r, m, ν, λ (and e^ν if a consumer uses it) |
| Coordinates | The 2,049 original production support coordinates, including the terminal surface node at the production cut (an X2 segment endpoint, never extrapolated) and the pinned neutron-edge node under decision C's convention (physical edge versus pinned cell) |
| Stars | X2 at targets 0, 10, 20, interpolated in x. Withheld targets 5 and 15 serve as prospective nodewise checks. |
| Implementation | Passive collocation dense output. Stage solves, step sequence and integration are untouched. X2 primary outputs must be byte-identical with and without the observer. |
| Required interpolation accuracy | Nodewise reference error ≤ 1 % of the cap, i.e. ≤ 1.3×10⁻⁷ scaled. GL4 (s = 2) collocation output is uniformly O(h³) in the interior, with no fourth-order claim. At N = 2048 per segment the expected error is ~10⁻⁹–10⁻¹⁰. |
| Acceptance conditions | (i) Dense-output order test at withheld interior coordinates (stage midpoints) across N = 512/1024/2048: observed order in [2.5, 3.5], signed and consistent at all three targets; finest conditional tail ≤ 1 % of the cap. (ii) Dense output at step endpoints equals the stepped solution to ≤ 10⁻¹⁵ relative. (iii) Independent closure at the observed nodes: mass equation and metric relation residuals ≤ the tail. (iv) Nodewise P − I_X at all 21 samples is computed and charged in the joint inequality. |
| Prospective test (predeclared numbers before execution) | At targets 5 and 15: \|X_obs − I_X\| ≤ K_node h²/8 + bands at every node, with K_node from targets 0/10/20 times 10 |
| Computational limit | ≤ 120 CPU s total (benchmark plus 9 reference runs plus 2 holdouts) on one pinned core. Stop on any failed gate. No N = 4096, no second estimator, no further development. On failure, fall back to the consumer route. |

---

## 7. Final recommendation

**Disposition: B — ready for a narrower, explicitly conditional acceptance.** Full BA13 (A) is one bounded step away: the §6 observer. The consumer route keeps the claim at B, as an explicit narrowing of the §11 table.

Why the other dispositions do not apply:
- **C is not the right letter.** For the narrower claim, no further evidence is needed. What is missing is model statements and acceptance.
- **D is not indicated.** No production-method correction is needed. Nothing measured approaches a threshold.
- **E is not justified.** The evidence suffices for a qualification once the owner accepts the standard and conditions.

**The smallest remaining gap.** Write down the six error models from existing evidence (Appendix A), name the cut law and tangent object in decision C, and record the owner's itemized acceptance. Everything else is already measured or estimated, with multipliers of 10³–10¹².

**Can we establish an adequately justified frozen-background approximation for the controlled BNV experiment without a mathematical proof of every numerical operation?**

**Yes.** The evidence standard is the established solution-verification standard of computational science (code verification, solution verification and independent-implementation comparison, as in ASME V&V 20 practice), applied term by term:

1. **Code verification.** Manufactured-solution order (4.00025), the predicted centre orders, and thermodynamic, charge and β residuals.
2. **Solution verification.** Each referenced quantity has ≥ 3 resolutions in the asymptotic range with observed order within ±0.1 of formal at three targets. The estimate is the finest Richardson tail. Iterative and roundoff effects are shown subordinate by quantity-specific order cleanliness and precision ablations.
3. **Parameter convergence.** For the EOS quadrature, centre expansion and arithmetic precision: a measured sensitivity, plus a stated convergence law (exponential, known order, or ε-scaling).
4. **Independent-implementation agreement.** Production, X1 and X2 share only the physical constants.
5. **Prospective falsification.** Withheld targets with predictions committed before execution. This has been done and passed quantitatively.
6. **Named physical-regularity conditions** with explicit quantitative content and failure multipliers: C-CONT, C-SUPPORT and the slope hypothesis.
7. **Itemized owner acceptance** of each conditional in a registry committed separately from the assessor.
8. **Fixed qualification wording**, so that "qualified" never reads as "rigorously enclosed".

What makes this standard sufficient here is the drift structure (§3.1). The only consequential errors this standard cannot exclude are shared O(1) errors in a physical B-derivative. Those are model-specification errors, which no numerical proof could exclude either.

**Recommended path:**
1. **Owner adjudication** of Option A, interpretation (2) and A–E, with the §4 amendments, and of each Appendix A condition: accepted, rejected or amended.
2. **One mechanical static completion.** Populate the six null term types per Appendix A with values, classes and provenance. Commit the registry separately. Re-run the assessor, tests and mutation controls, ≤ 60 CPU s with no stellar evaluation. Verify against a predeclared expectation: every non-metric row PASS under the registry, utilization changes ≤ 10⁻³.
3. **Metric route.** Commission a consumer inventory (code reading only). If it confirms consumer coverage, record the owner's D-narrowing. Otherwise run the §6 observer under its 120 CPU-s cap.
4. **Short acceptance check** of the completed sidecar. Then the owner decides whether to authorize the already scoped certificate-dependent controlled-BNV pilots.

**Stop rules:**
- If the owner rejects validated numerical estimates, no qualification is attainable within scope. Only ADR-0021 Option B, historical-method reproduction, remains.
- If the owner rejects C-CONT, only a 21-point sampled claim remains, which is not usable by a continuous runtime monitor. Stop.
- If the consumer inventory finds a direct nodewise metric consumer, the observer becomes mandatory for the experiment.
- In no case: another solver campaign, holdout, tolerance ladder or root-refinement campaign.

---

## 8. Reviewer artifacts and reproduction

| File | Purpose |
|---|---|
| `scripts/authenticate.py` → `outputs/authentication.json` | 120 identity, history, boundary, manifest, receipt, mutation, handoff, certificate and registry checks |
| `scripts/quantitative_review.py` → `outputs/quantitative-review.json`, `outputs/quantitative-review.txt` | A: caps and fit norms. B: holdout curvature-prediction accuracy. B2: common-abscissa explanation of the residual. C: raw-output tangent recomputation (60 digits). D: evidence scales, requirements and multipliers. E: muon sector and onset structure. F: endpoint sliver. G: cut law. |
| `outputs/rerun-test_ba13_{recovery,headroom,closure,evidence,final_decision}.out` | Re-run of the 55 unit tests with toolchain Python `-I -B` |
| `outputs/no-write-check.txt` | `find -newer` over the candidate worktree and every evidence and handoff tree (empty), plus the worktree status including ignored files |
| `SHA256SUMS`, `SHA256SUMS.check.txt` | Manifest of this directory, and the `sha256sum -c` output |

Re-run from the review directory:

```
cd /mnt/sdd/zaki/CompactStar/linux-transition-20261010/ba13-final-decision-independent-review-20261030
../toolchain/python/bin/python -I -B scripts/authenticate.py outputs/authentication.json
../toolchain/python/bin/python -I -B scripts/quantitative_review.py outputs/quantitative-review.json > outputs/quantitative-review.txt
```

Both scripts take about 5 s on a login node.

**Changed files in any repository: none.** No ADR, certificate, registry, branch or evidence tree was touched.

---

## Appendix A — Proposed quantitative statements for owner consideration (not accepted)

These are a reviewer's proposals. Each needs explicit owner acceptance, rejection or amendment. Values are per row, from retained files. "Multiplier" is the reviewer's robustness figure, not part of the condition.

**S — Validated-numerical-estimate standard** (amends decision B). An estimate is class B if all of the following hold:

- (a) It derives from ≥ 3 resolutions or parameter levels with signed, consistent increments and observed order within ±0.5 of formal at every reported target.
- (b) Its value is the finest Richardson tail, or the measured sensitivity under a stated convergence law.
- (c) Non-asymptotic iterative and roundoff effects are shown subordinate, by quantity-specific order cleanliness or a precision ablation.
- (d) An independent implementation or a prospective prediction corroborates it where one exists.

Decision-value stop rule (an investigation criterion, not an admissibility condition, consistent with §4): stop refining any estimate once it is ≤ 10 % of its single-axis requirement.

**C-REF — Reference error.**
- X2's error in each referenced quantity equals its charged GL4, EOS and centre estimates, plus a stage and arithmetic remainder no larger than that quantity's finest Richardson increment.
- No numerical error is shared between X2 and production beyond the physical constants and the model specification.
- Any shared residual is ≤ 10 × the observed production–X2 agreement for that quantity. For Ltilde_Me, the X1–X2 primary-output agreement of 9.1×10⁻¹² is used instead.
- Multiplier ≥ 4×10⁴ (arbitrary-error treatment). Shared B-derivative errors must reach 85 % of the physical sensitivity (N_μ) to matter.

**C-TAN — Tangents.**
- |t_true − t_X2| ≤ (finest GL4 tail on the ratio) + 2 × (narrowest-width stencil tail) + (propagated EOS and centre sensitivities) + (cross-geometry moving-cut residual, if decision C certifies the path tangent).
- All terms are summed in absolute value, with no cancellation. The denominator relative error is ≤ 10⁻¹⁴.
- t_μ total ≈ 1.6×10⁻¹¹ absolute. Multiplier 1.8×10³.

**C-CONT — Continuous response.**
- For each group, |Q_xx| ≤ K_group on [−10⁻⁶, 0], with K_group as recorded in `conditional-intervals.json`: 10 × the maximum of the target and stencil curvatures.
- Fixed species topology; B_H ≥ 0.999 × the minimum sampled secant (observed variation 5.8×10⁻⁶).
- The endpoint sliver is closed by the secant rate plus K·sliver²/2.
- Validated by both holdouts (prediction within bands for all 15 groups). For N_μ the curvature would have to exceed K by ≥ 7×10⁶, i.e. the measured curvature by 7×10⁷, to consume the joint right-hand side.

**C-INPUT — Inputs.**
- The abscissa design errors are as recorded per gap in `ba13-closure` m3, with endpoint sums ≤ 5.1×10⁻¹² in x, and the B₀ normalization error is as charged.
- The fit response is (I − H)e with f = 5.185 (all rows), plus a second-order design remainder ≤ (δx/Δx)² × the maximum residual.
- The historical weights are part of the frozen fitter, not error inputs.

**C-SUPPORT — Support.** No species threshold is crossed, given C-CONT's B_H condition: H_c − H_μ = 0.0197 against a window of 3.97×10⁻⁷.

**C-CUT — Surface law** (amends decision C).
- Q is defined with ln P_s(x) linear through the retained cuts (slope 8.3839 per unit x).
- The historical cuts deviate from it by ≤ 2.1×10⁻¹², charged through the surface sensitivity.
- The certified tangent is the one named in decision C.

**C-MUN — μ_n_profile.** The equilibrium identity μ_n^∞ = μ_B^∞ at η = 0 supplies the reference for all 2,049 nodes, within the equilibrium/reference definition only. The η-dependent actual potential is evaluated separately, as §8 requires.
