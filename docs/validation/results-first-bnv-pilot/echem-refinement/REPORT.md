# Echem reconstruction study — full handoff to Fable

Date: 2026-11-09. **Status: UNQUALIFIED numerical candidate; investigation stopped
for owner-requested handoff. No Stage A retry or matched control was executed.**

## Finding and decision needed

The original failure is a **chemical-state tolerance floor propagated into
whole-star chemical free energy, Echem**, not a measured factor-eight integration
error. The bounded study reproduces it. Tightening the test-only O2 chemical atol
to **5e-27 MeV** makes the existing diagnostic budgets pass at all 13 sampled
checkpoints, with maximum `U/(0.20 F) = 0.7330768176357471`. This is the largest
predeclared tested tolerance meeting the study's 0.80 utilization margin.

**That is not a feasibility PASS.** Every one of the six tested tolerances
produced the same serialized binary64 O2 state at each checkpoint. The reduced U
therefore reflects the configured floor, not a demonstrated improvement in the
integrated state. Independent local witnesses agree at roughly 2e-15 relative
Echem, but a complete coupled accumulated-roundoff bound and the corresponding
convergence assessment remain unestablished. The conditional chemical-arithmetic
audit does not fill that gap. No production tolerance or acceptance policy was
changed, and the conditional Stage A retry authorization was not exercised.

Fable's immediate task is to assess whether the retained traces and independent
witnesses justify the smaller floor under the existing accuracy policy, including
coupled thermal/chemical roundoff. There is no evidence here that a new solver or
relaxed energy criterion is necessary. Nor is there a measured impossibility
limit that would justify declaring binary64 inadequate.

## Provenance and scope

Repository: `/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot`.
Branch: `candidate/results-first-bnv-pilot-20261108`.
The clean starting source was `56dd1959a945bc13de574e83ba664998ae87e3b7`, a descendant
of `a8d48cf`; all 15 local branch tips were ancestors and nine worktrees existed.
The executed study binary was built from
**`76b1907598881dab5fe72039f037a8e206e3f9d8`**. Earlier study commits were `cbc8b03`
(predeclaration/harness) and `e6ac1d1` (build isolation). The report and offline
analysis are committed separately; obtain their exact containing commit with
`git log -1 --format=%H -- docs/validation/results-first-bnv-pilot/echem-refinement/REPORT.md`.

Authority: owner's bounded reconstruction/refinement request, then explicit
request to stop and finalize the handoff; GOVERNANCE §§2–3 and §5, ADR-0017 §§5–7
and its breakpoint amendment, INV-11f and INV-15. Changes are a numerical-method
**test candidate**, engineering harness/analysis, documentation and generated
evidence. AI-authored analysis remains unverified pending domain-owner review.
The observed baseline was preserved; no pre-baseline scientific-refactoring
exception is invoked. No alternative scientific policy is selected or proposed
for implementation by this report.

The study reused job 19174's source accepted history and five completed original
reference observation tables. It created local reconstructions from retained
accepted left states, not new source/control trajectories from t=0. Original
numerical archives, physical inputs, owner/runtime guards and snapshots were
retained. Snapshot SHA-256 identities:

- Source: `c6a29044e8330f079a1593be166238e11eafb66745186c5322c04c272358ec18`.
- Control: `0588ddf4518b3dec0f1cf103dcb8a3b356cf4631f2394df41b6a2d6baadf4ceb`.

No stellar solve, BA13 reinvestigation, physical-ledger change, master merge,
certificate promotion, Stage B or longer campaign occurred. Static owner
assembly in the local harness is not a stellar solve. No Mac artifact was made.

## Original refusal, reproduced

The full prior diagnosis, reporter regression tests and diagnostic-only replay
are in [checkpoint-diagnosis/REPORT.md](../checkpoint-diagnosis/REPORT.md).
The failed channel is **Echem**, whole-star chemical free energy in MeV.
Time is **131,490 s**, checkpoint **1 on the 240-interval grid**, equivalently
**8 on the 1920-interval grid**. This distinction prevents an index ambiguity.

| Quantity | Original and reproduced baseline |
|---|---:|
| Echem, both O1 and O2 | 2.8395321373652751e28 MeV |
| U | 4.5873384329419162e17 MeV |
| F | 2.8395321373652749e17 MeV |
| Allowed 0.20 F | 5.6790642747305504e16 MeV |
| U/F | 1.6155261539664711 |
| U/(0.20 F) | 8.0776307698323553 |
| O1/O2 energy disagreement | 0 MeV |
| Propagated state contribution | 4.581708933407703e17 MeV |
| Explicit disagreement/ulp contribution | 5.62949953421312e14 MeV |

The propagated state contribution is 99.8773% of U. All three state channels
pass; only Echem fails among 18 diagnostic channels. The five retained global
references differ from O2 by at most 3.0786325577728e13 MeV at this time, 0.0005421
of the permitted budget. This is measured agreement, not certified reference
uncertainty. The previous replay reproduced the main accepted history byte for
byte; the new study reproduces the failing local diagnostic and its utilization.

## Exact diagnostic being tested

For Echem, the existing implementation is:

```text
M = max(abs(Echem(O1)), abs(Echem(O2)))
F = max(1e-11 * max(1, M), 64 ulp(M))
u_i = 2 max(abs(O2_i-O1_i),
            atol_O2,i + rtol_O2 max(abs(O1_i), abs(O2_i)),
            64 ulp(max(abs(O1_i), abs(O2_i)))) + U_event,i
P = 2 max over eight state-box corners abs(Echem(O2 +/- u)-Echem(O2))
U = 2 max(abs(Echem(O2)-Echem(O1)), 64 ulp(M)) + P
Require U <= 0.20 F.
```

The local witness/event radius is the original doubled maximum of half-span O2
and independent-witness disagreements, plus its 64-ulp term. The fixed quadratic
energy algebra is

```text
Echem = [Z11 eta_e^2 - (Z01+Z10) eta_e eta_mu + Z00 eta_mu^2]
        / [2 (Z00 Z11 - Z01 Z10)]  MeV.
```

Echem uses energy normalization, not the gross-power scale G. At these nonzero
states the accepted relative energy uncertainty is `2e-12`. The tolerance floor
is a conservative uncertainty construction; it is not a measured local RK error.
See the prior report's code citations and the serialized per-channel budgets.
Neither the formula, safety factors, solver algorithms nor thresholds changed.

## Predeclared numerical matrix and actual execution

[PREDECLARATION.md](PREDECLARATION.md) was committed before numerical execution.
The 13 finest-grid indices were
`1,2,4,8,16,32,64,128,256,512,960,1440,1919`, at `31557600 * index/1920` seconds.
The six O2 chemical atols were `1e-24,1e-26,5e-27,1e-27,1e-28,1e-29` MeV.
O2 rtol=1e-13 and x atol=1e-18, O1 settings, built-in witnesses and original
algorithms were unchanged. Test-only configuration injection changed the stored
O2 atol on a nonconst reconstructor; no production numerical library was edited.

There are **78 completed candidate reconstructions**. Each checkpoint also has
six independent local witnesses: rk8pd and RKF45 at maximum steps span/4, /8,
/16, rtol=1e-15 and atols `(1e-20,1e-30,1e-30)`, for **78 witnesses**. These
restart from the same retained left state. They do not independently validate
the main prefix. All witness rejection counts are zero. Raw traces preserve RHS
stages, diagnostic evaluations, witness accepted steps and local error estimates.

The prospective checks required original gates, maximum Echem utilization ≤0.80,
observed discrepancy plus twice independent refinement spread within the existing
uncertainty envelopes, and contraction ≤0.5 above an independently established
roundoff bound. Missing roundoff/reference evidence was explicitly UNQUALIFIED.

| Job | Actual outcome |
|---|---|
| 19177, build | Failed before execution: `redefinition of int main`; helper extraction fixed the include/macro collision. |
| 19178, build | Succeeded; 23.24 s, child maximum RSS 495948 KiB. |
| 19179, harness | Failed before first reconstruction: `vector::_M_range_check`, index 18446744073709551615, size 12. The observation schedule omitted its required t=0 sentinel. Output tables contain headers only. |
| 19180, build | Succeeded after sentinel fix; 22.90 s, 495072 KiB. `--check-history` passed 13 strict-interior brackets plus initial sentinel, with retained local errors loaded and no solver calls. |
| 19181, study | Completed the declared matrix, process exit 0; 231.99 s, 434688 KiB. Collection completion is not qualification. |

The failed build and pre-reconstruction harness attempt are preserved. They were
harness repairs, not additional tolerance searches or Stage A retries. Jobs used
Slurm with 4 GiB, 30 minutes and one logical execution CPU; recorded process
affinity is `[0]`. Slurm reports two reserved CPU threads but the process was
pinned to one logical CPU. Resource reports are Linux maximum child RSS, not a
simultaneous aggregate-memory measurement. Final `scontrol` lookups at handoff
returned `Invalid job id specified`; the retained start-time allocation records,
completion logs and measured resource receipts are the available evidence, not
fresh scheduler completion records.

## Results and uncertainty budgets

| O2 chemical atol (MeV) | Maximum U/(0.20 F) over 13 points | Configured gates and ≤0.80 margin |
|---:|---:|---|
| 1e-24 | 61.64009742541691 | Fail |
| 1e-26 | 1.0388223570426696 | Fail |
| 5e-27 | 0.7330768176357471 | Pass, conditional evidence only |
| 1e-27 | 0.4878918182025361 | Pass, conditional evidence only |
| 1e-28 | 0.43290718472256773 | Pass, conditional evidence only |
| 1e-29 | 0.43000168663530697 | Pass, conditional evidence only |

The worst point is index 1, t=16,436.25 s, except for 1e-29 where it is index
1440. Thus fixing only the original index-8 failure would miss a stricter early
checkpoint. All six tolerances give **one identical O2 state per checkpoint**.
The table must not be described as convergence of the solution with tolerance.

For the smallest prospective change, atol=5e-27 MeV:

| Quantity | Earliest tested point, 16,436.25 s | Original failed point, 131,490 s |
|---|---:|---:|
| Echem (MeV) | 4.4367694567187725e26 | 2.8395321373652751e28 |
| U (MeV) | 6.50498566782976e14 | 2.6370686880579584e16 |
| F (MeV) | 4.4367694567187725e15 | 2.8395321373652749e17 |
| 0.20 F (MeV) | 8.873538913437545e14 | 5.6790642747305504e16 |
| U/(0.20 F) | 0.7330768176357471 | 0.46434915339694355 |
| U/Echem | 1.466153635271494e-12 | 9.286983067938871e-13 |
| Propagated state term (MeV) | 6.41702473760768e14 | 2.580773692715827e16 |
| Explicit disagreement/ulp term (MeV) | 8.796093022208e12 | 5.62949953421312e14 |
| O1/O2 Echem disagreement (MeV) | 0 | 0 |

These totals include the unchanged relative-tolerance, witness/event and ulp
contributions through the existing box formula. Nothing was removed to obtain
these utilizations. The roughly 0.43 limiting utilization as atol is reduced
illustrates that other contributions remain. At the earliest point, reevaluating
the recorded corners at 100 decimal digits gives a propagated term
6.4136849356101597e14 MeV versus 6.41702473760768e14 from the binary64 diagnostic;
this local subtraction/algebra rounding is much smaller than the allocation.

Across the 13 checkpoints, the measured independent local-witness Echem
disagreement is at most **2.0135258411715982e-15 relative**. The measured envelope
(disagreement plus twice maximum adjacent/finest-method spread) uses at most
**0.0030202887617573975 of the allowed energy budget**, and at most
**0.007044797687861272 of candidate U**. The corresponding maximum state-envelope
fraction of its propagated radius is **0.014387835356780667**. These are observed
numbers, not certified bounds on the physical trajectory's error.

## Roundoff, cancellation and what is still missing

The offline assessment reevaluated Echem at 100-digit precision from the exact
serialized binary64 states and frozen Z coefficients. Maximum relative center
algebra discrepancy is **2.900224229261538e-16**. The earlier original-failure
audit found numerator cancellation condition 1.04052 and determinant condition
1.00114. Catastrophic energy-algebra cancellation is not supported as the cause.

The exploratory chemical roundoff script processed **390 RHS-bearing contexts**
(all supported by its trace parser). Context IDs also include diagnostic-only
contexts, so they are not a contiguous 390-ID sequence. It uses recorded RK stage
residuals, GSL rational tableau coefficients, a high-precision chemical RHS and
an enlarged chemical Jacobian box to propagate conditional arithmetic estimates.
Its maximum relative Echem enclosure is **2.262268871997437e-15** across contexts,
and **5.882253107010768e-16** for O2 contexts.

**This audit is an unverified candidate, not a full roundoff certificate.** It is
conditional on recorded thermal-stage inputs. It excludes thermal-stage
accumulated roundoff, main-prefix error, truncation and reference interpolation.
It freezes the serialized input constants and does not certify their physical
accuracy. Its inferred step arithmetic and Jacobian enclosure also require
independent review. A 100-digit residual evaluation is not a 100-digit independent
coupled trajectory. The scripts deliberately retain false Stage A gate flags;
these are explicit missing-evidence status, not a discovered solver failure.

Observed witness differences do not uniformly contract by 0.5: for example at
index 8 the rk8pd adjacent Echem differences are 1.3194139533312e13 and
3.5184372088832e13 MeV. They are small and irregular, compatible with arithmetic
noise, but the predeclared exemption below a justified roundoff bound has not
been established for the full coupled problem. Agreement alone cannot supply it.
The actual reliable lower tolerance limit has therefore **not been determined**.

The 13-point local study does not certify all checkpoints, all main tiers/caps,
global-reference uncertainty, R18/R20, main accumulated roundoff, checkpoint/
restart, or the matched source/control signal. It contains no qualified one-year
temperature figure. Presenting a temperature plot as qualified would overstate
what was executed.

## Scientific interpretation and smallest next action

Classification remains **E: conservative state-floor/energy-budget incompatibility
near equilibrium**. A genuine nonlinear integration defect (A), reference-driven
refusal (B), dominant cancellation (C), or an incorrect normalization/formula (D)
is not demonstrated. The missing failure-channel observability defect was fixed
in the prior diagnostic work; it is distinct from this numerical question.

Because Echem is quadratic, a fixed imbalance error floor produces a roughly
`atol/abs(eta)` relative energy estimate away from zero. Near exact zero the
corner contribution becomes quadratic in atol. The implemented criterion already
has an absolute part: `F >= 1e-11 MeV`, hence permitted uncertainty at least
`2e-12 MeV`. It is not a purely relative condition all the way to zero. The
owner-certified exact-zero control exception does not justify suppressing the
floor in a nonzero driven source. This study supplies no quantitative necessity
for a different absolute-plus-relative policy, and none was implemented.

Recommended focused handoff sequence:

1. Review the retained local traces and conditional audit, specifically the
   thermal-to-chemical error propagation and the independent convergence/roundoff
   requirement. Determine whether the remaining evidence can be completed from
   these artifacts. Do not infer a pass from unchanged states or smaller atol.
2. If that verification succeeds, the smallest identified numerical candidate is
   **only O2 chemical atol 1e-24 → 5e-27 MeV**. Preserve rtol, x atol, all algorithms,
   snapshots and `U <= 0.20 F`. This setting was tested, not adopted in production.
3. Only after publishing the missing verification, use the owner's one conditional
   Stage A retry: qualify source first, then matched zero-source control, under
   the original [Stage A plan](../native-regeneration/RUN_PLAN.md). Every reference,
   Echem/Tinf/imbalance, R18/R20, restart and accumulated-roundoff requirement remains.
   Each stage retains one logical CPU, 4 GiB and 30 minutes. No Stage B.
4. If verification instead reveals a real precision obstacle, report its measured
   size and affected terms. Any scientifically necessary accuracy-policy change
   needs a quantitative energy-accounting proposal and explicit owner approval.

This report does not authorize new experiments beyond the owner's existing
conditional scope. Work stopped for Fable's review rather than expanding the
investigation into another campaign.

## Files, evidence and validation

Study files committed in the preceding three commits: `PREDECLARATION.md`,
`tests/diagnostics/checkpoint_refinement.cpp`, and
`tools/linux/native_qualification/{build_refinement.py,run_refinement.py,refinement.sbatch}`.
This handoff adds this report, `REPRODUCE.md`, compact evidence/checksums, and the
already-executed offline scripts `assess_refinement.py` and
`audit_chemical_roundoff.py`. There is no production solver or policy edit.

Start with [evidence/summary.json](evidence/summary.json),
[evidence/assessment.json](evidence/assessment.json), and
[evidence/candidates.tsv](evidence/candidates.tsv). The assessment's external
filename was `partial-assessment.json`; its final content has `complete=true`.
The committed `assessment.json` is a byte-identical copy, not a new experiment.
Representative complete checkpoint records include baseline index 8 and candidate
indices 1 and 8. All points' raw records remain externally retained.

External root:
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/echem-refinement-20261109-r1`.
[external-SHA256SUMS](evidence/external-SHA256SUMS) indexes 354 files (13,608,659
bytes), including both failed harness attempts, all traces, binary and manifests.
Historical failed-attempt manifests can name paths subsequently reused by a
successful attempt; use this content inventory and archived attempt directories
rather than interpreting those old paths as current identity assertions.
Original job 19174 and prior diagnostic replay evidence remain untouched in
their separate existing directories.

Validation actually performed: compiled successful harness, no-solver bracket
regression, complete bounded matrix, offline high-precision assessment and
conditional audit. At handoff, all 25 files in the successful build manifest
were rehashed successfully; both offline scripts passed Python AST syntax checks;
78 candidate and 78 witness rows and the reproduced baseline were checked.
No independent unit/regression validation of the new conditional-roundoff
algorithm was completed. The earlier reporter regression evidence remains in
the prior diagnosis. Handoff checksum and whitespace checks cover packaging only. The two raw Slurm
allocation receipts retain their original trailing whitespace; authored text
passes the whitespace check with those byte-preserved receipts excluded.
No additional numerical experiment was run to finalize this report.

**Unresolved outcome:** a promising smaller floor with substantial observed
agreement, but incomplete justification for full numerical qualification.

## 2026-11-09 addendum — bounded D1 convergence resolution

The preceding report and its findings are historical. The original offline
D1 global-grid calculation remains [FAIL at C1](d1/REPORT.md#L1); its receipt,
implementation and evidence were preserved byte-for-byte. The owner then
authorized one bounded diagnostic to resolve that reference failure, starting
at commit `871dbd77230f9e21b54ea0129473fe658236fd56`.

The [resolution receipt](d1-resolution/REPORT.md#L1) records
**PASS_CONDITIONAL for the 13 existing local reconstructions**. The method and
unchanged acceptance criteria were [predeclared](d1-resolution/PREDECLARATION.md#L1)
and published at `7842938` before job 19185. Aligning RK4 steps to every retained
Hermite segment restored fourth-order convergence: maximum C1 ratio
0.06250000005471215. No thermal interpolant, input, production solver or tolerance
was changed. All original inputs and exact initial-state conversions passed
authentication/fidelity checks; all 17 analytic tests passed before the fixed
156 local reference integrations.

Original C2/C3 pass. A stronger comparison, including a conservative bound for
time-dependent thermal perturbations within the assumed band, also passes:
maximum state ratio 0.004303188584056 and energy ratio 0.001660771997586. Actual
errors, both reference contributions and original budgets are reported
separately in [the tables](d1-resolution/evidence/TABLES.md#L1). This qualification
is conditional on a uniform thermal band and continuation of the observed
refinement contraction; finite sampled agreement is not an unconditional IVP
bound. C4 remains report-only. Main-prefix error and thermal-channel roundoff
obligations remain open, as do the excluded full-run qualifications.

The only proposed future production change is O2 chemical absolute tolerances
`(1e-24,1e-24) -> (5e-27,5e-27) MeV`, pending the owner's go-ahead. No such change,
Stage A retry, Stage B, stellar solve, master merge or certificate promotion
occurred. Any later authorized Stage A attempt must apply the diagnostic to its
actual output checkpoints. See the receipt for complete scope and provenance.
