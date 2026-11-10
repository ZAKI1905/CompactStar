# BNV-THERMO-01-D1 — refinement-report addendum

2026-11-09. **AUTHENTICATION FAIL; chemical reconstruction unassessed.**
This is a failure receipt, not D1 qualification or permission for a production
tolerance change. The original refinement report and all original numerical
inputs remain unchanged.

The maintained 110-digit Decimal diagnostic was implemented and invoked on an
EKU compute allocation. It stopped because
`BNV-THERMO-01_FABLE_ADJUDICATION.md` and its independent identity were not
available. No replacement adjudication was fabricated from the task prompt.
The exact failure is retained in
[job-19183/failure-details.json](evidence/job-19183/failure-details.json#L1).
Every one of the 13 point records explicitly says `NOT_RUN`; missing reference
values are JSON nulls, not zeros or fabricated 40-digit answers.

## Authentication and preserved evidence

The clean starting worktree was
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot`,
branch `candidate/results-first-bnv-pilot-20261108`, exactly at reviewed commit
`ad4e91aee3837b8b07e7ed723afdd40b765011fe`. All 15 local branch tips were ancestors;
nine worktrees existed. The committed records identify the actual retained study
as `/mnt/sdd/zaki/CompactStar/linux-transition-20261010/echem-refinement-20261109-r1/study/output`.
See [preflight.json](evidence/preflight.json#L1) and the original
[REPORT.md](../REPORT.md#L30).

Job 19183 freshly hashed **all 354 external files** against the original,
reviewed-commit `evidence/external-SHA256SUMS`. They all match. This includes the
four summary/trace tables and all 13 sets of `checkpoint.json`, `state.tsv`, and
`diagnostics.tsv` at `atol-5e-27`. This is a new cluster check, not reliance on the
reviewer's inaccessible traces. The source audit, conditional-roundoff output,
external inventory and physical-input authentication record were also bound to
their exact reviewed-commit blobs. The coefficient and tangent tables are
**outside** the 354-file inventory; their hashes were instead obtained from the
reviewed commit's `native-regeneration/evidence/authentication.json` and verified
against the actual cluster bytes. The complete chain is in
[results.json](evidence/job-19183/results.json#L1), under `authentication`.

The missing review was absent from tracked paths at the reviewed commit and
searches under `/mnt/sdd/zaki/CompactStar`, `/home/zaki`, and `/tmp`. A narrow
connected Google Drive search returned no match. No access to the office Mac or
Scriptorium is claimed. The owner was asked for the actual file/path and trusted
digest; none was available during execution. Consequently the review itself
could not be authenticated or preserved as durable evidence.

A post-job rehash confirmed all **360 checked file identities** remained intact:
[post-execution-input-check.json](evidence/post-execution-input-check.json#L1).
Original trace contents, modes, physical inputs, audit code, and audit results
were not modified. Original checksum rows remain intact; this addendum's
manifests add the new artifacts.

## Implementation and validation actually performed

[verify_chemical_reference.py](../../../../../tools/linux/native_qualification/verify_chemical_reference.py#L1)
contains an independent classical RK4 implementation with rational coefficients,
the fixed `N=64,128,256,512` matrix and offsets `-1e-12,0,+1e-12`. It reconstructs
the mathematical reaction/RHS definition from the authenticated audit, without
importing its recursion or any GSL stepping code. Retained numbers are decoded
to their original binary64 values and lifted with `Decimal.from_float`.
Calculation precision is 110 digits. Time is in seconds, chemical state and
energy are in MeV, and `x=ln(Tinf/1e8 K)`. Each local integration uses
`tau=t-t_left`, with the exact retained left state, rather than starting at t=0.

The implementation checks all local RHS rows with the fixed two-ulp limit;
constructs the accepted rk8pd/8 Hermite path using exact node/derivative matches;
checks rk8pd/16, RKF45/16 and O2 residuals against `rho<=1e-13`; retains both
uncertainty contributions; and implements C1–C3 and report-only C4. It writes
failure details before nonzero exit. **Those retained-data numerical paths were
not executed**, because authentication must precede them. Their implementation
is not evidence that the 13 intervals pass.

[test_chemical_reference.py](../../../../../tools/linux/native_qualification/test_chemical_reference.py#L1)
passed **11/11 analytic tests** on the compute node: exact binary64 lifting and
ulp spacing; constant and cubic-forcing integrals; exponential-solution RK4
order; quadratic energy with a cross term; four-corner propagation; cubic
Hermite reproduction; zero-difference convergence and the small-difference
floor; reaction symmetry; and the exact two-ulp fidelity boundary. See the
[test log](evidence/job-19183/tests.log#L1). AST parsing, shell syntax and authored
whitespace were also checked. These tests use synthetic data and do not qualify
the retained scientific outputs.

| Job | Outcome | Resource evidence |
| --- | --- | --- |
| 19182 | 11 tests passed; launch failed before authentication because compute-node PATH lacked Git | [receipt](evidence/job-19182/results.json#L1), [completion](evidence/job-19182/completion.txt#L1) |
| 19183 | Retained toolchain added to PATH; 11 tests passed; authentication failed only on missing adjudication | [receipt](evidence/job-19183/results.json#L1), [completion](evidence/job-19183/completion.txt#L1) |

Both jobs used `phyn004`, one logical execution CPU (affinity `[0]`), 4 GiB and
a 30-minute ceiling. Slurm reserved two hardware threads; the execution process
was pinned to one. Each allocation completed in one second. Job 19183's Python
process reported 22,272 KiB maximum RSS and 0.557 seconds elapsed; this is not an
aggregate measurement of all child-process memory. Slurm accounting is disabled;
the preserved `scontrol` allocation/completion records are the resource evidence.
The only retry repaired a launch-environment defect before authentication. No
numerical setting changed and no retained-interval integration ran in either job.

## Disposition and outstanding obligations

| Gate | Disposition |
| --- | --- |
| Reviewed commit and supplied external/physical inputs | PASS |
| Governing adjudication authentication | **FAIL — unavailable file and trusted identity** |
| RHS fidelity, thermal residual, C1, C2, C3 | NOT_RUN for all 13 points |
| C4 audit validation | NOT_RUN; report-only when executed |
| D1 local reconstruction qualification | **NOT ESTABLISHED** |

No point/component magnitude is reported for the failed authentication gate:
it concerns a missing governing input, not a measured chemical error. The
machine-readable records retain that distinction. No radius, allowance,
criterion or contribution was changed; `U <= 0.20 F` remains in force.

For a later properly authenticated execution, interpretation remains conditional
on a continuous thermal band and an adequate chemical-sensitivity estimate.
Agreement at sampled witness nodes alone does not prove a uniform enclosure.
The two constant-offset runs are sensitivity probes, not demonstrated extrema
over all time-dependent thermal perturbations. The implementation retains these
limitations explicitly. A C4 ratio above one invalidates that audit enclosure
as a bound; C4 is not an independent D1 acceptance gate.

Even a future D1 PASS would apply only to these 13 local intervals under the
stated assumptions. The main-prefix error from t=0 and thermal-channel roundoff
obligations remain open. It would not establish every Stage A checkpoint,
complete paired-run R18/R20 or restart qualification, chemical heating in the
one-year temperature observable, or Stage B readiness. A subsequent authorized
Stage A run must apply the diagnostic to its actual output checkpoints.

The conditional next-step proposal remains O2 chemical absolute tolerances
`(1e-24,1e-24) -> (5e-27,5e-27) MeV`, with all other numerical settings unchanged.
**This failure receipt does not justify that change or request its go-ahead.**
No production change, Stage A retry, thermal retry, Stage B, master merge or
certificate promotion occurred.

Change classes: numerical-method **test candidate**, engineering, documentation,
and generated evidence. Authority: the owner's bounded D1 task; GOVERNANCE
§§2–3 and §5; ADR-0017 §§5–7; INV-11f and INV-15. No production algorithm is
refactored and no pre-baseline exception is invoked. AI-authored code remains an
unverified scientific candidate pending domain-owner review. Files added are the
three maintained diagnostic/test/Slurm files, this addendum, reproduction notes,
compact receipts and manifests; the parent refinement manifest gains new rows.
