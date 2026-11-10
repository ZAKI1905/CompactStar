# BNV-THERMO-01-D1 — completed verification receipt

2026-11-09. **D1 FAIL — C1 reference convergence fails at points 1, 64 and 1440.**
Authentication, RHS fidelity, thermal residuals and the computed C2/C3 inequalities
pass. These do not override C1. The local chemical reconstruction box is therefore
**not qualified by D1**, and this receipt does not support a Stage A go-ahead.

The fixed calculation was completed once in compute job **19184**, using Python
Decimal at **110 digits**, classical rational-coefficient RK4, exactly
`N=64,128,256,512` and offsets `b=-1e-12,0,+1e-12`. All **156 integrations** across
the 13 retained intervals are preserved. No numerical setting was changed after
seeing results, and no further numerical execution was attempted.

## Exact failures

Both chemical components fail C1 at each of these three points for all three
thermal offsets. There are 30 failed contraction comparisons. The table gives
the central `b=0` results, with `d1=|eta128-eta64|`,
`d2=|eta256-eta128|`, `d3=|eta512-eta256|`. Ratios must be at most 0.5 unless the
new difference is below `1e-60*|eta|`. None of these differences is below that
floor. Full 110-digit values and all offsets are in
[evidence/summary.json](evidence/summary.json#L1) and
[evidence/job-19184/failure-details.json](evidence/job-19184/failure-details.json#L1).

| Point | Time (s) | Component | d2/d1 | d3/d2 | d3 (MeV) | Floor (MeV) |
| --- | ---: | --- | ---: | ---: | ---: | ---: |
| 1 | 16,436.25 | eta_e | 0.5000017962425375 | 0.5000002245295105 | 2.509125069314903e-48 | 5.674142684618745e-74 |
| 1 | 16,436.25 | eta_mu | 0.5000025736801995 | 0.5000003217083688 | 7.915114588944062e-48 | 1.439316317796586e-73 |
| 64 | 1,051,920 | eta_e | 0.5067672101636893 | 0.5008346053624186 | 5.057726929099665e-43 | 3.631449596826166e-72 |
| 64 | 1,051,920 | eta_mu | 0.5079193205286363 | 0.5009744806070687 | 1.595412351470938e-42 | 9.211619003908873e-72 |
| 1440 | 23,668,200 | eta_e | 0.007465670868019653 | 3.623325227348866 | 3.368247284566130e-39 | 8.170676983626176e-71 |
| 1440 | 23,668,200 | eta_mu | 0.01627237687504062 | 1.357931987344474 | 1.033291333377201e-38 | 2.072587585674555e-70 |

The largest contraction ratio over all offsets is **3.6233252274994592**.
The small absolute differences do not authorize a replacement convergence
threshold. The other ten points pass C1. All retained witness nodes were kept,
including the very short final rk8pd/8 segments at these three points
(3.637978807091713e-12, 2.328306436538696e-10 and 3.725290298461914e-9 seconds).
Their presence is recorded, not asserted as a proven cause or removed to force
convergence. No altered interpolation or step partition was tested.

## Remaining measured gates

| Check | Actual result |
| --- | --- |
| Authentication | PASS: 354 external files, source-bound audit/inventories, coefficients, tangent and supplied review |
| RHS fidelity | PASS: all 14,796 rows / 29,592 chemical components; maximum errors 0.7404581109 ulp (e), 0.8022895880 ulp (mu), below the fixed 2-ulp limit |
| Thermal band sample check | PASS: maximum rho = 3.1101367210353021e-20, below 1e-13; delta_x remains exactly 1e-12 |
| C2 computed state inequalities | PASS at all points; maximum ratio 0.001370086143057981 (point 16, eta_mu) |
| C3 computed energy inequalities | PASS at all points; maximum ratio 0.0004982358353235805 (point 4) |
| C4 report-only comparisons | No ratio above 1; maximum 0.9865931756010835 (point 128, RKF45/4, eta_e) |

The largest O2 error against the proposed reference is 1.9680116805582353 ulp
(point 1, eta_e). Maximum reference-uncertainty contributions, over all points
and components, are 6.4428377934800896e-38 MeV from discretization and
2.6757870102213857e-26 MeV from the thermal offsets. These are reported
separately; no contribution was removed and `U <= 0.20 F` is unchanged.
C2/C3 are arithmetic comparisons against the calculated reference and specified
uncertainty estimate; because C1 fails they do **not** certify the true state.
C4 is report-only: a ratio above one would void that audit enclosure as a bound,
not independently fail D1. Here no enclosure is falsified against this reference;
that is not a proof of the audit bound or a cure for C1.

[evidence/job-19184/results.json](evidence/job-19184/results.json#L1) is the complete
machine-readable result: input hashes, 110-digit eta_ref, every integration,
contraction ratios and zero-difference handling, both u_ref contributions,
thermal nodes/residuals, recorded radii/energy allowances, O2 errors in MeV and
ulp, C2/C3/C4 ratios and explicit dispositions. Matching per-point JSON files
and [summary.tsv](evidence/summary.tsv#L1) support direct inspection.

## Authentication, execution and changes

The initial worktree was clean at reviewed commit
`ad4e91aee3837b8b07e7ed723afdd40b765011fe`, on
`candidate/results-first-bnv-pilot-20261108`, at
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot`.
All 15 local branch tips were ancestors and nine worktrees existed. The actual
study directory was recovered from committed records. The complete history and
initial stopped jobs are retained in [AUTHENTICATION_STOP.md](AUTHENTICATION_STOP.md#L1).

The owner subsequently supplied the original
[adjudication](review/BNV-THERMO-01_FABLE_ADJUDICATION.md#L1). Its exact bytes were
preserved and compared with the supplied attachment: SHA-256
`46201727575270055d237178b2ad0433ddc22837c0cc1fa646a8d34eddcb0066`.
This authenticates the copy against the owner-provided artifact; a separately
signed reviewer digest was not supplied and is not claimed. Its provenance and
Section D comparison are recorded in [review/provenance.json](review/provenance.json#L1).
The review was committed before execution, at `06066be`. It agrees with the
implemented fixed matrix and gates. Its arrival resolved the earlier missing
input without changing numerical settings.

Job 19184 rehashed all 354 original external files on the cluster rather than
assuming the reviewer had checked them. The coefficient and tangent tables are
outside that inventory and were authenticated through the exact reviewed
commit's `native-regeneration/evidence/authentication.json`. Source audit and
conditional-roundoff output were also checked against reviewed-commit blobs.
All 361 checked identities passed before and after the diagnostic; all original
inputs remained unchanged. The original refinement report is unchanged.

The maintained files are
[verify_chemical_reference.py](../../../../../tools/linux/native_qualification/verify_chemical_reference.py#L1),
[test_chemical_reference.py](../../../../../tools/linux/native_qualification/test_chemical_reference.py#L1)
and [chemical_reference.sbatch](../../../../../tools/linux/native_qualification/chemical_reference.sbatch#L1).
The implementation reuses the mathematical RHS definition, without importing
the audit recursion or GSL. Printed retained records are first decoded as
binary64, then lifted exactly into Decimal. Units are explicit: eta and Echem
in MeV, time in seconds, x=ln(Tinf/1e8 K). Each local origin is its exact retained
t_left; this is not a new trajectory from t=0.

All **11 analytic tests passed before the retained calculations**: binary64/ulp
handling, constant and cubic-forcing integrals, exponential RK4 order, quadratic
energy and four corners, cubic Hermite reproduction, convergence zero/floor
cases, reaction symmetry and the two-ulp fidelity boundary. The source script
SHA-256 is unchanged across all three jobs:
`92e1fb0a5c7e50a66dbd3ee7dd54cf084d4dc1705c74225978a6a6022e76b058`.
Job 19182 needed a PATH correction before authentication; job 19183 stopped on
the then-missing review; **19184 is the sole retained-interval numerical run**.

Job 19184 ran on `phyn004` with process affinity `[0]`, one logical execution
CPU, 4 GiB and a 30-minute ceiling. Its Python process used 9.452 seconds and
42,344 KiB maximum RSS (not aggregate child memory). Slurm reserved two hardware
threads while execution was pinned to one. Accounting storage is disabled;
[allocation](evidence/job-19184/allocation.txt#L1) and
[completion](evidence/job-19184/completion.txt#L1) receipts are preserved.

Change classes: numerical-method **test candidate**, engineering, documentation
and generated evidence. Governing authority: owner's D1 task and the supplied
review Section D; GOVERNANCE §§2–3 and §5; ADR-0017 §§5–7; INV-11f and INV-15.
No production refactor or pre-baseline exception is involved. AI-authored code
remains an unverified scientific candidate pending domain-owner review.

## Qualification scope and owner disposition

The thermal interpretation remains conditional: sampled-node agreement is not
a proof of a uniform continuous thermal enclosure. The constant-offset runs
are sensitivity probes, not proven extrema over arbitrary time-dependent
thermal perturbations. Neither that assumption nor the C1 threshold is silently
weakened by this receipt.

Main-prefix error from t=0 and thermal RHS roundoff remain open obligations.
This study does not qualify every Stage A checkpoint, complete paired-run
R18/R20 or restart behavior, chemical heating in the one-year temperature
observable, or Stage B. Even a later authorized Stage A run must apply the
specified diagnostic to its actual output checkpoints.

**No production change or Stage A retry is proposed for go-ahead on this FAIL.**
The sole conditional future change named by the review remains O2 chemical
absolute tolerances `(1e-24,1e-24) -> (5e-27,5e-27) MeV`, with every other setting
unchanged, but its D1 prerequisite is unmet. No production change, Stage A,
additional thermal retry, Stage B, master merge or certificate promotion occurred.
