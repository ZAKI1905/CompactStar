# Bounded Echem reconstruction feasibility study

Owner-authorized numerical-method investigation, not qualification or a policy
change. Starting source: clean `56dd1959a945bc13de574e83ba664998ae87e3b7`,
branch `candidate/results-first-bnv-pilot-20261108`, at
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot`.
All 15 local branch tips are ancestors; nine worktrees are present. Authority:
owner's Echem recovery request; GOVERNANCE sections 2–3; ADR-0017 sections 5–7
and its breakpoint amendment; INV-11f and INV-15. The owner's explicit request
authorizes the bounded tolerance study and supersedes the original prohibition
on tolerance changes for this investigation only. Acceptance allocations stay fixed.

## Before numerical execution

Reuse job 19174's immutable BASELINE source accepted history and five reference
observation tables. No new main, reference-from-t0, control, stellar solve or
Stage B trajectory is part of this study. Fresh isolated local reconstructions
start at the retained main's left accepted state, as in the original algorithm.
Authenticate the original source/control snapshots and retained files with their
existing manifests; preserve all bytes. Link the original numerical archives.

Test finest-grid indices **1,2,4,8,16,32,64,128,256,512,960,1440,1919** on the
1920-interval one-year grid. This includes the first point on every nested
scientific grid, the original failure, and later intervals with larger states.
Uniform O2 chemical absolute tolerances, in MeV:
**1e-24 (baseline), 1e-26, 5e-27, 1e-27, 1e-28, 1e-29**.
Keep O2 rtol=1e-13, x atol=1e-18, O1, built-in independent witnesses, algorithms,
eight-corner propagation and all acceptance criteria unchanged. A test-only
configuration override on a nonconst reconstructor changes its stored O2 setting;
no authenticated production source or library is edited. Serialize actual settings.

At each point independently integrate the same local IVP with rk8pd and RKF45,
rtol=1e-15, atols=(1e-20,1e-30,1e-30), maximum steps span/4, /8, /16.
Record actual steps, rejection counts, local errors, RHS values, terminal states,
Echem and its corner values. These witnesses are diagnostic only, never a fallback
reported state. No tolerances or caps may be added in response to results.

## Prospective acceptance and stopping rules

The baseline must reproduce the original checkpoint-8/1920 failure. All runs
must reach the target, remain finite and within original scientific guards, and
have at most 10% rejected steps. All original state and diagnostic gates must
pass for a candidate. Require **Echem U/(0.20F) <=0.80** at every study point,
reserving 20% margin; select the **largest** listed chemical atol that passes
all requirements, without changing any other setting. Serialize every candidate,
including failures. The study is one bounded batch, not an adaptive search.

Independent observed discrepancies plus twice the maximum finest method or
adjacent refinement discrepancy must fit inside the candidate's unchanged state
and Echem uncertainty envelopes. Refinement must contract by <=0.5 when above
an independently justified roundoff bound. A 64-ulp floor or agreement alone
does not bound accumulated roundoff. Check recorded Echem algebra at 100-digit
precision and audit chemical RHS/step arithmetic; use an analytic bound or
higher-precision evidence only where its assumptions are explicitly established.
Missing accumulated-roundoff or reference uncertainty evidence is **UNQUALIFIED**,
not zero. Smaller configured floors alone are not a feasibility PASS.

If these checks pass, publish the candidate setting and verification receipt
before implementing the minimal adjustment. One conditional Stage A retry may
then use the original matrix/criteria in `../native-regeneration/RUN_PLAN.md`,
with source fully assessed before control. All original R18/R20, restart, reference,
roundoff and matched-signal requirements remain mandatory. Source refusal prevents
control execution. No second Stage A retry or longer campaign is authorized.

Study allocation: phyn004, **one logical execution CPU, 4 GiB, 30 minutes**;
build separately with the same ceiling. Fresh external directory under
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/echem-refinement-20261109-r1`.
All substantial artifacts stay there. No Mac-side research artifacts are needed.
Timeout or unexpected refusal is INCOMPLETE/FAIL, never automatic resubmission.

If no candidate meets the evidence requirements, stop without Stage A and report
the measured precision and the missing or exceeded bound. A proposed change to
energy accuracy policy must remain unimplemented pending explicit owner approval.
