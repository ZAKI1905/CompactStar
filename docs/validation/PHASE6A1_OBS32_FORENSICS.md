# Observation-32 reconstruction forensics

## Pre-execution declaration (2026-09-30)

Owner-authorized bounded local investigation from
`7d2777561a17b56b466fb024ef13104da6c0e86a` on
`analysis/phase6a1-obs32-reconstruction-forensics`, worktree
`/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase6a1-obs32-reconstruction-forensics`.
Canonical master remains `9a5c1eca40758e4c39f9f80db7324e12ad23e505`.
The source branch, local upstream and live remote authenticate at entry. Historical
BA12, BA12R and the observation-32 campaign STOP are immutable failures.

Classes: diagnostic numerical-method experiments, engineering instrumentation,
documentation and generated evidence. Governing authority is the owner's bounded
request and accepted ADR-0017 sections 5.2–5.4; no production method change is
made or inferred. Phase-5/Cstar source, cache values, physics, source card, grid and
frozen domain are unchanged. No full trajectory or campaign control is executed.

The executable first calls the unmodified production reconstructor and requires
exact equality of all six archived O1/O2 endpoint values and the unresolved
status. A discrepancy halts all subsequent experiments. Instrumented copies of
those two local solves must also reproduce the same bits. GSL step-apply wrappers
record all attempted/rejected/accepted steps without changing controller logic;
RHS logging reads the actual cache cell selected by production. A test-only
member-pointer accessor reads the private cache without changing production
headers/layout, values, or lookup code.

The preserved target and three smooth controls are in
`evidence/obs32-forensics/cases.tsv`. Controls are selected before new solves:
three nearest distinct knot-free accepted brackets whose widths are 0.3–3 times
the target width; select an existing scientific observation nearest the target's
fractional position in each bracket. This yields observations 26, 36 and 18.

A: unsplit rk8pd O1/O2/O3, plus one O4 attempt; retain any GSL pathology rather
than forcing O4. B: find the knot event using a bounded bisection on local O3 IVP
solutions from the exact saved left state, at most 64 event evaluations; no
interpolation of archived endpoint solutions. Use the same event time for O1/O2/O3
splits, preserve each level's state continuously without projection, and restart
numerical state at the split. Record event residuals and bracket uncertainty.
C: independent GSL RKF45 O2/O3, unsplit and split at that same event.
D: the three preselected smooth controls at rk8pd O1/O2/O3.

Full currentness checks bracket each ordinary diagnostic solve. The root-finding
experiment has full checks before/after the bounded root sequence; every trial
has fresh state/GSL objects and retains cheap currentness on every RHS evaluation.
No production currentness rule is changed. The diagnostic root search never feeds
an authoritative trajectory. Any successful split remains a diagnostic prototype.

Inspect actual Cstar values at all 160 cache points, compare all 32 archived
qualified prefix Cstar values exactly, and rebuild a cold thermal cache to test
payload identity. Audit both endpoint-ownership histories at the knot and complete
RHS values around it. Compare direct main-style RHS evaluation against disposable
replay contexts at identical prescribed states. Error-scale arithmetic is checked
against ADR-0017 without relaxing the failed criterion.

Current source finding: `StarContext.cpp:795-802` interpolates Cstar linearly in
**log(T)**. `Bracket` at lines 41–53 retains an inclusive cached interval, otherwise
uses `upper_bound`. Neither fact by itself establishes the failure's cause.

The original cache was runtime-only; no complete standalone cache payload was
archived with the failure. Authentication therefore binds original input/source
bytes, exact reproduced O1/O2 and archived Cstar samples, plus a newly exported
and independently cold-rebuilt cache payload. This limitation is explicit.

Results follow below after execution. No candidate, campaign restart, physical
BNV model/rate, A18, superfluidity, Regime-II/MixedStar, or sliding background is
authorized. Any remedy changing ADR-0017's uniform rule requires owner acceptance.
