# Stage A recovery — stopped on required restart coverage

**FAIL_REQUIRED_RESTART_COVERAGE.** The approved O2 change was implemented and
the one authorized source retry was executed. The actual source accepted history
cannot supply the prescribed positive-length restart tail. The retry was stopped
after this was established; **control was not run and no paired figure is
qualified**. No further tolerance, cap, restart rule or acceptance threshold was
changed. [Stop receipt](evidence/stop-receipt.json#L1).

This is a failure of mandatory coverage for the short Stage A horizon, not an
observed chemical-state accuracy failure. The new tolerance clears the previous
native Echem refusal at 131,490 s. All 129 serialized native checkpoint packets
pass their existing state/diagnostic self-qualification checks. They do not yet
have independent high-precision actual-output qualification, and are not promoted
to such a claim. [Actual-output assessment](evidence/assessment.json#L1).

## Exact obstacle and why execution stopped early

The original plan requires the quarter-cap main's **first accepted endpoint at
or beyond half-time**, with that endpoint strictly before year-end, as the start
of a restart tail. The existing collector enforces
`it != end && it->t_right_s < end` with refusal `no interior restart endpoint`
([collector](../../../../tests/native/stage.cpp#L107),
[original plan](../native-assembly/NONLINEAR_PLAN.md#L65)).

The newly executed BASELINE uncapped source main has 12 accepted steps, no
rejections and no actual Cstar event. Its accepted history and endpoint-only
schedule repeat are byte-identical. The decisive times are:

| Quantity | Seconds |
| --- | ---: |
| Half-year | 15,778,800 |
| Last accepted endpoint before half-year | 10,025,181.919443566 |
| First accepted endpoint at or after half-year | 31,557,600 |
| Final time | 31,557,600 |
| Available restart-tail duration under the rule | **0** |
| Smallest planned cap (quarter-cap) | 9,680,951,939.25 |
| Largest actual next-step proposal | 39,098,210.485829912 |

The last main step spans 21,532,418.080556434 s, from approximately 0.317679 yr
to 1 yr. There are **no accepted endpoints in [half-year, year-end)**. Passive
observation checkpoints inside that step are reconstructed states; they are not
accepted main endpoints eligible under the existing restart rule.

The quarter-cap member was **not executed**. Its inability to supply the restart
is established by an explicit inactive-cap proof using the actual history:
`BreakpointStepper::Advance` uses its cap numerically only in `min(h,max_step_)`
and in its GSL target clamp. For every actual incoming proposal, all three finite
planned caps leave both the proposal and target unchanged. The smallest cap is
about 306.77 years, well above the one-year horizon and every proposal. Starting
from the same state with the same deterministic RHS, induction therefore gives
the same numerical calls, accepted endpoints and next proposals for the
quarter-cap policy. No rejected trials or events are omitted from this proof.
The endpoint-schedule repeat independently checks determinism for the observed
history. This reasoning is **program/path equivalence**, not a claim that an
unexecuted matrix member was measured.

The per-step cap checks, exact selected endpoint, input/source hashes and proof
are in [restart-gate.json](evidence/restart-gate.json#L1), reproduced by
[assess_recovery_restart.py](../../../../tools/linux/native_qualification/assess_recovery_restart.py#L1).
The relevant cap semantics are
[BreakpointStepper.hpp](../../../../CompactStar/Physics/BNV/BreakpointStepper.hpp#L61).
The failure receipt was written **before** `scancel 19188`. Slurm reports
**CANCELLED**, not a timeout and not a native exception from the as-yet-unreached
quarter-cap branch. Continuing the remaining collection could not repair the
missing restart interval under the fixed settings, so the owner's first-failure
instruction was applied. The authorized retry is consumed.

## What completed and what the O2 change accomplished

Job **19188** completed five original-t0 one-year source references, the
BASELINE uncapped one-year source main and its endpoint-schedule repeat. All five
reference histories and observation files are byte-identical to their
authenticated prior counterparts; the new main history also matches the prior
main exactly. This confirms the O2 change did not change those main/reference
trajectories. [Regression/context record](evidence/scientific-context.json#L1).

The collector saved checkpoints 0–128 on the first 240-interval scientific grid:
129 packets, including the initial sentinel, through t=16,830,720 s. No complete
scientific grid was reached, and no restart tail or control trajectory was run.
No partial checkpoint JSON was found. The aggregate uncertainty TSV is a partial
buffered stream; the 129 complete per-point packets are authoritative. All
complete packets and original state
radii/energy allowances remain external and hashed; their compact table is
[checkpoints.tsv](evidence/checkpoints.tsv#L1).

At checkpoint 1 (131,490 s), the previous failed Echem channel now has:

| Quantity | Value |
| --- | ---: |
| Echem | 2.8395321373652751e28 MeV |
| Recorded U_Echem | 2.6370686880579584e16 MeV |
| Original threshold 0.20 F | 5.6790642747305504e16 MeV |
| U_Echem / threshold | **0.46434915339694355 — PASS** |

That is also the largest Echem utilization among the 129 packets. The largest
native state utilization is 0.1032127688845073 (x, checkpoint 76). Every native
state and diagnostic channel in the saved packets passes. The actual
[first packet](evidence/checkpoint-1/checkpoint.json#L1),
[state allowances](evidence/checkpoint-1/state.tsv#L1) and
[diagnostic decomposition](evidence/checkpoint-1/diagnostics.tsv#L1) preserve all
terms, including the propagated state contribution; U<=0.20F is unchanged.

These are the production self-qualification inequalities, **not** the D1 C2/C3
comparisons against an independent reference. The earlier 13-point D1
conditional PASS remains limited to its own intervals. The actual-output
high-precision diagnostic was not reached because restart coverage failed
first; no additional local-reference calculation was launched after the stop.

## Scientific relevance and outstanding accuracy obligations

The restart obstacle does not demonstrate an inaccurate temperature or chemical
state. It prevents the required comparison of a positive-length restarted tail
with the uninterrupted reference. An empty tail at year-end cannot establish
restart compatibility or carry a tested accumulated-error budget. Creating an
eligible endpoint would require a separately authorized change to the restart
test or run settings; none was made here.

The completed source references give the following **diagnostic-only** context
(reference 4, the predeclared finest independent-method RKF45 member):

| Quantity at one year | Retained value |
| --- | ---: |
| Redshifted temperature | 100,007,404.51388617 K |
| Change from source initial temperature | +7,404.51388617 K |
| eta_e | -1.0894196614137575e-10 MeV |
| eta_mu | -2.7634376967724932e-10 MeV |
| Actual direct BNV deposition power | 3.3583239602973686e34 erg/s |
| Reaction heating LH | 5.054307851959247e15 erg/s |
| Photon luminosity | 2.7754493815996762e33 erg/s |
| Neutrino luminosity | 1.9125473233321554e32 erg/s |
| Chemical free energy | 1.6355209485677405e33 MeV |

The original 1920-interval reference-grid trapezoid gives direct deposition
1.059806442094805e42 erg, outgoing neutrino/photon/other energy
9.361242830453821e40 erg, and integrated LH 5.316178661932201e22 erg. These are
descriptive integrals of retained data, **not** a completed R20 or energy-ledger
qualification. No alternative quadrature or trajectory was used. The
temperature increase above is relative to the source's initial state, not a
source-minus-control difference; that paired difference is unavailable.

At year-end LH/direct-deposition power is only 1.5050090199e-19. This short
mathematical P2 experiment primarily tests **direct BNV deposition**; it does
not establish a measurable chemical-heating contribution to temperature.

Maximum raw reference spreads are 3.2526065175e-19 in x,
2.7141863849e-25 MeV in eta_e and 5.9453606527e-25 MeV in eta_mu. Main endpoint
differences from reference 4 are 2.0328790734e-19, 1.8094575900e-25 MeV and
5.1698788285e-26 MeV respectively. These small discrepancies are **not certified
bounds** and do not supply missing uncertainty terms.

| Required assessment | Disposition |
| --- | --- |
| Actual native checkpoint state/Echem inequalities | PASS for 129 saved packets only |
| Independent chemical reference on actual outputs | NOT REACHED after mandatory restart failure |
| Uniform continuous thermal enclosure | Conditional assumption retained, not proved by sampled agreement |
| Continued high-precision refinement contraction | Conditional assumption retained; actual-output test not reached |
| Main-prefix accuracy / accumulated error | UNQUALIFIED; raw reference agreement is insufficient |
| Thermal RHS and accumulated thermal roundoff | UNQUALIFIED; no new bound claimed |
| Global reference convergence / interpolation allowance | UNQUALIFIED pending full actual-output assessment |
| R18 | Static contract passed; dynamic grid packet was buffered and not published before stop |
| R20 with all uncertainties | NOT EVALUATED; original R20 file is empty; buffered output was not published |
| Checkpoint schedule consistency | Actual main and endpoint-only repeat byte-identical |
| Restart consistency | FAIL required positive-length restart coverage |
| Actual breakpoint coverage | No events on completed main; no crossing qualification inferred |
| Matched control / signal resolution | NOT RUN / UNAVAILABLE |

## Implementation, authentication and execution

Starting worktree was clean at requested `4387e56`, branch
`candidate/results-first-bnv-pilot-20261108`, at
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot`.
All 15 local branch tips were ancestors; nine worktrees existed. Original
runtime snapshots and D1 resolution manifests authenticated before edits;
[preflight.json](preflight.json#L1) records this.

The **only production-code change** is the two O2 chemical absolute-tolerance
literals in `PassiveCheckpointOutput.cpp`, from 1e-24 to 5e-27 MeV. O1, O2 rtol,
O2 thermal atol, all main/reference tiers, caps, physical inputs and criteria are
unchanged. The collector gained explicit source/control role selection and
durable actual checkpoint packets before refusal. A compiled configuration
contract and standalone build/publication/assessment scripts were added.
These reporting/execution changes do not change the solver equations.

Build commit `c5d59b1` passed all **seven** contracts (six retained tests plus
the O1/O2 configuration contract). A wrapper typo appended `+` to the static
publisher path and stopped job 19186 **before static verification or thermal
execution**. Commit `a909285` corrected this setup error; the failed log is
preserved. Job 19187 then passed static owner reconstruction, existing ownership/
certificate/currentness/source-control refusals, and byte-identical comparison
of all twelve static RHS states. No stellar regeneration occurred. Replacement
snapshots were committed at `509053d` before source job 19188.

| Identity | SHA-256 |
| --- | --- |
| Previous source snapshot, preserved | c6a29044e8330f079a1593be166238e11eafb66745186c5322c04c272358ec18 |
| Previous control snapshot, preserved | 0588ddf4518b3dec0f1cf103dcb8a3b356cf4631f2394df41b6a2d6baadf4ceb |
| New source snapshot | 7ed42a3ef02d989caa06246940de2cfe75fe7bd554b2650f6314f778fb77c175 |
| New control snapshot | 9226ea1f2ed0542b53b25012fc4cc9d42221efa4093a20174db3b297fd5daf5a |
| Unchanged full-owner archive | 26dd2b060ad8ed10d26ec46b67a9b957e3517fe6718708ba91d09c5f204158cf |

The publisher checked the production source against 4387e56, permitting exactly
the two approved literals. Its explicit snapshot diff verifies unchanged
physical/common configuration and identical runtime linkage, while binding the
new binaries, source and O1/O2 configuration. Both roles share the same common
snapshot; only the original source object differs. Both new snapshots rechecked
successfully after the stop. Original snapshots, binaries, failures, D1 FAIL
and conditional PASS receipts remain preserved.

Source job 19188 ran on phyn004 for **6 minutes 3 seconds** before deliberate
cancellation, within its 30-minute ceiling, 4 GiB limit and one logical execution
CPU (affinity [0]). Slurm reserved two hardware threads; execution was pinned to
one. The canceled measuring process did not publish a final peak-RSS receipt,
so none is claimed. Build/static jobs also used one logical CPU and 4 GiB.
See [source allocation](evidence/source-allocation.txt#L1),
[completion](evidence/source-completion.txt#L1),
[source log](evidence/source-stage.log#L1),
[contracts](evidence/contracts.log#L1) and
[snapshot diff](evidence/snapshot-diff.json#L1).

Large artifacts remain under
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/stage-a-recovery-20261109-r1`.
Git contains maintained code, source/control snapshots, compact actual-output
receipts, the cap/restart proof, this report and SHA-256 manifests. Analysis is
reproducible without another solver call; see [REPRODUCE.md](REPRODUCE.md#L1).

Change classes: numerical-method (explicitly approved tolerances), engineering,
dependency/build, documentation and generated evidence. Authority: owner's
Stage A recovery instruction, GOVERNANCE §§2–3/5, ADR-0015/0016/0017, INV-11f/
INV-15 and the original Stage A plan. This remains a scientific candidate;
there is no certificate promotion or claim of human ratification of new analysis.

No control, Stage B, 14,648-year benchmark, 500 kyr campaign, stellar
regeneration, BA13 reopening, master merge or additional thermal retry occurred.
The task stops with the exact restart obstacle preserved. No paired thermal
result is presented as trustworthy while the required qualification remains
incomplete.
