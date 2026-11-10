# Stage A source retry — INCOMPLETE at wall-time limit

**Source is not qualified.** The exact owner-authorized command from preparation
commit `37a9a21a85dfb03d9913d8e39739b294bb1acdd5` ran once, as Slurm job **19191**.
Slurm stopped it with `TIMEOUT / TimeLimit`, exit `0:15`, while collecting the
BASELINE uncapped 960-interval checkpoint grid. There was no native checkpoint
refusal or observed necessary-accuracy-test failure before termination. The
obstacle is incomplete qualification coverage within the fixed runtime budget;
this attempt does not demonstrate numerical inaccuracy or complete accuracy.
No resubmission, tolerance change, control trajectory or Stage B followed.

## Execution and authentication

The clean starting worktree was
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot`,
branch `candidate/results-first-bnv-pilot-20261108`. All 15 local branch tips were
ancestors; nine worktrees were inventoried. Preparation manifests authenticated
31 compact and 1,349 external files, and the output path was absent before the
single submission ([preflight](evidence/preflight.json)). The exact executed
command is retained in [REPRODUCE.md](REPRODUCE.md).

The replacement source snapshot remains
`4fa012bda09d5022df7dce6ad7278978b7f638c3d15786d84750a96da256c4d5`;
the authenticated, **unexecuted** matched control snapshot remains
`584251e18fb131cb8f141e7972de6fcfaecf52a53f818d5f5f6a3f3c8d19fb27`.
Final snapshot checks pass. Production sources, harness, solver settings,
physical inputs, acceptance thresholds and old evidence are unchanged.
The approved latest-noninitial-accepted-endpoint-at-or-before-half-time selector
is bound into the executed binary, but this attempt did not reach a restart tail.

The allocation was phyn004, one task, one execution CPU (affinity `[0]`), 4 GiB,
and a 30-minute limit. Slurm reserved two hardware threads while execution was
pinned to one. Start: **2026-11-09 22:43:36 EST**; end:
**2026-11-09 23:13:42 EST**. Slurm records **00:30:06** elapsed: its timeout
termination occurred six seconds after the configured deadline; no extension
was requested. The termination signal prevented the wrapper's final resource
receipt. Exact process CPU time and peak RSS are unavailable: Slurm accounting
is disabled, and `sstat` returned an invalid CPU value and blank RSS. These are
not replaced with estimates ([completion](evidence/completion.txt),
[resource counters](evidence/resource-counter-query.txt)).

## Coverage and gate disposition

| Required evidence or gate | Actual result |
| --- | --- |
| Five original-t0, one-year independent reference histories | 5/5 collected; byte-identical to authenticated prior outputs; uncertainty/convergence bounds remain unqualified |
| Main trajectories / endpoint-only schedule repeats | 1/12 each collected; BASELINE uncapped, 12 accepted steps; repeat byte-identical |
| Four scientific grids per main | 2/48 complete: 240 and 480 intervals; 960 grid has points 0–356 of 0–960; no 1920 grid |
| Actual native checkpoint packets | 1,079 serialized, all native state and diagnostic checks pass |
| Echem local native gate | Maximum utilization 0.5792747394 of the unchanged `U <= 0.20 F` allowance |
| Actual-output high-precision chemical reference | NOT REACHED before runtime stop; no checkpoint certified from the prior 13-point study |
| Main-prefix accuracy / accumulated roundoff | UNQUALIFIED; raw agreement and native local allowances are insufficient |
| Thermal RHS roundoff / continuous thermal enclosure | UNQUALIFIED |
| R18 | All 722 published samples pass; max normalized residual 1.63307347e-16 versus 1e-10 |
| R20 | Two original grid residuals retained; necessary partial tests pass; total uncertainty and quadrature contraction UNQUALIFIED |
| Scientific domain | Published samples pass; maximum frozen utilization 3.47484089e-6, maximum absolute fractional depletion 1.00002518e-12 |
| Actual event crossing coverage | None in completed main; remains UNQUALIFIED |
| Quarter-cap / restart tails | 0/3 tails reached; restart numerical consistency NOT TESTED |
| Control / paired signal | NOT RUN / NOT AVAILABLE |

The final complete packet is 960-grid point **356**, at **11,702,610 s**.
Collection stopped before point 357 was serialized. The complete 240/480 grids
reach the one-year endpoint, but the source matrix has not completed. The
collector's fixed order processes all four grids before advancing to the next
cap, and only tests a restart after the quarter-cap grids (`tests/native/stage.cpp`).
No suitable endpoint was rejected in this run; the selector was not reached.
No passive output packet was substituted as an accepted restart state.

At the worst native Echem packet (960-grid point 1, **32,872.5 s**), recorded
Echem is **1.7747077545683595e27 MeV**, U is **2.05608674394112e15 MeV**, and
`0.20 F` is **3.549415509136719e15 MeV**. Native acceptance is conditional on
its state/error premises, which still require the on-output independent check.
Maximum native state utilization is **0.1054382381** (thermal x, 240-grid point
239). Compared to the finest retained reference, the largest necessary partial
state-budget ratios are about **2.16e-6 (x), 1.92e-6 (eta_e), 2.07e-6 (eta_mu)**.
These omit unknown reference/prefix terms only to test necessary conditions;
the missing terms are never taken as zero for qualification.

## Energy accounting: retained residual and a bounded arithmetic component

For the 240 grid, original R20 residual is **-2.6681674277032373e37 erg**,
normalizer **9.6619401379124758e41 erg**, and recorded reconstruction allowance
**3.0333357575090393e32 erg**. The residual plus that allowance divided by N is
**2.7615548461e-5**, below the **2e-4** total criterion. The 480-grid residual is
**-2.6681674273570009e37 erg**, with allowance **5.2888729365761848e32 erg**.
Neither is a qualified total ledger result: main, event, thermal and quadrature
uncertainties and contraction assessment are incomplete.

One arithmetic contribution can be isolated without changing a trajectory.
The retained source evaluates `B(t) = fl(B0 + fl(Bdot*t))`; R20 evaluates
`fl(fl(MeVToErg*mu_B) * (B_end-B0))`. Subtraction of these nearby represented
endpoint counts is exact (verified explicitly), but their rounding loses
precision in the small count change. Source expressions/constants and hashes
are bound in the maintained assessor and receipt. Exact lifts of binary64
inputs, with 200-digit Decimal arithmetic, give a signed equilibrium-energy
arithmetic error **-2.6681675260272892e37 erg** relative to the exact linear count
change. It accounts for essentially all of the displayed raw residual.

A nearest-rounding bound for this component is **6.060344943817435e37 erg**.
Writing `c = exact(MeVToErg)*exact(mu_B)` and `s = fl(MeVToErg*mu_B)`, it is

```
U_B = (ulp(B_end) + ulp(fl(Bdot*t))) / 2
U_eq = |c| U_B + |B_end-B0| ulp(s)/2 + ulp(fl(s*(B_end-B0)))/2
```

This assumes the retained binary64 nearest rounding and `-ffp-contract=off`,
normal finite arithmetic, exact input doubles, and the verified constant Bdot
and mu_B. The actual exact error is checked against the bound. The original
residual is **not corrected or discarded**. Including this bound and the
recorded reconstruction allowance gives **9.03394409e-5 (240)** and
**9.03396744e-5 (480)** of N. These still leave room beneath 2e-4, but bound only
the stated component. They do not bound thermal RHS arithmetic, accumulated
trajectory/reference error, thermal quadrature, or the entire ledger. Small
raw differences do not prove those bounds.

## Descriptive scientific output — not a qualified paired result

The finest retained one-year reference gives Tinf **100,007,404.51388617 K**,
starting at **100,000,000 K**: a source-only increase of **7,404.51388617 K**.
Endpoint eta_e is **-1.0894196614137575e-10 MeV**, eta_mu
**-2.7634376967724932e-10 MeV**, and Echem **1.6355209485677405e33 MeV**.
These outputs carry no completed total uncertainty qualification.

Endpoint direct deposition is **3.3583239602973686e34 erg/s**, net thermal
power **3.0616535489040794e34 erg/s**, photon luminosity
**2.7754493815996762e33 erg/s**, and neutrino luminosity
**1.9125473233321554e32 erg/s**. The retained 1920-interval reference trapezoid
for direct deposition is **1.059806442094805e42 erg**; outgoing radiation totals
**9.361242830453821e40 erg**. Main-grid thermal energy change is about
**9.66194013791248e41 erg**, chemical free-energy change **2.6203934482e27 erg**.

Endpoint chemical heating is **5.054307851959247e15 erg/s**, about
**1.50500902e-19** of direct deposition; its reference-grid integral is
**5.316178661932201e22 erg**. The short experiment primarily tests direct BNV
deposition. It does not establish a measurable chemical-heating effect or a
source-minus-control temperature signal. No paired figure is published.

## Evidence, changes and scope

The maintained read-only `tools/linux/native_qualification/assess_source_retry.py`
produces the complete external assessment, compact [summary](evidence/summary.json)
and compressed [packet inventory](evidence/checkpoint-inventory.tsv.gz).
It verifies old trajectory hashes before comparison, uses exact binary64 lifts,
checks native failures, necessary state/energy conditions and domain guards,
and reports every missing obligation. All 12 compared reference/main files are
byte-identical to the authenticated prior records. Postprocessing asserts the
source count expression, exact endpoint subtraction, and component bound for
both completed grids. Snapshot reauthentication and source-diff checks pass. Final manifest checks
authenticate all 31 compact and 1,349 external preparation entries and 4,360
new external files. Authored code/prose whitespace checks pass; original Slurm
allocation/completion trailing whitespace is retained byte-for-byte.

All raw outputs, logs and partial files remain under
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/stage-a-restart-retry-20261109-r1/source-stage`.
They are covered by `evidence/external-SHA256SUMS`; compact records and assessor
sources are covered by `SHA256SUMS`. Previous D1 FAIL/PASS and Stage A failure
records remain unchanged. New files are only the assessor and this evidence
package; no production or harness numerical result is changed by this commit.
Classification: analysis/engineering, documentation, generated evidence.
Authority: explicit owner instruction, GOVERNANCE §§2–3 and 5, ADR-0017 §§4–6
and Amendment 1, approved restart amendment and the unchanged fixed run plan.
The analysis remains an unratified scientific candidate.

Conditional uniform thermal enclosure is still an assumption; agreement at
sampled nodes is not a continuous enclosure. Chemical reference error also
requires the declared refinement assumption. Neither assumption is discharged
by this incomplete run or by the older 13 intervals. Accumulated main-prefix
error, thermal roundoff, reference convergence, full R18/R20 coverage and actual
restart accuracy remain open. This attempt consumes the authorized source retry
and reaches the stage runtime ceiling. Any further numerical execution requires
an explicit decision; control remains withheld. No automatic continuation,
stellar solve, BA13 reopening, master merge, certificate promotion or longer
campaign occurred.
