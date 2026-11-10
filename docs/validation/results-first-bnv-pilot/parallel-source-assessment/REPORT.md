# Stage A source collection complete; full qualification UNQUALIFIED

Job **19195 finished the complete source matrix with exit code 0** in
**7,953.050998695195 seconds (2h12m33s)**. All 43,248 native checkpoint packets
pass. There is no observed failure of the necessary accuracy tests performed
here. The source is nevertheless **not fully qualified**: required independent
actual-output chemical, thermal and accumulated-error bounds are missing.
No control, retry or Stage B was launched. This is a scientific candidate,
not a certificate (`evidence/summary.json:1`, `evidence/completed-tests.json:1`).

## Authority, authentication and scope

The clean entry was d21a9d4 on candidate/results-first-bnv-pilot-20261108 in
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot`.
All local branch tips remain ancestors; nine worktrees are recorded
(`evidence/preflight.json:1`). The owner requested status and continued analysis,
and conditionally approved public publication after sensitive-data review.
This change is engineering/read-only analysis, documentation and generated
evidence under GOVERNANCE §2; no production source, tolerance, physical input,
integration method or acceptance criterion changes. ADR-0017 §§4–6 and Amendment
1 govern the accepted histories and isolated reconstructions. The unchanged
numerical gates are `../native-assembly/NONLINEAR_PLAN.md:95-134`; missing terms
stop qualification. The earlier D1 review and knot-aligned resolution remain
local to their 13 intervals and are not promoted to universal coverage.

The replacement pair passed a fresh snapshot check. Source identity remains
`cbe2cd7f4042b7586149b28428d6ea86bfe245997a520c5ca82829734b5b8a9f`.
The unexecuted control remains
`880dee60847576fe829dc0c8fd854862c664768a101453d7430f3c1f4bf126e6`.
Full original source bytes are inventoried in the compressed external manifest;
compact evidence has SHA256SUMS. Original D1 FAIL, previous restart failures,
serial TIMEOUT and preparation failures remain intact. Files added here are the
read-only `tools/linux/native_qualification/assess_completed_source.py`, this
receipt, reproduction instructions and evidence; no numerical solver was edited.

## Completed coverage and measured accuracy indicators

| Item | Result |
| --- | --- |
| Original-t0 reference histories | 5/5, each with 1,921 reference observations |
| Main trajectories / endpoint-only repeats | 12/12 each; every main has 12 accepted steps; every repeat is byte-identical |
| Scientific grids | 48/48 (240, 480, 960, 1920 for each main) |
| Native checkpoints | 43,248/43,248 pass; no native refusal |
| Maximum native state utilization | 0.1058747896, thermal x at 1920-grid point 1918 |
| Maximum native Echem utilization | 0.7330768176 at 1920-grid point 1, t=16,436.25 s |
| R18 residual / gross power | Maximum 1.7040699035e-16, below 1e-10 |
| Scientific domain | All 43,248 published samples pass; frozen utilization <=3.4748408919e-6; depletion <=1.0000251760e-12 |
| Quarter-cap restart tails | 3/3 collected from ordinal 11 at t=10,025,181.919443566 s; exact binary64 selection/round-trip verified |
| Actual event crossing | None; crossing coverage remains UNQUALIFIED |

These counts and extrema are in `evidence/summary.json:1`. At the worst native
Echem point, Echem=4.4367694567187725e26 MeV,
U=6.50498566782976e14 MeV and the original 0.20 F allowance is
8.873538913437545e14 MeV. Native acceptance still depends on its error-box
premise; it is not an independent verification of that premise.

Each restart retains a 21,532,418.080556434-second tail. At the final endpoint,
restart minus uninterrupted state is 0 in x, -1.2924697071141057e-26 MeV in
eta_e and +5.169878828456423e-26 MeV in eta_mu (one ulp in each chemical component).
Against the finest independent RKF45 reference, restart discrepancies are
2.032879073410321e-19 in x (15 ulp), 1.6802106192483375e-25 MeV in eta_e
(13 ulp), and zero in eta_mu. These are raw comparisons, not prefix-inclusive
uncertainty bounds (`evidence/completed-tests.json:1`). The old coverage defect
is resolved; full numerical restart qualification remains conditional on the
missing reference/prefix bounds.

## Reference convergence and remaining numerical obstacle

The maximum necessary state ratios, including recorded local output radii but
excluding unknown positive reference terms, are 0.021575 (x), 0.0191995 (eta_e),
and 0.0207461 (eta_mu), all in ULTRA. Reference allocations also pass their
necessary lower-bound tests: twice the largest required difference plus the
64-ulp minimum uses at most 3.656%, 2.779% and 3.547% of the respective 5%-of-ULTRA
allocations. This is budget headroom, not proof that accumulated roundoff fits
inside it (`evidence/summary.json:1`, `evidence/completed-tests.json:1`).

Global reference refinement does not show clean contraction at every sample.
For example, at t=31,130,257.5 s, successive rk8pd x differences are
1.3552527156068805e-20 and 4.065758146820642e-20 (ratio 3), while the required
64-ulp minimum is 8.673617379884035e-19. At this point the necessary reference
allowance is 1.517883041479706e-18 against a budget of 4.151991183819577e-17.
The differences are a few binary64 ulps; a ratio above 0.5 alone does not prove
physical inaccuracy. Equally, the minimum 64-ulp floor is not an independent
accumulated-roundoff bound and cannot automatically excuse noncontraction.
The x/e/mu noncontracting counts are 1570/1248/1605 of 1921 samples. Zero
predecessor differences are reported explicitly without undefined division
(`evidence/completed-tests.json:1`).

The decisive evidence gap is concrete: **zero `witness-steps.tsv`, `local-rhs.tsv`
or `witnesses.tsv` files are present for the actual output**. The collector saves
O1/O2 states, allowances and aggregate diagnostic packets, but not accepted
witness nodes or their derivatives (`tests/diagnostics/checkpoint_evidence.hpp:43`,
`tests/native/checkpoint_transport.hpp:53`). Consequently, the approved
knot-aligned thermal path and two-ulp RHS-fidelity check cannot be reconstructed
from these packets alone. The actual-output high-precision diagnostic was
**not executed**; neither its thermal-band premise nor chemical-state box can be
claimed verified. The prior 13-point PASS is not substituted.

Independently supported thermal RHS evaluation error, accumulated main/reference
roundoff, reference interpolation at accepted endpoints, and complete R20
U_main/U_event/rounding terms are also not supplied. These are mandatory
qualification obligations, not demonstrated excess physical error. No tolerance
adjustment or additional trajectory has been made to hide them. The present
qualification stops under the existing missing-uncertainty rule.

A useful next evidence task would capture bounded local witness/RHS traces
against these immutable actual accepted histories, then apply the unchanged
knot-aligned diagnostic to the actual outputs. It must also establish the
independent thermal/main-prefix bounds needed to excuse arithmetic-floor
noncontraction. That requires no new main trajectory merely to recover coverage;
it is not authorized here as a new solver retry or an amended acceptance rule.
Uniform continuous thermal enclosure and continuation of refinement contraction
must remain explicit assumptions. Constant-offset thermal runs alone are not
extrema over arbitrary time-dependent perturbations.

## Energy and temperature: descriptive, not yet qualified

The finest reference ends at **Tinf=100,007,404.51388617 K**, a rise of
**7,404.51388617 K from the source initial state**. This is not a source-minus-control
temperature difference. Final chemical imbalances are
eta_e=-1.0894196614137575e-10 MeV and eta_mu=-2.7634376967724932e-10 MeV
(`evidence/summary.json:1`).

Original 1920-grid reference trapezoids give direct deposited energy
1.059806442094805e42 erg, photon losses 8.757867510379589e40 erg,
neutrino losses 6.033753200742326e39 erg and net thermal input
9.661940137902668e41 erg. Recorded main thermal-energy increase is
9.661940137912547e41 erg and chemical-energy increase is
2.6203934482127427e27 erg. Chemical heating integrates to only
5.316178661932201e22 erg; endpoint chemical heating/direct power is
1.5050090199e-19. This short experiment primarily tests direct BNV deposition
and establishes no measurable chemical-heating contribution.

The original R20 residual is -2.6681674272428783e37 erg at N=9.661940137912547e41
erg. Its magnitude plus recorded reconstruction uncertainty is
2.7617440174e-5 N, below the 2e-4 ceiling. A separately supported equilibrium-count
arithmetic component has bound 6.060344943817435e37 erg (6.2723892482e-5 N).
Adding that bound leaves a known subtotal about 9.0341e-5 N. The residual is
retained, not corrected away; unknown positive terms are not set to zero.
The final nested-grid contraction ratios are 0.297619 for the residual,
0.142857 for thermal-energy change and 0.321429 for outgoing energy; Echem
change is identical across grids. Full quadrature/thermal error qualification
still requires its independent rounding evidence
(`evidence/summary.json:1`, `evidence/completed-tests.json:1`).

## Execution, validation and publication

Source allocation: 24 physical worker cores, 96 GiB, five-hour ceiling.
Measured user CPU is 134,510.877333 s and system CPU 45,452.539897 s;
maximum child RSS is 436,492 KiB (not simultaneous aggregate RSS).
Slurm accounting is disabled and the completed job record has expired, so
exact allocation end time and aggregate peak memory are unavailable. The
successful wrapper resource receipt supplies measured process wall time and
exit status (`evidence/resource-usage.json:1`). The inherited summary's generic
resource-counter caveat refers to the earlier serial attempt; this paragraph
and the actual resource receipt govern this completed run.

Read-only EKU jobs 19196/19197/19198 assessed existing files; none integrated a
trajectory. The first supplemental report reversed the two restart comparison
row labels. The corrected assessor checks them against the independently saved
main/tail endpoints, and job 19198 passes these assertions. Original report and
script are preserved externally; `evidence/ASSESSOR_ERRATUM.md:1` records the
correction. No scientific criterion changed. Full assessment and compressed
per-point inventory remain under /mnt/sdd and are externally hashed.

Publication review scanned all 59 changed historical blobs since the public
branch tip 88da801, and inspected the new code/logs/results. No credentials,
private keys, personal contact records or authentication configuration contents
were identified. Paths, author username, node labels and job metadata are
reproducibility metadata; their presence was explicitly raised with the owner.
The owner approved publication after this review. Review limitations and exact
scope are retained in `evidence/publication-review.json:1`; it is not a claim
that pattern scanning can guarantee the absence of every possible secret.
