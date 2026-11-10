# Stage B — 2,048-year paired thermal evolution

**Working report: preliminary native experiment; new conditional verification
candidates, pending completion of the original matrix and numerical gates.**
Do not interpret this document as an owner acceptance of Stage B.

Monitoring has stopped at the owner's request. See [STATUS.md](STATUS.md) for
the saved job snapshot and remaining final-report, plot, and manifest work.

The scoped Stage A CONDITIONAL PASS is recorded in
[OWNER_ACCEPTANCE.md](../preliminary-pair/OWNER_ACCEPTANCE.md). Its accepted matched
one-year endpoint is +8,121.86954011 K. Historical failures, original thresholds,
and immutable Stage A evidence remain unchanged. The Stage A report now explicitly
records identical histories across every tier/cap, the thermal RHS systematic
offset, and R20 count quantization. Those one-year conditions are not assumed for B.

## Physical result within the fixed-input mathematical experiment

Both original-t0 roles reached exactly 64,629,964,800 s (2,048 Julian years).
All five independent native references per role completed. The independent
piecewise high-precision reference gives:

| Endpoint quantity | Source | Zero-source control |
|---|---:|---:|
| Redshifted temperature T∞ | 113,875,356.00171844 K | 98,549,552.79750898 K |
| Electron imbalance | −2.1199137019465853e−7 MeV | exactly zero |
| Muon imbalance | −5.311266839549596e−7 MeV | exactly zero |

Matched temperature difference: **+15,325,803.20420946 K**. The finest native
RKF45 reference difference is +15,325,803.20420949 K. These are different numerical
representations, not interchangeable last-digit claims. The new fixed-input
primary-formula high-precision endpoint enclosure for the matched difference is
7.47e−9 K; this is a candidate verification bound, not a complete native Stage B
qualification or an absolute stellar error bar.

The mathematical source's direct deposition dominates the thermal separation.
At the endpoint, direct deposition is 3.3583239603e34 erg/s, photon cooling
3.8002790062e33 erg/s, and full neutrino cooling 5.4049654469e32 erg/s. Source
chemical heating is only 4.1450912065e22 erg/s, while excess neutrino losses are
6.2176368096e22 erg/s. Chemical energy stored is 9.8695948979e33 erg. Control
endpoint photon and neutrino powers are 2.6785500459e33 and 1.7005609542e32 erg/s.
The observed excess-neutrino/heating ratio is approximately 1.5, so the net
reaction term cools despite the positive heating channel. Heating/direct deposition
is only about 1.23e−12 at the endpoint.
The equilibrium-potential deposition is 3.3583239602973852e34 erg/s; the actual
deposition is lower by approximately 3.223e23 erg/s (9.60e−12 of deposition),
consistent with the source contribution to chemical-energy storage. The isolated
LH and net-reaction temperature bounds below do not represent the total effect
of every chemical coupling, including this actual-potential correction.
The muon imbalance is only about 5.4e−5 of kB T∞ at the endpoint, consistent with
the very small thermal contribution from reactions in this benchmark.
Thermal insignificance does not imply absent chemical relaxation. At the endpoint,
the reaction terms oppose about **12.3% of electron and 15.1% of muon source
driving** in the imbalance derivatives. Reaction losses remove about 12.9% of the
instantaneous source input to the chemical-energy reservoir. The imbalances
continue growing in magnitude; they have not reached stationary balance. These
are instantaneous term ratios, not integrated suppression fractions or a new
counterfactual trajectory. Both chemical derivative components are included in
the retained source/control curves.

The retained 1,920-interval trapezoids give the following descriptive integrated
channels. These channel integrals are not assigned the total-ledger error bar.

| Integrated channel [erg] | Source | Control |
|---|---:|---:|
| Direct BNV deposition | 2.17048359e45 | 0 |
| Chemical heating | 7.79464697e32 | 0 |
| Excess neutrino loss | 1.16919705e33 | 0 |
| Full neutrino cooling | 2.22889132e43 | 1.16575820e43 |
| Photon cooling | 2.12429843e44 | 1.76205797e44 |
| Net thermal power | 1.93576484e45 | −1.87863379e44 |

The new analytic coupled comparison in [CHEMICAL_RESPONSE.md](CHEMICAL_RESPONSE.md)
bounds the temperature contribution of LH alone by **+4.936 to +5.615 microkelvin**
(display bounds rounded outward). The net reaction term LH−DeltaLnu contributes
**−2.812 to −2.464 microkelvin**.
These bounds include the response of both chemical variables, thermal damping,
cell changes, and the high-precision/normalization uncertainty. They compare
vector fields analytically; no extra counterfactual trajectory was run. The
relative half-width is about 6.43% for LH and 6.58% for the net term, below the
existing 10% resolved-signal criterion. Both intervals explicitly add an 8.12e−9 K
allowance on each side for normalizing both comparison flows, without cancellation.
The earlier narrower preview remains preserved; this additional allowance changes
no trajectory or threshold.
Thus this contribution is conditionally numerically resolved by the new comparison
verification, but physically negligible in this experiment. It is not resolved
merely by comparing the baseline native trajectory or by subtracting two rounded
temperatures. The omitted finite-temperature source scale alone is about 6.26 K;
accepted stellar inputs are held fixed, with no new absolute-accuracy claim.
The total R20 uncertainty is much larger than the chemical energy terms. The
microkelvin conclusion comes from the separate coupled comparison enclosure;
a passing total-energy balance alone cannot resolve this contribution.

## Actual interpolation and events

The source encounters Cstar knots 94 and 95 at approximately **676.579678 years**
and **1,865.377507 years**, respectively. The control stays in cell 93. The source
ends in cell 95. The lower knot 93 and upper knot 96 are not crossed. These are
actual retained trajectories and high-precision sign brackets, not extrapolated
crossing times. Every source tier has two native events. The high-precision
brackets are 0.01 s wide; native localization width and full trajectory timing
error are reported separately. A zero rounded native residual/bracket is not a
zero physical timing uncertainty.

The cache is continuous at its knots but has derivative jumps. The verification
uses fresh cellwise majorants and a Lipschitz crossing bridge with propagated
prefix error, described in [PIECEWISE.md](PIECEWISE.md). It does not apply the
single-cell Stage A analytic enclosure across knots. The uncapped source native
histories contain 21, 26, and 40 accepted steps for BASELINE, REFINED, and ULTRA;
all have zero rejected steps. The control uncapped counts are 17, 17, and 18.
Therefore Stage A's identical-tier observation does not extend to Stage B.
The BASELINE source quarter cap also changes its history; the observed direct
state comparison improves relative to its three looser caps. Identical subsets
of cap histories are not interpreted as cap convergence.

## Qualification and preserved evidence

Status remains preliminary until every original state, power, reference,
reconstruction, event, restart, and R20 condition is accounted for. New independent
bounds are explicitly candidates under GOVERNANCE §5. Missing terms are not zero,
and observed agreement alone does not waive a contraction or arithmetic test.

The first completed source BASELINE/uncapped configuration has a conditional
piecewise R20 bound of 4.794e−5 N (limit 2e−4), including the original residual,
main/reconstruction/event errors, independent piecewise quadrature bound, original
nested-grid differences, and count quantization. Its maximum R18/G is 1.965e−16.
All 1,921 local coupled reconstruction comparisons fit their stored state and
chemical-energy radii. These initial checks do not by themselves qualify the
other configurations. Final full-matrix evidence will supersede this progress
paragraph without altering the underlying evidence.

Historical BA12 remains **FAIL**, with maximum scaled difference
1.7255917120989046. Historical BA12R c1/c2/c3 remain separately **all FAIL**:
state 32.26041661869904 / 0.9329162771682704 / 43.54669142406804;
ledger 7.611500166686796e−9 / 0.9329163323161163 / 94.83035685919317;
matched stability/contraction/floor 21.923678867554578 / 4.101136780587133 /
41.03343563659818. Their original definitions and limits are retained in
`../preliminary-pair/REPORT.md:203` and
[the historical BA12R report](../../PHASE6A1_CONTROLLED_BNV_BA12R.md).
Neither Stage A acceptance nor a new Stage B result changes those findings.

Large artifacts stay on EKU under:

* `/mnt/sdd/zaki/CompactStar/linux-transition-20261010/stage-b-20261110-r1`
* `/mnt/sdd/zaki/CompactStar/linux-transition-20261010/stage-b-20261110-r2`
* `/mnt/sdd/zaki/CompactStar/linux-transition-20261010/stage-b-20261110-r3`

The first directory retains the original ten references, six uncapped histories,
and all initial launch failures. The second authenticates/copies those retained
histories and completes the first source tiers. The third migrates diagnostic
work interrupted by NFS stalls and provides the unified matrix. Completed source
configuration directories link back to their immutable original r2 evidence. Rehydration restores
accepted records and event metadata; it does not claim bitwise adaptive process
restart or introduce an extra original-t0 evolution.

## Execution, scope, and governance

Requested source: `c90b2205e553afcd102c7acd5f84b253e0a67f1b`, branch
`candidate/results-first-bnv-pilot-20261108`. Execution changes are scheduling,
read-only numerical verification, generated evidence, and documentation.
Production equations, stellar inputs, RK methods, tolerances, interpolation,
breakpoint rules, and scientific guards are unchanged. Governing authority:
GOVERNANCE §§1–5; ADR-0014, ADR-0015, ADR-0016, ADR-0017 and Amendment 1; INV-15.

Independent native configurations use Slurm with **48 hardware threads providing
24 distinct physical reconstruction cores and 48 GiB per task**, unlimited
partition walltime, and single-thread numerical libraries. The initial allocation
of 24 hardware threads provided only 12 physical cores and was correctly refused
by the existing guard after retaining the uncapped histories. Reference jobs
use one core each. Node-local startup delays required explicit placement and a
producer-authenticated native dispatcher. Native workers on phyn001 subsequently
spent most of their time blocked in NFS openat calls (about 6% average CPU, versus
98% on phyn004); phyn003 also stalled during assembly. Their incomplete diagnostic
outputs were preserved before migrating the remaining work to phyn004. The
migration array initially refused a missing RF_MODE before dispatch; its logs
remain intact and the corrected launch supplied RF_MODE=tier. All operational
failures are preserved.
The original one-CPU/30-minute proposal was superseded by the owner's explicit
resource authorization, not treated as a scientific threshold.

This mathematical P2 benchmark does not specify a realistic microscopic BNV rate,
a realistic decay prediction, or an astrophysical constraint. No new stellar-solver
campaign, master merge, historical 500 kyr calculation, or automatic horizon
extension was performed. Stage A acceptance does not accept Stage B in advance.
