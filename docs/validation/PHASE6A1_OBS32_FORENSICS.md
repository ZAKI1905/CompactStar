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

### Supplemental diagnostic declaration (before execution)

After exact reproduction and the initial split comparisons, two targeted checks
will complete the requested boundary/main-crossing audit. In a second freshly
assembled context, evaluate the complete RHS at the nearest representable packed
states below/at/above the exact knot, warmed from both directions; report if no
packed state represents the knot exactly. Extend each already-computed O3 split
left leg from its own event state to the preserved `t_R`, once with rk8pd and once
with independent RKF45. These are **two local legs inside the same saved bracket**,
not main trajectory reruns. They quantify the saved BASELINE main step's local
endpoint discrepancy. No original result is replaced, and no production fix is
made. Primary O3 event states are authenticated inputs to these extensions.

## Completed bounded investigation

**Disposition: A — KNOT-CROSSING ORDER-REDUCTION / PIECEWISE-SMOOTH NUMERICAL EFFECT.**
This is a numerical inference supported by multiple controlled comparisons, not
an empirical measurement of a formal RK order exponent. The archived failure
reproduces exactly. No production implementation defect was found. The existing
fail-closed gate behaved correctly. The successful event-split diagnostic changes
ADR-0017 semantics and is **not** a governed replacement or permission to resume.

The authenticated lookup is linear in **log(T)**, not linear in T. Its value is
continuous and its slope changes at the knot. Tight unsplit rk8pd approaches the
split value; splitting removes the O1/O2 discrepancy without changing the RHS;
three preselected smooth controls pass. Context, cache payload, actual cell
selection and scale audits do not explain the failure. The independent RKF45
method supports the split result but shows a nonmonotonic unsplit tolerance trend.

### Authentication, scope and reproducibility

Entry: `7d2777561a17b56b466fb024ef13104da6c0e86a`, evidence branch
`physics/phase6a1-controlled-bnv-resume`. New work is confined to diagnostics,
manual test targets, evidence and documentation. All CompactStar production bytes,
Phase-5, original archives/raw files, protected artifacts and baselines were
rechecked unchanged; counts are in `evidence/obs32-forensics/exit-verification.json`.
The authenticated dependency heads are Zaki
`7aca56cfdf3875cf2a8c0acea4acdac44f69c9f5` and CONFIND
`b0cbd510fd3fd0c772fa50499cd749287cb39e7b`. The same qualified local Debug packages,
AppleClang 21.0.0, arm64 and GSL 2.7.1 are recorded with artifact hashes in
`packages.json`. No dependency modernization or qualification was performed.

`authentication.json` binds every saved campaign artifact, the original manifest,
source payload and input trees. Original failure manifest SHA256:
`79c1ee93b767f4a55313499a5a6a56678381acf66c879b9d3e41523f43d6a855`.
The original executable hash is
`a0a55b3280ade8016e8d0dc502f12b56acfc9fab1f3c788653b0c07414388fea`.
The new 71-file deterministic-gzip run manifest is
`evidence/obs32-forensics/runs/manifest.json`, SHA256
`71cb3e84039d683a4da0be3f643c9b8530bcdb199fbe79a73ff0ff10dac6f271`.
It binds every trial/RHS trace, root evaluation, result, console log and execution
receipt. The receipts include commands, executable hashes and source commits.
The absent original full cache dump remains the authentication limitation stated
above; all 32 original Cstar samples and all six O1/O2 values reproduce exactly,
and the new independent cold cache rebuild matches all 160 entries bit for bit.

The unchanged source card is `CPL-P2-LINEAR-QSS-v1`: prescribed neutron sink,
`Bdot/B0=-1e-12 yr^-1`, `Bdot=-2.4136520263641375e37 counts/s`, P2 terminal fate,
`5e5 yr`, 8193 uniform passive observation times, initial packed state `(0,0,0)`
with `Tinf=1e8 K`. Structure-1 uses `rho_c=1.10e15 g/cm^3`, radial 80000, EOS8192,
Me/Mmu on, De/Dmu off, StaticZeroSpinHistory/AnalyticControl. No rate/model was
selected. These are authenticated inputs, not a new campaign. No full trajectory,
matched campaign control, REFINED or ULTRA was executed. Maximum prescribed
campaign depletion remains the archived `4.9999999996026683e-7`, below `1e-6`.

Raw ignored output roots: `build/obs32/run-1` and `build/obs32/boundary-1`.
The committed compressed data can reproduce the entire reduction without an ODE:

```sh
python3 tests/bnv/analyze_obs32_forensics.py \
  --root docs/validation/evidence/obs32-forensics/runs/primary/output \
  --evidence docs/validation/evidence/obs32-forensics \
  --supplement docs/validation/evidence/obs32-forensics/runs/supplement/output \
  --output /tmp/obs32-analysis.json
```

This reduction was run on both raw and compressed evidence; outputs are identical.
It independently checks every archived budget and endpoint exactly, all logged
trial-cell containment, trace step counts, smooth-bracket knot absence and
warm-history RHS equality. `python3 tests/bnv/check_obs32_reduction.py` also
rejects a deliberately wrong archived budget and a noncontaining trial cell in
temporary copies, without any ODE. Results are in `reduction-checks.json`.
C++ diagnostics compiled successfully. They are manual
targets and are not added as automatic campaign reruns to CTest. No general
regression/campaign acceptance is claimed from these local checks.

### Exact local problem and event

Packed state is `(x=ln(Tinf/1e8 K), eta_e [MeV], eta_mu [MeV])`; time is seconds.
Observation 32 is `1953.125 yr`, strictly between accepted ordinals 20 and 21.

| Quantity | Exact saved decimal |
|---|---|
| t_L | 56152132796.889648 |
| t_obs | 61635937500 |
| t_R | 64268647225.388481 |
| y_L | (0.11500605945683190, -1.8594970078177559e-7, -4.6695379642316489e-7) |
| y_R | (0.12930946174191599, -2.1089622937816103e-7, -5.2843718654675489e-7) |
| T_L [K] | 112188023.55362001 |
| T_obs O1 [K] | 113283677.74343133 |
| T_obs O2 [K] | 113283677.74347705 |
| T_R [K] | 113804225.02152768 |
| actual cache knot index | 95 (zero based) |
| exact knot T [MeV] | 0.0097145208258697432 |
| exact knot T [K] | 112732332.96598108 |
| Cstar at knot [erg/K] | 1.4700798497005163e38 |

The knot lies strictly inside `(t_L,t_obs)`, not between observation and right
endpoint. A bounded event solve uses repeated independent O3 IVPs from the saved
left state and the condition `Tinf(t)-T_knot=0`; no endpoint interpolation or state
projection is used. Fifty bisection evaluations give the **numerical-IVP** bracket
`[58866820737.305344, 58866820737.305351] s`, width
`7.62939453125e-6 s`; the common chosen split is `58866820737.305344 s`.
This floating-point bracket is **not a certified exact-IVP event-error bound**.
The split O1/O2/O3 first-leg states agree and end one T-MeV ulp below the knot
(residual `-1.7347234759768071e-18 MeV`). Each second leg starts from its own first
leg without projection and with fresh numerical state. Common event-time error
must be treated separately in any future production method.

### A/B: rk8pd hierarchy and split prototype

O1/O2 retain governed relative tolerances `1e-12/1e-13` and absolute vectors
`(1e-17,1e-23,1e-23)/(1e-18,1e-24,1e-24)`. Diagnostic O3/O4 use
`1e-14/(1e-19,1e-25,1e-25)` and `1e-15/(1e-20,1e-26,1e-26)` respectively.
All start from the same left state and full local-interval initial trial step.

| Method | Endpoint (x, eta_e, eta_mu) | Accepted / rejected | RHS calls |
|---|---|---|---|
| unsplit-O1 | `(0.12472490938574241, -2.0287605965623057e-07, -5.0871537810150439e-07)` | 11 / 15 | 339 |
| unsplit-O2 | `(0.12472490938614615, -2.0287605965622805e-07, -5.0871537810149634e-07)` | 16 / 22 | 495 |
| unsplit-O3 | `(0.12472490938613169, -2.0287605965622813e-07, -5.0871537810149666e-07)` | 20 / 19 | 508 |
| unsplit-O4 | `(0.12472490938613229, -2.0287605965622808e-07, -5.0871537810149655e-07)` | 23 / 17 | 521 |
| split-right-O1 | `(0.12472490938613233, -2.0287605965622813e-07, -5.0871537810149666e-07)` | 1 / 0 + 1 / 0 | 28 |
| split-right-O2 | `(0.12472490938613233, -2.0287605965622813e-07, -5.0871537810149666e-07)` | 1 / 0 + 1 / 0 | 28 |
| split-right-O3 | `(0.12472490938613233, -2.0287605965622813e-07, -5.0871537810149666e-07)` | 1 / 0 + 1 / 0 | 28 |

All GSL statuses are zero, including O4: no reported tolerance/roundoff pathology.
O4 is already near the binary64 precision floor and was not tightened further.
All three split endpoints are bit-identical. Their O1/O2 `d` is zero, but the
nonzero `F_O` floor remains: thermal `U/(0.20 F_i)=0.09645486279170176` and both
chemical ratios approximately 0.0962. This is a **diagnostic criterion evaluation**,
not a governed checkpoint pass or a bound on shared event-location error.

Absolute thermal differences from split rk8pd O3 decrease:
O1 `3.899242040361628e-13`, O2 `1.3822276656583199e-14`,
O3 `6.3837823915946501e-16`, O4 `4.163336342344337e-17`.
The supported local value is approximately `x=0.1247249093861323`; it is not a
reinterpreted campaign final state. No fitted formal convergence order is claimed.

Both O1 and O2 first attempt the full `5483804703.1103516 s` interval from t_L;
that attempt crosses the knot and is rejected, with thermal embedded error
`1.4164053722969252e-7`. First accepted knot-crossing attempts follow:

| Level | Attempt | Start t [s] | Trial h [s] | x before | x after | embedded x error |
|---|---|---|---|---|---|---|
| unsplit-O1 | 22 | 58866141304.730392 | 69270965.383422658 | 0.11984488367529816 | 0.11996765306154851 | 1.3058978833803524e-13 |
| unsplit-O2 | 32 | 58865934946.203499 | 5276773.8541865023 | 0.11984451788895019 | 0.11985387127282893 | 3.6891518864879692e-15 |
| unsplit-O3 | 30 | 58865928186.397911 | 1062825.0839249322 | 0.11984450590666981 | 0.11984638984242106 | 3.7782772031809173e-16 |
| unsplit-O4 | 28 | 58866546776.825699 | 302932.94581294747 | 0.11984560240481074 | 0.11984613937522641 | -3.1637776619705901e-17 |

Every trial, including rejected steps, is retained in `*.steps.tsv.gz`; every RHS
stage and actual cache cell is in `*.rhs.tsv.gz`. The reduction also retains all
near-knot attempts. O1 uses cell94/cell95 on 131/208 RHS calls; O2 on 218/277.
Every logged trial temperature belongs to its selected cell. Stages on opposite
sides correctly choose opposite cells; this is not an ownership defect.

### C: independent local solver

GSL RKF45 has a distinct embedded tableau and error estimator from rk8pd; it shares
GSL's controller and the identical physical RHS. It is numerical-method
independence, not independent physics or an independent GSL implementation.

| Run | Endpoint (x, eta_e, eta_mu) | Accepted / rejected | RHS calls |
|---|---|---|---|
| independent-unsplit-O2 | `(0.12472490938612608, -2.028760596562285e-07, -5.0871537810149793e-07)` | 11 / 10 | 127 |
| independent-unsplit-O3 | `(0.12472490938609974, -2.028760596562284e-07, -5.0871537810149719e-07)` | 17 / 18 | 211 |
| independent-right-O2 | `(0.12472490938615989, -2.0287605965622848e-07, -5.0871537810149772e-07)` | 3 / 1 + 3 / 1 | 50 |
| independent-right-O3 | `(0.12472490938613265, -2.0287605965622816e-07, -5.0871537810149687e-07)` | 5 / 2 + 4 / 1 | 74 |

Split RKF45 O3 differs from split rk8pd O3 by
`(3.191891195797325e-16, 2.6469779601696886e-23, 2.1175823681357508e-22)`.
Split RKF45 O2's thermal discrepancy is `2.756128658631951e-14`; O3 improves it.
Unsplit RKF45 is **nonmonotonic**: its O2/O3 thermal discrepancies from split rk8pd
are `6.2450045135165055e-15` and `3.2585045772748344e-14`. This behavior supports
caution about ordinary unsplit error estimation at the kink; it does not establish
an exact reference value from one tight solve.

### D: preselected nearby smooth controls

Full identities and states are in `cases.tsv` and `analysis.json`. All have zero
actual Cstar knots in the entire accepted bracket, checked against the dumped
160-point cache. Their longer/shorter local durations avoid selecting only tiny
or unusually easy replay intervals.

| Observation | Replay duration / failing replay | O1/O2 d (all components) | O2/O3 absolute differences (x, eta_e, eta_mu) | Max O1/O2 U/(0.20 F_i) |
|---|---|---|---|---|
| 26 | 1.947604871 | (0,0,0) | (2.7755575615628914e-17, 2.6469779601696886e-23, 1.0587911840678754e-22) | 0.09040538801 |
| 36 | 0.9248656246 | (0,0,0) | (0.0, 0.0, 0.0) | 0.09641886794 |
| 18 | 1.598069126 | (0,0,0) | (0.0, 0.0, 0.0) | 0.08898386153 |

O1/O2 each use one accepted step without rejection for all three controls.
Observations 18 and 36 are bit-identical at all three levels; observation 26 O3
differs by roundoff-scale amounts. All controls pass the unchanged criterion.
This corrects an intermediate chat claim that all three controls were identical
across all levels. Three controls discriminate this failure; they do not prove
O1/O2 sufficiency for every future smooth bracket.

### Cell ownership, cache and complete RHS

`StarContext.cpp:41-53` retains the hinted interval when
`grid[i] <= T <= grid[i+1]`; otherwise `upper_bound` chooses the lower cell index.
At exact equality an inclusive lower-cell hint can remain 94 while a warm upper
hint gives 95. Both produce the exact same Cstar. At the nearest distinct
representable temperatures below/above, cells are 94/95. Log rounding causes a
one-ulp Cstar plateau across these adjacent inputs, not a jump or wrong cell.

The payload is unchanged, but the mutable `last_i` accelerator **does persist
through shared ordinary context owners**. Disposable packed/GSL state is isolated.
No value/RHS dependence on the hint is observed: all 15 main-style/fresh-replay/
oppositely-warmed RHS components agree exactly, both warm directions give identical
complete adjacent RHS rows, and a separate cold StarContext rebuild matches all
160 T/C entries. `context-route.tsv`'s `replay_cold` column denotes fresh state,
not a cold cache; cold cache rebuilding is a separate explicit test.

For `C(T)=C_i+(C_{i+1}-C_i) log(T/T_i)/log(T_{i+1}/T_i)`, at the knot:

- Cstar is continuous at `1.4700798497005163e38 erg/K`.
- dC/dT below: `1.257951884727429e30 erg/K^2`.
- dC/dT above: `1.3524168435956085e30 erg/K^2`.
- Slope jump: `9.44649588681795e28 erg/K^2`, or `7.509425441073052%`.

Fixed chemical state for this RHS audit is the O3 event state
`(-1.9436528131428633e-7,-4.877399106566095e-7) MeV`; the sampled packed x values
are the nearest ones giving distinct temperatures below/at/above the knot:
`0.11984608801947824`, `0.11984608801947835`, `0.11984608801947845`.

| Quantity | Below | At | Above |
|---|---|---|---|
| T_MeV | 0.0097145208258697414 | 0.0097145208258697432 | 0.0097145208258697449 |
| xdot | 1.7725709358350855e-12 | 1.7725709358350855e-12 | 1.7725709358350849e-12 |
| eta_dot_e | -3.0870345128172983e-18 | -3.0870345128172979e-18 | -3.0870345128172971e-18 |
| eta_dot_mu | -7.616737809961695e-18 | -7.6167378099616935e-18 | -7.6167378099616919e-18 |
| Cstar | 1.4700798497005163e+38 | 1.4700798497005163e+38 | 1.4700798497005163e+38 |
| Pnet | 2.9376025975904325e+34 | 2.9376025975904325e+34 | 2.9376025975904316e+34 |
| Pdir | 3.3583239602677939e+34 | 3.3583239602677939e+34 | 3.3583239602677939e+34 |
| LH | 3.2824551920485368e+22 | 3.2824551920485397e+22 | 3.2824551920485452e+22 |
| DeltaLnu | 4.9236827879519451e+22 | 4.9236827879519468e+22 | 4.923682787951956e+22 |
| Lnu_eq | 4.9858923559729203e+32 | 4.9858923559729246e+32 | 4.9858923559729354e+32 |
| Lnu_full | 4.9858923564652884e+32 | 4.9858923564652927e+32 | 4.9858923564653035e+32 |
| Lgamma | 3.7086243911599089e+33 | 3.7086243911599089e+33 | 3.7086243911599135e+33 |
| sigma_e | 7.4466405596551908e+35 | 7.4466405596551908e+35 | 7.4466405596551908e+35 |
| sigma_mu | 8.1526982270859939e+34 | 8.1526982270859939e+34 | 8.1526982270859939e+34 |
| R_e | -7.85357676745236e+34 | -7.8535767674523665e+34 | -7.8535767674523794e+34 |
| R_mu | -1.0708263105552357e+34 | -1.0708263105552366e+34 | -1.0708263105552385e+34 |

Units: xdot s^-1; eta dots MeV/s; Cstar erg/K; powers/luminosities erg/s;
sigma and R counts/s. The unchanged prescribed ordinary source is
`S=(Bdot,0,0,0)` in (n,p,e,mu), with the Bdot above; Sigma and sigma retain the
accepted moving-reference tangent projection. `DeltaPbeta=LH-DeltaLnu`; no
second heat source is inserted. Source contributions are constant across this
fixed-state temperature scan. The complete physical RHS itself is continuous.

The verified discontinuity is in **partial xdot/partial x**, holding eta fixed:
finite differences at `h=1e-8` give left `-4.264720913338948e-12`, right
`-4.393125860231013e-12`, jump `-1.284049468920649e-13 s^-1`.
The independently computed Cstar-slope prediction
`-xdot * T_knot * (C'_above-C'_below)/Cstar`
is `-1.2840499959160055e-13 s^-1`. The time derivative of xdot along the solution
therefore has the corresponding jump (approximately `-2.2761e-25 s^-2`).
Chemical derivative differences tend to zero with h (zero at h=1e-9 in this
floating-point scan). Pnet, direct power, beta terms and luminosities show no
additional finite jump. The h-sequence and all contribution derivatives remain
in `analysis.json`; cancellation at binary64 limits is not interpreted as exact
analytic differentiability proof for unrelated physics.

### Main-trajectory crossing and scale audit

Only archived main evidence was read. Accepted step21 spans `[t_L,t_R]` in one
`8116514428.4988327 s` step, cell94 to95. Cumulative rejections rise 3→7 (four
rejections during this evolve call); the preceding next-step proposal was
`17813012587.102097 s`. The sole t1 ceiling remained `15778800000000 s`, with zero
intermediate observation ceilings. The main solver adapted but did **not** split
at, stop at or explicitly resolve the knot. Its rejected trial endpoints were
not archived, so their individual positions/errors cannot be reconstructed
without a prohibited main rerun.

The two authorized supplemental local O3 continuations to saved t_R yield:

| Component | Saved main right | Split rk8pd O3 right | Split RKF45 O3 right | Main minus rk8pd |
|---|---|---|---|---|
| 0 | 0.12930946174191599 | 0.12930922373678003 | 0.12930922373678039 | 2.3800513596072825e-07 |
| 1 | -2.1089622937816103e-07 | -2.1089623236974881e-07 | -2.1089623236974895e-07 | 2.9915877825761639e-15 |
| 2 | -5.2843718654675489e-07 | -5.2843719586873925e-07 | -5.2843719586873947e-07 | 9.3219843609531152e-15 |

The thermal discrepancy is `2.3800513596072825e-7` while the methods agree within
`3.6082248300317588e-16`. BASELINE main uses rtol1e-7/atol1e-12 in x; the local
thermal discrepancy is about 18.4 times its simple endpoint tolerance scale.
That scale is not a rigorous global-error bound and this is **not** a retroactive
new gate verdict. It is a material local main-step finding for independent review.
No main defect or campaign convergence is adjudicated without the deferred
hierarchy. Reconstruction is intentionally judged against an ULTRA-derived
uncertainty budget; main-step acceptance does not waive that independent budget.

Production `PassiveCheckpointOutput.cpp:399-415` uses:
`M_O=max(abs(O1),abs(O2))`, `D_O1=atol_O1+rtol_O1*M_O`,
`F_O=max(atol_O2+rtol_O2*M_O,64 ulp(M_O))`,
`M_i=max(abs(y_L),abs(y_R),abs(O1),abs(O2))`,
`F_i=max(atol_ULTRA+rtol_ULTRA*M_i,64 ulp(M_i))`,
`U=2 max(d,F_O)`. ULTRA here is the fixed **budget configuration**, not an executed
trajectory. The independent Python reduction uses explicit libc fma to match the
qualified AppleClang arithmetic and reproduces all archived values exactly.

| Component | d | D_O1 | F_i | U | d/D_O1 | U/(0.20 F_i) |
|---|---|---|---|---|---|---|
| x | 4.0374648069274599e-13 | 1.2473490938614616e-13 | 1.2931946174191598e-12 | 8.0749296138549198e-13 | 3.2368362848836014 | 3.1220859973768409 |
| eta_e | 2.5146290621612041e-21 | 2.0288605965623057e-19 | 2.1090622937816101e-18 | 4.0577211931246114e-20 | 0.012394291980543082 | 0.096197281727724571 |
| eta_mu | 8.0468129989158532e-21 | 5.0872537810150436e-19 | 5.2844718654675484e-18 | 1.0174507562030088e-19 | 0.015817596969401236 | 0.096267969827954492 |

There is no wrong absolute component, unit mismatch, cross-solution magnitude
substitution, stale BA12R tolerance or double application of 0.20. The failed
thermal gate is real. The chemical coordinates pass. Z, tangent, actual-potential
ledger, Cstar table values, source, grid and physics remain unchanged.

### Hypothesis disposition, limitations and next decision

| Hypothesis / classification | Bounded evidence and disposition |
|---|---|
| H1 / A | Supported: continuous RHS with a verified thermal derivative kink; tight unsplit convergence, stable split hierarchy, independent split agreement, smooth controls pass. |
| H2 / B | Not supported as a general smooth-segment insufficiency; O1/O2 are insufficient for this unsplit crossing. Three smooth controls cannot prove universality. |
| H3 / C | No defect found: exact saved failure reproduction, unchanged source/input/package hashes, exact route/warm equality and cold payload equality. |
| H4 / D | No defect found: all actual trial cells contain T; exact-knot ownership differs by valid inclusive hint but Cstar/RHS do not. |
| H5 / E | No broader RHS defect demonstrated; smooth controls and complete RHS audit localize the effect. Main-step discrepancy is a separate review finding. |
| F | Not the primary classification. Exact-IVP event-error certification and general production qualification remain unresolved implementation prerequisites, not an unexplained failure reproduction. |

The clean prototype is a uniform, event-aware piecewise-smooth local reconstructor.
It must be proposed separately; accepted ADR-0017 is untouched. A common O3 root
used by both tiers can hide common error, so identical split O1/O2 is not sufficient
production evidence. Event error propagation, multiple/nonmonotone crossings,
endpoint ownership, schedule passivity and broad regression are future requirements.
Uniform tighter tolerances are an alternative to assess, but O4 is near roundoff
and one bracket does not authorize a global replacement. No fallback/retry has
been added to production. No candidate is eligible.

The primary executable took `557.0286429999396 s` externally (`556.528403792 s`
internally; initial context ready at `271.791057625 s`). The supplemental executable
took `316.4879429170396 s` externally (`315.965795583 s` internally).
Total external execution time was `873.5165859169792 s` (~14.56 min), sequential.
There were **80 local integrations**: 2 unmodified-production reproductions,
4 traced unsplit levels, 50 event-root evaluations, 6 split legs, 6 independent
method legs, 9 smooth controls, 1 event-state audit solve, and 2 right-endpoint
extensions. No full source/control trajectory. Per-solve integration-only timings
and counts are in solutions.tsv; context construction/full-currentness overhead
is included only in the executable totals. These timings do not qualify cluster
execution or motivate changing currentness checks.

### Requested 38-item disposition

1. Entry SHA: `7d2777561a17b56b466fb024ef13104da6c0e86a`.
2. Branch/worktree: `analysis/phase6a1-obs32-reconstruction-forensics` at the absolute path in the declaration.
3. Identity: source BASELINE observation32, accepted endpoints20/21, strict interior, one knot.
4. Exact t_L/t_obs/t_R: `56152132796.889648 / 61635937500 / 64268647225.388481 s`.
5. Exact left/right states: local-problem table above; authenticated cases.tsv.
6. Knot: `0.0097145208258697432 MeV =112732332.96598108 K`, cache95.
7. Location: strictly between t_L and t_obs; numerical root bracket above.
8. Cstar: continuous; dC/dT jumps by 7.509425441073052%, slopes above.
9. Archived O1/O2 failure: **exactly reproduced**, all six states and budgets.
10. O1: 11 accepted /15 rejected,339 RHS evaluations.
11. O2: 16 accepted /22 rejected,495 RHS evaluations.
12. O3: x=`0.12472490938613169`,20/19, converges toward split value.
13. O4: x=`0.12472490938613229`,23/17, GSL success; near binary64 floor.
14. Unsplit trend: thermal distance to split decreases 3.90e-13→1.38e-14→6.38e-16→4.16e-17.
15. Split O1: x=`0.12472490938613233`, two one-step legs, no rejections.
16. Split O2: exactly the same endpoint as split O1.
17. Split O3: exactly the same endpoint as split O1/O2.
18. Split convergence: all component differences zero; nonzero uncertainty floor remains; common event error not independently bounded.
19. Split–unsplit: exact vector differences in analysis.json; thermal values in item14.
20. Independent solver: split RKF45 O3 agrees within 3.19e-16 in x; unsplit nonmonotonic trend retained.
21. Smooth controls: observations26,36,18 all pass; full differences above.
22. Cell ownership: inclusive-hint/upper_bound rule valid; no floating-point ownership defect found.
23. RHS: complete below/at/above table above; RHS continuous, thermal derivative kink verified.
24. Main crossing: one accepted knot-spanning step after four rejections; local O3 right-end discrepancy disclosed above.
25. D_O1/F_i: exact independent reproduction; no scale/unit/0.20 defect.
26. Classification: **A**, with H2 limited to this unsplit crossing; no evidence requiring B/C/D/E classification.
27. Production defect: none identified under existing ADR; the fail-closed refusal is correct.
28. ADR-0017 change required: **YES** to adopt event-split reconstruction; not made here.
29. Proposed fix: uniform event-aware piecewise-smooth local reconstruction, subject to separate owner acceptance.
30. Prototype: positive local split evidence; not production qualified.
31. Campaign rerun authorized: **NO**.
32. Historical BA12: **FAIL**, immutable.
33. Historical BA12R: **FAIL**, immutable.
34. Physical BNV model/rate selected: **NO / NO**.
35. Canonical master changed: **NO**, remains `9a5c1eca40758e4c39f9f80db7324e12ad23e505`.
36. Blockers: owner architecture decision and event-error qualification before production adoption; campaign remains stopped.
37. Nonblocking findings for completing this investigation: log(T) interpolation clarification; harmless shared cache hint; missing original full cache dump explicitly bounded by reproduction; finite control sample. Material main-step discrepancy is routed to the same owner review, not silently dismissed.
38. Recommended next action: **one independent scientific/numerical review of this evidence and the separately proposed architecture change**, then an explicit owner decision. No campaign or physical specialization follows automatically.
