# Bounded actual-path qualification proposal

> Stage A restart-selection amendment accepted by the owner after proposal
> `3557da9`: for each role/tier's quarter-cap main, use the latest accepted
> **noninitial** endpoint at or before half-time, with a valid nonzero tail.
> Fail explicitly if none exists; no passive checkpoint or initial-state fallback.
> This supersedes only the Stage A selection sentence below; its original wording
> is retained as historical provenance. Stage B and all numerical criteria are
> unchanged. Preparation is authorized; another Stage A trajectory is not.
> See `../restart-preparation/OWNER_APPROVAL.md` and the exact applied patch at
> `../restart-coverage/proposed-selector.patch`.


**NOT AUTHORIZED, NOT RUN, NOT LAUNCH-READY.** This is a future run plan, not a
submission script. Native pair assembly is blocked as recorded in `REPORT.md:1`.
Commit the eventual runner, snapshot hashes, input manifest and completed
non-evolution receipt before requesting owner authorization. A permission to run
cannot substitute for the missing assembly evidence.

Reuse `UninterruptedBnvTrajectory`, `BreakpointStepper`,
`PassiveObservationSchedule`, `Rk8pdCheckpointReconstructor` and
`ComputeCheckpointR20`; do not alter their equations or acceptance logic.
The existing APIs are in `CompactStar/Physics/BNV/PassiveCheckpointOutput.hpp:21`.
The numerical criteria below retain `../../PHASE6_CAMPAIGN_ACCURACY_PREDECLARATION.md:10`
and `../EXPERIMENT.md:40`. The shorter horizons are a new proposed staging choice,
not a replacement acceptance test for the 14,648.4375-year experiment.

## Prerequisites without a trajectory

Restore the owned central/sequence/chemical/tangent dependencies; authenticate
all pinned input bytes and the conditional sidecar. Serialize the pair and load
it twice with separate mutable caches/tokens. Require exact semantic identity,
identical common configuration and deterministic reconstructed values. The only
physical difference is P2 neutron sink versus exact-zero source/product sector.

Evaluate each role at Tinf=1e8 K, eta=(0,0), (1e-7,-2e-7), (-1e-7,2e-7) MeV,
at t=0 and t=462269531250 s: 12 RHS snapshots, with no integrator calls.
Verify source/control ordinary terms agree at the same supplied state; the
control reduces exactly to the ordinary RHS. Check moving-source signs,
actual-potential correction, all eleven power channels and R18/G_P <=1e-10.
Require missing/foreign/modified certificate, source hash, owner, token, thermal
table, role, and serialized-state mismatches to refuse before integration.
These are the predeclared states in `../OWNER_APPROVAL.md:57`.

Inspect the initial Cstar cell/knots and round-trip state/units, solver settings,
role and dependency identities. Distinguish passive checkpoint replay from
integrator restart: the present API exposes `Integrate(start,end,state)` and
`Replay`, but no full process-state serializer. Do not claim bitwise continuation
of an adaptive solver from a three-component state alone
(`CompactStar/Physics/BNV/PassiveCheckpointOutput.hpp:83`, `:137`).

## Exact proposed execution bounds

Each stage starts BOTH roles independently at original t=0, x=0, eta=(0,0),
Tref=Tinf=1e8 K, spin=0. B0=7.616906518718875e56, source Bdot=-2.413652026364133e37
count/s; control Bdot=0. P2 is the existing mathematical benchmark, not a decay
prediction (`../EXPERIMENT.md:11`).

| Stage | Final time | Purpose | Initial reference cap H=duration/256 |
|---|---:|---|---:|
| A | 31,557,600 s (1 Julian year) | Startup, matched reduction, state restart and prefix accounting | 123,271.875 s |
| B | 64,629,964,800 s (2048 Julian years) | Coupled response and actual breakpoint/reconstruction coverage | 252,460,800 s |

Advance to B only after A's capped policies and references pass. B targets the
vicinity of the historical source observation-32 failure (1953.125 years;
`../../PHASE6A1_CONTROLLED_BNV_RESUME.md:1`), but does **not** assume a crossing.
No crossing means actual event coverage remains UNQUALIFIED; do not lengthen the
run automatically. Neither stage qualifies the full 14,648.4375-year interval.
No warm historical window or 500 kyr campaign is included.

For each horizon and role, run the following fixed matrix:

* Main RKF45: BASELINE rtol=1e-7, atol=(1e-12,1e-18,1e-18);
  REFINED 1e-9/(1e-14,1e-20,1e-20); ULTRA 1e-11/(1e-16,1e-22,1e-22).
  Initial step 1 s. Each tier compares uncapped, 38,723,807,757 s,
  19,361,903,878.5 s, and 9,680,951,939.25 s caps. Preserve all failed policies.
* Independent references: rk8pd at 1e-14/(1e-19,1e-25,1e-25), cap H;
  rk8pd at 1e-15/(1e-20,1e-26,1e-26), caps H/2 and H/4;
  independent RKF45 at the latter tolerances, caps H/4 and H/8.
  Every reference owns its events/caches and starts at original t=0.
* Passive scientific grids: uniform 240, 480, 960 and 1920 intervals, including
  both endpoints. Evaluate the nested grids from the same accepted history.
  Rerun each main policy/tier with only the two endpoint observations and require
  bit-identical main accepted states/steps. Events are not quadrature nodes.
* For the quarter-cap main at each tier/role, branch one restart tail at its first
  accepted endpoint at or beyond half the horizon. Round-trip its exact state and
  dependency identities, rebuild independent caches, and compare against the
  uninterrupted reference. Carry all prefix errors and cumulative energies;
  retain original t=0 for source B(t). Restart compatibility is a numerical
  comparison, not a claim of bitwise adaptive continuation.

Maximum per stage: 24 candidate main runs, 10 reference runs, 24 schedule repeats,
6 restart tails. No unlisted retries, horizon changes or refinement extensions.
Identical capped histories do not establish cap convergence.

## Norms, acceptance and stops

At all scientific times and accepted endpoints use the componentwise infinity
norm of (|y-y_ref|+U_ref+U_output)/D, where
D_i=atol_i+rtol*max(|y_i|,|y_ref_i|); require <=1 for each role/tier. Require
U_ref<=0.05 D_ULTRA and reconstruction<=0.20 F_i using the existing F_i definition.
Temperature uncertainty uses exp(U_x)-1. U_ref is twice the largest adjacent
refinement or finest independent-method difference plus actual-path roundoff,
at least 64 ulps. Above independently bounded roundoff, refinement discrepancy
must contract by <=0.5. Operation-count or higher-precision evidence must bound
accumulated roundoff; missing evidence is UNQUALIFIED, not zero.

For the eleven existing powers, use noncancelling
G_P=max(1,|Pdir_actual|+|LH|+|DeltaLnu|+|Lnu_eq|+|Lgamma|+|Lother|+|Lout_fluid|).
Error+reference+output uncertainty must be <=1e-7 G_P for BASELINE and <=1e-9 G_P
for REFINED/ULTRA; reference allocation <=5e-11 G_P. R18/G_P<=1e-10.
See `tests/bnv/breakpoint_qualification.cpp:145` for the exact channel list.

R20 uses the original scientific-grid trapezoids and
N=max(1 erg,|DeltaU|,|DeltaEchem|,integral(|Lnu_full|+|Lgamma|+|Lother|)dt),
with DeltaEchem in erg. Require
(|residual|+U_main+U_reconstruction+U_event+U_quad)/N<=2e-4,
and combined reconstruction/event allocation <=5e-6 N. U_quad is twice the
largest adjacent nested-grid discrepancy; final contraction <=0.5 unless all
differences are below independent rounding bounds. Publish every uncertainty
term and avoid cancellation (`CompactStar/Physics/BNV/src/PassiveCheckpointOutput.cpp:109`).

Form matched source-minus-control uncertainties by addition, with no assumed
cancellation. Require discrepancy+uncertainty <=D_source+D_control for states;
use summed corresponding budgets for powers/energies. A resolved signal also
needs U_delta<=0.10|delta_ref|; otherwise label the observable/window SIGN UNRESOLVED.
Exact-zero control chemistry must use the existing owner-checked invariant.

At every actual event check localization, exact-endpoint bypass, both nextafter
neighbors, continuity/transversality and event uncertainty. Preserve full-current
checks at construction/restart/integration boundaries and cheap checks during RHS.
Require rejection <=10%, finite progress, exact final time, Tinf in [1e6,1e9] K,
historical and conditional utilization <=1, depletion<=1e-6, finite diagnostics,
and all source/fate/currentness guards. Stop the affected run on any violation;
keep its partial output. Missing error terms stop qualification. Historical
BA12/BA12R remain FAIL; report their original c1/c2/c3 separately.

## Resources and authorization boundary

Proposed Slurm allocation per stage: phyall/phyn004, one node/task/CPU, 4 GiB,
30-minute **aggregate stage wall ceiling**, single-thread libraries. Timeout is
INCOMPLETE, with preserved outputs; no automatic resubmission. Build separately
with the existing Linux lock, Release and -ffp-contract=off, four CPUs, 8 GiB,
15 minutes. These are ceilings, not measured runtime predictions.

All large outputs/builds/tmp stay on /mnt/sdd; Git contains manifests, code,
compact receipts and the fixed plan. An eventual authorization request must name
the exact runner commit, source/control snapshot SHA-256s, input/sidecar hashes,
these two horizons, matrix, ceilings and stopping rules. Those snapshot identities
are currently absent, so no executable launch request is issued by this delivery.
