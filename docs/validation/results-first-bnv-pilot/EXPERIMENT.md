# RF-BNV-P2-NATIVE-BOUNDED-v1 — pre-execution declaration

**Defined, not authorized or executed.** Commit this declaration before any
numerical execution. Release depends on `OWNER_DECISION.md`, completed conditional
sidecar acceptance, exact native context binding and affected prerequisite tests.
No claimed result or baseline exists. Authority: ADR-0015 §§9,14,15; accepted
ADR-0017 Amendment 1; `PHASE6_CAMPAIGN_ACCURACY_PREDECLARATION.md:1` (candidate
accuracy policy adopted for scoring this separately named diagnostic, not ratified
for the full campaign).

## Physical and mathematical definition

Use the existing mathematical P2 neutron sink (`campaign_fixture.hpp:90`), same
spin-off whole-star free-gas Structure-1, rho_c=1.10e15 g/cm^3, radial80000,
EOS parameter8192, MU-only SMe=1e-51 and SMmu=2e-51 in existing units. No physical
decay rate is selected. The native target-0 profile, chemical certificate, thermal
table, coefficient table and native 21-star certificate are pinned by
`evidence/preflight.json`. No Mac certificate is attached to a native star.

| Parameter | Source and matched control |
|---|---|
| Common initial state | t=0; Tinf=1e8 K; eta_e=eta_mu=0 MeV; x=ln(Tinf/1e8 K)=0 |
| Interval | 0 to 462269531250 s = 14648.4375 Julian years, the previously declared bounded-from-zero interval |
| Source | Bdot/B0=-1e-12 yr^-1, constant; uniform proper neutron sink proportional to n_n; P2 Eesc_fluid=Eesc_star=EX=0 |
| Native normalization | B0=7.616906518718875e56 counts; Bdot=-2.413652026364133e37 count/s (evaluate using the pinned binary64 inputs) |
| Expected depletion | 1.46484375e-8; 0.0146484375 of the 1e-6 hard ceiling. This arithmetic is not measured certificate utilization |
| Control | Same B0, T, eta, ordinary microphysics, thermal/envelope/spin/solver configuration and observations; ExactZeroHistory + ExactZeroPartition, no BNV events or products |
| Numerical start | initial step 1 s; three original tolerance tiers, unchanged |
| Main cap comparison | unchanged adaptive policy, 3.8723807757e10 s, half and quarter; report all outcomes separately |
| Scientific times | Uniform grid with 240 intervals; nested 480/960 and witness1920 for R20. Also report accepted endpoints; do not force integration onto these times |

This is a new native diagnostic card derived from the existing benchmark and
bounded window. It is **not** the original CPL-P2-LINEAR-QSS-v1 500 kyr/8193-point
campaign. Do not edit that card, replace its absolute Bdot literal, claim reached
QSS, or extend the duration after inspecting the answer. Existing paired assembly
is in `tests/bnv/breakpoint_fixture.hpp:9`; its native binding gate remains pending.

## Fixed numerical acceptance and stops

Use the complete equations, reference construction, budgets and stop rules in
`../PHASE6_CAMPAIGN_ACCURACY_PREDECLARATION.md:10`, restricted to this interval:

- BASELINE rtol1e-7/atol(1e-12,1e-18,1e-18); REFINED1e-9/(1e-14,1e-20,1e-20);
  ULTRA1e-11/(1e-16,1e-22,1e-22). State discrepancy + reference + output
  uncertainty <= each tier's component budget. Reference allocation <=0.05 of
  ULTRA, output reconstruction <=0.20 F_i. Exact exponential temperature transform.
- Independent rk8pd and RKF45 references start at the same original t=0 state,
  with independent events/caches and the prescribed tolerance/cap refinements.
  No warm restart can erase prefix error. Quantify actual-path roundoff; 64 ulps
  alone is not an accumulated-error proof. Missing evidence stays UNQUALIFIED.
- Eleven power-channel errors <=1e-7 gross G_P (BASELINE) or 1e-9 G_P (REFINED,
  ULTRA); reference allocation <=5e-11 G_P. R18/G_P<=1e-10.
- R20 total residual plus main, reconstruction, event and quadrature uncertainty,
  divided by the original N_R20, <=2e-4; reconstruction/event <=5e-6. Scientific
  grid trapezoids only, with nested-grid contraction <=0.5 above roundoff and
  twice the largest adjacent quadrature discrepancy. Publish every contribution.
- Matched differences sum source and control uncertainty without assumed
  cancellation. A resolved nonzero signal needs U_delta<=0.10*abs(delta_reference).
  Otherwise report SIGN UNRESOLVED for that observable/window.
- All-tier passive-schedule bit identity, endpoints, event/reconstruction bounds,
  <=10% rejection, finite progress, exact final time; runtime historical and
  sidecar conditional utilization <=1; Tinf in [1e6,1e9] K. Every currentness,
  source/fate, finite-value and certificate refusal remains effective.

Score historical BA12/BA12R criteria separately; their failures never become
successor passes. A failed policy remains a failed result even if another passes.
No tolerance relaxation or post-result horizon reduction. The qualified Figure 1
requires the predeclared capped pilot and its refinements to pass the applicable
criteria; an unchanged-policy failure is disclosed. If a budget cannot be scored,
stop with its exact missing term. Do not commission broader numerical work.

## Mandatory output and resources

For both roles and every scored time retain: time_s/year, x, Tinf_K, both eta_MeV,
xi, signed reaction rates; Lnu_eq/full, DeltaLnu, Lgamma, Lother; Pdir_eq/actual;
LH plus each C_MeV_to_erg*eta_l*R_l contribution and DeltaPbeta; Lout_fluid,
Lesc_star, JX; Cstar, Echem, DeltaEeq, DeltaU; R18 and cumulative R20 with its
normalizer and all uncertainty terms; finite-T omitted floor; depletion and
per-group historical/conditional certificate utilization; accepted/rejected step,
breakpoint, reconstruction, currentness and validity diagnostics. Missing values
cannot be serialized as zeros. See existing `BnvDiagnostics.hpp` and preflight §16.

Only after passing, publish Figure 1: log10(Tinf/K) against log10(time/yr), source
and control, visibly labelled **CONTROLLED MATHEMATICAL BNV BENCHMARK**. Retain
t=0 in data and explicitly omit it from the logarithmic axis. Include uncertainties,
maximum absolute/fractional DeltaT and its time, first resolved departure on the
declared grid, direct/reaction integrated energies, ledger accuracy and maximum
utilization. Deliver numerical data and the plotting script together. No placeholder
curve is produced while blocked.

Proposed first execution allocation: Slurm, one task/one CPU, 4 GiB, 30-minute
hard wall cap, single-thread libraries. This is a resource ceiling, not a runtime
prediction. Retain partial output on timeout; request a separately declared
extension if needed. Any projected >5–6-hour extension must use Slurm. Large
outputs/builds/tmp stay in the external directory named in REPORT.md, with SHA-256
manifests in Git. **Actual allocations/runtime in this task: none/zero trajectory
seconds.** No full-campaign submission script is created.
