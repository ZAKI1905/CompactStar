# Phase-6A-1 BA12R segmentation diagnostic

Status: **NUMERICAL FORENSICS / DIAGNOSTIC EVIDENCE**

Classification: **NOT A BNV CANDIDATE; NOT A GOVERNED BASELINE; NOT A
PHYSICAL RESULT**

Disposition: **SEGMENTATION EFFECT CONFIRMED — ENDPOINT-ONLY AND
CHECKPOINT-SEGMENTED RKF45 DIFFER BY >10F — RETURN TO OWNER FOR
NUMERICAL-ARCHITECTURE DECISION**

## Scope and owner authorization

The human owner accepted the BA12R numerical-forensics result at
`0aab1c2b5b748585326de35499d1830d0c788fae` and authorized exactly two
local P2-source integrations on the same Mac/toolchain:

1. Arm S: the current checkpoint-segmented architecture through archived
   checkpoint index 240; and
2. Arm E: the identical source problem exposed to GSL through one final
   endpoint ceiling only.

Arm S was executed first and passed its exact archived-reproduction gate
before Arm E was started. No control, tolerance sweep, other trajectory, full
P2 run, full test suite, cluster job, candidate operation, or merge was
performed. Production source, tests, CMake, baselines, `ScaledRKF45`, Cstar,
and Phase-5D were not changed.

Historical BA12 remains **FAIL**, with maximum scaled difference
`1.7255917120989046`. Historical BA12R remains **FAIL**. This bounded
diagnostic does not requalify REFINED and does not change either result.

## Authenticated entry and provenance

| Item | Authenticated value |
| --- | --- |
| BA12R final SHA | `971b2bac47c320bd262d8686841e679f703e5ccb` |
| Forensics entry SHA | `0aab1c2b5b748585326de35499d1830d0c788fae` |
| Entry branch parity | local = upstream = live = forensics entry SHA; clean |
| Canonical master | local = origin = live = `bd697ffdc474863d7a39f42e17ad8e8dbf105e5d`; clean |
| Diagnostic branch | `analysis/phase6a1-ba12r-segmentation-diagnostic` |
| Diagnostic worktree | `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase6a1-ba12r-segmentation-diagnostic` |
| Historical raw-evidence aggregate SHA-256 | `81314ddceba77eaaf29f493928560e93f3a6b45d59f302c11fbec0e799da0945` |
| Historical BASELINE source SHA-256 | `f331b5bbfd6bdfb0f8e93095052cf1edc712219180370d23383bed1c18b31093` |
| Historical REFINED source SHA-256 | `1889336755a87377336e1b688e34032dab428304861bd2a63e142e27555b7a3a` |
| Archived ULTRA source SHA-256 | `584799310627541259cc1164c8a18f8950afb8dddd7f0f14a1f0b3e9dc557b52` |
| Archived ULTRA source steps SHA-256 | `80dff6f0620c156d942a16ab1f42a724ced156d1216fe44dab17bac8163ba085` |
| BA12R comparison JSON SHA-256 | `ca122ad9ad92911fa604ae2ec1ae5208c107a4b6e3dd96ed4aca7c9187e4cce4` |

The historical manifest and the per-file hashes above were authenticated
before execution. The first 241 archived rows supplied the checkpoint times
without regeneration or numerical alteration.

## Platform and build identity

The diagnostic matched the BA12R local numerical platform:

| Property | Value | Match? |
| --- | --- | --- |
| OS | macOS 26.6.2, build 25G83; Darwin 25.6.0 | yes |
| Architecture | arm64 | yes |
| Compiler | Apple clang 21.0.0 (`clang-2100.3.34.2`) | yes |
| Compiler target | `arm64-apple-darwin25.6.0` | yes |
| GSL | 2.7.1, `/opt/local/lib/libgsl.27.dylib` | yes |
| Build type | Debug | yes |
| Effective C++ flags | `-g -std=c++17 -arch arm64 -pthread -Xclang -fopenmp` | yes |

The untracked scratch executable linked the authenticated BA12R
`libCompactStar.a`, the same GSL/GSL CBLAS, OpenMP, Python, Zaki, and Confind
libraries, and the unchanged source at the forensics entry. Its SHA-256 was
`2631c5756f4b29201d4422ef6452687ef272a7de6966db16535e67f3023b2dd3`.
The scratch harness source SHA-256 was
`f4444792215e5b4d9b7b3a21ad4bde4f362159f1600fb87749638a00c4d4d732`.
The harness copied the public `ScaledRKF45::Integrate` call order into an
untracked audit adapter solely to serialize every accepted step; it made no
tracked or production change.

## Endpoint authentication

The endpoint artifact contains 241 rows, indices `0 ... 240`, comprising the
initial state and 240 positive-time ceilings. Its SHA-256 is
`43ec23ada72bfa59c4e89672ae2987e9927164db765f05afb7b45c924cc0c67e`.
Index 240 is exactly `462269531250 s`. Arm S used all positive-time entries in
their archived order. Arm E used only index 240 as its positive-time ceiling.

## Common mathematical and numerical configuration

Both arms used fresh independent process and context state with:

- fixture `CPL-P2-LINEAR-QSS-v1`, source ON;
- `Bdot/B0 = -1.0e-12 yr^-1` and
  `Bdot = -2.4136520263641375e37 count/s`;
- P2 partition, spin OFF, Me/Mmu enabled, De/Dmu disabled;
- initial `(x_state, eta_npe, eta_npmu) = (0, 0, 0)`, equivalent to
  `Tinf = 1e8 K`;
- identical frozen background, tangent, Z, Cstar cache, Ltilde, metric,
  envelope, source normalization, direct-energy semantics, currentness, and
  validity rules;
- GSL RKF45, `rtol = 1e-11`,
  `atol = (1e-16, 1e-22, 1e-22)`, initial `h = 1 s`;
- start `t = 0` and final `t = 462269531250 s`.

The same compiled binary and authenticated input bytes were used. Separate
work/output paths were required for process isolation and were not scientific
inputs. The sole scientific/numerical difference was:

| Property | Arm S | Arm E |
| --- | --- | --- |
| GSL `t1` ceilings | archived indices 1 ... 240 | archived index 240 only |
| Positive-time ceiling count | 240 | 1 |
| GSL step/control/evolve lifetime | continuous for the complete arm | continuous for the complete arm |
| Initial `h` | exactly 1 s, assigned once | exactly 1 s, assigned once |
| Intermediate observer samples | 239 before the common final endpoint | none |

Thus this experiment changes the integration ceilings/checkpoint segmentation,
not the equation, source, controller, tolerance, initial state, or final time.

## Arm S exact archived reproduction gate

Arm S completed successfully. It reproduced the archived ULTRA source through
index 240 exactly:

| Check | Result |
| --- | ---: |
| Trajectory rows | 241 |
| Binary64 state values compared | 723 |
| Unequal state values | 0 |
| Maximum absolute state difference | 0 |
| Unequal times | 0 |
| Deterministic serialized rows compared | 241 |
| Unequal deterministic serialized rows | 0 |
| Archived/Arm-S checkpoint step-summary rows equal | yes, all 240 |
| Currentness/validity | pass; every row valid |

Arm S recorded 289 accepted steps, 31 rejected steps, and 2210 RHS
evaluations. Its minimum and maximum accepted steps were `1 s` and
`1926123046.875 s`. Wall/user/sys timing was
`301.67 / 300.09 / 0.91 s`.

This exact gate passed before Arm E was executed.

## Arm E execution

Arm E completed successfully and emitted two trajectory rows: the initial
state and the single final endpoint. It recorded 232 accepted steps, 60
rejected steps, and 1985 RHS evaluations. Its minimum and maximum accepted
steps were `1 s` and `2883572725.0536194 s`. Wall/user/sys timing was
`300.80 / 299.71 / 0.90 s`.

The integration executable emitted `SEGMENTATION_DIAGNOSTIC PASS arm E` and
complete outputs. After the child completed, macOS `/usr/bin/time -lp`
reported a sandbox denial for `sysctl kern.clockrate` and returned nonzero.
This was a timing-wrapper limitation, not an integration failure; no retry was
made. Every Arm E row passed currentness/frozen validity. Both arms had maximum
frozen utilization `0.007921712322551693` and maximum
`|DeltaB|/B0 = 1.46484375124617e-08` on the bounded window.

## Predeclared endpoint discrimination

The BA12R scale was evaluated at archived index 240 using only archived
BASELINE, REFINED, and ULTRA values. `Q = 0`, as authenticated for the exact
binary64 round trip. Arm E was not included in `M`, `D_U`, or `F`.

| Component | Arm S final | Arm E final | `Delta_seg` | `F` | `10F` | `Delta_seg/F` | `>10F`? |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| `x_state` | 0.49240008843840127 | 0.49240008824076903 | 1.9763224390345613e-10 | 4.92410091788327e-12 | 4.92410091788327e-11 | 40.135701359348374 | yes |
| `eta_npe` / `eta_e` (MeV) | -2.5123474220450743e-07 | -2.5123474256442210e-07 | 3.599146754419561e-16 | 2.512447441357819e-18 | 2.5124474413578192e-17 | 143.25261874829306 | yes |
| `eta_npmu` / `eta_mu` (MeV) | -4.7906772975623889e-07 | -4.7906773046561003e-07 | 7.093711409632817e-16 | 4.790777327717226e-18 | 4.790777327717226e-17 | 148.07015489097935 | yes |

All three components exceed the predeclared `10F` threshold. Therefore:

> **CHECKPOINT SEGMENTATION EFFECT CONFIRMED AT RESOLVABLE LEVEL on this
> bounded window.**

## Secondary comparison with the failed hierarchy

All `d_BR` and `d_RU` denominators at index 240 were resolvable. These ratios
are descriptive only:

| Component | `d_BR` | `d_RU` | `Delta_seg/d_BR` | `Delta_seg/d_RU` |
| --- | ---: | ---: | ---: | ---: |
| `x_state` | 3.778579504842838e-09 | 3.3499257257041393e-09 | 0.05230331759598008 | 0.05899600769862285 |
| `eta_e` | 9.327498968956214e-15 | 7.396224462184362e-15 | 0.03858640742173483 | 0.048661946008012415 |
| `eta_mu` | 1.7253838577663207e-14 | 1.4238354804504917e-14 | 0.041113815790628325 | 0.04982114511845439 |

The segmentation effect is decisively above the resolution floor, but at this
checkpoint is approximately 3.9% to 5.9% of the historical adjacent-tier
differences. The experiment establishes a resolvable architecture effect; it
does not prove that segmentation alone accounts for the full failed B/R/U
hierarchy.

## Step-history comparison

| Statistic | Arm S | Arm E |
| --- | ---: | ---: |
| Accepted steps | 289 | 232 |
| Rejected steps | 31 | 60 |
| Minimum step (s) | 1 | 1 |
| 10th percentile step (s) | 206325523.0661621 | 195142011.38098997 |
| 25th percentile step (s) | 1926123046.875 | 1915657007.4191284 |
| Median step (s) | 1926123046.875 | 2212334854.861145 |
| 75th percentile step (s) | 1926123046.875 | 2526398054.752228 |
| 90th percentile step (s) | 1926123046.875 | 2788408529.1384125 |
| Maximum step (s) | 1926123046.875 | 2883572725.0536194 |
| Final accepted step (s) | 1926123046.875 | 844795280.7041626 |
| Exactly 1-s accepted steps | 1 | 1 |
| Accepted steps ending at a requested ceiling | 240 | 1 |

For Arm S, 222/240 intervals used one accepted step, 11 used two, and 7 used
three or more. All 240 intervals necessarily ended an accepted step at their
requested ceiling. Arm E had one interval and 232 adaptive steps; only its
final accepted step ended at a requested ceiling. Arm E used 57 fewer accepted
steps but 29 more rejected trials. The 1-s step in each arm was the common
initial `h`, imposed once at construction of the arm's integration, not a
late-time collapse.

The current architecture preserves GSL step/control/evolve state and `h`
across Arm S checkpoints, but each call to `gsl_odeiv2_evolve_apply` still
receives the next observation time as a hard `t1` ceiling. The equality of the
Arm S step-summary rows to archived ULTRA and the 240 ceiling-ending steps show
that observation times materially constrain the accepted-step history even
without checkpoint-level solver reinitialization.

## Cstar-knot observation

Both arms crossed the same seven Cstar interpolation-knot indices, 94 through
100, once each. They crossed them at different adaptive step endpoints:

| Knot | Arm S crossing endpoint (s) | Arm E crossing endpoint (s) |
| ---: | ---: | ---: |
| 94 | 21717796199.715218 | 21366289325.392822 |
| 95 | 58880042881.52441 | 58924195795.75158 |
| 96 | 103681156811.36728 | 103822991792.349 |
| 97 | 157942089843.75 | 158525294275.20718 |
| 98 | 225288641840.51474 | 225437606411.4205 |
| 99 | 313958056640.625 | 312922590327.025 |
| 100 | 439156054687.5 | 439323538805.5538 |

No arm crossed an additional knot. The differing crossing steps are secondary
evidence that the two ceiling structures generate different adaptive paths;
no independent cache-interpolation pass/fail criterion is inferred.

## Scratch evidence hashes

All raw artifacts remain under the ignored scratch build root and are not
tracked repository content.

| Artifact | SHA-256 |
| --- | --- |
| Scratch harness source | `f4444792215e5b4d9b7b3a21ad4bde4f362159f1600fb87749638a00c4d4d732` |
| Scratch executable | `2631c5756f4b29201d4422ef6452687ef272a7de6966db16535e67f3023b2dd3` |
| Analysis script | `e56ae96654fcb8992fe6f5281b6cd808f8733e33239d7cc0b4bf9c88cd62f66a` |
| Endpoint list | `43ec23ada72bfa59c4e89672ae2987e9927164db765f05afb7b45c924cc0c67e` |
| Arm S trajectory | `7d714ce4c41894c7f272c8a6b15f87f94fcf602244aa9045f2d7b4253be5cbcd` |
| Arm S checkpoint steps | `41ea5b2821ff4e16fc2cae7e23864de14f9ca98f2f2de8f2735f35eadfac8fed` |
| Arm S internal steps | `6cad4c9e7970dab4c4f3a802a6fa351d1529295d76816984a7315bb562ff59fb` |
| Arm S log | `9ea3064c9d778902882e55ff5680e3a693c0a5f9d3d8e2ee38eb353158f17fa3` |
| Arm S reproduction JSON | `b8cb078589de66f851f58cf59cbe9705a25848838b15968db32e8ec51cd8ef32` |
| Arm E trajectory | `8c9531c87b53d8bd189e50802f3d1dd526db6f9b0d256b120c4d635e1e31a98c` |
| Arm E checkpoint steps | `912b0e6300c745f02d733da91fb1072e927aadcd4616b944cc8d02c8a26b2ecb` |
| Arm E internal steps | `fe9afd8c1ddfc7abefca3f1e57a76c62345be1026d7042f41692650f531742c8` |
| Arm E log | `fe465bec45acd2288de4713379c94cb1f13791982d91a3e6d351288257de4575` |
| Final comparison JSON | `a016861b593e1c66cb2bc9724829fefb71bd5d0bd6d2a0d0c05f3e450b757bf3` |

## Interpretation and bounded next decision

The bounded experiment confirms that using observation checkpoints as hard
GSL integration ceilings changes the endpoint solution by a resolvable amount
on the same platform, with the same RKF45 method, tolerances, initial state,
and final time. It does not identify a defect in GSL RKF45 or justify changing
the governed shared integrator. It also does not establish a platform effect:
the proximate mechanism is platform-independent algorithmic segmentation,
while the last bits and exact magnitude can remain platform-sensitive.

The next action requires separate owner authorization and a numerical-output
architecture decision. The smallest proposed follow-up is:

1. specify, before coding, a Phase-6-only passive observation mechanism that
   does not present output times as owner-integrator `t1` ceilings and does not
   mutate the owner path through observer/cache evaluation;
2. keep production `ScaledRKF45`, Cstar, Phase-5D, the physics card, ULTRA
   tolerances, initial `h`, and the archived 0...240 window unchanged;
3. authorize at most one new bounded P2-source integration to validate that
   mechanism, requiring its owner-path final binary64 state and accepted-step
   history to reproduce the existing Arm E evidence exactly before using any
   passively reconstructed checkpoint samples; and
4. predeclare and independently validate the checkpoint reconstruction/dense
   output error before considering a full BA12R retry.

This follow-up must not begin without explicit owner approval of the passive
sampling/interpolation semantics, its bounded implementation scope, its one-run
ceiling, and its validation criteria. A full P2 rerun or BA12R retry remains
unauthorized.
