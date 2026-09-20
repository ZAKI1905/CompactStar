# Phase-6A-1 controlled-BNV recovery: BA12R

Date: 2026-09-20

Scope: the separately owner-approved BA12R P2 ULTRA convergence witness only

Classification: **BA12R VALIDATION EVIDENCE; NOT BNV CANDIDATE; NOT GOVERNED BASELINE**

## Disposition

**BA12R STATE CONVERGENCE FAIL — RETURN TO OWNER.**

This is disposition **B**.  Gate A, gate B, and gate C each have at least one
state comparison failure.  Ledger convergence, matched source-minus-control
convergence, and the final-decade minimum-step requirement also fail.  The
endpoint-energy, R18/R-a/R-b/R-c, R20, frozen-validity, currentness, process,
grid, row-count, and run-card checks pass.  No result was retried or retuned.

The historical BA12 result remains immutable and remains **FAIL**.  BA12R is a
new qualification result and does not replace or reinterpret BA12.

No candidate was classified or created.  No physical BNV rate/model was
selected.  No merge or canonical integration was performed.

## Entry authentication

- Canonical governance authority and canonical `master`:
  `bd697ffdc474863d7a39f42e17ad8e8dbf105e5d`.
- R0-R3 branch: `physics/phase6a1-bnv-recovery-r0-r3`.
- R0-R3 local, upstream, and live SHA at entry:
  `d2822656aae8e264f4985f0d2f7228b90320a577`.
- Entry was clean in both the R0-R3 and canonical worktrees.
- Historical failed implementation SHA
  `e56e6e50040dbcd9dcbee1acecf58843f3dddf1c` is not an ancestor of R0-R3.
- BA12R branch: `physics/phase6a1-bnv-recovery-ba12r`.
- BA12R worktree:
  `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase6a1-bnv-recovery-ba12r`.
- BA12R branch was created from exact R0-R3 SHA
  `d2822656aae8e264f4985f0d2f7228b90320a577`.

The recovery preflight, acceptance record, implementation record, ADR-0015,
and ADR-0016 were read before execution.

## Platform and toolchain authentication

BA12R ran on the same local Mac numerical platform/toolchain used for the
historical BASELINE/REFINED campaign.  It was not moved to EKU.

- OS: macOS 26.6.2, build 25G83; Darwin 25.6.0.
- Architecture: arm64.
- C++ compiler: `/usr/bin/clang++`.
- Compiler version: Apple clang 21.0.0 (`clang-2100.3.34.2`).
- Compiler target: `arm64-apple-darwin25.6.0`.
- GSL version: 2.7.1; linked as `/opt/local/lib/libgsl.27.dylib` and
  `/opt/local/lib/libgslcblas.0.dylib`.
- CMake version: 4.2.1.
- Build type: Debug.
- Effective C++ flags for both the historical trajectory executable and the
  BA12R executable: `-g -std=c++17 -arch arm64 -pthread -Xclang -fopenmp`.
- Python headers/runtime used by both builds: Miniforge Python 3.12.
- BA12R executable SHA-256:
  `83d05403828cd589cbb7b24dd2317345ca6060b0a1a33ad1643a848efdd55f45`.

The first CMake configure attempt stopped before compilation because the
system Python 3.14 environment did not provide NumPy.  The same fresh build
root was configured with the historical Miniforge Python 3.12 interpreter.
This did not change any scientific setting.  No material platform or
toolchain mismatch was found.

## Historical evidence authentication

The historical evidence was read in place from the immutable failed-campaign
scratch root.  It was not regenerated and was not copied into a tracked path.
The exact historical 32-file aggregate command was reproduced, including its
historical relative path labels.

- Historical raw-evidence aggregate SHA-256:
  `81314ddceba77eaaf29f493928560e93f3a6b45d59f302c11fbec0e799da0945`
  — authenticated.
- P2 BASELINE source SHA-256:
  `f331b5bbfd6bdfb0f8e93095052cf1edc712219180370d23383bed1c18b31093`.
- P2 REFINED source SHA-256:
  `1889336755a87377336e1b688e34032dab428304861bd2a63e142e27555b7a3a`.
- P2 BASELINE matched control SHA-256:
  `be807f432c81b7f62594005b593bdd5bcfd301ddfd9ed2837178831fc2dbb9b8`.
- P2 REFINED matched control SHA-256:
  `be807f432c81b7f62594005b593bdd5bcfd301ddfd9ed2837178831fc2dbb9b8`.
- Entry manifest SHA-256:
  `5cfbf4b1d623b58fa2a94101631dce20c602950eec8eb514dc1fcd3c7a04620d`.
- Historical pretrajectory authority SHA-256:
  `1d4e77780fde474f24782879d6ac43fb5bd6d01f3919c359437c0c5f01d327a6`.
- Qualification certificate SHA-256:
  `7fc892b50bfbcdd2c5d963e0ff866777f3628f40fc282e1ce5a0b0364200a453`.
- Frozen certificate SHA-256:
  `9217e278eb0a4b06b193a2a02f4d59f285d74e4033427b87ccd51181c3f7beab`.
- Frozen coefficients SHA-256:
  `3efa060d99a7a92266679aadb31885abbed9c940404247f1d3009cafa94832d6`.

Historical evidence authentication: **PASS**.

## Immutable configurations

The three tiers were used without modification:

| Tier | `rtol` | `atol(x, eta_e, eta_mu)` |
| --- | ---: | --- |
| BASELINE | `1e-7` | `(1e-12, 1e-18, 1e-18)` |
| REFINED | `1e-9` | `(1e-14, 1e-20, 1e-20)` |
| ULTRA | `1e-11` | `(1e-16, 1e-22, 1e-22)` |

Exactly two new trajectories were executed, concurrently as separate
processes with distinct scratch roots, output files, logs, and retained exit
statuses:

1. `CPL-P2-LINEAR-QSS-v1`, ULTRA, controlled source ON.
2. Its exact matched no-BNV control, ULTRA, source OFF.

For the source, `Bdot/B0 = -1.0e-12 yr^-1` and serialized
`Bdot = -2.4136520263641375e37 count/s`.  The control serialized `Bdot = 0`.
Both used duration `5.0e5 yr`, 8193 output samples, initial `Tinf = 1e8 K`,
initial `eta_npe = eta_npmu = 0`, spin OFF, P2 partition for the source,
`Me/Mmu` ON, and `De/Dmu` OFF.  The source final
`DeltaB/B0 = -4.9999999996026683e-7`; the control remained exactly zero.

No run-card, tolerance, drive, duration, grid, or boundary value changed.

## Execution evidence

| Quantity | ULTRA source | ULTRA matched control |
| --- | ---: | ---: |
| process exit status | 0 | 0 |
| wall-clock time | 343.94 s | 342.95 s |
| CPU user time | 336.64 s | 335.74 s |
| CPU system time | 2.99 s | 2.94 s |
| accepted steps | 8355 | 8329 |
| rejected steps | 33 | 12 |
| overall rejection fraction | 0.0039341917024320458 | 0.0014386764176957199 |
| RHS evaluations | 58692 | 58384 |
| output rows | 8193 | 8193 |
| overall minimum accepted step | 1 s | 1 s |
| maximum accepted step | 1926123046.875 s | 1926123046.875 s |

There was no GSL failure, nonfinite state, currentness failure,
frozen-validity refusal, 100000-step checkpoint refusal, output-grid mismatch,
row-count mismatch, run-card mismatch, or external interruption.  No retry was
performed.

## State convergence

The comparison parsed the max-digits-10 tables back to their exact binary64
values.  Therefore declared serialization uncertainty `Q = 0`.  All 49158
source/control state component/checkpoint comparisons were tested individually.

- Maximum `d_RU / D_R = 32.26041661869904` — **FAIL** (limit 1).
  - Mode: source.
  - Component: `eta_mu_MeV`.
  - Index: 216.
  - Time: `416042578125 s`.
  - `d_BR = 2.4715452183430588e-14`.
  - `d_RU = 1.869212699188061e-14`.
  - `D_R = 5.794136886950843e-16`.
  - `D_U = F = 5.794136886950842e-18`.
- Resolvable comparisons: 2444.
- Maximum resolvable `d_RU / d_BR = 0.9329162771682704` — **FAIL**
  (limit 0.10).
  - Mode: source.
  - Component: `x_state`.
  - Index: 179.
  - Time: `344776025390.625 s`.
  - `d_BR = 5.750012188610043e-9`.
  - `d_RU = 5.36427996467026e-9`.
  - `F = 4.310955145946337e-12`.
- Floor-limited comparisons: 46714.
- Maximum floor-limited `d_RU/(10 F) = 43.54669142406804` — **FAIL**
  (limit 1).
  - Mode: control.
  - Component: `x_state`.
  - Index: 58.
  - Time: `111715136718.75 s`.
  - `d_BR = 0`.
  - `d_RU = 1.0959697385737321e-10`.
  - `F = 2.516769248669016e-13`.

State gates A, B, and C: **FAIL**.

## Ledger convergence

The gross-power scale is exactly the accepted `G_P` definition.  No new scale
was introduced.

- Maximum normalized REFINED/ULTRA difference:
  `7.611500166686796e-9` — **FAIL** (limit `1e-9`).
  - Source `Pnet`, index 168, time `323588671875 s`.
  - `d_BR = 3.799778380406335e26 erg/s`.
  - `d_RU = 3.544875091136552e26 erg/s`.
  - `G_P = 4.657262055450493e34 erg/s`.
- Maximum resolvable contraction ratio:
  `0.9329163323161163` — **FAIL** (limit 0.10).
  - Source `Lgamma`, index 216, time `416042578125 s`.
  - `d_BR = 9.401918303328511e25 erg/s`.
  - `d_RU = 8.771203140276998e25 erg/s`.
  - `G_P = 5.055898230134086e34 erg/s`.
  - Floor `= 5.055898230134086e23 erg/s`.
- Maximum floor utilization `d_RU/F = 94.83035685919317` — **FAIL**.
  - Control `Pnet`, index 600, time `1155673828125 s`.
  - `d_BR = 0`.
  - `d_RU = 1.4904755921078965e24 erg/s`.
  - `G_P = 1.5717283383431714e33 erg/s`.
  - Floor `= 1.5717283383431714e22 erg/s`.

Ledger convergence: **FAIL**.

## Endpoint energy and R20

The accepted `N_R20` definition was used without `DeltaEeq` or integrated
`|L_out_fluid|` in the normalizer.

Endpoint REFINED/ULTRA convergence: **PASS**.

- Source normalized differences:
  - `DeltaEeq`: 0.
  - `MeVToErg DeltaEchem`: 0.
  - `DeltaUth`: `1.698152813581614e-14`.
  - Combined endpoint-state propagation: `1.7011061228226257e-14`.
- Control normalized differences:
  - `DeltaEeq`: 0.
  - `MeVToErg DeltaEchem`: 0.
  - `DeltaUth`: `5.203901206568989e-12`.
  - Combined endpoint-state propagation: `5.203901206568989e-12`.

R20 results:

| Evidence | `|R20|/N_R20` | luminosity quadrature / N | thermal quadrature / N |
| --- | ---: | ---: | ---: |
| source REFINED | `1.3774406615827498e-8` | `3.483646054360719e-8` | `6.798480907221146e-10` |
| source ULTRA | `1.414393758909499e-8` | `3.4835964243017257e-8` | `6.79848853911549e-10` |
| control REFINED | `8.728196366555655e-8` | `2.7000514918058557e-7` | `2.9757566776308896e-9` |
| control ULTRA | `8.783182080952754e-8` | `2.700050231581484e-7` | `2.9757766029900346e-9` |

All R20 residual and quadrature criteria pass.  All endpoint-state propagation
ratios are below `5e-5`.

## R18, R-a/R-b/R-c, frozen validity, and currentness

- Source REFINED maximum R18 budget utilization:
  `0.01268636912376294` at index 784, time `1510080468750 s`.
- Source ULTRA maximum R18 budget utilization:
  `0.012686297455034145` at the same checkpoint.
- Control REFINED and ULTRA maximum R18 residual: exactly zero.
- All REFINED and ULTRA R-a/R-b and R-b/R-c residuals: exactly zero.
- Maximum frozen utilization across applicable REFINED/ULTRA source/control
  evidence: `0.27037946636901355` — **PASS**.
- Maximum `|DeltaB|/B0 = 4.999999999602668e-7` — **PASS**.
- Every serialized checkpoint has `valid_through_sample = 1`; full and cheap
  currentness checks did not refuse execution.

R18/R-a/R-b/R-c and frozen/currentness requirements: **PASS**.

## Solver trend

The final temporal decade is `[t_final/10, t_final]`.  Fixed output spacing is
`1926123046.875 s`, so the predeclared minimum-step floor is
`1926.123046875 s`.

| Quantity | ULTRA source | ULTRA control |
| --- | ---: | ---: |
| final-decade rejection fraction | 0 | `0.0004003736821032964` |
| first-half median steps/interval | 1 | 1 |
| last-half median steps/interval | 1 | 1 |
| first-half 95th percentile | 1 | 1 |
| last-half 95th percentile | 1 | 1 |
| final-decade minimum accepted step | 1 s | 1 s |

The rejection and count-trend limits pass for both processes.  The minimum
accepted step fails because `1 s < 1926.123046875 s` for both.  These 1 s steps
are the unchanged integrator initialization at each existing 999-checkpoint
batch boundary; they are nevertheless actual accepted steps and were not
excluded from the predeclared gate.

- Source solver trend: **FAIL**.
- Control solver trend: **FAIL**.

No scientific multithreading was added inside one evolution.

## Matched source-minus-control convergence

REFINED/ULTRA source-minus-control convergence: **FAIL**.

- Maximum stability utilization: `21.923678867554578` for `Tinf`, index 12,
  time `23113476562.5 s`.
- Maximum resolvable contraction ratio: `4.101136780587133` for
  `Tsurface_inf`, index 708, time `1363695117187.5 s`.
- Maximum floor utilization: `41.03343563659818` for `Tsurface_inf`, index
  709, time `1365621240234.375 s`.
- No heating/cooling classification was made.

## Historical BA12

Historical BA12 remains **FAIL** with immutable maximum scaled difference
`1.7255917120989046` at `23113476562.5 s` (`732.421875 yr`), where baseline
`x = 0.051132258160285306` and REFINED
`x = 0.05113226698535251`.

## Output and compact evidence hashes

Raw trajectories and logs remain under the ignored BA12R build/scratch root;
they are not committed.

- ULTRA source trajectory SHA-256:
  `584799310627541259cc1164c8a18f8950afb8dddd7f0f14a1f0b3e9dc557b52`.
- ULTRA source step evidence SHA-256:
  `80dff6f0620c156d942a16ab1f42a724ced156d1216fe44dab17bac8163ba085`.
- ULTRA control trajectory SHA-256:
  `2a742d0db6860e8606b94d3c7e0b9756d5656fdc3abe1a87452608c6ef4c55bb`.
- ULTRA control step evidence SHA-256:
  `69fa548c6d48983de14d01c7b48e68faa837887253f5ed111bdb21a41404f508`.
- Source execution log SHA-256:
  `b2487827cfe91604abaacc3e457a8bbb82ad7bffc726a665e171926841838a04`.
- Control execution log SHA-256:
  `fd541f744ca0bc716332c0546fc695c066060d1f353ca37c06fd04de0c45f918`.
- Compact comparison JSON SHA-256:
  `ca122ad9ad92911fa604ae2ec1ae5208c107a4b6e3dd96ed4aca7c9187e4cce4`.

## Scope attestation and next action

- Other ULTRA cards run: **NO**.
- New BASELINE or REFINED trajectories run: **NO**.
- Reaction-free trajectory run: **NO**.
- Full test suite run: **NO**.
- Unfinished BA11 completion beyond BA12R: **NO**.
- BA16, trajectory BA17, BA14 finalization, or full BA15 run: **NO**.
- Candidate producer or candidate artifact run: **NO**.
- Candidate created or ratified: **NO**.
- Canonical merge performed: **NO**.
- Physical BNV rate/model selected: **NO**.

Exact recommended next action: return this BA12R state-convergence failure and
the nonblocking ledger, matched-difference, and solver-trend findings to the
owner.  Do not begin the bounded unfinished-validation completion task.  Any
new numerical investigation, implementation change, retry, or altered gate
requires a separately owner-approved task.
