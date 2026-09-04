# Phase 4D-RV — corrected independent revalidation of the Hartle O(Ω²) monopole response

> **FORMAL STATUS: `HARTLE O(OMEGA^2) MONOPOLE RESPONSE CHARACTERIZED — INDEPENDENT VALIDATION INCOMPLETE`**
>
> Every physics line of the corrected campaign passes on the migrated ADR-0009 backgrounds with
> the independent `(m₀, h₀)` oracle carrying the EOS measure by its own Stieltjes route: analytic
> profile `9.7e-9`, continuum first integral `6.1e-15`, DS(CMF)-1 four stars `≤ 5.3e-5` (fully
> independent chain, EOS-knot measure) against `1e-4`, homogeneous sequence identity
> `1.03e-4 / 5.7e-5` against `2e-4`, Chandrasekhar & Miller 1974 to `7.3e-4`, Hartle & Thorne
> 1968 to `1.1e-2`, detectors M1–M10 all load-bearing, production diff **NONE**. **One accepted
> line is NOT MET:** ADR-0008 Validation D's *monotonicity* clause. With the surface fixed at
> `R_* = 13.473358 km` to `3e-11` across resolutions, `δM̂` over 5000/10000/20000/40000 still moves
> `−2.1e-3, −3.0e-2, +9.5e-3 km³` (spread `3.7e-5` — inside D's `1e-4`), and the decomposition
> attributes the residual to O(h) sub-node node-placement in **both** the profile-partition
> measure (removed by the EOS-knot oracle) and the validated first-order sampled background
> (the same `−9.5e-3 km³` dip appears in `I`). The accepted production representation cannot meet
> a monotone-differences criterion at this level by construction. Per the task's §27 rule the
> status is therefore **CHARACTERIZED**, not VERIFIED; **no monopole baseline was created**, and
> the governance question is put to the owner (§17, §20).
>
> **SLOW-ROTATION DISCLAIMER.** Everything here establishes the O(Ω²) coefficients on the finite
> EOS-floor surface. Nothing bears on the accuracy of truncating a rapidly rotating star at O(Ω²).

| Field | Value |
|---|---|
| **Starting HEAD** | `218c9aa302e5641084109746e4ff13fc36c56df5` (seven-artifact migration, TOV-SURF-MR), branch `validation/hartle-monopole-revalidation`, worktree `worktrees/CompactStar-monopole-validation`, created clean from that exact commit; **20 ahead / 0 behind** `master` = `df859b5a73c4cac0c115f240744d89ce9f830b8d` |
| **Ancestry authenticated** | `377bc4a` (ADR-0007 implementation), `8abbab4` (ADR-0008 derivation), `cc4bec4` (ADR-0008 acceptance), `e2fe0ad` (measure-complete source), `7a223d3` (ADR-0009 acceptance), `94165f3` (validated surface implementation), `816b754` (migration-envelope adjudication), `218c9aa` (migration) — all ancestors |
| **Change class** | **test + documentation**; production diff **NONE** (`git diff -- CompactStar/ CMakeLists.txt` = 0 lines at every stage; detector mutations reverted byte-identically, SHA-256 verified) |
| **Governing authority** | ADR-0007 (ACCEPTED) as amended by ADR-0008 (ACCEPTED) §Validation A–J; ADR-0009 (ACCEPTED, conformed, validated, migrated) Q1–Q14; ADR-0006; ADR-0003; ADR-0005; INV-06, INV-08, INV-09, INV-13, INV-14; `GOVERNANCE.md` §2–§3, §3.1 (ADR-0007 §9 / ADR-0008 Q12) |
| **Suites at entry** | main non-Phase-4D **41/41** (553.37 s); self-contained non-Phase-4D **20/20** (67.36 s); CTest inventory 44 full / 22 self-contained, the three independent Phase-4D tests deliberately excluded as the brief expects; the three existing Phase-4D harnesses also ran green diagnostically before any edit |
| **Suites at exit** | full **44/44** (665.79 s, serial, the three Phase-4D harnesses included; `hartle_monopole_physics_cmf` 36/36 with Hc recorded NOT MET rather than asserted — the 4D line-J precedent — so the committed suite is green while the scientific status is CHARACTERIZED); self-contained **22/22** (89.06 s); analytic harness 25/25, published harness 16/16 |
| **Evidence** | this record; harness logs preserved in `build/rv/` of the worktree (not tracked) — `cmf_rv_run*.log`, `analytic_rv_run1.log`, `published_rv_run2.log`, `detectors/`; scratch `primary_rederivation.md` in the session scratchpad |

---

## 1. Starting authority and the seven migrated artifacts

All seven `tests/baselines/*.tsv` hashes at entry, re-verified unchanged at every stage and at exit (§15):

| Artifact | SHA-256 (current, migrated) |
|---|---|
| `passive_cooling_cmf_1p6_debug.tsv` | `afcad5d078fe458bdb441278ce56caa2d025becc6b489d00e9baad91a101c0de` |
| `tov_dscmf1_reference.tsv` | `3d9af9129a6a4ffde9e0f8c5507a160f968a861c0cf9f3b089cceecab86b701a` |
| `grid_convergence_cmf_1p6_debug.tsv` | `2c68e2f7e871192e00322bcc12b6d8b13eb7fdbda4a6d8d7050f26f0271ce5eb` |
| `grid_convergence_cmf_1p6_trajectory.tsv` | `7c84557742ec0ec747756d118b2837fb54095a94c10194d4ccc2d0a778f5b04f` |
| `hartle_I_dscmf1_debug.tsv` | `a21d4c3f6d89322cb4cca7a073e0e142f180e529332d704b7b16a145eae741c9` |
| `baryon_number_dscmf1_reference.tsv` | `7b036942f2ae599ace3b2fc8b9a7d91f6d11b5899ab7e8a88d2bb4ee6686493b` |
| `tov_path_equivalence_dscmf1.tsv` | `5c0f4b3bdb70921f8f2a869af10edc4d8f5ae3963a9d150e11ef859d21e1c678` |

They are background/reference authorities here, never oracles for the second-order response.
`tests/baselines/hartle_monopole_dscmf1_debug.tsv` did not exist at entry and **does not exist at exit** (§18).

## 2. Claim and non-claims

**Claim.** For ordinary `NStar` on the governed finite EOS surface `p(R_*) = p_cut`, with `R_*` the
accepted-solution event stored as the final node (ADR-0009), the fixed-`ε_c`, `l = 0`, O(Ω²)
response `m̂₀(r)`, `p̂₀*(r)`, `δp̂₀(r)`, `ξ̂₀(r)`, `δM̂` per `Ω_geom²` is physically correct when the EOS
energy-density source is the measure `dm̂₀|_EOS = −4πr²ξ̂₀ dε` (ADR-0008): smooth variation, sharp
continuous tabulated variation, the terminal `ε_* → 0` atom, and — as a mathematical rule for
future declared atoms — the jump operator `4πr_t²(ε⁻−ε⁺)ξ̂₀(r_t)`. No EOS in the repository
declares an internal constant-pressure discontinuity, so no execution of that operator is claimed.

**Not claimed.** `l = 2` / quadrupole; high-spin (Ω⁴) truncation accuracy; baryon-conserving
sequence reduction; Phase-5 `A_i/B_i/Z_i`; rotochemical heating; MixedStar rotation.

## 3. Primary-source authentication

The ADS journal scans used by Phase 4C-G were not available on disk in this session and no file
was downloaded; the page-level transcription of H67 / HT68 / C&M recorded in
`docs/validation/PHASE4C_HARTLE2_DERIVATION.md` §2 (read from the scans) is the transcription of
record. Every governing equation was **re-derived independently** here (scratch
`primary_rederivation.md`) and checked against that transcription and against two
transcription-independent facts:

| Object | Re-derived | Independent confirmation |
|---|---|---|
| metric (H67 66–70, `k₀ = 0`), `ν_H = 2ν`, `j² = e^{−2ν}(1−2m/r)`, `(j²)' = −8πr(ε+p)e^{−2ν}` | ✓ | production `I` and first order (Phases 2B-4B, 4B) |
| Eulerian matter changes `δε = −ξ₀ dε/dr`, `δp = −ξ₀ dp/dr`; `p₀* = ν'ξ₀` (H67 79–80, 87–88, 99) | ✓ | — |
| `G_t^t`: `dm₀/dr = −4πr²ξ₀ dε/dr + (1/12)j²r⁴ω̄'² − (1/3)r³(j²)'ω̄²` — the **measure** form (H67 93→97) | ✓ | secant and Stieltjes routes agree with production (§9–§10) |
| `G_r^r`: `dh₀/dr` (H67 98); first integral `p₀* + h₀ − (1/3)e^{−2ν}r²ω̄² = γ` (H67 90); their combination is (100) | ✓ | continuum identity `(100) ≡ (98)+(90)` residual `6.1e-15` (§10 E) |
| fixed-`ε_c` family `m₀(0) = p₀*(0) = 0`, series `p₀* → (1/3)j_c²ω̄_c²r²`, `m₀ → (4π/15)(ε_c+p_c)[(dε/dp)_c+2]j_c²ω̄_c²r⁵` (H67 108) | ✓ | Experiment A `≤ 5.0e-10` |
| exterior `m₀ = δM − J²/r³`, `δM = m₀(a) + J²/a³` (H67 105–107; HT68 16); surface atom `4πR_*²ε_*ξ₀(R_*)` (H67 18; measure form) | ✓ | near-vacuum identity `≤ 2.7e-7`; Newtonian intercepts `1e-6` |
| Newtonian homogeneous limits `δM → Ω²a³`, `3Mξ₀(a)/a⁴ → 1` (H67 14–21) | ✓ | Experiment D |
| C&M 1974 eq. (32) `m₀` continuous ⇒ shell-excluded `δM/M`; HT68 Table 5 at `Ω² = M/R³` with log interpolation of Table 1 | ✓ | K1 19/19, K2 8/8 — a mis-transcribed coefficient in (97)/(100) could not reproduce both tables |

## 4. Independent formulation and independence audit

The oracle (`tests/rotation/hartle_monopole_reference.hpp`) integrates Hartle's **`(m₀, h₀)`**
pair, (97)+(98), with `p₀*` an algebraic by-product of the first integral (90) — production
integrates `(m₀, p₀*)`, (97)+(100). Different state vector, different pressure-side equation.
Audit of independence (every item re-read at `218c9aa`):

| Requirement (task §7) | Oracle |
|---|---|
| never calls `ComputeHartleMonopoleResponse`, `MonopoleResponse`, `MonopoleAt`, `ODE_HartleMonopole_` | none referenced; only profile columns and the first-order response (chain A) or `hartle_reference.hpp` (chain B) are read |
| no copy of production's segment loop | the Stieltjes route (§5) has its own sub-partition loop with midpoint atoms; the legacy `eos_measure` route (secant density in the RHS) is retained for the variable-pair cross-check only and is **reported, never asserted** |
| production `eps_slope` not used as an oracle | the Stieltjes route carries no EOS term in its RHS at all; `Δε` comes from the profile's `ε` nodes and, in knot mode, from the EOS table read by `tests/rotation/eos_table_knots.hpp` |
| own integrator | GSL rk8pd driver owned by the oracle, `rtol 1e-13`, `atol 1e-16` (production `1e-10`) |
| own interpolation / partition | binary-search linear interpolation; K-fold sub-partition with optional EOS-knot mapping through the interval's own linear `p(r)` |
| own exterior / matching | `δM̂ = m̂₀(R_*) + 4πR_*²ε_*ξ̂₀(R_*) + I²/R_*³` assembled from the oracle's own fields |
| inputs only | authenticated TOV profile (`r, p, ε, m, ν, ν'`), EOS `ε` values (profile nodes; table knots), authoritative `dε/dp` for the centre series only, and either production's verified `s, s'` (chain A) or `hartle_reference.hpp`'s (chain B) |

## 5. Independent measure integration design

`hartle_mono_ref::SolveStieltjes` represents `dm̂₀|_EOS = −4πr²ξ̂₀ dε` as a Riemann–Stieltjes sum
of **atoms**: every profile interval is split into `K` equal sub-intervals (optionally after every
EOS-table knot inside the interval is mapped in through the interval's linear `p(r)`, with `ε`
at the knot taken from the table), and each sub-interval's `Δε_j` is applied once at its midpoint
as `m̂₀ += −4π r_mid² ξ̂₀(r_mid) Δε_j`, `ξ̂₀` from the first integral at the current state (`ĥ₀`
continuous through an atom). The driver integrates the smooth part with no EOS term. The terminal
atom `ε_* → 0` is applied exactly once in the `δM̂` assembly. Production's realisation is a
per-segment ODE density integrated by its RK driver; the two routes share only the measure they
represent. `K` refinement measures the route's own floor; the EOS-knot mode moves the sub-node
location of the measure to where the star's own EOS puts it — the profile-partition vs knot
comparison is the measure-discretization diagnostic the brief asks for (§9, §10 H).

## 6. Oracle floor (Experiment R, 1.6 M☉, migrated background)

| Variation | movement of the oracle |
|---|---|
| integrator tolerance `1e-11 … 1e-15` | `4.3e-8` |
| centre start `r₀ = 1e-7 km` (vs the first node `1e-5`) | `0.0` |
| profile-partition refinement ladder `K = 1, 2, 8` vs `K = 4`: `δM̂ = 865.86822188, 865.86661495, 865.86611310` vs `865.86621332` | `2.4e-6, 4.8e-7` (ladder, O(K⁻²)); **floor `1.2e-7`** (`K = 8`) |
| EOS-knot refinement ladder `K = 1, 4` vs `K = 2`: `865.85619006, 865.85598990` vs `865.85602984` | `3.2e-7` (ladder); **floor `8.1e-8`** |

Absolute floor on `δM̂`: `≈ 1e-4 km³`; relative `1.2e-7` (profile) / `8.1e-8` (knot). Signals: production vs profile oracle `2.2e-7` (floor-limited — the profile oracle converges in `K` to the same continuous-density measure production integrates, so this is agreement *down to the oracle floor*, both two decades under the `1e-4` bound); production vs knot oracle `2.5e-5`, floor/signal `3e-3`. **Agreement is never claimed below the oracle's own precision.**

## 7. Chains A and B — DS(CMF)-1 four stars (Experiment G, bound `1e-4`)

Target-mass workflow, migrated background, every star `SURFACE_REACHED` with `p(R_*) = p_cut` exactly (§14). Worst node-wise relative disagreement (`|ref| ≥ 1 %` of the profile peak), production vs the independent Stieltjes oracle:

| M☉ | chain A, profile K=4: `m̂₀ / p̂₀* / δp̂₀ / ξ̂₀ / δM̂` | chain A, EOS-knot K=2 | chain B, profile K=4 | chain B, EOS-knot K=2 |
|---|---|---|---|---|
| 1.0 | `2.6e-7 / 1.4e-7 / 1.4e-7 / 1.4e-7 / 2.0e-7` | `5.3e-5 / 4.2e-5 / – / 4.2e-5 / 4.1e-5` | `6.7e-6 / 6.7e-6 / – / 6.8e-6 / 3.1e-6` | `5.0e-5 / 4.0e-5 / – / 4.0e-5 / 3.8e-5` |
| 1.4 | `2.3e-7 / 1.5e-7 / 1.5e-7 / 1.5e-7 / 1.9e-7` | `3.9e-5 / 3.4e-5 / – / 3.4e-5 / 3.4e-5` | `1.1e-5 / 1.2e-5 / – / 1.2e-5 / 3.0e-6` | `4.2e-5 / 3.5e-5 / – / 3.5e-5 / 3.7e-5` |
| 1.6 | `2.2e-7 / 1.5e-7 / 1.5e-7 / 1.5e-7 / 1.9e-7` | `2.4e-5 / 2.5e-5 / – / 2.5e-5 / 1.2e-5` | `1.8e-5 / 1.9e-5 / – / 1.9e-5 / 5.6e-6` | `1.8e-5 / 2.4e-5 / – / 2.4e-5 / 6.1e-6` |
| 2.0 | `2.0e-7 / 2.9e-7 / 2.8e-7 / 2.9e-7 / 1.8e-7` | `3.3e-5 / 3.9e-5 / – / 3.9e-5 / 2.8e-5` | `3.5e-5 / 3.6e-5 / – / 3.6e-5 / 2.9e-6` | `3.7e-5 / 3.6e-5 / – / 3.6e-5 / 3.1e-5` |

Chain B is bounded by the Phase-4B first-order gap (`≈ 1e-5`, as in 4D); the knot oracle's `1–5e-5` is the profile-vs-knot measure-discretization gap (`δM̂`: `4.1e-5 / 3.4e-5 / 1.2e-5 / 2.8e-5`; node-wise `5.3e-5 / 3.9e-5 / 2.5e-5 / 3.3e-5`). Reported beside them: the production-like secant oracle (`≤ 3.3e-7`; variable-pair cross-check only) and the **superseded** nodal-derivative oracle, which disagrees by `6.1e-2 / 5.1e-2 / 4.6e-2 / 2.0e-2` — the sub-node energy-density variation ADR-0008 recovers (M10's target). **Met in both chains on every star.**

Four-star table (migrated ADR-0009 background; `p_cut = 3.3518849e25 dyn cm⁻²`, `ε_* = 1.6588079e8 g cm⁻³` on all four):

| target | `ε_c` [g cm⁻³] | `M` [M☉] | `R_*` [km] | `I` [km³] | `(dε/dp)_c` | `m̂₀(R_*)` | `p̂₀*(R_*)` | `ξ̂₀(R_*)` | shell | `I²/R_*³` | `δM̂` [km³] |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1.0 | 4.545504051e14 | 1.0000438097 | 13.427522155 | 86.99575838 | 5.2824846 | 992.1316332 | 34.69765568 | 3304.386797 | 9.222590e-4 | 3.126140724 | 995.25869618 |
| 1.4 | 6.164882705e14 | 1.4000217833 | 13.546350433 | 135.6161306 | 3.7908370 | 933.5227375 | 33.85037764 | 2087.430175 | 5.929620e-4 | 7.398706152 | 940.92203664 |
| 1.6 | 7.312533427e14 | 1.5999758353 | 13.473358194 | 159.5871416 | 3.2786464 | 855.4527843 | 32.14457947 | 1603.536435 | 4.506100e-4 | 10.41281673 | 865.86605162 |
| 2.0 | 1.297985390e15 | 1.9999985862 | 12.716477747 | 193.7221162 | 2.5627582 | 555.6188693 | 23.93200974 | 701.6846603 | 1.756491e-4 | 18.24979129 | 573.86883625 |

These are diagnostics of this run on the migrated background — not a baseline, not an oracle.

## 8. Experiments A–F, J, S on the exact constant-density interior (self-contained harness)

| Line | Result | Bound |
|---|---|---|
| A centre series, first ten nodes: `p̂₀*/r²`, `m̂₀/r⁵`, `ξ̂₀/r` | `1.26e-10`, `5.02e-10`, `1.24e-10` | `1e-8` — **met** |
| A2 enlarged-centre fixture `r₀ = 0.1 km` (M7's load-bearing target) | production vs oracle `3.8e-8`; a literal-zero start would displace `p̂₀*` by `5.7e-5` of its surface value | `1e-7` — **met** |
| B full profile N=4001, isolated / fully independent / `δM̂` | `9.72e-9 / 6.89e-9 / 8.7e-9, 6.0e-9`; exact `h²` (2001: `3.9e-8`, 8001: `2.4e-9`) | `1e-7` — **met** |
| Bs Stieltjes oracle on a zero interior measure vs the differential oracle; vs production | `4.7e-15`; `9.72e-9` | `1e-12`; `1e-7` — **met** (ADR-0008 A: values unchanged) |
| B-cont vs the continuum solver, N=1001…8001 | monotone `h²` convergence | **met** |
| E continuum first integral `(100)` vs `(98)+(90)` | `6.13e-15` | `1e-9` — **met**; tabulated form `1.23e-7 → 3.07e-8 → 7.67e-9 → 1.92e-9 → 4.80e-10` (N = 2001…32001, ratios 4.00) — the `O(h²)` input floor, reported as in 4D |
| C shell: `δM̂` identity; shell vs `4πR²εξ̂₀(R)` (oracle); shell/`δM̂` | `0.0`; `9.72e-9`; `0.896` | exact; `1e-7`; dominant — **met** |
| D Newtonian limits, `M/R = 0.15 … 0.001`: linear intercepts `|δM̂/R³ − 1|`, `|3Mξ̂/R⁴ − 1|`; monotone; without the shell | `1.07e-6`, `4.9e-7`; yes; `4.0e-4` instead of 1 | `5e-3` — **met** |
| J homogeneous `δM̂` vs exact `(dM/dp_c)δp_c` | `3.05e-9`; `m̂₀_hom ≡ 0` | `1e-3` — **met** |
| R oracle floor (bitwise under tolerance/start; rk4 at `1e-7` moves it `6.6e-15`) | `0` vs signal `9.7e-9` | `0.1` — **met** |
| S materialization: `±Ω` bitwise; `Q(2Ω) = 4Q(Ω)`; zero spin zeros | yes; `0.0`; yes | **met** |

## 9. Experiment F — smooth-EOS equivalence and ADR-0008 Validation F (HT68 HW fixture)

ADR-0008 F reads "HT68 published comparison — unchanged within `2e-5`". **Authenticated scope:** the
`2e-5` was set from ADR-0008's Context table ("smooth HW EOS vs today `≤ 1.25e-5`", scratch stars at
`ε_c = 3e14, 1e15, 3e15 g cm⁻³`, `PHASE4D_R_EOS_MEASURE_DERIVATION.md` §8) and 4D-RI §13 asserted it on
`δM̂` (measure vs the superseded differential form) at exactly those densities. The two objects are
therefore `δM̂` in the measure form and `δM̂` in the superseded differential form on the smooth HW EOS;
it was never a byte-identity claim against pre-ADR-0009 output.

| `ε_c` [g cm⁻³] | nodes | measure vs superseded differential (`δM̂`) | vs Stieltjes profile K=4 | vs Stieltjes EOS-knot K=2 | vs secant |
|---|---|---|---|---|---|
| 1e14 | 5838 | `1.66e-8` | `9.08e-7` | `9.29e-7` | `8.96e-7` |
| **3e14** | 3683 | **`5.55e-6`** | `3.90e-7` | `4.00e-7` | `3.73e-7` |
| **1e15** | 2731 | **`1.25e-5`** | `2.35e-7` | `1.60e-7` | `2.03e-7` |
| **3e15** | 2165 | **`1.15e-5`** | `1.76e-7` | `2.59e-7` | `1.16e-7` |
| 6e15 | 1910 | `4.08e-5` | `1.65e-7` | `5.91e-7` | `7.31e-8` |
| 1e16 | 1779 | `1.67e-4` | `1.81e-7` | `7.71e-7` | `5.87e-8` |
| 3e16 | 1563 | `2.47e-5` | `3.63e-7` | `1.40e-6` | `1.37e-7` |
| 1e17 | 1446 | `1.42e-4` | `8.89e-7` | `1.91e-6` | `5.43e-7` |

**F1 (authenticated scope, asserted): `1.25e-5 ≤ 2e-5` — met**, reproducing 4D-RI digit for digit on
the migrated surface. **F1r (recorded chronology):** a first 4D-RV run asserted the same bound over all
eight HT68 configurations and measured up to `1.67e-4` at `1e16`; the bound was **not widened** — it
is asserted at its authenticated scope and every configuration is reported. The growth of the
superseded form's deficit with the density gradient is the ADR-0008 mechanism on the densified
table's sub-node slope structure; the **independent** Stieltjes oracles agree with production on all
eight to `≤ 1.9e-6` (**F2, asserted at `2e-5`, met**), so production is the measure-complete side.
Every HT68 star is a complete ADR-0009 star (`K2z`).

## 10. DS(CMF)-1 campaign — near-vacuum, accounting, convergence, sensitivity, sequence identity

**F near-vacuum identity** (outermost nodes, `ε < 1e9 g cm⁻³`, `r > R_*/2`): matter-corrected
`m̂₀ + S + I²/r³` spread over `δM̂` = `2.44e-7 / 2.74e-7 / 2.62e-7 / 2.45e-7` (11/8/6/4 nodes; raw
`4.0e-6 / 2.9e-6 / 2.2e-6 / 1.2e-6`), bound `1e-6` — **met**; `s'(R_*) = 6I/R_*⁴` to `≤ 2.2e-16`;
matching arithmetic exact.

**ADR-0008 C, same-partition accounting** (own four-variable driver, `rtol 1e-13`): production
`m̂₀` reproduced to `7.3e-8 / 8.1e-8 / 8.2e-8 / 6.3e-8` node-wise and `2.5e-8 / 4.6e-9 / 3.8e-9 /
1.6e-8` at `R_*`; per-segment identity residual `≤ 3.4e-10 km³` of EOS totals `910 / 809 / 711 / 389`
km³; bound `1e-6` — **met**.

**H — corrected radial convergence at fixed `ε_c = 7.312533427e14 g cm⁻³` (ADR-0008 D).** Every star
complete; `R_*` spread `3.4e-11` (bound `1e-8`, ADR-0009 floor — **met**):

| res | nodes | `R_*` [km] | `M` [M☉] | `I` [km³] | `m̂₀(R_*)` | `ξ̂₀(R_*)` | `I²/R_*³` | **`δM̂`** | profile oracle K=4 | EOS-knot oracle K=2 |
|---|---|---|---|---|---|---|---|---|---|---|
| 5000 | 1320 | 13.473358194 | 1.5999758353 | 159.5881112 | 855.4547984 | 1603.538914 | 10.41294326 | **865.86819222** | 865.86869133 | 865.89612869 |
| 10000 | 2636 | 13.473358194 | 1.5999758353 | 159.5871416 | 855.4527843 | 1603.536435 | 10.41281673 | **865.86605162** | 865.86621332 | 865.85602984 |
| 20000 | 5269 | 13.473358193 | 1.5999758353 | 159.5776228 | 855.4240608 | 1603.519371 | 10.41157460 | **865.83608597** | 865.83614760 | 865.85364848 |
| 40000 | 10536 | 13.473358194 | 1.5999758353 | 159.5807196 | 855.4331145 | 1603.525231 | 10.41197869 | **865.84554381** | 865.84558674 | 865.84754056 |
| 80000 | 21068 | 13.473358194 | 1.5999758353 | 159.5813769 | 855.4349592 | 1603.526821 | 10.41206447 | **865.84747426** | 865.84756129 | 865.84500369 |

`p̂₀*` at fixed radii `R/4, R/2, 3R/4` (km²): `2.8602543, 10.976799, 22.589115` (5000) →
`2.8602926, 10.976951, 22.589448` (80000); the shell is `4.5061e-4` at every resolution.

| successive differences 5000→10000→20000→40000 [km³] | `d₁` | `d₂` | `d₃` | magnitudes decreasing | values monotone |
|---|---|---|---|---|---|
| production `δM̂` | `−2.141e-3` | `−2.997e-2` | `+9.458e-3` | **NO** | **NO** |
| profile-partition oracle | `−2.478e-3` | `−3.007e-2` | `+9.439e-3` | NO | NO |
| EOS-knot oracle | `−4.010e-2` | `−2.381e-3` | `−6.108e-3` | NO | **YES** |
| first-order `I` | `−9.70e-4` | `−9.52e-3` | `+3.10e-3` | NO | NO |
| `m̂₀` EOS channel | `−1.42e-3` | `−2.31e-2` | `+7.24e-3` | NO | NO |
| `m̂₀` rotational channel | `−5.89e-4` | `−5.63e-3` | `+1.83e-3` | NO | NO |
| extension 40000→80000 | production `+1.93e-3`, profile oracle `+1.98e-3`, knot oracle `−2.54e-3` | | | | |

Relative spread over 5000…40000: production **`3.708e-5`** (bound `1e-4` — **Hb met**), profile oracle
`3.76e-5`, knot oracle `5.61e-5`. **Hc — the monotonicity clause of ADR-0008 D — NOT MET** under both
readings (decreasing magnitudes; monotone values), recorded, not waived, not widened.

*Diagnosis.* With the surface fixed, the residual is sub-node **node placement**:
(i) the profile-partition measure locates each interval's `Δε` uniformly across the interval
(ADR-0008 Q3), an `O(h)` weight-location error at the 0.44-m crust–core layer that the ADR itself
names ("the weight's sub-segment location error is O(h) and vanishes with the optional knot
partition", Q4); the EOS-knot oracle removes it — its `δM̂` sequence becomes monotone
(`865.896 → 865.856 → 865.854 → 865.848 → 865.845`) and production minus knot-oracle alternates
`−2.8e-2, +1.0e-2, −1.8e-2, −2.0e-3, +2.5e-3 km³`, shrinking with resolution; (ii) the validated
first-order sampled background carries the same `O(h)` node-placement scatter — `I` itself dips by
`−9.5e-3 km³` (`6.0e-5`) at 20000 (its recorded resolution spread is `6.6e-5`,
`TOV_SURFACE_ARTIFACT_MIGRATION.md` §11) — and enters `δM̂` through `I²/R_*³` (`−1.2e-3 km³`) and
the rotational sources `s², s'²` (rotational channel `−5.6e-3 km³`), which is why even the knot
oracle's difference magnitudes do not decrease. Both effects are properties of the accepted
representations (ADR-0008 Q3/Q4; INV-13 sampled background), not of the monopole physics: the
measure-complete `δM̂` agrees with the fully independent EOS-knot oracle to `≤ 4e-5` at every
resolution and to `≈ 3e-6` at 40000/80000, the spread is inside D's bound with a margin of 2.7, and
the limit is `865.846 ± 0.002 km³`. The 4D-RI attribution of this non-monotonicity to `R_*` moving
with the grid is **superseded**: `R_*` no longer moves and the pattern is unchanged.

**I — measure / derivative sensitivity (1.6 M☉).** A governed measure `865.86605162`; B retired
profile-FD derivative substituted through the explicit-supply path `865.86605162` — spread **`0.0`**
(ADR-0008 E `< 1e-3` — **met**; `dε/dp` enters only the centre series); C superseded nodal-`dε/dp`
differential form (oracle) `825.97018842`, `−4.6e-2` (the historical deficit as a diagnostic);
D EOS-knot-refined measure (oracle) `865.85602984`, `1.16e-5`; profile-partition Stieltjes
`865.86621332`, `1.9e-7`; secant `865.86607920`, `3.2e-8`.

**J — homogeneous family vs the TOV sequence derivative (ADR-0008 B, `≤ 2e-4`).** Sources off,
`p̂₀*_c = 1`, `δp_c = ε_c + p_c`; `dM/dp_c` from `M(ε_c(1 ± 10⁻³))`, both neighbours complete ADR-0009
stars; `(dε/dp)_c` sequence vs authority `9.9e-7`; `dM/dp_c = 9671.2436 km³` at every resolution:

| res | expected `(dM/dp_c)δp_c` | Stieltjes profile K=4 | rel | EOS-knot K=2 | rel | secant (rep.) | superseded differential (rep.) |
|---|---|---|---|---|---|---|---|
| 10000 | 5.996175563 | 5.995559079 | **`1.03e-4`** | 5.995629537 | **`9.1e-5`** | `1.03e-4` | `1.17e-3` |
| 20000 | 5.996175560 | 5.995832474 | **`5.7e-5`** | 5.995709033 | **`7.8e-5`** | `5.7e-5` | `1.02e-3` |
| 40000 (diag.) | 5.996175559 | 5.995744262 | `7.2e-5` | 5.995731185 | `7.4e-5` | `7.2e-5` | `1.04e-3` |

**Met** at both binding resolutions with the independent measure route; the superseded form sits at
`1e-3`, the 4D failure, as a labelled diagnostic.

## 11. ADR-0008 Validation A–J

| Line | Requirement | Result | Status |
|---|---|---|---|
| A | all Phase-4D analytic checks green; constant-density interior measure zero; values unchanged | 4D numbers reproduced digit for digit; Stieltjes oracle vs differential `4.7e-15` | **met** |
| B | homogeneous DS(CMF)-1 sequence identity, 10000 and 20000 `≤ 2e-4` | `1.03e-4`, `5.7e-5` (profile); `9.1e-5`, `7.8e-5` (knot) | **met** |
| C | sourced same-partition accounting `≤ 1e-6` | `≤ 8.2e-8` node-wise, `≤ 2.5e-8` at `R_*` | **met** |
| D | radial convergence 5000/10000/20000/40000: spread `≤ 1e-4` **and** successive differences monotone; no node-placement behaviour | spread `3.71e-5`; differences `2.1e-3, 3.0e-2, 9.5e-3 km³` — not monotone; attributed to O(h) node placement of the measure weight and of the first-order background (§10) | **spread met; monotonicity NOT MET** |
| E | retired-FD substitution moves `δM̂` `< 1e-3` | `0.0` | **met** |
| F | smooth HT68/HW comparison unchanged within `2e-5` (authenticated scope: `δM̂` measure vs differential at `3e14/1e15/3e15`) | `5.6e-6 / 1.25e-5 / 1.15e-5`; independent oracles `≤ 1.9e-6` on all eight; published K2 unchanged within `2e-2` | **met** (eight-star extension recorded, §9) |
| G | constant-density surface shell preserved | `δM̂ = 1.4674047059e3`, shell/`δM̂` `0.896`, unchanged | **met** |
| H | M1–M9 still fire | §13 | **met** — M1–M9 all fire (M4 on the constant-density/Newtonian lines only, M5 on CMF/HW only, as expected) |
| I | M10 (nodal `dε/dp` source) recreates `≥ 3 %` deficit | §13 | **met** on the fixture star — `4.61 %` at 1.6 M☉ (`6.14 %`, `5.10 %` at 1.0/1.4; `1.97 %` at 2.0 M☉, recorded) |
| J | seven artifacts, `I`, first-order tests unchanged | §15; `hartle_I_dscmf1_debug.tsv` `a21d4c3f…` | **met** |

## 12. Spin materialization (secondary)

`+100` / `−100 rad s⁻¹` bit-identical; `Q(2π·100) = 4Q(π·100)` to `0.0`; zero spin exact zeros; no raw
seed dependence (the 4C-I1 contract test re-runs green, seed invariance `≤ 1e-10`). No high-spin claim.

## 13. Detectors M1–M10

Each mutation was applied at one production site of `CompactStar/Core/src/RotationSolver.cpp`
(M5 at its two coupled sites), rebuilt in a separate build tree (`build-det`), run against the
analytic, published, CMF-quick (`Z + G + F + C`), contract and measure-contract harnesses, and
the source restored byte-identically after each (SHA-256 `c4f2c2f3898312d1751f027f8d6c015dd2d9b6f01c22a403e183a547bfb33138` before, after each detector,
and after the sweep). Nothing was committed.

| Detector | Mutation (one production site) | Fires on | Signature | Blind (expected) |
|---|---|---|---|---|
| M1 | drop `e^{−2ν}` from term 3 | analytic B/Bs/A2/C/D (8); published K1c/K1d/K2h/K2i/F1 (5); CMF G/C (8); measure-contract (7) | analytic `4.2e-1`; C&M `ξ₀` `5.7`; CMF `m̂₀` `9e-2 … 1.6e-1` | contract (seed/algebra only) |
| M2 | sign-flip term 3 | analytic (8); published (6); CMF (8); measure-contract (7) | analytic `2.0`; C&M `ξ₀` `13`; CMF `m̂₀` `0.5–0.65` | contract |
| M3 | omit `I²/R_*³` from `δM̂` | analytic Bc/Bs2/A2a/Ca (4); published F1/F2 (2); CMF matching arithmetic + G `δM̂` on all four stars (8); contract C13a; measure M4b/M2 (4) | `7.66e-3` (analytic, contract); CMF `δM̂` `3.1e-3 … 3.2e-2` | — (published K1/K2 inside `2e-2`, as in 4D) |
| M4 | omit the terminal surface atom | analytic Bc/Bs2/A2a/Cb/Cc/Da/Dc (7); contract C13b; measure M4a | `δM̂` `0.896` (analytic); `δM̂/R³` intercept `1.0` instead of 1 | published, CMF (shell `≲ 1e-6 δM̂`) — fires strongly on constant density / Newtonian as required |
| M5 | retired **profile-FD** derivative source `4πr²(ε+p)(Δε/Δp)_seg p̂₀*` (two coupled edit points) | CMF G/C on all four stars (8); published F2 (`2.6e-5`); measure M3b (`2.2e-6`, `4.8e-6`) | CMF `δM̂` `8.5e-3 / 2.8e-4 / 4.7e-3 / 1.0e-3`, `m̂₀` `≤ 1.0e-2` | analytic, contract (`dε/dp = 0`) |
| M6 | restore `p₀*(R_*) = 0` (shift the family so `p̂₀*` vanishes at `R_*`) | everything: analytic (13, `1.1e4`), published K1c/K2h (`1.0`), CMF (12), contract C5a (`1.7e12`)/C13b, measure M3b (3) | `ξ₀(R_*) = 0`, centre series `1.7e12` | — |
| M7 | literal-zero start | analytic Aa/Ab (`1.0`) and **A2a on the enlarged-centre fixture (`6.3e-1`)** (4); contract C5a/C5b | the initial value itself, and downstream displacement at `r₀ = 0.1 km` | published, CMF, measure (production `r₀ = 1e-5 km`: `≈ 6e-13`, below every bound — as ADR-0007 §8 predicted; A2 makes M7 load-bearing) |
| M8 | seed leak (`s ← stored ω̄`) | everything: analytic (15), published (6), CMF (8), contract Sa (`4.0e10` at seed `1e3`)/C5a/C5b, measure (7) | outputs `≈ 1.0` off; seed invariance `4e10` | — |
| M9 | drop `(1 + 8πr²p)` in term 4 | analytic B (`1.3e-3`) (7); published K1c (`1.2`)/F1/F2 (3); CMF G/C (`2.5e-3 … 8.3e-3`) (8); measure (6) | — | contract |
| M10 | superseded **nodal** `dε/dp` differential source | CMF G/C on all four stars (8); published F2 (`1.67e-4`); measure M3b (3) | **DS(CMF)-1 `δM̂` deficit `6.14 % / 5.10 % / 4.61 % / 1.97 %` at 1.0/1.4/1.6/2.0 M☉** (`δM̂ = 934.11 / 892.93 / 825.97 / 562.57` vs `995.26 / 940.92 / 865.87 / 573.87` km³) | analytic, contract (`dε/dp = 0`) |

ADR-0008 line I asks M10 to recreate a `≥ 3 %` deficit "under the accepted fixture": on the ADR-0008
fixture star (1.6 M☉) it is `4.61 %`, on 1.0 and 1.4 M☉ `6.14 %` and `5.10 %`; on the 2.0 M☉ target
star it is `1.97 %` (its selected `ε_c` moved under ADR-0009 V7b from `1.2983e15` to `1.2980e15 g cm⁻³`,
and the crust's share of `δM̂` is smallest there), still 200× the `1e-4` bound. Sweep totals: 10/10
build, 10/10 fire, 10/10 restored byte-identically.

M5 and M10 are distinct: M5 is the retired **profile-FD** derivative `(Δε/Δp)_seg` as the source,
M10 the **authoritative nodal** Steffen `dε/dp` differential form; both re-create the sub-node
deficit ADR-0008 corrected, and both are blind on constant density (`dε/dp = 0`) by construction.
M4 is blind on CMF and HT68 (shell `≲ 1e-6 δM̂`) and fires strongly on the constant-density and
Newtonian lines, as the brief requires. M7 fires on the enlarged-centre fixture A2 as well as on
the centre-series line.

## 14. ADR-0009 background checks

Every CMF star consumed (four target-mass stars; the H ladder at five resolutions; the J stars and
their sequence neighbours at three resolutions) and every HT68 star (eight) was checked: solve
status `SURFACE_REACHED`, last node `p == PressureCutoff()` **exactly**, finite `ε_*`, `m`, `dε/dp`
at the surface, strictly increasing radii with no sub-cutoff node, `n ≥ 4`. All pass. `R_*` is
partition-invariant to `3.4e-11` over 5000…80000 (ADR-0009 floor). No 2500 star and no truncated
star exists in any artifact or in this campaign.

## 15. Seven-artifact preservation and first-order protection

`shasum -a 256 tests/baselines/*.tsv` equals §1 at entry, after the detector sweep, and at exit; in
particular `hartle_I_dscmf1_debug.tsv = a21d4c3f6d89322cb4cca7a073e0e142f180e529332d704b7b16a145eae741c9`.
`hartle_normalization_*`, `hartle_first_order_physics_*`, `hartle_moment_inertia_*` pass in the
final suites (§17 field). No pre-ADR-0009 hash was used as authority.

## 16. Production diff

`git diff -- CompactStar/ CMakeLists.txt` = **0 lines** at every stage; `RotationSolver.cpp` SHA-256
identical before and after the detector sweep. Tracked changes: `tests/rotation/hartle_monopole_reference.hpp`
(the Stieltjes route added; historical routes byte-identical in behaviour), new
`tests/rotation/eos_table_knots.hpp`, `tests/rotation/hartle_monopole_physics_cmf.cpp` (rewritten
campaign), `tests/rotation/hartle_monopole_physics_analytic.cpp` (Bs, A2 added),
`tests/rotation/hartle_monopole_published.cpp` (K2z, F added), this record and the document sync of
§20. The stale in-source doc comment on `HartleMonopoleResponse` ("not yet independently validated
… no baseline may be created") remains true and was not edited (production diff NONE).

## 17. Scientific status

> **`HARTLE O(OMEGA^2) MONOPOLE RESPONSE CHARACTERIZED — INDEPENDENT VALIDATION INCOMPLETE`**

VERIFIED requires (task §27) every line including ADR-0008 D's accepted monotonicity; D's
monotonicity clause is NOT MET (§10). Everything else — independent solver admissibility,
analytic validation, ADR-0008 A/B/C/E/F/G/J, the sequence identity, the published comparisons,
M1–M10 coverage, production unchanged — holds. No substantive production disagreement was found
(the C status does not apply). `GOVERNANCE.md` §3.1 state: **CORRECTION IMPLEMENTED — INDEPENDENT
REVALIDATION EXECUTED — ONE ACCEPTED CONVERGENCE CRITERION NOT MET — NO BASELINE**.

## 18. Baseline chronology

**Not earned.** `tests/baselines/hartle_monopole_dscmf1_debug.tsv` was not created, no regression
test was registered, and no `test: revalidate …` / `test: establish …` commit exists. The durable
evidence is committed as `test: characterize corrected Hartle monopole response` (identified by Git history and the delivery report).

## 19. INV-08 / INV-09 disposition

**INV-08:** `GOVERNED (ADR-0007 + ADR-0008 ACCEPTED) — MEASURE-COMPLETE O(Ω²) MONOPOLE SOURCE CONFORMED
AND INDEPENDENTLY CHARACTERIZED ON THE ADR-0009 TOV SURFACE CONTRACT; ADR-0008 VALIDATION D
MONOTONICITY NOT MET (NODE-PLACEMENT FLOOR 3.7e-5) — OWNER ADJUDICATION REQUIRED; NO BASELINE`.
`l = 2` remains not implemented, not validated. No O(Ω²) number may be cited as a validated result
until the owner adjudicates §20.

**INV-09:** the fixed-`ε_c` structural response (the `A_i`-side source) is independently
characterized; the baryon-conserving reduction (`B_i`, `A_B/B_B`, `Z_i`), `dn_i/dp` and its measure
completeness (ADR-0008 Q11) remain unimplemented. **Not resolved.**

## 20. Phase-4E disposition and the owner question

The Phase-5 structural fields carried by `HartleMonopoleResponse` / `PhysicalHartleMonopole` are the
intended Phase-5 interface and are physically characterized to `≤ 5e-5` against a fully independent
route; 4E stays blocked only on the baseline, which stays blocked on ADR-0008 D. **Phase 5 does not
begin.**

**Smallest next action (one):** owner adjudication of ADR-0008 Validation D's monotonicity clause
against the measured node-placement floor — a post-validation clarification of ADR-0008 (the same
mechanism as ADR-0009's V7 clarification), choosing between (a) taking D's spread bound (`≤ 1e-4`,
met at `3.7e-5`) as the operative convergence criterion with the O(h) node-placement terms of
ADR-0008 Q3/Q4 and INV-13 explicitly recognised, after which the first monopole baseline may be
created from the already-committed characterization evidence, or (b) a governed increment adopting
the optional EOS-knot partition (ADR-0008 Q2b) in production, which removes the measure-location
term but, as the knot oracle shows, not the first-order sampled-background term.
