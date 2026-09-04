# TOV surface artifact migration — stopped candidate audit

**Status: TOV SURFACE ARTIFACT MIGRATION FAILED — OWNER REVIEW REQUIRED**

**Current status (TOV-SURF-MA, 2026-09-03): TOV SURFACE MIGRATION ENVELOPE ADJUDICATED — ARTIFACT MIGRATION MAY RESUME** — see §19. The historical stop status above and §§1–18 are preserved unchanged.

Date: 2026-09-03. Task: TOV-SURF-M. One agent; serial execution.

No candidate artifact was copied into `tests/baselines`. The seven historical files
remain the current regression authorities. ADR-0009 remains **ACCEPTED / SOURCE
CONFORMED / NUMERICALLY VALIDATED**; this stop concerns the migration impact limit
for a nondefault historical resolution, not a new production change.

The copy gate stopped at the nontruncated `A_fixed_ec / 5000` row: its same-central-density
mass change is **+1.6200244697598425e-9 relative**, above the task §7 limit `1e-9`.
Task §16 prohibits copying when an impact limit remains unresolved. The candidate
matches the previously validated seven-resolution dry run to serialization precision.
The owner must decide the scope of the impact bound before migration can proceed.

## 1. Authenticated starting state and authority

Worktree: `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-tov-surface`.
Branch: `governance/tov-surface-contract`.
Starting HEAD: `94165f359ccfadf34bf64ede62e7bef9c581a067`.
The starting tree was clean; local HEAD, upstream and live remote agreed.
Relative to master `df859b5a73c4cac0c115f240744d89ce9f830b8d`: **18 ahead / 0 behind**.
Other worktrees and ancestor relationships were inspected. These ancestors remain:

- `7a223d3390415e6b12737b5a9fbd9c87cfeadab6` — ADR-0009 acceptance.
- `96c1425f8503c2845e72f5755b02840798f6b958` — historical V7 stop.
- `03b7d585b8cc807d78c7d366b1f18340f6e40f5b` — V7 clarification.
- `94165f359ccfadf34bf64ede62e7bef9c581a067` — validated implementation.

Change classes: generated-artifact candidate, regression/test orchestration, and
documentation. No scientific-semantic or numerical-method production change; the
GOVERNANCE §3.1 pre-baseline exception was not invoked.

Read authorities: `AGENTS.md`; `GOVERNANCE.md`; ADR-0009 (Q11–Q14 and its accepted
post-decision V7 clarification); `SCIENTIFIC_INVARIANTS.md` INV-06/13/14;
`TOV_SURFACE_CONTRACT_DERIVATION.md` §§11–13;
`TOV_SURFACE_IMPLEMENTATION.md` §§6, 8, 9; `TOV_RADIAL_RES_2500_AUDIT.md`;
`GRID_CONVERGENCE.md`; `TOV_REFERENCE.md` §5; `CURRENT_ARCHITECTURE.md`;
`MODERNIZATION_ROADMAP.md`. The heuristic review also used
`HARTLE_MOMENT_INERTIA.md` and existing Phase-4A/4B validation authority.

ADR-0009's Decision, PressureCutoff, SolveToProfile's absolute `1e-4 Msun` tolerance,
first-acceptable bisection policy, coarse N=24, stable-branch criterion, and inherited
adaptive driver history are unchanged. `R_*` denotes the finite EOS-floor
`p=p_cut` event, never a zero-pressure physical surface.

## 2. Old authority authenticated before generation

All seven files were hashed and copied to the evidence archive before any producer
edit. The trajectory hash in the task text was incomplete and did not match. The
hash gate was paused and the complete value was authenticated from the repository
bytes and historical `TOV_SURFACE_IMPLEMENTATION.md` §6; no hash was guessed.

| Artifact | Class | OLD SHA256 — still current |
|---|---|---|
| baryon_number_dscmf1_reference.tsv | B | `8da5799d21da2017dd7dc49dfec8571ade6efba22846a652796118f248d4a646` |
| grid_convergence_cmf_1p6_debug.tsv | A | `61d84ddcb87645197c5406c880b648fdf3bb9b0ed8c58350800ca2f2d296ff40` |
| grid_convergence_cmf_1p6_trajectory.tsv | A | `ca32863dabaa28fad63d5c36b287a3b94e9b6b85f11980bf2be4e65499d9a0c6` |
| hartle_I_dscmf1_debug.tsv | B | `ddf018579364e9b3bf24f1e3a3e2577e70fb09b8b84eb718237d9da683aa9d15` |
| passive_cooling_cmf_1p6_debug.tsv | B | `831744b0a206541fd0e24adc67876cc1ee4d02d89a580942a9fb0c6749999453` |
| tov_dscmf1_reference.tsv | B | `ba9f6ee51e501e5e5a2133f72d3d16f351e5c721eb3f7a7c04e4d922fbc13e28` |
| tov_path_equivalence_dscmf1.tsv | D | `bbf61e5fddb4709500f22a1eb11b1e20554f7463376619e86e96ea0a2540d871` |

## 3. Producer authentication and retained evidence

Durable archive: `/Users/keeper/.codex/diagnostics/tov-surf-m-20260903`.

- `evidence/old/`: untouched historical artifact bytes.
- `run1/` and `run2/`: independent canonical producer outputs, **uninstalled**.
- `producer-binaries/`, `producer-sources/`, `producer-CMakeCache.txt`: generation provenance.
- `candidate-producers.patch`: exact four-file test/producer candidate, subsequently restored.
- `evidence/producer-binary-hashes.json`: fixed executable hashes across both runs.
- `evidence/external-input-hashes.json`: ten authenticated EOS/reference input hashes.
- `evidence/production-hashes.json`: authenticated production tree bytes.
- `evidence/generation.json`: all twelve exact argv vectors, exit codes, runtimes and log names.
- `authenticated-prior-grid-dry.log`: preserved TOV-SURF-I-R comparison evidence.
- `SHA256SUMS.json`: final archive file integrity manifest (excluding itself).

Candidate producer patch SHA256:
`f2674572f78b0d9371d45eaa41856ff5ef56d5c78b807f2061fb926086c66060`.

All artifact numbers came from the C++ canonical producer. Python only orchestrated
processes and compared/reported existing output; it did not synthesize artifact rows.
Both generation directories were created clean outside `tests/baselines`.

For each `RUN` in `/tmp/compactstar-tov-migration-run1` and
`/tmp/compactstar-tov-migration-run2`, these commands ran **sequentially**, from the
named worktree. `EOS=/Users/keeper/Documents/CompactStar/data/compose`:

```sh
build/tests/passive_cooling_regression "$EOS" --emit "$RUN/passive_cooling_cmf_1p6_debug.tsv"
build/tests/tov_reference_cmf "$EOS" --emit "$RUN/tov_dscmf1_reference.tsv"
build/tests/grid_convergence_cmf "$EOS" --emit-dir "$RUN"
build/tests/hartle_moment_inertia_cmf "$EOS" --emit "$RUN/hartle_I_dscmf1_debug.tsv"
build/tests/baryon_number_cmf "$EOS" --emit "$RUN/baryon_number_dscmf1_reference.tsv"
build/tests/tov_path_equivalence_cmf "$EOS" --emit "$RUN/tov_path_equivalence_dscmf1.tsv"
```

All twelve processes exited 0. Emit-mode success establishes their applicable
producer checks, **not** the post-copy legacy-regression pass required for completion.

## 4. Same-build repeatability and candidate hashes

**Seven of seven candidates were byte-identical between the two independent runs.**
There were no differing fields or formatting differences. No regression tolerance
was widened. The following are candidate hashes, not promoted current authorities.

| Artifact | UNINSTALLED candidate SHA256 | Producer | Repeatability |
|---|---|---|---|
| passive_cooling_cmf_1p6_debug.tsv | `afcad5d078fe458bdb441278ce56caa2d025becc6b489d00e9baad91a101c0de` | passive_cooling_regression | byte-identical |
| tov_dscmf1_reference.tsv | `3d9af9129a6a4ffde9e0f8c5507a160f968a861c0cf9f3b089cceecab86b701a` | tov_reference_cmf | byte-identical |
| grid_convergence_cmf_1p6_debug.tsv | `2c68e2f7e871192e00322bcc12b6d8b13eb7fdbda4a6d8d7050f26f0271ce5eb` | grid_convergence_cmf | byte-identical |
| grid_convergence_cmf_1p6_trajectory.tsv | `7c84557742ec0ec747756d118b2837fb54095a94c10194d4ccc2d0a778f5b04f` | grid_convergence_cmf | byte-identical |
| hartle_I_dscmf1_debug.tsv | `a21d4c3f6d89322cb4cca7a073e0e142f180e529332d704b7b16a145eae741c9` | hartle_moment_inertia_cmf | byte-identical |
| baryon_number_dscmf1_reference.tsv | `7b036942f2ae599ace3b2fc8b9a7d91f6d11b5899ab7e8a88d2bb4ee6686493b` | baryon_number_cmf | byte-identical |
| tov_path_equivalence_dscmf1.tsv | `5c0f4b3bdb70921f8f2a869af10edc4d8f5ae3963a9d150e11ef859d21e1c678` | tov_path_equivalence_cmf | byte-identical |

## 5. Full delta methodology and row classification

`evidence/row-column-deltas.csv` and `.json` contain **1,045 field records**, including
all unchanged, changed and added fields: artifact, class, row key, column, old/new
serialized values, signed absolute and relative delta, semantics and explanation.
Decimal arithmetic at precision 40 uses `new-old` and `(new-old)/abs(old)`; a zero
old denominator or nonnumeric value has no relative delta. Summaries below round
only for readability. `evidence/metadata-comments.diff` separately preserves all
comment/provenance differences. No data columns were changed and no data rows removed.
The generic physical explanations in the field audit are not approval of a bound:
`evidence/delta-gates.json` records the failed 5000-row copy gate explicitly.

Row semantics were assigned before impact acceptance:

- Grid A is fixed central density. The old 2500 row is class-A truncation recovery;
  its large delta is not a surface-only tail. Other common A rows are surface-tail
  comparisons at unchanged central density; the 5000 row triggers the stop.
- Grid B is the target-mass workflow. At 2500 its selected central density changes
  because the missing-crust compensation disappears. At the other common rows it
  happens to select the same central density as before and as corrected A.
- The grid trajectory belongs to B (fixed target mass), with nine epochs per grid.
- TOV reference, passive cooling, Hartle I and baryon reference use target-mass
  production calls. The 1.0/1.4/1.6 target calls happen to retain their old central
  densities. The 2.0 target changes central density under V7b; it is not subject to
  the fixed-central-density old/new `1e-9` comparison.
- Path-equivalence rows hold central density fixed; the artifact records comparisons
  between the two ordinary paths, not independent target-mass roots.
- 1250 and 80000 rows are **NEW ROW — no historical comparator**.

| Artifact | Old → candidate rows | Changed existing fields | New rows | Removed rows |
|---|---|---|---|---|
| baryon_number_dscmf1_reference.tsv | 4 → 4 | 13 | 0 | 0 |
| grid_convergence_cmf_1p6_debug.tsv | 10 → 14 | 111 | 4 | 0 |
| grid_convergence_cmf_1p6_trajectory.tsv | 45 → 63 | 175 | 18 | 0 |
| hartle_I_dscmf1_debug.tsv | 4 → 4 | 57 | 0 | 0 |
| passive_cooling_cmf_1p6_debug.tsv | 9 → 9 | 94 | 0 | 0 |
| tov_dscmf1_reference.tsv | 3 → 3 | 3 | 0 | 0 |
| tov_path_equivalence_dscmf1.tsv | 17 → 17 | 34 | 0 | 0 |

## 6. Class A — grid debug candidate and failed copy gate

The candidate uses the accepted seven-resolution set:
**1250, 2500, 5000, 10000, 20000, 40000, 80000** (implementation §9.7).
Both arms contain seven complete `SURFACE_REACHED` stars; all use the governed
surface and pass the unchanged target-mass residual. The arms select the same
central density at every corrected resolution. Profile sizes are respectively
332, 662, 1320, 2636, 5269, 10536, 21068. The corrected 2500 point is on the same
physical branch; it is not a truncated-star detector.

At 2500, A restores mass and crust at fixed central density; B removes the old
central-density compensation. In B, old `epsilon_c=7.328604318265e14` becomes
`7.312533426775e14 g/cm^3`. Common-resolution M/R/B changes follow; all other
structural and thermal columns are in the field audit.

| Row key | Column | Old | Candidate | Signed Δ | Signed Δ / abs(old) |
|---|---|---|---|---|---|
| A_fixed_ec / 2500 | achieved_M | 1.597681155885e+00 | 1.599975835285e+00 | +0.0022946794 | +0.001436256159 |
| A_fixed_ec / 2500 | R_km | 1.290422815654e+01 | 1.347335819259e+01 | +0.5691300361 | +0.04410415169 |
| A_fixed_ec / 2500 | B | 2.121057787217e+57 | 2.124747791571e+57 | +3.690004354e+54 | +0.001739700057 |
| A_fixed_ec / 5000 | achieved_M | 1.599975832701e+00 | 1.599975835293e+00 | +2.592e-09 | +1.62002447e-09 |
| A_fixed_ec / 5000 | R_km | 1.346345807665e+01 | 1.347335819373e+01 | +0.00990011708 | +0.0007353324104 |
| A_fixed_ec / 5000 | B | 2.124579914504e+57 | 2.124579918436e+57 | +3.932e+48 | +1.850718805e-09 |
| A_fixed_ec / 10000 | achieved_M | 1.599975834172e+00 | 1.599975835301e+00 | +1.129e-09 | +7.056356577e-10 |
| A_fixed_ec / 10000 | R_km | 1.346832307596e+01 | 1.347335819395e+01 | +0.00503511799 | +0.0003738489166 |
| A_fixed_ec / 10000 | B | 2.124575695480e+57 | 2.124575697169e+57 | +1.689e+48 | +7.949822657e-10 |
| A_fixed_ec / 20000 | achieved_M | 1.599975834821e+00 | 1.599975835302e+00 | +4.81e-10 | +3.006295405e-10 |
| A_fixed_ec / 20000 | R_km | 1.347101807557e+01 | 1.347335819349e+01 | +0.00234011792 | +0.0001737150011 |
| A_fixed_ec / 20000 | B | 2.124509533317e+57 | 2.124509534034e+57 | +7.17e+47 | +3.3748966e-10 |
| A_fixed_ec / 40000 | achieved_M | 1.599975835170e+00 | 1.599975835302e+00 | +1.32e-10 | +8.250124602e-11 |
| A_fixed_ec / 40000 | R_km | 1.347268057533e+01 | 1.347335819385e+01 | +0.00067761852 | +5.029574599e-05 |
| A_fixed_ec / 40000 | B | 2.124531386072e+57 | 2.124531386268e+57 | +1.96e+47 | +9.225563872e-11 |
| B_fixed_mass / 2500 | achieved_M | 1.600076364417e+00 | 1.599975835285e+00 | -0.000100529132 | -6.282770887e-05 |
| B_fixed_mass / 2500 | R_km | 1.290422815654e+01 | 1.347335819259e+01 | +0.5691300361 | +0.04410415169 |
| B_fixed_mass / 2500 | B | 2.124619525864e+57 | 2.124747791571e+57 | +1.28265707e+53 | +6.037114196e-05 |
| B_fixed_mass / 5000 | achieved_M | 1.599975832701e+00 | 1.599975835293e+00 | +2.592e-09 | +1.62002447e-09 |
| B_fixed_mass / 5000 | R_km | 1.346345807665e+01 | 1.347335819373e+01 | +0.00990011708 | +0.0007353324104 |
| B_fixed_mass / 5000 | B | 2.124579914504e+57 | 2.124579918436e+57 | +3.932e+48 | +1.850718805e-09 |
| B_fixed_mass / 10000 | achieved_M | 1.599975834172e+00 | 1.599975835301e+00 | +1.129e-09 | +7.056356577e-10 |
| B_fixed_mass / 10000 | R_km | 1.346832307596e+01 | 1.347335819395e+01 | +0.00503511799 | +0.0003738489166 |
| B_fixed_mass / 10000 | B | 2.124575695480e+57 | 2.124575697169e+57 | +1.689e+48 | +7.949822657e-10 |
| B_fixed_mass / 20000 | achieved_M | 1.599975834821e+00 | 1.599975835302e+00 | +4.81e-10 | +3.006295405e-10 |
| B_fixed_mass / 20000 | R_km | 1.347101807557e+01 | 1.347335819349e+01 | +0.00234011792 | +0.0001737150011 |
| B_fixed_mass / 20000 | B | 2.124509533317e+57 | 2.124509534034e+57 | +7.17e+47 | +3.3748966e-10 |
| B_fixed_mass / 40000 | achieved_M | 1.599975835170e+00 | 1.599975835302e+00 | +1.32e-10 | +8.250124602e-11 |
| B_fixed_mass / 40000 | R_km | 1.347268057533e+01 | 1.347335819385e+01 | +0.00067761852 | +5.029574599e-05 |
| B_fixed_mass / 40000 | B | 2.124531386072e+57 | 2.124531386268e+57 | +1.96e+47 | +9.225563872e-11 |

For the blocking A/5000 row, both old and candidate serialized central densities
are `7.312533426775e14 g/cm^3` (inherited input `731253342677476.125`). Its
`ΔM=+2.592e-9 Msun` gives `ΔM/M=+1.6200244697598425e-9` using the artifact decimals.
The radius moves **+9.90011708 m**, and `ΔB/B=+1.8507188048e-9`.
This is outside the task's universal-looking fixed-density `1e-9` mass wording;
it is not the known 2500 truncation and not a target-root shift.

Across the complete corrected matrix, comparison with the authenticated TOV-SURF-I-R
dry run yields maximum serialized relative difference **3.6889440369160784e-13**,
consistent with printed precision. Default 10000 has the expected +5.03511799 m
radius shift and `ΔM/M=7.05636e-10`. The longer old 5000 surface interval explains
the larger tail physically, but does not authorize relaxing or silently narrowing
an explicit owner limit. The candidate remains uninstalled.

## 7. Class A — full trajectory candidate

45 historical rows become 63: nine epochs for each of seven resolutions.
Eighteen rows at 1250/80000 are new; no old row is removed. The 2500 trajectory
changes because its background is now complete, while other resolutions acquire
the corrected surface factors. The terminal epoch (`1e6 yr`) comparisons are:

| Row key | Column | Old | Candidate | Signed Δ | Signed Δ / abs(old) |
|---|---|---|---|---|---|
| 2500 / 1.000000000000e+06 | Tinf_K | 4.575298532239e+05 | 4.663644064456e+05 | +8834.553222 | +0.01930923886 |
| 5000 / 1.000000000000e+06 | Tinf_K | 4.668964825035e+05 | 4.667991415157e+05 | -97.3409878 | -0.0002084851599 |
| 10000 / 1.000000000000e+06 | Tinf_K | 4.669660595508e+05 | 4.669165849205e+05 | -49.4746303 | -0.0001059490926 |
| 20000 / 1.000000000000e+06 | Tinf_K | 4.668539919058e+05 | 4.668310110647e+05 | -22.9808411 | -4.922490007e-05 |
| 40000 / 1.000000000000e+06 | Tinf_K | 4.668023052084e+05 | 4.667956522849e+05 | -6.6529235 | -1.425212221e-05 |
| 1250 / 1.000000000000e+06 | Tinf_K | — | 4.677466154262e+05 | — | — |
| 80000 / 1.000000000000e+06 | Tinf_K | — | 4.667775322991e+05 | — | — |

The full CSV gives Tinf, heat capacity, total/component neutrino luminosity,
photon luminosity and temperature derivative at every common epoch, with signed
deltas. The 175 changed common fields are not summarized solely by a file maximum.

## 8. Remeasured finite-resolution behavior

These diagnostics were recomputed from the actual canonical candidate bytes,
not copied from the previous dry run. The finest-three effective order uses the
producer's `dr_eff_km` spacing ratios and differences when signs and numerical
floor permit it. The trajectory norm is `max_epoch |ln(Tinf_a/Tinf_b)|` and its
order uses contraction of the finest adjacent-pair norms. This does not establish
a universal or asymptotic convergence order.

| Quantity | Finest-three diagnostic order / finest-pair contraction |
|---|---|
| achieved_M | not reliably measurable (solver floor / nonmonotone) |
| R_km | not reliably measurable (solver floor / nonmonotone) |
| B | 2.234619361 |
| Cstar_1e8 | 2.246126826 |
| Lnu_1e8 | 0.164427240 |
| Lgamma_1e8 | not reliably measurable (solver floor / nonmonotone) |
| dlnT_dt_1e8 | 0.147859455 |
| terminal_Tinf | 0.964651796 |
| trajectory_norm | 0.146235960 |

Adjacent trajectory norms, 1250→2500 through 40000→80000:
`2.959412654e-3, 9.317447140e-4, 4.535652739e-4, 2.381947703e-4,
1.434199071e-4, 1.295971832e-4`.
M and R are already at the solver/output floor; Lgamma is not monotonic on the
finest triple. B/Cstar and thermal trajectories retain sampling error.
Candidate G1–G5 pass. The thermal time-integration floor is `1.263e-8` versus
`5.111e-4` default-grid radial trajectory error; its original subdominance gate
is unchanged. The historical truncated-grid radius-error interpretation is
superseded for these candidates, but current baseline authority was not migrated.

## 9. Class B — TOV reference

All three target rows retain their achieved masses at the artifact's six-decimal
precision. Every external CompOSE reference value is unchanged. CompactStar moves
closer to that curve while the finite-floor-to-vacuum radius difference remains.

| Row key | Column | Old | Candidate | Signed Δ | Signed Δ / abs(old) |
|---|---|---|---|---|---|
| 1.000000 | achieved_M | 1.000044 | 1.000044 | +0 | +0 |
| 1.000000 | R_production_km | 13.426323 | 13.427522 | +0.001199 | +8.930218646e-05 |
| 1.000000 | R_official_km | 13.472825 | 13.472825 | +0 | +0 |
| 1.400000 | achieved_M | 1.400022 | 1.400022 | +0 | +0 |
| 1.400000 | R_production_km | 13.545323 | 13.546350 | +0.001027 | +7.581952826e-05 |
| 1.400000 | R_official_km | 13.574875 | 13.574875 | +0 | +0 |
| 1.600000 | achieved_M | 1.599976 | 1.599976 | +0 | +0 |
| 1.600000 | R_production_km | 13.468323 | 13.473358 | +0.005035 | +0.0003738401581 |
| 1.600000 | R_official_km | 13.494859 | 13.494859 | +0 | +0 |

The candidate B7 replacement retains the omitted-layer upper bound (each fraction
in `(0,1)`) and every independent reference comparison. The old mass-independent
fraction-spread assertion failed at `0.209690`; it depended partly on the old
one-step-short radius. It is replaced by direct complete-star, `p=p_cut`, same-εc
M/R consistency using ADR-0009's existing `1e-9`/`1e-8` bounds. The approximate
fraction spread remains reported, without inventing a wider threshold. Official
radius and maximum-mass budgets remain `0.005` and `0.01`. This test edit passed
in both candidate runs, but was archived and restored after the migration stop.

## 10. Class B — passive cooling

No thermal configuration, equation, stepper or regression tolerance changed.
The candidate's source provenance names validated science `94165f3` and retains
the original thermal-configuration provenance explicitly. Profile version stays 5;
node count moves 2635→2636. R moves `13.46832307596→13.47335819395 km`.
The same canonical default star has surface lapse
`0.8057091628456→0.8057905213345` (+`1.009774900e-4` relative),
squared lapse +`2.019651764e-4`, and redshifted photon area +`9.499538097e-4`.
The latter two are diagnostic evaluations of existing formulas, not synthesized
baseline rows; see `evidence/surface-factor-deltas.json`.

All nine Tinf epoch changes are recorded below; the full field audit includes
every luminosity, heat capacity, surface/base temperature and derivative field.

| Row key | Column | Old | Candidate | Signed Δ | Signed Δ / abs(old) |
|---|---|---|---|---|---|
| 1.000000000000e+02 | Tinf_K | 1.000000000000e+09 | 1.000000000000e+09 | +0 | +0 |
| 3.019951720402e+02 | Tinf_K | 1.208661442309e+07 | 1.208661387872e+07 | -0.54437 | -4.50390805e-08 |
| 1.000000000000e+03 | Tinf_K | 8.307535178790e+06 | 8.307533560605e+06 | -1.618185 | -1.947852119e-07 |
| 3.019951720402e+03 | Tinf_K | 6.168097604078e+06 | 6.168094181117e+06 | -3.422961 | -5.549459849e-07 |
| 1.000000000000e+04 | Tinf_K | 4.496566670363e+06 | 4.496559209229e+06 | -7.461134 | -1.659295758e-06 |
| 3.019951720402e+04 | Tinf_K | 3.310461884891e+06 | 3.310446731326e+06 | -15.153565 | -4.577477563e-06 |
| 1.000000000000e+05 | Tinf_K | 2.242118838269e+06 | 2.242087151311e+06 | -31.686958 | -1.413259523e-05 |
| 3.019951720402e+05 | Tinf_K | 1.329027695287e+06 | 1.328973673180e+06 | -54.022107 | -4.064784142e-05 |
| 1.000000000000e+06 | Tinf_K | 4.669660595508e+05 | 4.669165849205e+05 | -49.4746303 | -0.0001059490926 |

Terminal Tinf moves **−1.059490926e-4 relative**. Maximum absolute relative movements
are `6.354870176e-4` total neutrino luminosity, `8.472193628e-4` modified-Urca
luminosity, and `1.554014926e-4` photon luminosity. These agree with the validated
surface-impact preview; actual luminosity also depends on the evolving temperature.
The old/new surface shift was not used to set regression tolerance.

## 11. Class B — first-order Hartle I

All four rows are target-mass calls. The first three retain the old selected central
density; their I changes are surface-tail contributions. The 2.0 target changes
central density under V7b. No Hartle equation or first-order physics bound changed.

| Row key | Column | Old | Candidate | Signed Δ | Signed Δ / abs(old) |
|---|---|---|---|---|---|
| 1.000000000000e+00 | R_km | 1.342632308196e+01 | 1.342752215526e+01 | +0.0011990733 | +8.930764534e-05 |
| 1.000000000000e+00 | epsilon_c | 4.545504050785e+14 | 4.545504050785e+14 | +0 | +0 |
| 1.000000000000e+00 | I_production_km3 | 8.699575833405e+01 | 8.699575838205e+01 | +4.8e-08 | +5.517510384e-10 |
| 1.400000000000e+00 | R_km | 1.354532306496e+01 | 1.354635043274e+01 | +0.00102736778 | +7.58466797e-05 |
| 1.400000000000e+00 | epsilon_c | 6.164882705061e+14 | 6.164882705061e+14 | +0 | +0 |
| 1.400000000000e+00 | I_production_km3 | 1.356161305669e+02 | 1.356161306127e+02 | +4.58e-08 | +3.377179382e-10 |
| 1.600000000000e+00 | R_km | 1.346832307596e+01 | 1.347335819395e+01 | +0.00503511799 | +0.0003738489166 |
| 1.600000000000e+00 | epsilon_c | 7.312533426775e+14 | 7.312533426775e+14 | +0 | +0 |
| 1.600000000000e+00 | I_production_km3 | 1.595871413252e+02 | 1.595871415881e+02 | +2.629e-07 | +1.647375834e-09 |
| 2.000000000000e+00 | R_km | 1.271232318396e+01 | 1.271647774735e+01 | +0.00415456339 | +0.0003268138585 |
| 2.000000000000e+00 | epsilon_c | 1.298349261930e+15 | 1.297985389789e+15 | -3.63872141e+11 | -0.0002802575175 |
| 2.000000000000e+00 | I_production_km3 | 1.937231442198e+02 | 1.937221161501e+02 | -0.0010280697 | -5.306901786e-06 |

At 2.0, I moves `193.7231442198→193.7221161501 km^3` (`−5.306901786e-6` relative).
The selected central density moves `1.298349261930e15→1.297985389789e15 g/cm^3`.
Independent-reference I, volume I, residuals and fit diagnostics were recomputed
on the same corrected background; every changed field has its own audit record.
A large fractional movement of an already tiny residual is not a large change in I.

The candidate B3a replaces sampled-I spread versus sampled-R spread with I's own
contraction toward the finest computed value. Relative errors against resolution
40000 are `2.044e-4, 4.632e-5, 4.024e-5, 1.941e-5, 0` at 2500/5000/10000/20000/40000.
The finest four contract; no new continuum-accuracy tolerance was invented. All
five resolutions remain in independent B3b, whose unchanged `1e-3` gate has worst
error `4.027e-5`. Phase-4A/4B normalization and physics tolerances remain unchanged.
This test edit passed in both candidate runs and was then archived/restored.

## 12. Class B — baryon number

These are target-mass rows. The first three hold their historical central density
and satisfy the default-grid fixed-density tail expectation. The 2.0 row is V7b:
its larger B movement accompanies the permitted target-root change; it is not
compared against the fixed-density `1e-9` tail bound.

| Row key | Column | Old | Candidate | Signed Δ | Signed Δ / abs(old) |
|---|---|---|---|---|---|
| 1 | achieved_M_Msun | 1.0000438094972419 | 1.0000438097285522 | +2.313103e-10 | +2.313001669e-10 |
| 1 | R_km | 13.426323081955102 | 13.427522155261627 | +0.001199073307 | +8.930764582e-05 |
| 1 | ec_km^-2 | 454550405078491.75 | 454550405078491.75 | +0 | +0 |
| 1 | B | 1.2738887310935455e+57 | 1.2738887314075874e+57 | +3.140419e+47 | +2.465222373e-10 |
| 1.3999999999999999 | achieved_M_Msun | 1.4000217830781583 | 1.4000217832812236 | +2.030653e-10 | +1.450443861e-10 |
| 1.3999999999999999 | R_km | 13.545323064955092 | 13.546350432736507 | +0.001027367781 | +7.584667981e-05 |
| 1.3999999999999999 | ec_km^-2 | 616488270506054.5 | 616488270506054.5 | +0 | +0 |
| 1.3999999999999999 | B | 1.8321833625787515e+57 | 1.8321833628708972e+57 | +2.921457e+47 | +1.594522175e-10 |
| 1.6000000000000001 | achieved_M_Msun | 1.5999758341716293 | 1.5999758353006395 | +1.1290102e-09 | +7.056420328e-10 |
| 1.6000000000000001 | R_km | 13.468323075955098 | 13.473358193949251 | +0.005035117994 | +0.0003738489169 |
| 1.6000000000000001 | ec_km^-2 | 731253342677476.12 | 731253342677476.12 | +0 | +0 |
| 1.6000000000000001 | B | 2.1245756954797212e+57 | 2.124575697168995e+57 | +1.6892738e+48 | +7.951111385e-10 |
| 2 | achieved_M_Msun | 2.0000952962048859 | 1.9999985861663816 | -9.67100385e-05 | -4.835271534e-05 |
| 2 | R_km | 12.712323183955165 | 12.716477747345323 | +0.00415456339 | +0.0003268138585 |
| 2 | ec_km^-2 | 1298349261929558.8 | 1297985389788764.5 | -3.638721408e+11 | -0.0002802575173 |
| 2 | B | 2.7457630611447906e+57 | 2.7455916278968613e+57 | -1.714332479e+53 | -6.243555766e-05 |

The 2.0 selected density is `1298349261929558.8→1297985389788764.5 g/cm^3`.
Both its old mass `2.0000952962048859` and candidate mass `1.9999985861663816 Msun`
satisfy the unchanged absolute `1e-4 Msun` target tolerance. Implementation §8
already established that the old search consumed a truncated bisection member;
this task did not alter or re-adjudicate that search policy.

The legacy baryon artifact header `ec_km^-2` is mislabeled: the canonical producer
writes the CGS central-density value. The header and bytes were not hand-corrected;
the audit labels the physical interpretation explicitly. No baryon-current formula,
INV-14 or `1e-15` structural/B regression tolerance changed.

## 13. Class D — path equivalence

All 17 row keys remain. Only the two matching node-count fields change (34 fields).
Both paths append the same governed event. All M/R/B/I/profile comparison fields
remain zero, including ULP differences; every row is equivalent. Invalid/below-floor
and incomplete profiles fail closed in the unchanged implementation and its tests.
No second surface locator exists. The candidate passes both producer invocations.

The canonical producer's comment header omits a historical Phase-3E-I2 annotation
that had been present in the old file. `metadata-comments.diff` preserves this
text-only difference; it does not reopen the closed geometry gap or change B.
Historical geometry records remain untouched.

## 14. Tolerances and heuristic candidate disposition

Run1/run2 byte equality is the repeatability evidence. Existing tight comparisons
were retained: passive state `1e-5` and luminosity `1e-4`; baryon structural/B
`1e-15`; path equality/ULP requirements and all independent scientific bounds.
No tolerance was derived from the migration amplitude.

Only these four test/producer files were edited temporarily, then restored exactly
from the starting HEAD after archiving their patch, source and executable bytes:

- `tests/thermal/grid_convergence_cmf.cpp`: seven-resolution producer set, complete
  star/surface checks, G1 complete-branch and G3 event-partition checks; G2/G4/G5
  retained; failed runs prevented from emitting; provenance headers updated.
- `tests/core/tov_reference_cmf.cpp`: bounded-layer B7 retained, obsolete fraction
  heuristic replaced by direct governed-event consistency.
- `tests/rotation/hartle_moment_inertia_cmf.cpp`: obsolete radius-spread B3a replaced
  by I's own finite-grid contraction; independent physics gates unchanged.
- `tests/thermal/passive_cooling_regression.cpp`: provenance metadata only.

No production science file was edited. The restored build completed successfully.
Candidate tests are evidence in the archive, not installed changes or promoted
regression authority.

## 15. Hash occurrence audit and historical preservation

The audit found 60 occurrences of full old hashes in 14 documentation files;
`evidence/hash-occurrence-audit.json` records each path/line/context and classification.
They are historical validation snapshots, including Phase-3 closeout, Phase-4
validation records, TOV path equivalence and historical implementation sections.
No historical hash was globally replaced. No current authority pointer or hash
registry was promoted because the copy gate failed. The historical radial-2500
audit and surface derivation are unchanged.

## 16. Actual test state and exclusions

Starting source was built and tested serially, before the producer candidate edits:

```sh
cmake --build build -j6
ctest --test-dir build --output-on-failure -E 'hartle_monopole_physics|hartle_monopole_published|passive_cooling_regression|tov_reference_cmf|hartle_moment_inertia_cmf|baryon_number_cmf'
cmake --build build-selfcontained -j6
ctest --test-dir build-selfcontained --output-on-failure -E 'hartle_monopole_physics|hartle_monopole_published'
ctest --test-dir build --output-on-failure -R '^(passive_cooling_regression|tov_reference_cmf|hartle_moment_inertia_cmf|baryon_number_cmf)$'
```

| Execution | Included / excluded | Result | Wall time |
|---|---|---|---|
| Authenticated main behavior/scientific suite | 37 of 44; four expected legacy failures and three independent Phase-4D tests excluded from this invocation | 37/37 PASS | 463.76 s |
| Four known legacy regressions, separate serial invocation | 4 | 0/4 PASS; four expected failures | 24.05 s |
| Authenticated self-contained suite | 20 of 22; two independent Phase-4D tests excluded | 20/20 PASS | 65.40 s |
| Candidate producers | six per run, twelve processes total | 12/12 exit 0; 7/7 artifacts byte-identical | exact individual times in generation.json |
| Post-copy full non-Phase-4D suite | would include 41 of 44 | NOT RUN: copy gate failed | — |
| Post-copy self-contained suite | would include 20 of 22 | NOT RUN: copy gate failed | — |
| Post-copy producer reproduction | seven artifacts | NOT RUN: no copy occurred | — |

The combined starting main executions covered all 41 applicable non-Phase-4D tests:
37 passed, four known legacy failures, total serial wall time **487.81 s**. This is
not a successful final migrated suite. The 13 TOV surface tests, path equivalence,
first-order normalization/physics, thermal core, EOS and applicable baryon/geometry
checks in the 37-test suite passed. The shared heat-capacity fixture was never
used by simultaneous suite runs.

The independent exclusions are exactly `hartle_monopole_physics_analytic`,
`hartle_monopole_physics_cmf` and `hartle_monopole_published` in the full build;
the self-contained build has the analytic and published tests only. Existing
monopole API/measure contract checks included in the behavior suite are not the
reserved corrected Phase-4D independent scientific revalidation.

Expected legacy failures remain current after restoration:

- `passive_cooling_regression`: old trajectory/reference mismatch.
- `tov_reference_cmf`: obsolete B7 fraction spread `0.209690` fails `<0.15`.
- `hartle_moment_inertia_cmf`: obsolete B3a compares I spread `6.573e-5` with
  event-radius spread `9.477e-11`.
- `baryon_number_cmf`: old M/R/εc/B reference mismatches (eight field-group failures).

No new unexpected executable/test failure occurred. The unexpected **migration
limit conflict** is the A/5000 old/new mass bound above. All raw failure output is
retained in `evidence/ctest-legacy-before.log`; it has not been hidden or normalized.

## 17. PROPOSED validation-scope addendum — owner decision pending

*Adjudicated in §19 (TOV-SURF-MA, 2026-09-03): alternative 1, with the derived envelope of §19.5 Q3; no bound widened. This section is retained as the historical question.*

This proposed ADR-0009 validation addendum records alternatives under AGENTS §8 and
GOVERNANCE's ambiguity rule. **Status: PROPOSED; no decision accepted here.** It does
not edit ADR-0009's accepted scientific Decision.

Context: the canonical default-grid V7a impact experiment passes, and the corrected
seven-resolution matrix matches the accepted implementation dry run. A common,
nontruncated historical 5000-resolution row has a longer omitted surface interval,
causing `ΔM/M=1.62002447e-9`. TOV-SURF-M §7 states `<=1e-9` without explicitly
limiting that historical migration bound to the default resolution.

Alternatives for the owner:

1. Confirm that the existing `1e-9` V7a old/new bound governs the inherited canonical
   default-resolution experiment, while adjudicating historical nondefault rows
   under an explicitly recorded resolution-dependent migration envelope. The
   corrected-vs-corrected partition bounds remain unchanged.
2. Require `1e-9` for every nontruncated historical resolution. This candidate then
   fails migration; any further investigation or scientific change requires a
   separate governed task, because this task authorizes no production changes.

No alternative is selected and no bound is widened by this record. The archived
exact patch/candidates can support a resumed migration after the owner resolves
this scope question; they are not authorization to copy.

## 18. Final disposition and dependency gate

Seven baseline overwrites: **0**. New current artifact hashes: **none**.
Production diff: **none**. Four temporary producer edits: archived and restored.
Repository change remaining: this stopped-attempt record only, uncommitted.
No migration commit was made and nothing was pushed: the requested atomic
`test: migrate TOV surface artifacts` commit requires successful candidate gates.
HEAD and live upstream remain `94165f359ccfadf34bf64ede62e7bef9c581a067`; commit
histories are equal. The worktree contains this new untracked record, so it is
not described as clean. No merge, squash or history rewrite occurred.

The final integrity check re-authenticates all seven historical baseline hashes,
the full production hash map, accepted ADR bytes, preserved ancestors, and local/
remote equality; `git diff --check` and the restored build are clean.
See `final-integrity.json` in the archive.

Corrected Phase-4D independent revalidation performed: **NO**.
First monopole scientific baseline created: **NO**. Phase 5 begun: **NO**.
Current architecture/invariant/roadmap and ADR index remain unchanged: migration
is still required. INV-06/INV-13 validation is not revoked by this migration-scope stop.

Dependency gate: resolve this migration impact scope; complete the authorized
artifact audit/copy/reproduction and full non-Phase-4D green suites; then corrected
Phase-4D independent monopole revalidation against ADR-0007 + ADR-0008 may run on
the migrated ADR-0009 backgrounds. Only after that successful independent
revalidation may the first monopole scientific baseline be created.

**Exactly one recommended next task:** obtain owner adjudication of the
nondefault-resolution fixed-εc migration envelope, specifically the historical
5000 row, before resuming TOV-SURF-M from the preserved candidates.

## 19. OWNER-ADJUDICATED migration-scope clarification — TOV-SURF-MA (2026-09-03)

**Status: TOV SURFACE MIGRATION ENVELOPE ADJUDICATED — ARTIFACT MIGRATION MAY RESUME.**

This section resolves the §17 question. It is a validation-scope clarification of the
artifact-migration gate only. It changes no production source, no TOV physics, no
ADR-0009 accepted Decision (Q1–Q14), no corrected-vs-corrected scientific bound, and no
invariant; it copies no artifact. Change class: documentation. Governing authority:
ADR-0009 Q11/Q12 and its V7 post-acceptance clarification; INV-06/INV-13/INV-14;
`GOVERNANCE.md` §2–§3 (the §3.1 exception is not invoked).

### 19.1 Authentication

HEAD `94165f359ccfadf34bf64ede62e7bef9c581a067` on `governance/tov-surface-contract`;
the only untracked file was this record. The seven current baseline hashes equal §2.
The corrected library `build/libCompactStar.a` is newer than every production source and
`CompactStar/Core/src/TOVSolver.cpp`, `CompactStar/Core/TOVSolver.hpp`,
`CompactStar/Core/src/NStar.cpp` hash-match `evidence/production-hashes.json` of the
TOV-SURF-M archive. Recomputed independently from the archive bytes: run1 and run2 are
byte-identical for all seven candidates (§4 hashes reproduced), and the grid-debug
candidate agrees with the TOV-SURF-I-R validated seven-resolution dry run
(`authenticated-prior-grid-dry.log`, `SURFACE_GRID` rows) to a maximum relative
difference of `3.69e-13` over 14 rows × 9 columns — 13-significant-digit serialization
of the same numbers. **The candidates are the validated corrected solver's bytes**, not
new results.

### 19.2 Where the `1e-9` old→new mass bound came from, and what it covered

| Source | Text (verbatim scope) |
|---|---|
| `docs/validation/TOV_SURFACE_CONTRACT_DERIVATION.md:192` | impact map "Measured at `radial_res = 10000`, production last node → event"; four canonical εc; `ΔM/M` = `2.3e-10`, `1.5e-10`, `7.1e-10`, `3.8e-10`; summary "`M` `≤ 7e-10`" (`:209`) |
| `docs/validation/TOV_SURFACE_CONTRACT_DERIVATION.md:248` | V7: "1.0 / 1.4 / 1.6 / 2.0 M☉ canonical stars — `M` within `1e-9`, `R_*` within the impact map (+1.0/+1.2/+5.0/+3.7 m ± 0.5 m)" — `1e-9` is the rounded-up envelope of those four 10000-resolution measurements (margin 1.4×) |
| `docs/validation/TOV_SURFACE_CONTRACT_DERIVATION.md:245` | V4: corrected-vs-corrected partition invariance "`M` `≤ 1e-9`, `R_*` `≤ 1e-8`" — a different quantity (two decades of margin) |
| `docs/adr/ADR-0009-tov-surface-event-and-termination.md:172` | Alternative C "Migration risk — every durable artifact carrying `R_*` moves by the impact map (+1–5 m **at the canonical masses**; … `M` `≤ 7e-10`)"; `:206` "any delta outside the impact map is a STOP" |
| `docs/adr/ADR-0009-tov-surface-event-and-termination.md:239` | V7a: "compare old and corrected solves at identical inherited εc. The existing `|ΔM/M| <= 1e-9` and radius impact map apply here. The four inherited εc values are …"; `:249` "The `1e-9` bound belongs only to V7a" |
| `docs/validation/TOV_SURFACE_IMPLEMENTATION.md:98` | the old comparator rows came from "the production `SolveToProfile` at the default resolution 10000" |
| `tests/core/tov_surface_contract.cpp:104`, `:109`, `:110` | V7a solves `Solve(s, densities[i], 10000)` and requires `|dm| <= 1e-9` and `|dr − expectedR| <= 0.5 m` |

The bound was therefore derived from, measured at, and tested on **one output partition
(the default `radial_res = 10000`) and four fixed canonical central densities**, comparing
the old node-quantized locator with the corrected event. No document derives, measures or
claims `1e-9` for any other partition. Its companion radius clause (+5.0 ± 0.5 m at
1.6 M☉) is itself a 10000-only statement: the same star's physical shift at 5000 is
+9.9 m, so applying V7a to a 5000 row would fail its radius clause before its mass clause —
V7a was never a resolution-independent contract. The §10 impact map has no entry for
nontruncated nondefault-resolution rows; the migration procedure
(`docs/validation/TOV_SURFACE_CONTRACT_DERIVATION.md:232`) refers to "the impact map
above" and sends class-A rows to §13 (retain 2500, extend the set, re-measure). The
TOV-SURF-M brief's "fixed-εc `|ΔM/M| <= 1e-9`" restated V7a's number without its
resolution qualifier (§6, §17). **Adjudication: the `1e-9` old→new bound is V7a's,
i.e. default-resolution, canonical, fixed-εc; it is not a universal migration bound.**

### 19.3 The three quantities, kept apart

| | comparison | bound | status |
|---|---|---|---|
| **A** | corrected vs corrected, across output partitions (ADR-0009 Q6/Q12, V4) | `M ≤ 1e-9`, `R_* ≤ 1e-8` | validated on both EOSs (`TOV_SURFACE_IMPLEMENTATION.md:428`); **unchanged** |
| **B** | old locator vs corrected locator, same inherited εc, default 10000 (V7a) | `|ΔM/M| ≤ 1e-9`, radius map ± 0.5 m | passed for 1.0/1.4/1.6/2.0 M☉ (`TOV_SURFACE_IMPLEMENTATION.md:379`); **unchanged** |
| **C** | historical artifact row (old output-grid-quantized surface) vs corrected event, at **each** historical resolution | none predeclared; adjudicated here | C at 10000 coincides with B; C at other resolutions is resolution-dependent by construction (§19.4) |

### 19.4 Physical derivation and scratch verification

**Structure of the old rows.** The old and corrected right-hand sides are identical for
every trial state with `p ≥ p_cut` (`TOVSolver.cpp:1490`); the old solve died only when a
trial stage of the step *toward the first node below the cutoff* dipped under `p_cut`
(`TOV_RADIAL_RES_2500_AUDIT.md:233-248`). Every accepted step up to the last node above
the cutoff is therefore bit-identical between the old and corrected solvers, and a
nontruncated old star is **the corrected star with its final (event) node removed**. The
artifacts already say so: old `N_profile` = corrected `N_profile` − 1 at every
nontruncated resolution (1319/1320, 2635/2636, 5268/5269, 10535/10536), and the
omitted thickness `R_* − R_old` lies inside one outer output step
`Δr_outer(N) = (r_max − r_min)/N` (`TOVSolver.cpp:2532`, `:2600-2601`; the artifact
column `dr_max_km`).

**Physical relation.** The old→new mass change of such a row is the mass of the single
omitted layer of one solution, `ΔM = ∫_{R_old}^{R_*} 4πr²ε dr = m(R_*) − m(R_old)`. With
`dm/dp = f_m/f_p` (the code's own radial right-hand side, the same relation the canonical
locator integrates, `TOVSolver.cpp:2479`) and `p ≪ ε`, the leading order is
`ΔM ≈ 4πR⁴(1 − 2M/R)(p_old,last − p_cut)/M` (G = c = 1): the layer's weight equals its
pressure deficit times area over gravity. The density cancels at leading order, so the tail
is set by how far above `p_cut` the old last node sat — a **phase** of the old partition,
bounded by the weight of one full outer step, `E(N) = m(R_*) − m(R_* − Δr_outer(N))`.

**Scratch probe** (`/Users/keeper/.codex/diagnostics/tov-surf-ma-20260903/tail_check.cpp`,
output `tail_check.out`, `SHA256SUMS.json`): linked against the corrected
`build/libCompactStar.a`, it calls the production primitive at the inherited εc
`731253342677476.125 g/cm³` for `radial_res ∈ {1250, 2500, 5000, 10000, 20000, 40000, 80000}`,
takes the penultimate node `k` and the event node `*`, evaluates `TOVSolver::ODE` at both,
and builds `NStar` on the full profile and on the profile without its event node. No
production, test, CMake or baseline file was touched.

*Table 19-1 — old artifact row vs corrected profile without its event node (relative differences).*

| `radial_res` | `R` | `M` | `B` (NStar on nodes 0…k) | `R_* − R_old` | `Δr_outer` | fraction | class |
|---|---|---|---|---|---|---|---|
| 2500 | `4.34e-2` | `1.44e-3` | `1.74e-3` | 569.130 m | 28.000 m | 20.3 steps | **truncated star** (class A recovery, §6) |
| 5000 | `1.2e-14` | `2.1e-13` | `2.1e-13` | 9.900 m | 14.000 m | 0.707 | surface tail |
| 10000 | `3.6e-13` | `2.3e-13` | `1.3e-13` | 5.035 m | 7.000 m | 0.719 | surface tail (= V7a star) |
| 20000 | `3.4e-14` | `4.0e-14` | `1.6e-13` | 2.340 m | 3.500 m | 0.669 | surface tail |
| 40000 | `3.2e-13` | `2.6e-14` | `2.1e-14` | 0.678 m | 1.750 m | 0.387 | surface tail |

The old 5000/10000/20000/40000 rows are reproduced to 13-digit serialization in `R`, `M`
**and** `B`; the old 2500 row is not a node of the corrected star at any precision and
its omitted interval is twenty outer steps — it is the audited truncated star, not a
surface-tail case, and is excluded from the envelope derivation (§9 of the brief).

*Table 19-2 — omitted-layer mass tail on the corrected solution (`ΔM = m_* − m_k`).*

| `radial_res` | `p_k/p_cut` | `ΔM` (M☉) | `ΔM/M` | actual / trapezoid-in-`p` (code RHS) | actual / leading-order identity | endpoint bracket `[min,max](|dm/dp|)·Δp` | `E(N)/M` | fraction of step |
|---|---|---|---|---|---|---|---|---|
| 1250 | 8.098 | `1.5044e-8` | `9.403e-9` | 0.99870 | 1.00362 | inside (width `1.1e-2`) | `5.01e-8` | 0.551 |
| 2500 | 2.096 | `2.3336e-9` | `1.4585e-9` | 0.99986 | 1.00115 | inside (`3.1e-3`) | `9.54e-9` | 0.326 |
| **5000** | **2.218** | **`2.5927e-9`** | **`1.6205e-9`** | **0.99983** | **1.00125** | **inside (`3.4e-3`)** | **`2.898e-9`** | **0.707** |
| 10000 | 1.530 | `1.1290e-9` | `7.0564e-10` | 0.99996 | 1.00057 | inside (`1.7e-3`) | `1.099e-9` | 0.719 |
| 20000 | 1.225 | `4.8075e-10` | `3.0047e-10` | 0.99999 | 1.00015 | inside (`8.0e-4`) | `4.79e-10` | 0.669 |
| 40000 | 1.062 | `1.3183e-10` | `8.2396e-11` | 1.00000 | 0.99989 | inside (`2.3e-4`) | `2.23e-10` | 0.387 |
| 80000 | 1.006 | `1.3252e-11` | `8.283e-12` | 0.99999 | 0.99978 | inside (`2.4e-5`) | `1.07e-10` | 0.079 |

The 2500/1250 rows are the corrected star's own one-step tails (no historical comparator).
The record's artifact-decimal `1.62002e-9` for 5000 (§6) is the same quantity as the
probe's `1.62047e-9` after 13-digit rounding of `M`.

Findings against the brief's §4 tests:

- **Signs:** `M`, `R_*`, `B` all increase at every nontruncated resolution (the corrected
  star adds the omitted layer; nothing is removed).
- **Consistency with the omitted thickness:** the actual tail equals the trapezoid of the
  code's own `dm/dp` over the omitted pressure interval to `≤ 1.7e-4` of the tail for the
  historical 5000–40000 rows (`≤ 1.3e-3` for the leading-order identity; the thick 1250
  one-step tail reaches `1.3e-3` / `3.6e-3`) and lies inside the parameter-free two-sided
  bracket at all seven resolutions. The **5000 value is quantitatively expected**: its
  layer is 1.97× thicker than the 10000 layer (9.90 vs 5.04 m) *and* denser
  (`ε_k` `3.00e8` vs `2.27e8 g cm⁻³`, `p_k` `2.22` vs `1.53 p_cut`), giving 2.30× the
  10000 tail. No multiplier was fitted; the identity has no free parameter.
- **Convergence:** the envelope `E(N)` falls as `1/N` (`5.0e-8 → 1.1e-10` over
  1250 → 80000); the individual tail is `E(N)` times a phase fraction (measured
  0.08–0.72, scattered, as the audit's step-phase mechanism predicts), so the tail
  vanishes with refinement in envelope, not pointwise. **`E(10000) = 1.099e-9 > 1e-9`**: even
  at the default resolution a surface phase above 0.93 of the outer step would give a
  fixed-εc tail above `1e-9`. V7a's `1e-9` holds for the four canonical stars because their
  phases happen to be ≤ 0.72, not because `1e-9` is a physical bound on the old locator's
  error. A universal `1e-9` would demand of the OLD locator a property that ADR-0009 was
  written to remove and that the audit measured it never had (33 m / 0.24 % scatter,
  `TOV_RADIAL_RES_2500_AUDIT.md:243`); it contradicts the accepted resolution-dependent
  surface error.

*Table 19-3 — baryon tail (INV-14 integrand `4πr² n_B e^λ · 1e54`, INV-13 linear interpolant).*

| `radial_res` | `ΔB/B` (NStar, corrected minus corrected-without-event) | actual / ½(f_k + f_*)·ΔR | `(ΔB/B)/(ΔM/M)` | `e^λ (n_B/ε) (M/B)` at the surface |
|---|---|---|---|---|
| 1250 | `1.1855e-8` | 1.000000 | 1.2608 | 1.1205 |
| 2500 | `1.6615e-9` | 1.000000 | 1.1392 | 1.1205 |
| **5000** | **`1.8510e-9`** | **1.000000** | **1.1423** | **1.1206** |
| 10000 | `7.9511e-10` | 1.000000 | 1.1268 | 1.1206 |
| 20000 | `3.3714e-10` | 1.000000 | 1.1220 | 1.1206 |
| 40000 | `9.2341e-11` | 1.000000 | 1.1207 | 1.1206 |
| 80000 | `9.2814e-12` | 1.000000 | 1.1206 | 1.1206 |

`ΔB` is exactly the one added trapezoid interval of the sampled integral
(`CompactStar/Core/src/NStar.cpp:310-336`; `dependencies/include/Zaki/Vector/DataSet.hpp:572`,
`gsl_interp_linear`), and the reconstructed old `B` equals the artifact to `2e-13`. The
ratio `(ΔB/B)/(ΔM/M)` converges to `e^λ(R_*)·(n_B/ε)_surface·(M/B) = 1.1206`; its excess at
coarse sampling (1.9 % at 5000) is the trapezoid's endpoint weighting of a curved integrand
over a thicker interval — a property of the sampled `B`, absent from `M`, which the event
locator integrates exactly. **The larger 5000 `ΔB/B` is likewise the expected omitted
low-density tail at coarser historical sampling.** No baryon-current physics and no
INV-14 element is touched.

### 19.5 Adjudication (the owner's §12 prior, tested)

The prior is confirmed on every point, with one sharpening: `E(10000) > 1e-9` shows the
V7a number is not a physical bound even at the default resolution (§19.4). Decisions:

- **Q1 — scope of the `1e-9` V7a bound.** V7a only: default `radial_res = 10000`, the four
  inherited canonical εc, old locator vs corrected locator. It is not a universal bound on
  historical output-grid resolutions. ADR-0009's corrected-vs-corrected bounds (A) and the
  V7a result (B) are unchanged and not weakened.
- **Q2 — acceptance rule for nondefault, nontruncated fixed-εc rows.** A historical row with
  unchanged εc is a *surface-tail row* and migrates when all of the following hold:
  - **R1 (screen)** `0 < R_* − R_old ≤ Δr_outer(N)`, with `Δr_outer(N) = (r_max − r_min)/N`
    (the artifact's `dr_max_km`);
  - **R2 (identity)** the old `(R, M)` coincide with a node of the corrected profile at the
    same εc and `N`, and the old sampled integral(s) (`B`) are reproduced by the corrected
    profile without its event node, both to serialization precision (`≤ 1e-12` relative);
  - **R3 (mass tail)** `M_new − M_old` lies inside the endpoint bracket
    `[min, max](|dm/dp|_k, |dm/dp|_*)·(p_k − p_cut)` and equals the trapezoid-in-`p` value
    to better than the bracket width, with positive sign;
  - **R4 (baryon tail)** `B_new − B_old = ½(f_k + f_*)(R_* − R_old)` of the INV-14 integrand
    on the INV-13 interpolant (ratio `1 ± 1e-5`), and `(ΔB/B)/(ΔM/M)` is explained by
    `e^λ(n_B/ε)(M/B)` plus the trapezoid endpoint term;
  - **R5 (validation of the new numbers)** corrected-vs-corrected partition invariance of the
    migrated rows within ADR-0009 V4 (`M ≤ 1e-9`, `R_* ≤ 1e-8`);
  - **R6 (producer)** byte-identical independent repeat and `≤ 1e-12` agreement with the
    validated dry run (§3–§4).
  No universal old→new mass number is imposed. Rows meeting R1–R2 have their delta
  *defined* by the omitted layer; R3–R4 verify that physics; R5–R6 validate the candidate.
  A row failing R1 or R2 is not a surface-tail row and needs its own classification
  (the 2500 rows: class-A truncation recovery, §6).
- **Q3 — envelope.** Required as a *derived* quantity, never a free number:
  `E(N) = m_new(R_*) − m_new(R_* − Δr_outer(N))`, the mass of one full outer-step surface
  layer of the corrected solution; every R1–R2 row satisfies `ΔM ≤ E(N)` automatically
  because its layer is a sub-layer. Reviewer's closed form:
  `ΔM ≈ 4πR⁴(1 − 2M/R)(p_old,last − p_cut)/M`. Values: Table 19-2. Combination of the
  brief's options 1 + 2; option 3 is rejected for the reason in §19.4.
- **Q4 — `B`.** Judged by R4. The `baryon_number_cmf` `1e-15` comparison remains a same-build
  repeatability tolerance on the migrated artifact, not a bound on the old→new tail.
- **Q5 — the 5000 row.** **Scientifically acceptable.** R1: 9.900 m ≤ 14.000 m. R2:
  `1.2e-14 / 2.1e-13 / 2.1e-13`. R3: inside the bracket, trapezoid ratio 0.99983, positive.
  R4: ratio 1.000000, `(ΔB/B)/(ΔM/M) = 1.1423` explained. R5: corrected 5000 `M`, `R_*` within
  `4.5e-12` and `1.6e-11` of the corrected 10000 values (both-EOS V4 sweep
  `TOV_SURFACE_IMPLEMENTATION.md:428`). R6: byte-identical; `3.69e-13`. The 20000 and
  40000 rows pass the same rule (they also happen to be below `1e-9`, which is not the
  criterion). The `B_fixed_mass` 5000/10000/20000/40000 rows carry the identical numbers
  (same selected εc) and are accepted by the same rule; as target-mass rows they are also
  V7b rows.
- **Q6 — resume from the preserved candidates.** **Yes, without regeneration.** run1/run2 are
  byte-identical, were produced by the authenticated HEAD build of the validated solver
  (`evidence/producer-binary-hashes.json`, `evidence/production-hashes.json`), and match the
  validated dry run to `3.69e-13`; regenerating them would add nothing. The migration task's
  own post-copy reproduction and complete non-Phase-4D suites remain its gates.
- **Q7 — new production/TOV validation before migration.** **None.** No production change;
  ADR-0009 remains ACCEPTED / SOURCE CONFORMED / NUMERICALLY VALIDATED; this adjudication
  analysed the validated solver's outputs only.
- **Q8 — mechanism.** An artifact-migration record clarification (this section) plus a
  post-validation clarification note in ADR-0009 — the same mechanism as the V7 clarification
  (commit `03b7d58`). Not a new ADR: the accepted Decision, the scientific bounds, the
  invariants, `PressureCutoff`, and production are untouched; what is clarified is the scope
  of one impact-experiment number and the gate for rows the §10 impact map never covered.
  The §3.1 exception is not involved.

### 19.6 Preserved special cases

- **2500 rows:** truncated different star (R1 fails by a factor 20; R2 fails at `1e-3`);
  class-A defect recovery justified in §6 and `docs/validation/TOV_RADIAL_RES_2500_AUDIT.md`;
  not used in deriving the envelope.
- **Target-mass rows (V7b):** unchanged. `SolveToProfile` may select a different εc under
  its existing absolute `1e-4 M☉` tolerance; the 2.0-M☉ `I` and `B` changes remain
  target-root movement (§11, §12) and no fixed-εc envelope applies to changed-root rows.
- **Thermal columns and trajectories:** governed by the surface-impact scale already
  recorded (§7, §10); not contested by the stopped gate and not re-adjudicated here.

### 19.7 What this adjudication did not do

No artifact copied. No production, test, CMake or baseline change. No ADR Decision change.
No corrected Phase-4D revalidation, no monopole baseline, no Phase 5. No test suite was
re-run: the repository change is documentation only (this section, the top status line, a
§17 pointer, and the ADR-0009 clarification note), and the seven baseline hashes were
re-verified unchanged.

**Exactly one recommended next task:** Codex resumes TOV-SURF-M from the preserved
run1/run2 byte-identical candidates, performs the authorized copy, updates the four legacy
regressions, runs the complete non-Phase-4D suites, records new hashes, and pushes one
atomic migration commit.
