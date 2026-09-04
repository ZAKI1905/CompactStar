# Phase 4D-BL — first verified Hartle monopole regression baseline

> **FIRST VERIFIED HARTLE MONOPOLE BASELINE ESTABLISHED — PHASE-4E CLOSEOUT READY**
> Scientific status remains **HARTLE O(OMEGA^2) MONOPOLE RESPONSE VERIFIED**.
> Production science and all seven prior artifacts are unchanged.

## 1. Entry authentication and chronology

Worktree `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-monopole-validation`,
branch `validation/hartle-monopole-revalidation`; starting SHA
`eccbfa6951ec7ed489e7dfde1fc93c2759d57e2a`, clean, upstream and live remote equal.
Master `df859b5a73c4cac0c115f240744d89ce9f830b8d`; **24 ahead / 0 behind**.

The verified state precedes this baseline commit:

| Commit | Role |
|---|---|
| `218c9aa302e5641084109746e4ff13fc36c56df5` | ADR-0009 artifact migration |
| `42b34ac17ffb3aecd4470d5df54d448f457c30db` | corrected independent Phase-4D characterization |
| `2d85ad5178af33ad9605594d3f49d27a0d2bf655` | owner Validation-D adjudication |
| `c87fb09508f5db03f89bcb77fd97fead4c784bd5` | ADR-0008 convergence clarification / VERIFIED sync |
| `eccbfa6951ec7ed489e7dfde1fc93c2759d57e2a` | documentation correction; final VERIFIED entry authority |

Each was authenticated by `git merge-base --is-ancestor <commit> HEAD`; none is rewritten.
ADR-0008's post-validation status is `CORRECTION IMPLEMENTED — INDEPENDENTLY VERIFIED
(Phase 4D-RV + 4D-DA) — FIRST MONOPOLE BASELINE NEXT`
(`docs/adr/ADR-0008-measure-complete-eos-energy-density-source.md:312`).
The scientific VERIFIED disposition is `docs/validation/PHASE4D_CORRECTED_MONOPOLE_REVALIDATION.md:611`.

Before code or baseline work, both unmodified configurations built and their complete suites
passed serially: **44/44 full, 674.02 s; 22/22 self-contained, 87.85 s; zero skips/failures**.
The raw CTest inventories, JUnit results and command logs are in the evidence directory (§11).

**Proof of first baseline:** filesystem absence and
`git cat-file -e eccbfa6:tests/baselines/hartle_monopole_dscmf1_debug.tsv` returning nonzero;
`git ls-tree -r --name-only eccbfa6 tests/baselines` contained exactly the seven artifacts in §9.
No monopole regression was registered. The filename search found only test sources, including
`tests/rotation/hartle_monopole_reference.hpp`, which is independent solver code, not a baseline.
Historical candidate/debug tables remain diagnostics, never retrospectively adopted as authority
(`docs/validation/PHASE4D_CORRECTED_MONOPOLE_REVALIDATION.md:59`, `:158`).

## 2. Authority and limits of this artifact

Change classes: **generated artifact, engineering/test, dependency/build registration, documentation**
(`GOVERNANCE.md:47`). Production scientific behavior is unchanged. The owner authorized creation
under ADR-0007/0008 and the completed ADR-0009 background contract; this discharges the deferred
baseline chronology of `GOVERNANCE.md:109` (§3.1 condition 7), rather than performing a new correction.

The artifact detects changes in verified production output. **It is not a physics oracle, does
not establish equation correctness, and does not elevate scientific status.** Physical authority
remains the analytic constant-density solution, independent `(m0,h0)` solver, independent
Stieltjes/profile and EOS-knot oracles, homogeneous sequence derivative, C&M 1974, HT68, and
M1–M10 detector evidence in
`docs/validation/PHASE4D_CORRECTED_MONOPOLE_REVALIDATION.md:111`, `:134`, `:165`, `:279`, `:305`, `:611`.
The superseded 4D failure and CHARACTERIZED chronology remain intact. No new detector campaign,
scientific adjudication, representation change or production EOS-knot partition was performed.
Some production API comments retain their pre-revalidation status wording; the later accepted
ADR-0008 status governs. Production files, including those comments, were left byte-identical.

## 3. Canonical producer and star definition

Producer/regression: `tests/rotation/hartle_monopole_regression.cpp:101`.
Registry: `tests/CMakeLists.txt:266`; labels
`rotation;hartle;monopole;scientific;regression;external-data`, without an independent-oracle label.

Targets **1.0, 1.4, 1.6, 2.0 M☉**, production/default **radial_res=10000**
(`CompactStar/Core/TOVSolver.hpp:508`), via `TOVSolver::SolveToProfile`.
No frozen historical central densities, no search override, no alternative solver.
The existing absolute `mass_tol=1e-4 M☉` and first-acceptable stopping policy remain
(`CompactStar/Core/src/TOVSolver.cpp:2742`, `:2771`). These are **V7b target-mass backgrounds**:
selected central density can change under a future governed workflow change; each monopole solve
then computes the fixed-central-density response of that returned background.

The exact EOS is `DS-CMF-1-with-crust/DS(CMF)-1_with_crust.eos`, external root
`/Users/keeper/Documents/CompactStar/data/compose`, SHA256
`5747dd73256c0c28bc56be337cbb96d0918a54bc9ed9fc40984c5befd47ae5dd`, authenticated against
`docs/validation/TOV_REFERENCE.md:132`. Debug build, Apple clang 17.0.0
(clang-1700.6.4.2), CMake 4.2.1, arm64 macOS; no cross-build equivalence claimed.

For each star the producer requires a positive return and `SURFACE_REACHED`, exact final
`TOVPoint::p == PressureCutoff()`, and target residual `<1e-4 M☉`. A test-only subclass exposes
the existing protected cutoff for inspection; it introduces no cutoff formula or locator.
`NStar` construction is required to leave the monopole response absent; its verified first-order
response must exist. Only the explicit `ComputeHartleMonopoleResponse()` supplies coefficients
(`tests/rotation/hartle_monopole_regression.cpp:112`, `:124`).

## 4. Schema, units and normalization audit

Schema **1**, 14 named tab-separated columns, four rows. Header comments carry the verification
anchor, EOS identity, target/default-resolution policy, finite-surface semantics and normalization.
17 significant digits (`max_digits10`), classic locale, lossless round-trip doubles.
No raw first-order seed, private state, l=2 or Phase-5 coefficient is stored.

| Column, in serialized order | Unit | Source |
|---|---|---|
| `target_M_Msun` | M☉ | requested target passed to `SolveToProfile` |
| `epsilon_c_gcm3` | g cm⁻³ | `const NStar::GetSequence().ec`; central energy density in mass-equivalent cgs units |
| `achieved_M_Msun` | M☉ | `const NStar::GetSequence().m` |
| `R_km` | km | `HartleMonopoleResponse::R_surface` |
| `p_cut_dyn_cm2` | dyn cm⁻² | existing `TOVSolver::PressureCutoff()`; exactly equal to final `TOVPoint::p` |
| `epsilon_surface_gcm3` | g cm⁻³ | final `TOVPoint::e`, energy density in mass-equivalent cgs units |
| `I_km3` | km³ | `HartleFirstOrderResponse::I`, required equal to monopole response `I` |
| `m0_surface_over_Omega2_km3` | km³ | last `m0_over_Omega2` node, interior `R_*⁻` value |
| `p0star_surface_over_Omega2_km2` | km² | last `p0star_over_Omega2` node |
| `delta_p0_surface_over_Omega2_dimensionless` | 1 | last existing public `delta_p0_over_Omega2` node |
| `xi0_surface_over_Omega2_km3` | km³ | `surface_xi0_over_Omega2`, equal to last `xi0_over_Omega2` node |
| `surface_shell_mass_over_Omega2_km3` | km³ | `surface_shell_mass_over_Omega2` |
| `I2_over_R3_km3` | km³ | test-side arithmetic `I*I/(R*R*R)` using the verified first-order `I` |
| `deltaM_over_Omega2_km3` | km³ | `deltaM_over_Omega2` |

Unit authority: `CompactStar/Core/RotationSolver.hpp:283` through `:311` and ADR-0007 definitions
(`docs/adr/ADR-0007-hartle-second-order-monopole-response.md:84`).
`Omega_geom = Omega_phys/c`, where physical Ω is rad s⁻¹ and `c` is km s⁻¹, hence geometric Ω
is km⁻¹ (`CompactStar/AngularVelocity.hpp:34`). Thus division by Ω_geom² adds km²:
`m0`/`xi0`/shell/`deltaM` [km] become km³, dimensionless `p0star` becomes km², and `delta_p0`
[km⁻²] becomes dimensionless. **No coefficient is labelled per physical s⁻².**
`I²/R³` has km⁶/km³ = km³.

`TOVPoint` carries `p,e` in cgs (`CompactStar/Core/TOVSolver.hpp:332`);
`NStar::BuildFromTOV` converts pressure and density to km⁻² (`CompactStar/Core/src/NStar.cpp:189`).
Its sequence central density is converted back to g cm⁻³ (`:344`); its mass is M☉ (`:356`).
The shell check uses the **geometric profile density**, not the cgs value stored for external audit.
The existing public `delta_p0_over_Omega2` is retained, not invented: `(eps+p)*p0star_over_Omega2`
(`CompactStar/Core/RotationSolver.hpp:289`).

All surface columns refer to the **finite `p=p_cut` event `R_*`**, never the vacuum `p=0` surface.
For all four rows `p_cut=3.3518848999999998e25 dyn cm⁻²` and
`epsilon_surface=165880787 g cm⁻³` (`tests/baselines/hartle_monopole_dscmf1_debug.tsv:12`).

## 5. Four canonical rows and generation sanity

The installed artifact contains every field at full serialization precision. The following tables
are transcribed from the producer bytes; no numbers were used as generation targets.

| target_M_Msun | epsilon_c_gcm3 | achieved_M_Msun | R_km | I_km3 |
|---|---|---|---|---|
| 1 | 454550405078491.75 | 1.0000438097285522 | 13.427522155261627 | 86.99575838204774 |
| 1.3999999999999999 | 616488270506054.5 | 1.4000217832812236 | 13.546350432736507 | 135.61613061268929 |
| 1.6000000000000001 | 731253342677476.12 | 1.5999758353006395 | 13.473358193949251 | 159.587141588083 |
| 2 | 1297985389788764.5 | 1.9999985861663816 | 12.716477747345323 | 193.72211615013879 |

| target_M_Msun | m0_surface_over_Omega2_km3 | p0star_surface_over_Omega2_km2 | delta_p0_surface_over_Omega2_dimensionless | xi0_surface_over_Omega2_km3 |
|---|---|---|---|---|
| 1 | 992.13163319988598 | 34.697655682172176 | 4.2752170537193241e-09 | 3304.3867974489731 |
| 1.3999999999999999 | 933.52273752584563 | 33.850377637372169 | 4.1708210224845076e-09 | 2087.4301753782556 |
| 1.6000000000000001 | 855.452784276741 | 32.144579465178531 | 3.9606437844957026e-09 | 1603.5364346092667 |
| 2 | 555.61886931325898 | 23.932009736539932 | 2.9487449265341039e-09 | 701.68466031308253 |

| target_M_Msun | surface_shell_mass_over_Omega2_km3 | I2_over_R3_km3 | deltaM_over_Omega2_km3 |
|---|---|---|---|
| 1 | 0.0009222590010069141 | 3.1261407235308112 | 995.25869618241779 |
| 1.3999999999999999 | 0.00059296196755143278 | 7.3987061522737028 | 940.92203664008684 |
| 1.6000000000000001 | 0.00045060998685392392 | 10.412816728421443 | 865.86605161514933 |
| 2 | 0.00017564909748357254 | 18.249791288764293 | 573.8688362511208 |

All 44 comparable fields (11 per row) agree with the already committed Phase-4D-RV table
at its printed rounding precision (`docs/validation/PHASE4D_CORRECTED_MONOPOLE_REVALIDATION.md:154`).
All signs are positive and every serialized value finite. At 1.6 M☉,
`R=13.473358193949251 km`, `I=159.587141588083 km³`, `deltaM_hat=865.86605161514933 km³`;
the verified table prints `13.473358194`, `159.5871416`, `865.86605162` respectively.
This records the accepted finite-resolution response; the separate continuum diagnostic does
not supply the baseline or its tolerance. No unexplained target-root difference was found.

The four absolute target residuals are `4.38097285522e-5`, `2.17832812237e-5`,
`2.41646993606e-5`, `1.41383361840e-6 M☉`; all satisfy the unchanged V7b contract.
Producer logs include the existing 18 rejected low-density coarse samples per target
(72 logged errors per four-star run). These samples are excluded by the fail-closed production
scan; none supplies a baseline row. Each returned star is complete. These are expected
scan diagnostics, not failed returned profiles or regression failures
(`CompactStar/Core/src/TOVSolver.cpp:2674`, `:2791`).

## 6. Decomposition and provenance

For every row, independently form
`deltaM_hat = m0_hat(R_*^-) + surface_shell_hat + I²/R_*³`.
The C++ arithmetic residual is zero on this build. Re-evaluating the serialized decimals with
60-digit arithmetic gives signed residuals (km³):
`-8.1141e-15`, `-4.423278e-14`, `+3.307608e-14`, `+4.342746e-14`;
maximum relative residual `7.568e-17`. Independently recomputing `I²/R³` from serialized
`I,R` differs by at most `1.070e-15 km³`. The shell and Eulerian pressure identities also pass.
The arithmetic acceptance bound `1e-14` is inherited from the existing conformance identity
(`tests/rotation/hartle_monopole_cmf.cpp:57`, `:180`); it is not the baseline comparison tolerance.

The producer reads through `const NStar&`, calls `MonopoleResponse()`, requires validity through
`MatchesSource(&profile, profile.Version())`, and checks again immediately before storing the row
(`tests/rotation/hartle_monopole_regression.cpp:128`, `:167`). This reuses the existing
identity/version authority, without a second cache implementation or production API mutation.
The production accessor refuses stale responses (`CompactStar/Core/src/RotationSolver.cpp:1335`).
Existing `hartle_monopole_contract` C6/C7 proves version-bump rejection, refusal of stale
materialization, and explicit recomputation (`tests/rotation/hartle_monopole_contract.cpp:368`).
That test passed in both starting suites; final results are recorded in §8.

## 7. Repeatability, installation and regression precision

Three separate producer processes used fresh paths in the external evidence directory (§11).
The exact commands are in `producer-runs.json` and `final-runs.json`; equivalent invocation:

```sh
build/tests/hartle_monopole_regression /Users/keeper/Documents/CompactStar/data/compose --emit /absolute/external/fresh.tsv
build/tests/hartle_monopole_regression /Users/keeper/Documents/CompactStar/data/compose
```

Run1 and run2 were generated before baseline creation. They are byte-identical: **2,145 bytes**.
After the sanity and arithmetic audit, exact run1 bytes were copied once into
`tests/baselines/hartle_monopole_dscmf1_debug.tsv`; its immediate hash and bytes matched run1.
Run3, generated after installation, equals run1, run2 and installed bytes.
Every one has SHA256:

`bd49e5a091ebcc59f7c4899422200181d4e71ecf552284840454d01aac4b8d52`

Measured same-build repeatability is zero differing bytes/fields. **The regression requires exact
serialized-byte equality**, including schema/header, using the repository's 17-significant-digit
round-trip convention (`tests/core/baryon_number_cmf.cpp:130`). There is no adjustable numerical
regression tolerance and no allowance derived from scientific validation bounds, D′, the old 4.6%
deficit, migration shifts, or the continuum estimate. A future toolchain discrepancy must be
investigated, not silently tolerated. Independent scientific bounds remain untouched.

The regression only reads its expected file; only explicit emit mode writes, and emit rejects
existing paths and any resolved path under `tests/baselines` (`tests/rotation/hartle_monopole_regression.cpp:84`).
Two expected negative checks passed: emit into the baseline directory is refused before solving;
a one-last-digit corruption of an external comparison copy is rejected after recomputation.
The corrupt copy is clearly labelled scratch and is never an authoritative candidate.

## 8. Complete test inventory and final integrity

| Configuration | Before baseline | Final | Skipped / excluded / failed |
|---|---|---|---|
| Full Debug, external data | 44/44, 674.02 s | **45/45, 679.30 s** | **0 / 0 / 0** |
| Self-contained Debug | 22/22, 87.85 s | **22/22, 91.38 s** | **0 / 0 / 0** |

`cmake --build build -j6` and `cmake --build build-selfcontained -j6` passed.
CTest used `--test-dir <build> -j1 --output-on-failure --output-junit <external-path>`.
`hartle_monopole_regression` is green; the three independent monopole physics/published tests
and existing provenance/staleness tests are included and green. The only deliberately failing
executions were the two negative checks in §7. No unexpected test failure occurred.

`git diff eccbfa6 -- CompactStar/` is empty; all 197 production source hashes still match entry.
The seven previous baseline hashes match §9; run1/run2/run3/installed bytes are identical.
`git diff --check` is clean. Only the new producer, test registration, one new baseline and
current-state documentation are included; no scratch files or other scientific artifacts.
The atomic commit is `test: establish Hartle monopole baseline`, after the verified chronology
in §1. Final clean-tree/live-remote equality is recorded in external `delivery.json` (§11).

Build/test suites run serially, including all independent monopole analytic/CMF/published tests,
first-order, TOV, thermal, EOS and baryon tests. No suite runs concurrently with the other.
No exclusions, expected test failures or skips are permitted. The new regression is external-data
only; it does not enter the self-contained inventory. Existing scientific tolerances and tests
are unchanged. No new scientific adjudication or detector mutation campaign is part of this task.

## 9. Current eight-artifact authority

The seven ADR-0009 migrated artifacts were authenticated before any edit and remain byte-identical.
The eighth is the first monopole regression. Historical seven-hash tables retain their original
experiment-time role and values; this table is the current complete registry after Phase 4D-BL.
Hashes below are recomputed from repository bytes, not supplied as expected values by the producer.

| Artifact under `tests/baselines/` | Current SHA256 |
|---|---|
| `passive_cooling_cmf_1p6_debug.tsv` | `afcad5d078fe458bdb441278ce56caa2d025becc6b489d00e9baad91a101c0de` |
| `tov_dscmf1_reference.tsv` | `3d9af9129a6a4ffde9e0f8c5507a160f968a861c0cf9f3b089cceecab86b701a` |
| `grid_convergence_cmf_1p6_debug.tsv` | `2c68e2f7e871192e00322bcc12b6d8b13eb7fdbda4a6d8d7050f26f0271ce5eb` |
| `grid_convergence_cmf_1p6_trajectory.tsv` | `7c84557742ec0ec747756d118b2837fb54095a94c10194d4ccc2d0a778f5b04f` |
| `hartle_I_dscmf1_debug.tsv` | `a21d4c3f6d89322cb4cca7a073e0e142f180e529332d704b7b16a145eae741c9` |
| `baryon_number_dscmf1_reference.tsv` | `7b036942f2ae599ace3b2fc8b9a7d91f6d11b5899ab7e8a88d2bb4ee6686493b` |
| `tov_path_equivalence_dscmf1.tsv` | `5c0f4b3bdb70921f8f2a869af10edc4d8f5ae3963a9d150e11ef859d21e1c678` |
| `hartle_monopole_dscmf1_debug.tsv` | `bd49e5a091ebcc59f7c4899422200181d4e71ecf552284840454d01aac4b8d52` |

Producer for the eighth artifact: `tests/rotation/hartle_monopole_regression.cpp`; schema 1;
verified source anchor `eccbfa6951ec7ed489e7dfde1fc93c2759d57e2a`.
Prior seven producer/migration provenance remains `docs/validation/TOV_SURFACE_ARTIFACT_MIGRATION.md:889`.

## 10. Status and Phase-4E dependency

**FIRST VERIFIED HARTLE MONOPOLE BASELINE ESTABLISHED — PHASE-4E CLOSEOUT READY**

Scientific status remains **HARTLE O(OMEGA^2) MONOPOLE RESPONSE VERIFIED**. INV-08 gains a trusted
regression artifact, not a new physical claim. `l=2` is not implemented or validated and no
high-spin accuracy is claimed. INV-09 remains unresolved: baryon-conserving reduction,
`A_i/B_i/Z_i`, and particle-number measure completeness (ADR-0008 Q11) are not supplied by this
baseline. No accepted ADR-0007/0008/0009 Decision changes.

The chronology is correction implemented → independent revalidation → owner convergence
adjudication → VERIFIED → first monopole baseline. This discharges GOVERNANCE §3.1 condition 7.
The **scientific blocker for Phase 4E is cleared**: the complete final suites in §8 pass.
Phase-4E closeout/interface ratification is not performed here; Phase 5 is not begun.

**Exactly one next action:** perform only the smallest Phase-4E / Phase-4 closeout needed to ratify
the verified `HartleMonopoleResponse` fields as the Phase-5 structural interface and close Phase 4.
Do not begin Phase 5 automatically.

## 11. External reproducibility evidence

`/Users/keeper/.codex/diagnostics/phase4d-bl-20260904/` contains starting authentication
(197 production-file hashes, seven artifact hashes, ancestry), EOS/toolchain authentication,
initial/final CTest JSON inventories, JUnit and complete logs, all three exact producer outputs,
producer command/runtime records, 44-field sanity comparison and decimal decomposition audit,
negative-check records, installation record, and final integrity evidence. These scratch files
are not committed. The final evidence manifest authenticates them individually.
