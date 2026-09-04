# Phase 4E — rotation correctness closeout and structural-interface ratification

> **PHASE 4 ROTATION CORRECTNESS COMPLETE — PHASE-5 STRUCTURAL INTERFACE RATIFIED**
>
> **PHASE 5 NOT YET BEGUN**
>
> Scientific status remains **HARTLE O(OMEGA^2) MONOPOLE RESPONSE VERIFIED**.
> This closeout records the owner's Phase-4E ratification of existing verified objects;
> it adds no scientific derivation, numerical method, API, or downstream implementation.

## 1. Entry authority and change class

Starting SHA: `2e2f01682705edd95f81992f149279c9d00a4602`, branch
`validation/hartle-monopole-revalidation`, worktree
`/Users/keeper/Documents/CompactStar/worktrees/CompactStar-monopole-validation`.
The tree was clean, local HEAD equalled upstream and live remote, and ancestry was
**25 ahead / 0 behind** master `df859b5a73c4cac0c115f240744d89ce9f830b8d`.
Other worktrees were inventoried; this branch contains the prior rotation and surface
work and is the owner-designated authority. No history is rewritten.
Entry commands and file hashes: external `starting-authentication.json` and `entry-git.txt` (§10).

**Change class: documentation**, including the explicitly permitted source-comment correction.
This is the owner's closeout instruction applied to the accepted ADR-0006/0007/0008/0009
contracts, not a new structural abstraction or scientific adjudication. Authority:
`GOVERNANCE.md:43`, `GOVERNANCE.md:225`, `AGENTS.md:42`,
`docs/adr/ADR-0006-hartle-first-order-physical-normalization.md:243`,
`docs/adr/ADR-0007-hartle-second-order-monopole-response.md:418`,
`docs/adr/ADR-0008-measure-complete-eos-energy-density-source.md:61`,
`docs/adr/ADR-0009-tov-surface-event-and-termination.md:47`.
All accepted Decisions remain byte-identical. The following chronology retains the
failed and CHARACTERIZED experiments rather than retroactively relabelling them.

## 2. Complete Phase-4 chronology and accepted scope

| Increment / commits | Completed work and durable evidence |
|---|---|
| Entry `59714a3` | Provenance, normalization and candidate audit; `docs/validation/PHASE4_ROTATION_ENTRY.md:32` |
| 4A `a9c423c` → `36ea3b3` | ADR-0006 acceptance then physical normalization; `docs/validation/PHASE4A_FIRST_ORDER_NORMALIZATION.md:17`, `docs/validation/PHASE4A_FIRST_ORDER_NORMALIZATION.md:86` |
| 4B `bb073c8` | First-order shape and I independently verified; `docs/validation/PHASE4B_FIRST_ORDER_PHYSICS.md:3` |
| 4C-G `bcef5b5` → acceptance `13111de` | Governed l=0 derivation then ADR-0007 acceptance; `docs/validation/PHASE4C_HARTLE2_DERIVATION.md:224`, `docs/validation/PHASE4C_I0_EOS_DERIVATIVE.md:16` |
| 4C-I0 `a97f9c5` | EOS-owned derivative authority; `docs/validation/PHASE4C_I0_EOS_DERIVATIVE.md:67` |
| 4C-I1 `377bc4a` | Atomic candidate replacement, fixed-ε_c response; `docs/validation/PHASE4C_I1_MONOPOLE_IMPLEMENTATION.md:185` |
| Original 4D `246f3f2` | Implementation verified but stepped-crust physical contract failed; no baseline; `docs/validation/PHASE4D_MONOPOLE_VALIDATION.md:3` |
| 4D-RG `8abbab4` → ADR-0008 `cc4bec4` → RI `e2fe0ad` | Measure derivation, owner acceptance, measure-complete source correction; `docs/adr/ADR-0008-measure-complete-eos-energy-density-source.md:61`, `docs/adr/ADR-0008-measure-complete-eos-energy-density-source.md:221` |
| TOV audit `b93be5f` → proposal `8034a9c` → ADR-0009 `7a223d3` | Surface defect characterized and contract accepted; `docs/adr/ADR-0009-tov-surface-event-and-termination.md:18`, `docs/adr/ADR-0009-tov-surface-event-and-termination.md:47` |
| `96c1425` → `03b7d58` → `94165f3` | Historical V7 STOP, owner clarification, validated accepted-state event; `docs/adr/ADR-0009-tov-surface-event-and-termination.md:74`, `docs/adr/ADR-0009-tov-surface-event-and-termination.md:234`, `docs/adr/ADR-0009-tov-surface-event-and-termination.md:256` |
| `816b754` → `218c9aa` | Migration-envelope adjudication then seven-artifact migration; `docs/adr/ADR-0009-tov-surface-event-and-termination.md:269`, `docs/adr/ADR-0009-tov-surface-event-and-termination.md:306` |
| 4D-RV `42b34ac` | Corrected independent characterization, all physics lines and M1–M10 passed; D monotonicity remained open; `docs/validation/PHASE4D_CORRECTED_MONOPOLE_REVALIDATION.md:3`, `docs/validation/PHASE4D_CORRECTED_MONOPOLE_REVALIDATION.md:303` |
| 4D-DA `2d85ad5` → `c87fb09` → `eccbfa6` | Owner D′ convergence adjudication and VERIFIED authority; `docs/validation/PHASE4D_CORRECTED_MONOPOLE_REVALIDATION.md:411`, `docs/validation/PHASE4D_CORRECTED_MONOPOLE_REVALIDATION.md:602` |
| 4D-BL `2e2f016` | First trusted monopole regression baseline after verification; `docs/validation/PHASE4D_MONOPOLE_BASELINE.md:193`, `docs/validation/PHASE4D_MONOPOLE_BASELINE.md:224` |
| 4E, this closeout | Existing response objects ratified as the Phase-5 structural input boundary; §§3–8 below |

Thus 4A (normalization), 4B (independent first order), 4C (governed l=0 implementation),
4D (independent validation, corrected EOS/TOV contracts and first baseline), and
4E (interface ratification) are **COMPLETE for the accepted ordinary-NStar scope**.
The second-order claim is the fixed-central-energy-density l=0 O(Ω²) coefficient
response on complete ADR-0009 TOV backgrounds. It makes no high-spin/O(Ω⁴) accuracy
claim. `l=2` is not implemented or validated; it remains a future independent rotation
extension and is **nonblocking for Phase-4 completion, Phase 5 and BNV thermal work**.
MixedStar rotation is not covered. Scope authority:
`docs/adr/ADR-0007-hartle-second-order-monopole-response.md:428`,
`docs/validation/PHASE4D_CORRECTED_MONOPOLE_REVALIDATION.md:61`.

## 3. Ratified first-order structural interface

The existing `HartleFirstOrderResponse` is the supported normalized first-order input.
Its fields are ratified with their current public names and units:

| Field | Unit / meaning | Source |
|---|---|---|
| `I` | km³, moment of inertia | `CompactStar/Core/RotationSolver.hpp:162` |
| `omega_bar_over_Omega` | dimensionless, s(r) | `CompactStar/Core/RotationSolver.hpp:165` |
| `domega_bar_over_Omega_dr` | km⁻¹, s′(r) | `CompactStar/Core/RotationSolver.hpp:168` |
| `r_grid` | non-owning radius column, km | `CompactStar/Core/RotationSolver.hpp:171` |
| `valid` | boolean, complete usable solve | `CompactStar/Core/RotationSolver.hpp:174` |

This is seed-free, scale-free, a property of the nonrotating background, normalized
to angular velocity at infinity. It needs no supplied physical spin. The supported
verified scope is ordinary NStar. `NStar::RotationResponse()` exposes the owned
response; ordinary construction computes first order without conferring a spin.
Authority: `CompactStar/Core/RotationSolver.hpp:146`, `CompactStar/Core/NStar.hpp:370`,
`docs/adr/ADR-0006-hartle-first-order-physical-normalization.md:256`.
No raw internal omega_bar seed is part of this interface.

## 4. Ratified monopole structural interface and units

The existing `HartleMonopoleResponse` is the supported second-order structural input.
No `Phase5RotationInterface`, downstream wrapper or replacement API is introduced.

| Field | Unit / meaning | Source |
|---|---|---|
| `m0_over_Omega2` | km³, m̂₀(r) | `CompactStar/Core/RotationSolver.hpp:283` |
| `p0star_over_Omega2` | km², p̂₀*(r) | `CompactStar/Core/RotationSolver.hpp:286` |
| `delta_p0_over_Omega2` | dimensionless, Eulerian δp̂₀(r) | `CompactStar/Core/RotationSolver.hpp:289` |
| `xi0_over_Omega2` | km³, isobar displacement ξ̂₀(r) | `CompactStar/Core/RotationSolver.hpp:292` |
| `deltaM_over_Omega2` | km³, m̂₀(R_*⁻) + shell + I²/R_*³ | `CompactStar/Core/RotationSolver.hpp:295` |
| `surface_shell_mass_over_Omega2` | km³, terminal surface atom | `CompactStar/Core/RotationSolver.hpp:298` |
| `surface_xi0_over_Omega2` | km³, displacement coefficient at R_* | `CompactStar/Core/RotationSolver.hpp:301` |
| `R_surface` | km, governed finite p=p_cut event radius R_* | `CompactStar/Core/RotationSolver.hpp:305` |
| `I` | km³, verified first-order moment used by this solve | `CompactStar/Core/RotationSolver.hpp:309` |
| `r_grid` | non-owning radius column, km | `CompactStar/Core/RotationSolver.hpp:313` |
| `source_profile` | non-owning profile identity, `const void*`; no physical unit | `CompactStar/Core/RotationSolver.hpp:317` |
| `source_version` | `std::uint64_t` provenance version; no physical unit | `CompactStar/Core/RotationSolver.hpp:321` |
| `valid` | boolean, complete all-finite atomically published solve | `CompactStar/Core/RotationSolver.hpp:323` |

All hatted second-order coefficients are **per Ω_geom²**, where
`Ω_geom [km⁻¹] = Ω_phys [rad s⁻¹] / c [km s⁻¹]`, normalized at infinity.
They are not coefficients per physical s⁻². The conversion owner is
`CompactStar/AngularVelocity.hpp:34`, `CompactStar/AngularVelocity.hpp:150`.
The current source defines δp̂₀=(ε+p)p̂₀* and ξ̂₀=p̂₀*/ν′;
the unit audit already exists at `docs/validation/PHASE4D_MONOPOLE_BASELINE.md:87`.
R_* is the EOS-floor event surface, **never the exact p=0 physical surface**
(`docs/adr/ADR-0009-tov-surface-event-and-termination.md:54`).

## 5. Runtime consumption, explicit spin and provenance

Phase 5 must consume **HartleFirstOrderResponse + HartleMonopoleResponse**.
It must not independently rerun Hartle ODEs, access raw seed profiles, reconstruct
m₀/p₀*/ξ₀ from private solver state, infer normalization from a physical spin,
or consume regression TSV files as computational inputs. Existing production
`NStar::ComputeHartleMonopoleResponse()` remains the explicit operation that obtains
or refreshes the response; it is never implicit during construction. Existing
`NStar::MonopoleResponse()` is a cheap current-response accessor.
Source: `CompactStar/Core/NStar.hpp:411`, `CompactStar/Core/NStar.hpp:438`.

`PhysicalFirstOrderRotation` and `PhysicalHartleMonopole` remain the explicit-spin
materialization layer via `At(AngularVelocity)`, `NStar::RotationAt` and
`NStar::MonopoleAt`. They scale existing coefficients without another integration.
First-order physical units are Ω_geom km⁻¹, J km², I km³, omega_bar km⁻¹ and its
derivative km⁻². Physical monopole units are m₀/δM/shell/ξ₀ km, p₀* dimensionless,
and δp₀ km⁻². These materializations are for explicit-spin observables, not for
recovering the normalized Phase-5 coefficients.
Source: `CompactStar/Core/RotationSolver.hpp:111`, `CompactStar/Core/RotationSolver.hpp:176`,
`CompactStar/Core/RotationSolver.hpp:200`, `CompactStar/Core/RotationSolver.hpp:331`.

A monopole response is usable only when **valid == true** and
**MatchesSource(&profile, profile.Version()) == true** (the actual API accepts
profile identity as a pointer). `MatchesSource` itself includes the valid check.
A retained response must be rejected after source identity/version changes;
the supported NStar accessor already refuses absent or stale data.
Source: `CompactStar/Core/RotationSolver.hpp:317`, `CompactStar/Core/RotationSolver.hpp:325`,
`CompactStar/Core/src/RotationSolver.cpp:1335`.
First-order references and both non-owning grids follow the owning StarProfile
lifetime; they cannot outlive that profile. Consumers must use the owning star's
current response rather than detached stale copies. No second cache/provenance
mechanism is created. Source: `CompactStar/Core/RotationSolver.hpp:131`,
`CompactStar/Core/RotationSolver.hpp:171`, `CompactStar/Core/NStar.hpp:370`.

The baseline is a **test regression authority only**, established after verification.
Physical authority remains the analytic constant-density solution, independent
(m₀,h₀) solver, Stieltjes/profile and EOS-knot oracles, homogeneous sequence
identity, C&M 1974, HT68 and M1–M10 campaign. This closeout does not elevate the
scientific status through agreement with the baseline.
Evidence: `docs/validation/PHASE4D_MONOPOLE_BASELINE.md:42`,
`docs/validation/PHASE4D_CORRECTED_MONOPOLE_REVALIDATION.md:94`,
`docs/validation/PHASE4D_CORRECTED_MONOPOLE_REVALIDATION.md:303`.

## 6. Fixed-ε_c semantics and invariant disposition

The monopole is the **fixed-central-energy-density particular solution**, with no
homogeneous admixture. It is **not** the baryon-conserving rotating sequence.
In particular `deltaM_over_Omega2` is not a fixed-baryon-number mass correction.
The response alone supplies no B_i, Z_i or sequence-reduction coefficient.
Authority: `CompactStar/Core/RotationSolver.hpp:263`,
`docs/adr/ADR-0007-hartle-second-order-monopole-response.md:428`.

**INV-08: CLOSED / VERIFIED for the governed scope** — GOVERNED (ADR-0007 +
ADR-0008), CONFORMED, INDEPENDENTLY PHYSICALLY VERIFIED, REGRESSION-PROTECTED.
Scope is ordinary NStar, fixed ε_c, l=0, O(Ω²), complete ADR-0009 backgrounds.
No l=2, MixedStar or high-spin/O(Ω⁴) claim is added. Supporting verification and
baseline: `docs/validation/PHASE4D_CORRECTED_MONOPOLE_REVALIDATION.md:602`,
`docs/validation/PHASE4D_MONOPOLE_BASELINE.md:272`.

**INV-09 remains unresolved.** The A_i-side *structural rotation inputs* are now
verified and ratified: m̂₀, p̂₀*/δp̂₀, ξ̂₀, s, s′, I and surface/provenance data.
This does not verify or implement A_i. Remaining work includes measure-complete
particle-number response, composition/species measure dn_i, B_i, Z_i,
baryon-conserving sequence reduction and exact ownership of the homogeneous /
sequence-derivative response. A nodal dn_i/dp column alone is not a permitted
particle-number source. Authority: `docs/adr/ADR-0008-measure-complete-eos-energy-density-source.md:78`,
`docs/adr/ADR-0007-hartle-second-order-monopole-response.md:434`.

## 7. Current eight-artifact authority — unchanged

These are the same eight authorities established by Phase 4D-BL, rehashed before
editing and at exit. Historical seven-artifact tables retain their original hashes.
No artifact is generated, edited or added by this closeout; no ninth artifact exists.
Baseline authority: `docs/validation/PHASE4D_MONOPOLE_BASELINE.md:250`.

| Artifact in tests/baselines | SHA256 |
|---|---|
| `passive_cooling_cmf_1p6_debug.tsv` | `afcad5d078fe458bdb441278ce56caa2d025becc6b489d00e9baad91a101c0de` |
| `tov_dscmf1_reference.tsv` | `3d9af9129a6a4ffde9e0f8c5507a160f968a861c0cf9f3b089cceecab86b701a` |
| `grid_convergence_cmf_1p6_debug.tsv` | `2c68e2f7e871192e00322bcc12b6d8b13eb7fdbda4a6d8d7050f26f0271ce5eb` |
| `grid_convergence_cmf_1p6_trajectory.tsv` | `7c84557742ec0ec747756d118b2837fb54095a94c10194d4ccc2d0a778f5b04f` |
| `hartle_I_dscmf1_debug.tsv` | `a21d4c3f6d89322cb4cca7a073e0e142f180e529332d704b7b16a145eae741c9` |
| `baryon_number_dscmf1_reference.tsv` | `7b036942f2ae599ace3b2fc8b9a7d91f6d11b5899ab7e8a88d2bb4ee6686493b` |
| `tov_path_equivalence_dscmf1.tsv` | `5c0f4b3bdb70921f8f2a869af10edc4d8f5ae3963a9d150e11ef859d21e1c678` |
| `hartle_monopole_dscmf1_debug.tsv` | `bd49e5a091ebcc59f7c4899422200181d4e71ecf552284840454d01aac4b8d52` |

## 8. Phase-5 entry gate

The Phase-4 prerequisite is discharged. Phase 5 **may begin in a subsequent task**
after normal reviewed integration of this completed branch; **PHASE 5 NOT YET BEGUN**.
Available verified/governed inputs are the EOS/StarProfile background, ADR-0001
baryon/species semantics, heat-capacity/thermal infrastructure, first-order response,
fixed-ε_c monopole response, explicit-spin materialization and ADR-0009 complete
TOV surface. Authority: `docs/adr/ADR-0001-species-profile-semantics.md:1`,
`docs/adr/ADR-0002-thermal-heat-capacity-ownership.md:1`, §§3–6 above.

Still to build and govern in Phase 5: standard rotochemical structural particle-number
response, A_i, B_i, Z_i, baryon-conserving reduction, chemical-imbalance evolution,
and modified reaction/neutrino rates. INV-09 and the chemical-convention gate INV-11
remain open; entry into later Phase-5 governance work does not authorize bypassing
those scientific prerequisites. No such design or implementation occurs here.
Authority: `docs/adr/ADR-0008-measure-complete-eos-energy-density-source.md:78`,
`docs/adr/README.md` (current deferred ownership and chemical-convention entries).

## 9. Complete tests and production/artifact protection

Starting complete authority was authenticated from Phase 4D-BL's individually hashed
38-file evidence archive: **45/45 full, 679.30 s; 22/22 self-contained, 91.38 s**;
zero failures, disabled or skipped tests. Actual current CTest inventories agree.
These starting suites were not redundantly rerun before the comment edit.
Evidence: `docs/validation/PHASE4D_MONOPOLE_BASELINE.md:224`, external
`starting-authentication.json:1` (§10).

Both configurations were rebuilt with `cmake --build build -j6` and
`cmake --build build-selfcontained -j6`. Complete final suites run serially with
`ctest --test-dir build -j1 --output-on-failure` and then the same command for
`build-selfcontained`, with JUnit output recorded externally. No suite concurrency,
exclusions, expected failures or scientific tolerance changes are introduced.

| Final configuration | Result | CTest wall time | Failures / skipped / disabled |
|---|---|---|---|
| `build` | **45/45 PASS** | **680.12 s** | **0 / 0 / 0** |
| `build-selfcontained` | **22/22 PASS** | **90.78 s** | **0 / 0 / 0** |

All current tests were included, including the analytic/CMF/published monopole physics,
monopole regression, first-order, TOV, thermal, EOS and baryon tests. Exact inventories
and full outputs: external `build-final-inventory.json:1`,
`build-selfcontained-final-inventory.json:1`, `final-full.log:1`, `final-self.log:1` (§10).

The only production-file difference from the starting SHA is the four-line status
comment at `CompactStar/Core/RotationSolver.hpp:277`. Every other byte of that
header and every other production file is unchanged; there is **no executable-token,
declaration, signature, ABI or scientific behavior change**. Existing tests, producers,
CMake files and all baseline bytes are unchanged. No new scientific test was needed.
The analogous historical warning at `CompactStar/Core/NStar.hpp:446` is recorded
as stale documentation outside this task's one permitted source-comment block;
ADR-0008's verified authority and this closeout govern its current interpretation
(`GOVERNANCE.md:14`). It is not a physics blocker.

Accepted ADR Decisions and all historical validation records remain unchanged;
ADR-0007/0008 receive append-only closeout pointers. No EOS-knot production
partitioning, measure change, l=2, particle-number response, coefficients or Phase-5
work is introduced. Evidence: external `comment-only-proof.json:1`, final integrity
record (§10), and the reviewed diff against the starting SHA.

## 10. Evidence, branch delivery and final disposition

External evidence directory (not committed):
`/Users/keeper/.codex/diagnostics/phase4e-20260904/`.
It contains entry Git/remote authentication, tracked-file and eight-artifact hashes,
authenticated prior-suite evidence, CTest inventories, build logs, complete serial
CTest logs/JUnit, exact comment-only proof, final integrity and delivery records.
`SHA256SUMS.json` authenticates the evidence files individually.

One separate commit: **docs: close Phase 4 rotation correctness**. Non-force push
on `validation/hartle-monopole-revalidation`; no merge and no master fast-forward.
Master remains `df859b5a73c4cac0c115f240744d89ce9f830b8d`. The commit is identified
by the Git history containing this record; post-commit clean-tree and local/upstream/
live-remote equality are recorded in external `delivery.json:1`.

**PHASE 4 ROTATION CORRECTNESS COMPLETE — PHASE-5 STRUCTURAL INTERFACE RATIFIED**

**PHASE 5 NOT YET BEGUN**

**Exactly one next action:** integrate the completed Phase-4 branch into the canonical
development line through the normal reviewed fast-forward/merge process before
starting Phase 5.
