# TOV surface implementation and validation history

**Current TOV-SURF-I-R status: TOV SURFACE EVENT IMPLEMENTED AND VALIDATED — ARTIFACT MIGRATION REQUIRED.**

The current evidence is in §9. Sections 1–7 preserve the prior stop; §8 records the
owner clarification and archived-binary diagnosis. Their historical dispositions
are superseded by §9, not erased.

## Historical TOV-SURF-I stop

**TOV SURFACE IMPLEMENTATION FAILED — GOVERNANCE REVIEW REQUIRED**

ADR-0009 was accepted before source mutation. A candidate compiled, but the early V7
canonical-star check exceeded its predeclared mass-impact bound at the 2.0 M☉ target.
The owner request §21 requires **STOP** on an out-of-envelope delta. Validation stopped;
the candidate source was preserved outside the worktree and restored byte-identically
to the authenticated starting source. This record does **not** authorize the candidate,
artifact migration, corrected Phase-4D revalidation, or Phase 5.

## 1. Authentication and authority

- Worktree: `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-tov-surface`.
- Branch: `governance/tov-surface-contract`.
- Starting HEAD, upstream, and live remote: `8034a9c7d70a208a8c049a937675660785ac5f30`; clean tracked/untracked state.
- Audit ancestor: `b93be5f1df73eec0df1db1fc0ec39c524551e0b0`.
- Rotation/science ancestor: `e2fe0adf53d7975d3a5c57e84a8d3481ffd2ce41`.
- Master: `df859b5a73c4cac0c115f240744d89ce9f830b8d`; both named ancestry checks passed.
- ADR acceptance commit: `7a223d3390415e6b12737b5a9fbd9c87cfeadab6`, message `docs: accept TOV surface-event contract`.
- Exact owner Q1–Q14 decisions: `docs/adr/ADR-0009-tov-surface-event-and-termination.md:49`.
- Classification: scientific-semantic + numerical-method + structural, strictest governs.
  Normal migration, not GOVERNANCE §3.1; see `GOVERNANCE.md:44` and ADR-0009 Decision.
- One agent. No additional implementation branch, no merge, no artifact rewrite.

The starting Debug build passed **31/31** authenticated CTests (235.93 s) and the
self-contained build passed **17/17** (20.30 s), before ADR acceptance. This authentication
included the pre-existing suite; no corrected Phase-4D revalidation was launched.
Toolchain: AppleClang 17, CMake 4.2, GSL 2.7.1, macOS arm64. External EOS:
`/Users/keeper/Documents/CompactStar/data/compose/DS-CMF-1-with-crust/DS(CMF)-1_with_crust.eos`,
SHA256 `5747dd73256c0c28bc56be337cbb96d0918a54bc9ed9fc40984c5befd47ae5dd`.

Evidence location (local, not a scientific baseline):
`/Users/keeper/.codex/diagnostics/tov-surf-i-20260903/`.
`SHA256.json` authenticates the 26 preserved files, including old/new probe source,
old/new statically linked probe binaries, output TSVs, build/test logs, candidate source,
and `candidate.patch` (SHA256
`4f4f28cdc14d50bda499e26a891b455e10ae8e068834b4814479ca12002552dc`).
The patch is an **unvalidated stopped candidate**, never a production authority.
Original execution paths were `/tmp/compactstar-tov-surf-i/`; the preserved sources
retain those paths. Key raw evidence is `scratch/old.tsv:1`, `scratch/new.tsv:1`,
`scratch/old_probe.cpp:1`, and `scratch/new_probe.cpp:1` under the evidence location.

## 2. Candidate implementation and restoration

The preserved patch changed only `CompactStar/Core/TOVSolver.hpp`,
`CompactStar/Core/src/TOVSolver.cpp`, and `CompactStar/Core/src/NStar.cpp`:

- Added `TOVSolveStatus` and `LastSolveStatus()`, with the six accepted completion
  states plus `MASS_TARGET_NOT_REACHED` for an unsuccessful mass search.
- Removed the ordinary RHS's fatal **strict `<`** trial-pressure guard (the exact
  original differs from the request's schematic `<=`). Used `GetEDens(max(p,p_floor))`
  only for the energy lookup, retaining trial `p` in the force. Finite-state, metric,
  and EOS-error failures remained fatal. No unrelated EOS evaluator changed.
- Added private `PressureCoordinateODE` and `LocateSurface`: reuse the radial RHS
  with `dr/dp=1/f_p`, `dm/dp=f_m/f_p`; integrate from the last accepted state above
  the cutoff to exactly `p_cut`; require the result inside the radial bracket.
- Kept `PressureCutoff()` byte-identical, RK8PD and `1e-10` tolerances, inherited
  adaptive history, and the existing radial sampling ladder. No driver reset.
- Appended one final node at the event with EOS values/derivative, TOV `ν′`, and
  `ν=0` for reconstruction. Rejected nonfinite fields and invalid partitions.
- Cleared primitive output on every failure; positive return only at the event.
  `SolveToProfile` excluded failed coarse samples, removed nearest-sample fallback,
  rejected failed bisection members, and retained its existing `1e-4 M☉` tolerance.
- Guarded `Solve`, `GenTestSequence`, and `NStar::SolveTOV_Profile` publication.

These changes **compiled but did not complete validation**. They were restored after
V7 failed. Current tracked source still has the original guard and partial-return
behavior (`CompactStar/Core/src/TOVSolver.cpp:1496`, `:2572`, `:2653`); no completion
API or surface locator is installed. The restored source matches starting commit
`8034a9c` byte-for-byte. The authenticated build was rebuilt after restoration so
its binaries do not silently retain the stopped candidate.

Read-only caller inventory covered all compiled ordinary uses of the primitive:
`Solve`, `SolveToProfile`, `GenTestSequence`, the NStar wrapper, core/EOS/rotation/
thermal test builders and sequence-derivative helpers. `StarBuilder` reads existing
profiles rather than invoking this primitive. Test-builder adaptation and focused
failure-propagation tests were not completed before STOP. Thus there is no claim
that V12 covers every orchestrator.

Archived candidate source anchors (relative to the evidence directory):
`candidate-source/CompactStar/Core/TOVSolver.hpp:479`, `:1011`;
`candidate-source/CompactStar/Core/src/TOVSolver.cpp:1489`, `:2479`, `:2493`, `:2516`, `:2618`;
`candidate-source/CompactStar/Core/src/NStar.cpp:1297`. They describe only the stopped
candidate. The durable record carries the measurements needed to review the stop;
the local archive carries reproduction code and the exact uninstalled patch.

## 3. V7 stop evidence: canonical mass workflow

Both probes called the production `SolveToProfile` at the default resolution 10000,
then constructed `NStar` to obtain `I` and `B`. They used identical input masses,
EOS, compiler settings, and the unchanged mass-search tolerance. Values below are
raw diagnostic output, not fitted to the impact map. The signed relative difference
is `(new-old)/old`. Bounds: `|ΔM/M| <= 1e-9`; expected radius shifts approximately
+1.2/+1.0/+5.0/+3.7 m, each ±0.5 m (owner request §21; derivation §10/§12).

| Target M☉ | old εc (g/cm³) | candidate εc (g/cm³) | old M (M☉) | candidate M (M☉) | ΔM/M | ΔR (m) | V7 |
|---|---|---|---|---|---|---|---|
| 1 | 454550405078491.75 | 454550405078491.75 | 1.0000438094972419 | 1.0000438097285522 | 2.313000902e-10 | 1.199073307 | within map |
| 1.4 | 616488270506054.5 | 616488270506054.5 | 1.4000217830781583 | 1.4000217832812236 | 1.450444209e-10 | 1.027367781 | within map |
| 1.6 | 731253342677476.12 | 731253342677476.12 | 1.5999758341716293 | 1.5999758353006395 | 7.056419893e-10 | 5.035117994 | within map |
| 2 | 1298349261929558.8 | 1297985389788764.5 | 2.0000952962048859 | 1.9999985861663816 | -4.835271534e-05 | 4.154563390 | STOP |

At 2.0 M☉, the central density changes by approximately −0.0280%; both masses are
within the existing absolute `1e-4 M☉` target tolerance, but their separation is
`9.6710038504e-5 M☉` (`4.8352715337e-5` relative), outside V7 by about 48,353 times.
Its +4.15456 m radius change is inside +3.7±0.5 m; the **mass** bound triggers STOP.

This is a target-mass workflow discrepancy. It does **not** establish that the
pressure-coordinate solver changes mass outside the envelope at a *fixed* εc:
the two inputs actually used by the primitive differ. The precise change in coarse
bracketing/bisection was not isolated after STOP. No tolerance was tightened or
loosened, no central density was forced, and no fixed-density check was substituted
for this failed target-workflow comparison. Owner review must resolve the V7
comparison and mass-search contract before work resumes.

| Target M☉ | old R (km) | candidate R (km) | old I (geometric km³) | candidate I | ΔI/I | ΔB/B |
|---|---|---|---|---|---|---|
| 1 | 13.426323081955102 | 13.427522155261627 | 86.995758334054358 | 86.99575838204774 | 5.516749280e-10 | 2.465221360e-10 |
| 1.4 | 13.545323064955092 | 13.546350432736507 | 135.61613056691036 | 135.61613061268929 | 3.375626445e-10 | 1.594522292e-10 |
| 1.6 | 13.468323075955098 | 13.473358193949251 | 159.58714132523306 | 159.587141588083 | 1.647062042e-09 | 7.951110881e-10 |
| 2 | 12.712323183955165 | 12.716477747345323 | 193.72314421978143 | 193.72211615013879 | -5.306901490e-06 | -6.243555766e-05 |

The first three stars also satisfy the scratch tail estimates `|ΔI/I| <= 1e-8`,
`|ΔB/B| <= 1e-9`; the changed central density at 2.0 M☉ does not. Compactness,
`z_surf`, reconstructed surface `ν`, and photon factors were not independently
collected or checked before STOP. No downstream consistency claim is made.

## 4. Preliminary DS(CMF)-1 diagnostics, not full V1–V12 validation

The same probe completed 150 logarithmic densities `3e14*(1.6e15/3e14)^(i/149)`,
`i=0..149`, at resolution 10000; four canonical target-mass outputs; and 14
fixed-density resolution outputs. **168 nonempty outputs** all had final pressure
exactly `3.3518848999999998e25 dyn/cm²`, the EOS floor: the observed final-pressure
relative residual is **0**. Explicit status values and internal failed-solve counts
were not emitted by this smoke probe. Candidate control flow returns positive only
for `SURFACE_REACHED`; this is narrower evidence than a completed status-audited V1.
The `SolveToProfile` log records 72 rejected coarse-scan calls across the four searches;
their error categories were not counted. A claim of zero GSL failures over the
governed sweep would therefore be unsupported.

The five previously truncated default-resolution densities are exact points in that
150-density logarithmic scan. The frozen 1.6-M☉ audit density was also tested at 2500:

| εc (g/cm³) | resolution | old R (km) | candidate R (km) | old M (M☉) | candidate M (M☉) | ΔR (km) | ΔM/M |
|---|---|---|---|---|---|---|---|
| 1321824566394199.8 | 10000 | 12.348323235955199 | 12.684161211893699 | 2.0048458924293584 | 2.0061168267276099 | 0.335837976 | 6.339311680e-04 |
| 1351861473269815.2 | 10000 | 12.313323240955201 | 12.643765719058809 | 2.0119975360089364 | 2.0132269580548847 | 0.330442478 | 6.110455028e-04 |
| 1413998457266518 | 10000 | 12.243323250955209 | 12.56139017499819 | 2.0250080431301027 | 2.0260234087128306 | 0.318066924 | 5.014131110e-04 |
| 1429973936092014 | 10000 | 12.22232325395521 | 12.540482931206592 | 2.0277686568438091 | 2.0289339359820877 | 0.318159677 | 5.746607900e-04 |
| 1529689333192825.5 | 10000 | 12.11032326995522 | 12.412653398947755 | 2.0431255520766043 | 2.044075025098516 | 0.302330129 | 4.647159451e-04 |
| 731253300000000 | 2500 | 12.904228156539933 | 13.473358234380582 | 1.597681092159376 | 1.5999757719735503 | 0.569130078 | 1.436256475e-03 |

These outputs reach the floor and exhibit the expected restored outer layers. They
are promising diagnostics, not ratified fixes: the source was restored and the
full five-resolution/both-EOS sweep and missing-transition checks remain unperformed.

The 14 diagnostic resolutions were
`1250,2000,2450,2500,2510,2525,2550,3000,4000,5000,10000,20000,40000,80000`, at
`εc=7.312533e14 g/cm³`. Spread `(max-min)/value_at_10000` was
**1.218091310e-10 for M** and **3.378939567e-10 for R**. This subset is within
V4's bounds; it is not all 29 audit resolutions or the required uniform-grid/EOS sweep.

### Structural subset of the requested convergence dry run

| resolution | candidate R (km) | candidate M (M☉) |
|---|---|---|
| 1250 | 13.473358233159281 | 1.59997577184728 |
| 2500 | 13.473358234380582 | 1.5999757719735503 |
| 5000 | 13.473358234845024 | 1.5999757719822179 |
| 10000 | 13.473358234883204 | 1.5999757719894794 |
| 20000 | 13.473358234785138 | 1.5999757719905292 |
| 40000 | 13.47335823490849 | 1.5999757719906371 |
| 80000 | 13.473358235018576 | 1.599975771990775 |

The structural subset places 2500 with the other resolutions. The existing
convergence test's full configuration, fixed-mass arm, thermal quantities,
trajectories and observed orders were **not run** on the candidate after STOP.
The required complete grid-convergence gate is therefore **not passed**.

## 5. Validation disposition

| Line | Outcome |
|---|---|
| V1 | Incomplete: DS(CMF)-1 150-density smoke at 10000 only; no full status/GSL-failure accounting, no HW or three uniform grids. |
| V2/V3 | 168 diagnostic outputs nonempty, final `p/p_cut-1 = 0`; exhaustive status and node-property tests pending. |
| V4 | 14-resolution structural subset within bounds; full sweep and independent radial locator unperformed. |
| V5 | Five former 10000 failures and frozen 2500 case reach the floor in diagnostics; all-29 certification pending. |
| V6 | Not run: no claim concerning 9734/20361/38041 or the three perturbation steps. |
| V7 | **FAILED — STOP** at the 2.0-M☉ target-mass workflow; other three within map. |
| V8 | Not independently checked; no monotonicity, redshift, geometry, or photon-factor certification. |
| V9 | Preview only below; migration is blocked. |
| V10 | Starting 31/31 and 17/17 green; candidate compiled. Post-change first-order/contract/thermal/baryon suites not run after STOP. |
| V11 | Not run: the prerequisite corrected validation was not green. |
| V12 | Not run: no forced-failure/publication mutation performed. |

The independent bracketed radial cross-check was not implemented before STOP.
Driver-reset sensitivity was not measured; **DRIVER RESET NOT REQUIRED is not certified**.
No claim is made that the candidate has zero unexpected failures or that any final
baseline mismatch is an expected migration-only discrepancy. No post-change golden
comparison was run; the V7 failure is the recorded unexpected scientific gate failure.
The pre-change green suite remains evidence for the restored byte-identical source,
not for the candidate. No new tests or CMake changes were committed.

## 6. Migration preview and artifact preservation

All seven hashes below were authenticated before source work and rechecked unchanged
after restoration. Producers are the existing test executables (`tests/CMakeLists.txt:214`).
Classes are those already declared in `TOV_SURFACE_CONTRACT_DERIVATION.md:216`.
Dry-run deltas in this record are diagnostic star outputs, **not generated artifact
files or complete column-by-column artifact comparisons**.

| Artifact (`tests/baselines/`) | SHA256 (before = after) | Class | Producer | Expected movement / available diagnostic |
|---|---|---|---|---|
| `baryon_number_dscmf1_reference.tsv` | `8da5799d21da2017dd7dc49dfec8571ade6efba22846a652796118f248d4a646` | B | `baryon_number_cmf` | B upper limit; first three canonical ΔB/B ≤7.96e-10; 2.0 workflow outside map |
| `grid_convergence_cmf_1p6_debug.tsv` | `61d84ddcb87645197c5406c880b648fdf3bb9b0ed8c58350800ca2f2d296ff40` | A | `grid_convergence_cmf` | 2500 contamination; frozen-density R +0.569130 km, M +0.143626%; thermal/fixed-mass arms pending |
| `grid_convergence_cmf_1p6_trajectory.tsv` | `ca32863dabaa28fad63d5c36b287a3b94e9b6b85f11980bf2be4e65499d9a0c6` | A | `grid_convergence_cmf` | 2500 thermal trajectory and surface factors; no trajectory dry run |
| `hartle_I_dscmf1_debug.tsv` | `ddf018579364e9b3bf24f1e3a3e2577e70fb09b8b84eb718237d9da683aa9d15` | B | `hartle_moment_inertia_cmf` | R and I; first three ΔI/I ≤1.65e-9; 2.0 workflow outside map |
| `passive_cooling_cmf_1p6_debug.tsv` | `831744b0a206541fd0e24adc67876cc1ee4d02d89a580942a9fb0c6749999453` | B | `passive_cooling_cmf` | Tinf trajectory, R, lapse and photon area; no thermal dry run |
| `tov_dscmf1_reference.tsv` | `ba9f6ee51e501e5e5a2133f72d3d16f351e5c721eb3f7a7c04e4d922fbc13e28` | B | `tov_reference_cmf` | canonical R/M; §3 shows the stopped comparison |
| `tov_path_equivalence_dscmf1.tsv` | `bbf61e5fddb4709500f22a1eb11b1e20554f7463376619e86e96ea0a2540d871` | D | `tov_path_equivalence_cmf` | both ordinary paths and surface-dependent fields; exact outcome unmeasured |

No class C artifact is expected by the accepted impact map. **Artifacts regenerated:
NO.** The two historical audit/derivation records and historical GRID_CONVERGENCE
results are unchanged. No regression tolerance or baseline was changed.

## 7. Final governed state and next gate

- ADR-0009: **ACCEPTED**, Q1–Q14 recorded; **source conformance and validation not achieved**.
- INV-06: accepted event locator, unchanged cutoff/EOS-floor semantics;
  current source remains nonconformant. Do not label it numerically validated or `p=0`.
- INV-13: sampling-only target partition accepted; full numerical validation pending.
- ADR-0008 physics unchanged; corrected Phase-4D revalidation and first monopole
  baseline remain blocked. No Phase 5 work.
- MixedStar/dark core: **SEPARATE UNVALIDATED TRACK**, scientifically untouched.
  Its existing fatal guards remain (`CompactStar/Core/src/TOVSolver.cpp:1359`, `:1429`).
- All production, test, CMake and baseline files restored/unchanged relative to
  starting HEAD. Only ADR acceptance, this stopped record and status documents change.
- No implementation-success commit is made. The second commit records the stop
  instead of using the proposed `fix:` title. Non-force push only; no merge.

**Exactly one recommended next task:** owner-reviewed diagnosis and adjudication of
the 2.0-M☉ V7 target-mass discrepancy and its comparison/termination contract before
resuming TOV-SURF-I. Artifact migration remains downstream of successful V1–V12;
this stopped attempt does not authorize it.

Changed documentation files are this record, `docs/adr/ADR-0009-tov-surface-event-and-termination.md`,
`docs/adr/README.md`, `docs/SCIENTIFIC_INVARIANTS.md`,
`docs/architecture/CURRENT_ARCHITECTURE.md`, and `docs/MODERNIZATION_ROADMAP.md`.
The acceptance commit is listed in §1; the stopped-record commit is identified by
Git history and the final delivery report. The final report also records verified
local/remote equality after the non-force push.

## 8. TOV-SURF-I-R owner clarification and pre-restoration diagnosis

**Resume starting HEAD:** `96c1425f8503c2845e72f5755b02840798f6b958`.
Clean worktree, upstream and live remote equal; accepted ancestor `7a223d3` and
proposal ancestor `8034a9c` authenticated. Production/test bytes equal the pre-candidate
source. All seven artifact hashes in §6 remain unchanged. Candidate patch exists
and matches SHA256 `4f4f28cdc14d50bda499e26a891b455e10ae8e068834b4814479ca12002552dc`.
The complete archived manifest and starting 31/31 + 17/17 logs were authenticated;
these recent suites are reused because source and build state were restored and
no intervening source change exists. The original stopped record above remains
historical evidence. Classification remains scientific-semantic + numerical-method
+ structural under the normal migration path, not GOVERNANCE §3.1.

The owner clarified **V7a** as the identical-εc impact experiment, retaining the
`1e-9` relative mass bound and the existing radius envelope. **V7b** is the complete,
stable-branch, deterministic target-mass workflow with its unchanged absolute
`1e-4 M☉` tolerance; historical εc equality is not required. ADR-0009's accepted
scientific Decision is not rewritten. See its post-acceptance clarification.

### 8.1 Exact archived-binary bisection diagnosis

Before restoring the patch, LLDB traced the **preserved original executables**
`scratch/old_probe` and `scratch/new_probe` from the authenticated TOV-SURF-I archive.
Breakpoints read local variables only; no numerical variable or source was modified.
Commands/logs and their SHA256 manifest are preserved at
`/Users/keeper/.codex/diagnostics/tov-surf-ir-20260903/` (`old.lldb:1`, `new.lldb:1`,
`old-profile.lldb:1`). Both binaries reproduce their archived final masses exactly.

Both selected the same positive-slope coarse bracket:
`[806030255434967.75, 1551240399781547.5] g/cm³`.
Old endpoint masses `[1.7001029821130762, 2.0467345600684186] M☉`;
candidate `[1.700102983637509, 2.0467345605021521] M☉`.
The candidate excludes 18 incomplete low-density coarse samples, changing compacted
indices from `(20,21)` to `(2,3)` but **not** these endpoints or stable-branch selection.

Zero-based iteration trace; accepted means `|M-2| < 1e-4 M☉`:

| Behavior | iteration | εc (g/cm³) | M (M☉) | absolute residual (M☉) | accepted? |
|---|---|---|---|---|---|
| old | 0 | 1178635327608257.5 | 1.9617227924250455 | 0.038277207574954453 | no |
| old | 1 | 1364937863694902.5 | 2.0161260692509688 | 0.016126069250968822 | no |
| old | 2 | 1271786595651580 | 1.9927512688280649 | 0.007248731171935141 | no |
| old | 3 | 1318362229673241.3 | 2.0052551548445821 | 0.0052551548445820551 | no |
| old | 4 | 1295074412662410.5 | 1.9992211126858292 | 0.00077888731417075086 | no |
| old | 5 | 1306718321167826 | 2.0022907812973161 | 0.0022907812973160802 | no |
| old | 6 | 1300896366915118.3 | 2.0007693277192287 | 0.00076932771922866294 | no |
| old | 7 | 1297985389788764.5 | 1.9987052318295861 | 0.0012947681704138514 | no |
| old | 8 | 1299440878351941.5 | 2.0003847970065332 | 0.00038479700653315518 | no |
| old | 9 | 1298713134070353 | 2.0001919017217573 | 0.0001919017217573149 | no |
| old | 10 | 1298349261929558.8 | 2.0000952962048859 | 0.000095296204885908508 | yes |
| new | 0 | 1178635327608257.5 | 1.9617227927722212 | 0.038277207227778831 | no |
| new | 1 | 1364937863694902.5 | 2.016126070354213 | 0.016126070354212985 | no |
| new | 2 | 1271786595651580 | 1.9927512699268111 | 0.0072487300731889359 | no |
| new | 3 | 1318362229673241.3 | 2.0052551558256626 | 0.0052551558256626052 | no |
| new | 4 | 1295074412662410.5 | 1.999221112887285 | 0.00077888711271500988 | no |
| new | 5 | 1306718321167826 | 2.0022907828021772 | 0.0022907828021772048 | no |
| new | 6 | 1300896366915118.3 | 2.0007693277528409 | 0.00076932775284088706 | no |
| new | 7 | 1297985389788764.5 | 1.9999985861663816 | 0.0000014138336184021938 | yes |

**Cause: E — restoration of a previously truncated bisection member**, followed by
the unchanged first-acceptable stopping rule. It is not simply neighboring iterates
crossing a nearly equal threshold (A), a changed coarse bracket (B), removal of a
coarse sample affecting that bracket (C), or a stable-branch selection change (D).

At iteration 7 both evaluate `εc=1297985389788764.5 g/cm³`. LLDB inspection of the old
`tmp` shows **2479 nodes**, last `R=12.376323231955196 km`,
`M=1.9987052318295861 M☉`, `p=4.4009432910546427e31 dyn/cm²`,
`ε=70785195087364.766 g/cm³`: about `1.31e6` times the governed pressure cutoff,
a crust-truncated profile. The old root finder incorrectly treats its mass as a
sample and continues to iteration 10. The corrected member completes and returns
`1.9999985861663816 M☉` at iteration 7, satisfying the unchanged tolerance.

The candidate residual is `1.4138336184021938e-6 M☉`, versus the old accepted
`9.5296204885908508e-5 M☉`: the candidate is about **67.4 times closer** to the target.
This diagnosis requires no mass-tolerance or bisection-policy change. V7a and V7b
must still pass before the remaining V1–V12 campaign resumes.

## 9. TOV-SURF-I-R implementation and numerical validation

**TOV SURFACE EVENT IMPLEMENTED AND VALIDATED — ARTIFACT MIGRATION REQUIRED**

This section supersedes the implementation disposition in historical §§1–7 under
the owner's §8 clarification; it does not erase that stop. Clarification commit:
`03b7d58` (`docs: clarify TOV surface V7 validation semantics`). The accepted ADR-0009
Decision remains byte-identical. The exact archived candidate was applied only after
`git apply --check`; all three production files matched their archived SHA256 values.
No scientific implementation was reconstructed or opportunistically modified.

### 9.1 Production conformance and caller scope

The implementation is the candidate described in historical §2, now validated:

| Contract | Source authority / implementation |
|---|---|
| Completion enum and read API | `CompactStar/Core/TOVSolver.hpp:479`, `:1011`: six accepted states plus `MASS_TARGET_NOT_REACHED`; `LastSolveStatus()` |
| Trial-only continuation | `CompactStar/Core/src/TOVSolver.cpp:1490`: ordinary RHS evaluates only energy density at `max(p,p_floor)`; trial p retained in force; nonfinite state, invalid metric, missing/invalid EOS remain fatal |
| Canonical terminal locator | `CompactStar/Core/src/TOVSolver.cpp:2479`, `:2493`: `dr/dp=1/f_p`, `dm/dp=f_m/f_p` by calling the existing radial RHS; start at last accepted state above cut, land exactly on cut inside the bracket |
| Output and partition | `CompactStar/Core/src/TOVSolver.cpp:2516`: same RK8PD, abs/rel `1e-10`, initial radial step and sampling ladder; inherited adaptive history; invalid partitions fail closed |
| Final node | `CompactStar/Core/src/TOVSolver.cpp:2572`: exactly one event point, EOS fields and derivative at cutoff, event mass and TOV ν′, initial ν=0; strictly increasing r, no published sub-cutoff point |
| Failed primitive | `CompactStar/Core/src/TOVSolver.cpp:2522`: clears output and returns 0 with explicit status; positive result only `SURFACE_REACHED` |
| Target-mass search | `CompactStar/Core/src/TOVSolver.cpp:2613`: only complete coarse/bisection samples; no nearest fallback; failed bisection aborts. `N=24`, positive-slope criterion, `mass_tol=1e-4`, first-acceptable stopping unchanged (`:2654`, `:2707`, `:2742`, `:2771`) |
| Other ordinary callers | `CompactStar/Core/src/TOVSolver.cpp:1820`, `:3232`; `CompactStar/Core/src/NStar.cpp:1297`: refuse incomplete primitive results before append/finalization/analysis/export |

`PressureCutoff()` is byte-identical to the restored source, including its value
`max(1e-15 p_c,eos_tab.pre[0])`. No other EOS continuation, RotationSolver equation,
Geometry, thermal, BNV or rotochemical physics changed. Test builders already
checking positive return inherit the primitive's cleared-output guarantee. The
previously unchecked reference/label/sequence builders now explicitly guard
completion. `StarBuilder` imports profiles; it introduces no second ordinary TOV
integration path. The pre-existing Path-1 mirror surface-scalar asymmetry remains
separate postprocessing debt, characterized by `tov_path_equivalence_cmf`; the TOV
nodes and reconstructed surface boundary agree between paths.

MixedStar/dark-core equations and guards are unchanged: **SEPARATE UNVALIDATED TRACK**.
No new public uniform-grid or second production locator was introduced. Uniform-grid
experiments and the radial oracle are test-side only.

### 9.2 V7a — identical central conditions

`tests/core/tov_surface_contract.cpp:101` calls the primitive at the four inherited
central densities, comparing against authenticated old diagnostic rows; it does not
call independent target-mass workflows. Every bound passes, unchanged.

| Label M☉ | εc (g/cm³) | old M | corrected M | ΔM/M | old R (km) | corrected R (km) | ΔR (m) |
|---|---|---|---|---|---|---|---|
| 1 | 454550405078491.75 | 1.0000438094972419 | 1.0000438097299555 | 2.3270341209524759e-10 | 13.426323081955102 | 13.427522155210758 | 1.1990732556554917 |
| 1.4 | 616488270506054.5 | 1.4000217830781583 | 1.4000217832812236 | 1.450444209183388e-10 | 13.545323064955092 | 13.546350432736507 | 1.0273677814147675 |
| 1.6 | 731253342677476.12 | 1.5999758341716293 | 1.5999758353006395 | 7.0564198928479982e-10 | 13.468323075955098 | 13.473358193949251 | 5.0351179941525714 |
| 2 | 1298349261929558.8 | 2.0000952962048859 | 2.0000952969645365 | 3.7980707467966113e-10 | 12.712323183955165 | 12.715982774062567 | 3.659590107401911 |

At the inherited 2.0-M☉ density the surface-only relative mass change is
`3.7980707468e-10`, with radius +3.6595901074 m. This meets `1e-9` and +3.7±0.5 m.

### 9.3 V7b — complete stable-branch target solves

`tests/core/tov_surface_contract.cpp:113`: each target is called twice, and every
TOVPoint field and species-label vector compares bit-identically. All statuses are
`SURFACE_REACHED`, cutoff exact, no fallback; **18 failed coarse samples excluded,
zero failed bisection samples consumed**, for each call. Successful production
bisection cannot consume a failed member because it immediately returns failure.

| Target M☉ | Returned M☉ | absolute residual M☉ | selected εc (g/cm³) | complete bracket εc low / high |
|---|---|---|---|---|
| 1 | 1.0000438097285522 | 4.3809728552224314e-05 | 454550405078491.75 | 418816305176200.12 / 806030255434967.75 |
| 1.4 | 1.4000217832812236 | 2.1783281223708428e-05 | 616488270506054.5 | 418816305176200.12 / 806030255434967.75 |
| 1.6 | 1.5999758353006395 | 2.4164699360618158e-05 | 731253342677476.12 | 418816305176200.12 / 806030255434967.75 |
| 2 | 1.9999985861663816 | 1.4138336184021938e-06 | 1297985389788764.5 | 806030255434967.75 / 1551240399781547.5 |

The 1.0/1.4/1.6 bracket masses are `[0.89254069392360491,1.700102983637509] M☉`;
the 2.0 bracket masses are `[1.700102983637509,2.0467345605021521] M☉`. Each has
positive slope and contains its target. Coarse exclusion changes neither the
historical 2.0 bracket endpoints nor the stable-branch criterion (§8.1).

Test-only tighter 2.0 root (`tests/core/tov_surface_contract.cpp:437`):
`εc=1297990706760487.5 g/cm³`, `M=2.0000000000831948 M☉`, absolute residual
`8.319478439489103e-11 M☉`; final bracket
`[1297990706066456.5,1297990706760487.5] g/cm³`.
The production εc lies below this bracket by about `5.317e9 g/cm³` and is accepted
with residual `1.4138336184021938e-6 M☉`, well inside `1e-4`. Production tolerance,
coarse N, branch rule and stopping semantics were not changed.

**SOLVETOPROFILE TARGET-ROOT PRECISION / DETERMINISM DEBT:** future owner-governed
work may choose a tighter target contract; this is not needed for surface conformance.
Repeated identical calls are deterministic now. Cross-version selected εc equality
is not promised by the existing first-acceptable stopping rule.

### 9.4 V1–V4 — both EOS families and independent locator

`tests/core/tov_surface_contract.cpp:167`, `:210`. Each EOS receives 150 logarithmic
central densities over `[3e14,1.6e15]` plus the five exact inherited failure densities
(155 entries; repeated densities intentionally retained), at production resolutions
2500/5000/10000/20000/40000. **775 production solves per EOS, 1550 total, zero GSL
failures, every status SURFACE_REACHED.** Every published node is validated, not only
its last pressure: increasing radius, no sub-cutoff state, exactly one cutoff node,
finite derivative and fields, surface EOS values and TOV ν′ consistent.

For each production solve, an independent diagnostic integrates the same radial RHS
and bisects by fresh radial re-integration from the last above-cut state until the
bracket is ≤`1e-6 cm`. It never calls the pressure-coordinate locator. Three uniform
grids (1.75/7/28 m) are also integrated this way. **1240 independent solves per EOS,
2480 total**. The eight-partition spread is `(max-min)/abs(value_at_production_10000)`.
The independent uniform results are diagnostic, not additional production outputs.

| Metric | CMF | HW | Accepted bound |
|---|---|---|---|
| production final `abs(p/p_cut-1)` | 0 | 0 | 1e-8 |
| eight-partition M spread | 1.7684480612e-10 | 4.4847630666e-10 | 1e-9 |
| eight-partition R spread | 2.1580025491e-9 | 6.5252812003e-9 | 1e-8 |
| compactness M/R spread | 2.1554609700e-9 | 6.5099498277e-9 | 1e-8 |
| surface lapse spread | 3.0356459584e-10 | 5.1869391140e-10 | 1e-8 |
| maximum independent-locator ΔM/M | 5.4944937489e-13 | 3.9601655288e-13 | 1e-9 |
| maximum independent-locator ΔR/R | 6.4100824648e-10 | 6.4765219854e-9 | 1e-8 |

HW uses the existing governed HT68/Harrison-Wheeler densified fixture
(`tests/rotation/hartle_thorne_1968_hw_eos.hpp:136`), unchanged. CMF's cutoff is
floor dominated; HW's is `1e-15 p_c` dominated. The HW radius margin is smaller than
the CMF margin and smaller than the earlier single-star scratch margin; the measured
full-sweep bound passes without changing production tolerances or the locator.

### 9.5 V5/V6 — former failures and sequence derivative

All **29** audit resolutions reach the event at `εc=7.312533e14`, on both EOS families
(`tests/core/tov_surface_contract.cpp:264`). The four formerly truncating resolutions
and five exact formerly failing 10000 densities also pass the dedicated detector
control (9/9). The old/new CMF comparisons in historical §4 remain the measured
surface correction; no central density was adjusted for those comparisons.

At fixed audit εc, resolution 2500 changes from
`R=12.904228156539933 km, M=1.597681092159376 M☉` to
`R=13.473358234380582 km, M=1.5999757719735503 M☉`:
+0.569130077841 km and +0.1436256475% M. All five default failures recover the
floor endpoint, restoring +0.30233…0.33584 km (historical §4 table).

V6 tests ten representative densities (four inherited canonical densities, the
frozen audit density and five former failure densities) × perturbations
5e-4/1e-3/2e-3 × eight production resolutions (including 2000/2500/2550/4000)
and three uniform grids, for **both EOSs**. Derivatives use the actual central
pressure difference and existing geometric-unit conversion, not nominal ε increments
(`tests/core/tov_surface_contract.cpp:283`).

| Derivative check | CMF | HW | Bound |
|---|---|---|---|
| maximum fixed-step partition spread | 1.3971355937e-6 | 2.1696618402e-6 | 1e-5 |
| maximum 5e-4 vs 1e-3 step-pair relative difference | 2.4714621830e-6 | 4.2407260603e-7 | 1e-5 |

At the frozen CMF audit density and step 5e-4:

| resolution | corrected dM/dpc (km³) |
|---|---|
| 2000 | 9671.2418632523331 |
| 2500 | 9671.2399645384176 |
| 2550 | 9671.2394091534643 |
| 10000 | 9671.2398954312903 |

The former 9734/20361/38041 pathologies are absent. Step 2e-3 is reported in the raw
diagnostic table; the accepted step-pair gate remains 5e-4 versus 1e-3.

### 9.6 V8 and inherited adaptive history

The built NStar surface uses the event's r and m; reconstructed `ν(R)=0.5 ln(1-2M/R)`,
`ExpNuSurface`, metric lambda, and compactness agree (`tests/core/tov_surface_contract.cpp:415`).
Here `z_surf` means the **lapse** `exp(ν)`; physical gravitational redshift is
`1/z_surf-1`. No claim that this finite-cutoff surface is the `p=0` physical surface.

| CMF label | M/R | surface lapse | physical redshift | I (km³) | B |
|---|---|---|---|---|---|
| 1 | 0.10997484953106074 | 0.88320456347206366 | 0.13224052655342811 | 86.995758381598748 | 1.2738887314048103e+57 |
| 1.4 | 0.15260990251756115 | 0.83353475930214072 | 0.19971001669711863 | 135.61613061268929 | 1.8321833628708972e+57 |
| 1.6 | 0.17535081786373727 | 0.80579052133450002 | 0.24101732835459799 | 159.587141588083 | 2.124575697168995e+57 |
| 2 | 0.23225816253921916 | 0.73176750059124773 | 0.36655426647402067 | 193.72314439588695 | 2.7457630623997459e+57 |

At the same fixed εc, all four I/B shifts meet the surface-tail impact map
(I ≤1e-8, B ≤1e-9). The 2.0 target-mass workflow's larger I/B movement is separately
attributed to its permitted εc change (§8.1/§9.3), not confused with a surface-only tail.
HW's four representative built profiles also pass these boundary checks.

The extended thermal dry run checks the **actual** photon-derived radius,
`exp(2ν_surface)` and emitting area `4π R_cm² exp(2ν_surface)` against the same final
profile node, to 1e-13 relative (`tests/thermal/grid_convergence_cmf.cpp:256`).
No downstream formula was repaired. Radius/lapse values are invariant within the
governed bounds; their residual numerical scatter is not strictly monotonic at the
solver floor, and no meaningful convergence order is claimed for these surface scalars.

**DRIVER RESET NOT REQUIRED.** Test-only full driver reallocation per output segment
(discarding both RK state and inherited step size) versus inherited history, at four
representative densities and resolutions 2500/3000/10000:

| EOS | maximum ΔM/M | maximum ΔR/R |
|---|---|---|
| CMF | 4.1772807435e-11 | 2.5732260767e-10 |
| HW | 2.0258439370e-10 | 1.8686086012e-10 |

Both satisfy the accepted partition bounds. Production never resets per target.

### 9.7 V11/V12 mutation detectors

The exact script and logs are in the resume diagnostic archive (`detectors.py:1`).
Every temporary source edit was restored byte-identically, including on error via
`finally`; restored TOVSolver.cpp SHA256
`d5da96e512453ae9ce18d7e6b3334b82bf17dd1def3230dd24ad267bc1fc7ad8` equals the candidate.

- **V11:** normal control 9/9 complete, exit 0. Restoring `p <= p_cut -> GSL_EBADFUNC`
  at the unique ordinary RHS site makes all nine known failures recur, detector exit 1.
  Diagnostic-only last-accepted logging reproduces the four audit truncation radii
  12.904228/12.908594/12.894285/12.904391 km and all five default-resolution truncations
  byte-for-byte in the recorded r/m/p values. No trial state was published.
- **V12:** force GSL failure beyond 12.9 km. With governed publication, primitive,
  target search, Solve, GenTestSequence, NStar wrapper and derivative builder all
  reject; primitive status GSL_FAILURE, outputs cleared; focused test exit 0.
  Temporarily mislabeling the partial profile complete and bypassing clearing makes
  the primitive and target search return **2553 points**, Solve export 2 stars,
  GenTestSequence export 16 stars, and NStar finalize a partial star; detector exit 1.
  The test-side derivative builder still rejects through its independent terminal-point
  validation. Both mutation and forced-failure injection were removed.

Invalid-input contracts also pass: zero/nonfinite/nonincreasing partition,
nonfinite central density, missing EOS, radius exhaustion, unbracketed target,
cleared species labels and status reset success→failure→success
(`tests/core/tov_surface_contract.cpp:326`).

### 9.8 Full corrected convergence dry run

Existing 1.6-M☉ configuration, both fixed εc and fixed target-mass arms, resolutions
1250/2500/5000/10000/20000/40000/80000; static temperature 1e8 K; cooling 100 yr→1 Myr,
RKF45, rtol 1e-6, atol 1e-10. The dedicated `--surface-dry-run` mode cannot be combined
with `--emit-dir`; it writes no durable artifact. It retains the original thermal
model and reports actual sampled-integral convergence (`grid_convergence_cmf.cpp:477`).

Both arms select `εc=731253342677476.125 g/cm³` at every resolution and have the
same structure. **2500 is on the same complete physical branch: YES**.

| resolution | M (M☉) | R (km) | Cstar (erg/K) | Lnu (erg/s) | Lgamma (erg/s) | Tinf at 1 Myr (K) |
|---|---|---|---|---|---|---|
| 1250 | 1.59997583515575 | 13.4733581921643 | 2.4311040554e+38 | 4.4124361239e+39 | 8.7024922545e+33 | 467746.6154261609 |
| 2500 | 1.59997583528467 | 13.4733581925929 | 2.4279138914e+38 | 4.4462589593e+39 | 8.7052265908e+33 | 466364.4064455564 |
| 5000 | 1.59997583529338 | 13.4733581937279 | 2.4283250949e+38 | 4.4450569927e+39 | 8.6994880009e+33 | 466799.1415156969 |
| 10000 | 1.59997583530064 | 13.4733581939493 | 2.4280622085e+38 | 4.4526532394e+39 | 8.6943571297e+33 | 466916.5849204648 |
| 20000 | 1.59997583530169 | 13.4733581934893 | 2.4280434984e+38 | 4.4568634339e+39 | 8.6947047675e+33 | 466831.0110647354 |
| 40000 | 1.59997583530179 | 13.473358193845 | 2.4280645886e+38 | 4.4594614189e+39 | 8.6947425571e+33 | 466795.6522849235 |
| 80000 | 1.59997583530193 | 13.4733581938797 | 2.4280690353e+38 | 4.4617796011e+39 | 8.6944799205e+33 | 466777.5322990529 |

Raw `SURFACE_GRID` rows also contain B, lapse and dlnT/dt; `SURFACE_TRAJ` rows contain
all nine epochs and Lnu/Lgamma/Cstar along each trajectory. Measured finest-three
orders: B **2.235**, Cstar **2.246**, Lnu **0.164**, dlnT/dt **0.148**, terminal Tinf
**0.965**; trajectory-norm contraction order **0.146**. Lgamma's order is not reliable;
M/R/lapse are at the numerical floor. These are measured finite-resolution orders,
not a claim of universal asymptotic order. No final regression tolerance was chosen.

### 9.9 V9/V10 — artifact migration preview and legacy-test review

All seven original hashes in historical §6 were reverified unchanged. The preview
below supplements those hashes; no expected file was edited or regenerated.

| Artifact | Class | Producer | Expected / measured movement |
|---|---|---|---|
| grid_convergence_cmf_1p6_debug.tsv | A | grid_convergence_cmf | contaminated 2500 rows replaced by full star; R +0.569130 km; fixed-mass εc no longer compensates missing crust; complete extended dry run above |
| grid_convergence_cmf_1p6_trajectory.tsv | A | grid_convergence_cmf | 2500 trajectory changes; nine epochs × seven resolutions measured in scratch |
| tov_dscmf1_reference.tsv | B | tov_reference_cmf | fixed-εc ΔR +1.199/+1.027/+5.035/+3.660 m; mass within 1e-9; target workflows interpreted under V7b |
| passive_cooling_cmf_1p6_debug.tsv | B | passive_cooling_regression | surface factors and thermal trajectory; terminal Tinf −1.05949e-4 relative, max reported luminosity change 8.47219e-4 |
| hartle_I_dscmf1_debug.tsv | B | hartle_moment_inertia_cmf | R columns move; fixed-εc I tails <1e-8; permitted target-root movement separated |
| baryon_number_dscmf1_reference.tsv | B | baryon_number_cmf | fixed-εc B tails <1e-9; target 2.0 root changes B by −6.24356e-5 relative because εc changes |
| tov_path_equivalence_dscmf1.tsv | D | tov_path_equivalence_cmf | both physical paths agree; emission not performed, exact byte outcome deferred |

Four **expected legacy/reference failures** are retained, not hidden by adjusted
expectations (`ctest-impact.log` in the archive):

1. `passive_cooling_regression`: 29 comparisons against old surface-dependent values;
   changes 1e-5…8.47e-4, within the declared thermal/surface-impact scale.
2. `tov_reference_cmf`: B7's approximate omitted-layer attribution spread moves to
   0.209690; individual radius shifts are independently within V7a. The old heuristic
   threshold is retained for migration review.
3. `hartle_moment_inertia_cmf`: B3a compares sampled I spread (`6.573e-5`) with R spread
   (`9.477e-11`). R now converges as an event while I remains a sampled integral;
   this old radius-quantization assumption is obsolete. First-order independent
   physics and normalization tests pass; no I equation was changed.
4. `baryon_number_cmf`: old exact structure/B expectations move by the surface tail;
   the 2.0 target-root component is explicitly allowed by V7b, while its fixed-density
   surface-only result meets the original tail bound.

The initial `tov_path_equivalence_cmf` below-floor expectation also failed because it
expected an incomplete clamped star to be published. This API-contract assertion was
updated narrowly to require both paths to reject; no numerical reference or baseline
was relaxed. The final path test passes.

Final behavior suite: **37/37 passed** (468.21 s), including all 13 new surface tests
and 24 tests carrying the scientific label. Four legacy/reference failures above
are classified separately; this is not a claim that the full golden suite is green.
Final self-contained suite: **20/20 passed** (65.23 s), including five HW surface tests.
The three independent Phase-4D physics/published tests (two in the self-contained
configuration) were excluded because corrected Phase-4D revalidation is not this task.
A concurrent launch initially caused `heat_capacity_v1` in the self-contained build
to lose `compactstar_hcv1_fixture/linear2x/eos.t`: both copies delete the same fixture
root (`tests/thermal/heat_capacity_v1.cpp:273`, `:532`). That runner collision is retained
in `ctest-self-concurrent.log`; the successful serial rerun is in
`ctest-self-serial.log`. No heat-capacity source or expectation changed.
There are no unresolved unexpected failures or unexplained scientific discrepancies.

### 9.10 Final status and dependency boundary

ADR-0009: **ACCEPTED / SOURCE CONFORMED / NUMERICALLY VALIDATED** for ordinary NStar.
INV-06: **GOVERNED / CONFORMED / VALIDATED**, at the finite `p=p_cut` event, never a
claim of a zero-pressure physical surface. INV-13: radial output partition is
sampling only within the governed bounds; sampled integrals retain resolution error.
The accepted Decision, PressureCutoff, mass_tol, stable-branch criterion and coarse
N are unchanged. The prior acceptance, stop and clarification commits are preserved.

ADR-0008's physics is not reopened. Corrected Phase-4D revalidation and the first
monopole baseline require the **separate artifact migration first**. No Phase 5,
MixedStar correction, merge, or artifact regeneration occurred. Documentation of
source conformance does not constitute new scientific ratification of unrelated code.

**Exactly one recommended next task:** migrate the seven governed artifacts once,
recording every old→new delta and new hash, before corrected Phase-4D revalidation.

### 9.11 Reproduction, retained evidence and changed files

The resume archive is `/Users/keeper/.codex/diagnostics/tov-surf-ir-20260903/`.
Its original root `SHA256.json` authenticates the pre-restoration bisection traces.
The added `validation/SHA256.json` authenticates all retained scratch reports, logs,
mutation scripts, final build reports, source snapshots and final validation binary.
Original scratch execution paths remain `/tmp/compactstar-tov-surf-ir/`; these are
validation evidence, never governed golden artifacts. In particular,
`validation/scratch/detectors.py:1` reproduces V11/V12 and restores source in `finally`.
`validation/final-integrity.json:1` records unchanged artifact hashes and exact
candidate-source and accepted-Decision checks.

Both existing Debug builds were rebuilt. Equivalent commands from this worktree
(run the suites serially to avoid their pre-existing shared fixture root):

```sh
cmake --build build -j6
cmake --build build-selfcontained -j6
ctest --test-dir build --output-on-failure -E 'hartle_monopole_physics|hartle_monopole_published|passive_cooling_regression|tov_reference_cmf|hartle_moment_inertia_cmf|baryon_number_cmf'
ctest --test-dir build-selfcontained --output-on-failure -E 'hartle_monopole_physics|hartle_monopole_published'
```

The impact run retained in `validation/scratch/ctest-impact.log:1` excludes only the
independent Phase-4D tests and includes the legacy comparisons. The final path-test
contract correction and new surface tests then pass in
`validation/scratch/ctest-behavior-final.log:1`. The final scientific source is the
authenticated candidate, not either temporary detector mutation.

Files changed by the implementation commit (16):

- Production: `CompactStar/Core/TOVSolver.hpp`, `CompactStar/Core/src/TOVSolver.cpp`,
  `CompactStar/Core/src/NStar.cpp`.
- Validation: `tests/CMakeLists.txt`, `tests/core/tov_surface_contract.cpp`,
  `tests/core/tov_gen_test_sequence_cmf.cpp`, `tests/core/tov_path_equivalence_cmf.cpp`,
  `tests/core/tov_reference_cmf.cpp`, `tests/rotation/hartle_moment_inertia_cmf.cpp`,
  `tests/thermal/grid_convergence_cmf.cpp`.
- Documentation: this record, `docs/adr/ADR-0009-tov-surface-event-and-termination.md`,
  `docs/adr/README.md`, `docs/SCIENTIFIC_INVARIANTS.md`,
  `docs/architecture/CURRENT_ARCHITECTURE.md`, `docs/MODERNIZATION_ROADMAP.md`.

The implementation commit and post-push local/remote equality are recorded by Git
history and the delivery report, preserving `7a223d3`, `96c1425` and `03b7d58`.
