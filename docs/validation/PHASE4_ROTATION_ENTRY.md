# Phase 4 rotation entry — audit, first-order normalization derivation, second-order equation map

> **FORMAL STATUS: `PHASE-4 ROTATION ENTRY AUDIT COMPLETE — ADR-0006 OWNER ADJUDICATION REQUIRED`**
>
> **SUPERSEDED IN PART — 2026-09-02.** ADR-0006 was **ACCEPTED** and implemented in Phase 4A.
> This record remains the durable audit of the state **at Phase-4 entry** (`df859b5`) and its
> derivations, provenance map and second-order equation map stand unchanged. But its
> descriptions of *current* first-order behavior — the API map (§4), the call/reachability table
> for the first-order entry points (§5), the unit-annotation verdicts (§9) and the source line
> numbers throughout — describe the pre-4A source. For current first-order behavior see
> `docs/validation/PHASE4A_FIRST_ORDER_NORMALIZATION.md` and INV-07. **The O(Ω²) sections
> (§10–§13) are NOT superseded:** the candidate is byte-identical and every defect recorded there
> is still present.
>
> This is an **audit / scientific-derivation / governance** record (increment 4A-0). **No
> production source, test, baseline or CMake file was changed.** No candidate was activated,
> repaired, normalized or baselined. The one scratch measurement program used in §12–§13 lives
> outside the repository tree and is not tracked. `ADR-0006` is drafted **PROPOSED**; nothing here
> accepts it, implements it, or certifies any O(Ω²) quantity.

| Field | Value |
|---|---|
| **Starting `master`** | `df859b5a73c4cac0c115f240744d89ce9f830b8d` — `master` = `origin/master` = `HEAD` of the canonical checkout, working tree clean |
| **Branch / worktree** | `physics/rotation-correctness` at `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-rotation`, created from that exact SHA (no force operations; the Phase-3 consolidation worktree was not reused) |
| **Change class** | documentation / governance (`GOVERNANCE.md` §2) — every factual claim below carries a `path:line` citation at `df859b5` |
| **Authority read, in order** | `GOVERNANCE.md`; ADR-0001…ADR-0005 (all ACCEPTED); `docs/SCIENTIFIC_INVARIANTS.md` (INV-05, INV-07, INV-08, INV-09); `docs/MODERNIZATION_ROADMAP.md` Phase 4; `docs/architecture/CURRENT_ARCHITECTURE.md`; `docs/architecture/TARGET_ARCHITECTURE.md`; `docs/validation/PHASE3_CLOSEOUT.md`; `docs/validation/HARTLE_MOMENT_INERTIA.md`; `docs/validation/GRID_CONVERGENCE.md`; `docs/validation/PHASE2B_CLOSURE.md`; `AGENTS.md` |
| **Literature authority** | J. B. Hartle (1967), ApJ **150**, 1005 — first-order frame-dragging equation (already validated, 2B-4B) and the l = 0 second-order equations for `(m₀, p₀*)`; Hartle & Thorne (1968), ApJ **153**, 807; the same monopole system as reproduced in Paschalidis & Stergioulas, *Living Rev. Relativ.* **20**, 7 (2017). The transcription used here is written out in §10 so that it is itself auditable |
| **Build / toolchain** | `Debug` (default), AppleClang 17, CMake 4.2, GSL 2.7.1 at `/opt/local`, macOS arm64 — as `docs/build/MACOS_BUILD.md` |

---

## 1. Authentication

```
cd /Users/keeper/Documents/CompactStar/repo/CompactStar && git fetch origin
master        = df859b5a73c4cac0c115f240744d89ce9f830b8d
origin/master = df859b5a73c4cac0c115f240744d89ce9f830b8d
HEAD          = df859b5a73c4cac0c115f240744d89ce9f830b8d     working tree: clean
```

Existing worktrees: the canonical checkout (`master`) and
`worktrees/CompactStar-consolidation` (`modernization/behavior-preserving-consolidation`, also at
`df859b5`). Neither carried `physics/rotation-correctness`; it did not exist locally or on
`origin`. It was created with `git worktree add -b physics/rotation-correctness … df859b5…`.

## 2. Inherited Phase-3 baseline — re-authenticated in the new worktree at `df859b5`

| Configuration | Registered | Result | Wall time |
|---|---|---|---|
| Full — `-DCOMPACTSTAR_EOS_DATA_ROOT=/Users/keeper/Documents/CompactStar/data/compose` | 19 | **19/19 PASS** | 203.14 s |
| Self-contained — no data root (`build-selfcontained/`) | 10 | **10/10 PASS** | 14.82 s |

`git status --porcelain` empty before and after both builds (both build directories are
gitignored: `.gitignore:6-7`). Nothing was regenerated or re-baselined.

Seven durable artifacts, `shasum -a 256` at `df859b5`, identical to `PHASE3_CLOSEOUT.md` §1:

| Artifact | SHA-256 |
|---|---|
| `passive_cooling_cmf_1p6_debug.tsv` | `831744b0a206541fd0e24adc67876cc1ee4d02d89a580942a9fb0c6749999453` |
| `tov_dscmf1_reference.tsv` | `ba9f6ee51e501e5e5a2133f72d3d16f351e5c721eb3f7a7c04e4d922fbc13e28` |
| `grid_convergence_cmf_1p6_debug.tsv` | `61d84ddcb87645197c5406c880b648fdf3bb9b0ed8c58350800ca2f2d296ff40` |
| `grid_convergence_cmf_1p6_trajectory.tsv` | `ca32863dabaa28fad63d5c36b287a3b94e9b6b85f11980bf2be4e65499d9a0c6` |
| **`hartle_I_dscmf1_debug.tsv`** | **`ddf018579364e9b3bf24f1e3a3e2577e70fb09b8b84eb718237d9da683aa9d15`** |
| `baryon_number_dscmf1_reference.tsv` | `8da5799d21da2017dd7dc49dfec8571ade6efba22846a652796118f248d4a646` |
| `tov_path_equivalence_dscmf1.tsv` | `bbf61e5fddb4709500f22a1eb11b1e20554f7463376619e86e96ea0a2540d871` |

## 3. Governing-document read — statements found stale at `df859b5`

The invariant register's rotation entries carry `d91c31b`-era line numbers and one reachability
statement that the current header contradicts. Every item below was checked against the actual
source, not assumed.

| Statement | Where | Finding |
|---|---|---|
| "Zero call sites; `rot_solver` is private with no accessor — **structurally unreachable**" | `SCIENTIFIC_INVARIANTS.md` INV-08 | **STALE / INCORRECT.** `RotationSolver` is a public class whose `AttachNStar`, `FindNMomInertia`, `SolveHartle2_N` and `GetHartleResult` are public and defined. An external instance attached to any `NStar` executes the candidate. Demonstrated §13. *Zero repository callers* is true; *unreachable* is not |
| `RotationSolver — O(Ω²)`: "UNREACHABLE SCAFFOLDING · CANDIDATE — `rot_solver` private, no accessor" | `CURRENT_ARCHITECTURE.md:95` | same — corrected in this task |
| "It is unreachable and unverified" | `CURRENT_ARCHITECTURE.md:487` | half stale — corrected |
| "SCAFFOLDED · CANDIDATE — unreachable" | `TARGET_ARCHITECTURE.md:46` | stale (non-normative) — minimal correction |
| "UNREACHABLE SCIENTIFIC CANDIDATE — PROTECTED" | `PHASE3_CLOSEOUT.md:171`, `:268` | historical closeout record — **left as written**, superseded by this record |
| "The rotation half of this entry must be re-audited against `9f70f14`" | INV-05 | discharged here, §7.3 |
| `RotationSolver.cpp:1240`, `:1276-1288` | INV-07 evidence | `d91c31b` numbering; current: seed `:701`, extraction `:737-739` |
| `:1554-1700`, `:1591-1678`, `:1502`, `:1504-1507`, `:1514`, `:1696`, `:1444`, `:1518`, `:1542-1548`, `:1702-1705` | INV-08 evidence | `d91c31b` numbering; current: `SolveHartle2_N` `:1126-1270`, superposition `:1233-1250`, `j²` comment `:1076-1082`, `S_m` `:1088`, `delta_M` `:1268`, `[FIX…]` `:1019`, `???` `:1091`, MixedStar stubs `:1114-1120` and `:1274-1277` |
| "(incompressible is dε/dp → ∞)" | INV-08 | **wrong**: incompressible matter has `dε/dp → 0`; the candidate's `1.0` fallback is the causal-limit `c_s² = 1`, not an incompressible limit (§11C) |
| "The superposition construction is correctly implemented" | INV-08 verified sub-claim | the **arithmetic** is correct; the **condition it imposes** (`δp(R) = 0`) is not Hartle's fixed-central-density condition and forces `ξ₀(R) = 0` (§11E) |
| "Make `HartleResult` reachable from `NStar`" | roadmap Phase 4 item | premise partly stale: a raw `HartleResult` is already reachable through an external solver; the real gap is that no **normalized** product is exposed (§16, ADR-0006 Q4) |

## 4. Current rotation API map — `CompactStar/Core/RotationSolver.hpp` at `df859b5`

| Member | Declared | Defined in `RotationSolver.cpp` | Symbol in `libCompactStar.a` | Access |
|---|---|---|---|---|
| `RotationSolver()`, `~RotationSolver()` | `:293-294` | `:45-52` | yes | public |
| `AttachNStar(NStar*)` / `AttachMixedStar(MixedStar*)` | `:306`, `:308` | `:55-66`, `:69-80` | yes | public |
| `GetNStar()`, `GetMixedStar()` | `:310`, `:312` | `:652-661` | yes | public |
| `GetMass(const double&)`, `GetPress(double)` | `:315`, `:318` | **none** | **absent** | public — **DECLARED, UNDEFINED** |
| `GetEDens`, `GetNu` | `:324`, `:334` | `:216-231` | yes | public |
| `GetInitOmegaBar()` | `:327` | `:354-357` | yes | public |
| `GetOmegaSeq()` | `:330` | `:646-649` | yes | public |
| `GetHartle{Omega,DOmega}Coeff_{N_Fast,Mixed,Mixed_Fast}` | `:352-367` | `:300-350` | yes | public |
| `Solve(const Zaki::Math::Axis&, dir)` | `:370-371` | **none** | **absent** | public — **DECLARED, UNDEFINED** |
| `Solve(const double&, dir)` | `:374-375` | **none** | **absent** | public — **DECLARED, UNDEFINED** |
| `Solve_Mixed(const double&, dir)` | `:378-379` | `:363-627` | yes | public |
| `FindNMomInertia()` | `:382` | `:673-752` | yes | public |
| `FindMixedMomInertia()` | `:385` | `:905-999` | yes | public |
| `SolveHartle2_N()` | `:390` | `:1126-1270` | yes | public |
| `SolveHartle2_Mixed()` | `:394` | `:1274-1277` (empty body) | yes | public |
| `GetHartleResult()` | `:397` | `:1003-1006` | yes | public |
| `ExportResults(dir) const` | `:401` | `:632-642` | yes | public |
| `static ODE(...)` | `:403` | **none** | **absent** | public — **DECLARED, UNDEFINED** |
| `static ODE_Mixed`, `ODE_Mixed_Out`, `ODE_Mixed_Fast`, `ODE_N_Fast` | `:404-407` | `:236-296` | yes | public |
| `static ODE_Hartle2_N_Fast`, `ODE_Hartle2_Mixed_Fast` | `:410`, `:413` | `:1024-1110`, `:1114-1120` (stub) | yes | public |
| `Reset()` | `:416` | `:665-667` (empty) | yes | public |
| `SetFastProfilePtrs_`, `SetFastMixedPtrs_`, `EvalFastPEM_`, `EvalFastMixedPEM_` | `:255-266` | `:83-212` | yes | private |
| data: `init_omega_bar` `:212`; `fast_p/e/m` `:230-232`; `fast_nu`, `fast_nu_prime`, `fast_dEdP`, `fast_dEdP_v`, `fast_dEdP_d` `:269-273`; `stored_omega_bar_`, `stored_domega_bar_` `:277-278`; `fast_omega_bar`, `fast_domega_bar` `:281-282`; `hartle_result_` `:285`; `include_m0p0_source_` `:289` | | | | private |

Value types: `TOVTable` `:56-68` (unused), `OmegaSeqPoint` `:79-98`, `HartleResult` `:102-124`.
`fast_nu` (`:269`), `fast_dEdP_v`, `fast_dEdP_d` (`:272-273`) are declared and **never assigned
or read** anywhere in the file.

**Finding — declared-undefined public API.** `Solve(Axis,…)`, `Solve(double,…)`, `ODE`,
`GetMass` and `GetPress` have no definition in any translation unit (`nm -C libCompactStar.a`
lists no such symbol). Any caller fails at link time. Their only reference,
`main/Examples/rotating_ns.cpp:62` (`r_solver.Solve({…}, …)`), is not built:
`main/Examples/CMakeLists.txt:1-5` lists only `sig_omega_rho_nstar` and `Fermi_gas_many`. These
are the owner's *historical* first-order sequence entry points (defined at the common base
`f6bbfea`, removed by the owner's rework `3639d71`, declarations retained).

## 5. Current call / reachability map

| Entry point | Repository callers at `df859b5` | Classification |
|---|---|---|
| `RotationSolver::Solve(Axis,…)` | `main/Examples/rotating_ns.cpp:62` only — file **not compiled**, call would **not link** | **PUBLIC — DECLARED, UNDEFINED (LINK-DEAD) — LEGACY** |
| `RotationSolver::Solve(double,…)` | none | **PUBLIC — DECLARED, UNDEFINED — LEGACY** |
| `RotationSolver::Solve_Mixed(double,…)` | none (the only rotation-solver references in `MixedStar::Find_MomInertia`, `MixedStar.cpp:1288-1290`, are commented out and name the obsolete `Solve_Dark`/`GetOmegaSeq`) | **PUBLIC — ZERO CALLERS — LEGACY (two-fluid)** |
| `RotationSolver::FindNMomInertia()` | `NStar::Find_MomInertia` (`NStar.cpp:1105-1111`, call `:1109`) ← `NStar::BuildFromTOV` `:335` and `NStar::FinalizeSurface` `:659` — i.e. **every** `NStar` construction | **PUBLIC AND USED — LIVE** (all Phase-2/3 harnesses; `SeqPoint::I`) |
| `RotationSolver::FindMixedMomInertia()` | `MixedStar::Find_MomInertia` (`MixedStar.cpp:1280-1300`, call `:1293`) ← `MixedStar.cpp:289`, `:699` | **PUBLIC AND USED — COMPILED, UNEXERCISED** (no live `main/` or test builds a `MixedStar`; `PHASE3_CLOSEOUT.md` §10) |
| `RotationSolver::SolveHartle2_N()` | **zero** in `CompactStar/`, `main/`, `tests/` | **PUBLIC SCIENTIFIC CANDIDATE — ZERO REPOSITORY CALLERS — EXECUTABLE FROM USER CODE — UNVALIDATED** (§13) |
| `RotationSolver::SolveHartle2_Mixed()` | zero | **PUBLIC STUB — CANDIDATE** (empty body `:1274-1277`) |
| `RotationSolver::GetHartleResult()` | zero. `RotochemicalCache::ComputeStructuralDerivative` takes a `HartleResult&` (`RotochemicalCache.cpp:51`, reads `hartle.p0` `:66`) but is **NOT COMPILED** and never calls the accessor | **PUBLIC — ZERO CALLERS** |
| `RotationSolver::GetInitOmegaBar()`, `GetOmegaSeq()`, `ExportResults()` | zero compiled callers (`rotating_ns.cpp:67` uses `ExportResults`, not built) | **PUBLIC — ZERO CALLERS — LEGACY** |
| `NStar::Find_MomInertia()` (`NStar.hpp:362`, `[[nodiscard]]`) | `NStar.cpp:335`, `:659` (internal); no test calls it directly — tests read `GetSequence().I` (`hartle_moment_inertia_analytic.cpp:281,326`; `…_cmf.cpp:209,314`) | **PUBLIC AND USED — LIVE** |
| `MixedStar::Find_MomInertia()` (`MixedStar.hpp:607`) | `MixedStar.cpp:289`, `:699` | COMPILED, UNEXERCISED |

**Construction sites of `RotationSolver`.** `NStar::rot_solver` (`NStar.hpp:100`, in the
`private:` block opened at `:86`; attached in all three constructors, `NStar.cpp:32,40,50`);
`MixedStar::rot_solver` (`MixedStar.hpp:253`; attached `MixedStar.cpp:37,427`);
`main/Examples/rotating_ns.cpp:41` (not compiled). **No test constructs one.** `NStar` declares
`friend class RotationSolver` (`NStar.hpp:82`), which is what lets *any* `RotationSolver`
instance — not only the member — read `nstar_ptr->prof_` (`RotationSolver.cpp:693-695`) and write
`nstar_ptr->MomI` (`:741`). The same friendship exists for `MixedStar` (`MixedStar.hpp:179`).

**What is and is not reachable.** `NStar`'s **own** `rot_solver` and therefore **its**
`hartle_result_` are unreachable (private member, no accessor — the true part of the old
wording). A **second** `RotationSolver` constructed by user code and attached to the same
`NStar` reaches everything, and its `FindNMomInertia` overwrites `NStar::MomI` as a side effect
(`:741`) after `SeqPoint::I` has already been captured (`NStar.cpp:335`).

**Other current consumers of first-order output.** `SeqPoint::I` is exported by the TOV
sequence workflow with the label `I(km^3)` (`TOVSolver.cpp:119`, `:179`). `DriverContext`
exposes no Hartle quantity and no thermal driver reads one (`PHASE2B_CLOSURE.md` §4 —
re-checked: unchanged). `SpinState::Omega()` is a separate, physical spin variable in **rad/s**
(`Physics/State/SpinState.hpp:65,117`; set as `spin.Omega() = 100.0; // rad/s` in
`main/Test/spin_therm_evol_2_main.cpp:205`) that never touches `RotationSolver`.

## 6. Provenance / lineage map

```
f6bbfea (2026-01-06, owner)  ── common base: RotationSolver.cpp 1410 lines
  ├── 675b4a9 (2026-04-05, author "Claude")   AI CANDIDATE: +HartleResult, stored_omega_bar_,
  │      │                                     ODE_Hartle2_N_Fast, SolveHartle2_N, stubs,
  │      │                                     RotochemicalCache/Rotochemical (+964 lines, no tests,
  │      │                                     no CMake change)                       cpp 1707 lines
  │      └── e60e656 (2026-04-05, owner)      "Merge pull request #1 …/claude/hartle-rotochemical-heating-FwAHP"
  │                                            (branch still on origin)
  └── 3639d71 (2026-04-07, owner, "updates")  OWNER FIRST-ORDER REWORK: kR_EPS_KM/SafeR0,
         │                                     profile-backed EvalFastPEM_ at the true RHS radius,
         │                                     cached brackets, master-grid FindMixedMomInertia;
         │                                     deleted the legacy Solve()/ODE()/GetMass/GetPress
         │                                     definitions (declarations kept)         cpp  981 lines
         └── 9f70f14 (2026-04-07, owner)  MERGE of 3639d71 + e60e656 — malformed (candidate side
                │                          taken for three hunks; brace imbalance)       cpp 1265 lines
                └── 57334d8 (2026-08-31, Phase 1B)  MERGE REPAIR: restored the owner's ten
                       │                             interpolation members / four declarations /
                       │                             FindNMomInertia setup; kept the candidate's
                       │                             storage + HartleResult population; re-commented
                       │                             the resurrected dead function     cpp 1279 lines
                       └── … Phase 0.5 / 1 / 2 / 3 … ── df859b5 (current master)
```

Ancestry, verified with `git merge-base --is-ancestor`: `675b4a9 ⊂ e60e656 ⊂ 9f70f14`;
`3639d71 ⊂ 9f70f14`; `9f70f14 ⊂ 57334d8 ⊂ df859b5`; `675b4a9` and `3639d71` are **not** ancestors
of each other (merge base `f6bbfea`). **`RotationSolver.hpp/.cpp` have not changed since
`57334d8`** (`git log` on both files). In particular `git diff a24fe95 df859b5` and
`git diff 10c397d df859b5` on both files are **empty**: the source validated in Phase 2B-4B is
byte-identical to the current source (§7).

Per-block provenance of the **current** file (`git blame -w`, line ranges at `df859b5`).
Provenance is recorded to locate authorship; **it confers no scientific authority either way**
(`GOVERNANCE.md` §1 rule 7, §5).

| Current lines (`RotationSolver.cpp`) | Content | Provenance |
|---|---|---|
| `:28-40` | `kR_EPS_KM`, `SafeR0` | **OWNER FIRST-ORDER** (`3639d71`) |
| `:83-212` | `SetFastProfilePtrs_`, `SetFastMixedPtrs_`, `Lerp_`, `EvalFastPEM_`, `EvalFastMixedPEM_` | OWNER FIRST-ORDER (`3639d71`) |
| `:236-296` | `ODE_Mixed`, `ODE_N_Fast`, `ODE_Mixed_Fast`, `ODE_Mixed_Out` (`r_safe` lines `:244,260,276,292` from `3639d71`) | OWNER FIRST-ORDER (`7fd0132` base + `3639d71`) |
| `:300-350` | coefficient helpers; `_N_Fast`/`_Mixed_Fast` bodies route through `EvalFast*` (`:316-319,326-329,336-339,346-349`) | OWNER FIRST-ORDER |
| `:354-357` | `GetInitOmegaBar` | OWNER (`7fd0132`) |
| `:363-627` | `Solve_Mixed` (two-fluid, `× LIGHT_C_KM_S` exports at `:515,540-542,556,591-592`) | OWNER FIRST-ORDER — LEGACY (`7fd0132` + `3639d71` edits) |
| `:632-667` | `ExportResults`, `GetOmegaSeq`, `GetNStar`, `GetMixedStar`, `Reset` | OWNER — LEGACY (`7fd0132`) |
| `:671-672` | comment "Also stores the omega_bar profile for use by the second-order solver" | **AI CANDIDATE** (`675b4a9`) |
| `:675-699` | `FindNMomInertia` setup — `i0` scan, `r0`, `P/E/M`, `SetFastProfilePtrs_` | **MERGE REPAIR** (`57334d8`, restoring `3639d71`) |
| `:701-713` | seed `init_omega_bar = 5e-3`, `y[]`, GSL driver | OWNER FIRST-ORDER (`7fd0132`; `:702` `3639d71`) |
| `:715-718`, `:729-732` | `stored_omega_bar_` / `stored_domega_bar_` allocation and per-node storage | **AI CANDIDATE** (`675b4a9`; `N→n` by `57334d8`) |
| `:720-728`, `:733-742` | grid loop, `J`, `Ω`, `I` extraction, `nstar_ptr->MomI` | OWNER FIRST-ORDER (`7fd0132`, `:720` `57334d8`) |
| `:743-749` | `HartleResult` population (`Omega`, `J`, `I`, `omega_bar`, `domega_bar`, `r_grid`) | **AI CANDIDATE** (`675b4a9`) |
| `:754-757` | re-commented header of the obsolete `FindMixedMomInertia` | **MERGE REPAIR** (`57334d8`) |
| `:758-902` | fully commented-out obsolete mixed implementation | OWNER — DEAD (commented) |
| `:905-999` | `FindMixedMomInertia` (master-grid, "Added on Jan 7, 2026") | OWNER FIRST-ORDER (`3639d71`) — two-fluid |
| `:1003-1006` | `GetHartleResult` | **AI CANDIDATE** |
| `:1024-1110` | `ODE_Hartle2_N_Fast` | **AI CANDIDATE** — three long lines (`:1068`, `:1088`, `:1100`) are blamed to the merge `9f70f14` only because the merge re-joined wrapped lines; `diff -w` against `675b4a9` is **empty** for the whole function |
| `:1114-1120`, `:1274-1277` | MixedStar second-order stubs | **AI CANDIDATE** |
| `:1126-1270` | `SolveHartle2_N` | **AI CANDIDATE** — `diff -w` against `675b4a9` empty |
| `RotationSolver.hpp:100-125`, `:268-290`, `:387-398`, `:409-414`, `:45` | `HartleResult`; second-order members; declarations; `DataColumn` include | **AI CANDIDATE** |
| `RotationSolver.hpp:227-266` | comment + `fast_p/e/m` declarations (kept for the candidate) + the owner's ten interpolation members and four private declarations | **MERGE REPAIR** (`57334d8`) |
| everything else in the header | | OWNER (`7fd0132`, `dfb4443`, `eb7608f`) |
| `tests/rotation/hartle_reference.hpp`, `hartle_moment_inertia_{analytic,cmf}.cpp` | Phase 2B-4B harnesses (`10c397d`) | **LATER VALIDATION-ONLY** — no production code |
| `MixedStar.hpp/.cpp` rotation hooks | `Find_MomInertia`, master-grid totals | OWNER (`3639d71`) — out of core scope (§14) |

**UNKNOWN: none.** Every line of both files resolves to one of the four lineages above.

## 7. First-order equation — re-authenticated, not re-derived

**7.1 Source identity.** `RotationSolver.{hpp,cpp}` are byte-identical to the source validated
in Phase 2B-4B (`git diff a24fe95 df859b5 -- …` empty; last change `57334d8`, 2026-08-31).
`NStar.cpp` changed in Phase 3 (3D, 3E-I2), but the Hartle golden
`hartle_I_dscmf1_debug.tsv` is hash-identical (§2) and both Hartle CTests pass at `df859b5`, so
the profile the solver consumes is equivalent. **`HARTLE_MOMENT_INERTIA.md` §3–§4
(`EQUATION MATCH`) stands unchanged**; nothing in it was re-run except what §12–§13 below
needed anyway.

**7.2 What is confirmed against the current source.**

- `ω̄ = Ω − ω(r)` — `RotationSolver.hpp:72`; state `y[0] = ω̄ [km⁻¹]`, `y[1] = dω̄/dr [km⁻²]`
  (`ODE_N_Fast`, `RotationSolver.cpp:252-264`), coefficients `16π(p+ε)/(1−2m/r)` and
  `4π(p+ε)r/(1−2m/r)` (`:314-330`) with `p, ε, m` interpolated at the true RHS radius
  (`EvalFastPEM_`, `:114-164`).
- Surface extraction (`:737-739`): `J_raw = R⁴ ω̄'(R)/6`, `Ω_raw = ω̄(R) + R ω̄'(R)/3`,
  `I = J_raw/Ω_raw`; `R` = last profile radius (`:681`).
- Linear homogeneity: under `ω̄ → A ω̄`, `J → AJ`, `Ω → AΩ`, `I` unchanged
  (`HARTLE_MOMENT_INERTIA.md` §3.1, §7). Re-measured here with the independent reference solver
  on the analytic star: seed × 37 gives `J` and `Ω` proportional to `6.3e-15` and `6.1e-15`,
  `I` changed by `1.8e-16` (§12, scratch run).
- The raw seed is `init_omega_bar = 5e-3` (`:701`), assigned as `y[0]` at the first
  strictly-positive grid radius (`:684-706`), `y[1] = 0`.

**7.3 INV-05 rotation half — re-audited (discharges the ⚠).** `kR_EPS_KM = 1e-6 km` (`:31`);
`SafeR0()` (`:33-38`) used by `Solve_Mixed` (`:375`) and `FindMixedMomInertia` (`:947`);
`r_safe` guards in the four first-order RHS functions (`:244`, `:260`, `:276`, `:292`);
`FindNMomInertia` starts at the first strictly-positive grid radius, falling back to
`kR_EPS_KM` (`:684-690`) — on production profiles that is `r₀ = 1e-5 km` (INV-05 TOV half;
measured `r0 = 1.000000e-05 km` on both audited stars). Centre condition `ω̄'(r₀) = 0`; **no
series expansion**. The regular series is `ω̄ = ω̄_c[1 + (8π/5)(ε_c+p_c) r² + …]`, so starting at
`r₀` with the centre value carries a relative truncation `O((ε_c+p_c) r₀²) ≲ 1e-12`, and
`HARTLE_MOMENT_INERTIA.md` §15.1 shows `I` is immune to the centre condition by construction.
`SolveHartle2_N` starts at `r_min = r[0]` with literal zeros (`:1135`, `:1170-1171`,
`:1204-1205`) — assessed in §11D. **INV-05 rotation half: VERIFIED CURRENT BEHAVIOR at
`df859b5`.**

## 8. INV-07 — the physical normalization contract, derived

**8.1 Derivation.** Let `L[ω̄] = 0` denote the first-order equation in its production form
(`ω̄'' = −4ω̄'/r + [16π(ε+p)/(1−2m/r)] ω̄ + [4π(ε+p)r/(1−2m/r)] ω̄'`, `EQUATION MATCH` with
Hartle's `(r⁴ j ω̄')' /r⁴ + (4/r) j' ω̄ = 0`).

1. **One-dimensional solution space.** `L` is linear and homogeneous; with the regularity
   condition `ω̄'(0) = 0` the regular solutions form a one-dimensional space. Any regular
   solution is `A·ω̄_raw(r)` for the production solution `ω̄_raw` obtained from the arbitrary
   seed. (The irregular solution `∝ r⁻³` is excluded by regularity; 2B-4B §15.1 shows it decays
   across the star even when excited.)
2. **Exterior matching is linear.** Outside the star `ε = p = 0`, `j = 1`, so `(r⁴ω̄')' = 0`;
   the asymptotically flat exterior is `ω̄ = Ω − 2J/r³`. Matching value and derivative at `R`
   gives `J = R⁴ω̄'(R)/6` and `Ω = ω̄(R) + Rω̄'(R)/3` — both **linear functionals** of the
   solution. Hence `J[Aω̄_raw] = A J_raw` and `Ω[Aω̄_raw] = A Ω_raw`.
3. **Physical target fixes `A` uniquely.** A requested physical angular velocity `Ω_phys [s⁻¹]`
   converts to geometric units by `Ω_target_geom = Ω_phys / c`, `c = LIGHT_C_KM_S = 2.99792458e5
   km s⁻¹` (exact by definition; the value is what `Zaki::Physics::LIGHT_C_KM_S` holds, printed
   in the scratch run). Requiring `Ω[Aω̄_raw] = Ω_target_geom` gives
   **`A = Ω_target_geom / Ω_raw`**, well defined because `Ω_raw = ω̄(R) + 2J_raw/R³ > 0` for a
   positive seed (both terms positive: `ω̄ > 0` throughout and `ω̄' > 0`).
4. **Consequences.** `ω̄_phys(r) = A ω̄_raw(r)`, `dω̄_phys/dr = A dω̄_raw/dr`, `J_phys = A J_raw`,
   `Ω_geom = A Ω_raw = Ω_target_geom`, and `I = J_phys/Ω_geom = J_raw/Ω_raw` — **unchanged**,
   which is why `I` was validatable before this contract existed. The frame-dragging angular
   velocity `ω(r) = Ω − ω̄(r) = A(Ω_raw − ω̄_raw(r))` scales too, and the ratio `ω(r)/Ω` is
   seed-independent (a property of the background star).
5. **Zero spin.** `Ω_phys = 0 ⇒ A = 0 ⇒ ω̄_phys ≡ 0`, `J_phys = 0`, while `I = J_raw/Ω_raw`
   remains defined (it is the `Ω → 0` limit of `J/Ω`). An implementation must therefore form `I`
   from the raw solution, never as `J_phys/Ω_phys`.
6. **Regime.** The contract is exact for the linear first-order problem at any `A`; physical
   validity of the *slow-rotation expansion* additionally requires `Ω ≪ Ω_K`. In geometric
   units `Ω_K ≈ (M/R³)^{1/2} [km⁻¹]` needs no `G` — a useful diagnostic (ADR-0006 §7 item 9).

**8.2 Dimensions.** `ω̄, Ω [km⁻¹]`; `ω̄' [km⁻²]`; `J = R⁴ω̄'/6 [km²]`; `I [km³]`; `A`
dimensionless; `Ω_phys = c Ω_geom [km s⁻¹ · km⁻¹ = s⁻¹]` ✓. Only `c` enters the Ω/ω̄
conversion. Converting `J` or `I` to cgs additionally needs `G` (`J_cgs = J_km² · 10¹⁰ · c³/G`,
`I_cgs = I_km³ · 10¹⁵ · c²/G`), which is the still-unadjudicated solar-mass/`G` authority
(`PHASE3_CLOSEOUT.md` §12) — ADR-0006 therefore fixes the `c`-only part and defers cgs `J`, `I`
to that ADR.

**8.3 Numerical illustration (scratch run, production `FindNMomInertia` through the public API).**

| Star | `Ω_raw` [km⁻¹] | `Ω_raw·c` [s⁻¹] | `J_raw` [km²] | `I` [km³] | `ω_c/Ω` |
|---|---|---|---|---|---|
| analytic uniform, `M = 2 km`, `R = 13 km`, `N = 4001` | `7.310742328681e-3` | `2.1917e3` | `1.148757927677` | `1.571328705116e2` (= `SeqPoint::I` bitwise) | `0.3160749` |
| DS(CMF)-1 `1.400022 M☉`, `R = 13.54532 km`, `N = 2646` | `8.064116072884e-3` | `2.4176e3` | `1.093624218247` | `1.356161305669e2` (= `SeqPoint::I` bitwise) | `0.3799692` |

So the hard-coded seed corresponds, physically, to a star spinning at **≈ 2.2–2.4 × 10³ s⁻¹
(≈ 350–385 Hz)** — a fast millisecond pulsar — which is why the candidate's O(Ω²) outputs are of
the magnitude they are (§12) and why they must never be read as physical.

| Requested `Ω_phys` | `Ω_geom = Ω_phys/c` | `A` (analytic / CMF) | `J_phys` (analytic / CMF) [km²] | `J_phys/Ω_geom − I` (rel) |
|---|---|---|---|---|
| `100 s⁻¹` | `3.335641e-4 km⁻¹` | `4.562657e-2` / `4.136400e-2` | `5.241388e-2` / `4.523667e-2` | `−1.8e-16` / `0` |
| `2π·716 Hz = 4498.761 s⁻¹` | `1.500625e-2 km⁻¹` | `2.052630` / `1.860867` | `2.357975` / `2.035090` | `−1.8e-16` / `−2.1e-16` |

The contract is arithmetic once `ω̄_raw` exists; its only numerical hazard is indirect — the
production driver uses a **fixed absolute tolerance `1e-10`** (`:712-713`), so the *accuracy* of
`ω̄_raw` depends on the seed's magnitude relative to that tolerance. That is a reason to keep
the seed internal and stable, not a reason to expose it (ADR-0006 §5-Q2, §7 item 1).

## 9. Unit-annotation audit

| Quantity | Documented unit | Actually stored | Dimension | Consumers | Verdict |
|---|---|---|---|---|---|
| `init_omega_bar` (`hpp:212`) | "Initial bar{omega} (at r = r(i=0))" — **no unit** | `5e-3` assigned to `y[0]` (`cpp:701-705`); `Solve_Mixed` takes it from the caller (`:389`) | km⁻¹ | `FindNMomInertia`, `FindMixedMomInertia` (`:954`), `Solve_Mixed`, `GetInitOmegaBar` (zero callers) | **UNDOCUMENTED UNIT; arbitrary numerical seed** — must be documented and must not escape |
| `omega_bar`, `domega_bar` (`HartleResult`, `hpp:109-110`) | `[km⁻¹]`, `[km⁻²]` | raw geometric, seed-normalized (`:747-748`) | km⁻¹, km⁻² | none compiled | annotation **correct**; normalization arbitrary (INV-07) |
| `OmegaSeqPoint::omega_bar_c` (`hpp:81`) | struct comment: value at r = 0, no unit; export header `omega_bar_c (1/s)` (`cpp:637-638`) | populated **only** by `Solve_Mixed`: `init_omega_bar × LIGHT_C_KM_S` (`:540`) | s⁻¹ | `ExportResults` | label consistent with the one populating path; **value is the caller's seed × c, not a physical spin**; never populated on the `NStar` path |
| `OmegaSeqPoint::J` | header `J` — **no unit** | `ang_mom_J` geometric (`:542`) | km² | `ExportResults` | **unit missing in export header** |
| `OmegaSeqPoint::Omega` | header `Omega (1/s)` | `ang_vel_Omega × LIGHT_C_KM_S` (`:542`) | s⁻¹ | `ExportResults` | consistent for the mixed path; seed-normalized |
| `OmegaSeqPoint::m`, `::r` | header `M`, `R` — no units | `sequence.v.m` (M☉), `r_surface` (km) | M☉, km | `ExportResults` | **units missing in export header** |
| `HartleResult::Omega` (`hpp:105`) | `[s^-1]` | `ang_vel_Omega` **without** `× c` (`:744`) | km⁻¹ | none | **WRONG** (INV-07 secondary defect, re-confirmed: value `7.31e-3` and `8.06e-3` on the audited stars) |
| `HartleResult::J` (`hpp:106`) | `[km^2]` | `R⁴ω̄'/6` | km² | none | correct |
| `HartleResult::I` (`hpp:107`) | `[km^3]` | `J/Ω` (`:746`); also `NStar::MomI` → `SeqPoint::I` | km³ | tests; `_Sequence.tsv` `I(km^3)` | correct |
| `stored_omega_bar_`, `stored_domega_bar_` (`hpp:277-278`) | "Cached omega_bar profile … for interpolation" | raw, seed-normalized (`:731-732`) | km⁻¹, km⁻² | `SolveHartle2_N` (`:1184-1185`) | **the candidate's only source; quadratic seed dependence** (§12) |
| `m0` (`hpp:113`) | `delta_m(r) [km]` | `y[0]` of `ODE_Hartle2_N_Fast` | km (from `dm0/dr = 4πr²(dε/dp)p0`, `:1054`) | — | intended unit consistent; but one source term is dimensionally km⁻¹ (§10) |
| `p0` (`hpp:114`) | `delta_p(r) [km^-2]` | `y[1]` | intended km⁻² (Eulerian δp) | `RotochemicalCache.cpp:66` (not compiled) | **intended variable = Eulerian pressure perturbation δp**, not Hartle's dimensionless `p₀*`; its ODE is dimensionally inconsistent with either (§10) |
| `xi0` (`hpp:115`) | `[km]` | `−p0/(dp/dr)` (`:1258-1261`) | km | — | relation correct in form for `p0 ≡ δp` |
| `delta_M` (`hpp:118`) | `[km]` | `m0[-1]` (`:1268`) | km | — | incomplete (no `J²/R³`); sign/magnitude unphysical (§12) |
| `p0_c` (`hpp:117`) | "Central pressure perturbation (shooting result)" — no unit | coefficient of the homogeneous solution with `p0_hom(r₀) = 1` (`:1205`, `:1239`) = `δp(r₀)` | km⁻² | — | **unit missing**; nonzero by construction ⇒ central density changes (§11E) |
| `fast_nu`, `fast_dEdP_v`, `fast_dEdP_d` (`hpp:269,272-273`) | "Metric nu at current grid point", … | **never assigned** | — | none | dead members |
| profile `ν'` column consumed by the candidate (`:1133`) | `StarProfile::Column::MetricNuPrime` "dν/dr" | converted `1/cm → 1/km` at `NStar.cpp:187`, `:721` | km⁻¹ | `SolveHartle2_N` | correct |
| `ExportResults` header (`:637-638`) | `omega_bar_c (1/s)  M  R  J  Omega (1/s)` | see rows above | mixed | `rotating_ns.cpp:67` (not built) | **advertises physical units for a seed-normalized quantity**; `M`, `R`, `J` unlabeled; exports nothing on the `NStar` path (`omega_seq_pts` never populated there) |

No label was changed. The audit is recorded for ADR-0006 (`HartleResult` and export semantics).

## 10. Second-order candidate — equation and variable map (not fixed)

**10.1 Reference system, in the repository's metric convention** (`g_tt = −e^{2ν}`,
`g_rr = e^{2λ}`, INV-03; Hartle's own `e^{ν_H}` is `e^{2ν}` here), with
`j ≡ e^{−(ν+λ)} = e^{−ν}(1 − 2m/r)^{1/2}`, so `j² = e^{−2ν}(1 − 2m/r)` and, using the identity
validated in `HARTLE_MOMENT_INERTIA.md` §4 (`j'/j = −(ν'+λ') = −4πr(ε+p)/(1−2m/r)`),

```
dj²/dr = −8π r² (ε+p) j² / (r − 2m)                                              (J)
```

Hartle's l = 0 sector is written in the variables `m₀(r)` [km] (mass perturbation) and the
dimensionless **pressure perturbation factor** `p₀*(r)`, defined through the Eulerian pressure
perturbation `δp(r,θ) = (ε+p)[p₀*(r) + p₂*(r)P₂(cosθ)]`, i.e. **`δp₀ = (ε+p) p₀*`**:

```
dm₀/dr  = 4π r² (ε+p)(dε/dp) p₀*  +  (1/12) j² r⁴ (dω̄/dr)²  −  (1/3) r³ (dj²/dr) ω̄²          (H1)

dp₀*/dr = − m₀ (1 + 8π r² p)/(r − 2m)²  −  4π (ε+p) r² p₀*/(r − 2m)
          + (1/12) r⁴ j² (dω̄/dr)²/(r − 2m)  +  (1/3) d/dr[ r³ j² ω̄² /(r − 2m) ]              (H2)
```

with, by (J), `−(1/3) r³ (dj²/dr) ω̄² = +(8π/3) r⁴ (ε+p) e^{−2ν} ω̄²`. Dimensions: every term of
(H1) is dimensionless, every term of (H2) is km⁻¹ ✓. Regular centre (derived from (H1)–(H2), with
`m ∝ r³`, `dω̄/dr = (16π/5)(ε_c+p_c) ω̄_c r + …`):

```
p₀* = (1/3) j_c² ω̄_c² r² + O(r⁴),     m₀ = (4π/15)(ε_c+p_c)[(dε/dp)_c + 2] j_c² ω̄_c² r⁵ + O(r⁷)
```

Boundary/normalization for the **fixed-central-density** family (the one INV-09 and
`RotochemicalCache` presuppose): `m₀(0) = 0`, `p₀*(0) = 0` — **no** homogeneous admixture; the
homogeneous solution (`p₀*(0) = const`) is a shift of central density along the nonrotating
sequence and is added only when comparing at fixed baryon mass. Derived quantities:
`ξ₀ = −p₀*(ε+p)/(dp/dr) = −δp₀/(dp/dr) = p₀*/ν'` (using the TOV relation `dp/dr = −(ε+p)ν'`);
**`δM = m₀(R) + J²/R³`** (exterior solution `m₀(r) = δM − J²/r³`); `dε/dp` is the barotropic
EOS derivative `1/c_s²`.

> **Transcription caveat.** (H1)–(H2) are transcribed from the literature listed in the header,
> in this repository's conventions. The candidate's own comment cites "Hartle (1967),
> Eqs. (30)–(33)" (`RotationSolver.cpp:1013`, `:1087`); those equation numbers were **not**
> verified against the primary source here. Owner confirmation of the transcription against
> the primary source is part of the second-order ADR's validation (§15), not of this record.

**10.2 What the code integrates** (`ODE_Hartle2_N_Fast`, `:1024-1110`; state `y[0] = m0`,
`y[1] = p0`; `p, ε, m, ν', dε/dp, ω̄, ω̄'` are per-node scalars set in `SolveHartle2_N`,
`:1179-1185`):

```
dm0/dr = 4π r² (dε/dp) p0 + S_m,                                                       (:1054, :1102)
   S_m  = (1/12) r⁴ (ω̄')² /(1 − 2m/r)  −  (4π/3) r³ (ε+p) ω̄² /(1 − 2m/r)                 (:1088)
dp0/dr = − p0 [4π(ε+p) r + m/r²]/(r − 2m)  −  m0 (ε+p)/(r² (r − 2m)) + S_p,             (:1068, :1103)
   S_p  = (1/12) r (1 − 2m/r) (ω̄')²  +  (1/3) r ν' ω̄²                                    (:1100)
```

**10.3 Term-by-term map.** To compare, (H1)–(H2) are rewritten for `δp ≡ (ε+p)p₀*` using
`(ε+p)' = −(1 + dε/dp)(ε+p)ν'`:

```
dm₀/dr = 4π r² (dε/dp) δp + (1/12) j² r⁴ (ω̄')² + (8π/3) r⁴ (ε+p) e^{−2ν} ω̄²                          (T1)
dδp/dr = −(1+dε/dp) ν' δp − 4π(ε+p) r² δp/(r−2m) − (ε+p) m₀ (1+8πr²p)/(r−2m)²
         + (1/12)(ε+p) r³ e^{−2ν} (ω̄')² + (1/3)(ε+p) d/dr[ r² e^{−2ν} ω̄² ]                          (T2)
```

| Code term | Dimension (code) | Required (for `p0 ≡ δp`) | Verdict |
|---|---|---|---|
| `4π r² (dε/dp) p0` | — (ok) | `4π r² (dε/dp) δp` | **matches** (only term of the homogeneous system that does) |
| `S_m` term 1: `(1/12) r⁴ (ω̄')²/(1−2m/r)` | — | `(1/12) e^{−2ν}(1−2m/r) r⁴ (ω̄')²` | **`1/(1−2m/r)` is the reciprocal of the required factor; `e^{−2ν}` missing** (INV-08 statement re-confirmed) |
| `S_m` term 2: `−(4π/3) r³ (ε+p) ω̄²/(1−2m/r)` | **km⁻¹** ✗ | `+(8π/3) r⁴ (ε+p) e^{−2ν} ω̄²` (dimensionless) | **wrong sign, factor 2, one power of `r` short, wrong metric factor** |
| `dp0/dr` homogeneous `p0` coefficient `[4π(ε+p)r + m/r²]/(r−2m)` | **km⁻²** ✗ (needs km⁻¹) | `(1+dε/dp)ν' + 4π(ε+p)r²/(r−2m)` | **dimensionally inconsistent by `1/r`** — the comment at `:1057` writes `r²(r−2m)` where the TOV factor is `r²(1−2m/r) = r(r−2m)`; the `(dε/dp) ν' δp` piece is absent |
| `dp0/dr` homogeneous `m0` coefficient `(ε+p)/(r²(r−2m))` | **km⁻⁵** ✗ (needs km⁻⁴) | `(ε+p)(1+8πr²p)/(r−2m)²` | **dimensionally inconsistent by `1/r`; GR factor `(1+8πr²p)/(1−2m/r)` absent** |
| `S_p` term 1: `(1/12) r (1−2m/r)(ω̄')²` | km⁻³ (ok) | `(1/12)(ε+p) r³ e^{−2ν} (ω̄')²` | **different structure: no `(ε+p)`, no `r³e^{−2ν}`** |
| `S_p` term 2: `(1/3) r ν' ω̄²` | **km⁻²** ✗ | `(1/3)(ε+p) d/dr[r² e^{−2ν} ω̄²]` | **dimensionally inconsistent; not the derivative term** |
| `ξ₀ = −p0/(dp/dr)`, `dp/dr = −(ε+p)ν'` (`:1252-1261`) | km | `−δp/(dp/dr)` | **correct in form** for `p0 ≡ δp` |
| `delta_M = m0[-1]` (`:1268`) | km | `m₀(R) + J²/R³` | **exterior term missing** |
| shooting `p0(R) = 0` (`:1233-1239`) | — | `p₀*(0) = 0`, no admixture (fixed ε_c) | **wrong boundary condition** (§11E) |

**10.4 Conclusion on the variable `p0`.** By the header annotation (`[km⁻²]`, `hpp:114`), the
one correct homogeneous term (`:1054`), the `ξ₀` relation (`:1258-1261`) and the way
`RotochemicalCache` consumes it (`dn_i/dp × p0`, `RotochemicalCache.cpp:98-99`), **the code's
`p0` is intended to be the Eulerian monopole pressure perturbation `δp₀ = (ε+p) p₀*` in km⁻²
— not Hartle's dimensionless `p₀*`.** Its evolution equation, however, is not the correct
equation for `δp` (T2) *or* for `p₀*` (H2): the homogeneous part is dimensionally inconsistent
under either reading. **`p0` is therefore "δp by intent, with an invalid ODE"**; the Phase-4C
ADR must choose the variable explicitly (recommendation: Hartle's `p₀*`, with `δp` derived) and
must not inherit the candidate's equations.

**10.5 A false premise in the candidate.** The comment at `:1078-1080` says `e^{−2ν}` "is not
available in the fast cache". The profile **does** carry `ν` (`StarProfile::Column::MetricNu`,
`StarProfile.hpp:246`; read by the tests via `Profile().GetMetricNu()`), so the `j²` omission had
no data-availability justification.

## 11. Known INV-08 defects — re-authenticated against current source and literature

| | Item | Finding at `df859b5` |
|---|---|---|
| **A** | `j²` factor | Required: `j² = e^{−2ν}(1−2m/r)` multiplying `r⁴(ω̄')²` in (H1), and `dj²/dr` (equivalently `e^{−2ν}`) in the `ω̄²` source. Production uses **only `1/(1−2m/r)`** in `S_m` (`:1088`) — the reciprocal of the `(1−2m/r)` part, with `e^{−2ν}` absent; `S_p` uses `(1−2m/r)` in one term and nothing in the other (`:1100`). Neither `j²` nor `1/j²` is used consistently anywhere |
| **B** | exterior mass correction | Required in the code's own convention (`m₀ [km]`, `J [km²]`): `δM = m₀(R) + J²/R³`. Production: `m0[-1]` only (`:1268`). With the raw seed, `J²/R³` would also be seed²-dependent; the correction is only meaningful after INV-07 or per `Ω²` |
| **C** | `dε/dp` | Correct authority: the **EOS** barotropic derivative `dε/dp = 1/c_s²` (the TOV/EOS spline already owned by `TOVSolver`; not a profile derivative). The profile finite difference `(dε/dr)/(dp/dr)` (`:1140-1161`) equals it by the chain rule wherever `dp/dr ≠ 0` (the profile is a barotrope), so it is **mathematically equivalent but numerically inferior**: on the CMF star it ranges `3.79 – 4.87e3`, amplifying discretization noise across the crust; on the constant-density star it is identically `0` (correct: incompressible ⇒ `dε/dp = 0`). The `1.0` fallback (`:1148`, `:1154`, `:1160`) is `c_s² = 1` — the **causal limit**, not an incompressible limit; INV-08's parenthetical "(incompressible is dε/dp → ∞)" is itself wrong and is corrected. No fallback was triggered on either audited star |
| **D** | centre behaviour | Starts at `r₀ = r[0] = 1e-5 km` with `{m0, p0} = {0, 0}` (particular) and `{0, 1}` (homogeneous) (`:1170-1171`, `:1204-1205`). Against the regular series (§10.1: `p₀* ∝ r²`, `m₀ ∝ r⁵`; the homogeneous solution with `p₀*(0) = 1` is itself regular) the literal-zero start has relative truncation `O((r₀/R)²) ≈ 6e-13` — **an acceptable documented approximation for the correct equations**, downgraded from "defect" to "no series expansion; bound recorded". The `r < 1e-10` guard (`:1042`) is never reached |
| **E** | surface condition / superposition | The superposition (`:1233-1250`) is arithmetically correct **for the condition it imposes**: `p0(R) = 0` is achieved (measured `−1.06e-21`, `5.0e-22`). That condition is **not Hartle's**: it forces `ξ₀(R) = −δp(R)/(dp/dr) = 0` (measured `−5e-10 km`, `1e-16 km`) — no monopole surface displacement — and yields `p0_c = δp(r₀) ≠ 0` (measured `−0.47 p_c` and `−0.52 p_c` at the raw seed), i.e. **the central density changes**. Hartle's fixed-central-density solution has `p₀*(0) = 0`, `p₀*(R) ≠ 0`, `ξ₀(R) = p₀*(R)/ν'(R) ≠ 0`. The imposed condition confuses the Eulerian perturbation `δp` with the Lagrangian one `Δp = δp + ξ₀ dp/dr`, whose vanishing is already the *definition* of `ξ₀`. Consequently the candidate does **not** produce the fixed-`ε_c` derivatives INV-09 requires |
| **F** | `ξ₀` | `ξ₀ = −p0/(dp/dr)` with `dp/dr = −(ε+p)ν'` from the profile `ν'` column (km⁻¹, `NStar.cpp:187,721`): **correct in form** for `p0 ≡ δp`; but because E forces `δp(R) = 0`, the physical monopole radius change is lost |
| **G** | `δM` | dimension km ✓; matching incomplete (B); **measured sign and magnitude unphysical**: `δM/M = −1.595` on the CMF star and `−3.4e-6` on the analytic star at the raw seed (`Ω ≈ 2.2–2.4e3 s⁻¹`), whereas the physical fixed-`ε_c` monopole mass change is **positive** and small (percent level or below for `Ω ≲ 2.5e3 s⁻¹`) |
| **H** | stored first-order source | **Confirmed**: `SolveHartle2_N` reads `stored_omega_bar_[i]`, `stored_domega_bar_[i]` (`:1184-1185`) written by `FindNMomInertia` at seed `5e-3` (`:729-732`). Under `ω̄ → Aω̄`: all sources `∝ A²`; the homogeneous pass is `A`-independent; `p0_c = −p0_part(R)/p0_hom(R) ∝ A²`; hence `m0, p0, ξ₀, p0_c, δM ∝ A²` (measured, §12) |

Also re-confirmed: `[FIX: confirm exact from textbook]` (`:1019`) and `???` (`:1091`) are still in
shipped source; `ODE_Hartle2_Mixed_Fast` returns zeros (`:1114-1120`) and `SolveHartle2_Mixed` is
empty (`:1274-1277`).

## 12. Candidate scaling evidence — scratch only, nothing committed

Scratch program (outside the tree; links the built `libCompactStar.a`; production source
untouched): constructs an `NStar`, attaches a **second, external** `RotationSolver` through the
public API, runs `FindNMomInertia()`, then rescales the private `stored_omega_bar_` /
`stored_domega_bar_` by `A` (scratch-only access override), calls `SolveHartle2_N()` and reads
`GetHartleResult()`. Seeds tested: `A ∈ {1, 0.1, 10, 10³, 10⁻³}` relative to the production seed.

**Executability.** Both stars: `valid = true`, all of `m0, p0, xi0, delta_M, p0_c` finite,
`p0(R) = 0` satisfied to `≤ 1e-21`. The candidate is **not** BLOCKED numerically; it runs and
returns numbers. (That is the hazard: they look like results.)

**Raw-seed values (A = 1; `Ω_raw` ≈ 2.2e3 / 2.4e3 s⁻¹).**

| Star | `δM` [km] | `δM/M` | `p0_c` [km⁻²] | `p0_c/p_c` | `ξ₀(N/2)` [km] | `m0(N/2)` [km] |
|---|---|---|---|---|---|---|
| analytic | `−6.806915e-6` | `−3.4e-6` | `−1.256840e-5` | `−0.515` | `−5.006` | `−1.15e-5` |
| CMF 1.4 M☉ | `−3.297389` | **`−1.595`** | `−2.468505e-5` | `−0.469` | `−4.618` | `−3.08e-2` |

**Quadratic-scaling residuals `Q(A)/Q(1)/A² − 1`** (0 = exact `A²`):

| Star | `A` | `δM` | `p0_c` | `m0(N/2)` | `p0(N/2)` | `ξ₀(N/2)` |
|---|---|---|---|---|---|---|
| analytic | `0.1` | `−1.8e-15` | `−1.6e-15` | `2.2e-16` | `−1.9e-15` | `−2.1e-15` |
| analytic | `10` | `−9.4e-4` | `7.4e-6` | `−1.1e-16` | `8.2e-6` | `8.2e-6` |
| analytic | `10³` | `−1.5e-3` | `2.3e-5` | `1.9e-6` | `2.2e-5` | `2.2e-5` |
| analytic | `10⁻³` | `−1.3e-14` | `−6.7e-16` | `2.2e-16` | `−8.9e-16` | `−8.9e-16` |
| CMF | `0.1` | `5.3e-6` | `−1.3e-5` | `−1.3e-5` | `−1.4e-5` | `−1.4e-5` |
| CMF | `10` | `5.2e-10` | `5.4e-6` | `4.9e-6` | `4.3e-6` | `4.3e-6` |
| CMF | `10³` | `3.7e-9` | `6.9e-6` | `5.3e-6` | `4.4e-6` | `4.4e-6` |
| CMF | `10⁻³` | **`4.3e-3`** | `−1.5e-5` | `−1.5e-5` | `−1.5e-5` | `−1.5e-5` |

**Reading.** Every candidate output scales as `A²` to within `1.5e-3` over six decades, exactly
(`≤ 2e-15`) where the driver's **fixed absolute tolerance `1e-10`** (`:1168`, `:1202`) is not
engaged and at the `10⁻⁵ – 10⁻³` level where it is. So (i) the candidate's outputs are pure
functions of the arbitrary seed squared — **normalization dependence confirmed**; (ii) even the
seed-normalized ratios `Q/Ω_raw²` are not seed-independent beyond `~1e-5 – 1e-3`, because the
integrator's absolute tolerance couples to the seed magnitude — a second, numerical reason the
seed must be governed. **Successful `A²` scaling is not correctness**: the values themselves are
unphysical (§11G). No candidate output was written to any baseline.

For the record, the candidate's own "response coefficients" `Q/Ω_raw²` at `A = 1`:
`δM/Ω²` = `−1.273583e-1 km³` (analytic), `−5.070569e4 km³` (CMF); `p0_c/Ω²` = `−0.2351565`,
`−0.379595`. These numbers are **diagnostic only** and must not be cited.

## 13. Reachability correction

**Claim tested:** "the O(Ω²) candidate is structurally unreachable." **Result: false.** From user
code, with public API only:

```cpp
CompactStar::Core::NStar ns(points);              // any profile construction (public ctor)
CompactStar::Core::RotationSolver rs;             // public class, public ctor
rs.AttachNStar(&ns);                              // public   RotationSolver.hpp:306
rs.FindNMomInertia();                             // public   :382   (raw first order, seed 5e-3)
rs.SolveHartle2_N();                              // public   :390   (the O(Omega^2) candidate)
const auto &h = rs.GetHartleResult();             // public   :397   -> m0, p0, xi0, p0_c, delta_M
```

This compiled, linked against the unmodified `libCompactStar.a`, and executed on both audited
stars (§12). No private state is needed: `FindNMomInertia` reaches `NStar` internals through the
`friend` declaration (`NStar.hpp:82`), and `SolveHartle2_N` reads only `Profile()` (public,
`NStar.hpp:249`) plus the solver's own members.

**Correct classification:** `PUBLIC SCIENTIFIC CANDIDATE — ZERO REPOSITORY CALLERS —
EXECUTABLE FROM USER CODE — UNVALIDATED` (GOVERNANCE.md §5 unchanged). What *is* unreachable:
`NStar`'s own `rot_solver.hartle_result_`, and any normalized/physical product.
`SCIENTIFIC_INVARIANTS.md` INV-08, `CURRENT_ARCHITECTURE.md` and `TARGET_ARCHITECTURE.md` are
corrected in this task (§20); the INV-08 **status** is unchanged (`INTENDED BUT UNVERIFIED` /
`UNVERIFIED SCIENTIFIC CANDIDATE`): public reachability is a reason for more caution, not for a
status change.

## 14. MixedStar / two-fluid rotation — classified separately, out of Phase-4 core scope

Ordinary `NStar` Hartle: `FindNMomInertia` on the single `StarProfile` (§5, §7) — **the Phase-4
subject.** `MixedStar` rotation: `FindMixedMomInertia` (`:905-999`) integrates the same
first-order equation over master-grid **total** `p, ε, m` (`SetFastMixedPtrs_`, `:960`;
`MixedStar.hpp:307-339`, `:471-505`) with the same seed `5e-3` (`:954`) and the same extraction
(`:992-994`); `Solve_Mixed` (`:363-627`) is the legacy sequence path with core/mantle two-driver
integration; second order is a stub (`:1114-1120`, `:1274-1277`). `MixedStar` is
`COMPILED — UNEXERCISED` with zero tests (`PHASE3_CLOSEOUT.md` §10) and its geometry migration
is deferred by ADR-0004 §0-Q2. **Classification: `MIXED-STAR ROTATION — COMPILED, UNEXERCISED,
UNVALIDATED — SEPARATE TRACK`.** INV-07's normalization contract applies to it by the same
mathematics, but its validation waits on `MixedStar` coverage. **No ordinary-`NStar` Phase-4
item is blocked by it.**

## 15. Phase-4 internal dependency plan (proposed; ratified item list unchanged)

| Increment | Content | Class | Gate |
|---|---|---|---|
| **4A-0** | this audit; ADR-0006 PROPOSED | documentation | done |
| **4A** | ADR-0006 **accepted**, then implemented: physical spin input, `A = Ω_geom/Ω_raw` rescaling, seed internalized, `HartleResult`/export unit contract, normalized first-order response exposed (`I`, `ω̄/Ω`, `J_phys`, `Ω`) | scientific-semantic + structural (API) | ADR-0006 acceptance; `I` and the seven goldens byte-identical (ADR-0005 §13) |
| **4B** | physical first-order validation — the nine tests and four detectors of §18, reusing `hartle_moment_inertia_{analytic,cmf}` and `hartle_reference.hpp` | validation | 4A |
| **4C-G** | second-order governance: a **separate ADR** fixing the variable (`p₀*` vs `δp`), (H1)–(H2) with `j²`, the fixed-`ε_c` boundary condition, the centre start, `δM = m₀(R) + J²/R³`, EOS `dε/dp` authority, exposure as `Ω²`-normalized coefficients; expected to invoke `GOVERNANCE.md` §3.1 because the candidate's output cannot serve as a baseline (§12) | governance | may be drafted in parallel with 4A; accepted before 4C-I |
| **4C-I** | implement the corrected l = 0 system (replace, do not patch, the candidate) | scientific-semantic | 4C-G accepted; 4A for the physical-units layer (the `Ω²`-normalized coefficients themselves are seed-free and can be computed from the raw solution, exactly as `I` is) |
| **4D** | independent O(Ω²) validation: regular-centre series, Newtonian limit, an independent test-side solver in a different formulation, quadratic scaling and seed independence at a predeclared bound, grid convergence, comparison with published slow-rotation results (constant-density: Chandrasekhar & Miller 1974; polytropes / realistic EOS: Hartle & Thorne 1968 and later tabulations) | validation | 4C-I |
| **4E** | fixed-`ε_c` structural response products for Phase 5: `m₀/Ω²`, `δp₀/Ω²` (or `p₀*/Ω²`), `ξ₀/Ω²`, `δM/Ω²`, `ω̄/Ω`, `I`, at fixed central density; `RotochemicalCache` consumes these, never raw profiles | structural | 4D |
| — | `MixedStar` rotation | separate track | `MixedStar` coverage (ADR-0004 §0-Q2) |

Roadmap items re-read against this plan: "replace `init_omega_bar = 5e-3` and scale to a physical
Ω; correct the `[s^-1]` annotation" → 4A; "verify true RHS-radius interpolation and cached
bracket indices" → largely discharged by 2B-4B (production agrees with the independent solver to
`2e-5`, resolution-independent) — a targeted check remains in 4B; "verify the second-order
equations against a cited derivation" → 4C/4D; "make `HartleResult` reachable from `NStar`" →
re-scoped to 4E (expose normalized response, not the raw struct); "convergence tests and
validation against published `I` and `δM`" → 4B/4D. **Phase 5 must not begin before 4E.**

## 16. What Phase 5 genuinely needs from Phase 4

From `RotochemicalCache` (`RotochemicalCache.hpp:13-25`, `.cpp:47-104`, `:107-177`) and the
`Rotochemical` driver (`Rotochemical.cpp:22-25`, `:62-122`), both NOT COMPILED, and from
Fernández & Reisenegger (2005):

| Needed | Source | Normalization requirement |
|---|---|---|
| physical `Ω(t)`, `Ω̇(t)` | `SpinState::Omega()` [rad/s] and the Spin RHS — **not** Hartle | already physical |
| `I` (for the spin-down torque / energy budget; `MagneticDipole::use_moment_of_inertia` is today a logged no-op) | first order | scale-free; physical `I_cgs` needs the `G` authority |
| `A_i = (∂N_i/∂Ω²)|_{ε_c}` — needs the **monopole** Eulerian pressure perturbation **per unit `Ω²`**, `δp₀/Ω²` (dimensionless in geometric units), plus the O(Ω²) metric/volume corrections that involve `m₀/Ω²` and `ω̄²/Ω²` | second order, **fixed-`ε_c` family** (`p₀*(0) = 0`) | seed-free `Ω²`-response coefficients; the l = 2 sector integrates out of `N_i` (`∫P₂ dΩ = 0`), so only l = 0 is required |
| `B_i = (∂N_i/∂ε_c)|_Ω` | TOV sequence (not Hartle) | — |
| `ξ₀/Ω²`, `δM/Ω²` | second order | diagnostics / cross-checks of the same solution |

Therefore the Phase-4 output contract should be **rotation response coefficients** —
`ω̄(r)/Ω`, `I`, and (after 4C/4D) `m₀(r)/Ω²`, `δp₀(r)/Ω²`, `ξ₀(r)/Ω²`, `δM/Ω²` at fixed `ε_c`,
in geometric units with `Ω_geom = Ω_phys/c` — rather than any seed-dependent profile. Phase 5
needs derivatives with respect to `Ω²`, not a solution at one arbitrary amplitude (INV-09: `A_i`
"never divided by Ω²" is exactly the symptom of exposing the wrong product).

## 17. Owner questions (real decisions only)

| | Question | Recommendation |
|---|---|---|
| **Q1** | Public spin input: `Ω [s⁻¹]`, `f [Hz]`, `Ω_geom [km⁻¹]`, or a typed quantity? | **`Ω` in rad s⁻¹ (= s⁻¹)**, the quantity `SpinState::Omega()` already evolves (`SpinState.hpp:65`), passed through an explicitly named parameter/type with `FromHz`/`FromPeriod` helpers; geometric `km⁻¹` strictly internal via one conversion `Ω/c` |
| **Q2** | Does the arbitrary seed remain strictly internal? | **Yes.** It is a numerical normalization (§8); it must not appear in `HartleResult`, exports or accessors as a physical quantity; its magnitude is governed only through the integrator-tolerance requirement (§12) |
| **Q3** | `HartleResult`: store both `Ω_geom` and `Ω_phys`, or one canonical quantity + conversion accessors? | **One canonical geometric representation** (consistent with `J [km²]`, `I [km³]`, `ω̄ [km⁻¹]`), annotation corrected, with explicitly named conversion accessors (`c` only for Ω/ω̄; cgs `J`, `I` deferred to the `G`/solar-mass ADR). No duplicate state |
| **Q4** | How is first-order (and later second-order) response exposed to Phase 5? | **Explicit normalized response fields** (`ω̄/Ω`, `I`, later the `Ω²` coefficients at fixed `ε_c`) — never raw seed-dependent profiles and never `NStar`'s private struct |

## 18. Validation and detector plan for ADR-0006 (plan only; tolerances predeclared now)

Bounds are chosen **before** implementation from the measurements in §8.3 and §12 and from
`HARTLE_MOMENT_INERTIA.md`; they are not to be revised after the fact.

| # | Requirement | Predeclared bound | Basis |
|---|---|---|---|
| 1 | changing the internal seed over `[10⁻³, 10³]` leaves `ω̄_phys(R)`, `J_phys`, `Ω` unchanged | rel `≤ 1e-10` | exact arithmetic + production `abstol = 1e-10`; the reference solver gives `6e-15`. A failure is an integrator-tolerance finding, not a reason to widen |
| 2 | requested `Ω_phys` is recovered at the surface: `c·(ω̄_phys(R) + Rω̄'_phys(R)/3)` | rel `≤ 1e-13` | pure arithmetic |
| 3 | `J_phys = I·Ω_phys/c` | rel `≤ 1e-13` | measured `2e-16` |
| 4 | `I` bit-identical (golden `hartle_I_dscmf1_debug.tsv` hash unchanged; both Hartle CTests unchanged) | **bitwise** | ADR-0006 authorizes no change to the ODE, the extraction or the seed value; any seed change needs its own predeclared bound |
| 5 | `ω̄_phys(r)` linear in requested `Ω`: `ω̄_phys(r; 2Ω) = 2 ω̄_phys(r; Ω)` node by node | rel `≤ 1e-14` | arithmetic |
| 6 | conversion checked against an **independent literal** `c = 299792.458 km/s` (exact by definition): `Ω_geom·c = Ω_phys` | rel `≤ 1e-15` | the test must not import `LIGHT_C_KM_S` as its oracle |
| 7 | every unit annotation and exported header token matches the stored quantity (parse the export; assert on `HartleResult` accessor units) | exact | §9 |
| 8 | zero spin: `Ω_phys = 0 ⇒ ω̄_phys ≡ 0`, `J_phys = 0`, `I` finite and equal to the raw `I` | exact | §8.1 item 5 |
| 9 | slow-rotation regime reported: `Ω/Ω_K` with `Ω_K = (M/R³)^{1/2}` geometric; validation set `1.0–2.0 M☉` at `Ω_phys ∈ {100, 2π·100, 2π·716} s⁻¹` all below a documented threshold | diagnostic, threshold to be governed | §8.1 item 6 |

Detectors (each applied, measured, reverted, byte-identical tree; expected failing test in
brackets):

| | Mutation | Expected failure |
|---|---|---|
| **D1** | omit the rescaling factor `A` (return raw `ω̄`) | requested-Ω surface test [2], `J = IΩ` [3] |
| **D2** | multiply by `c` instead of dividing on physical → geometric | independent-constant test [6] fails by `c² ≈ 9e10`; [2], [3] |
| **D3** | rescale `ω̄` but not `dω̄/dr` / `J` | [3] and [2] (surface matching inconsistent) |
| **D4** | let `init_omega_bar` leak into the normalized output (e.g. return `ω̄_raw` for `ω̄_phys`) | multi-seed invariance [1] |

No second-order detector is defined until the second-order equations are governed (§15, 4C-G).

## 19. Exact non-scope of this record

Not done, deliberately: no `RotationSolver` production change; no repair of any O(Ω²) equation;
no rotochemical work; no BNV / Decay / DarkCore change; no `MixedStar` rotation validation; no
ADR acceptance; no baseline of any candidate output; no unit-label edit in source; no change to
`tests/`, `tests/baselines/`, `CompactStar/**` or any `CMakeLists.txt`; no adjudication of the
solar-mass/`G` authority; no relocation of the Phase-2B convergence scope (`PHASE2B_CLOSURE.md`
§6 still awaits owner adjudication).

## 20. Documents synchronized in this task (documentation class)

| File | Change |
|---|---|
| `docs/validation/PHASE4_ROTATION_ENTRY.md` | **new** — this record |
| `docs/adr/ADR-0006-hartle-first-order-physical-normalization.md` | **new — PROPOSED** |
| `docs/SCIENTIFIC_INVARIANTS.md` | INV-05 rotation half re-audited (⚠ discharged); INV-07 current line references, re-audit note, ADR-0006 PROPOSED (status unchanged: UNRESOLVED); INV-08 current-behaviour portion corrected — reachability, variable definition, boundary condition, `dε/dp` parenthetical, dimensional findings, line references; historical `d91c31b` findings retained and labelled; INV-09 cross-reference; summary table |
| `docs/architecture/CURRENT_ARCHITECTURE.md` | evidence-scope row; new label; `RotationSolver` rows (O(Ω) note; O(Ω²) reclassified; declared-undefined API row); §6 non-claim corrected |
| `docs/MODERNIZATION_ROADMAP.md` | Phase-3 merge status (merged at `df859b5`); Phase-4 entry status and sub-phase plan; dependency summary |
| `docs/adr/README.md` | ADR-0006 index row; anticipated second-order ADR; MixedStar row extended to rotation |
| `docs/architecture/TARGET_ARCHITECTURE.md` | three status cells in §2 (non-normative; minimal) |

No accepted ADR decision text was altered.

## 21. Status

> **`PHASE-4 ROTATION ENTRY AUDIT COMPLETE — ADR-0006 OWNER ADJUDICATION REQUIRED`**

No fundamental contradiction blocks the normalization proposal: the first-order equation is
linear-homogeneous with a one-dimensional regular solution space, the exterior matching is
linear, and the contract of §8 follows uniquely. The second-order candidate is a separate matter
(§10–§12) and is not certified, touched, or baselined here.

**Exact next action.** The owner reviews and adjudicates ADR-0006 (Q1–Q4) **before any
`RotationSolver` production change**. Neither the normalization implementation (4A) nor any
O(Ω²) repair (4C) begins before that adjudication.
