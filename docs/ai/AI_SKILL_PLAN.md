# AI Skill Plan — CompactStar

> **STATUS: RATIFIED — 2026-08-31.** Contracts only. **No tool-specific skill files are created
> by this plan**, and none exist yet.
>
> The generation gate is now **open**. Claude skills (`.claude/skills/`), Codex prompts
> (`.codex/`), or equivalents were withheld until `GOVERNANCE.md`,
> `docs/SCIENTIFIC_INVARIANTS.md`, and the accepted ADRs (**ADR-0001**, **ADR-0002**) were
> ratified — generating adapters against draft governance would have baked in conventions the
> owner had not approved, and ADR-0001 in particular could have inverted the meaning of every
> per-species calculation a skill might touch, while ADR-0002 fixes the denominator of every
> thermal-evolution term. That ratification occurred at commit `617bb0e` (`GOVERNANCE.md` §7).
>
> Generating the skill files is therefore now authorized work, governed by the contracts below.
> It is a separate task and has not been performed.

Five skills. The taxonomy is derived from recurring CompactStar workflows, not from directory
structure. Prefer a few strong reusable skills over many bespoke prompts.

All five inherit `AGENTS.md` in full. The contracts below add only what is specific to each.

---

## 1. `repo-authenticate`

- **Purpose.** Establish ground truth before any other work: source of truth, branch ancestry,
  build reality, and which claims in documentation currently hold.
- **Inputs.** Repository path; optionally a subsystem to focus on.
- **Allowed mutations.** **None. Strictly read-only.**
- **Preflight.** None — this *is* the preflight for everything else.
- **Validation.** Every claim carries a `path:line`. Every documentation claim is labeled
  VERIFIED / PARTIALLY VERIFIED / STALE / CONTRADICTED / UNKNOWN.
- **Provenance.** Emits a dated snapshot under `docs/reconnaissance/`, marked
  DESCRIPTIVE / NON-NORMATIVE and scoped to the audited commit.
- **Stop conditions.** Always terminates in a report; never edits. Escalates immediately on
  branch/worktree divergence.
- **Governing documents.** `GOVERNANCE.md` §1 (to label authority correctly).

## 2. `scientific-change`

- **Purpose.** Any change that alters a number the code produces.
- **Inputs.** Change description; governing invariant or ADR; baseline reference.
- **Allowed mutations.** Source, plus an ADR, plus validation artifacts. **Never** governance
  documents in the same change.
- **Preflight.** Cite the governing authority · confirm the build works · confirm a baseline
  exists · confirm no open fail-closed condition covers the area.
- **Validation.** Regression against baseline; convergence evidence if a numerical method
  changed; explicit statement of which previously generated outputs are invalidated.
- **Provenance.** ADR recording what changed, why, under which equations, with what validation.
- **Stop conditions.** Ambiguous units, state meaning, ownership, cache validity, or authoritative
  path → **halt, draft a PROPOSED ADR, do not choose.** Missing baseline → halt; creating the
  baseline becomes the task — **except** under the narrow pre-baseline correctness exception
  defined by `GOVERNANCE.md` §3.1, which requires an ACCEPTED ADR authorizing it for a named
  defect. This skill does **not** define that exception and must not extend it; it implements the
  workflow §3.1 permits, including §3.1's four required report items. First and currently only
  case: **ADR-0002** (`PhotonCooling` heat-capacity conformance, roadmap Phase 2A).
- **Governing documents.** All of them.

## 3. `numerical-audit`

- **Purpose.** Verify a solver against its governing equations and expected order of accuracy.
- **Inputs.** Target solver; cited derivation or reference.
- **Allowed mutations.** **Read-only**, plus a report.
- **Preflight.** Identify the governing equations, conventions, and normalization.
- **Validation.** Term-by-term derivation trace; dimensional check; order-of-accuracy assessment
  against the actual interpolation and integration schemes in use (INV-13 — linear, not spline).
- **Provenance.** Report; a PROPOSED ADR if a defect implies a convention change.
- **Stop conditions.** Cannot locate an authoritative derivation → halt and report. **Never
  "fix" physics inside an audit.**
- **Governing documents.** `SCIENTIFIC_INVARIANTS.md`; `NUMERICAL_POLICY.md` once written — until then, record the numerical rationale in the report (`GOVERNANCE.md` §1).
- **First targets.** Hartle O(Ω) normalization (INV-07); Hartle O(Ω²) source terms and the
  dropped j² factor (INV-08); `A_i` Ω² normalization (INV-09).

## 4. `governance-sync`

- **Purpose.** Keep governance and architecture documents true to the code as it changes.
- **Inputs.** A commit range, or a document to re-verify.
- **Allowed mutations.** **Governance and architecture documents only. Never source, never
  CMake, never tests.**
- **Preflight.** Diff documents against code; re-check every factual claim.
- **Validation.** Every surviving claim carries a `path:line`. Unverifiable claims are marked
  UNKNOWN — **never silently deleted**.
- **Provenance.** Documentation-class commit; no ADR unless a convention actually changed.
- **Stop conditions.** A claim cannot be verified → mark UNKNOWN and report. Code and governance
  genuinely conflict → halt; that is a scientific finding, not a documentation defect.
- **Governing documents.** `GOVERNANCE.md` §1.

## 5. `build-validation`

- **Purpose.** Establish and maintain a reproducible build and the validation baseline the other
  skills depend on.
- **Inputs.** Target platform; dependency versions.
- **Allowed mutations.** CMake, CI configuration, test files, baseline fixtures. **Never
  scientific source. Never numerical constants.**
- **Preflight.** Confirm no build step writes into the source tree; record exact dependency
  versions and platform.
- **Validation.** Clean-clone configure and build succeeds; tests run; results reproducible
  across two runs.
- **Provenance.** Recorded toolchain and dependency versions; ADR for any platform-support change.
- **Stop conditions.** A build step mutates tracked scientific products → halt. A test appears
  destructive → halt. Required external data unavailable → halt and report precisely what is
  missing.
- **Governing documents.** `VALIDATION_POLICY.md` once written — until then, roadmap Phases 1, 2A, and 2B (`GOVERNANCE.md` §1).

---

## Sequencing

```
repo-authenticate     ← available immediately; needs no ratified governance
build-validation      ← Phases 1–2B; unblocks everything else
numerical-audit       ← after governance ratification
governance-sync       ← after the first documents are ratified
scientific-change     ← last; requires a baseline, except the Phase-2A carve-out below
```

`scientific-change` is deliberately last. Until a validation baseline exists, no scientific
change can be shown correct — so a skill that performs one would generally be operating outside
the evidence standard `GOVERNANCE.md` §2 requires.

The sole exception is the one `GOVERNANCE.md` §3.1 defines, invoked in **Phase 2A** under
**ADR-0002**. Its conditions, its evidence standard, and its reporting duties live in §3.1 and
nowhere else.

## Deliberately excluded

One skill per subsystem — EOS-change, TOV-change, Hartle-change, BNV-change, cache-change. Each
would restate `scientific-change` with a different noun, and the maintenance cost would exceed
the benefit. Subsystem specificity belongs in the *invariants* a skill cites, not in the skill.
