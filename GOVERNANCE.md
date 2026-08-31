# CompactStar Governance

> **STATUS: RATIFIED — 2026-08-31.**
> Ratified by the project owner at commit `617bb0e`, subject to the status distinctions and
> unresolved items recorded in the governing documents. **This document is normative.** See §7
> for the ratification record, including what ratification explicitly does *not* settle.

CompactStar is a scientific instrument. Its outputs are intended to support publishable
physics. Governance exists so that a result can be traced to the assumptions, equations,
numerical methods, and validation that produced it.

---

## 1. Authority hierarchy

When two sources disagree, the **higher-numbered rule loses**.

| # | Authority | Status |
|---|---|---|
| 1 | `GOVERNANCE.md` (this document) | **Normative** — ratified 2026-08-31 |
| 2 | Accepted decision records — `docs/adr/ADR-*.md` with status ACCEPTED | Normative |
| 3 | `docs/SCIENTIFIC_INVARIANTS.md` | **Normative** — ratified 2026-08-31 |
| 4 | `docs/NUMERICAL_POLICY.md` | Normative once ratified *(not yet written)* |
| 5 | `docs/VALIDATION_POLICY.md` | Normative once ratified *(not yet written)* |
| 6 | `docs/architecture/CURRENT_ARCHITECTURE.md` | Descriptive; authoritative for component boundaries and ownership |
| 7 | Implementation source and code comments | Descriptive |
| 8 | Historical or target-design documents — `docs/architecture/TARGET_ARCHITECTURE.md`, `notes.txt`, `README.md`, `docs/reconnaissance/*` | Non-normative |

Two consequences worth stating plainly:

- **Code is not self-justifying.** That the source does something is evidence of current
  behavior (rule 7), not evidence that the behavior is correct or intended.
- **Aspirational documents bind nothing.** `README.md` and `TARGET_ARCHITECTURE.md` describe
  intent. Where they disagree with `CURRENT_ARCHITECTURE.md`, they are wrong about the present.
- **Unwritten authorities bind nothing.** Ranks 4 and 5 are reserved slots. `NUMERICAL_POLICY.md`
  and `VALIDATION_POLICY.md` do not exist, and ratifying this document does **not** make them
  normative or imply they are pending. They will be written when Phase-1/Phase-2 work supplies
  enough evidence to write them usefully — not before. Until then, no requirement anywhere may be
  discharged *only* by citing them, and no work may be blocked *solely* because they are absent.

---

## 2. Change classes

Every change is classified before work begins. The class determines required evidence.

| Class | Definition | Required evidence |
|---|---|---|
| **Scientific-semantic** | Changes physical meaning: equations, conventions, units, state definitions, conserved quantities, reaction channels | Governing authority cited · ADR · validation before and after · provenance record. **Narrow exception:** §3.1 governs the one case where a *prior* baseline cannot legitimately exist |
| **Numerical-method** | Changes numerical results without changing physical meaning: integrator, tolerance, interpolation order, finite-difference scheme, grid | Convergence evidence · regression against baseline · `NUMERICAL_POLICY.md` citation **once that document exists**; until then the numerical rationale is recorded in the change itself |
| **Structural / architecture** | Moves ownership, changes boundaries, promotes or retires a code path | ADR · `CURRENT_ARCHITECTURE.md` updated in the same change |
| **Engineering** | Behavior-preserving: refactor, rename, dead-code removal, warning fixes | Proof of behavior preservation — bit-identical output, or a documented tolerance |
| **Dependency / build** | Toolchain, third-party libraries, build configuration, platform support | Recorded versions · reproducible build instructions · ADR if it changes platform support |
| **Documentation** | Governance and descriptive documents only | Every factual claim carries a `path:line` citation |
| **Generated artifact** | Doxygen output, run results, build directories, profiler traces | Never hand-edited. Governed by an ADR on what may be tracked |

A change that spans classes takes the **strictest** applicable requirement.

---

## 3. Fail-closed conditions

Work **stops and reports** — it does not proceed on a best guess — when any of the following
is true:

1. **Ambiguous units.** The unit of a quantity cannot be established from a normative document
   or unambiguous code.
2. **Ambiguous state meaning.** What a state variable represents — including its redshift frame
   and ordering — is undefined.
3. **Uncertain authoritative code path.** Two or more implementations are live and no document
   establishes which is canonical.
4. **Uncertain cache validity.** It cannot be shown when a cached quantity is invalidated.
5. **Uncertain ownership.** Two components compute or own the same quantity and no document
   says which is authoritative.
6. **Absent validation for a scientific-semantic change.** No baseline exists against which the
   change could be shown correct. **Unless §3.1 applies** — and §3.1 applies only when an ACCEPTED
   ADR says it does.
7. **Source-of-truth disagreement.** Branch, worktree, or checkout states disagree about the
   content of the files being changed.

Encountering a fail-closed condition is a **successful outcome** of a task. The deliverable is
the report and, where appropriate, a PROPOSED decision record — not a guess.

### 3.1 Narrow exception — pre-baseline correctness work

**This is the only exception to condition 6, and this section is its only definition.** Lower-
ranked documents may require compliance with it or implement a workflow around it; none may
restate, extend, or define its own version of it.

There is one situation in which condition 6 would defeat its own purpose. A baseline exists to
detect unintended change. If the behavior about to be captured has already been adjudicated as
scientifically invalid, freezing it makes it the reference against which every later correction
registers as a *regression* — converting a known defect into a durable obligation. Requiring a
baseline first would then not protect the science; it would entrench the error.

A scientific-semantic change may therefore proceed without a prior regression baseline **only
when every one of the following holds:**

1. An **ACCEPTED ADR** identifies a specific current behavior as scientifically invalid or
   internally inconsistent.
2. Capturing that behavior as the project's golden/reference baseline would knowingly enshrine
   the rejected behavior.
3. That ADR **explicitly identifies the minimum correction** that must precede the baseline.
4. **Independent verification exists** that does not depend on agreement with the superseded
   output.
5. The change is **narrowly scoped** to the defect the ADR governs — no adjacent cleanup, no
   unrelated corrections.
6. The change **records which historical outputs are no longer suitable as reference results.**
7. A regression baseline is **created immediately after** the pre-baseline correction is
   validated — the exception defers the baseline, it does not waive it.

The evidence standard is **substituted, not relaxed.** Condition 4 is satisfied by evidence that
would stand even if no prior run of this code existed. Depending on the quantity, that may include:

- dimensional analysis;
- exact or asymptotic analytical limits;
- thermodynamic identities;
- convergence behavior under grid refinement;
- comparison with independently published or reference calculations;
- internal conservation identities;
- cross-check against an independent implementation.

**This is not a general license to change untested scientific code.** The repository is full of
code that is unvalidated, and that alone is never grounds for invoking §3.1. Absent an ACCEPTED
ADR that explicitly authorizes the exception for a named defect, condition 6 remains fully
active and the normal answer stands: **if no baseline exists, creating one is the task.**

An agent invoking §3.1 must state, in its report: which ADR authorizes it; why the current output
cannot serve as a baseline; what independent verification substitutes for regression; and what
baseline is established immediately afterward.

**First and currently only use:** `docs/adr/ADR-0002-thermal-heat-capacity-ownership.md`, which
rejects `PhotonCooling`'s constant heat-capacity denominator and names the minimum correction that
must precede any passive-cooling baseline (roadmap Phase 2A).

---

## 4. Provenance for scientific changes

A scientific-semantic change is incomplete without a durable record of what changed, why, under
which governing equations or authority, what validation was run, and what impact it has on
previously generated outputs. Decision records under `docs/adr/` carry this until a dedicated
provenance policy exists.

---

## 5. AI-authored contributions

**AI-authored code is not scientifically ratified merely because it compiles, passes review, or
has been committed and merged.**

This is not hypothetical. Commit `675b4a9` (2026-04-05) added 964 lines implementing Hartle
O(Ω²) structural perturbations and rotochemical heating coefficients. It was AI-authored,
landed during a gap in active human development, added no tests, modified no build file — so
its two new source files are still not compiled — and its own message describes parts of it as
stubs. It was merged to `master` via PR #1. Subsequent reconnaissance found dropped physics
terms and an unverified derivation carrying `[FIX: confirm exact from textbook]` markers in
shipped source.

Accordingly:

- AI-authored scientific code is an **UNVERIFIED SCIENTIFIC CANDIDATE** until a human with
  domain authority validates it against a cited derivation *and* a numerical baseline.
- Merging does not confer ratification. Neither does the passage of time.
- Candidate status must be visible where the code is described — see
  `docs/architecture/CURRENT_ARCHITECTURE.md`.
- Governance work performed on a branch containing candidate code — including the branch that
  introduced this document — **does not ratify that code.**

---

## 6. Scope of this document

Deliberately excluded for now, to keep governance maintainable: a separate contribution policy
(folded into this document), a separate provenance policy (carried by ADRs), and a release
policy (premature — see `docs/MODERNIZATION_ROADMAP.md` Phase 1).

Governance grows only when a real recurring decision demands it.

---

## 7. Ratification record

**Ratified by the project owner on 2026-08-31**, at commit
`617bb0e78ea2da22c194a5de32126abb628fc50a` — the exact tree reviewed. Ratification is of the
package as a set, subject to the status distinctions and unresolved items recorded within the
documents themselves.

### Ratified

| Document | Rank | Effect |
|---|---|---|
| `GOVERNANCE.md` | 1 | Authority hierarchy, change classes, fail-closed conditions, and the narrowly governed pre-baseline correctness exception (§3.1) are normative |
| `docs/SCIENTIFIC_INVARIANTS.md` | 3 | The invariant register, its status vocabulary, and the accepted ADR-backed contracts are the project's governing record |
| `docs/MODERNIZATION_ROADMAP.md` | — | The modernization phase ordering and prerequisites are accepted |
| `docs/ai/AI_SKILL_PLAN.md` | — | The AI-skill contracts and sequencing are accepted as the basis for creating repository-specific AI skills |

`docs/adr/ADR-0001` and `docs/adr/ADR-0002` remain **independently normative at rank 2**; their
authority derives from their own ACCEPTED status, not from this ratification.

### Explicitly not ratified

This decision does **not**:

- resolve **INV-07** (Hartle normalization) or **INV-11** (chemical-imbalance convention);
- resolve the unresolved sub-items of **INV-06** or **INV-16**;
- validate any `INTENDED BUT UNVERIFIED` item;
- ratify the Hartle O(Ω²) or rotochemical candidate code from commit `675b4a9` (§5 continues to
  govern it);
- certify the numerical correctness of `StarContext::HeatCapacityStar_Tinf`;
- excuse known implementation nonconformance — `PhotonCooling` (INV-15) and `RotochemicalCache`
  (INV-01) must be corrected, not tolerated;
- make the nonexistent `NUMERICAL_POLICY.md` or `VALIDATION_POLICY.md` normative (§1);
- authorize scientific source changes outside the governed roadmap and change-control
  requirements.

**Ratification settles what the conventions are. It does not certify that the code meets them,
and it does not certify that any number the code produces is correct.**

### Consequence

The owner-review criterion of roadmap **Phase 0.5 is satisfied; Phase 0.5 is complete.** Phase 1
— reproducible macOS build and minimal validation plumbing — becomes the active phase.

### Amending a ratified document

A ratified document is changed the way any governed artifact is: by a classified change citing
its authority. A change to §1, §2, §3, or §3.1 is **structural** and requires an ADR, because it
alters what evidence the project demands. Corrections of fact in the invariant register are
**documentation** class. Re-ratification is not required for a documentation-class correction; it
is required for any change that would alter what a previously ratified rule permits or forbids.
