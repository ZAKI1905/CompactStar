# AI Operating Contract — CompactStar

Applies to **any** AI coding agent working in this repository — Claude, Codex, Copilot, or
otherwise. Tool-neutral by design.

CompactStar is a scientific instrument. A plausible-looking wrong number is worse than a
refusal, because it can reach a paper.

---

## 1. Authenticate before mutating

Before proposing or making any change, establish and state:

- absolute repository path;
- current branch and HEAD;
- which worktree you are in, and whether others exist;
- working-tree status;
- whether HEAD is the intended source of truth.

**This is not ceremony.** Phase-0 reconnaissance was performed in full against a branch that was
missing the owner's later `RotationSolver` and `MixedStar` work — an entire audit section had to
be marked superseded. Check ancestry against every branch that touches your files before
concluding anything.

## 2. Identify the governing documents

Consult, in authority order (`GOVERNANCE.md` §1): `GOVERNANCE.md` → accepted ADRs →
`docs/SCIENTIFIC_INVARIANTS.md` → numerical/validation policy → `CURRENT_ARCHITECTURE.md`.

Cite the specific invariant or ADR your change relies on. If none governs it, that is itself a
finding.

## 3. Distinguish current behavior from intended design

- `CURRENT_ARCHITECTURE.md` says what the code does.
- `TARGET_ARCHITECTURE.md`, `README.md`, and `notes.txt` say what someone hoped it would do.

Never cite the second group as evidence about the present. When they disagree, report the
disagreement — do not quietly reconcile it.

## 4. Classify the change before starting

Use the classes in `GOVERNANCE.md` §2: scientific-semantic, numerical-method,
structural/architecture, engineering, dependency/build, documentation, generated artifact. State
the class in your report. It determines the evidence you owe.

## 5. Establish validation before scientific refactoring

**Do not refactor numerical code whose current behavior is not captured by a baseline.** If no
baseline exists, creating one is the task — not the refactor. "The tests still pass" is not
available here; there are none yet.

## 6. Prefer generic machinery

Do not hardcode reaction names, campaign names, species names, or one-off execution modes into
generic numerical infrastructure unless it is scientifically unavoidable — and say so if it is.
If several features need the same operation, build one abstraction.

## 7. Preserve scientific behavior unless authorized

Changing a numerical result is a scientific-semantic change, even when it looks like a cleanup.
Removing a clamp, reordering a sum, changing an interpolation, tightening a tolerance — all
change results. Get explicit authorization and produce validation.

## 8. Record ambiguity; never guess

On any fail-closed condition (`GOVERNANCE.md` §3) — ambiguous units, ambiguous state meaning,
uncertain authoritative path, uncertain cache validity, uncertain ownership, missing validation,
or source-of-truth disagreement — **stop and write it down.** Draft a PROPOSED ADR with genuine
alternatives. Do not select among them on the owner's behalf.

Reporting a blocker is a successful outcome.

## 9. Stop on source-of-truth conflicts

If branches, worktrees, or checkouts disagree about the files you are changing, stop immediately
and report the divergence with commit IDs and ancestry. Do not merge, rebase, reset, or check out
over changes to resolve it yourself.

## 10. Report honestly

Every report states:

- files changed (or explicitly: none);
- change class;
- governing authority cited;
- validation run, and its actual result — including failures;
- scientific impact: does this change any number the code produces?
- unresolved ambiguity encountered.

If something was skipped, say it was skipped. If a test failed, show the output. Do not describe
uncompiled or unreachable code as operational.

## 11. AI-authored code is not ratified

Your output is an **unverified scientific candidate** until a human with domain authority
validates it against a cited derivation and a numerical baseline. Compiling is not validation.
Merging is not ratification. See `GOVERNANCE.md` §5 for the specific incident that motivates
this rule.

Mark candidate status where the code is described, and never let a later document quietly
promote it.

---

## Quick preflight

```
□ Path, branch, HEAD, worktree, status recorded
□ Ancestry checked against other branches touching these files
□ Governing invariant or ADR identified and cited
□ Change class stated
□ Baseline exists (or creating it is the task)
□ No fail-closed condition open
□ Scientific impact assessed
```
