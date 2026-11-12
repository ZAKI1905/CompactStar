# Stage B documentation closeout and next-study handoff — 2026-11-12

**Documentation response complete; owner ratification pending.** This record
covers the owner's request to address all four corrections in the complete
independent review, including §9, and formulate a separate physics study.
Authority: `GOVERNANCE.md:42` (documentation/provenance), `GOVERNANCE.md:147`
(owner ratification), ADR-0015 §§3–11 (energy and composition), INV-15
(`docs/SCIENTIFIC_INVARIANTS.md:1256`). No implemented scientific-semantic or
numerical-method change is part of this submission.

## Authenticated source and exact evidence

Repository: `/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot`.
Branch: `candidate/results-first-bnv-pilot-20261108`.
Reviewed and entry HEAD: **`49253200a801ad9b336789dffb4e93d535a8e903`**.
Accepted Stage A ancestor: **`c90b2205e553afcd102c7acd5f84b253e0a67f1b`**.
All 15 local branch tips were ancestors of entry HEAD; nine worktrees existed.
No competing relevant branch changes or tracked worktree changes existed.
Three pre-existing untracked Python bytecode directories were left untouched.
Evidence: `review-evidence/IMPORT.json:1`; the reviewed preflight worktree list
is in `evidence/preflight.json:1`. This is the intended candidate source, not
an authorization to merge into master.

The original review root is
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/stage-b-independent-review-20261111`.
Its `REVIEW.md`, nine numerical/provenance check outputs plus `longterm.json`
are copied byte for byte under [review-evidence](review-evidence/IMPORT.json).
All 11 imported files were rehashed against the original `SHA256SUMS`; original
paths, full digests and that manifest's digest are in `review-evidence/IMPORT.json:1`.
The original review, its clean export and all sealed campaign roots are unmodified.
This selected-file check is not a new full campaign-manifest verification.

Campaign roots remain `/mnt/sdd/zaki/CompactStar/linux-transition-20261010/`:
`stage-b-20261110-r1`, `stage-b-20261110-r2`, `stage-b-20261110-r3`.
Full manifest: `stage-b-finalization-20261110-r1/sealed/manifest.jsonl.gz`;
SHA-256 `4eae8a0132354450457f6db79460f30f768f1d9af2123def691a331688af25e0`.
Original compact records remain under `evidence/`, unchanged. The independent
review's full rehash and 58/58 assessment rerun are historical retained evidence,
not work repeated in this closeout (`review-evidence/REVIEW.md:96`, `:134`).
For the exact delivery commit use `git log -1 --format=%H -- docs/validation/results-first-bnv-pilot/stage-b/CLOSEOUT.md`
in this checkout; the handoff message supplies its full hash. A commit cannot
embed its own hash without changing that hash.

## Changes and correction map

| File | Change |
|---|---|
| REPORT.md | D1 full chemical-channel table and diagnostic provenance; D3 conservative R20 decomposition and Stage A/B totals; D4 three distinct uncertainty classes; current independent-review status and exact acceptance limits |
| PIECEWISE.md | D2 native knot offset (15 ulp at knot 94), old-cell RKF45 stage overshoot, and distinction between reference audits and main-history evidence |
| CHEMICAL_RESPONSE.md | Distinguish original analytic comparison from subsequent reviewer term-suppressed integrations; do not extend the reaction enclosure to storage/total |
| CONDITIONAL_ACCEPTANCE.md | Draft disposition with C1–C7, exact numerical scope and review limits; no owner ratification |
| NEXT_PHYSICS_STUDY.md | Governing equations, dimensionless powers, regime scalings, conservative retention choices and physical limitations |
| docs/adr/ADR-0025-bnv-retained-energy-study.md and docs/adr/README.md | Proposed alternatives and validation obligations; Decision remains empty |
| review-evidence/ | Eleven unchanged source copies and IMPORT.json provenance |
| CLOSEOUT.md and closeout-checks.json | This change record and document/evidence validation receipt |

The four endpoint diagnostics are +5.19, −2.59, −70.5 and −73.1 µK. R20 totals
are 4.793855e−5 / 4.036164e−5 N for Stage B; 9.034e−5 / 6.99e−9 N belong to
Stage A. Original analytic bounds and reviewer numerical estimates retain their
different evidence classes (`review-evidence/REVIEW.md:281`, `:339`, `:397`).
Historical failures, certificates, thresholds and scientific output bytes are unchanged.

## Validation and remaining decisions

`closeout-checks.json:1` records selected-file hash equality, channel-value
transcription, existing 58-check status, local document links, unmodified original
evidence/physics/tests/tools, and git whitespace checks. These are documentation
checks, not a solver qualification. Governing identities and dimensions in the
proposal were checked algebraically against ADR-0015 §§6, 8–10; no predictor
solve or integration was used.

No build, solver test suite, native evolution, new numerical qualification campaign,
cluster submission, full campaign rehash, or historical failed run was repeated.
The task needs only read-only evidence inspection and document checks on EKU;
no compute allocation or CPU limit was needed. Prior tests and review reruns are
reported only with their original provenance.

No source-of-truth conflict remains. The proposal's retention convention is an
explicit owner decision: actual-budget fraction is recommended, while fixed
equilibrium-budget retention and a microscopic process are distinct alternatives.
The finite-T floor at low retained energy is a physical interpretation limit,
not a repair to Stage B. Owner Stage B ratification and approval of any next-study
implementation remain outstanding; this documentation task supplies both reviewable
records without taking those decisions.
