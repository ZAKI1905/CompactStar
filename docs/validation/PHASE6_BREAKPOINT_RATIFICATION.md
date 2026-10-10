# ADR-0017 Amendment 1 — owner ratification

Status: **ACCEPTED / HUMAN-RATIFIED — bounded numerical architecture only**.
Date: 2026-10-10. Authority: the human owner's explicit instruction in the
CompactStar Phase-6 BNV Campaign Readiness task of this date. This record consumes
that acceptance; it does not infer acceptance from an AI review or a merge.
Change class: documentation/governance. Scientific implementation bytes unchanged.

## Exact accepted object

Accepted branch: `physics/phase6-breakpoint-aware-evolution`.
Accepted candidate: `0c3f9269c884161c9d89163660a47b141a93ace6`.
Canonical master: `9a5c1eca40758e4c39f9f80db7324e12ad23e505`.
Both local and live origin refs matched these identities at task entry, and both
worktrees were clean. Readiness starts from the accepted commit in fresh branch
`physics/phase6-bnv-campaign-readiness`; master is not modified.

The owner reports independent Fable 5.1 disposition B: PASS with nonblocking
caveats, zero BLOCKING and zero MATERIAL findings. This is attributed to the
owner's instruction; a separate byte-authenticated Fable report was not supplied
with that instruction. The older obs32 review archive is a different review.

## Accepted scope and evidence

Ratification covers breakpoint-aware main RKF45 and local rk8pd, frozen-cell RHS,
independently located events, correct restoration of both GSL caches, passive
scientific observations, event uncertainty propagation, the bounded successor
accuracy assessment, and accepted fail-closed behavior. The exact contract is
`PHASE6A1_PIECEWISE_RECONSTRUCTION_PROPOSAL.md:10`; bounded criteria are in
`PHASE6_BREAKPOINT_VALIDATION_SPEC.md:1`; recorded evidence and limits are in
`PHASE6_BREAKPOINT_QUALIFICATION.md:1`. Their numerical content is unchanged.
This acceptance follows `GOVERNANCE.md:164` (scientific provenance) and
`GOVERNANCE.md:171` (human ratification), retaining ADR-0015/0016 physics and
INV-11's restricted controlled/frozen scope (`docs/SCIENTIFIC_INVARIANTS.md:1046`).

It does **not** establish campaign-global accuracy; replace or retroactively pass
BA12/BA12R; accept a BNV trajectory or numerical baseline; select physical BNV
rates/models; change Phase-5D mathematics; or qualify Linux/cluster execution.
The full 500 kyr campaign remains unauthorized. Later readiness experiments and
solver policies are candidates requiring their own qualification and acceptance.

## Ancestry and integration boundary

Canonical master is an ancestor of the accepted candidate. There are 24 commits
in `master..0c3f9269`, including prior campaign execution, certificate-domain
changes, forensics, imported evidence, and the breakpoint correction. Acceptance
of this amendment does not ratify every earlier Phase-6 commit. Canonical
integration is DEFERRED pending an explicit commit-by-commit authority review and
owner integration instruction. No merge, rebase, or master write occurs here.
The readiness branch preserves that history without promoting it.

Historical qualification documents retain their original dated conclusions.
This record supersedes only their amendment/implementation acceptance status;
their numerical caveats and failures stand.
