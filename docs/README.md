# Repository and research guide

Keep the existing source and evidence layout. New study proposals belong beside
their scientific context, with consequential choices recorded in `docs/adr/`.
Do not move completed campaigns or rewrite historical status records to make
them look current. Authority and change classes remain those of
[GOVERNANCE.md](../GOVERNANCE.md) §§1–5 (`GOVERNANCE.md:18`, `:42`, `:147`).

## Where to look

```text
CompactStar/                         production C++ headers and implementation
tests/                              tests, governed baselines, bounded fixtures
cmake/, CMakeLists.txt               build configuration
main/                               existing examples and executable entry points
tools/linux/                        cluster reproduction and verification tools
docs/
  architecture/                     current and target architecture records
  adr/                              accepted decisions and proposed alternatives
  references/                       seven tracked PDFs, catalog.tsv, SHA256SUMS
  validation/
    results-first-bnv-pilot/
      preliminary-pair/             Stage A report, acceptance and compact evidence
      stage-b/                      Stage B report, closeout, proposed next study
        evidence/                   compact numerical records and seal receipt
        review-evidence/            imported independent review and provenance
    ...                             earlier campaigns, failures and validation records
  maintenance/                      dated repository audit and verification receipt
```

This is a navigation map, not a new architecture. Source/test ownership is described
in [CURRENT_ARCHITECTURE.md](architecture/CURRENT_ARCHITECTURE.md); scientific
rules are in the [invariant register](SCIENTIFIC_INVARIANTS.md) and
[ADR index](adr/README.md). Build instructions remain in [docs/build](build/EXTERNAL_PACKAGES_MAC.md).

## Current research handoff

| Record | Purpose and status |
|---|---|
| [Stage A report](validation/results-first-bnv-pilot/preliminary-pair/REPORT.md) and [owner acceptance](validation/results-first-bnv-pilot/preliminary-pair/OWNER_ACCEPTANCE.md) | Accepted bounded preliminary pair; retain its exact scope (`OWNER_ACCEPTANCE.md:1`) |
| [Stage B report](validation/results-first-bnv-pilot/stage-b/REPORT.md) | Numerical results and qualification limits (`REPORT.md:1`) |
| [Independent review](validation/results-first-bnv-pilot/stage-b/review-evidence/REVIEW.md) and [import receipt](validation/results-first-bnv-pilot/stage-b/review-evidence/IMPORT.json) | Original review copied byte for byte, with source paths and digests (`IMPORT.json:1`) |
| [Conditional acceptance](validation/results-first-bnv-pilot/stage-b/CONDITIONAL_ACCEPTANCE.md) | Owner-ratified restricted conditional acceptance, 2026-11-12; C1–C7 retained |
| [Closeout](validation/results-first-bnv-pilot/stage-b/CLOSEOUT.md) | Historical review-correction record at `0e88ec7`; its next-study description refers to the earlier proposal (`CLOSEOUT.md:1`) |
| [BNV–rotochemical competition study](validation/results-first-bnv-pilot/stage-b/NEXT_PHYSICS_STUDY.md) | Historical proposal; approved S1/S2 results are linked below |
| [ADR-0025](adr/ADR-0025-bnv-retained-energy-study.md) | ACCEPTED for restricted S1/S2: A, P-A, S-obs, H-new, V-proportional, G-short |
| [Competition results](research/bnv-rotochemical-competition/README.md) | Completed S1 map and S2 short tests; independently cross-checked scientific candidates, awaiting scientific review |
| [Tracked references](references/README.md) | Exact PDFs and provenance; distinguish missing R2006 from RFJ2006 (`references/README.md:1`) |

The owner explicitly ratified restricted Stage B and authorized S1/S2 on
2026-11-12. The new results remain candidates for scientific review; no S3,
certificate extension, microscopic mechanism or retention campaign is authorized.
Historical failures and the original proposal remain available.

## External numerical evidence

The cluster workspace root is `/mnt/sdd/zaki/CompactStar/linux-transition-20261010`.
It is a workspace containing Git worktrees, builds, and external evidence; it is
not itself the repository root. The active candidate checkout is
`worktrees/results-first-bnv-pilot`; shared Git metadata live in
`sources/CompactStar/.git`. The [dated audit](maintenance/REPOSITORY_AUDIT_2026-11-12.md)
lists all nine worktrees and their starting commits.

Large histories, checkpoint packets, native builds and full manifests remain
outside Git. Their authoritative locations are recorded in the committed receipts:

| Evidence | Location relative to cluster workspace | Committed authentication record |
|---|---|---|
| Stage A | `preliminary-pair-20261110-r1/` | [external receipt](validation/results-first-bnv-pilot/preliminary-pair/evidence/external-receipt.json) |
| Stage A manifest | `preliminary-pair-evidence-20261110-r1/external-SHA256SUMS.gz` | Same receipt's `manifest_sha256` |
| Stage B | `stage-b-20261110-r1/`, `stage-b-20261110-r2/`, `stage-b-20261110-r3/` | [manifest receipt](validation/results-first-bnv-pilot/stage-b/evidence/manifest-receipt.json) |
| Stage B manifest | `stage-b-finalization-20261110-r1/sealed/manifest.jsonl.gz` | Same receipt's `manifest_sha256` and ordered directory roots |
| Original independent reviews | `stage-a-preliminary-pair-independent-review-20261110/`, `stage-b-independent-review-20261111/` | Original `SHA256SUMS`; Stage B [import receipt](validation/results-first-bnv-pilot/stage-b/review-evidence/IMPORT.json) |

Authenticate the compressed manifest itself against the committed receipt before
using its entries to authenticate individual external files. A manifest digest
match alone does **not** establish that every payload file still matches. Full
payload verification must read the files and compare their SHA-256 values; the
historical sealing passes and independent review retain that broader evidence.

Use the existing [Stage A reproduction instructions](validation/results-first-bnv-pilot/preliminary-pair/REPRODUCE.md)
and [Stage B reproduction instructions](validation/results-first-bnv-pilot/stage-b/REPRODUCE.md).
Scripts remain in `tools/linux/preliminary_pair/` and `tools/linux/stage_b/`.
Verification outputs must go to a fresh external destination, never into sealed
roots. A repository housekeeping audit does not require a numerical rerun.

## Git and sparse checkouts

The existing remote is `git@github.com:ZAKI1905/CompactStar.git`. Use configured
SSH authentication; do not copy keys or embed credentials. The primary candidate
tracks `origin/candidate/results-first-bnv-pilot-20261108`. Configuration is local
to this clone and is documented in the [audit](maintenance/REPOSITORY_AUDIT_2026-11-12.md).

Before working, inspect `git status --short --branch`, `git worktree list`,
`git branch -vv`, and `git sparse-checkout list`. Check Git's tree before
concluding that a path absent from disk is missing from history:

```bash
git ls-tree -r --name-only HEAD -- docs/references docs/adr docs/validation/results-first-bnv-pilot
```

For the existing EKU **non-cone** sparse checkout, extend the current patterns
without replacing its source/test selections:

```bash
git sparse-checkout add --stdin <<'PATTERNS'
/docs/adr/
/docs/references/
/docs/validation/results-first-bnv-pilot/
/docs/README.md
/docs/maintenance/
/tests/native/
/tests/diagnostics/
/CompactStar/Physics/Rotochemical/QualifiedBackgroundInputs.hpp
PATTERNS
```

Only apply this in the intended worktree after preserving pending changes.
Older worktree branches predate the reference import; sparse rules cannot make
files appear in a commit that never contained them. Do not switch their branches
or materialize large unrelated directories just to read references available in
the candidate. The clone remains partial (`blob:none`); omitted historical blobs
can be fetched through `origin` when needed.

## Commit identity policy

The owner should use their own explicitly chosen name and verified GitHub email
for owner-authored commits. This audit does not infer those values from the
cluster login or historical authors and does not set a global identity.

Agents should use an explicit agent identity per commit, following the existing
`Codex <codex@localhost>` and `Claude <claude@localhost>` convention. For example:

```bash
git -c user.name=Codex -c user.email=codex@localhost commit
```

Do not configure an agent as the shared repository's default author. This clone
sets `user.useConfigOnly=true` to refuse an accidental host-derived identity.
Explicit owner or agent settings are needed for both author and committer;
check `git var GIT_AUTHOR_IDENT` and `git var GIT_COMMITTER_IDENT` with the same
configuration before committing. Preserve historical author/committer records.
SSH access identifies who may push; it does not choose commit authorship or
confer scientific approval (`GOVERNANCE.md:147`).
