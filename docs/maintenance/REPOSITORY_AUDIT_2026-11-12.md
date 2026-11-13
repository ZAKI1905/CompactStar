# CompactStar repository audit — 2026-11-12

**Class: documentation and local Git housekeeping. Scientific impact: none.**
Authority: `GOVERNANCE.md:18` (hierarchy), `:42` (change classes), `:147`
(AI candidates), and `AGENTS.md:12` (authenticated worktree entry).
Measured Git state and verification results below are retained in
[repository-audit-checks-2026-11-12.json](repository-audit-checks-2026-11-12.json).
No production code, test, numerical tool, equation, acceptance criterion,
certificate, baseline or historical failure was changed by this task.

## Entry and preserved worktrees

Cluster workspace: `/mnt/sdd/zaki/CompactStar/linux-transition-20261010`.
Intended checkout: `worktrees/results-first-bnv-pilot`, branch
`candidate/results-first-bnv-pilot-20261108`, entry HEAD
`129c72dc01e21c6a046c185d24987d696e183e13`. All 15 local branch tips were
ancestors of that candidate. No detached HEAD, stale worktree registration,
invalid pointer, ownership mismatch or competing local research revision was found.
All nine worktrees and shared metadata are owned by `zaki:zaki`.

Paths below are relative to the cluster workspace. These are **entry** SHAs;
the candidate advances in this task, while every other tip remains unchanged.

| Worktree | Branch | Entry HEAD |
|---|---|---|
| `sources/CompactStar` | `candidate/linux-ba13-recovery-20261012` | `8e6ac589d8d34e78447ee893e43e5673490eb728` |
| `worktrees/ba13-closure` | `candidate/linux-ba13-closure-20261019` | `093d121135e35bc50ca2a807470906bcb28ef5c3` |
| `worktrees/ba13-evidence` | `candidate/linux-ba13-evidence-20261026` | `84e1062db78b91e4fb0988074941b1aea6019917` |
| `worktrees/ba13-final-decision` | `candidate/linux-ba13-final-decision-20261030` | `bb0f9c5af8ebaf86b5a59aece5551bb9dafbd264` |
| `worktrees/ba13-onset-pilot` | `candidate/linux-ba13-onset-pilot-20261013` | `dad05b17056f3ccb1e3871e6e4bb24ebe5db5df4` |
| `worktrees/ba13-uncertainty` | `candidate/linux-ba13-uncertainty-20261012` | `4f3cf97462ca030f715e26808938f781a2a40945` |
| `worktrees/results-first-bnv-pilot` | `candidate/results-first-bnv-pilot-20261108` | `129c72dc01e21c6a046c185d24987d696e183e13` |
| `worktrees/x1-oracle-anchor` | `candidate/linux-x1-oracle-anchor-20261016` | `6cc7e20713e3c0210faa9036316650816f0290d5` |
| `worktrees/x2-fixed-step` | `candidate/linux-x2-fixed-step-20261016` | `286aafb63002093b3052e23a3972903011f36035` |

Shared Git metadata live in `sources/CompactStar/.git`; the eight linked
worktrees point into its `worktrees/` directory. The Stage A and Stage B
independent-review `clean-tree/` directories are source exports without `.git`,
not detached or unregistered Git worktrees. They remain untouched.

Eight worktrees were clean. The candidate contained three modified research
files and 19 untracked `.pyc` files in three `__pycache__/` directories; there
were no ignored files in the worktree inventory at entry. The three research
files were preserved byte for byte in commit **`58638994845774fc91e9f5cb6feda520186ebefb`**:

- `docs/adr/ADR-0025-bnv-retained-energy-study.md`
- `docs/adr/README.md`
- `docs/validation/results-first-bnv-pilot/stage-b/NEXT_PHYSICS_STUDY.md`

The caches remain on disk and are now ignored by `.gitignore`. No file,
worktree or branch was deleted or moved. No pruning was needed or performed.

## Remote and ancestry findings

`origin` already named the correct existing public repository:
`git@github.com:ZAKI1905/CompactStar.git`. SSH access was verified from EKU.
Network access required the execution environment's network-enabled permission;
no SSH key or credential was read or modified.

GitHub's candidate started at **`191fac300ec3770bc79f505e3afd606075fc65c3`**.
The local entry candidate was **27 commits ahead, zero behind**. Thus the
Stage A/B reporting, Stage B independent-review import and review corrections
were committed locally but had not yet been pushed. In particular, preserve:

- `0e88ec7e5ee981fcc2d93c85cf8a17d53d6d15fd`: Stage B review corrections,
  conditional-acceptance draft and original retained-energy proposal.
- `129c72dc01e21c6a046c185d24987d696e183e13`: seven owner-supplied PDFs,
  reference README, catalog and manifest.
- `58638994845774fc91e9f5cb6feda520186ebefb`: incoming revised competition study and ADR,
  committed unchanged during this audit.

The full intervening history is preserved; no amend, rebase, squash, reset,
force push or master merge is used. Local and remote master both remain
`9a5c1eca40758e4c39f9f80db7324e12ad23e505`.

Fetching 55 remote branch tips identified these historical branches with commits
outside the candidate ancestry. They remain separate; they are not missing
candidate work and are not merged or rewritten by housekeeping.

| Remote branch | Remote SHA at audit | Commits outside candidate |
|---|---|---|
| `analysis/phase6a1-ba12r-numerical-forensics` | `0aab1c2b5b748585326de35499d1830d0c788fae` | 6 |
| `analysis/phase6a1-ba12r-segmentation-diagnostic` | `6057eb92339e5a0596baab6a652c6290d0930658` | 7 |
| `analysis/phase6a1-bnv-recovery-preflight` | `987d78f17906851dbd91342a2008ba65cf02bba1` | 12 |
| `analysis/phase6a1-checkpoint-reconstruction-candidate-recovery` | `ae0d5607fe8889cec70604747dfe3bfea7e5fd97` | 17 |
| `analysis/phase6a1-checkpoint-reconstruction-preflight` | `93e93c7f91a3cd8fced2f7a0961eda9c469c43fe` | 11 |
| `analysis/phase6a1-checkpoint-reconstruction-validation` | `8b783dbe73cc504b8aa00a5e9aeb48677bcf1ead` | 14 |
| `analysis/phase6a1-passive-observation-proof` | `9af8912107ea77a2e2ea51517c1776f26a4b7b49` | 10 |
| `claude/hartle-rotochemical-heating-FwAHP` | `d91c31bfcc9919915ebb0e745946f57c7f31f141` | 1 |
| `physics/adr0018-mac-source-equivalence` | `abf2b628a35e56542c4f8fc42612a4fe7fd70c88` | 3 |
| `physics/phase6a1-bnv-recovery-ba12r` | `971b2bac47c320bd262d8686841e679f703e5ccb` | 5 |
| `physics/phase6a1-bnv-recovery-r0-r3` | `d2822656aae8e264f4985f0d2f7228b90320a577` | 3 |
| `physics/phase6a1-controlled-bnv-implementation` | `e56e6e50040dbcd9dcbee1acecf58843f3dddf1c` | 11 |

None of those branch-only histories changes the research/reference documents
being preserved here. The old Claude branch's sole extra commit, `d91c31b`,
adds `ARCHITECTURE.md` to `.gitignore`; it is unrelated to the new Python-cache
rules. The candidate continues to track its architecture documents. No
source-of-truth conflict requiring a history resolution was found.

## Configuration corrections

- Removed the invalid `omitted` remote pointing to
  `./OMITTED_OBJECTS_NOT_SUPPLIED`. Successful filtered fetching established
  `origin` as the working promisor remote (`blob:none`). A missing historical
  `.gitignore` blob was then successfully retrieved from GitHub.
- The partial clone had 4,938 intentionally omitted objects after fetching
  branch metadata. This is not evidence of missing commits. No bulk blob fetch
  or large-directory materialization was attempted. `git fsck
  --connectivity-only --no-dangling` passed.
- Set same-name upstreams for the candidate, `candidate/linux-transition-20261010`,
  `physics/phase6-bnv-campaign-readiness`, and
  `physics/phase6-breakpoint-aware-evolution`. Other local historical branches
  lack same-name GitHub branches; they were left local rather than assigned an
  unrelated upstream. Master configuration and ref were left unchanged.
- Restored explicit non-cone sparse settings for the main worktree; each
  worktree now owns its sparse settings. Extended only the candidate's existing
  patterns to include research documents, ADRs, references and navigation.
  Added the native/diagnostics test directories and the qualified-background
  header because those tracked files were materialized but omitted from its
  old patterns. All previously materialized tracked files are retained.
- Global and repository `user.name`/`user.email` were unset. No personal identity
  was inferred. Repository `user.useConfigOnly=true` prevents a host-derived
  fallback. These commits use existing `Codex <codex@localhost>` attribution
  through per-command configuration; no agent default is imposed on other
  worktrees. Historical authors remain unchanged.

See [operating instructions](../README.md#git-and-sparse-checkouts) and
[identity policy](../README.md#commit-identity-policy). Local Git configuration
is not replicated by a push; these instructions let other clones reproduce it.

## Organization and literature

The existing layout is retained. Root README now links to the
[repository/research guide](../README.md), which connects production source,
tests, validation, ADRs, references, study proposals and reproduction scripts.
The guide records external artifact locations and authentication procedures.
`AGENTS.md` distinguishes the tracked collection, larger Mac library and
not-yet-imported sources. The reference README's incorrect private-repository
claim is corrected using the owner's explicit public-repository instruction.

All **seven PDFs**, `catalog.tsv`, `README.md` and `SHA256SUMS` were present and
tracked at entry. All PDF bytes match both the committed Git blobs and manifest;
catalog SHA-256 values and sizes agree. No PDF or manifest was changed, and no
catalog entry is missing for an available PDF.

ADR-0014's distinct **R2006** has expected SHA-256
`a286f15e083e52becd95b3000cbb5ec3ed97148681cf10a43f1a1cc5c4d23ae8`
(`docs/adr/ADR-0014-secular-rotochemical-evolution-contract.md:45`). A hash search
of 43 accessible PDFs under the CompactStar workspace, `/home/zaki` and `/tmp`,
plus the supplied reference ZIP, found no match. The Mac library was not mounted;
inaccessible system temporary directories were not inspected. The original
library filename is abbreviated in the historical source table, so no exact
filename or metadata was invented. RFJ2006 is explicitly **not a substitute**.
The [missing-source record](../references/README.md#not-yet-here) remains actionable.

## Scientific evidence checks and limits

The [verification receipt](repository-audit-checks-2026-11-12.json) records:

- Seven PDF hashes, sizes, catalog entries and committed blobs verified.
- Stage A compact manifest: **19/19** files match; original independent-review
  manifest: **26/26** files match.
- Stage B import: **11/11** copies and originals match recorded hashes; the
  original review manifest matches its committed digest.
- Stage B display/source bindings: **30/30** local display copies and external
  source files match their respective hashes; the seal receipt digest matches.
- Stage A external compressed manifest matches
  `4ff01cd3172a3bf92b7c5c2df5b69d5b138c6cd80910eeb53e7d22c6177c3992`;
  Stage B external compressed manifest matches
  `4eae8a0132354450457f6db79460f30f768f1d9af2123def691a331688af25e0`.
- Existing local links in the research/reference handoff resolve against Git's
  tree, including sparse paths. New navigation links are checked before commit.
- The 27 previously unpublished commits introduce 233 blobs (22,500,900 bytes).
  Filename and text credential-pattern scans found no matches. Large external
  payloads are not among those blobs; the largest is the 7,368,656-byte BGZ2022
  reference PDF. Only intended documents and ignore rules are added by this task.

This is **not** a fresh full-payload verification of the 296,997 Stage A and
663,232 Stage B files. The compressed manifest digests and selected source files
were rehashed; full byte passes remain the historical sealing/review evidence.
No calculations, numerical reruns, builds, test suites or cluster jobs were run.

Stage A/B reports, acceptance records, compact numerical evidence, historical
failures, original review exports and sealed roots remain intact. The revised
study makes BNV versus spin-down driving the primary question, with retention
secondary; ADR-0025 remains **PROPOSED** and owner selection is pending
(`docs/adr/ADR-0025-bnv-retained-energy-study.md:1`). Stage B's conditional
acceptance remains a **draft**, not ratification
(`docs/validation/results-first-bnv-pilot/stage-b/CONDITIONAL_ACCEPTANCE.md:1`).

## Publication and remaining owner decisions

Publish only `candidate/results-first-bnv-pilot-20261108` by ordinary fast-forward
push with upstream tracking. This audit file cannot embed its own delivery SHA;
identify it with:

```bash
git log -1 --format=%H -- docs/maintenance/REPOSITORY_AUDIT_2026-11-12.md
git rev-parse HEAD origin/candidate/results-first-bnv-pilot-20261108
git ls-remote origin refs/heads/candidate/results-first-bnv-pilot-20261108
```

The final handoff reports those independently checked local/GitHub tips.
Remaining owner decisions are Stage B ratification, ADR-0025 study alternatives,
which personal name/email to use for owner commits, and supplying the original
missing references if needed. Existing reference authority roles remain unchanged.
No history rewrite or discarded work is needed for this housekeeping task.
