# CompactStar Architecture — Index

Architecture documentation is split in two, because mixing observed behavior with intended
design is how a codebase comes to describe scaffolding as if it were operational.

| Document | Question it answers | Authority |
|---|---|---|
| [`docs/architecture/CURRENT_ARCHITECTURE.md`](docs/architecture/CURRENT_ARCHITECTURE.md) | What does the code **do** today? | Descriptive — authoritative for component boundaries and ownership (`GOVERNANCE.md` rank 6) |
| [`docs/architecture/TARGET_ARCHITECTURE.md`](docs/architecture/TARGET_ARCHITECTURE.md) | Where is the design **going**? | Non-normative (rank 8) |

**On disagreement about the present, `CURRENT_ARCHITECTURE.md` wins.**

## Related governance

| Document | Purpose |
|---|---|
| [`GOVERNANCE.md`](GOVERNANCE.md) | Authority hierarchy, change classes, fail-closed conditions |
| [`docs/SCIENTIFIC_INVARIANTS.md`](docs/SCIENTIFIC_INVARIANTS.md) | Physical and numerical conventions, with status labels |
| [`docs/adr/`](docs/adr/) | Decision records |
| [`docs/MODERNIZATION_ROADMAP.md`](docs/MODERNIZATION_ROADMAP.md) | Dependency-ordered plan |
| [`docs/reconnaissance/`](docs/reconnaissance/) | Point-in-time audit snapshots (non-normative) |
| [`AGENTS.md`](AGENTS.md) | Operating contract for AI coding agents |

---

### Note on the previous root `ARCHITECTURE.md`

An earlier root-level `ARCHITECTURE.md` was produced during Phase-0 reconnaissance and added to
`.gitignore` on branch `claude/hartle-rotochemical-heating-FwAHP` (commit `d91c31b`), so it was
never tracked. That `.gitignore` entry does **not** exist on `master`, and the untracked local
file remains untouched in its own worktree.

Its content has not been discarded: the observed-behavior half was authenticated into
`CURRENT_ARCHITECTURE.md`, and the intended-design half into `TARGET_ARCHITECTURE.md` with
explicit status labels. Several of its claims did not survive authentication — see
`docs/reconnaissance/2026-08-31-phase-0-reconnaissance.md` §B.
