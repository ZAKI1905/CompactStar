# Approved 24-worker implementation

The owner instructed: “proceed with the 24 workers implementation”, after the
proposal to preserve each sequential integrator while running independent
checkpoint calculations in isolated processes, with deterministic ordered
aggregation and coordinated stopping on failure. The preceding instruction
also authorizes longer runs. This preparation implements the proposed 24-worker
harness and tests it against retained authenticated source outputs.

The starting source of truth is `88da80156ae2be37bd8ad57a50034f133eaff8ab` on
`candidate/results-first-bnv-pilot-20261108`, at
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot`.
Preflight establishes a clean entry, all local branch ancestry, nine worktrees,
and authentication of the previous preparation and timed-out source receipts.

Class: engineering, test harness scheduling, build, documentation and generated
evidence under GOVERNANCE §2. ADR-0017 §§4–6 and Amendment 1 govern immutable
main histories and isolated local checkpoint reconstruction. There is no new
production owner or integration method; physical inputs, numerical settings,
point lists, restart selector and all acceptance conditions remain unchanged.
The intended numerical impact is zero, requiring exact regression evidence.

Preparation uses an EKU allocation of 24 physical execution cores, 96 GiB
(4 GiB per worker), and 30 minutes. It builds with four compiler processes,
runs scheduling/transport contracts and reconstructs the retained 240-interval
grid with one and 24 workers. This is a local reconstruction regression on an
existing main history, not a fresh Stage A trajectory or independent global
reference. No new main integration, restart tail or stellar solve is required.
The prepared source wrapper requests 24 workers/96 GiB/two hours, consistent
with the owner's longer-runtime authorization. The complete source matrix is
not launched as part of this implementation/preparation step.

Control execution still requires separate authorization after complete source
qualification. No Stage B, full campaign, BA13 reopening, master merge or
certificate promotion. Previous failures and raw outputs remain read-only.
Snapshot replacement records old identities without rewriting historical receipts.
