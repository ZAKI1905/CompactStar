# Reproduction and transfer

The numerical run is complete. These commands document reproduction; they do
not authorize another experiment. The retained runner refuses to overwrite
`runs.json`. No further stellar solve or tolerance was run for this report.

The cluster root is `/mnt/sdd/zaki/CompactStar/linux-transition-20261010`.
The candidate worktree is `worktrees/x1-oracle-anchor`; external evidence and
binaries are under `x1-oracle-anchor`. The exact prediction and execution gates
are committed as described in [REPORT.md](REPORT.md). The final evidence contains
all outputs needed to inspect or independently recompute the scores without
loading any large stellar profiles.

Build recipe: `bash tools/linux/x1_oracle/build.sh`. It uses the authenticated
GCC 14.1.0 wrapper, GSL 2.7.1 archives and Zaki constants; `-O3 -DNDEBUG
-std=c++17 -ffp-contract=off`. No production CompactStar library is linked.
Build CPU receipts are retained separately from numerical work. The binaries
and logs are hashed in `evidence/external-manifest.json`, and the executed
sources, dependencies and frozen references in `execution-gate.json`.

Numerical entry point was `sbatch tools/linux/x1_oracle/run.sbatch`, after
committing the source/input hash gate. It ran precisely the 18 prescribed
solves, the fixed mathematical checks, and `score.py` within one allocation.
The score script accepts the external result directory as its argument and
reads only frozen compact reference records. Its A2 decisions respect the
global gate even when individual residuals look favorable. A separately
committed `audit.py` checks all classifications with 60-digit Decimal arithmetic
and authenticates every raw output; it adds no stellar solve or tolerance.

GitHub dry-run push failed with `fatal: unable to get password from user`.
A verified incremental bundle is supplied under the cluster root:
`x1-oracle-anchor-candidate.bundle`, with SHA256 and verification receipts.
It requires the existing starting commit
`dad05b17056f3ccb1e3871e6e4bb24ebe5db5df4` on the office Mac.

After copying the bundle and its checksum to the Mac, from the repository:

```sh
shasum -a 256 -c x1-oracle-anchor-candidate.bundle.sha256
git cat-file -e dad05b17056f3ccb1e3871e6e4bb24ebe5db5df4^{commit}
git bundle verify /path/to/x1-oracle-anchor-candidate.bundle
git fetch /path/to/x1-oracle-anchor-candidate.bundle \
  refs/heads/candidate/linux-x1-oracle-anchor-20261016:refs/remotes/eku/x1-oracle-anchor
git worktree add -b review/x1-oracle-anchor /path/to/new-review-worktree \
  refs/remotes/eku/x1-oracle-anchor
```

This imports the candidate for review. It does not merge master, accept a
certificate, modify production owners, or authorize a continuation. Absolute
cluster paths in the provenance are historical identities; do not rewrite
them merely to make a Mac path look local. Compact output hashes can be
verified directly in the imported worktree. Full binary/source-gate checks
require the authenticated cluster toolchain and its recorded paths.
