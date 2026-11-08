# Reproduce the stop decision and deliver directly to Scriptorium

This package runs **no trajectory**. From the new worktree, create a fresh output
name under the external cluster root (the audit refuses to overwrite an output):

```bash
export TMPDIR=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/results-first-bnv-pilot-20261108/tmp
export PYTHONDONTWRITEBYTECODE=1
/mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain/python/bin/python -I -B \
  tools/linux/results_first/preflight.py \
  --output /mnt/sdd/zaki/CompactStar/linux-transition-20261010/results-first-bnv-pilot-20261108/preflight-recheck.json
```

It verifies the independent review manifest, exact native certificate identities,
branch ancestry and the native-versus-Mac context pins. It hashes the candidate
inputs without constructing a star. New audit time/HEAD/free-space fields may
differ; scientific source remains pinned to bb0f9c5. The review's 12 manifest
entries are authenticated, not silently rerun as a new independent review.

The fresh 12-test static refusal suite is reproducible with:

```bash
/mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain/python/bin/python -I -B \
  tests/bnv/test_ba13_final_decision.py
git diff bb0f9c5 -- CompactStar tests cmake GOVERNANCE.md AGENTS.md
```

The test suite must pass and the protected-source diff must be empty. This does
not execute or qualify a thermal trajectory. Source/control native prerequisites
remain release gates rather than being inferred from this test.

Verify committed compact evidence with `sha256sum -c SHA256SUMS` from this
directory. `evidence/external-SHA256SUMS` records external receipt hashes. Large
native inputs remain under `ba13-recovery/run1`, already bound by the inherited
manifests and the new preflight identities. Original review files remain under
`ba13-final-decision-independent-review-20261030`.

The inherited worktree uses sparse checkout. For tracked documents absent on disk,
read `git show bb0f9c5:path` (notably ADR-0020 and linux-corrections). Their absence
from this sparse working tree is not absent Git evidence. No existing worktree's
sparse rules or files were changed.

On the **Mac**, first verify the actual volume rather than creating a directory
on the internal disk:

```bash
diskutil info /Volumes/Scriptorium
mount | grep ' on /Volumes/Scriptorium '
test -w /Volumes/Scriptorium/Research/CompactStar
df -h /Volumes/Scriptorium/Research/CompactStar
```

Proceed only if those checks show the intended mounted volume and free space.
Create the destination on that volume, then transfer directly from EKU:

```bash
mkdir -p /Volumes/Scriptorium/Research/CompactStar/results-first-bnv-pilot-20261108
rsync -a --checksum \
  zaki@ekucn1.eku.edu:/mnt/sdd/zaki/CompactStar/linux-transition-20261010/results-first-bnv-pilot-20261108/ \
  /Volumes/Scriptorium/Research/CompactStar/results-first-bnv-pilot-20261108/
```

The external directory contains a delivery Git bundle, SHA-256 receipt and compact
checks. Verify SHA-256 after transfer and `git bundle verify` against a repository
containing bb0f9c5. Import the review branch into a Scriptorium-hosted repository
and place any Mac worktree under that same volume. Do not update master. No
internal-SSD staging copy is needed. This transfer has **not** been performed.
