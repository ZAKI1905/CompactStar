# Restricted preparation reproduction

Use the existing EKU installation. Source:
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot`.
All outputs, builds, temporary files and failed attempts:
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/results-first-bnv-preparation-20261108`.
Choose a fresh external RF_OUTPUT for every reproduction. Never overwrite the
retained inputs or these receipts. No command here launches a BNV trajectory.

```bash
source /mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain/environment.sh
export RF_SOURCE=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot
export RF_OUTPUT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/results-first-bnv-preparation-reproduction
export RF_COMMIT=$(git -C "$RF_SOURCE" rev-parse HEAD)
mkdir -p "$RF_OUTPUT"
sbatch --export=ALL --output="$RF_OUTPUT/slurm-%j.out" --error="$RF_OUTPUT/slurm-%j.err" \
  "$RF_SOURCE/tools/linux/results_first/cref_completion.sbatch"
sbatch --export=ALL --output="$RF_OUTPUT/native-%j.out" --error="$RF_OUTPUT/native-%j.err" \
  "$RF_SOURCE/tools/linux/results_first/native_binding.sbatch"
"$TC/python/bin/python" -I -B "$RF_SOURCE/tools/linux/results_first/bdot_witness.py" \
  --output "$RF_OUTPUT/bdot-witness.json"
"$TC/python/bin/python" -I -B "$RF_SOURCE/tests/bnv/test_results_first_cref.py"
"$TC/python/bin/python" -I -B "$RF_SOURCE/tests/bnv/test_ba13_final_decision.py"
```

The C-REF assay intentionally exits 2 after writing the failed result. Exit 0 is
not its success criterion: it never issues a certificate. Compare its summary and
historical-row preservation, not the diagnostic sidecar gzip byte stream (which
also records the execution commit). Its CPU cap is 60 seconds. The native job
uses the exact retained dependency lock, Release configuration and
`-ffp-contract=off`. Its mutations target only a fresh copied fixture; canonical
paths and inode inequality enforce isolation from the archive. It constructs an
ordinary thermal context but makes no actual BNV RHS calls or trajectories.

The Bdot witness reproduces the historical literal from the historical B0 and
the native literal from the native B0, with binary64 input identities and 80-digit
arithmetic. It is not a fitted correction or a new microscopic rate.

The complete diagnostic sidecar and build artifacts remain external. Their hashes
are in `preparation-evidence/external-manifest.json`; compact summaries and logs
are committed with `PREPARATION_SHA256SUMS`. Verify from this document's directory
with `sha256sum -c PREPARATION_SHA256SUMS`.

This session runs on EKU Linux; it cannot establish a Mac's mount or writability.
No Mac worktree or internal-SSD staging copy was made. Once the cluster-access
Mac has verified Scriptorium is mounted, writable and has space, use the direct
rsync procedure in REPRODUCE.md with the preparation output path above. The
cluster copies remain authoritative until a verified direct transfer is possible.
