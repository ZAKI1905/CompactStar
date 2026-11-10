# Read-only reproduction

Use the existing source output directory; do not submit parallel_stage.sbatch.
No production configuration changed. The assessment is not a qualification
certificate. Numerical assessment runs on an EKU compute allocation, one core,
4 GiB, at most 30 minutes; keep outputs in a new /mnt/sdd directory.

From the candidate worktree, within that allocation:

```bash
python3 -B tools/linux/native_pair/snapshot.py check \
  docs/validation/results-first-bnv-pilot/parallel-preparation/snapshots
python3 -B tools/linux/native_qualification/assess_source_retry.py \
  /mnt/sdd/zaki/CompactStar/linux-transition-20261010/stage-a-parallel-20261110-r1/source-stage \
  /mnt/sdd/zaki/CompactStar/linux-transition-20261010/source-assessment-NEW
python3 -B tools/linux/native_qualification/assess_completed_source.py \
  /mnt/sdd/zaki/CompactStar/linux-transition-20261010/stage-a-parallel-20261110-r1/source-stage \
  /mnt/sdd/zaki/CompactStar/linux-transition-20261010/source-assessment-NEW
```

The corrected supplemental assessor verifies the restart-comparison row order
against both retained accepted endpoints before doing arithmetic. Every printed
floating record is decoded as binary64 and lifted exactly into 110-digit Decimal.
It reports necessary lower-bound tests and missing evidence, never substitutes
zero for a missing uncertainty, and does not call GSL or any integrator.

To authenticate preserved evidence from the repository root:

```bash
sha256sum -c docs/validation/results-first-bnv-pilot/parallel-source-assessment/SHA256SUMS
gzip -cd docs/validation/results-first-bnv-pilot/parallel-source-assessment/evidence/external-SHA256SUMS.gz | sha256sum -c -
sha256sum -c docs/validation/results-first-bnv-pilot/parallel-source-assessment/evidence/analysis-SHA256SUMS
```

The latter external paths are cluster-local. They remain private filesystem
locations even though their provenance names/hashes are published. Only compact
results/logs/manifests are tracked; raw checkpoint packets, binary IPC and the
full per-point assessment stay under /mnt/sdd. Original failure records remain
unchanged. Control, Stage B and another thermal retry are not part of reproduction.
