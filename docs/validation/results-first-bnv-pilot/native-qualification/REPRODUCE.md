# Reproduce the startup-refusal evidence

The thermal attempt at source commit
`3436063dec533d0e7be7b9f01b74e95c64fd0da2` was job 19171, not a completed stage.
Do not resubmit it automatically. Stage B is blocked by the Stage A result.
The commands in the first block only authenticate and inspect existing evidence.

```bash
cd /mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot
source /mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain/environment.sh
export RF_SOURCE="$PWD"
export RF_SNAPSHOTS="$RF_SOURCE/docs/validation/results-first-bnv-pilot/native-regeneration/snapshots"
export RF_OUTPUT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/native-p2-qualification-1yr-r1
export RF_EVIDENCE="$RF_SOURCE/docs/validation/results-first-bnv-pilot/native-qualification/evidence"
"$TC/python/bin/python" -B tools/linux/native_pair/snapshot.py check "$RF_SNAPSHOTS"
"$TC/python/bin/python" -B tools/linux/native_pair/test_snapshot.py "$RF_SNAPSHOTS"
sha256sum -c "$RF_EVIDENCE/external-SHA256SUMS"
"$TC/python/bin/python" -B tools/linux/native_qualification/collect_refusal.py "$RF_OUTPUT" "$RF_EVIDENCE"
cat "$RF_EVIDENCE/assessment.json"
```

The first check passes on ekucn1 and refuses on phyn004. The read-only audit can
inventory all files without stopping at the first mismatch. It does not call the
stellar or thermal runner. Full retained login/compute inventories, exact ELF
interpreter, and resolved linkage are in `login-runtime-audit.json` and
`compute-runtime-audit.json` under RF_OUTPUT. `assessment.json` retains compact
summaries. The original read-only audit invocation (job 19172) was:

```bash
"$TC/python/bin/python" -I -B tools/linux/native_qualification/runtime_audit.py \
  "$RF_SNAPSHOTS" > "$RF_OUTPUT/login-runtime-audit.json"
sbatch --parsable --export=ALL --output="$RF_OUTPUT/runtime-audit-%j.log" \
  tools/linux/native_qualification/runtime_audit.sbatch
```

Do not overwrite retained inventories when auditing again: use a fresh output
directory. Audit code and collection assertions were exercised against both
retained inventories; shell syntax checks passed. No numerical tests were rerun
after the startup refusal.

Raw Slurm receipts retain their original trailing spaces and final blank lines;
they are preserved byte-for-byte. Source/documentation whitespace checks exclude
those raw `.txt` receipts. Package and external SHA-256 verification both pass.

For provenance, the **already executed, failed** Stage A command was exactly the
following, after the environment setup above. This is a historical command
record, **not authorization to retry**:

```bash
export RF_YEARS=1
mkdir -p "$RF_OUTPUT"
sbatch --parsable --hint=nomultithread --export=ALL \
  --output="$RF_OUTPUT/slurm-%j.log" tools/linux/native_pair/stage.sbatch
# Returned 19171. FAILED / exit 1 before numerical runner invocation.
```

The fixed run matrix, error norms, resource ceiling, and stop conditions remain
`../native-regeneration/RUN_PLAN.md` and `../native-assembly/NONLINEAR_PLAN.md`.
Any approved repair must retain them and establish new snapshot identities and
compute-node authentication before another bounded thermal attempt. Stage B
requires a passing Stage A and must start independently at the original initial
state. Neither full campaign is authorized.
