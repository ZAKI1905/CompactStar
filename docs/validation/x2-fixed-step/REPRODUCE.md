# X2 reproduction and office-Mac handoff

Candidate branch: `candidate/linux-x2-fixed-step-20261016`; prerequisite commit:
`6cc7e20713e3c0210faa9036316650816f0290d5`. No merge or promotion is authorized.
The initial declaration is commit 75d8c02; primary execution source is d313cc3.
Subsequent commits preserve diagnostics, corrections, evidence and the report.

## Inspect and authenticate without solving

1. In this checkout, verify `docs/validation/x2-fixed-step/SHA256SUMS` with
   `sha256sum -c SHA256SUMS` from that directory. Verify the unchanged X1 sums
   from `docs/validation/x1-oracle-anchor` as well.
2. `evidence/external-manifest.json` gives absolute artifact root, relative names,
   byte counts and SHA-256. Every execution gate binds executable, source and input
   hashes. Diagnostic generated sources live outside Git but their generators and
   exact source hashes are committed. Original X1 sources were never changed.
3. `evidence/raw-decimal-text.json` preserves exact numerical strings, including
   long-double output precision. `scores.json` is a readable binary64 summary;
   use the raw decimal strings for a high-precision reanalysis.
4. `tools/linux/x2_oracle/analyze.py` produces the compact evidence from the raw
   external artifacts. `audit.py` independently uses 80-digit Decimal arithmetic,
   verifies all hashes, checks bit-identical X1 replays, convergence, support and
   headrooms. The preserved audit corresponds to commit fc8dc14 and job 17034.
   A new audit is a new receipt, not a replacement of that evidence.

## Numerical reproduction on Slurm

The original artifact directory is immutable evidence. Runners refuse to overwrite
existing numerical stdout files. Make a separate scratch cluster root with layout
`<root>/worktrees/x2-fixed-step`, `<root>/toolchain` and `<root>/x2-fixed-step`.
Check out this branch into its new worktree location and point toolchain at the
same authenticated Linux compiler, Zaki headers/archive and Python. The repo
includes X1 inputs and headers; no production EOS table is needed.

Build with `bash tools/linux/x2_oracle/build.sh`. This produces extended and
binary64 executables outside Git. Direct command syntax is:

    extended B Hc_guess surfaceP N center_factor eos_order raw plain fixed_H

For example the prospective target-0 pilot command was:

    extended 7.6169065187188748e56 0.13154890784587198 4.9156528789505745e19 32 1e-4 96 0 0 0

`extended --check` performs the mathematical tests. Do not run substantial stellar
solves on the login node. Use one pinned CPU under Slurm. `run.sbatch` sets the
thread environment and `run.py` pins affinity. Its absolute paths must point to
YOUR scratch root for a new run; preserve the old launcher bytes in this branch.
Freeze the new source/executable/input hashes and new launch paths as a separate
reproduction gate before running. Recompute each `files` hash in the gate using
SHA-256 rather than claiming the original absolute paths or binary identity.
The compiler flags are C++17, -O3, -ffp-contract=off. The original cluster Real
has 64 mantissa bits; do not assume an office Mac's `long double` has that format.

The exact execution schedule was:

    run.py pilot
    run.py ladder
    run.py 2048
    run.py 4096
    run.py ablations

The last command stopped at the declared binary64 stage gate. Preserve that failure.
`x1_diagnostics.py` builds separate instrumented old-source copies;
`run_x1.py` runs their frozen replays. `continue_diagnostics.py build` generates the
separate early-Newton-stop correction and its gate; its no-argument mode resumes
unaffected ablations and records the remaining binary64 Newton failure. Finally,
`shell_diagnostics.py build` generates the post-trace mechanism copies and gate;
its no-argument mode executes the declared historical-tolerance tests. Follow the
chronological predeclarations; do not silently replace the failed outputs.

Submit the auxiliary Python runners using the same Slurm partition/CPU/memory
settings and a PATH containing `<root>/toolchain/bin`. Compute nodes do NOT provide
Git through their default PATH or `/usr/bin/time`; Python getrusage is used for
receipts. The final successful audit launcher is retained in its Slurm allocation
and source record. Original launch commands are also in controller receipts.
All original numerical commands, arguments, environment changes, CPU costs, return
codes and stdout/stderr hashes are in `evidence/resources.json`. No extra tolerance,
resolution, target, certificate run or BNV trajectory is required for reproduction.

## Transfer to the office Mac

GitHub CLI and noninteractive GitHub credentials were unavailable. External handoff:

    /mnt/sdd/zaki/CompactStar/linux-transition-20261010/x2-fixed-step-candidate.bundle
    /mnt/sdd/zaki/CompactStar/linux-transition-20261010/x2-fixed-step-candidate.bundle.sha256
    /mnt/sdd/zaki/CompactStar/linux-transition-20261010/x2-fixed-step-handoff.json
    /mnt/sdd/zaki/CompactStar/linux-transition-20261010/x2-fixed-step-bundle-verification.txt

The bundle is incremental against 6cc7e20. Copy these four files to the Mac.
In a repository already containing that commit, run `shasum -a 256` on the bundle
and compare the sidecar, then `git bundle verify /path/to/x2-fixed-step-candidate.bundle`.
Fetch the candidate into a review ref without changing the active checkout:

    git fetch /path/to/x2-fixed-step-candidate.bundle refs/heads/candidate/linux-x2-fixed-step-20261016:refs/remotes/x2-handoff/candidate
    git worktree add -b review/x2-fixed-step-20261016 ../CompactStar-x2-review refs/remotes/x2-handoff/candidate

The handoff JSON records final commit and bundle hash. Verification includes bundle
prerequisite validation and importing its candidate tip into an isolated bare
repository with access to the prerequisite object store, then checking committed
X2 bytes against the source checkout. This is transfer verification, not scientific
ratification. Raw external numerical artifacts can be copied separately according
to the committed source-bound manifest; no executable is embedded in the bundle.
