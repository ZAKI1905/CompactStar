# Reproduce the bounded uncertainty pilot

This recipe reproduces measurements and refusals, not a qualified certificate.
Use the authenticated installation and the committed predeclaration. No accepted
owner or nominal artifact is edited. The actual study used the default OUT below;
for a repeat set BA13_STUDY_OUT to a **fresh** directory. Scripts now honor this
output-only override; it was added after measurement without changing numerical
formulas. Inputs remain the authenticated `ba13-recovery/run1/raw` tree.

```bash
cd /mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/ba13-uncertainty
source tools/linux/eku_environment.sh   # login node only
export BA13_STUDY_OUT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/ba13-uncertainty-repeat
mkdir "$BA13_STUDY_OUT"                # refuse if it already exists
sbatch --wait --partition=phyall --cpus-per-task=4 --mem=4G --time=00:15:00 \
  --output="$BA13_STUDY_OUT/build-%j.log" tools/linux/ba13_uncertainty/build.sh
sbatch --wait --cpus-per-task=1 --mem=2G --time=00:15:00 --array=0,10,20 \
  --output="$BA13_STUDY_OUT/nominal-%A_%a.log" tools/linux/ba13_uncertainty/job.sh nominal
sbatch --wait --cpus-per-task=1 --mem=2G --time=00:15:00 --array=0-35%7 \
  --output="$BA13_STUDY_OUT/sweep-%A_%a.log" tools/linux/ba13_uncertainty/sweep.sh
"$TC/python/bin/python" -B tools/linux/ba13_uncertainty/predict.py
```

Freeze the predictions in Git before the next block, without changing criteria
from observed results. The completed study's first prediction commit is 1733092.
Compare fresh predictions to the committed reference as a reproducibility check;
a disagreement is evidence requiring investigation, not permission to overwrite.

```bash
sbatch --wait --cpus-per-task=1 --mem=2G --time=00:20:00 --array=0-8%7 \
  --output="$BA13_STUDY_OUT/withheld-%A_%a.log" tools/linux/ba13_uncertainty/holdouts.sh
```

**Expected refusal:** EOS 32768 tasks fail `unordered EOS values` before a star
is published. Retain their `.partial` inputs, `.receipt.json` and logs. Do not
remove nodes or call that holdout PASS. Slurm array failure is an expected
scientific disposition here; inspect every receipt before continuing independent
work. Conditional radial checks can continue even though EOS is inconclusive.

```bash
sbatch --wait --cpus-per-task=1 --mem=2G --time=00:15:00 --array=0-7%6 \
  --output="$BA13_STUDY_OUT/targets-%A_%a.log" tools/linux/ba13_uncertainty/target_holdouts.sh
"$TC/python/bin/python" -B tools/linux/ba13_uncertainty/freeze_target_predictions.py
# Freeze these predictions before the finest target cases.
sbatch --wait --cpus-per-task=1 --mem=2G --time=00:15:00 --array=5,15 \
  --output="$BA13_STUDY_OUT/target-finest-%A_%a.log" tools/linux/ba13_uncertainty/job.sh r320000
"$TC/python/bin/python" -B tools/linux/ba13_uncertainty/assess_holdouts.py
sbatch --wait --cpus-per-task=1 --mem=2G --time=00:20:00 \
  --output="$BA13_STUDY_OUT/chemical-%j.log" tools/linux/ba13_uncertainty/chemical_job.sh
sbatch --wait --cpus-per-task=1 --mem=2G --time=00:20:00 \
  --output="$BA13_STUDY_OUT/integrals-%j.log" tools/linux/ba13_uncertainty/integral_job.sh
sbatch --wait --cpus-per-task=1 --mem=2G --time=00:10:00 \
  --output="$BA13_STUDY_OUT/eos-audit-%j.log" tools/linux/ba13_uncertainty/eos_audit_job.sh
sbatch --wait --cpus-per-task=1 --mem=2G --time=00:10:00 \
  --output="$BA13_STUDY_OUT/final-analysis-%j.log" tools/linux/ba13_uncertainty/final_analysis_job.sh
```

Run `thermal_diagnostic.py` and `fit_diagnostic.py` on Slurm too, using
`"$TC/python/bin/python" -B tools/linux/measure.py "$TC/python/bin/python" -B SCRIPT`
and redirect to `thermal_diagnostic.log` / `fit_diagnostic.log`. They authenticate
and read both platforms, preserve the fitter and write separate diagnostics.
Then run `measurements.py`, `summarize.py` and `publish.py` for compact tables,
all 462 dispositions and hash-bound publication. `authenticate.py` rechecks the
original archive, nominal manifests and source owners; it writes the original
study's reauthentication receipt. Read-only authentication can run on the login
node. Publication refuses missing required evidence or inconsistent manifests.
A completed publication means complete **pilot evidence**, not qualified BA13.

## Actual execution and deviations

Actual jobs: 14587 initial build; 14591 thermal/unit diagnostics; 14592 nominal;
14595 sweep; 14596 chemical; 14632 withheld; 14641 independent integrals; 14642
held-target ladder; 14651 EOS audit/first profile mapping; 14652 held-target finest;
14654 completed profile mapping/P1 arithmetic. Resources are recorded per job.
Slurm accounting storage is disabled; measured process resources and receipts are
retained. The 53-case first mapping was rerun after the last two held-target
cases completed; both resource logs are counted.

Earlier failures: 14579/14583 missing sparse-checkout cmake files; 14585 absent
explicit Debug qualification prefixes; 14586 missing sparse test file; 14588
missing predeclaration path after sparse reapplication (no solves began).
An attempted `compute` partition submission was rejected; the actual partition
is phyall. Job14643 requested a new CMake target before regeneration and failed;
14651 regenerated and built it. These failures did not alter numerical owners.
A large combined execution request timed out in automatic review; splitting
file preparation and short submission calls succeeded. No scientific approval
was requested or inferred from that tooling failure.

Held-target predictions were written and staged before job14652 submission;
the first commit attempt failed for absent Git author identity. Commit 5d4fa20
followed submission, before examining those withheld values. Thus those two
predictions have a pre-submission file/index freeze but not a pre-submission
commit. This procedural limitation is disclosed rather than describing that
commit as earlier than it was. The original predeclaration and first-resolution
predictions were committed before their respective runs.

## Artifact and source binding

Large outputs live under
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/ba13-uncertainty`.
Each completed star directory has an atomic COMPLETE manifest, parameter/job
receipt and executable hash. Refused stages have no COMPLETE marker. Final
`STUDY_COMPLETE.json` inventories all retained case files, failed artifacts,
analysis logs, raw tables and selected build provenance/binaries. Compact copies
in this directory's `evidence/` are matched to those hashes. Scripts and the
predeclaration are separately hashed. No historical certificate was written.

The final incremental Git bundle uses 8e6ac58 as prerequisite; it includes only
the descendant study commits. It excludes large numerical files and is intended
for the existing office-Mac non-force push workflow. Verify it in a repository
containing the prerequisite before fetching its candidate branch.
