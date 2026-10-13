# Reproduction and evidence layout

Entry 4f3cf97; standalone descendant candidate. Read PREDECLARATION.md and
EQUATIONS_AND_LIMITS.md first. No production owner or certificate is modified.

The scripts intentionally bind this pilot's root:
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/ba13-onset-pilot`.
A future repetition must use an isolated path and its own resource ledger;
do not rerun fresh-output commands over this evidence tree. Paths in scripts
must be changed and committed for a repeat, with the original predictions retained.

1. Run `python3 -B tools/linux/ba13_onset/authenticate.py` (read-only hashing).
2. Source `tools/linux/eku_environment.sh` on login node, then submit `build.sh`
   with partition phyall, one CPU, 2 GiB, 10 minutes. Reuses the authenticated
   original libCompactStar and dependency libraries; generates exact original
   helpers through the hash-checking extractor. Build CPU is reported separately.
3. Submit `freeze.py` and `discovery.py` through Slurm and `tools/linux/measure.py`.
   Discovery reads only prior authenticated profiles; all its CPU counts toward cap.
4. Commit PREDECLARATION.md, predictions.json, eos-predictions.json and tools
   before submitting `controller.py`. Actual gate commit: 650871e. The controller
   requires a clean tree at startup and a Slurm job ID.
5. Submit the controller with one CPU, 2 GiB, 59-minute allocation. This is one
   sequential allocation, not 31 separate jobs; each subprocess is limited to ten
   minutes. A conservative elapsed-time budget plus all discovery CPU prevents
   spending 3600 additional numerical CPU-seconds. BLAS/OMP threads are one.
6. After completion, submit `score.py` on Slurm using the reserved budget. Rules
   are those committed before outcomes; scores may not redefine any threshold.
7. Reauthenticate, publish compact evidence with hashes of all external files,
   and create an incremental Git bundle with prerequisite 4f3cf97. No direct push
   is inferred from credentials or remote configuration.

Actual initial job records: 14686 could not find /usr/bin/time (no numerical
work); 14687 froze predictions through the repository resource wrapper; 14688
built; 14689 ran discovery; 14690 ran the fixed pilot. An earlier sandboxed
submission could not resolve the Slurm controller and submitted no job. The
network-enabled submission succeeded. The first add omitted sparse paths;
8bb72e2 subsequently committed the complete predeclaration, before any withheld
calculation. e657efb committed radial/secondary numeric predictions; 650871e
committed EOS predictions and all executed numerical tools. No timestamp-only
or post-submission prediction is presented as a blind commitment.

External directories: `tables/` strict-ordering/identity tables; `runs/` original
owner stars/thermal output and per-case COMPLETE hashes; `analysis/` reconstruction
outputs; discovery.json; high-precision.json; oracle logs; resources.json;
scores.json; metrics.json. Partial/refused outputs remain in place. Compact Git
evidence references the external manifest; large profiles and tables stay outside Git.

Scoring job14691 returned FAIL under unchanged rules. Post-outcome jobs14693
(diagnose.py) and14695 (geometry_diagnostic.py) analyzed only already-declared
profiles; neither adds a case nor changes a score. Their CPU is included in the
resource summary. The geometry subtraction is explicitly post-hoc discovery,
not a new production implementation or repaired holdout result. A prospective
successor must use new withheld data and its own predeclaration.
