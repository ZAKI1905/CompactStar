# Before-results launcher correction

Job 17021 failed before `execute()` and before any numerical test or solve:
compute nodes lack Git on the default PATH. Sources passed authentication.
The original Slurm output and allocation receipt are retained externally.
The launcher now prepends the existing authenticated toolchain Git wrapper.
It measures children with Python getrusage, avoiding a second unverified
compute-node dependency on GNU time. This is the repository's established
measurement mechanism (`tools/linux/measure.py`). Child peak RSS is explicitly
a cumulative upper bound, not falsely labelled a per-child measurement.

Slurm allocated two hardware threads for the first job's one requested core.
The replacement requests nomultithread and pins Python and every child to one
logical CPU using sched_setaffinity, independently of Slurm allocation rounding.
The initial job had zero numerical work and zero displayed elapsed seconds;
its unmeasured subsecond authentication/launcher CPU is bookkeeping, not hidden
numerical work. No numerical source, prediction, target, tolerance, scoring
definition or scorer changed. Regenerate the execution hash gate before the
first numerical run. This correction was made without observing numerical data.
