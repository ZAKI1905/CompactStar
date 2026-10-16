# Preserved audit launcher failure

Slurm 17032 completed all arithmetic, raw replay, support and headroom assertions,
then failed to locate Git in the compute-node PATH during the scope check.
It wrote no success receipt and is NOT an audit PASS. Its traceback and Slurm
FAILED state are retained. No new stellar solve occurred. Retry the identical
committed audit with the authenticated toolchain/bin added to PATH and an
external /usr/bin/time process receipt. No audit assertion or criterion changes.

Final Slurm records for jobs 17024–17031 had already expired when queried;
the preserved allocations, controller receipts, individual return codes, raw
outputs and stderr remain the evidence. Do not invent completed Slurm receipts.

Slurm 17033 then failed before Python because compute nodes lack /usr/bin/time.
Its FAILED/127 receipt is also preserved; numerical CPU is zero for that attempt.
The next launcher uses only Python getrusage and the authenticated PATH. Both
failed allocations lasted one second; conservatively charge their full two-CPU
allocation wall-times (four CPU-seconds total) in the resource upper accounting.
