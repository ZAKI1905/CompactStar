# Preserved binary64 diagnostic failure

Slurm 17028 completed five target-0 ablations at N=512, then double-n512 failed
`stage residual`. The loop stopped at 32*epsilon, about 7.1e-15 in binary64,
while the predeclared acceptance gate is 1e-15. This is an implementation
inconsistency in the binary64 ablation. No extended-precision result failed.

Preserve its empty stdout, stderr, failed Slurm receipt and FAIL disposition.
A separately named generated binary64 copy stops Newton at min(32*epsilon,5e-16),
without changing the 1e-15 gate, method, resolution or physical equations.
Freeze that copy and continuation runner before its results. Resume the already
specified, unaffected center/EOS/plain/onset/fixed-H diagnostics. The original
GL4 source, executable, gate and all ladder outputs remain unchanged.
