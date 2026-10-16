# Target-0 pilot and ladder release

Slurm 17024, source d313cc3, numerical CPU 2.041637 seconds, affinity [4].
All analytic checks PASS: GL4 manufactured order 4.00025; center u residual
order 2.00004 and v residual order 0.999982, matching the declared expansion.
Momentum EOS thermodynamic residual <4.62e-17; finite-difference dP check 4.47e-10.
At N=256 independent off-stage mass and baryon closure residuals are 6.37e-12
and -6.09e-11; independent metric drift 8.06e-14. Residuals decrease.
All seven physical scalar outputs and Hc have both measured orders in
[3.965,3.993]; B alignment is at the root/representation floor and has no
convergence order claim. No domain/Newton/regularity failure is observed.

Release the predeclared N<=1024 ladder for targets 0,10,20. Pilot cost is
approximately 0.0043 CPU seconds per N per common-B solve. Remaining base
ladder forecast <30 CPU seconds; even 3x reserve <90, well below 900.
This is an empirical correctness/convergence PASS for the pilot, not a rigorous
reference qualification. No historical X1 gate is changed.

Sandbox-local sbatch initially could not resolve ekucn1; submission outside the
network sandbox succeeded. The unsuccessful submission produced no allocation
or numerical work. The actual numerical work ran only on Slurm.
