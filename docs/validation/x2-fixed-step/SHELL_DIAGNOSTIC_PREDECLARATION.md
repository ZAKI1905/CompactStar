# Prospective mechanism test after X1 trace observation

This is explicitly a post-observation hypothesis, not part of the original blind
X2 convergence forecast. No new resolution or tolerance is introduced, and no
original score is changed. Before these diagnostic runs, we observed that the
X1 regularized fixed-H trace adds JL=1.04239e-7 BELOW neutron support at 2e-13.
The exact integral there is zero. C's final-shell increment also dominates its
remaining integration error. Source review shows that H=Hc-(Hc-Hn) can round
above Hn at the first shell stage, enabling the neutron-only L channel.

Test target 0 only, at the original four tolerances, both common-B and fixed at
the original regularized finest Hc. Generate three distinct X1 diagnostic copies:
1. support: use the segment's neutron-support flag for L, no other RHS change;
2. endpoint: additionally remove the rounded neutron density at lower-side EOS
   endpoint evaluations, rebuilding EOS from the remaining densities;
3. shell-s: additionally integrate the sub-neutron shell with t=-sqrt(H),
   applying its Jacobian to every equation and integral (the X2 shell map).

Before results: support should make the below-neutron JL increment EXACTLY zero.
If not, FAIL the mechanism test. C need not improve when only L is corrected:
changing which channel controls adaptive steps may change its cancellation.
Comparing support versus endpoint isolates phantom-neutron C. Comparing endpoint
versus shell-s tests shell parametrization/integration of the existing smooth
sub-neutron EOS. Trace contributions must distinguish accumulated interior errors
from shell errors. Center/roundoff ablations already observed cannot explain the
~1e-12 residual by themselves. Do not assert a unique cause of every C ratio
unless the isolated changes account for its scale; otherwise INCONCLUSIVE.

These are diagnostic candidates, never a replacement X1 ladder. Record source
hashes, all failures and CPU usage. Forecast <15 CPU seconds; same 900-second
cumulative ceiling, one pinned CPU via Slurm. No other target or tolerance.
