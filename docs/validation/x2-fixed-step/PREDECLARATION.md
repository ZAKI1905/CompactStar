# X2 prospective fixed-step experiment (UNVERIFIED SCIENTIFIC CANDIDATE)

Frozen before any X2 numerical run. Parent: 6cc7e20713e3c0210faa9036316650816f0290d5.
Authority: explicit X2 owner instruction; GOVERNANCE §§2–5; INV-02/03/04/06/14/15;
accepted ADR-0004 (proper volume), ADR-0009 (retained finite-pressure surface).
Class: validation-only numerical-method, documentation and generated evidence.
No production owner, fitter, frozen threshold, certificate or old result changes.
The Mac literature catalog is unavailable; equations are checked directly against
X1 MATHEMATICS.md and the authenticated test-side physical model. No new physics.

## Prior knowledge versus prospective claims

X1 A1 C contraction failed at all three targets; its mass cross-formulation gate
failed at all three; all A2 outcomes are INCONCLUSIVE. These facts and all old
bands remain fixed. We have read these results, not observed any X2 results.
X2 tests below are prospective; comparison to X1 is cross-method evidence, not a
blind prediction of previously unseen physical data.

## Method and regularity

Implement independent two-stage implicit Gauss–Legendre collocation (GL4),
c=(1/2±sqrt(3)/6), b=(1/2,1/2), A=[[1/4,1/4-sqrt(3)/6],
[1/4+sqrt(3)/6,1/4]]. Nominal global order four on smooth regular segments.
Reference: Hairer, https://www.unige.ch/~hairer/preprints/coimbra.pdf (collocation).
No adaptive RKF45 code or production EOS tables/interpolation in the new solver.
Fresh long-double momentum quadrature (96 points, independently generated nodes),
charge-neutral beta equilibrium inversion, physical constants identical to X1.
Use long double for solver arithmetic, compensated sums for integral increments.
Solve implicit structure stages by Newton with analytic Jacobian; reject on
nonconvergence, nonfinite/domain error, or stage residual >1e-15 scaled.

x=Hc-H, u=r², v=m/r³, D=1-2uv, A=D/(v+4πP):
u_x=2A; v_x=(4πε-3v)A/u; dV/dx=4πsqrt(u)A/sqrt(D).
B_x=1e-3*nB*dV/dx (B stored /1e57); JL_x=exp(6H)dV/dx on neutron
support; JC_x=exp(H)*sum(k_i mu_i/[3(hbar*c)^3])*dV/dx; Nmu_x=nmu*dV/dx.
Mass and metric identities must be evaluated independently of cancellation in v_x.

Split the central-to-muon interval at (Hc+Hmu)/2. Above that midpoint use
z=sqrt(Hc-H), dx/dz=2z. Then t=-sqrt(H-Hmu), dx/dt=-2t to Hmu.
Below Hmu use t=-sqrt(H-Hn), split exactly at Hb and Hn. Below Hn use
t=-sqrt(H), ending at the original positive Hs. All structure AND integral
RHS receive the Jacobian. Onsets Hn,Hmu are analytic; Hb is the retained
rho=1e10 g/cm³ cut. No step straddles a species onset or support discontinuity.
A species has k~sqrt(H-Hon), n~(H-Hon)^(3/2). In s these are one-sided
integer-power series; C's sqrt term becomes a smooth s² measured integrand.
Below an onset the absent species is identically zero. Center r~sqrt(x) makes
volume functionals nonsmooth in x; z makes them smooth. The positive surface
avoids the H=0 vacuum degeneracy. Verify these assumptions and actual convergence;
the singular u,v center means a uniform textbook Lipschitz theorem is NOT assumed.

Center startup x0=Hc*1e-4/N². With eH=dε/dH, pH=ε+P, v0=4πε/3,
a=2/(v0+4πP), c=-4πeH/5,
b=(a/2)*[-2*a*v0-(c-4πpH)/(v0+4πP)], initialize
u=a*x0+b*x0², v=v0+c*x0. Independently substitute the series into both TOV
equations. Integrals include their leading central-volume term INCLUDING C;
next omitted term is O(x0^(5/2)), not zero central C. u remainder O(x0³),
v remainder O(x0²) has a decaying outward homogeneous mode; quantify startup
sensitivity, do not claim a certified center remainder. Central metric increment
starts at x0. Alternative startup factors 1e-2 and 1e-6 are diagnostic only.

## Fixed schedule and stopping rules

Target 0 pilot: checks, then N=32,64,128,256 per segment (five segments).
If the correct equations/identities fail, Newton fails, persistent increasing
errors occur, or regularity assumptions fail: stop expansion, preserve results.
Non-asymptotic but decreasing increments permit the declared ladder to continue.
After pilot cost estimate, run N=32,64,128,256,512,1024 for targets 0,10,20,
reusing completed target-0 outputs. Withheld N=2048 and 4096 at all three targets
are allowed only if decreasing increments have not reached a stable order-four
regime or the cross-method mass/C question remains unresolved, and projected
cost fits the ceiling. No other resolution, no 21-target study.

Common-B secant uses original target strings and individual surface pressures;
max 10 shots, |B/Btarget-1|<=2e-16 (long-double computational alignment, not a
physical bound). Fixed-H shots at both retained finest X1 Hc values separate
shooting from fixed-H formulation differences. Shot sensitivities are diagnostics.

Predeclared ablations on target 0 only at N=512,1024 (and 2048 if already released):
center factors above; EOS order 192; binary64 arithmetic version; plain rather
than compensated sums; replace onset s by H while retaining exact segmentation.
The unregularized ablation is expected to lose order for C (~3/2 quadrature error)
and is not eligible for GL4 nominal-order claims. These diagnose components,
not substitutes selected to pass. X1 diagnostic replays may use ONLY its already
observed four tolerances, with center changes, accumulation changes, or passive
mass/B identity diagnostics; no smaller tolerance. Preserve distinct identities.

## Scoring, falsification, uncertainty

Report signed increments, p=log2(|d_N/d_2N|) and Richardson |d|/15 for all
M,R,r_neutron,C,L,muinf,Hc and achieved B. Order-four empirical PASS requires
two consecutive p in [3.5,4.5] with same-sign increments above a disclosed
arithmetic floor; otherwise INCONCLUSIVE (roundoff/crossover) or FAIL (persistent
nonconvergence, wrong equations). Do not fit an order after observing results.
Richardson is a conditional estimate, eligible only after that gate; it is not
a global enclosure. Report absolute and target/B0 paired errors separately.
Cross-method descriptive agreement: differences <= sum of eligible empirical
estimates plus measured component sensitivities. No arbitrary roundoff allowance;
if estimates are unavailable the classification is INCONCLUSIVE. X1 mass gate
is NEVER rescored with these estimates. Anchor residuals use unchanged bands.

Correctness checks: EOS epsilon+P=nB*mu <=1e-13 relative, charge/beta <=1e-13;
dP/dH finite difference <=1e-6; analytic transformed power moments <=1e-14;
center series residual decreases at its predicted power; independent off-stage
mass and B integral residuals decrease, and at finest level are <1e-7 relative;
independently integrated metric yields redshifted-mu drift <1e-12. Manufactured
smooth GL4 test must measure p in [3.8,4.2]. A failed criterion is retained.
Residuals reconstructed through the collocation polynomial at additional Gauss
nodes test integrated balance, not just the algebra of the same stage RHS.

Attribution requires reproducible component isolation and effect magnitude;
otherwise report cause INCONCLUSIVE, not a speculative unique explanation.
Compute historical headroom T_X-drift_bound for EVERY BA13 component and all
retained samples before small-residual pursuit. This is conditional remaining
numerical/nonlinearity budget, NOT uncertainty. For nonzero ratios, equal endpoint
relative errors e give ratio uncertainty <=2e*|ratio|/(1-e); solve this inequality
for allowed e using headroom. Profiles/matrices use their historical norms;
species support needs exact topology, not a scalar relative allowance.
No absolute physical accuracy threshold has been supplied: do not invent one.
Certificate accuracy also requires bounded EOS, stability, center, arithmetic,
root sensitivity and propagation, absent from mere cross-method agreement.
Assess a posteriori defects with stability bounds or interval integration as a
future justified enclosure; explicitly name every conditional bound.

## Resources and artifacts

All substantial numerical work on Slurm, one pinned logical CPU, 1 GiB memory,
900 CPU seconds total initial numerical ceiling. Pilot first; cumulative getrusage
(parent+children), per-process RLIMIT_CPU and wall timeout, forecast before more
jobs. Stop at 850 CPU seconds to reserve receipts; hard total 900. No BNV runs.
Builds/read-only analysis are recorded separately. No automatic budget extension.
Commit source and execution hashes before runs; store raw large artifacts outside
Git with source-bound SHA256 manifests, compact results/receipts/report in Git.
No merge, promotion or historical failure changes. If GitHub auth is unavailable,
produce and verify an incremental bundle against 6cc7e20.
