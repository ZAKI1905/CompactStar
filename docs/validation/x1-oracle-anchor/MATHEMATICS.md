# X1 equations and independent mathematical audit

Candidate only. Governing meanings: INV-02/03/04/06/14/15. This is a coordinate
change to the established validation oracle, not a new physical EOS or owner.
The Mac-only literature paths in AGENTS.md are unavailable on this cluster.
The derivation below checks the existing test-side equations directly; it does
not claim a new literature-based physical model. Source: frozen starting
`tests/eos/structure1/{local_oracle,tov_oracle}.hpp` and
`tests/bnv/onset_structure_oracle.cpp`. GSL's local `ode-initval2/cstd.c` and
`rkf45.c` define the unchanged error controller and RKF45 tableau.

Write x=Hc-H, u=r², v=m/r³, D=1-2uv, A=D/(v+4πP), with geometric P,epsilon.
The TOV chain rule gives du/dx=2A, dv/dx=(4πε-3v)A/u and dr/dx=A/r.
Thus dV/dx=4πrA/sqrt(D), denoted W. Checks independently reconstruct
dV/dx from 4πr² dr/dx/sqrt(1-2m/r), and dm/dx from both the mass equation
and u,v. These algebraic residuals are consistency checks, not global errors.

H=ln(mu_B/h0), h0=m_p+m_e. Thermodynamic consistency requires
epsilon+P=n_B mu_B and dP/dH=epsilon+P. Both are independently tested against
the direct momentum quadratures; the latter uses a fixed finite difference.
Charge neutrality is n_p=n_e+n_mu; beta equilibrium is mu_n=mu_p+mu_e on
neutron support. Muons start when mu_e=m_mu. Explicit cuts are
H_n=ln(m_n/h0), H_mu=ln((sqrt(m_p²+m_mu²-m_e²)+m_mu)/h0).

At each high-density onset side choose t=-s=-sqrt(H-H_on). Then H=H_on+t²,
dx/dt=-2t=2s>=0. ALL six right sides multiply by 2s. In particular
dB/dt=1e54*n_B*W*2s (the stored B channel is B/1e57),
dJ_L/dt=1_(H>H_n)*exp(6H)*W*2s,
dJ_C/dt=exp(H)*[sum_i k_i mu_i/(3(hbar c)^3)]*W*2s,
dN_mu/dt=n_mu*W*2s. At t=0 the transformed measure vanishes exactly;
the L support endpoint has no finite atom. The oriented identity is
integral_(H_on)^(H_high) F(H)dH = integral_0^S 2s F(H_on+s²)ds.
Analytic power moments p=0..5, both actual onset RHS transformations,
species absence/presence and the zero endpoint are tested independently.
Below neutron onset use the original x coordinate. Above H_mu use its own
s; between H_mu and H_n use neutron s, retaining the envelope cut H_b.
There is no claim that this also regularizes the *lower* side of a higher
species boundary or the stellar center.

At the finite-pressure surface nu_s=ln(1-2M/R)/2. Hydrostatic equilibrium gives
nu(H)=nu_s+H_s-H, and mu_B exp(nu)=h0 exp(nu_s+H_s). This identity is
checked against EOS mu_B, but nu is constructed from H; this is not an
independently integrated metric-conservation test. Output L is
1e-36*exp(-6(nu_s+H_s))*J_L and C/Tinf is
1e15*k_B²*1.602176634e33*exp(-(nu_s+H_s))*J_C. The inherited label Cstar
denotes this linear heat-capacity coefficient; a common fixed Tinf multiplies
all values and cancels in relative tests. M is geometric km, not solar masses.

The four integrated functionals do not feed back into the structure RHS.
Each RKF45 trial therefore starts their local bases at zero, yielding direct
increments. On acceptance, Kahan adds each increment to its running total;
rejected trials do not modify sums or compensation. GSL's y-controller sees
the full proposed totals, with unchanged eps_abs=tol*1e-5, eps_rel=tol.
This avoids inaccurately obtaining increments by subtracting large totals.
Stage sums and the shared 96-point EOS quadrature remain ordinary binary64;
compensation is not an all-operations roundoff enclosure. The initial small
central B,L,N_mu contributions match the old source; initial C remains zero.
The omitted central C scales as radius cubed but has no certified bound here.

## Hermite order at onset

For a massive species appearing continuously, k~sqrt(delta H), n~delta H^(3/2),
epsilon contributes O(delta H^(3/2)), and P contributes O(delta H^(5/2)).
There is no finite jump in the exact first derivative d epsilon/dH from that
species at the boundary, although rapid nearby variation can look like a kink
on sampled tables. This observation limits the review's language about a jump.
For a simple onset with dH/dr nonzero, m acquires an O(delta r^(5/2)) term.
The TOV nu' inherits such a term, and nu contains O(delta r^(7/2)). Its fourth
derivative can diverge like |delta r|^(-1/2). Hence the classical bounded-nu''''
Hermite O(h^4) remainder is not established on a cell straddling onset; a
fractional h^(7/2) contribution is possible. Exact nodal derivatives still
remove the linear interpolation O(h²) term. A tiny prefactor or apparent h^4
at sampled grids is not a uniform proof. One-sided smooth segments require
their own derivative bounds. This audit does not adjust any anchor band.

## Limits of absolute-reference claims

Both formulations share EOS quadrature, binary64 constants and inversion,
central asymptotics, finite surface convention and RKF45 tableau. Agreement
can miss those common errors. Event splitting and the coordinate algebra are
mathematically correct for the specified ordering, but do not supply a global
error estimate. Empirical U requires indefinitely continuing geometric decay;
four tolerances cannot establish that premise. A rigorous absolute reference
would need an integration remainder/stability estimate (or validated ODE
enclosure), EOS/quadrature/inversion and roundoff control, center truncation,
and propagation through the common-B root with a bounded nonzero derivative.
A predeclared independent fixed-step Gauss-Legendre method in s, with center
series control and step halving, is the smallest useful cross-method next
test, and still requires justified remainders before claiming rigorous bounds.
