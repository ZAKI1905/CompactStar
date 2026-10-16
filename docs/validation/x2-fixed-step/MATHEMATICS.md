# X2 formulation and regularity audit

UNVERIFIED SCIENTIFIC CANDIDATE. Physical meanings are inherited from authenticated
X1 MATHEMATICS.md: geometric r,m in km; epsilon,P in km^-2 after the same binary64
unit conversion; n in fm^-3; B a count; Cstar is C/Tinf, not a finite-temperature
trajectory. Governing contracts: INV-02/03/04/06/14/15, ADR-0004 and ADR-0009.
The Mac literature catalog is unavailable on this cluster. No physical EOS change
or production ownership change is inferred from this test-side derivation.

## Independent momentum EOS

For each species k_i=(3*pi²*n_i)^(1/3)*hbar*c and mu_i=sqrt(m_i²+k_i²).
Evaluate epsilon_i=3*n_i*integral_0^1 z² sqrt(m_i²+k_i² z²) dz and
P_i=n_i*k_i²*integral_0^1 z⁴/sqrt(m_i²+k_i² z²) dz. A fresh Legendre recurrence
and Newton root solve construct the 96-node rule in long double; no GSL quadrature,
old LocalOracle, table or interpolation implementation is included.
Charge neutrality gives np=ne+nmu; beta equilibrium gives muB=mup+mue.
Inverting these monotone relations yields the same cold free-gas physical model.
Below neutron onset np=ne and the two-body k² relation is evaluated directly.
Constants are linked from the same authenticated Zaki archive; this shared physical
input is disclosed. G and c match the authenticated GSL constants exactly.

Thermodynamics gives dP/dH=epsilon+P=nB*muB and
(d epsilon/dH)=muB²*[S_n+S_p*(S_e+S_mu)/(S_p+S_e+S_mu)],
where S_i=mu_i*k_i/[pi²*(hbar*c)³]. This derivative enters only the analytic center
startup. Local quadrature consistency, charge, beta equilibrium and finite-difference
pressure checks independently test these identities. EOS192 sensitivity is measured,
not promoted into a quadrature remainder enclosure.

## TOV and all functionals

The radial equations are dm/dr=4*pi*r²*epsilon and
-dH/dr=(m+4*pi*r³*P)/[r²*(1-2m/r)]. Define x=Hc-H, u=r², v=m/r³,
D=1-2uv and A=D/(v+4*pi*P). Chain rule yields

    du/dx = 2A
    dv/dx = (4*pi*epsilon-3v)*A/u
    dV/dx = 4*pi*sqrt(u)*A/sqrt(D) = W
    d(B/1e57)/dx = 1e-3*nB*W
    dJL/dx = 1_(H>Hn)*exp(6H)*W
    dJC/dx = exp(H)*sum_i[k_i*mu_i/(3*(hbar*c)^3)]*W
    dNmu/dx = nmu*W

The independent stage mass channel is dMint/dx=4*pi*epsilon*sqrt(u)*A.
No integral feeds the structure equations. Two implicit Gauss stages solve only
u,v with the analytic Jacobian. Accepted integral increments use Kahan sums.
Structure updates and stage algebra are ordinary long-double operations, not
interval arithmetic. The pilot records 64 mantissa bits on this cluster.

For each map multiply EVERY RHS by dx/dt. Center: t=z=sqrt(Hc-H), dx/dz=2z.
High-side onset: t=-s=-sqrt(H-Hon), dx/dt=-2t>=0. Surface shell: same map with
Hon=0, ending at t=-sqrt(Hs)<0. No high-side endpoint atom is present.
At Hn, JL is stopped by segment topology and interior Gauss nodes; no endpoint
sample enters its quadrature. Hb is the original rho=1e10 g/cm³ event, and Hs
is the original target-specific finite-pressure cutoff. Hn and Hmu are analytic.

At the surface nu_s=log(1-2M/R)/2 and K=nu_s+Hs.
Ltilde_Me=1e-36*exp(-6K)*JL, C/Tinf=1e15*kB²*1.602176634e33*exp(-K)*JC,
and mu_inf=h0*exp(K), exactly the X1 conventions. Independent radial metric
integration below tests the enthalpy-derived metric, rather than reconstructing
both sides from H and calling that an integration test.

## Independently checked center coefficients

Let epsilon(x)=ec-eH*x+O(x²), P(x)=pc-(ec+pc)*x+O(x²),
u=a*x+b*x²+O(x³), v=v0+c*x+O(x²). The zeroth mass balance gives
v0=4*pi*ec/3. The u equation gives a=2/Q, Q=v0+4*pi*pc.
The v equation then requires c=(-4*pi*eH-3c)/2, hence c=-4*pi*eH/5.
Expanding A to first order gives

    b=(a/2)*[-2*a*v0-(c-4*pi*(ec+pc))/Q].

These formulas are independently substituted in the ORIGINAL x equations in
--check. Residuals scale as x² in u and x in v: observed orders 2.00004 and
0.999982. Leading central volumes are included for all integrals, INCLUDING C.
The omitted functional term is O(x0^(5/2)); x0=Hc*1e-4/N² makes this O(N^-5).
The neglected v term is O(x0²). The leading homogeneous v error decays outward
as x^-3/2, but no certified constant or coupled global stability bound is supplied.
Changing the center factor by 100 in either direction tests startup sensitivity.
This is stronger evidence than assuming zero central C, but not a center enclosure.

## Actual smoothness and order limitations

Two-stage Gauss collocation has classical order four on smooth, regular problems;
see Hairer, https://www.unige.ch/~hairer/preprints/coimbra.pdf, collocation discussion.
That theorem is not applied across a species boundary or to the singular center.
For neutron onset k_n²=muB²-mn² has a simple positive zero in H-Hn.
At muon onset the equilibrium map has positive finite derivative dh/dmue, so
mue-mmu has a simple zero as well. In s>=0 the high-side species functions have
one-sided integer-power expansions: k~s, n~s³, and the measured C integrand ~s².
The implicit charge-neutral root remains locally invertible. Below each onset the
absent species is identically zero; those are separate branches. No assertion of
global C^4 regularity in H is made. In raw H, a heat-capacity sqrt onset can limit
quadrature order to 3/2; the three raw-H ablations have not reached that asymptote.

The central volume integrand scales as sqrt(x) but is smooth in z. The structure
system is still regular-singular at z=0, which is excluded by the analytic startup.
No uniform Lipschitz bound as x0 shrinks is proved. Positive Hs excludes the vacuum
endpoint; derivative sizes near Hs may still matter numerically. These statements
justify a fourth-order design and the prospective order test, not a rigorous global
order constant. The observed physical-scalar orders remain near four throughout
N=32..4096. Raw-H ablations have large C errors and changing orders.

## Independent mass, baryon and metric residuals

On each step construct the degree-two collocation polynomial u_h(t),v_h(t).
At FOUR additional Gauss nodes (none is a collocation stage), compute
r_h=sqrt(u_h), m_h=v_h*r_h³, and dr_h/dt=u_h'/(2*r_h). Integrate separately

    dM_check/dt = 4*pi*r_h²*epsilon(H(t))*dr_h/dt
    dB_check/dt = 1e-3*4*pi*r_h²*nB(H(t))*dr_h/dt/sqrt(1-2*m_h/r_h)
    dnu_check/dt = (m_h+4*pi*r_h³*P)/(r_h²*(1-2*m_h/r_h))*dr_h/dt.

Compare accumulated M_check to m_h at the surface, B_check to the primary B
channel, and nu_check to Hc-H at segment endpoints. The metric difference is
also the logarithmic drift of redshifted chemical potential. These integrals use
polynomial derivatives, not the same A/r RHS, and can detect a nonzero defect.
They still share an approximate trajectory and finite quadrature: small closures
are not a proof of proximity to the exact solution. A sampled u-equation defect
is also recorded; its maximum ~7.16e-8 at N=4096 is not a global error bound.

The X1 passive diagnostic applies its own RKF45 fifth-order weights to the radial
mass integrand at accepted stages, without adding a controlled state channel or
altering acceptance. Eight baseline replays reproduce ALL nine retained outputs
bit-for-bit. This matters because nonlinear invariants such as m=v*u^(3/2) need
not be preserved by Runge–Kutta discretization even when their continuous chain
rules are exact. The observed integrated residual detects precisely this distinction.
