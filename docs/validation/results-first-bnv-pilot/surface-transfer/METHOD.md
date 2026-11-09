# Finite local surface transfer

Scientific candidate; conditional on the already accepted numerical input bands.
This is the calculation specified in [NEXT_CALCULATION.md](../cut-assessment/NEXT_CALCULATION.md).
The source equations, literal EOS conventions, units and observable derivatives
are in [DERIVATION.md](../cut-assessment/DERIVATION.md). No physical definition,
surface law, fitted coefficient, threshold or production solver is changed.

## Terminal flow and its finite remainder

Use s=ln(Ps), H=Hs, D=1-2M/R, V=M/R^3+4*pi*pgeom, K=nu_s+H,
alpha=P/P_H, W=4*pi*R*D/(V*sqrt(D)). P_H differentiates the literal retained
EOS, including its separately widened binary64 masses and h0. The small literal
thermodynamic convention defect is preserved. Integrate at fixed central Hc:

    H_s = alpha
    R_s = -D*alpha/(R*V)
    M_s = 4*pi*R^2*egeom*R_s
    B_s = -W*alpha*nB*1e54 = (Ne)_s
    K_s = kappa = 4*pi*(egeom+pgeom)*alpha/V
    C_s = -kappa*C - W*alpha*exp(H-K)*cv*cvfactor
    G_s = -kappa*G + shell*e_e*e_e^T
    shell = -W*alpha*exp(H-K)*chi*1e54.

G0 has only the redshift term. Ne includes the shell; Nn and Nmu have exactly
zero terminal shell integrand because the entire tube is below their onsets.
Their **fixed-B** responses are nonzero. Ltilde transforms by exp(-6*Delta K),
muB/P0/P1 by exp(Delta K), and mu_n_profile follows the same equilibrium muB
change in every one of its original 2049 component scales. Radius and gravity
are reconstructed in their historical units. The nonlinear envelope and
redshifted surface temperatures are evaluated from the actual new g, H and K.
Their log rates are g_s/(4g)-0.605*kappa and
g_s/(4g)+0.395*kappa-alpha respectively; finite exponential bounds retain all
higher powers. Cstar's shell and redshift terms are both retained.

The terminal field and its first derivative use interval arithmetic at 90
digits. A 32-term convergent series for each positive Fermi integral has an
explicit geometric tail and derivative-tail bound. These bounds avoid the
cancellation in a direct interval evaluation of the closed pressure formula.
Independent momentum/closed-form EOS identities and numerical differentiation
are tested. Bound endpoints are serialized as exact decimals of the computed
binary numbers; they are not rounded to a short printed upper limit.

Original pointwise reference bands and C-CONT curvature enclose the on-law
starting states on [-1e-6,0]. R and g bound M through the exact Schwarzschild
surface-gravity relation. The pressure law bounds H by monotone EOS inversion.
Retained offset extrema are also included in the terminal box, for validation
of those local maps. Including sampled extrema is **not** proof of a uniform
equilibrium-family continuation outside that interval.

A proposed relative 2*epsilon tube is checked for strict Picard containment;
it is a domain to test, not a reported error multiplier. With scaled structural
Jacobian J, the calculation verifies tau*||J||<1. It then uses the supremum
field and Jacobian to bound the complete finite flow and its variation. The
structural linearization remainder is bounded componentwise by
tau^2*sum_j |J_ij|*sup|F_j|/2. Variation is bounded by exp(tau*||J||).
No small first-order product is substituted for this remainder.

For the chemical inverse put U=G_l^-1 and Z=U+11^T/G0. With
ell=integral(exp(Delta K)*shell ds), Sherman-Morrison gives exactly

    Z_new = exp(Delta K) * [Z - ell/(1+ell*Uee) * Ue*Ue^T].

The retained susceptibility construction has positive definite cooperative
lepton block (nonpositive off-diagonal), G0>0, and 0<=Uij<=Zij. Thus the
accepted Z bands bound the rank-one response without inventing an uncertainty
for an unbounded raw inverse. The calculation verifies
1-|ell|*sup(Zee)>0 for both signs, uses this lower denominator, and retains
the finite redshift exponential. Direct matrix inversion independently tests
the identity. These are conditional input enclosures, not a new rigorous
certification of the underlying stellar reference.

## Fixed-baryon matching and tangents

Let x=(B-B0)/B0, s0(x)=a*x+b, and w(x',tau) be the terminal change in B/B0
from an on-law initial state x'. To reach final x and pressure s0(x)+delta,

    x = x' + w(x',tau),       tau = delta + a*w(x',tau).

The source uses the original production B0, original binary64 fitted law and
epsilon=2.1e-12. If A bounds |B_s|/B0, then

    |tau| <= epsilon/(1-|a|*A),       |x-x'| <= eta=A*|tau|.

The accepted positive B_Hc condition gives a valid family coordinate in its
stated domain. At first order the mandatory central correction is
Delta Hc=-(B_s/B_Hc)*delta, about 1.2e-26 at the cut limit. The code propagates
the local map and its variation at high precision; it does not attempt a
binary64 central re-shoot. This matching is the finite form of
Q_s|B=F_Q-Q_Hc*B_s/B_Hc, previously independently verified.

For normalized q, the finite bound is terminal change plus
sup|q_x|*eta. The on-law derivative bound follows from retained secants,
their existing pointwise bands and accepted curvature. Those existing error
terms enlarge the domain on which the C-CUT response is bounded; they are not
added a second time to the historical error ledger.

For species c_n=c_mu=0, c_e=1, define T=N_x/B0 and let w_x denote the initial
family derivative of the terminal change at fixed flow length. Differentiating
the matching equations gives the finite identity

    dx'/dx = (1-a*w_tau)/(1-a*w_tau+w_x)
    T_new(x) = c + [T0(x')-c]*dx'/dx.

The first variation bounds w_x by
|tau|*||grad(B_s/B0)||_1*exp(||J||*|tau|)*||Y0_x||_inf. The denominator is
separated. Therefore

    |T_new-T0(x)| <= sup|T0_x|*eta
                    + sup|c-T0|*|w_x|/(1-|a|*A-|w_x|).

Use the accepted count curvature identity sup|T0_x|<=S_N*K_N/B0. This avoids
the inconclusive mixed-count Richardson remainder. It bounds tangents for
constant residual offsets delta, not derivatives of an arbitrary function
delta(x); no smoothness of the measured cut residuals is invented. Existing
C-TAN accounts separately for the historical fixed-cut/path distinction and
correlated quotient remainders.

## Domain obstruction and ledger propagation

The bounds above apply uniformly when both x and x' are in [-1e-6,0]. A
sufficient final domain is [-1e-6+eta,-eta]. At either endpoint, one residual
sign requires an on-law state outside the accepted domain. The actual retained
offset stars permit evaluating the local terminal maps there; their finite
samples do not bound the unsampled interior-family derivatives. In particular,
all six endpoint count-curvature grid convergence checks fail standard S.

The missing evidence is a uniform continuation of the existing normalized
point/slope bands and count curvatures, with positive B_Hc and fixed topology,
on [-1e-6-eta,-1e-6] and [0,eta]. It is not another derivative of the terminal
pressure flow. No endpoint extension is accepted in this delivery.

The successor ledger preserves every old key/value/type, term, threshold and
qualification. It records supported interior U separately and leaves the full
C-CUT bound, relative baseline contribution and full utilization null. Tangent
rows have no observable baseline denominator but still require the original
baseline tangent in their paired differences and all-row fit. Consequently
none of the 399 approximate rows is promoted.

A separately marked **unauthorized conditional scenario** applies the same
bands on the two missing endcaps. For original fit norm f, add
(f+2)*U to the numerator (f*U at row zero), U to each relative baseline, and
2*U*sliver/h at the final attained-coordinate endpoint. Retain the original
known numerator N and baseline b, and calculate

    scenario = (N + added_numerator + added_endpoint)/(1-b-U).

For absolute tangents omit U from this denominator. The existing endpoint
curvature/sliver terms and the ADR-0022 total reference fit responses remain
unchanged. The new term is solely the finite transverse-cut perturbation,
propagated through the historical linear fit and nonlinear normalization.
No threshold is inferred from available headroom. An independent Decimal audit
checks this arithmetic and refuses omission of the baseline, endpoint or
qualification restrictions.
