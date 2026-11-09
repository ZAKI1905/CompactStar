# Fixed-B surface-cut response

UNVERIFIED SCIENTIFIC CANDIDATE. This derives the historical observables, without
changing them. Governing definitions: accepted `../OWNER_DECISION.md:22–37`,
`../evidence/INDEPENDENT_REVIEW.md:536–579`, `../static-completion/MODEL.md:5–47`.
Implementation sources are `tests/bnv/x2_fixed_step.cpp:106–118,183–210`,
`tests/bnv/x2_uncertainty.cpp:55–89`, and
`tools/linux/ba13_closure/assess.py:58–69` (paths relative to repository root).

## Coordinates, EOS and terminal identities

Write h=Hc, s=ln(Ps), x=(B-B0)/B0, and s0(x)=a*x+b. B0 is the original
production baseline, not the slightly different X2 count. The original reviewer
calculation is binary64 `polyfit(achieved_x,log(Pgeom),1)`; multiplying pressure
by the fixed geometric/CGS conversion changes the intercept only. It is replayed
exactly, including the original 2.0676793610618915e-12 maximum residual. The
high-precision least-squares slope is a diagnostic and is NOT substituted for a.
The historical frozen-observable fit is an entirely different fit and is unchanged.

Use R and M in km, D=1-2M/R, p and e in geometric units in TOV, and
V=M/R^3+4*pi*p. The surface lies in the proton–electron branch, below neutron
and muon onset. At fixed h, a changed cutoff moves the endpoint on the same
interior equilibrium. Define alpha=dHs/ds and W=4*pi*R*D/(V*sqrt(D)). Then

    R_s = -D*alpha/(R*V)
    M_s = 4*pi*R^2*e*R_s
    B_s = -W*alpha*nB*1e54
    nu_s,s = (-M_s/R + M*R_s/R^2)/D
    kappa = K_s = nu_s,s + alpha,       K=nu_s+Hs.

nB is in fm^-3. W is in km^3 per unit enthalpy; 1e54 converts km^3 to fm^3.
For a thermodynamically exact EOS, alpha=p/(e+p) and
`kappa=4*pi*p/V`. The large radial term in nu_s,s cancels alpha in K_s;
discarding either term gives an incorrect thermal/chemical normalization.

There is a small literal-source detail. X2 widens its binary64 h0=mp+me
separately from the two binary64 masses. Its explicit k(H) therefore leaves a
nonzero thermodynamic convention defect at high precision. This sidecar
differentiates that exact k(H): alpha=P/P_H and
`kappa=4*pi*(e+p)*alpha/V`. The factor `(e+p)/P_H` is approximately
0.999999999999959 at these surfaces. It is recorded, not silently set to one.
This is differentiation of retained source semantics, not a new EOS correction,
pressure allowance or reopening of C-REF. The closed Fermi integrals and P_H
are independently tested against momentum quadrature and differentiation.

## Observable-specific terminal responses

In this table F_Q means Q_s at **fixed central enthalpy**, before the mandatory
fixed-B adjustment. All symbols have their historical units. C is the Cstar
coefficient before the common temperature/cache factor; J is the corresponding
thermal volume integral. A common factor at the same Tinf knot cancels only
where the historical definition already permits it.

| Observable | F_Q at fixed h |
|---|---|
| R, M, B | Terminal formulas above |
| N_n, N_mu | 0 terminal shell term; central adjustment is nonzero |
| N_e | B_s, since n_e=n_p=nB in the surface shell |
| Ltilde_Me, Ltilde_Mmu | -6*kappa*L; their emitting cores are inside the cutoff |
| mu_B, P0_average | kappa*Q |
| P1_average | kappa*P1; both neutron-weighted integrals have zero surface integrand |
| Cstar | -kappa*C - W*alpha*exp(Hs-K)*cv*cvfactor |
| surface gravity g | g*(M_s/M-2*R_s/R-kappa+alpha) |
| envelope Tlocal | Tlocal*(g_s/(4*g)-0.605*kappa) |
| Tsurface_inf | Tsurface_inf*(g_s/(4*g)+0.395*kappa-alpha) |
| mu_n_profile | Same equilibrium mu_B response at every node, in each original node scale |
| P2_average | Exactly zero, by the retained identity |
| species_support | No change while the accepted topology condition holds |
| metric_structure | No nodewise qualification; retain the accepted consumer restriction |

The two temperature rows are derivatives of the actual nonlinear temperature
outputs at fixed Tinf, not derivatives of a logarithmic temperature state.
Their second and higher responses must still enter a finite-cut estimate.

For chemistry let G0 be the neutron susceptibility integral, G_l the 2x2 lepton
block, U=G_l^-1, and Z=U+(1/G0)*11^T, exactly as in the retained reconstruction.
At the surface chi=chi_p*chi_e/(chi_p+chi_e). Put

    L_s = -W*alpha*exp(Hs-K)*chi*1e54
    (G0)_s = -kappa*G0
    (G_l)_s = -kappa*G_l + L_s*e_e*e_e^T
    F_Z = kappa*Z - L_s*(U*e_e)*(e_e^T*U).

This accounts for both the redshift change of the entire susceptibility matrix
and the rank-one terminal shell contribution. A small shell mass does not bound
the chemical inverse response. The inverse identity is independently tested.
The old full exterior-shell bounds are useful terminal comparisons, but they
are neither fixed-B bounds nor permission to scale a Ps-to-zero tail by 2.1e-12.

## Central adjustment and path tangents

For every scalar/component, with derivatives at fixed s,

    h_s|B = -B_s/B_h
    D_Q := Q_s|B = F_Q - (Q_h/B_h)*B_s.

This is the implicit-function result, including the required central adjustment.
The retained h stencils supply Q_h and B_h. B_h is positive under the existing
C-CONT condition. At the three anchors h_s is about 5.7126e-15 per unit s;
the first-order adjustment at the allowed residual is about 1.20e-26 in h.
It is not zero, and direct subtraction of new ordinary-precision stellar counts
at such a displacement would be a poor test.

For each species put c_n=c_mu=0, c_e=1, t=N_h/B_h, A=B_s, and u=a/B0.
The historical fixed-cut ratio and the certified path tangent are distinct:

    N_s|B = (c-t)*A
    T = t + u*(c-t)*A
    t_h = (N_hh-t*B_hh)/B_h
    t_s|B = (c-t)*A_h/B_h - A*t_h/B_h
    A_s|B = A_s|h - A*A_h/B_h
    D_T = T_s|B = (1-u*A)*t_s|B + u*(c-t)*A_s|B.

A_h and A_s|h follow by differentiating the explicit terminal formula with R,
M and Hs, using retained R_h, M_h. The code keeps the correlated combination
N_hh-t*B_hh; it does not divide unrelated count-error allowances by a tiny span.

There is also an exact on-law alternative. Since the constant-slope derivative
operator commutes with partial_s at fixed B, and
`N_s|B=(c-T)*A/(1-u*A)`,

    D_T = -T_x*A/[B0*(1-u*A)]
          + (c-T)*A_x/[B0*(1-u*A)^2],
    A_x = A_R*R_x + A_M*M_x + A_Hs*alpha*a.

Here all x derivatives are along the accepted law. C-CONT's count curvature
can bound T_x: `|T_x| <= S_N*K_N/B0`. Thus the failed mixed-count stencil
diagnostic does NOT prove that fresh interior solves are necessary. The
on-law identity is checked both against the other formula and by an independent
implicit-root toy problem. It provides an analytic route to a conservative
tangent bound at the law; it does not, by itself, enclose a finite transverse
strip. Differentiating arbitrary sample cut residuals as if they formed a smooth
function is not allowed: a bound on delta_s gives no bound on delta_s'.

## Finite residual, historical drift, fit and normalization

For each row/component the necessary input is a supported

    U_Q(B) >= sup_{|v|<=epsilon} |Q(B,s0(B)+v)-Q(B,s0(B))|/S_Q,
    epsilon=2.1e-12.

Point sensitivities give only `epsilon*|D_Q|/S_Q`. A Taylor realization also
requires a bound/validated estimate for
`epsilon^2*sup|partial_s D_Q|/(2*S_Q)` and transfer across B. Tangents require
the corresponding finite response of T, not merely N. An analytic finite
terminal map can supply these without estimating a noisy third mixed derivative.

For signed relative drifts, q=z+r gives exactly

    d_q-d_z = (r_i-r_0-d_z*r_0)/(z_0+r_0).

Thus baseline error remains. In the existing ledger conventions, add the
worst-case paired numerator `U_i+U_0` (zero at the origin), fit response
`||(I-H)r|| <= f*max U`, baseline U_0, and any supported endpoint contribution
before recomputing the full denominator. Never refit weights or replace the
old historical thresholds. Absolute tangent rows use S=1 and no observable
baseline denominator, but their B_h division and fit still require propagation.

C-TAN already charges fixed-cut/path discrepancy and its accepted ratio,
stencil and geometry terms. C-INPUT charges achieved coordinates and design;
C-CONT charges along-law interpolation/curvature. A new C-CUT term must cover
only the residual transverse cut displacement plus its finite/nonlinear transfer,
not charge those existing contributions a second time. ADR-0022's total
reference fit bound remains unchanged; no shared B-error model is introduced.
