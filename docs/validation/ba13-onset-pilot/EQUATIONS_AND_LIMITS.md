# Independent equation audit and limits

Validation-only, unratified. Production owners unchanged.

## Neutral equilibrium coordinate

Let nB=nn+np, neutrality np=ne+nmu, and beta equilibrium
μn=μp+μe=μp+μmu wherever the respective species is active.
The zero-temperature first law dε=Σ μs dns reduces to dε=μB dnB.
The Euler identity gives ε+P=μB nB; differentiating it gives dP=nB dμB.
Below neutron onset, vary the neutral pe pair: μB=μp+μe, not mn.
The inequality μp+μe<mn means adding a neutron costs energy. Equality is
precisely the onset; it does not justify setting μB=mn throughout the pe crust.

With gtt=−exp(2ν), the ordinary-star TOV equations are
m'=4πr²ε, ν'=(m+4πr³P)/(r(r−2m)), P'=−(ε+P)ν'.
Then dlnμB/dr=−ν', so μB exp(ν) is constant. This applies to the neutral
baryon combination; charged particle potentials individually include electric
potential terms. FR2005 §2.1, equations (1) and (3), establishes the temperature
and heat-capacity conventions; R2006 corrects the electrostatic response treatment.
Primary references: https://arxiv.org/pdf/astro-ph/0502116 and
https://arxiv.org/pdf/astro-ph/0606322. Cluster lacks the Mac-only literature catalog.
The above first-law/TOV derivation is independently given here, not inferred from
a small discrepancy between two profile reconstructions.

Source evidence: `TrackRFreeGasThermodynamics.cpp:221` onset formula;
`TOVSolver.cpp:1516` mass/pressure ODE; `:1544` ν derivative;
`NStar.cpp:841` inward ν reconstruction and surface matching.
The analytic inversion of μp+μe is algebraic in validation `equilibrium.py`,
independent of the review's iterative momentum inversion. Neutron inverse in the
prospective comparator intentionally matches the review; high precision uses the
factorized inverse. The independent structure oracle uses direct momentum
quadrature, enthalpy as the ODE coordinate, and its own common-B shooting.

## Equation consistency is not numerical exactness

The canonical EOS imports independent Steffen splines ε(P) and n(P), so
(ε+P)/n and ∫dP/(ε+P) need not obey the continuum thermodynamic identities between
nodes. NStar integrates ν' on stored radial nodes, introducing another common
error. Anchoring μ∞ at a pe node and imposing equilibrium elsewhere repairs a
reconstruction assumption; it does not prove that this approximate star has the
right physical μ∞, baryon number or metric.

Measure μ∞ from the independently evaluated species energies and enthalpy,
not just the imposed reconstruction. Check pressure onset using the independent
P_on; its log-linear spatial interpolation has its own truncation error.
The same caution applies to the envelope identity Tb/T∞=μb/μ∞. Agreement of
that formula with itself is tautological. The independent analytic EOS solves
εb and Pb; an independently integrated star supplies a separate μ∞ and rb.

Unit source hierarchy: the retained INV-02 text at
`docs/SCIENTIFIC_INVARIANTS.md:252` says ADR-0012 implementation is outstanding.
That statement is stale against the higher-authority accepted ADR-0012 header,
which records complete A1 ratification, and against the live implementation in
`CompactStar/RelativityUnits.hpp:1` / `NStar.cpp:185`. This investigation uses
ADR-0012's matched GSL G=6.673e-8 and c=2.99792458e10 conversion, confirmed in
source, and records the documentary inconsistency without editing either owner
or governance files. No modern-G replacement or mixed solar-mass convention.

## Heat capacity and envelope

Cstar=∫cV(Tlocal)dV. The energy redshift and dTlocal/dT∞ cancel; an extra lapse
must not be appended to cV. For this cold free-gas entropy adapter cV is linear
in Tlocal=T∞exp(−ν), producing the single exp(−ν) factor in Cstar/T∞.
This matches INV-15/ADR-0002, `StarContext.cpp`, and FR2005 equation (3).
The thermal cache linearly interpolates C against log T, so F(T)>1 between
knots for a T-linear capacity. Its common multiplicative factor cancels exactly
from equal-knot ratios. Node quadrature and composition reconstruction do not
necessarily cancel, and are reported separately at all 13 knots.

`TbDefinition.cpp:24` selects the outermost profile-energy node at or above
rho_b; validation therefore uses profile column21, not provider column10.
The independent base satisfies ε/c²=1e10 g/cm³. `EnvelopePotekhin1997.cpp:80`
uses Ts^4=1e24 g14(1.81 Tb/1e8)^2.42: a Tb fractional error transmits with
logarithmic factor 0.605. Local and redshifted surface temperatures share that
base-location bias on a fixed profile. Envelope-model physical uncertainty is
outside this numerical investigation.

## Correlated errors and nonlinear inverses

For Xhat_i=X_i+ε_i, the paired error is ε_i−ε_0. A common additive error
cancels in differences but changes ratio normalization. A common multiplicative
thermal factor cancels in ratios exactly. Neither rule licenses treating other
systematics as independent random noise or replacing them by fit covariance.

In particular, if Ghat=G+E, the exact resolvent is
Ghat^−1−G^−1=−Ghat^−1 E G^−1. The difference of these errors at i and 0
contains differences of inverse backgrounds multiplying common E, even when
E_i=E_0. Therefore the review's suggestion to feed only E_i−E_0 into an
unpaired Neumann majorant is insufficient in general. A paired Z enclosure must
also retain common-background Jacobian variation and inverse nonlinear remainder.
Existing `chemical_trackr_budget.py:zerror` remains unchanged; it is a valid
absolute perturbation majorant with appropriate input enclosures, not an
automatic physical paired-error bound. Retain the inherited chemical diagnostics.

Geometric-series tails U=2q|d|/(1−q) require continued contraction beyond the
last computed level. The finite pilot cannot establish that hypothesis globally,
nor can it rule out a shared deterministic bias. All reported such U are
conditional estimator values, not certificate-ready physical bounds. The
inter-sample/inner-half response remainder and full tangent/background budgets
remain separate obligations. Historical fitted covariance is not added to them.

## Rounding interpretation and limited audit coverage

The original depth13 grid's four refused pairs have distinct exact physical
pressures/energies, also distinct after correct rounding to binary64. This does
not contradict the review's more detailed §6.2: its “representational saturation”
includes the cube-root and equilibrium intermediates. It must not be read as a
proof that the output values themselves are intrinsically unrepresentable.
The capped construction changes validation sampling only. It neither repairs nor
reclassifies the historical 32768 refusal.

This pilot verifies every adjacent pair in each generated table in the provider
and serialized cgs arithmetic. Its high-precision oracle covers the four historical
refused pairs plus nearest pairs on both sides of the neutron onset in each capped
table. It is not the review's suggested every-node extended-precision error audit,
and supplies no maximum provider-ULP error over the entire 32768c table. Nor is
it an interval proof of the review's piecewise Gamma lower bound. Those are
remaining bounded evidence if an absolute EOS enclosure is required. The pilot
uses the explicitly predeclared N=32768 depth11/depth12 probe; the review's §6.3
N=8192 depth12 variant is not represented as having been run.
