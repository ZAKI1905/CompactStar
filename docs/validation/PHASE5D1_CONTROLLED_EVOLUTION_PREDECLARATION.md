# Phase-5D-1 controlled evolution predeclaration

Status: PRE-RESULT / IMMUTABLE BENCHMARK INPUTS. Scientific-semantic and architectural candidate under accepted ADR-0014 §§3–5, its preflight §§25,31–33, and GOVERNANCE §4. No governed baseline is installed. This record is committed before any production coupled trajectory is generated.

Canonical entry: `d019ae390be4f5e3daba05039903485cb497e397`.
Branch: `physics/phase5d-controlled-rotochemical-evolution`.
Worktree: `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase5d-controlled-evolution`.
Sole writer: lead / Agent D. Read-only specialists A/B/C: source, architecture, validation, all Astra HIGH. Their audits found no scientific contradiction with ADR-0014. Required upstream authority is ADR-0010/0011/0013; global INV-11 remains unresolved.

## Fixture and frozen normalizations

Use fresh canonical ordinary NStar Track-R Structure-1, central mass-energy density `1.10e15 g cm^-3`, baseline radial resolution 10000, fresh matched free-gas table and the governed Phase-5C production certificates/acceptance goals. Construct Z and W through production semantic objects; never read baseline JSON or copy coefficient numbers into evolution source.

Typed static process selection enables Me/Mmu and disables De/Dmu. D is exactly the associated G_y whole-star chemical domain, with the inherited center/onset/refusal/tail semantics. D_Me is the part with neutrons, protons and electrons present, ending at the neutron onset; D_Mmu ends at muon onset. Continuous endpoints have measure zero. The pe tail contains no enabled MU support. No density cutoff or temperature-dependent support is introduced. Generic support is an ordered union of disjoint intervals, not a last-index mask.

Local normalizations, fixed before results:

| Process | S(r) on its declared support | Units | Classification |
|---|---:|---|---|
| Me | `1.0e-51` | erg cm^-3 s^-1 K^-8 | declared mathematical/architecture input |
| Mmu | `2.0e-51` | erg cm^-3 s^-1 K^-8 | declared mathematical/architecture input |
| De, Dmu | disabled | erg cm^-3 s^-1 K^-6 | no benchmark contribution |

These constants are positive, lepton-separated aggregate modified-Urca normalizations. They are not FR2005/Yakovlev realistic rates or a neutron/proton branch model. Their order of magnitude gives an ordinary stellar coefficient near `10^-32 erg s^-1 K^-8`; that dimensional estimate is not a result target. No output-driven coefficient tuning is permitted.

Compute each Ltilde with independent production machinery, `10^15 integral_Da 4*pi*r_km^2 exp(lambda) S exp((2-q)*nu) dr_km`. The exact `10^15` length-volume conversion is owned by Units::KM3_TO_CM3; there is no particle-count factor `10^54`. Use the same immutable coefficient object for equilibrium cooling, reaction response, nonequilibrium neutrinos, and heating. Historical placeholder equilibrium coefficients are excluded from this controlled configuration, without altering default passive cooling.

## Thermal adapter and spin history

`Tinf(0)=1e8 K`; `(eta_npe,eta_npmu)(0)=(0,0) MeV`; state `(ln(Tinf/1e8 K),eta_npe,eta_npmu)`.

Use the existing iron Potekhin1997 envelope/photon cooling with no equation changes. Reuse StarContext::HeatCapacityStar_Tinf. Supply a generated scratch CompOSE-format **fixed-background** entropy adapter, `Q2=a(nB)*T_MeV`, with `a=sum_i pF_i*mu_i/[3*(hbar*c)^3*nB]`, using total free-gas potentials and g=2 already included. FR2005 §3.4 eq.(50) supplies the degenerate coefficient; values outside degeneracy are a declared mathematical extension, not finite-temperature EOS authority. Species absent at a node contribute zero. Use positive profile density nodes (no division at vacuum), two identical Yq planes `{0,1}`, and T planes `{0,1,2,4} MeV` covering all cache queries. No off-equilibrium composition derivative is claimed. Own an actually const thermal table with no mutable alias and preserve input bytes. The existing 160-point log-temperature cache and interpolation are retained; reject the coupled benchmark if any accepted sample reaches its `Tinf_MeV<=1e-5` floor or `>=1` ceiling. No claim of exact C proportional T between cache knots.

External prescribed dipole history (preflight §25.3): `B=1e8 G`, `P0=1e-3 s`, `PPdot=(B/3.2e19)^2`, `P(t)=sqrt(P0^2+2*PPdot*t)`, `Omega=2*pi/P`, `OmegaDot=-2*pi*PPdot/P^3`. Record and test the derivative. The rotochemical module consumes the history and owns no torque law. Run `t=0...1e10 yr`, `1 yr=365.25*86400 s`. A generic state-coupled adapter must share the same external torque evaluator with its spin driver.

## Numerical and pass/fail contract

Retain GSL RKF45. Baseline rtol `1e-7`; absolute tolerances `(1e-12,1e-18,1e-18)` in the above flat order; refinement rtol `1e-9`, absolute tolerances divided by 100. Generic per-component tolerances may be added; empty vector preserves legacy scalar behavior. Output initial state plus 40 samples per decade from 1 yr through 1e10 yr and exact endpoint. Maximum internal steps per output 100000; sample cap 1000. No eta clamp.

Global integration: Gauss-Legendre order 16 per explicit profile/support segment; refine to 32 and split every segment in two. Relative coefficient difference must be <=`1e-6`. This is numerical stability of the declared representation, not a continuum error bound. Reuse the governed Phase-5C background characterization and separately compare a radial-resolution 20000 Ltilde variant against baseline to <=`5e-3` relative.

ODE comparison at common output times: max relative T difference <=`2e-4`; each eta difference divided by `max(abs(eta_ref),1e-10 MeV)` <=`2e-4`. Report actual maxima; no tolerance relaxation. All values must be finite, T positive, both eta positive after spin-down begins, R of same sign as eta, and LH/DeltaLnu nonnegative. Initial passive cooling, imbalance growth, positive chemical power, enhanced neutrinos, and both negative and positive incremental beta power must be observed. Exact eta=0 incremental zeros and same-owner equilibrium reduction are separate analytic tests.

Quasi-steady predeclared fit window: `1e9...1e10 yr`; use samples only with both abs(xi)>=100 (print eligible count; insufficient eligibility is a failure, not permission to change window). Require >=10 samples. Compare each eta to `[2*kB_erg*(kB_MeV)^7*I_l*Omega*OmegaDot/(C_H*Ltilde_l)]^(1/7)`, `C_H=24/(11513*pi^8)`: max relative discrepancy <=`0.05`; log eta versus log abs(Omega*OmegaDot) slope within `0.015` of `1/7`. Report photon power versus large-xi `5/8 sum eta R` and thermal residual distinctly. Do not claim asymptotic behavior outside the eligible regime. Initial-condition variants `T0={1e7,1e9} K` and `eta0={kBT0,20*kBT0}` with T0=1e8 for chemical variants must approach endpoint within 1% of baseline; these are validation variants, not replacement initial conditions.

Stiffness diagnostic: report RHS evaluation count and local numerical Jacobian relaxation times/eigenvalue ratio at preselected `1e6,1e8,1e10 yr`; require successful explicit integration under the step budget. If major solver redesign is needed, stop.

## Independent RE10b and analytic gates

RE10b has **relative tolerance 1e-10**, predeclared here. Synthetic radius `x=r/km`, D=`[0,1]`, `nu=-0.4+0.1*x^2`, `lambda=0.2*x^2`, `S=S0*(2+x^2)`, `S0=1e-40` in the appropriate q units. Test q=6 and q=8; D_a=`[0.2,0.45] union [0.7,0.9]`. S is nonzero outside D_a, including the closed inner region. The independent expected value is `4*pi*1e15*S0*exp(0.4*(q-2))*integral_Da x^2*(2+x^2)*exp((0.4-0.1*q)*x^2) dx`, evaluated by closed-form/power-series integral and independent high-precision quadrature, never production. Metric inputs are analytic callbacks; assert ordered radial partitions.

Kill missing/extra/wrong-sign lapse, G_y inverse-lapse substitution, missing/inverted proper volume, wrong q, wrong domain and last-index sweep. Separately refuse chemical-domain identity mismatch. Triangle-open DU with positive input outside current nB_min must remain disabled by selection; test independently of the historical mask.

RE1–RE14 and RE16–RE18 use the accepted numbering (RE10 spin-only, RE12 two-temperature coupled linear oracle). Source polynomial tests include independent Fermi-convolution integrals, derivative identity F'=3H, parity, small/large limits and all four full/incremental roots. Linear small-eta comparisons use xi<=1e-4, relative tolerance 1e-6, two distinct temperatures. Spin-only analytic comparison tolerance 1e-6 relative with eta0 retained. Lyapunov nonincrease includes a dead channel and Z cross-coupling. M6 is the transpose of the complete semantic linear map Z*diag(rate slopes); a literal transpose of symmetric Z is algebraically identical and receives no mutation credit. M10/M11 is one integrated-lapse mutation. All mutation results identify actual production or transformed-input path and aliases; no inflated total.

## Frozen boundaries

Permitted refinements: solver tolerances, quadrature order/subdivision, radial resolution as explicitly above; repair implementation defects against these frozen oracles. Forbidden: changing normalization, support, initial conditions, spin law, interval, fit window, pass criteria, source authority, historical baseline or EOS/literature bytes after results. A pre-result implementation impossibility stops for an explicit record; no silent amendment. No realistic A18/FR2005 absolute-normalization claim and no BNV. Global INV-11 remains unresolved pending independent review and owner ratification.
