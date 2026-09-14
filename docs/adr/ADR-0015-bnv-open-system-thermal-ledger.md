# ADR-0015: Moving-equilibrium BNV thermal first-law contract

## 1. Status, scope and non-scope

**PROPOSED — NOT ACCEPTED — NOT OWNER-RATIFIED — NOT CANONICALLY INTEGRATED.**
Date:2026-09-13. Phase-6A-0: **SCIENTIFICALLY HARDENED DRAFT / READY FOR FINAL INDEPENDENT REVIEW.**
Canonical master `0a7418aecb7314cfa472a78f1faf477be8456a94`; draft entry
`5a6bf7cb9455d684ddb6fccb22ad2b9fec940b3a`. This is a docs-only scientific contract proposal,
not BNV production implementation, physical-rate selection, trajectory, EOS or A18 work.

The first controlled scope is spin OFF, whole-star diffusive non-superfluid free gas, cold
chemical response, frozen standard coefficients and an abstract charge-consistent source.
Phase-5B/C/D remain unchanged. Phase-5C is correct for its declared R2006/Cowling contract;
this ADR does not supersede it or its baseline. Global INV-11 and realistic A18 retain their
separate limits (`docs/SCIENTIFIC_INVARIANTS.md:1004`).

## 2. Evidence hierarchy

GOVERNANCE (`GOVERNANCE.md:14`) and accepted ADR-0011/0013/0014 remain project authority.
Catalog role defines published authority. Incoming local PDFs are supporting material, not
catalog-promoted sources. Author notes are unpublished/non-authoritative. The initial Opus
consultation is working history; the owner-supplied five-agent consolidated hardening is
independent evidence, not a governed contract. This proposal derives and checks its claims.

The [R1 preflight](../validation/PHASE6A0_BNV_THERMAL_FIRST_LAW_PREFLIGHT.md) contains exact
hashes, G1-G20 derivations, published and author-note audits, R0 corrections, scope restrictions
and validation. The [Cowling diagnostic](../validation/PHASE6A0_COWLING_BARYON_DIRECTION_DIAGNOSTIC.md)
is non-governed scientific evidence. Published omissions never establish the author's intent;
the supplied notes do not necessarily exhaust the author's intended theory.

## 3. Species, source and energy conventions

Metric(-,+,+,+), c=1 in contractions; rho and mu include rest energy with one consistent zero.
D=u.grad is proper-time derivative; dots are infinity-coordinate time. Gamma_i is positive
creation per proper volume/time; q_E=-u.Q is total energy INTO the ordinary thermal fluid.
Local power density is erg cm^-3 proper-s^-1, global power erg/s.

Use governed y=(N_n,N_e,N_mu), N_p=N_e+N_mu, g_y=(mu_n,mu_p+mu_e,mu_p+mu_mu),
b=(1,1,1), L=[[-1,-1],[1,0],[0,1]], P selecting e/mu. Eta=-L^Tg_infinity in MeV, R positive
for neutron decay/lepton creation in count/s. Z has MeV/count, t count/count, sigma count/s.
Convert chemical MeV to erg exactly once using the existing C_M; equations below suppress C_M
when all energy terms use common units. Basis/code authority:
`CompactStar/Analysis/ChemicalResponse.hpp:16`, `CompactStar/Analysis/src/ChemicalResponse.cpp:705`.

## 4. Local open-system first law

The proposed law is derived by contracting the perfect-fluid stress equation and applying
Gibbs/Euler (preflight R1-R4):

```text
nabla_mu(n_i u^mu)=Gamma_i
D rho+(rho+p)theta=q_E
T nabla_mu(su^mu)=q_E-sum_i mu_i Gamma_i
T n_B D s_b=q_E-sum_i mu_i Gamma_i-T s_b Gamma_B
c_V DT=q_E-sum_i[mu_i+T(partial s/partial n_i)_T]Gamma_i
       -T[s-sum_i n_i(partial s/partial n_i)_T]theta.          (A1)
```

Compression work cancels from the entropy equation. Reversible PdV, gravitational/hydrostatic
readjustment and binding-energy changes are not extra heat. The s_b denominator term is not
an independent flow. Exact finite-T temperature evolution differs from the entropy-energy
residual by composition/adiabatic terms; representative free-gas corrections are eV/event
at 10^8K, versus MeV direct energy. Controlled frozen-v1 omits them with a declared sign floor.

For static/quasistatic slices, source counts use `int e^Phi Gamma dV`, power uses
`int e^(2Phi)q dV`. One lapse is time dilation; the second is energy redshift. Thorne's
Gamma_B=0 thermal balance and FR05/ADR-0014 beta energy balance follow explicitly in preflight
section 4. Arbitrary rotating metrics require the corresponding Killing-energy/current treatment.

## 5. Direct-event energy ledger

Each directed channel declares Gamma_i=nu_ia R_a, R_a>=0. E_esc,fluid means energy leaving the
**ordinary thermal fluid**, not necessarily the star. Without external input,

```text
q_dir,a=[-sum_i nu_ia mu_i-E_esc,fluid,a]R_a
q_dir,neutron=(mu_n,actual-E_esc,fluid)R.                      (A2)
```

External input, if independently supplied, is a distinct owner. Split fluid exit into actual
star escape and retained-X transfer; later X release is counted once. For a neutron occupied
state E_n=E_dep+E_esc,fluid, `mu_n-E_esc,fluid=(mu_n-E_n)+E_dep`.
**Fermi-hole heating is not an additional independent term.** Supporting MPR21 three-piece
collision decomposition reduces to the same residual (preflight section 5).

Exactly one representation: R-a {chemical term,total matter-energy transfer}; R-b {partial
hole, independently defined deposition}; R-c {complete mu-E_esc}. No hole or deposition addend
on top of an inclusive residual. n->chi gamma is illustrative only: absorbed photon gives
mu_n-E_chi; both products leave gives mu_n-E_n; partial deposition is accounted by actual
remaining fluid-exit energy. No rate or RMF kinematics is ratified.

For the declared cold spontaneous neutron-removal conditions, E_n<=mu_n, E_dep>=0 imply
0<=Q_dir and the rate-weighted hole lower bound. Actual mu supplies the general upper budget;
`P_dir<=|Eeq_dot|-L_esc,rest` additionally requires equilibrium or the demonstrated neutron-sink
sign delta mu_n<=0. Here L_esc,rest=sum_a int e^Phi R_a(sum_truly_escaped m_j c^2)dV
uses one lapse for event counts and the minimum Killing energy at infinity; it is not a
two-lapse integral of local rest mass. Uniform nonrelativistic occupation gives2E_F,kin/5; no universal positive
lower bound exists for an arbitrarily Fermi-surface-tuned process. Invisible is not generally
zero heat. The full conditions and relativistic average are preflight R9-R10.

## 6. Moving-equilibrium chemical response

The primary chemical variables and proposed evolution are

```text
delta N_y=N_y-N_y^eq(B,Omega)=L ell,  eta=-Z ell
S_y=int e^Phi Gamma_y^BNV dV,  Bdot=b^T S_y
 t=(partial N_y^eq/partial B)_Omega,  b^Tt=1
sigma=P(S_y-t Bdot),  S_y-t Bdot=L sigma
ell_dot=R+sigma-2I_Omega Omega dotOmega
eta_dot=-Z(R+sigma)+2W Omega dotOmega+Zdot Z^-1eta,
W=Z I_Omega.                                                (A3)
```

I_Omega here has two lepton components; L I_Omega is the three-axis target spin tangent.
Zdot=0 in the first toy. The future positive Zdot sign is already governed by
`docs/adr/ADR-0014-secular-rotochemical-evolution-contract.md:380`.
S=0 recovers Phase-5D. S=t Bdot, zero initial eta and no spin forcing gives eta=0 identically:
this is the **physical sliding null**. Baryon-neutrality is derived from the moving target,
not imposed by projecting raw baryon loss through a fixed-metric inverse.

## 7. t, k, G_y and Z ownership

At spin OFF, t_i=B_i/B_B uses the existing Phase-5B structural derivative owner, with identical
neighboring-star/domain/surface policy, denominator qualification, errors and currency.
It is a newly derived proposed view, not a previously governed BNV object.
Required b^Tt=1 and t_p=t_e+t_mu; t approximately(0.965770,0.0308522,0.00337774).
The propagated numerical budgets and chart extensions are preflight section 6.

Cowling k=G_yb/(b^TG_yb) approximately(0.992218,0.00748454,0.000297569) is diagnostic only;
t_e/k_e approximately 4.12, t_mu/k_mu approximately 11.35. G_y builds governed Z and diagnostic k;
it MUST NOT map raw BNV S to physical eta_dot or reconstruct actual changing-B individual mu.
The rejected route equals -Z(S_l-k_l Bdot), not A3. Its old -3.442789339881122e-56 MeV/count
both-channel value is a **COWLING k-ROUTE NEGATIVE ORACLE**.

Physical neutron-only S=(Bdot,0,0), Bdot<0 gives sigma=t_l|Bdot| and
eta_dot/|Bdot| approximately(-1.4303e-55,-3.6281e-55) MeV/count. Both eta<0 imply capture R<0,
eta_lR_l>=0. The rejected route spuriously drives a physical slide by
(+1.086e-55,+3.284e-55)|Bdot|. Z+t suffice after moving-reference subtraction.
True-response findings remain non-governed diagnostic evidence; no Phase-5 coefficient changes.

## 8. Chemical free-energy reservoir and individual potentials

At fixed current B, `E_chem=ell^TZell/2=eta^TZ^-1eta/2>=0`, a state reservoir, not heat.
For spin-off first-order departures,

```text
delta g^infinity=-(I-bt^T)P^Teta+O(eta^2)
delta mu_n^infinity=eta^Tt_l+O(eta^2).                        (A4)
```

These follow from differentiating E_eq(B)+ell^TZell/2 at the moving reference. In the frozen
quadratic model they are exact; varying Z(B) adds b ell^TZ_,B ell/2 at quadratic order.

The correct derivative, including future coefficient variation, has two equivalent forms:

```text
Echem_dot=eta^TZ^-1eta_dot-(1/2)eta^TZ^-1Zdot Z^-1eta
        =-eta^T(R+sigma)+2Omega dotOmega eta^T I_Omega
         +(1/2)eta^TZ^-1Zdot Z^-1eta.                        (A5)
```

The hardening report's negative half **after substitution** is corrected here: A3 supplies
another positive full term. The negative half is valid in the first chain-rule line only.
Preflight R17 includes independent finite-difference counterexamples. Frozen coefficients
remove both metric terms and preserve the requested spin-filling sign. Old fixed-reference
alpha/a/alpha^2/(2a), E_2 and free F_ref constructions are rejected as primary physical storage.

## 9. Controlled thermal equation and global closure

Using actual potentials, P_dir=-g_actual^TS-L_out,fluid. The equilibrium decomposition satisfies

```text
P_dir(actual)=P_dir(eq)+eta^Tsigma,
P_dir(eq)=-mu_B^infinity Bdot-L_out,fluid
C_* T_infinity_dot=P_dir(actual)+eta^TR-Lnu_eq-DeltaLnu-Lgamma-Lother.
                                                               (A6)
```

For the neutron sink the plus sign follows directly from mu_n,actual=mu_B+eta^Tt_l.
At spin OFF/frozen coefficients, A5+A6 give
`Eeq_dot+Echem_dot+Udot=-L_out,fluid-Lnu_full-Lgamma-Lother`.
Add an explicit X reservoir to convert fluid-exit accounting to star-escape accounting.
There is no free reference-work addend, no dot E_chem heat, no generic PdV heat.
Changing coefficients/thermal structure or imposed spin requires its real state derivatives/
Omega dot J work; the spin-off identity cannot silently claim all-orders rotating ADM closure.

The static neutral stellar first law is independently derived by the TOV variation kernel in
preflight R23-R24; the cold limit gives d(M_eq c^2)/dB=mu_B^infinity with surface-work qualification.
Uniform-rotation extension is conditional on the stationary variational problem and proper
Killing conjugates. General quotable theorem authority remains source-limited, not a blocker
for the explicitly derived spin-off contract.

## 10. Beta-mediated thermal response

Retain governed `L_H=C_M sum eta_lR_l`, `DeltaLnu=sum Ltilde_a(F_a-1)T^q`,
`DeltaP_beta=L_H-DeltaLnu` with the same coefficients and signed rates. This is the actual-
potential decomposition. Relative to P_dir(eq), spin-off/frozen correction is
`-Echem_dot-DeltaLnu`, nonpositive while storage fills, exactly -DeltaLnu in reached QSS.
Do not confuse these decompositions or assert storage always fills.

For0<|xi|<4.9097100289 (M) or4.7870134733 (D), the incremental process power is negative.
Instantaneous cooling magnitudes are bounded by0.467659 Lnu_eq,M and0.527775 Lnu_eq,D;
DeltaP_beta/L_H lies in[-1/2,5/8] (pure D upper1/2). Preflight R21 proves these from the
polynomials. The approximately 0.724 keV/baryon source-driven cooling scale at 10^8K additionally
requires the capture/filling history assumptions stated there. It is not a bound on arbitrary
initial stored imbalance. Ordinary10-30 MeV direct sources dominate that scale.

QSS R=-sigma is Z-independent but requires reachability. R0 estimates3.4e12/9e9/3e6 yr at
abstract1e-17/1e-14/1e-10 yr^-1; these are attributed examples, not chosen physical rates.
Use eta=eta0-Z int sigma dt when reactions are negligible. QSS and thermal asymptotic formulas
must not be claimed for an unreached state.

## 11. Generic versus process-dependent inputs

Class A: total baryon loss plus sequence; B: ordinary stoichiometry/source; C: energy partition;
D: product fate; E: hidden-sector interactions/accumulated state. Sequence and chemistry are
generic conditional on A+B. Absolute thermal predictions require C-E. A single rate API must
not encode a hidden energy efficiency. Event identity links particle and energy owners.

## 12. Regime-I validity

Every terminal species must satisfy negligible production blocking/accumulation feedback,
stress/EOS and chemical-equilibrium influence, heat capacity, radiation/conduction/opacity,
and declared coefficient-drift budgets. N_X/B<<1 alone is insufficient. Near zero source,
eta or direct power use absolute tolerances, not impossible relative inequalities.
Frozen scope additionally bounds variations of Z,W,Ltilde,C_*,t,surface gravity,I and support.
No universal tolerance is chosen. Explicit testable criteria are preflight section 12.

## 13. Regime-II boundary

Mandatory fate flags: PROMPT_ESCAPE, SM_THERMALIZATION, BOUND_INERT, BOUND_INTERACTING,
including mixed branches and terminal products. Failure of X thermal/mechanical/chemical
conditions triggers the corresponding Regime II. Pure ordinary coefficient drift exits the
frozen contract but can lead to evolving ordinary Regime I, not necessarily hidden-sector II.
Conditional linear accumulation gives a rate-independent L_X/P_dir monitor; it is not universal.

Proposed progression: ordinary star -> MixedStar-lite test fluid with N_X/T_X ledger -> full
two-fluid state(B,B_X,J,S,S_X). Slaved T_X requires rapid relaxation and P_exchange=L_X,rad.
Goldman supports qualitative transfer through a mechanically dilute accumulated sector, with
quantitative source limitations. Cooling relative to Regime-I prediction differs from cooling
relative to passive control; constant positive deposition and a loss vanishing at zero T imply
a positive temperature floor. No Regime-II implementation or rate is authorized.

## 14. Matched control and sign classification

Match initial star/T/eta, standard microphysics, photon/neutrino models, spin, solver and
tolerances; remove only BNV source/product sector. First toy holds corresponding frozen B0.
Report DeltaT, DeltaT_s, DeltaLgamma, integrated direct/incremental powers and DeltaU;
`d DeltaU/dt=DeltaP_total` includes changed equilibrium neutrinos and photon feedback.

Report NET HEATING/NET COOLING/NEAR ZERO with SIGN UNRESOLVED when error spans zero. R0's
provisional1e-4 relative temperature floor must be remeasured for future runs and combined with
eV/event finite-T, partition, coefficient and t/Z errors. Thermal lag is C_*T/(a_cool L).
Generic Regime-I cooling cannot be claimed; a tuned small-direct-energy corner and applicable
history are required. No positive heating efficiency or arbitrary-history sign theorem.

## 15. First controlled toy proposal

Spin OFF, whole-star governed free gas, uniform abstract proper neutron sink; t from qualified
Phase-5B inputs, governed Z unchanged. P0: E_esc=actual current mu_n, artificial cold zero-direct
entropy oracle, not physical n->chi gamma; P1: E_esc=E_n, uniform-sea hole diagnostic; P2:
E_esc=0 maximal retention. P1gamma kinematics may follow only after source clarification.

R0 two-tier design examples are mathematical fixtures only: P1/P2 around 1e-17 yr^-1, P0 around
1e-13-1e-12 yr^-1 for 1e7-1e8 yr or transient analytic tests. This task selects no physical
BNV rate, implements no rate/trajectory and chooses no microscopic constraints/couplings.

## 16. Required validation oracles

Zero-source Phase-5D bit identity; physical sliding null; k-route refusal; t sum rules and
currency; physical neutron/capture signs; transient source-only solution; reached QSS balance;
actual/eq and R-a/b/c energy equality; chemical state never heat; variable-Z chain-rule check;
beta roots/B1/conditional B2; P0 actual-potential and P1 sea average; charge/fate closure;
controlled finite-interval conservation; matched-control sign/lag floors. G1-G20 proof/status
checklist is in preflight section 14. Tests are proposed, not implemented.

## 17. Source-limited items and author notes

The preflight preserves the four-paper published audit and separate Thermo_BNV/whiteboard
classifications. It includes MPR Eq 8 possible prefactor typo; Goldman quantitative transport/
normalization limitations; A24 AppendixE ambiguities; new author-note dimensional defects;
and corrections to unsupported R0 Goldman-volume claims. No source or note is edited/promoted.

Missing primary general first-law/transport/RMF sources are not claimed read. Static derived
closure stands; the general rotating theorem, realistic-EOS response, detailed escape/rates,
Regime-II transport, superfluidity and complete G_true numeric reproduction remain future
source-qualified work. A note's title/omission is not author intent or proof of false physics.

## 18. Consequences and implementation ownership

Candidate BNV-14 through BNV-25 in preflight section 15 are part of this proposal, not ratified
invariants. They cover the cold bracket, moving reference/t, diagnostic-only raw-G route,
individual-potential/state-energy identity, once-only event representation, product regimes,
beta bounds/QSS, actual-mu P0, Cowling scope and matched-control floors.

Separate future owners: process/charge source, energy partition, fate/X state, global integration,
structural t, governed Z/sigma, chemical storage, standard beta, thermal ledger and matched
comparison. Reuse generic Phase-5D machinery; no reaction-specific generic-driver branches.
No production source, tests, coefficients, baselines, EOS or data change under this ADR.

## 19. Rejected alternatives and double-count refusals

Rejected primary physics: raw D_BNV S; physical null k Bdot; old equal-channel sink oracle;
fixed-reference alpha^2/(2a)/E_2; free F_ref/Delta(delta g^T F) as physical thermal flow;
separate hole over complete mu-E_esc; dot E_chem or dot M_eff as heat; generic PdV/gravity heat;
duplicate weak bulk dissipation/radiation; bound-product positive heat counted twice;
P0 based on equilibrium mu in a disequilibrated star; generic Regime-I cooling; QSS without
reachability. Historical occurrences are negative tests only. The exact H1-H10 refusals are
preflight section 18, with boundary/energy-zero qualifications.

## 20. Ratification requirements

This is a scientifically hardened **proposal**, including explicit R1 corrections to the
independent report. Obtain a fresh-context final independent Opus scientific review of the
preflight, this ADR and the diagnostic, then explicit owner ratification before implementation.
The exact requested next-action text is in preflight section 19. Do not initiate it automatically.
No canonical merge, owner ratification, BNV rate, production trajectory or A18 work is authorized.
