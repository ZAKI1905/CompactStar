# ADR-0015: Moving-equilibrium BNV thermal first-law contract

## 1. Status, scope and non-scope

**Status:** ACCEPTED / HUMAN-RATIFIED.
Date accepted: 2026-09-18. Phase-6A-0: **PREFLIGHT COMPLETE / INDEPENDENTLY REVIEWED /
HUMAN-RATIFIED. CONTRACT ACCEPTED; PRODUCTION IMPLEMENTATION NOT BEGUN.**
Canonical master `0a7418aecb7314cfa472a78f1faf477be8456a94`; draft entry
`5a6bf7cb9455d684ddb6fccb22ad2b9fec940b3a`. R3 entry / PHASE6A0_R2_REVIEWED_SHA:
`7a31862f1e8a5cf046316882e05d76e4924e27d9`; final correction / reviewed branch SHA
`58b375631d6948ed254809567bd29050f9735089`. This is a docs-only accepted scientific contract,
not BNV production implementation, physical-rate selection, trajectory, EOS or A18 work.

The first controlled scope is spin OFF, whole-star diffusive non-superfluid free gas, cold
chemical response, frozen standard coefficients and an abstract charge-consistent source.
Phase-5B/C/D remain unchanged. Phase-5C is correct for its declared R2006/Cowling contract;
this ADR does not supersede ADR-0013/0014 coefficient mathematics, Q/Z/W ownership, governed
baselines, the Cowling source contract inside the declared fixed-baryon scope, or standard
rotochemical physics. It DOES narrow their forward-looking BNV-seam designation as explicitly
owner-ratified in section 7. Global INV-11 and realistic A18 retain their
separate limits (`docs/SCIENTIFIC_INVARIANTS.md:1004`).

## 2. Evidence hierarchy

GOVERNANCE (`GOVERNANCE.md:14`) and accepted ADR-0011/0013/0014 remain project authority.
Catalog role defines published authority. Incoming local PDFs are supporting material, not
catalog-promoted sources. Author notes are unpublished/non-authoritative. The initial Opus
consultation is working history; the owner-supplied five-agent consolidated hardening is
independent evidence, not a governed contract. This proposal derives and checks its claims.

The [preflight with R3 corrections](../validation/PHASE6A0_BNV_THERMAL_FIRST_LAW_PREFLIGHT.md) contains exact
hashes, G1-G20 derivations, published and author-note audits, R0 corrections, scope restrictions
and validation. The [Cowling diagnostic](../validation/PHASE6A0_COWLING_BARYON_DIRECTION_DIAGNOSTIC.md)
has status **DERIVED / NUMERICALLY SUPPORTED / NOT RATIFIED AS PRODUCTION PHYSICS**.
The R2 final report and Reviewer E were read completely; provenance and the E-1–E-27 correction
record are in preflight section 21. Published omissions never establish the author's intent;
the supplied notes do not necessarily exhaust the author's intended theory.

## 3. Species, source and energy conventions

Metric(-,+,+,+), c=1 in contractions; rho and mu include rest energy with one consistent zero.
D=u.grad is proper-time derivative; dots are infinity-coordinate time. Gamma_i is positive
creation per proper volume/time; q_E=-u.Q is total energy INTO the ordinary thermal fluid.
Local power density is erg cm^-3 proper-s^-1, global power erg/s.

Use governed y=(N_n,N_e,N_mu), N_p=N_e+N_mu, g_y=(mu_n,mu_p+mu_e,mu_p+mu_mu),
b=(1,1,1), L=[[-1,-1],[1,0],[0,1]], P selecting e/mu. Eta=-L^Tg_infinity in MeV, R positive
for neutron decay/lepton creation in count/s. Z has MeV/count, t count/count, sigma count/s.
Reserve R for the global beta-rate vector; R_a is local directed BNV event-rate density in
count cm^-3 proper-s^-1. C_(MeV->erg) denotes the existing `MeVToErg` conversion, owned by
`CompactStar/Physics/Rotochemical/ChemicalImbalanceState.hpp:14`; the exactly-once thermal
boundary is `CompactStar/Physics/Rotochemical/RotochemicalReactionResponse.hpp:48` and `:53`.
A1 uses rho, mu and q in consistent erg units; A2 and the chemical algebra use MeV potentials,
with conversion explicit at thermal boundaries. E_chem in A5 is MeV; C_(MeV->erg)E_chem is erg.
All global luminosities/powers are at infinity, including when the superscript is suppressed.
Basis/code authority:
`CompactStar/Analysis/ChemicalResponse.hpp:16`, `CompactStar/Analysis/src/ChemicalResponse.cpp:705`.

## 4. Local open-system first law

The accepted law is derived by contracting the perfect-fluid stress equation and applying
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
residual by composition/adiabatic terms. The finite-T floor is CHANNEL-WEIGHTING DEPENDENT:
smooth P0 gives (pi^2/6)(k_BT)^2/E_F,kin and uniform-sea P1 gives
(pi^2/3)(k_BT)^2/E_F,kin, approximately 1.2 and 2.4 eV/event at 10^8 K and 100 MeV.
A source concentrated within k_BT of the Fermi surface can instead have an O(k_BT) residual
of either sign, about 8.6 keV/event at 10^8 K. Every physical channel MUST redeclare its own
finite-T weighting floor. The compression coefficient (2/3)Ts is the NR-species example;
ultrarelativistic species have a different coefficient. P0's zero entropy residual leaves
+T s_,n R_a in the fixed-volume temperature equation (preflight R4).

Controlled frozen-v1 uses Udot=C_* T_infinity_dot at fixed B and the cold tangent
Eeq_dot=C_(MeV->erg)mu_B^infinity Bdot in erg/s. It explicitly omits finite-T B-dependent
state derivatives, including (partial U_th/partial B)_T Bdot and the finite-T correction to
mu_B. These belong to the same small finite-T class for the smooth toy; a future exact
physical conservation oracle must own them. The current P0 tier is a chemistry/transient
oracle, not a resolved thermal-sign experiment (section 15).

For static/quasistatic slices, source counts use `int e^Phi Gamma dV`, power uses
`int e^(2Phi)q dV`. One lapse is time dilation; the second is energy redshift. Thorne's
Gamma_B=0 thermal balance and FR05/ADR-0014 beta energy balance follow explicitly in preflight
section 4. Arbitrary rotating metrics require the corresponding Killing-energy/current treatment.

## 5. Direct-event energy ledger

Each directed channel declares Gamma_i=nu_ia R_a, R_a>=0. E_esc,fluid means energy leaving the
**ordinary thermal fluid**, not necessarily the star. Without external input,

```text
q_dir,a=C_(MeV->erg)[-sum_i nu_ia mu_i-E_esc,fluid,a]R_a
q_dir,neutron=C_(MeV->erg)(mu_n,actual-E_esc,fluid)R_a.                      (A2)
```

External input, if independently supplied, is a distinct owner. Split fluid exit into actual
star escape and retained-X transfer: E_esc,fluid=E_esc,star+E_X and
L_out,fluid^infinity=L_esc,star^infinity+J_X^infinity. Later X release is counted once.
For a neutron occupied
state E_n=E_dep+E_esc,fluid, `mu_n-E_esc,fluid=(mu_n-E_n)+E_dep`.
**Fermi-hole heating is not an additional independent term.** Supporting MPR21 three-piece
collision decomposition reduces to the same residual (preflight section 5).

Exactly one representation: R-a {chemical term,total matter-energy transfer}; R-b uses
`sum_removed(mu_i-E_i)+sum_created(E_j-mu_j)+other independently deposited energy`
per event, with stoichiometric multiplicities, one event rate and one common rest-inclusive
energy zero; R-c is the complete residual. The neutron R-b specialization is hole plus
independently defined deposition. No hole or deposition addend
on top of an inclusive residual. n->chi gamma is illustrative only: absorbed photon gives
mu_n-E_chi; both products leave gives mu_n-E_n; partial deposition is accounted by actual
remaining fluid-exit energy. No rate or RMF kinematics is ratified.

For the declared cold spontaneous neutron-removal conditions, E_n<=mu_n, E_dep>=0 imply
0<=Q_dir and the rate-weighted hole lower bound. Actual mu supplies the general upper budget;
`P_dir^infinity<=|Eeq_dot|-L_esc,rest,star^infinity` additionally requires equilibrium or the
demonstrated neutron-sink sign delta mu_n^infinity<=0. A truly escaping massive product must
satisfy `e^Phi E_local>=m c^2`; the sharper pointwise cold bound is
`Q_dir<=mu_n-e^-Phi sum_truly_escaped m_j c^2` in common local energy units.
Here L_esc,rest,star^infinity=sum_a int e^Phi R_a(sum_truly_escaped m_j c^2)dV
uses rest energies in erg (or converted exactly once from MeV) and one lapse for event counts
and the minimum Killing energy at infinity; it is not a
two-lapse integral of local rest mass. Uniform nonrelativistic occupation gives 2E_F,kin/5; no universal positive
lower bound exists for an arbitrarily Fermi-surface-tuned process. Invisible is not generally
zero heat. The full conditions and relativistic average are preflight R9-R10.

## 6. Moving-equilibrium chemical response

The primary chemical variables and accepted contract evolution are

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

### Normative reconciliation of the forward BNV seam (R2 E-1)

ADR-0013 Q1's designation of G_y as physical authority for later non-fixed-baryon sources
(`docs/adr/ADR-0013-corrected-rotochemical-chemical-coefficients.md:212`) and ADR-0014
section 3.17's forward-looking BNV seam description
(`docs/adr/ADR-0014-secular-rotochemical-evolution-contract.md:472`–`:478`) are narrowed by
this explicitly owner-ratified ADR. The coefficient mathematics, Q/Z/W ownership,
Cowling contract inside its declared Phase-5 scope, governed baselines and all Phase-5 results
remain unchanged. G_y remains the governed unreduced authority used to build Z, the owner
of that Cowling response, and the source of diagnostic k. It is NOT the physical raw-S_y
map to eta_dot for a changing-B star.

For BNV, the statement that a source cannot be represented in the fixed-baryon two-channel
space applies to the RAW source S_y. The physically relevant moving-reference source is
Sigma_y=S_y-t Bdot: b^T Sigma_y=0 by identity, so Sigma_y=L sigma for every charge-consistent
source, with sigma=P(S_y-t Bdot). The physical BNV seam is therefore {Bdot, sigma}, with
qualified structural t and eta_dot=-Z(R+sigma)+2W Omega dotOmega
[+Zdot Z^-1eta when applicable], not the raw-G_y source map. This is an explicit accepted
forward-contract narrowing, not a change to accepted ADR-0013/0014 coefficient mathematics or
Phase-5 numbers. The owner explicitly acknowledged and approved this narrowing on 2026-09-18.

At spin OFF, t_i=B_i/B_B uses the existing Phase-5B structural derivative owner, with identical
neighboring-star/domain/surface policy, denominator qualification, errors and currency.
It is the accepted primary structural object for this contract; no production BNV object is
implemented.
Required identities are b^Tt=1 and t_p=t_e+t_mu. On the Structure-1 free-gas fixture,
t approximately (0.965770,0.0308522,0.00337774).
The propagated numerical budgets and chart extensions are preflight section 6.

Cowling k=G_yb/(b^TG_yb) approximately(0.992218,0.00748454,0.000297569) is diagnostic only;
ON THE STRUCTURE-1 FREE-GAS FIXTURE, t_e/k_e approximately 4.12 and t_mu/k_mu approximately 11.35. G_y builds governed Z and diagnostic k;
it MUST NOT map raw BNV S to physical eta_dot or reconstruct actual changing-B individual mu.
The rejected route equals -Z(S_l-k_l Bdot), not A3. Its old -3.442789339881122e-56 MeV/count
both-channel value is a **COWLING k-ROUTE NEGATIVE ORACLE**.

On that fixture, physical neutron-only S=(Bdot,0,0), Bdot<0 gives sigma=t_l|Bdot| and
eta_dot/|Bdot| approximately(-1.4303e-55,-3.6281e-55) MeV/count. Both eta<0 imply capture R<0,
eta_lR_l>=0. For a physical slide, the rejected route spuriously gives
eta_dot/abs(Bdot)=(+1.086e-55,+3.284e-55) MeV/count. Z+t suffice after moving-reference subtraction.
Long digits are arithmetic-reproducibility oracles on governed bytes, not physical precision;
propagated t budgets are approximately (1.443e-7,4.353e-9,2.179e-9) absolute, not certified
intervals. The fixture's d ln N_mu^eq/d ln B approximately 54 requires a particularly strict
muon population/response drift check in the frozen-depletion budget (preflight section 12).
True-response findings remain non-governed diagnostic evidence; no Phase-5 coefficient changes.

## 8. Chemical free-energy reservoir and individual potentials

At fixed current B, `E_chem=ell^TZell/2=eta^TZ^-1eta/2>=0`, a state reservoir, not heat.
The Phase-5 Lyapunov quantity is V=eta^TZ^-1eta=2E_chem.
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

Using MeV actual potentials and erg/s powers, the exactly-once boundary gives
P_dir^infinity=-C_(MeV->erg)g_actual^{infinity T}S-L_out,fluid^infinity.
The equilibrium decomposition satisfies

```text
P_dir^infinity(actual)=P_dir^infinity(eq)+C_(MeV->erg)eta^Tsigma,
P_dir^infinity(eq)=-C_(MeV->erg)mu_B^infinity Bdot-L_out,fluid^infinity
C_* T_infinity_dot=P_dir^infinity(actual)+L_H^infinity
                  -Lnu_eq^infinity-DeltaLnu^infinity-Lgamma^infinity-Lother^infinity,
L_H^infinity=C_(MeV->erg)eta^TR.
                                                               (A6)
```

For the neutron sink the plus sign follows directly from mu_n,actual^infinity=mu_B^infinity+eta^Tt_l.
At spin OFF/frozen coefficients, A5+A6 give
`Eeq_dot+C_(MeV->erg)Echem_dot+Udot=-L_out,fluid^infinity-Lnu_full^infinity-Lgamma^infinity-Lother^infinity`,
where E_eq and U are in erg and E_chem in MeV. Preflight R24 plus R5 with actual potentials also
close the ledger directly; that is the strongest proof independent of the quadratic split.
Add an explicit X reservoir to convert fluid-exit accounting to star-escape accounting.
There is no free reference-work addend, no dot E_chem heat, no generic PdV heat.
Changing coefficients/thermal structure or imposed spin requires its real state derivatives/
Omega dot J work; the spin-off identity cannot silently claim all-orders rotating ADM closure.

The static neutral stellar first law is independently derived by the TOV variation kernel in
preflight R23-R24 around a hydrostatic/diffusive-equilibrium beta-disequilibrated Tolman star;
beta equilibrium is needed only for g_eq=mu_B b. The cold limit gives d(M_eq c^2)/dB=mu_B^infinity in common energy units
(C_(MeV->erg) multiplies MeV mu_B for Mc^2 in erg), with surface-work qualification.
Future finite-p_cut dM/dB oracles must declare the
-4pi R^2 p_R dR/dB surface term convention.
Uniform-rotation extension is conditional on the stationary variational problem and proper
Killing conjugates. General quotable theorem authority remains source-limited, not a blocker
for the explicitly derived spin-off contract.

## 10. Beta-mediated thermal response

Retain governed `L_H=C_(MeV->erg) sum eta_lR_l`, `DeltaLnu=sum Ltilde_a(F_a-1)T^q`,
`DeltaP_beta=L_H-DeltaLnu` with the same coefficients and signed rates. This is the actual-
potential decomposition. Relative to P_dir(eq), spin-off/frozen correction is
`-C_(MeV->erg)Echem_dot-DeltaLnu`, nonpositive while storage fills, exactly -DeltaLnu in reached QSS.
Do not confuse these decompositions or assert storage always fills.

The following roots, minima, B1/B2 and ratio bounds belong to the governed NON-SUPERFLUID
R1995/FR2005 polynomial model; they are not universal under superfluidity.
For 0<|xi|<4.9097100289 (M) or 4.7870134733 (D), the incremental process power is negative.
Instantaneous cooling magnitudes are bounded by 0.467659 Lnu_eq,M and 0.527775 Lnu_eq,D;
DeltaP_beta/L_H lies in [-1/2,5/8] (pure D upper 1/2); at eta=0 use the limit. Preflight R21 proves these from the
polynomials. The approximately 0.724 keV/baryon source-driven cooling scale at 10^8K additionally
requires the capture/filling history assumptions stated there. It is not a bound on arbitrary
initial stored imbalance. Ordinary 10-30 MeV direct sources dominate that scale.

QSS R=-sigma is Z-independent and needs nonzero normalization in every driven channel.
Every QSS claim requires tau_relax(T,xi) << relevant evolution time IN THE APPLICABLE REGIME.
R0's 3.4e12/9e9/3e6 yr build estimates at abstract 1e-17/1e-14/1e-10 yr^-1 concern the
LARGE-|xi| MODIFIED-URCA QSS. At 10^8 K the corresponding asymptotic xi estimates are only
0.43/1.15/4.27, not deep large-|xi| values. R2 E-17 instead finds fixture linear relaxation
scales approximately 3.5e4 yr (electron) and 2.9e4 yr (muon), proportional to T^-6.
Hot linear QSS can be reached rapidly; cooling lengthens relaxation, permitting freeze-out
and later drive-dominated evolution. Large-|xi| QSS is a distinct possible later regime.
Use eta=eta0-Z int sigma dt only while reactions are negligible (early t<<tau_lin or after
freeze-out with a verified small reaction term). These are fixture diagnostics, not physical rates.

## 11. Generic versus process-dependent inputs

Class A: total baryon loss plus sequence; B: ordinary stoichiometry/source; C: energy partition;
D: product fate; E: hidden-sector interactions/accumulated state. Sequence and chemistry are
generic conditional on A+B. Absolute thermal predictions require C-E and any applicable external inflow F. A single rate API must
not encode a hidden energy efficiency. Event identity links particle and energy owners.
Class F: external inflow owns incoming energy, conserved charges and angular momentum for
capture-induced or externally driven BNV. Source callbacks may depend on t and state, including
T, eta and structure: Sigma_y(t,state)=S_y(t,state)-t_sequence Bdot(t,state), the moving-reference
source of section 7. Preserve the distinct raw S_y and reduced sigma semantic outputs.

## 12. Regime-I validity

Every terminal species must satisfy negligible production blocking/accumulation feedback,
stress/EOS and chemical-equilibrium influence, heat capacity, radiation/conduction/opacity,
and declared coefficient-drift budgets. Separately test X-mediated modification of ORDINARY
weak-reaction kinematics/normalization (effective-mass shifts, new spectator/catalyst channels,
Urca thresholds/normalizations). Such feedback can survive DeltaB/B->0 and is not captured by
a depletion-proportional coefficient-drift estimator. N_X/B<<1 alone is insufficient. Near zero source,
eta or direct power use absolute tolerances, not impossible relative inequalities.
Frozen scope additionally bounds variations of Z,W,Ltilde,C_*,t,surface gravity,I and support.
No universal tolerance is chosen. Explicit testable criteria are preflight section 12.

## 13. Regime-II boundary

Mandatory fate flags: PROMPT_ESCAPE (massive terminal products must have e^Phi E_local>=mc^2),
SM_THERMALIZATION, BOUND_INERT, BOUND_INTERACTING,
including mixed branches and terminal products. Failure of X thermal/mechanical/chemical
conditions triggers the corresponding Regime II. Pure ordinary coefficient drift exits the
frozen contract but can lead to evolving ordinary Regime I, not necessarily hidden-sector II.
R25 is a LINEARIZED TOY: rate cancellation needs fixed f' and Q as well as linear accumulation.
Goldman's actual degenerate mirror fluid is not rate-independent: f' depends on N_X through
E_F'(N_X). Coupling, accumulation time and state-dependent transport control thermal importance.
X-induced ordinary weak-rate changes also require the corresponding chemical/thermal Regime-II
feedback owner, even if the depleted fraction tends to zero.

Proposed progression: ordinary star -> MixedStar-lite test fluid with N_X/T_X ledger -> full
two-fluid state(B,B_X,J,S,S_X). Slaved T_X requires rapid relaxation and P_exchange=L_X,rad.
Goldman supports qualitative transfer through a mechanically dilute accumulated sector, with
quantitative source limitations. Cooling relative to Regime-I prediction differs from cooling
relative to passive control; constant positive deposition and a loss vanishing at zero T imply
a positive temperature floor solving P_dir=L_std(T_floor)+L_X(T_floor).
If L_X=lambda T^p is only the added hidden loss, (P_dir/lambda)^(1/p) is an upper bound,
not the exact floor. No Regime-II implementation or rate is authorized.

## 14. Matched control and sign classification

The comparator is the **PASSIVE SAME-INITIAL-B0 NO-BNV STAR**. Match initial star/T/eta, standard microphysics, photon/neutrino models, spin, solver and
tolerances; remove only BNV source/product sector. First toy holds corresponding frozen B0.
Report DeltaT, DeltaT_s, DeltaLgamma, integrated direct/incremental powers and DeltaU;
`d DeltaU/dt=DeltaP_total` includes changed equilibrium neutrinos and photon feedback.

Every NET HEATING/NET COOLING/NEAR ZERO/SIGN UNRESOLVED label MUST name its observable and
time/window: e.g. DeltaT_infinity(t), DeltaT_s(t), DeltaLgamma(t), DeltaU_th(t), or integrated
DeltaP on [t0,t1]. Use SIGN UNRESOLVED when uncertainty spans zero. R0's provisional 1e-4
relative temperature floor must be remeasured and combined with channel-weighting-specific
finite-T, partition, coefficient and t/Z errors. Linear-response lag is C_*T/(a_cool L).
The quoted fixture range 1.4e5-1e6 yr instead denotes energy-content/luminosity time C_*T/L;
it inherits free-gas C_* and placeholder controlled Ltilde, not physical-star predictions.
In future sliding-background comparisons, separate radius, g_s, envelope and C_* changes from
genuine thermal changes when interpreting DeltaT_s.
Generic Regime-I cooling cannot be claimed; a tuned small-direct-energy corner and applicable
history are required. No positive heating efficiency or arbitrary-history sign theorem.

## 15. First controlled toy contract requirement

Spin OFF, whole-star governed free gas, uniform abstract proper neutron sink; t from qualified
Phase-5B inputs, governed Z unchanged. P0: E_esc,fluid=actual current local mu_n, artificial cold zero-direct
entropy oracle, not physical n->chi gamma; P1: E_esc,fluid=E_n, uniform-sea hole diagnostic; P2:
E_esc,fluid=0 maximal retention. P1gamma kinematics may follow only after source clarification.

R0 two-tier design examples are mathematical fixtures only: P1/P2 around 1e-17 yr^-1, P0 around
1e-13-1e-12 yr^-1 for 1e7-1e8 yr or transient analytic tests. This task selects no physical
BNV rate, implements no rate/trajectory and chooses no microscopic constraints/couplings.
THE CURRENT P0 TIER IS PRIMARILY AN ETA-EVOLUTION / MOVING-REFERENCE / TRANSIENT ORACLE.
At 10^8 K and those fixture rates, reachable beta thermal power is about 171/17 times smaller
than the omitted +T s_,n R_a contribution and far below the sign-resolution floor (R2 E-18).
It is NOT a resolved thermal-sign experiment. A future purely mathematical thermal-sign stress
test must include the finite-T terms and remeasure resolution; no new physical rate is chosen.

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
normalization limitations; A24 AppendixE ambiguities; author-note effective-mass notation and photon-luminosity dimensional defects;
and corrections to unsupported R0 Goldman-volume claims. No source or note is edited/promoted.

Missing primary general first-law/transport/RMF sources are not claimed read. Static derived
closure stands; the general rotating theorem, realistic-EOS response, detailed escape/rates,
Regime-II transport, superfluidity and complete G_true numeric reproduction remain future
source-qualified work. A note's title/omission is not author intent or proof of false physics.

## 18. Consequences and implementation ownership

BNV-14 through BNV-25 in preflight section 15 are accepted contract requirements under this ADR;
they are not separately numbered entries in the scientific-invariant register and are not an
implementation claim. They cover the cold bracket, moving reference/t, diagnostic-only raw-G route,
individual-potential/state-energy identity, once-only event representation, product regimes,
beta bounds/QSS, actual-mu P0, Cowling scope and matched-control floors.

Separate future owners: state-dependent process/charge source, external inflow of energy/charges/
angular momentum, energy partition, fate/X state, global integration,
structural t, governed Z/sigma, chemical storage, standard beta, thermal ledger and matched
comparison. Reuse generic Phase-5D machinery; no reaction-specific generic-driver branches.
No production source, tests, coefficients, baselines, EOS or data change under this ADR.

## 19. Rejected alternatives and double-count refusals

Rejected primary physics: raw D_BNV S; physical null k Bdot; old equal-channel sink oracle;
fixed-reference alpha^2/(2a)/E_2; free F_ref/Delta(delta g^T F) as physical thermal flow;
separate hole over complete mu-E_esc,fluid; dot E_chem or dot M_eff as heat; generic PdV/gravity heat;
duplicate weak bulk dissipation/radiation; bound-product positive heat counted twice;
P0 based on equilibrium mu in a disequilibrated star; generic Regime-I cooling; QSS without
reachability. Historical occurrences are negative tests only. The exact H1-H10 refusals are
preflight section 18, with boundary/energy-zero qualifications.

## 20. Ratification record

The full R2 independent review found **0 BLOCKING / 1 MATERIAL / 26 NONBLOCKING / 21 NOTE**.
The sole material issue was the previously implicit forward-contract narrowing now stated in
section 7. Correction commit `58b375631d6948ed254809567bd29050f9735089` applied E-1 through
E-27 without changing load-bearing physics. The bounded R4 re-review found
**0 BLOCKING / 0 MATERIAL / 0 NONBLOCKING / 3 NOTE** and concluded that this ADR was ready for
explicit owner ratification. Fable was not needed.

On 2026-09-18 the human owner explicitly ratified this ADR for the declared controlled
abstract-source Regime-I BNV thermal contract and knowingly approved the narrowing in section 7.
The complete owner decision, review provenance, exclusions, and status synchronization are recorded
in [the Phase-6A-0 ratification record](../validation/PHASE6A0_BNV_THERMAL_FIRST_LAW_RATIFICATION.md).

Acceptance governs the contract only. No physical BNV rate is selected; no `n->chi-gamma`,
realistic `E_esc`, realistic A18, superfluid, Regime-II, MixedStar thermal-evolution, production
source, test, trajectory, or numerical baseline is implemented or authorized by this ratification.
