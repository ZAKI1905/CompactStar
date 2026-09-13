# ADR-0015: BNV open-system thermal ledger and unreduced chemical response

**Status: PROPOSED — NOT ACCEPTED.**
**Date:** 2026-09-13.
**Owner ratification:** NOT REQUESTED OR CLAIMED by this document.
**Canonical entry:** `0a7418aecb7314cfa472a78f1faf477be8456a94`.
**Scope:** proposed physics/semantic architecture only; no production BNV implementation,
rate, test, baseline, EOS, A18 reconstruction or change to governed Phase-5 machinery.
**Scientific disposition:** D, source-limited for the complete requested authority audit.

## Context and authority

The governed standard problem uses a charge-neutral unreduced susceptibility G_y, beta
reduction Z, spin forcing W, signed beta reactions, chemical dissipation and enhanced neutrino
loss. Baryon-changing sources generally leave the fixed-baryon beta subspace; using only Z
would erase necessary information. A thermal source cannot be inferred from the particle
loss rate alone. These premises are authenticated and derived in the companion
[Phase-6A-0 preflight](../validation/PHASE6A0_BNV_THERMAL_FIRST_LAW_PREFLIGHT.md), sections1-12.
The source hierarchy is `GOVERNANCE.md:19`; governed definitions remain
`docs/adr/ADR-0013-corrected-rotochemical-chemical-coefficients.md:136` and
`docs/adr/ADR-0014-secular-rotochemical-evolution-contract.md:188`.

Phase-5D controlled frozen-v1 is canonically integrated and closed for its declared scope at
entry. Its baseline hash is
`2606916915b2da5c051b1a06637c0a371a74751c6e77b2775bc629a63bc9f6dd`, with eleven governed
baselines total. Global INV-11 remains unresolved and realistic A18/FR2005 remains
source-limited. This proposal does not reopen, supersede or modify those contracts.

**Evidence classes are separate:** published source claims; author working notes/hypotheses
(non-authoritative, not published); independent preflight derivation. The author reports
additional unfinished unpublished BNV thermal work. The published corpus is not the complete
intended theory, and omission is never evidence of intent to exclude a term. The preflight's
section17 audits available local planning fragments separately and records the unavailable
working manuscript. If notes become available, audit them against the completed proposed
ledger before independent Opus review / owner ratification. No note can override the first law.

## Proposed decision

### 1. Fix the matter boundary, signs and energy zero

Use metric (-,+,+,+), including rest energy in epsilon and mu. Define ordinary particle
creation positive, `nabla_mu(n_i u^mu)=Gamma_i`. For
`nabla_mu T_matter^{mu nu}=Q^nu`, define **total** local matter-energy transfer
`q_E=-u_nu Q^nu`, positive into matter. Do not call q_E a nonchemical thermal source.
Contracting the perfect-fluid equation and combining Gibbs/Euler gives

```text
D epsilon+(epsilon+p)theta=q_E
D n_i+n_i theta=Gamma_i
T nabla_mu(s u^mu)=q_E-sum_i mu_i Gamma_i.             (A1)
```

A1 is an open-matter entropy balance, not a demand that the matter entropy source be
nonnegative. Escape can remove entropy. At fixed volume the exact finite-T temperature law
instead has `c_V D T=q_E-sum_i[mu_i+T s_,ni]Gamma_i`. Equating the A1 entropy-energy residual
to the controlled frozen thermal RHS is an explicitly declared cold/degenerate approximation;
its neglected composition-entropy correction must be bounded near a cancelling residual.
Independent derivation: preflight P1-P5, including assumptions and units.

### 2. Keep event rate and energy partition separate

For directed channel a, `Gamma_i^a=nu_ia R_a`, R_a nonnegative per local proper volume/time.
With external input E_in, actual escaping E_esc and a separately retained nonthermal transfer
J_X, the prompt controlled residual is

```text
q_direct^a=[E_in,a-E_esc,a-J_X,a-sum_i mu_i nu_ia] R_a. (A2)
```

Delayed reservoir release is counted once at release. For one-particle disappearance,
`q_direct=(mu-E_occ+E_dep)R=(mu+E_in-E_esc-J_X)R` only if the event partition satisfies
`E_occ+E_in=E_dep+E_esc+J_X`. E_dep is actual returned product energy in the first expression.
If it instead names the entire residual, adding it to the second expression repeats energy.

A microscopic `(mu-E_occ)R` rearrangement term partially represents A2; its sum with true
product deposition is the same complete residual. Never add another hole term on top of A2.
The chemical term `-mu Gamma` alone is not universally the microscopic hole energy: the
occupied-state removal energy and matter boundary must also be translated. No standalone
named hole formula was found in the four authenticated BNV papers, and no published duplicate
thermal addend was established. Exact additional-primary attribution remains SOURCE-LIMITED.
See preflight sections5-6 and17, not an inference about unpublished author work.

### 3. Preserve local-to-infinity measures

For a static normalized lapse, `d tau=e^nu dt`, `E_infinity=e^nu E_local`:

```text
dot N_y^BNV=int e^nu Gamma_y dV_proper
P_direct^infinity=int e^(2nu) q_direct dV_proper
L_esc^infinity=int e^(2nu) sum_a E_esc,a R_a dV_proper. (A3)
```

The escape expression assumes truly escaping energy emitted locally, after retention and
transport have been resolved. Energy deposited at another radius uses that deposition site's
clock/frame. Rotating emission requires the appropriate Killing-energy/angular-momentum
ledger, not an unqualified static substitution. A global number vector alone cannot determine
radially varying energy loss. Preflight P9-P11.

### 4. Require the exact governed beta limit

For positive neutron decay, `Gamma_n=-DeltaGamma_l`, `Gamma_p=Gamma_l=+DeltaGamma_l`.
Then `-sum mu Gamma=eta_l DeltaGamma_l`. Beta neutrino escape has `q_E,beta=-Q_nu,beta`;
“no external source” does not mean that escaping neutrino energy vanishes. Consequently

```text
L_H^infinity=C_M sum_l eta_l^infinity R_l
DeltaP_beta=L_H-DeltaLnu
DeltaLnu=sum_a Ltilde_a [F_a(xi)-1] T^q.              (A4)
```

Use the existing unit conversion C_M and the same governed Ltilde for rates/emissivity.
This recovers the Phase-5D sign and normalization exactly within its controlled assumptions.
The incremental beta term can cool, including the established small-imbalance modified-Urca
limit. A proposed BNV path that fails this reduction must stop. Preflight section8; governed
`CompactStar/Physics/Rotochemical/RotochemicalReactionResponse.hpp:1` and ADR-0014.

### 5. Use the full neutral susceptibility before reduction

Authenticate and retain

```text
y=(N_n,N_e,N_mu), N_p=N_e+N_mu
g_y=(mu_n,mu_p+mu_e,mu_p+mu_mu)
b=(1,1,1), L=[[-1,-1],[1,0],[0,1]]
z=N_y-N_y,ref=G_y delta g_y^infinity
eta^infinity=-L^T delta g_y^infinity
D_BNV=-L^T G_y^-1,  eta_dot|BNV=D_BNV S_BNV
Z=L^T G_y^-1 L,  D_BNV L=-Z.                         (A5)
```

D has shape 2 by 3, units MeV/count, ordered eta_e/eta_mu by n/e/mu; it has no symmetry
property. Store G and provenance; D is a derived solve/view, never an independently fitted
coefficient. Actual number balance is `Ndot=L R+S_BNV+boundary_sources`. Reference/structural
forcing enters zdot as `-Ndot_ref`, not as fictitious particle production.

Since `range L=ker b^T`, S=(-R,0,0) has nonzero baryon loss and cannot be Lr. Z alone cannot
recover generic D. Free-gas neutron removal lowers mu_n and immediately drives both eta
components negative; the authenticated numeric oracle gives each derivative/R equal to
`-3.442789339881122e-56 MeV/count`. A neutral paired proton/electron removal gives both
positive in that oracle. Neither is a realistic rate model. Charge-changing processes must
close the ordinary/product/field currents before entering this neutral chart.
A mathematical null source is `S=G b dot N_B/(b^T G b)`; it changes baryon number with D S=0,
not necessarily along a physical hydrostatic sequence. Preflight sections9-10 and16.

### 6. Keep stored energy out of the thermal source

For frozen symmetric positive G, with a fixed cold tangent/reference,

```text
E_2=1/2 z^T G^-1 z
alpha=b^T z, a=b^T G b
E_2=alpha^2/(2a)+E_beta
E_beta=1/2 eta^T Z^-1 eta
E_2dot|beta=-eta^T R
E_2dot|BNV=delta g^T S_BNV.                           (A6)
```

E_2 is tangent-subtracted quadratic storage; only E_beta is beta-relaxable at fixed baryon
number. Absolute cold energy also has a linear reference term. Storage is not heat. For
changing G/reference, retain all derivative and parameter-work terms, notably
`eta_dot_extra=L^T G^-1 Gdot G^-1 z`, and in a fixed-baryon Z representation
`eta_dot_extra=Zdot Z^-1 eta`. The storage derivative includes the corresponding negative
half quadratic metric derivative. Frozen scope requires measured depletion/coefficient and
absolute residual error bounds; no numerical tolerance is set here. Preflight sections10-11.

### 7. Compare to a matched no-BNV trajectory

Use the same initial star, T, eta, spin history, standard microphysics, photon/neutrino model,
solver and tolerances, disabling only BNV sources in the control. Define Delta as BNV minus
control at equal coordinate time. Direct BNV source is A2 integrated by A3; beta-mediated
response is `Delta[L_H-DeltaLnu]`. Total thermal response also contains changed equilibrium
neutrino cooling, other cooling, photon feedback and separately owned reservoir/work terms.
Define Delta T, Delta Lgamma and Delta U with the same frozen U(T)=int C(T)dT. Neither
instantaneous thermal power nor integrated energy difference is a positive efficiency.

The controlled conservation target is preflight P27/P28. For constant g_0,G and
`zdot=L R+S+F`, it retains the cold tangent reservoir and reference-forcing work:

```text
Delta U+Delta E_2
 =int [P_in-L_esc-g_0^T S+Delta(delta g^T F)
       -Delta Lnu_full-Delta Lgamma-Delta Lother]dt.   (A7)
```

F acts on the deviation; it is not a particle source. Common F does not cancel its work.
Full neutrino difference includes both equilibrium-temperature feedback and disequilibrium
enhancement. This is a frozen-fixture identity, not a full evolving-star ADM closure.

### 8. Separate structural and rotational work from heat

No source-backed irreversible structural dissipation mechanism was found for the audited
quasistatic BNV scope. Reversible PdV/gravitational readjustment changes background/equilibrium
coefficients; it is not an installed heating term. H67/H70 support restricted first-integral
and fixed-baryon/entropy variational results. The independently derived cold static zero-
surface-pressure limit is `delta(Mc^2)=mu_B^infinity delta N_B`; the general finite-entropy,
multispecies rotating stellar first-law authority remains source-limited. Surface-cutoff work,
rotation and changing structure cannot be hidden in a heat residual. Preflight section14.

## Proposed architecture, invariants and verification

Adopt the ten proposed owners in preflight section18: process identity, local particle source,
energy partition, global integration, full-G response, storage diagnostics, standard beta
coupling, thermal ledger, matched comparison, and sign diagnostics. Candidate BNV-1 through
BNV-13 there form part of this proposal, not ratified additions to SCIENTIFIC_INVARIANTS.
Generic execution must have no reaction-specific hardcoded exceptions.

Future analytic oracles are predeclared in preflight section16: zero source exact standard
trajectory; nonzero source/zero direct residual; escape equals chemical release; complete
retention; no weak reactions; restored-beta Lyapunov decrease; null-eta baryon direction;
charge closure; matrix/storage/redshift identities; energy conservation; variable coefficients;
finite-T near-cancellation uncertainty. They are not implemented now.

The future sign map separates direct, beta-mediated and total effects over source direction,
energy partition, eta_e/kT, eta_mu/kT, T, depth and weak regime. Report NET COOLING, NET HEATING,
NEAR ZERO, with an uncertainty overlay when the sign cannot be resolved. No scan was run.

## Alternatives rejected by derivation

- Adding an unconditional positive Fermi-hole efficiency to an inclusive first-law residual.
- Passing arbitrary BNV sources through Z without their unreduced baryon component.
- Counting E_2dot as a thermal source in addition to the first-law terms.
- Treating emitted energy, equilibrium-sequence mass change or rotational work as heat.
- Calling q_E-mu Gamma an exact finite-T temperature RHS without entropy-composition terms.
- Inferring the author's unpublished intended theory from the published corpus's omissions.

## Source limits, consequences and acceptance gate

The exact authenticated manifest, four-paper/145-page audit and equation-level evidence are
in preflight sections2 and17. The prospective double-count rule is resolved. Additional named-
hole primary attribution and general stellar-first-law authority remain source-limited.
A24 Appendix E prints escape-frame and rate-normalization ambiguities; resolve these before
using its energy/escape prescription, without claiming its numerical results are wrong.
The additional unpublished thermal notes were not available/audited; supplied notes require
a separate non-authoritative hypothesis audit before review/ratification when available.

This proposal is therefore **NOT ACCEPTED** and does not declare source closure or owner
ratification. Next obtain/authenticate the missing authority and complete the specific
crosscheck, then independent scientific review and explicit owner ratification before any
implementation. Do not begin that next task automatically. No BNV rate, realistic A18, new EOS,
production trajectory or canonical merge is authorized by this ADR.
