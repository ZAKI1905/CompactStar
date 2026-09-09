# Phase-5D-0RAT — secular rotochemical evolution contract ratification

> **PHASE-5D SECULAR ROTOCHEMICAL EVOLUTION CONTRACT —
> SCIENTIFIC PREFLIGHT COMPLETE / INDEPENDENTLY REVIEWED /
> HUMAN-RATIFIED FOR CONTROLLED NON-SUPERFLUID V1 SCOPE —
> PRODUCTION IMPLEMENTATION NOT YET BEGUN.**

This is a documentation/governance ratification record. It changes no production source, test,
baseline, data, EOS, literature, build file, or computed result. It does not begin realistic A18
closure or BNV.

## 1. Authority and identity

| Authority | SHA / identity |
|---|---|
| Canonical `master` | `49ab2b8c2881b6ef7b9309307d18cea51d557f72` |
| Phase-5C human-ratified coefficient tip (`PHASE5C2_RATIFICATION_SHA`) | `27727016856a6a25a46e447c70e380722ea8ddbf` |
| Original Phase-5D preflight (`PHASE5D0_SHA`) | `080c5bcbb7c10242b6146da3d9fbee961b3d82e6` |
| Material-review revision (`PHASE5D0_REVISION_SHA`) | `109ebfbbb9543c5f8984f85b097a137f6cce8754` |
| Final review-finding closure (`PHASE5D0_CLOSURE_SHA`) | `5f04b5ef7cefc7ceb0d73fb0b3927bfbb28508be` |
| Branch | `analysis/phase5d-secular-rotochemical-evolution-preflight` |
| Worktree | `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase5d-rotochemical-preflight` |

The authenticated ancestry is exactly

```text
27727016856a6a25a46e447c70e380722ea8ddbf
  -> 080c5bcbb7c10242b6146da3d9fbee961b3d82e6
  -> 109ebfbbb9543c5f8984f85b097a137f6cce8754
  -> 5f04b5ef7cefc7ceb0d73fb0b3927bfbb28508be.
```

At entry, the worktree was clean and local/upstream/live branch refs all equalled
`5f04b5ef7cefc7ceb0d73fb0b3927bfbb28508be`. Local/origin/live `master` all equalled
`49ab2b8c2881b6ef7b9309307d18cea51d557f72`. The commit carrying this record defines
`PHASE5D0_RATIFICATION_SHA`; it is reported after commit and push rather than self-referenced here.

Authority: `GOVERNANCE.md:13`; `AGENTS.md:11`; the reviewed contract at
`docs/adr/ADR-0014-secular-rotochemical-evolution-contract.md:1`; and the complete preflight at
`docs/validation/PHASE5D0_SECULAR_ROTOCHEMICAL_EVOLUTION_PREFLIGHT.md:1`.

## 2. Scientific-review history and owner disposition

The original Phase-5D preflight was committed at `PHASE5D0_SHA`. The first independent review
returned a revision-required disposition with M-1 through M-8; the bounded documentation revision
at `PHASE5D0_REVISION_SHA` closed those findings by text. Phase-5D-0RR then returned R1 through R3;
the closure commit `PHASE5D0_CLOSURE_SHA` corrected the benchmark process-configuration basis,
added the independent production-facing `Ltilde` oracle, and restored the qualified INV-11a status.
The durable history is preserved in ADR-0014 sections 11–12 and preflight sections 34–35.

The owner-supplied final RR2 disposition governing this ratification is:

> **B — PHASE-5D FINAL BOUNDED RE-REVIEW PASS WITH NONBLOCKING FINDINGS —
> R1/R2/R3 CLOSED — ADR-0014 READY FOR HUMAN-OWNER RATIFICATION WITH EXPLICIT CAVEATS.**

Final review totals are **0 BLOCKING, 0 MATERIAL, 12 NONBLOCKING, and 7 NOTES**. Fable was not
needed. No unavailable external review report is reconstructed, and no absent finding detail is
invented. The durable closure text and the owner-supplied final disposition are the complete
review authority used here.

The human-owner disposition is:

> **RATIFY ADR-0014 WITH THE RETAINED REVIEW CAVEATS.**

ADR-0014 is therefore **ACCEPTED / HUMAN-RATIFIED** for the controlled non-superfluid v1 scope.
Acceptance governs the contract; it does not validate or claim an implementation.

## 3. Ratified state, sign, and secular equations

The evolved chemical state is

```text
y_chem = (eta_npe^infinity, eta_npmu^infinity) [MeV]
ordering = (Npe, NpMu)
eta^infinity = e^nu eta_local
xi = eta^infinity / (k_B^(MeV) T_infinity).
```

`DeltaGamma > 0` means the net reaction acts to destroy positive `eta`; in the accepted forward
direction it denotes net neutron decay and lepton creation. The thermodynamic sign is
`eta_l DeltaGamma_l >= 0`. For frozen coefficients,

```text
dot eta = -Z R + 2 W Omega Omega_dot.
```

The global reaction-rate integral carries `e^(+nu)`. Neutrino energy at infinity carries
`e^(2nu)`. Chemical heating is

```text
P_H,chem^infinity = sum_l eta_l^infinity R_l [MeV/s]
L_H^infinity = C_(MeV->erg) P_H,chem^infinity [erg/s].
```

`Z`, `W`, `I_Omega`, and `Ltilde` are frozen for v1. If `Z` becomes time-dependent, the accepted
equation requires the additional `+ dot(Z) Z^-1 eta` term. No spin-down law is hard-wired; the
secular layer consumes a supplied spin history. Authority: ADR-0014 sections 3.1–3.5 and 3.11–3.12;
preflight sections 5–8 and 18–19.

## 4. Ratified unit contract and energy ledger

For every channel `a`,

```text
Ltilde_a [erg s^-1 K^-q_a]
R_a = (Ltilde_a / k_B^(erg)) T_infinity^(q_a-1) H_a(xi) [count/s]
Z [MeV/count] R [count/s] = MeV/s.
```

`k_B^(erg) [erg K^-1]` is derived from the single repository Boltzmann authority and the governed
energy conversion. A duplicate literal is forbidden. Canonical chemical power remains in MeV/s
and is converted exactly once to erg/s at the thermal-luminosity/RHS boundary. Authority:
ADR-0014 sections 3.2, 3.4, 3.5, and 3.7; preflight sections 8.2, 10.1, and 13.

The accepted neutrino/heating ledger is

```text
L_nu,eq = Ltilde_a T^q
L_nu    = Ltilde_a F_a(xi) T^q
Delta L_nu = Ltilde_a [F_a(xi)-1] T^q
R_a     = (Ltilde_a/k_B^(erg)) H_a(xi) T^(q-1).
```

The same declared `Ltilde_a` feeds equilibrium cooling, the `F` correction, the `H` reaction rate,
and chemical-heating bookkeeping. At `eta=0`, `H=0`, `F=1`, `R=0`, `Delta L_nu=0`, and `L_H=0`.
This is the ratified RE9 **SAME-COEFFICIENT** identity. The historical placeholder `K_DU`/`K_MU`
cooling coefficients are not scientific normalization authority. Authority: ADR-0014 section 3.7;
preflight sections 10, 12–13, and RE9 in section 31.

## 5. Imbalance function, benchmark, and domain contract

The corrected modified-Urca `H_M` final term uses `pi^8`. Its exact classification is
**CONFIRMED PRINTED TYPO / INTERNAL SOURCE INCONSISTENCY; NO PUBLISHED ERRATUM LOCATED**. No
published erratum is claimed. Authority: ADR-0014 section 3.6; preflight section 9.3.

The controlled free-gas architecture benchmark is intentionally **MODIFIED-URCA-ONLY**:

```text
enabled = {Me, Mmu}
disabled = {De, Dmu}
D_De = empty
D_Dmu = empty.
```

Those empty domains follow from the benchmark process-configuration contract, not from a universal
claim about DU kinematics. The documented electron-DU triangle sliver at approximately
`7.36e-9` to `6.67e-8 fm^-3` is retained. Kinematic allowance, microphysical applicability, and
declared benchmark process configuration remain distinct. `nB_min` is not a DU threshold.
The benchmark coefficients are **MATHEMATICAL / ARCHITECTURE BENCHMARK NORMALIZATIONS**, not
realistic microphysical predictions. Authority: ADR-0014 sections 3.10 and 3.15; preflight
sections 15 and 25.

`D` is the declared chemical/reaction domain and `D_a subseteq D` is the process/channel support
domain. No implicit core boundary, saturation-density boundary, crust cutoff, or figure-derived
pressure boundary is accepted. Global coefficients integrate only over the declared `D_a`.
Authority: ADR-0014 sections 3.4, 3.10, and 3.14; preflight sections 6.1 and 11.

## 6. Ratified ownership architecture and future oracle

Five ownership layers remain separate:

1. The channel microphysics/normalization provider owns process, lepton, branch, effective masses,
   matrix-element and alpha/beta factors, local `S_a`, applicability/support, and provenance.
2. `GlobalUrcaChannelCoefficient` owns GR stellar integration over `D_a` producing `Ltilde_a`.
3. The imbalance-functions layer owns only `F_a(xi)` and `H_a(xi)`.
4. The reaction response combines `Ltilde`, `F/H`, `T_infinity`, and `eta_infinity` to produce
   `R`, `Delta L_nu`, and chemical heating.
5. The secular RHS owns state/source assembly only.

No monolithic owner is accepted. Authority: ADR-0014 section 3.13; preflight section 33.1.

Future implementation must provide RE10b, an oracle independent of the production integration
kernel, for

```text
Ltilde_a = integral_{D_a} 4 pi r^2 e^lambda S_a(r) e^[(2-q_a)nu] dr.
```

It must discriminate omitted/extra/wrong-sign lapse, omitted/inverted proper-volume factor, wrong
`q`, and wrong domain. Its expected value cannot be an injected/precomputed `Ltilde`. Authority:
ADR-0014 section 5; preflight RE10b in section 31 and mutations M10/M11, M12, M30, M33–M37 in
section 32.

## 7. Lyapunov qualification

For frozen symmetric SPD `Z` and zero spin drive,

```text
V = eta^T Z^-1 eta,
dot V = -2 sum_l eta_l R_l <= 0.
```

The correct general claim is **NONINCREASING**. Strict decay requires every relevant nonzero
imbalance direction to be dissipatively coupled to an active positive-normalization reaction
channel. An uncoupled dead imbalance may freeze. With `Z` cross-coupling, a dead individual
reaction channel does not imply that its `eta` component remains fixed. Authority: ADR-0014
section 3.3; preflight sections 7.3 and 21 limit B.

## 8. Source status and realistic blockers

- R1995 equation 34 is historical/supporting, npe-lumped, described by its source as somewhat
  uncertain, and is not definitive realistic FR2005 normalization.
- Yakovlev et al. 2001 is publicly available and sufficient **in form** for the planned
  normalization architecture. The exact `alpha_n` prescription remains unresolved.
- Bounded primary-source inspection establishes that APR/A18+delta-v+UIX* contains enough
  information **in form** to reconstruct the arbitrary-composition core functional. It is not
  installed as governed project authority.
- Realistic FR2005 reproduction remains source-limited and blocked on the durable ledger:
  authenticated APR/A18 authority; exact phase/Maxwell construction; crust treatment; rate
  normalization; effective masses; `alpha_n`; support authority; and benchmark arrays or governed
  digitization where required.

Authority: ADR-0014 section 6; preflight sections 10.4 and 26–27. This ratifies source status only,
not realistic physics implementation.

## 9. Retained final-review caveats and implementation requirements

The 12 nonblocking findings are retained through the following explicit requirements. These are
future implementation obligations, not work performed by this ratification:

A. Assign an explicit owner for the static declared enabled-process-set configuration.

B. Make the low-density DU-sliver negative-control fixture constructible even though the physical
sliver lies below the current `nB_min` guard.

C. Predeclare an actual numerical acceptance tolerance for RE10b.

D. Give the RE10b synthetic fixture sufficiently nontrivial `nu`/`lambda` amplitude to discriminate
lapse and proper-volume mutants.

E. Give the wrong-domain RE10b fixture nonzero `S_a` outside `D_a`, preventing a vacuous pass.

F. Preserve the domain-mismatch mutation's relation to `G_y`/domain consistency where applicable.

G. Cite the accepted upstream INV-11 authority at governance/status sites where appropriate.

H. Treat self-authored closure/status banners as evidence labels, not independent validation.

I. Preserve the pre-existing M10/M11 cross-reference/mislabel as non-new scientific debt and not a
ratification blocker; in the shared integrated-`Ltilde` representation they are one algebraic
mutation and receive no duplicate coverage credit.

J. Make innermost-first profile ordering an explicit future assertion wherever integration logic
relies on it.

K. Add a production negative control for disconnected or outer-region DU-support hazards.

L. Do not use the controlled free-gas benchmark as validation of realistic `Ltilde` construction or
normalization.

The closure text additionally fixes the following implementation details: RE10b declares every
ordinary unit prefactor and uses independent closed-form/high-precision quadrature; the rate-power
check uses two distinct temperatures; `D`/`D_a` belong in coefficient provenance; muon support is
represented separately; and frozen v1 introduces no temperature-dependent support predicate.
These details do not increase the final finding totals. Authority: ADR-0014 sections 3.10, 3.14,
and 5; preflight sections 15, 19, 31, 32.2, and 35.

The seven review notes are recorded only by total because their details are absent from the durable
evidence supplied for this task. No details are inferred.

## 10. INV-11 disposition after ratification

| Subpart | Ratified status |
|---|---|
| **INV-11a** — coefficient-object redshift semantics | **PARTIALLY RESOLVED upstream; accepted ADR-0014 extends/clarifies the secular contract** |
| **INV-11b** — evolved `eta` state ownership | **CONTRACT RESOLVED / IMPLEMENTATION PENDING** |
| **INV-11c** — reaction sign/index convention | **CONTRACT RESOLVED / IMPLEMENTATION PENDING** |
| **INV-11d** — thermal energy ledger and no double counting | **CONTRACT RESOLVED / IMPLEMENTATION PENDING** |
| **INV-11e** — frozen coefficient lifetime/update policy | **CONTRACT RESOLVED / IMPLEMENTATION PENDING** |
| **INV-11f** — ODE/source coupling | **UNRESOLVED / IMPLEMENTATION + VALIDATION PENDING** |

**Global INV-11 remains UNRESOLVED.** Contract resolution of subparts is not implementation or
validation. Authority: accepted ADR-0013 section 10; ADR-0014 section 7; preflight section 30.

## 11. Dependency and implementation boundary

Phase-5C corrected coefficients are **IMPLEMENTED / CANDIDATE-VALIDATED / INDEPENDENTLY REVIEWED /
HUMAN-RATIFIED** at `PHASE5C2_RATIFICATION_SHA`, but are **NOT CANONICALLY INTEGRATED**. The exact
status is recorded in
`docs/validation/PHASE5C2_CORRECTED_CHEMICAL_COEFFICIENT_RATIFICATION.md:298`.

Phase-5D production work cannot begin on canonical `master` until the complete accepted Phase-5C
coefficient history is canonically integrated. This ratification performs no merge and leaves
canonical `master` unchanged.

Production Phase-5D implementation has **NOT BEGUN**. Controlled free-gas implementation is now
authorized under accepted ADR-0014 only after its dependency is satisfied and under a separately
governed implementation task. Realistic A18 closure remains **SOURCE-LIMITED / BLOCKED**; no A18
implementation has begun. BNV has **NOT BEGUN** and remains sequenced after validated standard
rotochemical evolution (`docs/MODERNIZATION_ROADMAP.md:797`).

## 12. Ratification result

ADR-0014 is **ACCEPTED / HUMAN-RATIFIED**. The Phase-5D scientific contract is **PREFLIGHT
COMPLETE / INDEPENDENTLY REVIEWED / HUMAN-RATIFIED / NOT IMPLEMENTED**. The controlled free-gas
rotochemical evolution is **AUTHORIZED FOR IMPLEMENTATION / NOT YET IMPLEMENTED**. INV-09 remains
**VERIFIED / RESOLVED**. Global INV-11 remains **UNRESOLVED**. Realistic FR2005/A18 remains
**SOURCE-LIMITED / BLOCKED**. BNV is **NOT BEGUN**.
