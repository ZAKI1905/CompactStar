# Phase 5A-5R — owner ratification of Track-R free-gas local thermodynamics

> **HUMAN OWNER RATIFICATION:**
> **PHASE 5A-5 TRACK-R FREE-GAS LOCAL THERMODYNAMIC COVERAGE IS ACCEPTED.**

## 1. Ratification record and authority chain

| Item | Ratified value |
|---|---|
| Date | 2026-09-05 |
| Implementation SHA | `933494d86daf2cf8965079ece49fabd66d9390e5` |
| Independent-review SHA | `097546f441ec4497cc426a5bb7051d53c2d59da7` |
| Starting canonical SHA | `41ab66bd6e6f351691173a1e2a033c646ffd3772` |
| Change class | documentation; human-owner scientific ratification and current-state synchronization |

The authority chain is `GOVERNANCE.md` -> accepted ADR-0010 -> the Phase-5A-5
implementation and validation record -> the independent review -> this human-owner
ratification. ADR-0010 governs the cold, charge-neutral local provider boundary and
explicitly excludes global paper `B/Z/W`, evolved chemistry, superfluidity, and BNV
(`docs/adr/ADR-0010-rotochemical-off-equilibrium-thermodynamic-contract.md:370-379`).
The implementation record establishes the final local branch and its evidence
(`docs/validation/PHASE5A5_TRACKR_PE_BRANCH.md:63-80`, `:231-280`); the independent
review finds that evidence acceptable with nonblocking follow-up and ready for human
ratification and canonical integration
(`docs/validation/PHASE5A5_TRACKR_PE_INDEPENDENT_REVIEW.md:3-18`, `:804-828`).

The project owner now ratifies the code at implementation SHA
`933494d86daf2cf8965079ece49fabd66d9390e5` as scientific authority for the Track-R
local free-gas source model within the governed source domain described below. The
independent review at `097546f441ec4497cc426a5bb7051d53c2d59da7` remains review
evidence; this record supplies the separate human authority required by
`GOVERNANCE.md:152-170`.

## 2. Exact claim ratified

The ratified claim is:

> **TRACK-R FREE-GAS LOCAL THERMODYNAMIC COVERAGE HUMAN-RATIFIED AND COMPLETE —
> WHOLE-STAR STRUCTURE REPRODUCTION READY.**

Within the FR2005 noninteracting whole-star free-gas source model, the equilibrium
local EOS/thermodynamic provider is implemented and validated over the complete
density range from the zero-pressure vacuum surface through, but not including, the
Sigma-minus appearance ceiling. The independent review found no missing density
interval and no blocking local scientific or numerical defect
(`docs/validation/PHASE5A5_TRACKR_PE_INDEPENDENT_REVIEW.md:771-802`).

The complete active-set ladder is:

```text
vacuum -> p-e -> n-p-e -> n-p-e-mu -> Sigma-minus ceiling
```

Its active-response dimensions are respectively value-only, 1, 2, and 3; exact
species thresholds are value-only boundaries with no fabricated inactive-species
derivative. Vacuum has `epsilon=0`, the finite one-sided conjugate limit
`h_pe -> m_p+m_e`, and no finite Hessian
(`docs/validation/PHASE5A5_TRACKR_PE_INDEPENDENT_REVIEW.md:178-220`, `:283-324`,
`:497-520`). The independently verified threshold values are:

| Boundary | Value |
|---|---:|
| neutron appearance | `7.3567289037328326656351713e-9 fm^-3` |
| pressure at neutron appearance | `1.8964875026317865961252521e-9 MeV fm^-3` |
| muon appearance | `0.45698480541241986996 fm^-3` |
| Sigma-minus ceiling | `0.6173552079665349801 fm^-3` |

The neutron threshold and pressure were independently re-derived at 90 digits
(`docs/validation/PHASE5A5_TRACKR_PE_INDEPENDENT_REVIEW.md:222-281`); the muon and
Sigma-minus values remained bit-identical to the prior verified branch
(`docs/validation/PHASE5A5_TRACKR_PE_INDEPENDENT_REVIEW.md:471-481`). The p-e/npe
tangent restriction and inverse-response collapse are correct, including the exact
normalization
(`docs/validation/PHASE5A5_TRACKR_PE_INDEPENDENT_REVIEW.md:326-401`).

## 3. Validation evidence accepted

The owner accepts the Phase-5A-5 evidence: PE-V1 through PE-V13, R1-V1 through
R1-V10, RFG1 through RFG11, and V1 through V10 pass; the implementation run recorded
26/26 self-contained and 49/49 authenticated full-suite passes, and all eight governed
baselines remained byte-identical
(`docs/validation/PHASE5A5_TRACKR_PE_BRANCH.md:231-295`). The independent review
re-derived the active-set ladder, p-e thermodynamics, vacuum and neutron boundaries,
response-limit identities, and N-1 policy, repeated the focused evidence, and found no
blocking defect (`docs/validation/PHASE5A5_TRACKR_PE_INDEPENDENT_REVIEW.md:3-13`,
`:621-680`).

N-1 is accepted as closed by the explicit `2^30`-ULP fail-closed rule. The review
measured the worst allowed `H00` relative error at the policy boundary as
approximately `5.03e-7`; classification remains physical, no density floor or nearby
branch substitution occurs, and unresolved responses fail explicitly
(`docs/validation/PHASE5A5_TRACKR_PE_INDEPENDENT_REVIEW.md:403-469`).

## 4. Claims explicitly not ratified

This decision does **not** claim that the FR2005 whole-star benchmark, any TOV star,
stellar susceptibility integration, particle-number or equilibrium-sequence response,
paper `B/Z/W`, rotochemical evolution, APR/BPAL, DS(CMF) off-equilibrium physics,
superfluidity, or BNV has been reproduced, implemented, or validated. Those exclusions
are the implementation's explicit scientific boundary
(`docs/validation/PHASE5A5_TRACKR_PE_BRANCH.md:338-349`) and were independently
confirmed (`docs/validation/PHASE5A5_TRACKR_PE_INDEPENDENT_REVIEW.md:728-769`).

This ratification does **not** resolve INV-09 or INV-11 and does not authorize paper
`B/Z/W` or evolution work. INV-09 still lacks the measure-complete global
particle-number and baryon-conserving sequence reduction
(`docs/SCIENTIFIC_INVARIANTS.md:712-733`); INV-11 still lacks an evolved chemical-state
ordering and redshift convention (`docs/SCIENTIFIC_INVARIANTS.md:838-853`).

## 5. Nonblocking findings carried forward

The open hardening and design findings remain open: R2, R3, R5, R9, N-2, N-3,
muon-side N-4, N-6, and P-1 through P-7 where applicable. Their independently
authenticated dispositions and nonblocking status are preserved without cleanup or
reinterpretation
(`docs/validation/PHASE5A5_TRACKR_PE_INDEPENDENT_REVIEW.md:682-726`).

**P-7 is a mandatory input to the next whole-star planning task.** Within the N-1
fail-closed near-neutron-onset window, the smooth thermodynamic evaluation—including
`epsilon`—may be unavailable even though reliable equilibrium composition remains
available. The next whole-star structure task must decide how to obtain the equilibrium
barotropic structure through that numerically narrow region without weakening the
fail-closed derivative policy. This ratification does not select or implement that
policy (`docs/validation/PHASE5A5_TRACKR_PE_INDEPENDENT_REVIEW.md:688-694`,
`:721-726`, `:796-802`).

## 6. Next authorized scientific stage

The next authorized scientific stage is a **separate FR2005 whole-star free-gas
structure reproduction task**. It must begin by source-authenticating the FR2005
Table-1 structure targets—`M_max` before Sigma-minus appearance, central density,
`R`, and `R_infinity`—and determining the precise TOV pressure/energy-density
interface required from the ratified local provider while carrying P-7 explicitly.
This ratification does not begin that work and does not authorize global rotochemical
coefficients.
