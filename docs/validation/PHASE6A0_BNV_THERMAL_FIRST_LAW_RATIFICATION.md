# Phase-6A-0RAT — BNV thermal first-law contract ratification

> **ADR-0015 ACCEPTED / HUMAN-RATIFIED — PHASE-6A-0 PREFLIGHT COMPLETE /
> INDEPENDENTLY REVIEWED / HUMAN-RATIFIED — CONTRACT ACCEPTED;
> PRODUCTION IMPLEMENTATION NOT BEGUN.**

**Date:** 2026-09-18
**Change class:** scientific-semantic and structural/architecture contract ratification;
documentation/governance only.
**Accepted decision:** `docs/adr/ADR-0015-bnv-open-system-thermal-ledger.md`
**Preflight:** `docs/validation/PHASE6A0_BNV_THERMAL_FIRST_LAW_PREFLIGHT.md`
**Diagnostic:** `docs/validation/PHASE6A0_COWLING_BARYON_DIRECTION_DIAGNOSTIC.md`

This record accepts the controlled abstract-source Regime-I contract. It changes no scientific
equation, production source, test, numerical method, baseline, EOS/data, literature, or computed
result. It creates no Phase-6 numerical artifact and begins no BNV implementation.

## 1. Authenticated identity and ancestry

| Identity | SHA / value |
|---|---|
| Canonical `master` at entry | `0a7418aecb7314cfa472a78f1faf477be8456a94` |
| Phase-6 reviewed branch SHA at entry | `58b375631d6948ed254809567bd29050f9735089` |
| Final correction SHA | `58b375631d6948ed254809567bd29050f9735089` |
| R2-reviewed parent | `7a31862f1e8a5cf046316882e05d76e4924e27d9` |
| Original Phase-6A-0 draft | `5a6bf7cb9455d684ddb6fccb22ad2b9fec940b3a` |
| Branch | `analysis/phase6a0-bnv-thermal-first-law-preflight` |
| Worktree | `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase6a0-bnv-thermal-preflight` |

At entry the Phase-6 worktree was clean and local/upstream/live branch refs all equalled the
reviewed SHA. The canonical checkout was clean and local `master`, cached `origin/master`, and live
remote `master` all equalled the canonical entry SHA. Canonical `master` was authenticated as an
ancestor of the reviewed Phase-6 branch. The commit carrying this record defines
`PHASE6A0_RATIFICATION_SHA`; Git reports it after commit and push rather than self-referencing it
inside the commit.

## 2. Independent-review evidence

The full Phase-6A-0R2 independent scientific review reported:

- **0 BLOCKING**;
- **1 MATERIAL** before correction;
- **26 NONBLOCKING**;
- **21 NOTE**.

The single material issue was contract meaning: the earlier proposal narrowed ADR-0013 Q1 and
ADR-0014 section 3.17 for the future BNV seam without stating that narrowing explicitly. Correction
commit `58b375631d6948ed254809567bd29050f9735089` applied E-1 through E-27 and made the narrowing
explicit while preserving all load-bearing physics. The durable R2 provenance and correction
checklist are in `docs/validation/PHASE6A0_BNV_THERMAL_FIRST_LAW_PREFLIGHT.md:1225` and
`:1248`.

The bounded Phase-6A-0R4 fresh-context re-review of
`7a31862f1e8a5cf046316882e05d76e4924e27d9..58b375631d6948ed254809567bd29050f9735089`
reported:

- **0 BLOCKING**;
- **0 MATERIAL**;
- **0 NONBLOCKING**;
- **3 NOTE**.

R4 confirmed that all 27 mandatory corrections were applied, load-bearing equations were
unchanged, and Phase-5 source, test, coefficient, baseline, and data bytes remained unchanged. Its
final disposition was **PROPOSED ADR-0015 READY FOR EXPLICIT OWNER RATIFICATION**. The three notes
were optional symbol/citation/label observations and were not prerequisites to ratification. The R4
report was read completely (81 lines), SHA-256
`ed8f7c91d33a05b98d12240cd22cf4afb7b7526f7973f54745d546ddabc8554b`.

**Fable needed:** NO.

## 3. Human-owner decision — recorded faithfully

The owner states:

> I ratify PROPOSED ADR-0015 for the declared controlled abstract-source
> Regime-I BNV thermal contract.
>
> I explicitly acknowledge and approve the forward-contract narrowing of
> ADR-0013 Q1 and ADR-0014 §3.17:
>
> ADR-0013 Q1's designation of `G_y` as the unreduced physical authority for later
> non-fixed-baryon sources is narrowed for the BNV response seam.
>
> ADR-0014 §3.17's statement that a BNV source cannot be represented in the
> fixed-baryon two-channel space applies to the raw source `S_y`.
>
> After subtraction of the physical equilibrium-sequence motion,
>
> `S_y - t Bdot`,
>
> the source is baryon-neutral and is exactly representable as
>
> `S_y - t Bdot = L sigma`,
>
> `sigma = P(S_y - t Bdot)`.
>
> The BNV response seam is therefore `{Bdot, sigma}`, using the governed `Z`
> response.
>
> This does not supersede ADR-0013/0014 coefficient mathematics, Q/Z/W
> ownership, governed baselines, the Cowling fixed-baryon contract, or standard
> rotochemical results.
>
> This ratification is limited to the declared controlled abstract-source
> Regime-I framework.
>
> It does not ratify any physical BNV rate, `n->chi-gamma` implementation or rate
> model, realistic `E_esc`, realistic A18 physics, superfluid extensions,
> Regime-II transport, MixedStar thermal evolution, or production
> implementation.

ADR-0015 is therefore **ACCEPTED / HUMAN-RATIFIED** for precisely that scope. Acceptance governs
the contract; it does not claim that a production implementation exists or that any Phase-6 number
has been computed.

## 4. Exact forward-contract narrowing

The owner knowingly accepts all of the following together:

1. **ADR-0013 Q1.** `G_y` remains the governed unreduced coefficient authority inside its declared
   Phase-5/R2006/Cowling role. It builds the governed `Z` response and supplies diagnostic `k`.
   It is not the physical raw-BNV source-response map for changing total baryon number.
2. **ADR-0014 section 3.17.** “Cannot be represented in the fixed-baryon two-channel space”
   applies to the raw source `S_y`. ADR-0014's raw generic source symbol `Sigma_y(t,state)` is the
   object written `S_y` in ADR-0015; ADR-0015 reserves `Sigma_y` for the moving-reference
   baryon-neutral remainder.
3. **Moving reference.** With

   ```text
   Bdot = b^T S_y
   t = (partial N_y^eq / partial B)_Omega
   b^T t = 1,
   ```

   the remainder satisfies

   ```text
   Sigma_y = S_y - t Bdot
   b^T Sigma_y = 0
   Sigma_y = L sigma
   sigma = P(S_y - t Bdot).
   ```

4. **Physical BNV seam.** The accepted seam is `{Bdot, sigma}` with qualified structural `t` and
   governed `Z`; for the declared frozen contract,

   ```text
   eta_dot = -Z(R + sigma) + 2 W Omega Omega_dot
   ```

   with `+Zdot Z^-1 eta` when a future authorized changing-`Z` contract applies.
5. **No Phase-5 change.** No Phase-5 coefficient mathematics is invalidated. No Phase-5 baseline
   is superseded. No `Q/Z/W` ownership changes. The standard fixed-baryon rotochemical contract and
   all governed Phase-5 results remain intact.

Authority: `docs/adr/ADR-0015-bnv-open-system-thermal-ledger.md:158`; the predecessor clauses are
`docs/adr/ADR-0013-corrected-rotochemical-chemical-coefficients.md:203` and `:212`, and
`docs/adr/ADR-0014-secular-rotochemical-evolution-contract.md:470`.

## 5. Ratified local open-system first law and direct-energy ledger

The accepted local sign convention is creation-positive `Gamma_i` and energy-into-the-ordinary-
fluid positive `q_E=-u_nu Q^nu`. The local open-system identities are:

```text
nabla_mu(n_i u^mu) = Gamma_i
D rho + (rho + p) theta = q_E
T nabla_mu(s u^mu) = q_E - sum_i mu_i Gamma_i
T n_B D s_b = q_E - sum_i mu_i Gamma_i - T s_b Gamma_B.
```

The exact finite-temperature temperature equation retains its composition and adiabatic terms;
the cold/frozen controlled contract does not silently promote those omitted terms to zero.
Static/quasistatic number sources use one lapse, `int e^Phi Gamma dV`, while powers use two,
`int e^(2Phi) q dV`. Authority: ADR-0015 sections 3–4 and preflight sections 2–4.

For destruction of an ordinary neutron with event rate density `R_a > 0`, the accepted direct
energy per event is

```text
Q_dir = mu_n,actual - E_esc,fluid.
```

For an occupied-state energy `E_n=E_dep+E_esc,fluid`, this is identically

```text
mu_n,actual - E_esc,fluid = (mu_n,actual - E_n) + E_dep.
```

The Fermi-hole term is therefore **NOT INDEPENDENT**. It is the same chemical/free-energy residual
written microscopically and must not be added on top of an inclusive `mu-E_esc` ledger. Escaping
energy is subtracted exactly once; retained-X transfer is distinguished from energy that leaves the
star. No generic reversible `PdV`, gravitational readjustment, binding-energy change, reference
work, `E_chem_dot`, or effective-mass change is an independent heat term. Authority: ADR-0015
sections 4–5, 9, and 19.

## 6. Ratified structural, chemical, and sign semantics

`t=(partial N_y^eq/partial B)_Omega` is the primary structural object for moving-equilibrium BNV
chemistry. At spin OFF it is derived from the qualified Phase-5B equilibrium-sequence response as
`t_i=B_i/B_B`, with the same domain, surface, error, provenance, and currency requirements.

`k=G_y b/(b^T G_y b)` remains **diagnostic only**. It must not map raw `S_y` to a physical BNV
drive, define the physical sliding null, or reconstruct the changing-`B` individual potentials.
The accepted physical slide is

```text
S_y = t Bdot  =>  sigma = 0  =>  eta = 0
```

for zero initial imbalance and no spin forcing. For the controlled neutron sink, `Bdot<0`,
`sigma=t_l |Bdot|>0`, both initial imbalance drives are negative, the restoring beta response is
capture (`R_l<0`), and `eta_l R_l>=0`. These are sign-contract statements, not a physical BNV rate
or trajectory. Authority: ADR-0015 sections 6–7 and 16.

At fixed current `B`,

```text
E_chem = (1/2) eta^T Z^-1 eta = (1/2) ell^T Z ell >= 0
```

is a state reservoir and **NOT HEAT**. Standard beta chemical heating remains
`L_H=C_(MeV->erg) eta^T R`; the incremental beta thermal effect is
`DeltaP_beta=L_H-DeltaLnu` and may cool. The accepted contract makes no generic claim that
Regime-I BNV cools; the sign depends on the direct partition, beta response, history, comparator,
observable, time window, and uncertainty. Authority: ADR-0015 sections 8–10 and 14.

## 7. Regime-I scope, Regime-II boundary, and input ownership

The ratified scope is the declared controlled abstract-source Regime-I framework: spin OFF,
whole-star diffusive non-superfluid free gas, cold chemical response, frozen standard coefficients,
and an abstract charge-consistent source. Product blocking/accumulation, stress/EOS feedback,
chemical feedback, heat capacity, radiation/conduction/opacity, ordinary weak-rate modification,
and coefficient drift must remain within declared budgets.

Regime-II cooling is permitted **in principle**, but only through an explicitly owned physical
loss/transport model and a named comparator. No Regime-II transport or MixedStar thermal evolution
is accepted as implemented. Failure of the Regime-I hidden-sector mechanical, chemical, or thermal
conditions requires the corresponding Regime-II owner; ordinary coefficient drift alone may instead
require an evolving-background ordinary Regime-I model. Authority: ADR-0015 sections 12–14.

The accepted model-independent/process-dependent split is:

| Class | Accepted ownership |
|---|---|
| A | Equilibrium sequence and total `Bdot`; determines `E_eq(B)`, `mu_B`, `t`, and structural state |
| B | Ordinary stoichiometry/source; determines raw `S_y`, `sigma`, and the standard beta drive |
| C | Event energy partition; determines `E_esc,fluid`, deposition, and the direct thermal value |
| D | Product fate; distinguishes escape, ordinary thermalization, retained inert sector, and interacting sector |
| E | Hidden-sector state/interactions; owns accumulated-state, transport, radiation, chemical, and ordinary-weak feedback |
| F | External inflow; owns incoming energy, charges, and angular momentum |

Sequence and chemistry are generic conditional on A+B. Absolute thermal predictions require
process-dependent C–F inputs as applicable. A rate interface may not hide an energy efficiency.
Authority: ADR-0015 section 11 and preflight section 12.

## 8. Accepted contract requirements and source-limited items

BNV-14 through BNV-25 in the preflight section 15 are accepted as ADR-0015 contract requirements.
They govern the cold bracket and finite-temperature floor, moving-reference/t semantics, common
domain and provenance, explicit forward-seam narrowing, individual-potential identity, chemical
storage, once-only energy representation, product fate, non-superfluid beta bounds/QSS
reachability, common energy zero, Cowling scope, and matched-control/sign rules. They are not
separately numbered invariant-register entries and do not claim implementation.

Remaining source-limited items are: a general rotating multispecies stellar first-law theorem;
realistic-EOS `t` and Cowling-error closure; detailed physical channel kinematics, escape, rates,
and rate-weighted hole energy; realistic `E_esc`; authenticated realistic A18 response; superfluid
extensions; changing-background second-order potentials, heat-capacity and boundary terms;
quantitative hidden-sector transport/screening and Regime-II dynamics; and complete independent
reproduction of the diagnostic `G_true` matrix. These limitations do not block the declared
static, spin-off, abstract-source Regime-I contract and must not be represented as solved.
Authority: ADR-0015 section 17 and preflight sections 17–19.

## 9. Explicit implementation and artifact exclusions

This ratification performs and authorizes none of the following:

- physical BNV lifetime, cross section, operator coefficient, coupling, or density-dependent rate;
- `n->chi-gamma` implementation or rate model;
- realistic `E_esc` or detailed channel kinematics/transport;
- realistic A18 physics or FR2005 realistic normalization;
- superfluid extensions;
- Regime-II transport or MixedStar thermal evolution;
- production BNV source, energy partition, fate model, driver, or coupling;
- BNV trajectory or thermal-sign scan;
- Phase-6 numerical candidate or governed baseline.

The first controlled P0/P1/P2 definitions and all required oracles are contract requirements for a
future separately authorized implementation task. They are not executable products in this task.

## 10. Documentation-only validation boundary

The permanent change is restricted to ADR/status/ratification documentation. Validation must prove
that production source, tests, all 11 governed baselines, Phase-5B/C/D governed artifacts,
historical candidates, the 33 protected upstream paths, EOS/data, and literature are unchanged.
Repository-standard lightweight documentation/governance checks are sufficient; no multi-hour
Phase-5 trajectory rerun is required because no governed numerical input or implementation changes.

No numerical Phase-6 artifact exists. This is scientific-contract/preflight ratification only.

## 11. Final ratified status

| Item | Status |
|---|---|
| ADR-0015 | **ACCEPTED / HUMAN-RATIFIED** |
| Phase-6A-0 | **PREFLIGHT COMPLETE / INDEPENDENTLY REVIEWED / HUMAN-RATIFIED** |
| Open-system BNV thermal contract | **GOVERNED for declared controlled abstract-source Regime I** |
| BNV production implementation | **NOT BEGUN** |
| Physical BNV rate | **NOT SELECTED** |
| `n->chi-gamma` | **NOT IMPLEMENTED** |
| Realistic A18 | **NOT BEGUN** |
| Regime-II / MixedStar thermal evolution | **NOT BEGUN** |
| Phase-6 numerical baseline | **NONE** |

Canonical integration, if performed, must be a fast-forward of exact canonical entry
`0a7418aecb7314cfa472a78f1faf477be8456a94` to `PHASE6A0_RATIFICATION_SHA`, with no merge commit,
squash, cherry-pick, rebase, or force push.
