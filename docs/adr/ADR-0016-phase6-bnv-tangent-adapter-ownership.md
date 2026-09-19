# ADR-0016: Phase-6 BNV tangent-adapter ownership

## 1. Status and scope

**Status:** PROPOSED — OWNER RATIFICATION REQUIRED.

**Date proposed:** 2026-09-19.

This is a deliberately narrow architecture/ownership proposal. The human owner has accepted
the recovery plan and authorized the relocation concept, but has not ratified this exact ADR
text. No production relocation, R0 recovery step, trajectory, or numerical requalification
may begin from this proposal alone. The acceptance record is
`docs/validation/PHASE6A1_CONTROLLED_BNV_RECOVERY_ACCEPTANCE.md`.

This ADR changes no scientific equation. ADR-0015 physics, the Phase-5B structural-response
authority, Phase-5C/D mathematics, every existing Phase-5 source byte, the Phase-5D
provenance producer/comparator, and the governed Phase-5D baseline remain unchanged.

## 2. Context

The failed Phase-6A-1 implementation placed its typed equilibrium-baryon tangent adapter at:

```text
CompactStar/Analysis/EquilibriumBaryonTangent.hpp
CompactStar/Analysis/src/EquilibriumBaryonTangent.cpp
```

The governed Phase-5D producer recursively authenticates production files below
`CompactStar/Analysis`. The two new adapter files therefore changed fresh scientific
provenance from 91 to 93 source entries even though no existing governed source byte changed.
BA15 correctly remained **FAIL**
(`docs/validation/PHASE6A1_CONTROLLED_BNV_RECOVERY_PREFLIGHT.md:94-124`).

The conflict is architectural ownership/provenance, not a scientific change. Governance
requires an ADR when ownership or architecture boundaries move, plus a corresponding current-
architecture update (`GOVERNANCE.md:43-57`). Hiding the files from governed provenance would
weaken the Phase-5D contract and is not an acceptable repair.

## 3. Decision proposed

The typed Phase-6 BNV equilibrium-baryon tangent adapter is owned by the Phase-6 BNV module:

```text
header:    CompactStar/Physics/BNV/EquilibriumBaryonTangent.hpp
source:    CompactStar/Physics/BNV/src/EquilibriumBaryonTangent.cpp
type:      CompactStar::Physics::BNV::EquilibriumBaryonTangent
```

It is a typed consumer/view of the governed Phase-5B
`CompactStar::Analysis::EquilibriumSequenceNumberDerivative`. It is not and must not become
a new independent structural-response authority. It performs no new structural solve and
does not independently recompute the sequence derivative.

The proposed path follows the existing Phase-6 BNV public-header/source layout and makes the
adapter's Phase-6 ownership explicit
(`docs/validation/PHASE6A1_CONTROLLED_BNV_RECOVERY_PREFLIGHT.md:126-150`).

## 4. Scientific non-change

The governed tangent remains

```text
t = (partial N_y^eq / partial B)_Omega.
```

At `Omega = 0`, the components remain

```text
t_i = B_i / B_B,
B_B = B_n + B_p.
```

The reduced `B_n + B_e + B_mu` expression remains a charge-closure check only. ADR-0011
continues to own the complete-star baryon count, domain-qualified equilibrium sequence, and
complete-star number derivatives; ADR-0015 continues to own the moving-reference seam and
requires the existing Phase-5B derivative, domain/surface policy, numerical errors, and
currentness (`docs/adr/ADR-0011-particle-number-structural-response.md:44-84`,
`docs/adr/ADR-0015-bnv-open-system-thermal-ledger.md:137-203`).

Relocation must preserve:

- the same complete-star and declared-domain convention;
- the same Phase-5B derivative authority;
- the same component values and propagated numerical errors;
- the same `b^T t = 1` closure and certified representation;
- the same star, sequence, domain, provenance, revision, and currentness checks;
- the same BA2, BA3, BA4, and BA5 requirements; and
- the same M2, M17, and M18 falsifiers.

It is forbidden to substitute `k` for `t`, map the raw source through `G_y`, independently
recompute `t`, or add a structural solve to the Phase-6 adapter. ADR-0015 physics, equations,
energy ledger, scope, and exclusions are unchanged.

## 5. Provenance boundary

This decision is not authority to:

- edit the Phase-5D provenance producer;
- add an exclusion or ignore path for Phase-6 files;
- modify Phase-5D comparator semantics;
- regenerate, replace, or modify the governed Phase-5D baseline;
- change Phase-5B code or any existing governed Phase-5 source byte; or
- weaken source-provenance equality.

The architectural repair succeeds only if the relocated Phase-6 consumer leaves fresh
Phase-5D scientific provenance and governed artifact identity exactly equal to the governed
baseline contract. The existing contract requires scientific source-provenance equality
without a generic ignore-path escape
(`docs/validation/PHASE6A1_CONTROLLED_BNV_RECOVERY_PREFLIGHT.md:94-108`,
`docs/validation/PHASE6A1_CONTROLLED_BNV_RECOVERY_PREFLIGHT.md:191-205`).

## 6. Proof obligations before any recovered trajectory

After owner ratification and during a separately authorized recovery implementation, all of
the following are mandatory before any recovered trajectory:

1. Remove only the two failed-branch `CompactStar/Analysis` adapter files.
2. Create the semantically equivalent Phase-6 BNV-owned adapter at the paths in section 3.
3. Obtain focused BA2-BA5 PASS.
4. Obtain M2, M17, and M18 PASS.
5. Obtain the nonzero-eta R18 PASS.
6. Obtain focused BA10a and BA10b PASS.
7. Obtain retained-certificate parsing and currentness PASS.
8. Run the fresh empty-scratch Phase-5D regression and obtain exact governed scientific
   provenance and governed artifact identity.
9. Demonstrate no Phase-5 baseline, producer, or comparator change.

Fixture `t` values and source-to-`sigma` results must remain bit-identical where possible;
all component budgets, closure, and refusal semantics must remain within their already
accepted bounds. Any unexpected semantic value change is a stop condition, not permission to
widen a budget (`docs/validation/PHASE6A1_CONTROLLED_BNV_RECOVERY_PREFLIGHT.md:174-205`,
`docs/validation/PHASE6A1_CONTROLLED_BNV_RECOVERY_PREFLIGHT.md:440-460`).

## 7. Relation to BA12R

This ADR owns only the architectural location and type ownership of the Phase-6 adapter. It
does not redefine numerical convergence. The separately owner-accepted recovery plan owns
BA12R and retains historical BA12 as **FAIL**.

For scope identification only, the accepted hierarchy is BASELINE `rtol=1e-7`,
`atol=(1e-12,1e-18,1e-18)`; REFINED `rtol=1e-9`,
`atol=(1e-14,1e-20,1e-20)`; and ULTRA `rtol=1e-11`,
`atol=(1e-16,1e-22,1e-22)`. The detailed equations, thresholds, floor handling, ledger
requirements, solver refusals, nominal-tier choice, and restricted P2 source/control scope
remain in the accepted recovery plan
(`docs/validation/PHASE6A1_CONTROLLED_BNV_RECOVERY_PREFLIGHT.md:249-419`).

## 8. Rejected alternatives

- **Keep the adapter under `CompactStar/Analysis` and change provenance collection:** rejected;
  it would change the governed Phase-5D producer/contract.
- **Add an ignore path or relax comparator equality:** rejected; it would hide a scientific-
  provenance change rather than repair ownership.
- **Regenerate the Phase-5D baseline:** rejected; the Phase-5D result and source bytes did not
  change.
- **Move authority from Phase-5B into Phase-6:** rejected; the adapter is a typed consumer only.
- **Use `k` or raw `G_y` in place of `t`:** rejected by ADR-0015 and the recovery falsifiers.

## 9. Consequences and ratification gate

If ratified, ADR-0016 permits only the future R0 ownership relocation under the accepted
recovery plan. It does not accept the failed implementation, authorize a trajectory, select a
physical BNV rate/model, create a candidate artifact, or promote a governed baseline.

Until the human owner explicitly ratifies this exact ADR text:

- ADR-0016 remains **PROPOSED**;
- recovery implementation remains **NOT RESUMED**; and
- R0 must not begin.
