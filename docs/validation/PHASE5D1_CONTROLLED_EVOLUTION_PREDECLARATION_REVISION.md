# Phase-5D-1 controlled evolution predeclaration revision

**Status:** PRE-TRAJECTORY / QUALIFICATION-AUTHORIZED ADDENDUM. This document
preserves rather than rewrites the original predeclaration at
`697ec14610d63640866f3c0aed0a51d567cb037d`.

## Sole revision

Original predeclaration: **radial resolution 10000 — superseded after
pre-result numerical impossibility/refusal**.

Revised production radial resolution: **80000**.

Reason: the complete dependency chain was freshly reconstructed at the fixed
owner-selected 80000 resolution and passed the unchanged governed Phase-5C
G/Q/Z/W, numerical-error, validation-envelope, support/rank/conditioning,
refusal/tail, provenance/currentness, and upstream-regression gates. The
qualification is recorded in
`PHASE5D1_STRUCTURAL_RESOLUTION_QUALIFICATION.md` under controlling plan SHA
`9a74b245014ea1418b24c483dc066aefb7aca14d`.

This is not post-result tuning of a trajectory. The 80000 candidate was fixed
before admissible qualification output, and no trajectory existed when this
amendment was made.

## Everything else remains frozen

Only the production radial resolution changes. In particular:

- `SMe = 1e-51 erg cm^-3 s^-1 K^-8` and
  `SMmu = 2e-51 erg cm^-3 s^-1 K^-8`;
- Me/Mmu enabled; De/Dmu disabled;
- `T0 = 1e8 K` and `eta0 = (0,0) MeV`;
- `B=1e8 G`, `P0=1 ms`, `PPdot=(B/3.2e19)^2`,
  `P(t)=sqrt(P0^2+2 PPdot t)`, `Omega=2pi/P`, and
  `OmegaDot=-2pi PPdot/P^3`;
- run interval `0..1e10 yr` with `1 yr=365.25*86400 s`;
- GSL RKF45 intended `rtol=1e-7` and absolute tolerances
  `(1e-12,1e-18,1e-18)`, with the declared refinement by 100;
- the declared global integration refinement policy, output grid, step/sample
  limits, temperature-cache boundaries, no eta clamp, and stiffness checks;
- the `1e9..1e10 yr` quasi-steady window, eligibility and 1/7 criteria;
- source-function roots and all analytic/oracle/mutation criteria; and
- RE10b relative tolerance `1e-10`.

EOS resolution remains 8192 intervals; `rho_c` remains
`1.10e15 g cm^-3`; the ordinary NStar Track-R Structure-1 physical fixture,
whole-star chemical/structural domains, finite-cut/tail policies, provider and
source authority, benchmark normalizations, and every unchanged Phase-5C goal
remain exactly as originally declared.

## Boundary

This addendum does not resume coupled evolution. It does not restore
`build/phase5d-audit/uncommitted-coupling/`, enter an ODE, run a secular
trajectory, produce an evolution/candidate SHA, change Z/W physics, alter the
benchmark normalization, relax any goal, implement A18, begin BNV, resolve
global INV-11, or authorize a merge. Any archived coupling work may be used in
a later separately authorized task only as an unaudited patch source after the
six retained architecture guards are rechecked.
