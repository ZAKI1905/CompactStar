# Track-R Structure-1 — human-owner ratification

- **Date:** 2026-09-05
- **Implementation:** `3b0499e3cb5fb39db63d77fb136b4098b8084b4f`
- **Independent review:** `c671170aec5b4d0eded5d65b71332bb0c68dc0cf`

**Status:** **TRACK-R STRUCTURE-1 HUMAN-RATIFIED WITH QUALIFIED CLAIM.**

**HUMAN OWNER RATIFICATION: TRACK-R STRUCTURE-1 LEVELS 1 AND 2 ARE ACCEPTED.**

**THE SOURCE-QUALIFIED PRE-SIGMA "M_MAX" SEMANTICS ARE NOT RATIFIED AS A
CONTINUOUS-SUPREMUM REPRODUCTION.**

## Authority and scope

The scientific source authority is the authenticated FR2005 primary PDF,
`literature/rotochemical/2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf`,
SHA-256 `f184d7d1d7030b61a021eb5c7ac14b1f1b30c7ea69e9d53473d153cfb069ea88`.
Its Table 1 free-gas row, footnote c, whole-star free-gas statement, central-energy-density
meaning, and radius definitions are transcribed and authenticated in
`docs/validation/TRACKR_FREEGAS_WHOLESTAR_STRUCTURE1_INDEPENDENT_REVIEW.md:47-98`.

The owner authority is the human project owner's explicit 2026-09-05 instruction initiating this
record. The owner ratifies the qualified conclusion supported by the implementation and the
independent review; the review itself remained an AI-authored candidate and did not confer
ratification (`docs/validation/TRACKR_FREEGAS_WHOLESTAR_STRUCTURE1_INDEPENDENT_REVIEW.md:902-922`).
This is a documentation and canonical-integration decision. It changes no production result,
test, baseline, source model, numerical method, or scientific scope.

## Ratified claim

**LEVEL 1 — ACCEPTED:**

The Track-R FR2005 static free-gas structure implementation is numerically and scientifically
validated within its stated source domain. The independent review established this level using
analytic re-derivations, high-precision checks, convergence evidence, and three-way structural
agreement (`docs/validation/TRACKR_FREEGAS_WHOLESTAR_STRUCTURE1_INDEPENDENT_REVIEW.md:799-809`).

**LEVEL 2 — ACCEPTED:**

A source-compatible central state at the source's own printed
`rho_c = 1.10e15 g cm^-3` reproduces the printed FR2005 free-gas Table-1 values
`M = 0.62 M_sun`, `R = 12.77 km`, and `R_infinity = 13.80 km` simultaneously,
using one common central state, with no source bin widening, no mass fitting, no EOS tuning, and
independent TOV validation (`docs/validation/TRACKR_FREEGAS_WHOLESTAR_STRUCTURE1_INDEPENDENT_REVIEW.md:469-538`,
`:810-816`). The printed `rho_c` is the source-constrained central input, not a fourth independent
prediction.

**LEVEL 3 — NOT RATIFIED:**

The source-qualified meaning of the Table-1 label `M_max = "maximum mass before appearance of
Sigma-minus hyperons"` is not claimed to have been reproduced as a continuous pre-Sigma
supremum. Level 2 does not imply Level 3
(`docs/validation/TRACKR_FREEGAS_WHOLESTAR_STRUCTURE1_INDEPENDENT_REVIEW.md:817-827`).

## Common-state numerical evidence

| Quantity | Ratified midpoint | FR2005 printed value |
|---|---:|---:|
| `rho_c` | `1.10e15 g cm^-3` | `1.10e15 g cm^-3` |
| `n_B,c` | `~0.6049252413288511 fm^-3` | not printed |
| `M` | `~0.6236355692 M_sun` | `0.62 M_sun` |
| `R_0` | `~12.768154901 km` | `12.77 km` |
| `R_infinity,0` | `~13.8024 km` | `13.80 km` |

All three computed observables round into the published FR2005 row at this one common central
state. The independent central-state inversion agrees at about `4e-18` relative, the independent
direct-vacuum TOV result is `M = 0.623635569157 M_sun` and
`R_0 = 12.7681549010 km`, and both retained apparent-radius conventions round to `13.80 km`
(`docs/validation/TRACKR_FREEGAS_WHOLESTAR_STRUCTURE1_INDEPENDENT_REVIEW.md:469-503`). The
simultaneous matching interval is only about 22.4% of the printed density bin, with the printed
central density near its center; no target-mass inversion or separate state per observable was
used (`docs/validation/TRACKR_FREEGAS_WHOLESTAR_STRUCTURE1_INDEPENDENT_REVIEW.md:505-538`).

## Exact `M_max` exclusion and source ambiguity

The independently determined one-sided sub-Sigma limit of the reconstructed model is recorded
solely to explain why Level 3 is not claimed:

| Limit quantity | Approximate value | Value at source precision |
|---|---:|---:|
| `M` | `0.6259051 M_sun` | `0.63 M_sun` |
| `rho_c` | `1.123745e15 g cm^-3` | `1.12e15 g cm^-3` |
| `R_0` | `12.7051 km` | `12.71 km` |
| `R_infinity,0` | `13.744 km` | `13.74 km` |

These values do not correspond to the published `0.62 / 1.10 / 12.77 / 13.80` row. The mass is
still increasing through the sampled sub-Sigma branch; under a continuous-supremum reading the
source label is inconsistent with all four printed entries
(`docs/validation/TRACKR_FREEGAS_WHOLESTAR_STRUCTURE1_INDEPENDENT_REVIEW.md:540-573`). The exact
FR2005 retained mesh/full-precision input and the selection/retention semantics for its quoted
`M_max` model are unpublished. They remain **SOURCE AMBIGUITY**, not an invitation to infer or
fit a rule (`docs/validation/TRACKR_FREEGAS_WHOLESTAR_STRUCTURE1_INDEPENDENT_REVIEW.md:732-746`,
`:830-889`).

## Independent TOV and surface evidence

The repository oracle uses a separately coded enthalpy TOV formulation, direct phase-space EOS,
independent center expansion, and a distinct integrator. The independent review additionally
implemented a third TOV solve in `(r,m)` and obtained three-way agreement within
`9.4e-11 M_sun` and `2.4e-9 km`
(`docs/validation/TRACKR_FREEGAS_WHOLESTAR_STRUCTURE1_INDEPENDENT_REVIEW.md:345-419`).

The zero-pressure tail is kept distinct from the canonical finite-cutoff surface. At the accepted
cutoff the independently integrated vacuum radius lies within
`[12.7681545694, 12.7681549024] km`; the omitted tail is about 1.980 m and the bracket width is
`3.330e-7 km`. The bracket midpoint is biased low by about `1.65e-7 km`, but the true radius is
inside the stated uncertainty and the bias is immaterial to the ratified printed-bin claim
(`docs/validation/TRACKR_FREEGAS_WHOLESTAR_STRUCTURE1_INDEPENDENT_REVIEW.md:421-467`).

## Protected surfaces and accepted prior evidence

The protected TOV, `NStar`, and `StarProfile` surfaces are byte-identical to the starting master;
evolution and BNV code are untouched. All eight governed baselines are unchanged and no ninth
baseline exists (`docs/validation/TRACKR_FREEGAS_WHOLESTAR_STRUCTURE1_INDEPENDENT_REVIEW.md:698-712`).
Accepted prior evidence is focused 8/8 PASS, data-free 28/28 PASS, recorded full 51/51 PASS,
independent third-TOV and high-precision barotrope checks, and independent source review
(`docs/validation/TRACKR_FREEGAS_WHOLESTAR_STRUCTURE1_INDEPENDENT_REVIEW.md:714-730`). No suite is
rerun for this ratification-only change.

## Findings carried forward

**SOURCE AMBIGUITY**

1. Exact FR2005 retained mesh/full-precision input.
2. Source-qualified `M_max` selection semantics.

**MUST FIX / QUANTIFY BEFORE GLOBAL ROTOCHEMICAL RESPONSE**

1. Threshold-localized `d epsilon/dP` interpolation error, approximately 0.410%; harmless to the
   static result because the static TOV path does not consume it, but it must be quantified for
   future consumers.
2. INV-09 remains **INTENDED BUT UNVERIFIED**.
3. INV-11 remains **UNRESOLVED**.

**NONBLOCKING HARDENING**

1. Tail-bound midpoint is biased low, though well within uncertainty.
2. Missing executable oracle-inside-tail-bracket assertion.
3. S8 endpoint convergence assertion is tautological.
4. S12 coverage gap on floor-sweep and high-branch rows.
5. Pressure closed-form cancellation near `x ~ 0.1`.
6. Near-onset neutron-density error wording is overly optimistic/imprecise.
7. `finest_midpoint_oracle_disagreement` is an imprecise field name.
8. Adjacent-double threshold ambiguity.

**PRE-EXISTING TECHNICAL DEBT**

1. GSL/Zaki solar-mass convention offset in `R_infinity`.

This ranking and its numerical details are preserved from the independent review
(`docs/validation/TRACKR_FREEGAS_WHOLESTAR_STRUCTURE1_INDEPENDENT_REVIEW.md:732-797`). None of
these findings is fixed in this task.

## Next-stage restrictions

This ratification does not resolve INV-09 or INV-11 and does not authorize particle-number
spin-down response, `Btilde`, paper `Z/W`, global rotochemical coefficients, evolution, or BNV.
It does not convert the common-state Level-2 reproduction into a source-qualified Level-3
maximum-mass claim. Any later work must preserve the source ambiguity and independently resolve
the threshold-localized derivative requirement before relying on `dEdP` in a global response.

**Canonical scientific statement:**

> **TRACK-R STRUCTURE-1 HUMAN-RATIFIED WITH QUALIFIED CLAIM —
> STATIC FREE-GAS STRUCTURE AND COMMON-STATE FR2005 TABLE-1 NUMBERS VERIFIED.**
> FR2005 Table-1 common-state static numbers are reproduced; the source-qualified `M_max`
> selection semantics remain unresolved.
