# Relativistic Unit Boundary Surface Assertion Adjudication

**OWNER-ACCEPTED DISPOSITION C — SURFACE EXACT-EQUALITY SEMANTICS TOO STRONG;
REPLACE MASS BIT-EQUALITY WITH A DERIVED MACHINE-PRECISION CONVERSION CHECK.**

**THIS IS A TEST-SEMANTICS CORRECTION.
IT DOES NOT ALTER THE ADR-0012 A1 PRODUCTION IMPLEMENTATION.**

**THIS DECISION DOES NOT RATIFY UNIT-1, SUPERSEDE ANY GOVERNED BASELINE,
CLOSE ADR-0011 Q4, RESOLVE INV-09, OR RESOLVE INV-11.**

## 1. Authenticated scope and authority

Canonical SHA and candidate HEAD: `1d22dd1f5a0d1afa18c4cedebb36b28fdae49df4`.
Candidate branch: `fix/relativistic-unit-boundary-a1`.
Candidate worktree:
`/Users/keeper/Documents/CompactStar/worktrees/CompactStar-relativity-units-a1`.

The owner accepts this adjudication under ADR-0012's same-spacetime geometric mass rule
(`docs/adr/ADR-0012-relativistic-unit-boundary.md:56`) and preserved literal public-mass
semantics (`docs/adr/ADR-0012-relativistic-unit-boundary.md:85`). ADR-0012 explicitly requires
the geometric mass relation

```
m_km = G_TOV * m_grams / c_TOV^2 / 1e5
```

and does not require two independently ordered IEEE-754 implementations of that relation to be
bit-identical (`docs/adr/ADR-0012-relativistic-unit-boundary.md:61`). The accepted A1 scope also
forbids a TOV numerical change (`docs/adr/ADR-0012-relativistic-unit-boundary.md:119`).

Classification: owner-approved scientific-semantic test/validation correction and
documentation. There is no production, numerical-method, generated-artifact, or governed-
baseline change in this adjudication.

## 2. Original assertion and historical rationale

The historical combined assertion was:

```cpp
Require(r == p.back().r && m == p.back().m * relativity_fixture::solar_km,
        "surface profile coordinate");
```

It treated both fields as exact representation copies. That rationale remains valid for radius:
the surface profile radius is copied directly from `TOVPoint::r`. It is not valid for geometric
mass after ADR-0012 A1 because mass is a derived unit conversion rather than a direct copy.

Before A1, the fixture multiplication happened to reproduce the historical profile-mass bits.
A1 replaced the mixed-convention profile mass with the accepted same-spacetime conversion. The
arithmetic graph therefore changed even though the mathematical relation did not.

## 3. Production and independent fixture expressions

Production evaluates the public ratio through grams and then through the production conversion
factor:

```text
m_prod = (M_public * M_sun_GSL) * (G_TOV / c_TOV^2 / 1e5)
```

The production implementation is `LiteralSolarMassToMassKm`, composed from
`LiteralSolarMassToGrams` and `MassGramsToKm`
(`CompactStar/RelativityUnits.hpp:16`, `CompactStar/RelativityUnits.hpp:22`,
`CompactStar/RelativityUnits.hpp:35`).

The independently written fixture first forms a solar geometric length and then scales it by the
public ratio:

```text
solar_km_fixture = G_TOV * M_sun_GSL / c_TOV^2 / 1e5
m_expect = M_public * solar_km_fixture
```

The fixture definition is intentionally independent of production RelativityUnits
(`tests/relativity/fixture_units.hpp:1`, `tests/relativity/fixture_units.hpp:10`).

Both graphs evaluate the same real-number expression required by ADR-0012:

```text
G_TOV * (M_public * M_sun_GSL) / c_TOV^2 / 1e5
```

A decimal high-precision oracle using the accepted literals gives

```text
G_TOV * M_sun_GSL / c_TOV^2 / 1e5
= 1.476716181892116417835469812898496840892624040210026207253944... km
```

so the two double-precision graphs differ only by rounding order, not by physical convention.
Production is scientifically correct because it preserves `M_public`, reconstructs the literal
GSL solar mass in grams, and applies the accepted `G_TOV/c_TOV^2` conversion belonging to the
same solved spacetime. The already completed U3 checks independently close the geometric TOV
identities near `1e-15`; the surface lapse and metric-lambda identities also pass. Bit equality
between arithmetic graphs is neither an ADR-0012 invariant nor an independent scientific
oracle.

## 4. Independent four-star diagnosis

The owner-accepted independent review evaluated all four HW stars:

| Quantity | Result |
|---|---|
| Surface radius | Exact on 4/4 stars |
| Surface geometric mass | Bit-exact on 3/4 stars |
| Non-exact mass case | Exactly 1 ULP |
| Non-exact relative discrepancy | Approximately `1.648e-16` |
| Surface lapse identity | PASS on 4/4 stars |
| Surface metric-lambda identity | PASS on 4/4 stars |

This isolates the historical failure to an over-strong exact-equality test semantic. It does not
identify an A1 production defect.

## 5. Accepted tolerance derivation

The accepted relative tolerance is exactly:

```cpp
16.0 * std::numeric_limits<double>::epsilon()
```

Each correct expression is a chain of no more than approximately five rounded IEEE-754
operations. Two independently ordered graphs therefore have a first-order separation bound of
approximately

```text
10u = 10 * 2^-53 ~= 1.11e-15 ~= 5 epsilon.
```

The next conservative power-of-two operation-count bound is

```text
16 epsilon ~= 3.55e-15.
```

This is the existing Unit-1 relativistic-unit-test convention. The observed approximately
`1.65e-16` discrepancy was not used to select the bound. The bound must not be tightened,
widened, or replaced by an empirical tolerance.

## 6. Mutation detector authority

The owner-accepted independent mutations remain far outside the authorized gate:

| Mutation | Approximate relative separation | Result against 16 epsilon |
|---|---:|---|
| Historical Zaki `SUN_M_KM` | `6.172e-5` | Rejected |
| Modern CODATA G | `1.948e-4` | Rejected |
| Different solar-mass convention | `2.564e-4` | Rejected |
| Wrong power of c | `O(1)` | Rejected |
| Missing or doubled `1e5` | `O(1)` or larger | Rejected |
| Deliberate `1e-12` mass perturbation | `1e-12` | Rejected |

The retained test therefore detects the unique local defect class:

> TOVPoint public mass was not wired through the accepted ADR-0012 geometric mass conversion at
> the NStar surface node.

The mass check is not redundant and must not be removed.

## 7. Exact owner-approved replacement

Radius remains an exact representation invariant:

```cpp
Require(r == p.back().r,
        "surface profile radius");
```

Mass remains an independent conversion check:

```cpp
const double m_expect =
    p.back().m * relativity_fixture::solar_km;

Require(std::fabs(m / m_expect - 1.0) <=
            16.0 * std::numeric_limits<double>::epsilon(),
        "surface profile geometric mass");
```

The test retains the exact surface radius check, the independently derived coherent geometric
mass check, the surface lapse boundary, and the surface metric-lambda check
(`tests/core/tov_surface_contract.cpp:425`).

## 8. Corrected-test validation

The corrected `tov_surface_consumers_hw` was rebuilt alone and passed 1/1. Its four-star
diagnostics were:

| Target label | `M_public` | Surface radius [km] | Diagnostic mass relative error | Lapse |
|---:|---:|---:|---:|---:|
| 1.0 | `0.45615019716809196` | `18.097286376708052` | `-1.6481826152820721e-16` | `0.9620589463972562` |
| 1.4 | `0.49391427616962408` | `16.427339675185635` | `0` | `0.95456813584985656` |
| 1.6 | `0.5148743217274645` | `15.573380345199396` | `0` | `0.94992422999094217` |
| 2.0 | `0.58249289753018452` | `13.054030086879616` | `0` | `0.93177940174592577` |

The exact-radius, surface-lapse, and metric-lambda assertions pass on all four stars. The
complete authenticated focused selection was then rerun from the beginning, serially with
stop-on-failure: **42 PASS, 0 FAIL**, CTest rc=0, 503.92 s. No data-free or 53-test full suite was
started.

## 9. Preservation boundary

No production file, RelativityUnits definition, NStar implementation, TOVSolver equation,
Geometry equation, surface construction, lapse assertion, lambda assertion, numerical physics,
candidate artifact, or governed baseline is changed by this decision. The six prepared successor
hashes remain candidates only, and the two unchanged TOV candidates remain subject to the
existing byte-identity requirement (`docs/adr/ADR-0012-relativistic-unit-boundary.md:145`).

Passing the corrected surface test and focused set authorizes only resumption of the remaining
authenticated Unit-1 data-free/full-suite validation. It does not complete or ratify Unit-1,
authorize baseline supersession, close ADR-0011 Q4, or resolve INV-09 or INV-11.
