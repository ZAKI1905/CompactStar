# Relativistic unit boundary — B3c scientific-semantic adjudication

**OWNER ACCEPTED — DISPOSITION B**

**B3c HISTORICAL PREMISE SUPERSEDED BY ADR-0012 —**
**REPLACED WITH COMMON-LIMIT / ERROR-ENVELOPE VALIDATION.**

This record applies to canonical SHA
`1d22dd1f5a0d1afa18c4cedebb36b28fdae49df4`, candidate branch
`fix/relativistic-unit-boundary-a1`, and worktree
`/Users/keeper/Documents/CompactStar/worktrees/CompactStar-relativity-units-a1`.
The human project owner accepts the independent Unit-1B3c reviewer disposition B. Historical
validation documents remain preserved as records of what was believed when they were written.

## Original B3c and superseded premise

The original B3c required

```
worst(production/reference relative residual) /
best(production/reference relative residual) < 3
```

It interpreted an approximately constant residual across the radial-resolution ladder as proof
that the background, rather than the Hartle solver, owned the grid sensitivity. That premise was
characterized from the historical mixed-unit state; it was not derived from the Hartle equation.
It passes the known defective pre-A1 representation and fails the coherent A1 representation
because ADR-0012 removes the dominant resolution-independent offset and the corrected residual
approaches zero. Equivalent discretizations are expected to approach a common continuum limit;
they are not required to retain a nonzero, approximately constant disagreement.

The authenticated residual ladders at radial resolutions
`2500 / 5000 / 10000 / 20000 / 40000` are:

| Representation | 2500 | 5000 | 10000 | 20000 | 40000 | Old B3c |
|---|---:|---:|---:|---:|---:|---|
| Pre-A1 | `+4.0274e-5` | `+2.2101e-5` | `+2.1462e-5` | `+1.4589e-5` | `+1.6826e-5` | PASS |
| A1 | `+2.3567e-5` | `+5.4022e-6` | `+4.7504e-6` | `-2.1184e-6` | `+1.1908e-7` | FAIL |

A1 removes an approximately constant `S = 1.670623e-5` offset across the ladder. The exact
analytic-background production/reference result remains approximately `9.455e-9`, and the
independent reference numerical floor remains approximately `1e-15`. No first-order Hartle
physics defect was identified.

## Accepted B3c-prime definition

For every resolution `k`, define

```
res_k = abs(I_prod,k - I_ref,k) / abs(I_ref,k).
```

For every row except the finest, define the observed production background-resolution scale

```
s_k = abs(I_prod,k - I_prod,finest) / abs(I_prod,finest).
```

For the finest row, use the previous resolution:

```
s_finest = abs(I_prod,finest - I_prod,previous) / abs(I_prod,finest).
```

Every `s_k` must be finite and positive, and every resolution must satisfy

```
res_k <= alpha * s_k,    alpha = 1/2.
```

The owner accepts `alpha = 1/2` because subdominance requires a factor-of-two separation. It was
not selected by rounding the observed A1 ratio. The old worst/best residual assertion is removed
and is not retained as a hidden secondary gate. B3a, B3b's absolute `1e-3` accuracy gate, and B4
remain unchanged.

## Detector evidence

The detector below uses the authenticated residual ladders and their recorded production-I
values. Ratios reflect the precision retained in those records; they are evidence, not hard-coded
test targets.

| Resolution | Pre-A1 `res_k/s_k` | Pre result | A1 `res_k/s_k` | A1 result |
|---:|---:|---|---:|---|
| 2500 | `0.1970007944` | PASS | `0.1152762906` | PASS |
| 5000 | `0.4771875920` | PASS | `0.1166237775` | PASS |
| 10000 | `0.5333107151` | FAIL | `0.1180448810` | PASS |
| 20000 | `0.7517349448` | FAIL | `0.1091758706` | PASS |
| 40000 | `0.8670020002` | FAIL | `0.0061370198` | PASS |

Thus the detector reverses the defective historical behavior:

| Criterion | Pre-A1 | A1 | Maximum `res_k/s_k` |
|---|---|---|---:|
| Original B3c | PASS | FAIL | not applicable |
| Accepted B3c-prime | FAIL | PASS | pre-A1 `0.8670020002`; A1 `0.1180448810` |

## Focused U8 revalidation

On 2026-09-06, the six existing first-order Hartle tests were rebuilt and run serially from the
external Unit-1 build directory. No U9--U14 or full-suite test was run.

- `hartle_moment_inertia_analytic`: PASS;
- `hartle_normalization_contract`: PASS;
- `hartle_first_order_physics_analytic`: PASS;
- `hartle_moment_inertia_cmf`: PASS;
- `hartle_normalization_cmf`: PASS; and
- `hartle_first_order_physics_cmf`: PASS.

Raw result: **6/6 PASS, 0 FAIL**, total CTest time `14.17 s`. Every reported assertion passed.
The analytic-background production/reference relative difference remains `9.455e-9`.
B3a, B3b, the accepted B3c-prime, and B4 all pass. The live B3c-prime diagnostics are:

| Resolution | `res_k` | `s_k` | `res_k/s_k` |
|---:|---:|---:|---:|
| 2500 | `2.35659311e-5` | `2.04439902e-4` | `1.15270702e-1` |
| 5000 | `5.39846240e-6` | `4.63193476e-5` | `1.16548757e-1` |
| 10000 | `4.75026488e-6` | `4.02435597e-5` | `1.18037890e-1` |
| 20000 | `2.11912772e-6` | `1.94056205e-5` | `1.09201750e-1` |
| 40000 | `1.18050916e-7` | `1.94056205e-5` | `6.08333630e-3` |

**UNIT-1 U8 FIRST-ORDER HARTLE REVALIDATION PASSED —**
**A1 CANDIDATE MAY RESUME U9-U14.**

This status means only that the B3c semantic defect is resolved and the A1 first-order candidate
passes the accepted focused validation. Unit-1 remains incomplete; A1 is not ratified; baseline
supersession is not authorized; ADR-0011 Q4 remains pending.

## Scope and claim boundary

Only the B3c scientific meaning and its directly required diagnostic output change in
`tests/rotation/hartle_moment_inertia_cmf.cpp`. B1, B2, B3a, B3b, B4, the radial grid, EOS,
production solver, independent reference solver, and all production source remain unchanged.

**THE HARTLE EQUATION, PRODUCTION SOLVER, INDEPENDENT REFERENCE SOLVER,**
**AND B3b ABSOLUTE ACCURACY GATE ARE UNCHANGED.**

**THIS TEST-SEMANTICS CHANGE DOES NOT BY ITSELF VALIDATE UNIT-1, RATIFY A1, SUPERSEDE A**
**BASELINE, OR CLOSE ADR-0011 Q4.**

INV-09 remains unresolved. INV-11 remains unresolved. Phase-5B particle-number production,
Btilde, paper Z/W, evolution, BNV, and governed-baseline supersession remain outside this task.
