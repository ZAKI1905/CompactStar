# ADR-0001 — Species profile semantics: number densities or fractions?

| Field | Value |
|---|---|
| **Status** | **PROPOSED** — awaiting owner adjudication. Not accepted. |
| **Date** | 2026-08-31 |
| **Change class** | scientific-semantic |
| **Governing authority** | None yet — this ADR exists because none exists |
| **Affected invariants** | INV-01 (primary), INV-14, INV-16 |
| **Blocks** | Phase 5 (rotochemical heating); any trusted per-species integral |

## Context

`StarProfile` stores per-species radial columns alongside the eight fixed structural columns,
registered by label (`"n"`, `"p"`, `"e"`, `"mu"`, or CompOSE numeric codes such as `"10"`,
`"11"`, `"0"`) at column index `8 + j`.

**The profile documents these as number densities. Consumers read them as fractions.**

### Evidence that they are documented as number densities

- `CompactStar/Core/StarProfile.hpp:45` — species entries described as number densities.
- `CompactStar/Core/StarProfile.hpp:296` — same, at the `species_idx` declaration.
- The upstream label convention is inherited from the EOS file header in file-column order
  (`CompactStar/Core/src/TOVSolver.cpp:519-526`), where the source table's own units apply.
- `ARCHITECTURE.md` (untracked, `d91c31b`) line 109 repeats the density reading.

### Evidence that consumers treat them as fractions

- `CompactStar/Physics/Evolution/src/StarContext.cpp:544-546` — computes number densities by
  **multiplying the column by n_B**:
  ```
  nn = Yn * nB;
  ```
  This is only dimensionally correct if the column is a dimensionless fraction.
- `CompactStar/Physics/Evolution/src/StarContext.cpp:695` — charge fraction accumulated by
  summing `q_i × column_i` directly, with no `n_B` division:
  ```
  yq += t.q * (*(t.col))[i];
  ```
  A sum over number densities would yield a charge *density*, not `Y_q`.
- `CompactStar/Core/src/TOVSolver.cpp:1971` — forms a product
  `ds_vis[rho_i_v_idx[7]] * ds_vis[rho_idx]`, i.e. species-column × baryon-density column,
  again consistent with fraction × density.

### Why this cannot be resolved from the repository

Both readings are internally consistent within their own layer. The EOS ingestion path could
plausibly deliver either, depending on the CompOSE table variant and any normalization applied
during import — and no assertion, unit test, or unit contract anywhere constrains it. The
`UnitContract` machinery that could have recorded it exists but is never populated: every
producer returns an empty contract
(`CompactStar/Physics/Driver/Thermal/src/NeutrinoCooling.cpp:129-138`,
`.../PhotonCooling.cpp:155-164`).

There is no test suite, no assertion (`grep assert|static_assert` over `CompactStar/` returns no
files), and no reference output to arbitrate against.

### Why it matters

The distinction is a factor of `n_B ≈ 0.1 – 1 fm⁻³` — large enough to change results
substantially, small enough that plots still look physically plausible. It propagates into:

- **Enclosed species numbers** `N_i` (`Physics/Evolution/src/RotochemicalCache.cpp:25-44`) and
  therefore every `A_i`, `B_i`, and `Z_i` rotochemical coefficient.
- **The Direct-Urca threshold** (`Physics/Evolution/src/StarContext.cpp:564-569`), where
  `k_F = (3π² n)^{1/3}` requires a genuine number density. If the column is a fraction and is
  used unconverted, every Fermi momentum — and hence the DU boundary radius — is wrong.
- **Charge neutrality diagnostics** (`StarContext.cpp:636-695`).
- **Per-species baryon accounting** (INV-14).

The DU path already multiplies by `n_B` and so appears to assume fractions; the
`RotochemicalCache` integrand does **not**, and so appears to assume densities. **At least one
of these two is currently wrong**, whichever way the underlying convention actually goes.

## Decision

*Not decided. This ADR is PROPOSED.*

## Alternatives

### Alternative 1 — Species columns are canonical **fractions** `Y_i = n_i / n_B`

- **Statement.** The stored column is dimensionless. Number density is always derived as
  `n_i = Y_i · n_B` at the point of use.
- **Required code changes.** Correct the documentation at `StarProfile.hpp:45,296`. Audit and
  fix `RotochemicalCache::ComputeEnclosedNumber` (`RotochemicalCache.cpp:25-44`) to multiply by
  `n_B` before integrating. Verify the EOS import path actually normalizes
  (`TOVSolver.cpp:519-526`). Add a typed accessor `GetSpeciesFraction(label)` and retire the
  ambiguous `GetSpecies`.
- **Migration risk.** Low-to-moderate. Matches what the majority of live consumers already
  assume, so most compiled behavior is unchanged. The risk is that the EOS import does *not*
  normalize, in which case this alternative silently endorses existing incorrect DU thresholds.
- **Validation needed.** Confirm `Σ_i Y_i` behaves as expected and that `Σ_i q_i Y_i ≈ 0`
  (charge neutrality) on a reference profile. Confirm `Y_n + Y_p ≈ 1` for an npeμ EOS. Then
  re-derive the DU boundary radius and compare against a published threshold density.
- **Implications for existing outputs.** Previously generated DU boundaries and cooling curves
  remain as-is only if the import truly normalizes; otherwise all are invalidated.

### Alternative 2 — Species columns are canonical **number densities** `n_i` [fm⁻³]

- **Statement.** The stored column carries fm⁻³, matching the existing documentation and the
  CompOSE convention for tabulated composition.
- **Required code changes.** Fix `StarContext.cpp:544-546` (remove the `× n_B`),
  `StarContext.cpp:695` (divide by `n_B` to form `Y_q`), and `TOVSolver.cpp:1971`. Verify
  `RotochemicalCache.cpp:25-44` is already correct under this reading. Add
  `GetSpeciesDensity(label)`.
- **Migration risk.** **High.** This changes the numerical behavior of the Direct-Urca mask, the
  charge-fraction cache, and the heat-capacity integral — all of which are on the live thermal
  evolution path — with **no baseline to detect regression**.
- **Validation needed.** Everything in Alternative 1, plus a full re-validation of the passive
  cooling curve, since `Y_q` feeds the heat capacity (`StarContext.cpp:788-817`).
- **Implications for existing outputs.** Every committed result under `main/Test/results/` that
  depends on DU or heat capacity would be invalidated.

### Alternative 3 — Store **both**, behind distinct typed accessors

- **Statement.** `StarProfile` owns fractions `Y_i` as the stored primitive and exposes both
  `GetSpeciesFraction(label)` and `GetSpeciesDensity(label)`, the latter deriving `Y_i · n_B`
  internally. The ambiguous `GetSpecies` is removed so no caller can be unclear.
- **Required code changes.** Largest surface: new accessors on `StarProfile`, migration of all
  ~8 call sites, removal of the ambiguous accessor, and an update to
  `Diagnostics/UnitContract` so the distinction is machine-checkable rather than conventional.
- **Migration risk.** Moderate, but **fails safe** — the compiler locates every call site once
  `GetSpecies` is removed, converting a silent numerical error into a build error.
- **Validation needed.** As Alternative 1, plus confirmation that the derived-density path
  matches the direct one on a reference profile to round-off.
- **Implications for existing outputs.** Depends on which reading the stored primitive adopts;
  the accessor split itself is behavior-preserving once call sites are corrected.

## Consequences

Once accepted:

- INV-01 moves from **UNRESOLVED** to **VERIFIED CURRENT BEHAVIOR** or **PROPOSED**, as
  appropriate, and `docs/SCIENTIFIC_INVARIANTS.md` is updated in the same change.
- INV-14 (baryon-number definition) and INV-16 (Direct-Urca threshold) must be re-checked for
  consistency with the ratified reading.
- Phase 5 becomes unblocked; Phase 2 gains a concrete first validation target.
- Any alternative other than a pure documentation fix is a **scientific-semantic change** and
  therefore requires a validation baseline **before** implementation (`GOVERNANCE.md` §2).

## Validation

No baseline currently exists — this is itself a finding. Minimum evidence before implementing
any alternative:

1. A reference stellar profile committed as a fixture, with the EOS table version recorded.
2. Charge neutrality `Σ_i q_i n_i ≈ 0` demonstrated on that fixture under the ratified reading.
3. The Direct-Urca onset density compared against an independent published value for the same EOS.
4. A passive cooling curve captured as a regression baseline **before** any change lands.

## Provenance

Conflict identified during Phase-0 reconnaissance (2026-08-31) against commit `d91c31b`;
re-confirmed present in `9f70f14`. Alternatives drafted by an AI agent under
`GOVERNANCE.md` §5 and `AGENTS.md`. **No alternative has been selected. Selection requires the
project owner**, who alone holds the domain authority and the knowledge of what the EOS import
was intended to deliver.
