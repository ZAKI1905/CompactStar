# ADR-0001 — Species profile semantics: number densities or fractions?

| Field | Value |
|---|---|
| **Status** | **ACCEPTED** |
| **Date proposed** | 2026-08-31 |
| **Date accepted** | 2026-08-31 |
| **Authority** | **Project-owner adjudication.** The owner supplied the EOS schema contract as domain authority; this is not an inference from repository evidence. |
| **Selected alternative** | **Alternative 1 (as amended)** — stored species values are dimensionless fractions `Y_i` |
| **Change class** | scientific-semantic |
| **Governing authority** | This ADR is now the normative authority for species semantics |
| **Affected invariants** | INV-01 (primary), INV-14, INV-16 |
| **Unblocks** | The species-semantics prerequisite of Phase 5. Phase 5 remains blocked by its other prerequisites. |

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

### Why this could not be resolved from the repository alone

*(Historical: this is why the ADR was raised. Resolved by owner adjudication — see Decision.)*

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

**Alternative 1, as amended below.** The canonical CompactStar EOS and profile convention is:

- `n_B(r)` — **total baryon number density**, in fm⁻³. Stored in
  `StarProfile::Column::BaryonDensity`.
- `Y_i(r) = n_i(r) / n_B(r)` — **dimensionless species fraction**. This is what a
  `StarProfile` species column stores.
- `n_i(r) = Y_i(r) · n_B(r)` — **physical species number density**, in fm⁻³. Always derived
  explicitly at the point of use; never stored.

### Amendment to Alternative 1 — no normalization on import

Alternative 1 as originally drafted suggested that `TOVSolver` might need to normalize imported
species values. **That is wrong and is explicitly rejected.**

> The authoritative EOS schema supplies the extra species columns as dimensionless fractions
> `Y_i`. `TOVSolver` and `NStar` **preserve** those values. Consumers requiring physical number
> density must explicitly construct `n_i = Y_i · n_B`.

No division, rescaling, or normalization is to be introduced anywhere in the ingestion, TOV, or
profile-construction path. The input table already carries `n_B` in its baryon-density column and
`Y_i` in its extra composition columns; the pipeline's job is to carry them through unchanged.

### Verified data flow (authenticated at `9f70f14`)

| Stage | Evidence | Behavior |
|---|---|---|
| EOS read — primary | `Core/src/TOVSolver.cpp:545,547` | `rho_val` → `eos_tab.rho` (this is `n_B`) |
| EOS read — extras | `Core/src/TOVSolver.cpp:523-524` | `eos_tab.rho_i.push_back({})` + `AddExtraLabels(...)` |
| EOS read — extras | `Core/src/TOVSolver.cpp:552` | `eos_tab.rho_i[i-3].push_back(std::atof(...))` — **copied verbatim** |
| Normalization on import | `TOVSolver.cpp:515-560` | **None present. Correct — none is wanted.** |
| Profile — baryon density | `Core/src/NStar.cpp:199` | `radial[idx_nb].PushBack(tp.rho)` |
| Profile — species | `Core/src/NStar.cpp:236` | `tp.rho_i[j]` copied directly into the species column |

The pipeline already conforms to the accepted contract.

## Terminology (normative)

These meanings are binding for all future code, comments, and documentation.

| Term | Meaning | Units |
|---|---|---|
| `nB`, `n_B` | Total baryon number density | fm⁻³ |
| `Y_i` | Dimensionless species fraction, `Y_i = n_i / n_B` | — |
| `n_i` | Physical species number density, `n_i = Y_i · n_B` | fm⁻³ |
| **species column** | A `StarProfile` species column means **`Y_i`, never `n_i`** | — |

This applies identically whether the species label is human-readable (`"n"`, `"p"`, `"e"`,
`"mu"`) or a CompOSE numeric particle code (`"10"`, `"11"`, `"0"`, …).

The existing name `rho_i` is **semantically misleading** under this convention — it names a
density but holds a fraction. See Consequences.

## Implementation conformance (observed, not repaired here)

### Currently consistent with the accepted invariant

| Consumer | Evidence | Behavior |
|---|---|---|
| Direct-Urca number-density reconstruction | `Physics/Evolution/src/StarContext.cpp:544-546` | `nn = Yn*nB; np = Yp*nB; ne = Ye*nB;` — carries the explicit in-code comment *"Convert fractions to number densities in fm^-3."* |
| Charge-fraction construction | `Physics/Evolution/src/StarContext.cpp:691-696` | `yq += t.q * (*(t.col))[i];` — sums `q_i Y_i`, correctly yielding a dimensionless `Y_q` with no `n_B` division |
| Legacy per-species density reconstruction | `Core/src/TOVSolver.cpp:1971` | `b_den = ds_vis[rho_i_v_idx[7]] * ds_vis[rho_idx]` — species column × baryon density, correctly forming a density |

### Currently inconsistent

| Component | Evidence | Nature |
|---|---|---|
| `RotochemicalCache` | `Physics/Evolution/src/RotochemicalCache.cpp:147` passes `prof.GetSpeciesPtr(label)` — the raw stored `Y_i` column — directly into `ComputeEnclosedNumber` (`:25-44`) and `ComputeStructuralDerivative` (`:47-104`), both of which document the parameter as *"Species number density column (fm^-3)"* (`RotochemicalCache.hpp:116,133`) and integrate it as such. **No `× n_B` is applied anywhere.** | Uses `Y_i` where `n_i` is required, for `N_i`, `A_i`, and `B_i` |

**This is not a newly introduced regression.** It is pre-existing, uncompiled, unvalidated
candidate code from `675b4a9` — absent from every CMake source list, with `Build()` never called.
It is recorded here as an implementation nonconformance to be corrected in Phase 5, and is
deliberately **not repaired by this ADR**.

## Alternatives

*Retained for the record. Alternative 1, as amended in the Decision, was selected.*

### Alternative 1 — Species columns are canonical **fractions** `Y_i = n_i / n_B` ✅ **SELECTED**

- **Statement.** The stored column is dimensionless. Number density is always derived as
  `n_i = Y_i · n_B` at the point of use.
- **Required code changes.** Correct the Doxygen at `StarProfile.hpp:45,296`, which currently
  describes species columns as densities. Correct `RotochemicalCache` to construct
  `n_i = Y_i · n_B` before its species integrals. Add a typed accessor
  `GetSpeciesFraction(label)` and a derived `GetSpeciesNumberDensity(label)`, and retire the
  ambiguous `GetSpecies`.
  **Amended:** the original draft also proposed verifying that the EOS import normalizes.
  **This is rejected** — the schema supplies `Y_i` directly and the import must preserve it.
  Verification instead confirmed that no normalization is present, which is the correct behavior
  (`TOVSolver.cpp:552`).
- **Migration risk.** Low. Matches what every live compiled consumer already assumes, so
  **no currently compiled behavior changes.** The original draft's stated risk — that the import
  might not normalize — is resolved: it does not, and must not.
- **Validation needed.** On a reference profile, confirm `Σ_i q_i Y_i ≈ 0` (charge neutrality)
  and that `Y_n + Y_p ≈ 1` for npeμ matter. Then re-derive the DU boundary radius and compare
  against a published threshold density. These belong to the validation baseline (roadmap
  **Phase 2B**; renumbered from Phase 2 by ADR-0002).
- **Implications for existing outputs.** **None invalidated.** Every compiled consumer already
  conforms; the sole nonconformant component is not compiled and has never produced output.

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

Now in force:

- **INV-01 is resolved.** It moves from UNRESOLVED to **GOVERNED (ACCEPTED)** in
  `docs/SCIENTIFIC_INVARIANTS.md`, citing this ADR as the normative authority. The storage
  contract is settled; **implementation conformance is not** — see above.
- **INV-16 (Direct-Urca)** must document that its Fermi momenta are built from
  `n_n = Y_n n_B`, `n_p = Y_p n_B`, `n_e = Y_e n_B`. The existing implementation already does
  this correctly.
- **INV-14 (baryon number)** is unaffected in form: it integrates `n_B` directly, not a species
  column. Any *per-species* enclosed-number integral must apply `n_i = Y_i n_B`.
- **The species-semantics prerequisite of Phase 5 is satisfied.** Phase 5 remains blocked by all
  of its other prerequisites — INV-07, INV-11, out-of-equilibrium weak rates, and the CMake
  reachability gap.
- **No source change is authorized by this ADR.** The corrections below are recorded as future
  implementation tasks.

### Required future source corrections (recorded, not performed)

1. **`Core/StarProfile.hpp:45,296`** — Doxygen describes species columns as number densities.
   Must be corrected to dimensionless fractions `Y_i`. *(When edited, preserve Doxygen format.)*
2. **`rho_i` naming** — semantically ambiguous: the name says density, the content is a fraction.
   Should eventually be renamed or wrapped behind typed accessors.
3. **Typed accessors** — a preferred future API distinguishes `GetSpeciesFraction(label)` from
   a derived `GetSpeciesNumberDensity(label)`, retiring the ambiguous `GetSpecies`. Removing the
   ambiguous accessor converts any future misuse from a silent numerical error into a build error.
4. **`RotochemicalCache`** — must construct `n_i = Y_i · n_B` before its `N_i`, `A_i`, and `B_i`
   species number-density integrals (`RotochemicalCache.cpp:147`, `:25-44`, `:47-104`).
5. **Unit-contract machinery** — `Diagnostics/UnitContract` should eventually encode the
   fraction/density distinction machine-readably, so conformance is checked rather than
   conventional. It is currently populated by no producer.

## Validation

The decision itself required no numerical validation: it is a schema contract supplied by the
owner, and the compiled pipeline was verified to already conform.

Validation is still required before the **future source corrections** land, and belongs to the
validation baseline — roadmap **Phase 2B**, renumbered from Phase 2 by ADR-0002:

1. A reference stellar profile committed as a fixture, with the EOS table version recorded.
2. Charge neutrality `Σ_i q_i Y_i ≈ 0` demonstrated on that fixture.
3. `Y_n + Y_p ≈ 1` for npeμ matter, confirming the columns really are fractions in the tables in use.
4. The Direct-Urca onset density compared against an independent published value for the same EOS.
5. A passive cooling curve captured as a regression baseline **before** any ADR-0001 source
   correction lands. Those corrections are Phase-5 work, so the Phase-2B baseline precedes them.
   Note that the baseline itself is now preceded by the Phase-2A heat-capacity conformance work
   required by **ADR-0002**; it must not be captured while `PhotonCooling` still divides by a
   constant.

## Provenance

Conflict identified during Phase-0 reconnaissance (2026-08-31) against commit `d91c31b`;
re-confirmed present in `9f70f14`. Alternatives drafted by an AI agent under `GOVERNANCE.md` §5
and `AGENTS.md`.

**Adjudicated by the project owner on 2026-08-31**, who supplied the authoritative EOS schema
contract — domain authority that could not be recovered from repository evidence, since both
readings were internally consistent within their own layers. The agent then independently
verified the ingestion, profile-construction, and consumer data flow against `9f70f14` and
recorded conformance, including the single nonconformant component.

Per `GOVERNANCE.md` §5, acceptance of this ADR ratifies **only** the species-semantics contract.
It does not ratify the `RotochemicalCache` or Hartle O(Ω²) candidate code, and does not confer
accepted status on any other DRAFT governance document.
