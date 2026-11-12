# ADR-0025 — Retained-energy parameter for a controlled BNV competition study

| Field | Value |
|---|---|
| Status | **PROPOSED — owner selection and approval pending** |
| Date | 2026-11-12 |
| Change class | Scientific-semantic if implemented; this submission is documentation only |
| Governing authority | GOVERNANCE §§1–5; ADR-0015 §§3–11; ADR-0016 tangent ownership; ADR-0014 reaction/state contract; INV-15 |
| Affected invariants | INV-11 controlled state/reaction/ledger semantics; INV-15 common Cstar; no broader invariant closure |
| Blocks | Implementation and execution of a retained-energy BNV study |

## Context

Stage B at `49253200a801ad9b336789dffb4e93d535a8e903` has independent disposition
B — CONDITIONAL PASS, with owner ratification pending. It verifies fixed full-
retention P2 inputs only. The combined chemical response is −73.1 µK; gross
reaction heating alone is +5.19 µK. Existing term-suppressed reviewer diagnostics
are not new physical models (`docs/validation/results-first-bnv-pilot/stage-b/review-evidence/REVIEW.md:281`).

The code already separates ordinary source history, energy partition, actual
potentials and product fate. `IDirectBnvEnergyPartition` receives the actual
potential (`CompactStar/Physics/BNV/DirectEnergyLedger.hpp:33`); the ledger forms
q_actual=−nu·mu_actual−Eesc,fluid+Ein and separately q_eq from the same partition
(`CompactStar/Physics/BNV/src/DirectEnergyLedger.cpp:73`). Existing cold P0 and P2
fixtures use Eesc=mu_n,actual and Eesc=0, respectively
(`tests/bnv/energy_partition_fixtures.hpp:28`). The current P2 assembly has a fixed
card and fate (`tests/native/pair.hpp:7`, `:50`). Neither that card nor sealed
Stage B evidence may be repurposed as a variable-retention experiment.

The [scientific proposal](../validation/results-first-bnv-pilot/stage-b/NEXT_PHYSICS_STUDY.md)
derives the power ratios, relaxation criteria and finite-temperature limits from
ADR-0015:144, :211, :233. “Retention fraction” must identify its energy denominator;
scaling the thermal source alone does not define a conservative partition.

## Decision

## Alternatives

### A — Actual-potential residual fraction (recommended for owner consideration)

- **Statement:** choose constant f∈[0,1] for a neutron sink D=−Bdot, with
  Eesc,fluid=Eesc,star=(1−f)mu_n,actual and EX=Ein=0. Keep sigma=D t_l.
  Then Pactual=cD f mu_n,actual and P_eq=cD[mu_B−(1−f)mu_n,actual].
  The unchanged identity Pactual−P_eq=c eta·sigma closes
  (Eeq+U+Echem)dot=−Lesc,star−Lnu,full−Lgamma. The fraction refers to the
  actual chemical-energy budget; it is not a microscopic branching ratio.
- **Required code changes after approval:** add one parameterized implementation
  of `IDirectBnvEnergyPartition`, initially alongside controlled partition fixtures;
  add a distinct run card/assembly and event/fate identities through the existing
  `FrozenControlledBnvRunContext` constructor. Validate f and mu_n,actual and
  record units, source identity, f and finite-T weighting provenance. Use prompt
  escape for the unretained share; no X reservoir or new chemical ODE is needed.
  Preserve old P2 files. Inspect `BnvDirectEnergyLedger::Evaluate` cancellation at
  tiny f; if needed, supply an auditable algebraically factored residual through
  that owner, still checked against actual potential and the escape budget.
- **Migration risk:** mislabelling f as product kinetic-energy retention; silently
  using old P2 fate/card or omission floor; cancellation of two ~868 MeV numbers
  when retaining ~1e−7 MeV; changing P_eq independently of the chosen escape energy.
- **Validation needed:** equation identities, limiting fixtures, tiny-f arithmetic,
  independent reference and new scope-specific numerical evidence below.
- **Existing outputs:** f=1 recovers the old model; Stage B remains valid in its
  original scope. New f values require new qualification and owner acceptance.

### B — Fixed equilibrium-budget fraction

- **Statement:** Eesc=(1−f)mu_B; P_eq=cD f mu_B and
  Pactual=cD(f mu_B+t_l·eta). This keeps the explicit storage correction relative
  to a state-independent retained baseline. It is distinct from A at nonzero η.
- **Required code changes:** same partition interface/run-card work, different
  energy formula; an explicit admissibility or energy-extraction convention.
- **Migration risk:** for filling chemistry, f<−t_l·eta/mu_B makes Pactual negative.
  The f=0 limit cannot be called a nonnegative cold deposition model by default.
- **Validation needed:** all ledger identities and boundary tests, including
  negative-residual refusal or a separately approved extraction interpretation.
- **Existing outputs:** f=1 still recovers P2; no approval to reinterpret P0.

### C — Microscopic BNV process now

- **Statement:** derive rate, ordinary stoichiometry, occupied-state weighting,
  final-state spectra and product fate together. Compute source counts and energy
  transfer from the same event kernel; q_dep and D are then linked.
- **Required code changes:** process-specific source/partition providers, applicable
  blocking/transport/finite-T terms, and X reservoir/structural feedback if products
  accumulate. These exceed a one-parameter extension.
- **Migration risk:** unconstrained microscopic assumptions or mixing different
  event measures; double-counting hole energy; applying cold/frozen validity beyond
  its bounds. ADR-0015:93 and :295 define the separation of owners.
- **Validation needed:** cited channel derivation, thermodynamic and local/global
  energy checks, rates/kinematics, product escape conditions, independent validation.
- **Existing outputs:** Stage B remains a mathematical regression, not channel evidence.

### D — Keep P2 unchanged; analytic study only

- **Statement:** use the governing equations and existing review diagnostics to
  locate conditional competition scales without adding a partition.
- **Required code changes:** none.
- **Migration risk:** treating fixed-T scaling or unreachable stationary estimates
  as a coupled physical prediction.
- **Validation needed:** analytic identities, declared regimes and dimensional checks.
- **Existing outputs:** unchanged; no low-retention trajectory becomes available.

Multiplying only the direct term in the thermal RHS, adding “chemical work” on
top of Pactual+LH, or dropping the escape-energy ledger are inadmissible shortcuts,
not conservative alternatives (ADR-0015:233).

## Consequences if an option is accepted

The first study would remain a phenomenological cold, spin-off, frozen-background
experiment with variable energy partition. No microscopic BNV rate, invisible-
channel zero-heating claim, astrophysical constraint, or general event qualification
would follow. Document the chosen semantics in the architecture and run-card
provenance when implemented. A and B must not share an ambiguous parameter name.

For A, f=0 is the existing cold P0 cancellation control. It is not a claim that a
physical invisible decay deposits no heat. Finite-T source/state terms must be
owned or shown negligible before interpreting a small residual physically.
The proposed stationary competition scale (~0.065 eV per baryon at the reviewed
T_end) is below illustrative smooth-weight omission scales; changing f does not
justify scaling those omissions by f (ADR-0015:60;
`docs/validation/results-first-bnv-pilot/stage-b/review-evidence/outputs/longterm.json:1`).

## Validation required for an approved implementation

Before implementation, the owner must select the energy convention and scope;
this draft supplies the derivation and reviewed baseline, not permission to run.
After implementation and before interpreting new trajectories:

1. Reproduce f=1 against the reviewed Stage B model, f=0 against the cold P0
   residual, and exact-zero source. Preserve species driving independently of f.
2. Check per-event fate closure, Pactual−P_eq=c eta·sigma, Edot_chem=A−LH and
   total energy including escape. Use finite-difference chemical storage checks
   and independent high-precision residual oracles at small f; negative tests
   must detect a missing escape term or duplicated chemical-work term.
3. Report gross LH, extra neutrinos, source transfer A, net storage and direct
   deposition independently. R20 cannot validate tiny chemical channels by itself;
   require channel-level accuracy and correctly scaled non-cancelling budgets.
4. Use linear/nonlinear stationary predictors only in their regimes; predeclare
   a small matrix and horizon after checking relaxation, thermal change, depletion
   and frozen-population drift. Recompute boxes/majorants and qualify any downward
   or other newly encountered events. Do not reuse the 2,048-year bounds for a
   ~50,000-year study (`docs/validation/results-first-bnv-pilot/stage-b/review-evidence/REVIEW.md:417`).
5. Establish a new numerical baseline and independent scientific review for the
   approved scope. No revised physical sign claim below an unowned finite-T floor.

Execution may use EKU parallel resources according to measured work and scheduler
availability; this ADR proposes no arbitrary CPU or walltime ceiling.

## Provenance

Drafted by the AI assistant on the owner's 2026-11-12 request to close Stage B
review and propose the next physics study. The owner has **not selected an option,
ratified Stage B, accepted this ADR, or approved new BNV implementation** in this
record. Exact entry commit and evidence hashes are in
`docs/validation/results-first-bnv-pilot/stage-b/review-evidence/IMPORT.json:1`.
This submission changes documentation and copies retained review evidence only;
no predictor, new trajectory, numerical campaign or physics implementation ran.
