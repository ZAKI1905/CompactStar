# ADR-0025 — Controlled BNV–rotochemical competition study (retained-energy partition as a secondary axis)

| Field | Value |
|---|---|
| Status | **PROPOSED — owner selection and approval pending** |
| Date | 2026-11-12 (revised the same day from the draft at `0e88ec7e5ee981fcc2d93c85cf8a17d53d6d15fd`) |
| Change class | Scientific-semantic if implemented; this submission is documentation only |
| Governing authority | GOVERNANCE §§1–5; ADR-0014 §§3.1–3.12, 3.16; ADR-0015 §§3–14; ADR-0016 tangent ownership; INV-15 |
| Affected invariants | INV-11 controlled state/reaction/ledger semantics, extended to simultaneous spin and BNV driving; INV-15 common Cstar; no broader invariant closure |
| Blocks | Implementation and execution of any combined spin + BNV evolution, and of a variable-retention BNV study |

The file name is retained for link stability. The earlier draft proposed a
retained-energy parameter with spin off. This revision makes the competition
between BNV-induced and spin-down-induced chemical disequilibrium the object of
the decision.

## Context

Stage B at `49253200a801ad9b336789dffb4e93d535a8e903` has independent disposition
B — CONDITIONAL PASS, with owner ratification pending. It verifies fixed full-
retention P2 inputs with spin off. Its record is unchanged by this ADR
(`docs/validation/results-first-bnv-pilot/stage-b/CONDITIONAL_ACCEPTANCE.md:1`).

Three accepted facts frame the decision.

1. **The combined chemical equation is already contract.** ADR-0015 A3 reads
   `eta_dot = −Z(R + sigma) + 2 W Omega Omegadot`, and A5 gives the reservoir
   derivative with both terms (`docs/adr/ADR-0015-bnv-open-system-thermal-ledger.md:144`,
   `:219`). ADR-0014 §3.4 fixes the spin term's sign
   (`docs/adr/ADR-0014-secular-rotochemical-evolution-contract.md:149`).
2. **The code computes the sum but refuses to run it.** The ordinary evaluator
   carries `2 W Omega Omegadot`
   (`CompactStar/Physics/Rotochemical/FrozenRotochemicalRunContext.hpp:141`) and
   the BNV context subtracts `Z sigma`
   (`CompactStar/Physics/BNV/src/FrozenControlledBnvRunContext.cpp:139`). Spin
   with a nonzero BNV source is refused at construction and at every evaluation
   (`:35`, `:63`, `:114`).
3. **Global closure with imposed spin is explicitly not claimed.** ADR-0015 §9
   requires real `Omegadot J` work and state derivatives
   (`docs/adr/ADR-0015-bnv-open-system-thermal-ledger.md:255`).

The [scientific proposal](../validation/results-first-bnv-pilot/stage-b/NEXT_PHYSICS_STUDY.md)
derives the combined system, the ledger, the dimensionless measures and the
experiment. Its findings that bear on this decision:

- A neutron sink and spin-down drive the imbalances with opposite signs. On the
  governed Structure-1 coefficients the drive vectors are antiparallel to within
  16%. Heating depends on the net drive and is not additive.
- The controlling ratio is `beta_l = D t_l / (2|I_Omega,l| Omega|Omegadot|)`.
  On Structure-1 it is unity near `Gamma = 1e-12 yr^-1` for a 5 ms, 1e-20
  pulsar.
- The partition convention matters at order unity once spin makes `eta > 0`
  (see P-A and P-B below).
- Zero direct deposition is a diagnostic limit. For events whose products all
  leave the fluid, the deposit is the Fermi-hole energy, and competition needs
  it below about a keV per destroyed baryon.
- Reaching a BNV-driven quasi-steady state in time needs depletion near 1e-4,
  against the 1e-6 gate of the existing frozen certificate
  (`CompactStar/Physics/BNV/src/FrozenBnvValidity.cpp:66`). The Structure-1
  quasi-steady state for realistic drives may also lie below the 1e6 K
  certificate floor (`CompactStar/Physics/BNV/src/FrozenControlledBnvRunContext.cpp:108`).

## Decision

## Alternatives

The owner's choice has one primary part (study design) and four dependent parts
(partition, spin matching and history, gates, validation tier). The AI
assistant's recommendations are marked and are not selections.

### Study design

#### A — Staged competition study (recommended for owner consideration)

- **Statement:** study how a neutron-sink BNV drive changes a rotochemically
  heated star, in stages S1–S3 of the proposal:
  S1, a static quasi-steady-state map over spin drive and baryon-loss rate at
  zero direct deposition, whose corners are the four cases
  {neither, spin only, BNV only, both};
  S2, the minimal combined-drive enablement and four short matched coupled runs;
  S3, time-dependent approach to the quasi-steady state, with sweeps in `Gamma`,
  spin-down rate and deposit per event.
- **Governing equations:** proposal eqs. (2)–(5). The thermal equation is
  ADR-0015 A6 unchanged; spin adds no thermal term. The new conserved statement
  is the thermal + chemical subsystem ledger
  `d/dt[U + c Echem] = P_actual + A_B + A_Omega − DeltaLnu − Lnu,eq − Lgamma`,
  with `A_Omega = 2c Omega Omegadot eta·I_Omega`. It is a **derived, unratified
  candidate**; accepting this option accepts it for the controlled scope only.
- **Required code changes after approval:**
  S1 none in production source — a standalone tool reading governed coefficient
  artifacts.
  S2: (i) a combined run purpose admitting spin with nonzero BNV source in
  `FrozenControlledBnvRunContext`; (ii) a parameterized prescribed spin history,
  because `PrescribedDipoleHistory` is hard-coded and required by type
  (`CompactStar/Physics/Rotochemical/PrescribedSpinHistory.hpp:26`,
  `CompactStar/Physics/Rotochemical/FrozenRotochemicalRunContext.hpp:76`, `:105`)
  and native Structure-1 inputs are restricted to analytic controls (`:90`);
  (iii) the spin term in the reservoir diagnostic, absent from
  `Echem_dot_total` (`CompactStar/Physics/BNV/src/FrozenControlledBnvRunContext.cpp:164`),
  and the subsystem-ledger residual; (iv) a run-card zero-deposition partition,
  which exists only as a test fixture (`tests/bnv/energy_partition_fixtures.hpp:28`).
  S3: the parameterized partition of P-A or P-B, and whatever the gate choice
  below requires. No new ODE and no new chemical-work term.
- **Migration risk:** treating the two mechanisms as additive heaters; entering
  the spin drive twice (a BNV torque added on top of an observed `Omegadot`);
  adding `A_Omega`, `A_B` or `Omegadot J` to the thermal right-hand side;
  reading Structure-1 temperatures as predictions; relaxing existing gates
  without recording it; reusing Stage B run cards or evidence.
- **Validation needed:** the proportional validation below.
- **Existing outputs:** unchanged. Spin off recovers Stage B and the P0/P2
  fixtures; zero source recovers the governed Phase-5D baseline.

#### B — Retained-energy study with spin off (the earlier draft)

- **Statement:** vary the retained fraction at the existing interface with spin
  off, as in the draft at `0e88ec7`.
- **Required code changes:** one parameterized partition and run card.
- **Migration risk:** answers a narrower question. Without spin there is no
  rotochemical baseline to modify, and the competition scale found there
  (0.065 eV per baryon at 1e8 K) lies below the finite-temperature floor.
- **Validation needed:** as in the earlier draft: limiting fixtures at `f = 0, 1`,
  tiny-`f` arithmetic, ledger identities.
- **Existing outputs:** unchanged.

#### C — Microscopic BNV process first

- **Statement:** derive rate, stoichiometry, occupied-state weighting, in-medium
  kinematics, spectra and product fate for one process before any combined run.
- **Required code changes:** process-specific source and partition providers,
  finite-temperature terms, an X reservoir if products are retained.
- **Migration risk:** large scope before the coupling is characterized; the
  proposal's classification suggests most candidate processes fall outside the
  competition regime, which a cheap map would show first.
- **Validation needed:** cited channel derivation, local and global energy
  checks, escape conditions, independent review (ADR-0015 §§5, 11–13).
- **Existing outputs:** unchanged.

#### D — Analytic study only

- **Statement:** keep the proposal's equations and closed-form estimates; run
  nothing.
- **Required code changes:** none.
- **Migration risk:** the estimates assume the large-imbalance regime and a
  reached quasi-steady state; neither is checked.
- **Validation needed:** algebraic identities and dimensional checks.
- **Existing outputs:** unchanged.

### Energy-partition convention (needed for A and B)

Both are conservative. They coincide at `eta = 0` and at `f = 1`. With spin on
they differ by the full size of the effect under study, so the choice cannot be
deferred.

#### P-A — actual-potential fraction (recommended)

- **Statement:** `E_esc,fluid = E_esc,star = (1 − f) mu_n,actual`,
  `E_X = E_in = 0`. `P_actual = c D f mu_n,actual`, zero at `f = 0` for any
  `eta`. The escaping product carries the actual Fermi energy, as a Fermi-
  surface event with no further deposition would. Report the deposit as
  `q_dep = f mu_n,actual` in keV per event.
- **Required code changes:** one parameterized `IDirectBnvEnergyPartition`
  (`CompactStar/Physics/BNV/DirectEnergyLedger.hpp:33`) with its own run card,
  fate and finite-temperature weighting provenance; a factored residual in the
  ledger owner if cancellation at tiny `f` requires it
  (`CompactStar/Physics/BNV/src/DirectEnergyLedger.cpp:73`).
- **Migration risk:** reading `f` as a branching ratio or a kinetic-energy
  fraction; cancellation of two ~868 MeV numbers at `f ~ 1e-6`; reusing the
  Stage B omission floor.
- **Validation needed:** `f = 1` against Stage B, `f = 0` against P0,
  independent high-precision residual at small `f`.
- **Existing outputs:** unchanged; `f = 1` recovers P2.

#### P-B — equilibrium-budget fraction

- **Statement:** `E_esc = (1 − f) mu_B`.
  `P_actual = c D (f mu_B + t_l·eta)`. At `f = 0` this is
  `c D t_l·eta = −A_B`: positive when spin-down has made `eta > 0`, negative in
  a BNV-built state. The energy the sink drains from the chemical reservoir is
  deposited as heat instead of escaping.
- **Required code changes:** as P-A with a different energy formula, plus a
  negative-residual refusal or an approved energy-extraction convention.
- **Migration risk:** `f = 0` is not a zero-deposition control once `eta ≠ 0`;
  a negative direct residual in BNV-built states.
- **Validation needed:** all ledger identities and boundary tests, including the
  refusal.
- **Existing outputs:** unchanged; `f = 1` recovers P2; no reinterpretation of
  P0.

Neither may multiply only the direct term in the thermal right-hand side while
leaving the escape ledger unchanged
(`docs/adr/ADR-0015-bnv-open-system-thermal-ledger.md:233`).

### Spin matching and spin history (needed for A)

S1 needs neither choice: the map takes `Omega Omegadot` as a parameter.

| Option | Statement | Required code changes | Migration risk | Validation needed | Existing outputs |
|---|---|---|---|---|---|
| **S-obs** (recommended) | Spin-on cases share one prescribed `Omega(t)`: two stars with the same timing solution. Any BNV-induced spin change is already inside the prescribed `Omegadot`. | None beyond the history below. | Adding a separate BNV torque term on top. | Identical spin identity recorded in both run cards. | Unchanged. |
| S-torque | Spin-on cases share one external torque; `Omegadot` then differs by `kappa Gamma Omega`, `kappa = O(1)`. | `dI/dB`, product angular momentum, a state-coupled spin adapter (ADR-0014 §3.12). | On Structure-1 the effect is about `0.39 kappa (1 ms/P)^2` of the BNV drive, the order of terms the frozen model drops. | Angular-momentum ledger; second-order rotation estimate. | Unchanged. |
| **H-new** (recommended for S2) | Parameterized sub-Kepler prescribed history. | One class under `PrescribedSpinHistory`; admission in `CompactStar/Physics/Rotochemical/FrozenRotochemicalRunContext.hpp:76`, `:105`. | Weakening the type guard for other callers. | Spin-only oracle; sub-Kepler and `|I_Omega,l|Omega^2/N_l` checks. | Unchanged; governed history kept. |
| H-governed | Reuse the governed `P0 = 1 ms` dipole history. | None. | Super-Kepler on this star for about 2 Gyr; mathematical driver only (ADR-0014 §15). | Existing Phase-5D regression. | Unchanged. |

### Gates on depletion and temperature (needed for S3 of A)

No recommendation between G-extend and G-declared is made before S1 shows which
part of the plane matters.

| Option | Statement | Required code changes | Migration risk | Validation needed | Existing outputs |
|---|---|---|---|---|---|
| G-short | Keep the existing certificate: depletion ≤ 1e-6, 1e6–1e9 K. Sufficient for S2's slope, power and ledger checks. It yields no resolvable temperature comparison and cannot follow the approach to a BNV-driven quasi-steady state. | None. | Reading S2 as a thermal result. | Existing monitor. | Unchanged. |
| G-extend | Extend the frozen certificate to depletion of a few 1e-4 and state the thermal range needed. A bounded static task. | New certificate samples from the existing producer (`tests/bnv/produce_frozen_certificates.py:1`); the hard limit at `CompactStar/Physics/BNV/src/FrozenBnvValidity.cpp:66`. | Thresholds may not be met for the muon channel (`d ln N_mu/d ln B ≈ 54`). | Sensitivity samples and utilization at the new depletion. | Unchanged; old certificate kept. |
| G-declared | Exploratory tier with a declared drift estimate (about 1% at 2e-4 depletion) in place of a certificate. | An explicit, separately named run purpose that records the estimate. | Results mistaken for certified ones. | Drift estimate reported per run; status carried in provenance; no promotion. | Unchanged. |
| G-sliding | Evolving background and coefficients with the `+Zdot Z^-1 eta` term (ADR-0014 §3.11). | Time-dependent coefficient owners. | Large scope; mutation M18. | Variable-`Z` chain-rule oracle (ADR-0015 §16). | Unchanged. |

### Validation tier (needed for A)

| Option | Statement | Required code changes | Migration risk | Validation needed | Existing outputs |
|---|---|---|---|---|---|
| **V-proportional** (recommended) | Regression to existing governed corners, ledger closure, tolerance or step halving, independent re-derivation of the map. No high-precision enclosure campaign for capabilities Stage A/B and Phase-5D already validated. | Test and tool additions only. | Claiming certified status for exploratory results. | The list in "Validation required" below. | Unchanged. |
| V-certified | Stage-B-style enclosures and event qualification for every new trajectory. | New domain boxes, majorants and event coverage per scope. | Cost disproportionate to a first comparison. | As Stage B, per scope. | Unchanged. |

Multiplying only the direct term in the thermal equation, adding chemical or
spin work on top of `P_actual + L_H`, or dropping the escape-energy ledger are
inadmissible shortcuts, not alternatives (ADR-0015 §§9, 19).

## Consequences if an option is accepted

Under design A the first study remains a cold, frozen-background, non-superfluid,
modified-Urca experiment on a 0.62 M☉ free-gas star with benchmark reaction
normalization and an abstract neutron sink. Its results are ratios, exponents,
interaction contrasts and a regime map. No temperature for a real pulsar, no
BNV rate limit, no claim that an invisible decay deposits no heat, and no
general event qualification would follow.

Accepting A would extend the controlled INV-11 scope to simultaneous spin and
BNV driving and would accept the subsystem ledger for that scope. It would not
resolve global INV-11, the rotating global first law, superfluidity, realistic
normalization or A18.

The zero-deposition case is the existing cold cancellation control. The
proposal's classification finds no generic microscopic realization of it. A
physical statement needs a named process with its deposit per event shown to be
below the chemical work per event, and the finite-temperature terms owned
(ADR-0015 §4). Two sign-conditional statements of ADR-0015 (§5 upper bound,
§10 "nonpositive while storage fills") assume a neutron-sink-built state and do
not hold when spin-down makes `eta > 0`; any implementation must not test them
in spin-on cases. Whether ADR-0015 needs an explicit note is an owner decision.

Document the chosen partition, matching rule, gate route and tier in the run-
card provenance. P-A and P-B must not share a parameter name.

## Validation required for an approved implementation

Before implementation the owner selects the design and its dependent parts.
After implementation, and before interpreting new results:

1. **Corners.** Spin off, zero source: Stage B control. Spin on, zero source:
   governed Phase-5D baseline. Spin off, `f = 1`: Stage B source model. Spin
   off, `f = 0`: P0 fixture. The combined right-hand side must equal the sum of
   the two drive terms it already computes separately.
2. **Identities.** `P_actual − P_0 = c eta·sigma`;
   `c Echem_dot = A_B + A_Omega − L_H` against a finite-difference of `Echem`;
   the subsystem ledger over finite intervals; `b·t = 1`; per-event fate closure.
3. **Signs.** Neutron sink alone gives both `eta < 0`; spin-down alone gives both
   `eta > 0`; the electron null lies at `beta_e = 1` and the muon null at
   `beta_mu = 1`; the Cowling `k`-route negative oracle of ADR-0015 §7 still
   refuses.
4. **Negative tests.** A missing `A_Omega`, a duplicated chemical-work term in
   the thermal equation, a doubled spin drive and a missing escape term must each
   be detected.
5. **Map.** Independent re-derivation of the quasi-steady-state map; FR2005-form
   scaling at zero source; fitted exponents against 8/7 and 2/7; stability of
   every reported state; interaction contrast against the two reference curves.
6. **Gates.** Depletion, temperature range, `|I_Omega,l| Omega^2 / N_l` and
   relaxation-to-evolution ratios reported for every point and run. A result
   outside a certificate carries that status.
7. **Convergence and review.** Tolerance or step halving on every trajectory;
   an independent scientific review of the new scope. No physical sign claim
   below an unowned finite-temperature residual.

Execution may use EKU parallel resources according to measured work and
scheduler availability. This ADR proposes no CPU or walltime ceiling.

## Provenance

Drafted by the AI assistant on the owner's 2026-11-12 request to refocus the
next physics study on BNV–rotochemical competition, starting from
`0e88ec7e5ee981fcc2d93c85cf8a17d53d6d15fd` on
`candidate/results-first-bnv-pilot-20261108`. The owner has **not selected an
option, ratified Stage B, accepted this ADR, or approved new BNV or rotochemical
implementation** in this record.

This submission changes documentation only. No source, test, tool, run card,
certificate, baseline or evidence file changed. No predictor, map, trajectory or
numerical campaign ran. The fixture ratios quoted in the proposal are closed-form
arithmetic on existing governed values and are not retained evidence.

Published sources are tracked under `docs/references/` (owner-supplied
2026-11-12; hashes in `docs/references/SHA256SUMS`). Every equation, limit and
quoted statement cited in the proposal was checked against the PDF text at the
pages listed in its §11. `R2006` is not in the tracked set and is used only
through ADR-0014. The owner, a co-author of Berryman, Gardner & Zakeri 2023,
stated that its remarks on Urca cooling and net heating record an open question
and are not results; the proposal does not rely on them.
Stage B entry commit and evidence hashes remain in
`docs/validation/results-first-bnv-pilot/stage-b/review-evidence/IMPORT.json:1`.
