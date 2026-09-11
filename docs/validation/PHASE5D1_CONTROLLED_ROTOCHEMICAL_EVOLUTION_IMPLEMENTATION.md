# Phase-5D-1 response implementation and pre-result blocker

**Disposition B: PHASE-5D RESPONSE MACHINERY VALIDATED — COUPLED SECULAR EVOLUTION BLOCKED ON NUMERICAL ARCHITECTURE.** Here the numerical blocker is the frozen upstream coefficient-acceptance budget, before ODE construction. This is not a finding that RKF45 is unsuitable, or that accepted ADR-0014 is inconsistent. There is no complete controlled-evolution candidate.

## Identity and authority

Canonical entry `d019ae390be4f5e3daba05039903485cb497e397`; branch `physics/phase5d-controlled-rotochemical-evolution`; worktree `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase5d-controlled-evolution`. Entry was clean, branch/worktree absent, canonical master and live origin master matched the entry. Canonical master was not modified. Source/architecture/validation specialists ran concurrently, read-only, using Astra HIGH. Lead/root Agent D was the sole writer. Required governance, ADR-0010/0011/0013/0014, integration/preflight/ratification records and current status documents were read. All 22 authenticated shared literature manifest entries passed; no literature was downloaded.

`PHASE5D1_PREDECLARATION_SHA = 697ec14610d63640866f3c0aed0a51d567cb037d`.
`PHASE5D1_RESPONSE_SHA = 14b74cba1a0c7b6faec591aedf6f5c9a76fa42d6`.
`PHASE5D1_EVOLUTION_SHA = NOT CREATED`.
`PHASE5D1_CANDIDATE_SHA = NOT CREATED`.

The first two exact required commit messages were used. The successful-coupling/candidate commit messages would overstate this outcome and are not used. The scratch implementation map preceded production writes and remains at `build/phase5d-audit/implementation-map.md`.

## Exact blocker and responsibility

The lead selected radial resolution 10000 in the predeclaration before qualifying its W precision. This was a preparation error in this implementation attempt, not a regression of governed Phase-5C. The fresh Structure-1 assembly at rho_c=1.10e15 g/cm^3 computationally accepted G and Z, then unchanged production `RotochemicalSpinDrive::Compute` refused W with `AccuracyGoalUnmet` at `ChemicalResponse.cpp:804`. The unchanged upstream `chemical_production_fixture` reproduces raw rc=1 without the Phase-5D coupling wrapper. The separate structural-validation-envelope gate at line 805 was not the failing gate.

| Channel | Positive structural terms in attempted W numerical budget [MeV s^2] | Frozen numerical goal [MeV s^2] | Ratio |
|---|---:|---:|---:|
| Npe | 8.918515424724095e-12 | 1e-12 | 8.918515424724095 |
| NpMu | 3.263850588763438e-11 | 4e-12 | 8.159626471908595 |

These are `sum_j abs(Z_ij)*E_Ij`, already exceeding each goal before adding nonnegative Z-error/cross-error/arithmetic terms. They diagnose the attempted discrete calculation, not certified physical continuum uncertainty. No W result escaped its factory and no coupled trajectory was generated. The raw `STOP` and rc=1 are retained, irrespective of the unrelated logger footer saying zero logged errors.

An additional preparation defect: the scratch certificate reused background/provider/anchor errors characterized at radial 40000/80000, recomputing only the new tail for radial10000. Computational G/Z acceptance therefore does not scientifically qualify radial10000. Honest fresh background characterization is still required; reducing/relabeling error or weakening W goals cannot repair this defect. Read-only specialist comparison found unchanged coefficient arithmetic and correct certificate transport/order. A different baseline resolution would amend the frozen predeclaration; its 20000 variant is explicitly a Ltilde comparison, not an authorized replacement baseline. The predeclaration's pre-result-impossibility stop rule is obeyed. No ADR revision is proposed and no impossibility at every resolution is claimed.

Machine evidence: [phase5d1_w_assembly_refusal.json](phase5d1_w_assembly_refusal.json) and [raw refusal excerpt](phase5d1_refusal_excerpt.txt). Full raw logs, certificate inputs, generated profile/table, toolchain/build logs and diagnostic code remain under `build/phase5d-audit/`. The debugger attempt was terminated without a usable result; the recorded numerical bound comes from explicit test-side diagnostics, not the debugger. No acceptance goal was relaxed for diagnosis.

## Retained production components

All retained physics is under `CompactStar/Physics/Rotochemical/`:

- `ChemicalImbalanceState.hpp`: typed Npe/NpMu ordering, eta_inf in MeV, explicit ChemState read/store, local-redshift and xi helpers. It is not yet an integrated evolved chemical block.
- `UrcaImbalanceFunctions.hpp`: independent pure FD/HD/FM/HM and stable F-minus-one polynomials. HM ends in pi^8, consistent with the ratified primary-source reconciliation; no published-erratum claim.
- `UrcaChannelNormalization.hpp`: typed process selection, local positive normalization, explicit ordered disjoint support intervals, independently constructible DU triangle check.
- `GlobalUrcaChannelCoefficient.hpp` and `src/GlobalUrcaChannelCoefficient.cpp`: immutable independent stellar integral, shared governed chemical-domain owner and currentness. Integral is `1e15 integral_Da 4*pi*r_km^2 exp(lambda+(2-q)*nu) S dr_km`; units erg s^-1 K^-q. It does not reuse G_y's inverse lapse or 1e54 count-volume normalization.
- `RotochemicalReactionResponse.hpp`: xi=eta_inf/(kB_MeV*Tinf_K); rate=Ltilde*Tinf^(q-1)*H/kB_erg in count/s; equilibrium/full/incremental luminosities and internal eta dot R in MeV/s. `RotochemicalThermalPower::From` converts that power once to erg/s and distinguishes LH-DeltaLnu from LH-Lnu_full.

New CMake file registers the global integrator. Existing files modified for retained code: only `CompactStar/Physics/CMakeLists.txt` and `tests/CMakeLists.txt`. Tests are `tests/rotochemical/response.cpp` and `oracles.py`. All thermal, photon, heat-capacity, EvolutionSystem, GSLIntegrator, StarContext, Z/W and EOS production files remain unchanged from entry.

Single kB authority: Zaki::Physics::K_BOLTZ_EV, scaled to MeV/K; kB_erg is derived using the governed MeV-to-erg conversion. Thermal input is K. Intended future state remains `(ln(Tinf/1e8 K),eta_npe_inf,eta_npmu_inf)`, eta in MeV. No deltaN public evolved state, torque law, A18, superfluidity, DU benchmark activation or BNV was added.

## Frozen experiment and unfinished coupling

Predeclared local SMe=1e-51 and SMmu=2e-51 erg cm^-3 s^-1 K^-8, positive constants on their species supports, are mathematical/architecture inputs only. Me/Mmu enabled; De/Dmu explicitly disabled. T0=1e8 K, eta0=(0,0) MeV. External prescribed timing: B=1e8 G, P0=1ms, PPdot=(B/3.2e19)^2, P=sqrt(P0^2+2 PPdot t), Omega=2pi/P, OmegaDot=-2pi PPdot/P^3. Interval 0..1e10 yr, year=365.25*86400s. Intended RKF45 rtol=1e-7, absolute tolerances=(1e-12,1e-18,1e-18), refined by 100. These are declarations, not completed-run settings/results.

No final T/eta/xi, heating activation, thermal sign history, transient, stiffness diagnostic, ODE convergence, radial Ltilde convergence or quasi-steady/1/7 fit exists. No physical fixture Ltilde was produced because assembly stops before that operation. The exact same-Ltilde identity is tested at response/ledger level; actual equilibrium-driver replacement, no-double-counting in a coupled RHS, spin-only ODE, cross-Z reaction-only/Lyapunov/linear oracles and state-ordering independence remain pending.

Uncommitted coupling work was archived intact at `build/phase5d-audit/uncommitted-coupling/` (tracked patch, new files, SHA256 manifest) and removed from the active source tree before final testing. It includes draft spin seams, frozen Z/W bundle, thermal adapter, component tolerances and benchmark runner. It must not be treated as implemented/validated. The read-only architecture audit identified six unfinished guards: callback currentness, shared spin owner, semantic Z channel validation, empty-channel thermal denominator, dependency checks before RHS accumulation, and controlled radial/provenance diagnostics. The immutable thermal-table/source-byte boundary was also unfinished. No evolution commit was made.

## Analytic evidence and limitations

Two focused CTests cover state/units/redshift, source functions, parity, small/large limits, four ratified roots, reaction sign and eta R positivity, zero-imbalance rates/increments/heating, response-level same-owner ledger, support/process/domain/currentness controls. Independent high-precision R1995 Fermi convolutions check both DU and MU F/H, relative tolerance 3e-14. RE10b independently compares a power-series antiderivative and mpmath quadrature before checking actual production GL integration at the predeclared relative tolerance 1e-10.

RE10b uses nonconstant nu=-0.4+0.1x^2 and lambda=0.2x^2, S=1e-40(2+x^2), disjoint support [0.2,0.45] union [0.7,0.9], q=6 and8. S is nonzero outside support. Seventeen distinct transformed-input cases (nine integral mutation families across q=6 and q=8; the q=6 wrong-q/extra-two-lapse alias is counted once) are rejected by the oracle: missing lapse, altered one-/two-lapse powers, wrong-sign lapse, omitted/inverted proper volume, wrong q and wrong domain. These are actual-integrator transformed-input controls, not source-code mutation builds. A final test-only accounting correction removes that redundant q=6 case; the physics, oracle, benchmark inputs and tolerances are unchanged. Other requested coupling mutation families remain pending; no inflated mutation total or complete RE-ladder pass is claimed. Test banners are not independent evidence. The final banner explicitly limits the claim to focused response checks; RE12/RE16-18 references on individual assertions cover only their response/support subsets.

## Retained-caveat closure table

| Caveat | Disposition at this stop |
|---|---|
| 1 enabled-process owner | Implemented typed UrcaProcessSelection. |
| 2 constructible DU sliver/configuration control | Passed synthetic positive triangle below guard; disabled DU source must never execute. |
| 3 predeclared RE10b tolerance | 1e-10 committed before coupled output. |
| 4 nontrivial nu/lambda | Passed curved manufactured fixture. |
| 5 S nonzero outside D_a | Passed; wrong-domain control nonvacuous. |
| 6 G_y/domain semantic relationship | Shared current owner and identity, separate lapse/volume normalization; physical fixture qualification pending. |
| 7 upstream INV-11 authority | Preserved; global INV-11 unresolved. |
| 8 banners not independent evidence | Explicitly retained; actual assertions/oracles/results determine coverage. |
| 9 M10/M11 historical qualification | Historical record unchanged; one algebraic integrated-lapse family, not two credits. |
| 10 innermost-first assertion | Unordered partition explicitly rejected. |
| 11 disconnected/outer support | Disjoint integrals and closed-inner throwing-source test pass. |
| 12 benchmark not realistic Ltilde | Explicit mathematical classification; realistic closure remains blocked. |

These closures concern retained response components; they do not close the full Phase-5D candidate.

## Validation and protected state

Final results: focused analytic 2/2, data-free 47/47 (disjoint 25+22 serial selections), complete authenticated suite 70/70; every completed suite raw rc=0, no CTest failures/skips. The final focused oracle confirms 17 distinct mutation cases. Inspectable commands, counts, log/source hashes and protected-file hashes are in [phase5d1_response_validation.json](phase5d1_response_validation.json), the companion final report and `build/phase5d-audit/suite-return-codes.json`. The entry build and separately executed Phase-5B/5C governed regressions passed raw rc=0 (2/2). The W assembly diagnostics failed raw rc=1 as reported above; these are not hidden as skipped/passing CTests. No coupled CTest is registered or claimed complete.

All ten existing governed baselines, the Phase-5C candidate and all authenticated data/literature files were byte-compared against entry: 54 files checked, zero mismatches. Phase5B SHA256 remains `7588f0e9cd62f5b6be48bb725e4d0ba6e47b64f50d843117e70c17a78beeb5fa`; Phase5C remains `7027aa6179fe9111d76586d467d007c4223ec9d26ece0fee773327973c77bfe7`; Phase5C candidate remains `a6b430b2356987b54a7c82c87d0b00f3e1d9698f6f299c32c80efcc40c63833b`. No new baseline, no regenerated baseline, no EOS/data/literature modification, no merge.

INV-11b/c/d have partial response implementation with focused analytic validation; integrated evolved-state/reaction/thermal clauses remain pending. INV-11e run-level frozen Z/W semantics remain pending. INV-11f remains unresolved and blocked before coupling. Global INV-11 remains UNRESOLVED / AWAITING IMPLEMENTATION, INDEPENDENT REVIEW AND OWNER RATIFICATION. Realistic A18/FR2005 normalization remains SOURCE-LIMITED / BLOCKED; BNV NOT BEGUN.

## One recommended next action

Review this pre-result blocker and authorize a revised predeclaration only after freshly qualifying the background and structural resolution against the unchanged governed Phase-5C acceptance goals, before resuming coupled evolution. This next action has not begun.
