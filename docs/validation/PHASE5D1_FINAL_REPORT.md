# Phase-5D-1 complete requested-field report

Disposition B. This report distinguishes declarations, retained response implementation, and work not reached. See the implementation record for the exact pre-result refusal.

| # | Requested field | Result |
|---:|---|---|
| 1 | Canonical entry SHA | d019ae390be4f5e3daba05039903485cb497e397. |
| 2 | Branch/worktree | physics/phase5d-controlled-rotochemical-evolution; /Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase5d-controlled-evolution. |
| 3 | Entry authentication | Clean new worktree from exact canonical/master/live origin SHA; no conflicting branch; canonical untouched. |
| 4 | Read-only specialists | Source, architecture, validation ran concurrently using Astra HIGH; no writes. |
| 5 | Sole writer | Root/lead Agent D. |
| 6 | PHASE5D1_PREDECLARATION_SHA | 697ec14610d63640866f3c0aed0a51d567cb037d. |
| 7 | Benchmark normalization values | SMe=1e-51; SMmu=2e-51 erg cm^-3 s^-1 K^-8 on declared species support. |
| 8 | Normalization provenance | Declared mathematical/architecture inputs only. |
| 9 | Initial T | 1e8 K, predeclared; no coupled run. |
| 10 | Initial eta | (0,0) MeV, predeclared. |
| 11 | Spin history | External B=1e8 G, P0=1ms; PPdot=(B/3.2e19)^2; P(t)=sqrt(P0^2+2 PPdot t); Omega=2pi/P; OmegaDot=-2pi PPdot/P^3. |
| 12 | Run interval | Predeclared 0..1e10 yr, year=365.25*86400s; not executed. |
| 13 | Production files | Rotochemical/ChemicalImbalanceState.hpp, UrcaImbalanceFunctions.hpp, UrcaChannelNormalization.hpp, GlobalUrcaChannelCoefficient.hpp, src/GlobalUrcaChannelCoefficient.cpp, RotochemicalReactionResponse.hpp, CMakeLists.txt. |
| 14 | Modified existing files | Physics/CMakeLists.txt, tests/CMakeLists.txt; narrow branch-status addenda in docs/SCIENTIFIC_INVARIANTS.md, docs/MODERNIZATION_ROADMAP.md, docs/architecture/CURRENT_ARCHITECTURE.md. |
| 15 | State vector | Intended (ln(Tinf/1e8 K), eta_npe_inf, eta_npmu_inf); integrated vector pending. |
| 16 | Eta ordering | Typed Npe then NpMu. |
| 17 | Eta units | MeV. |
| 18 | T state convention | Thermal API K; log state ln(Tinf/Tref), intended Tref=1e8 K. |
| 19 | Spin interface | Draft prescribed/state-coupled ISpinHistory archived; not in retained production. |
| 20 | Process-selection owner | UrcaProcessSelection. |
| 21 | Enabled processes | Me and Mmu for declared controlled configuration. |
| 22 | Disabled processes | De and Dmu explicitly. |
| 23 | DU-sliver control | Passed synthetic triangle-open positive densities independent of nB_min; disabled DU callable must never execute. |
| 24 | Disconnected support | Ordered vector of disjoint UrcaSupportInterval; closed interior throwing-source control passed. |
| 25 | F_D implemented | Yes, pure function only; no DU benchmark activation. |
| 26 | H_D implemented | Yes, pure function only. |
| 27 | F_M implemented | Yes. |
| 28 | H_M implemented | Yes. |
| 29 | Pi^8 confirmed | Yes, ratified HM final denominator; independent source-function tests. No published-erratum claim. |
| 30 | Exact-source tests | Independent high-precision R1995 Fermi convolutions for DU/MU plus ratified roots, tolerance 3e-14. |
| 31 | Parity tests | Passed F even/H odd including positive/negative xi. |
| 32 | Small-xi tests | Passed analytic slopes/leading increments. |
| 33 | Large-xi tests | Passed implemented asymptote checks; independent convolutions/roots additionally constrain coefficients. |
| 34 | Ltilde object/API | GlobalUrcaChannelCoefficient::Compute(UrcaIntegrationRequest). |
| 35 | Ltilde units | erg s^-1 K^-q. |
| 36 | Local normalization API | UrcaChannelNormalization: process, S(r), explicit support, provenance identity. |
| 37 | Global integration | 1e15 integral_Da 4pi r_km^2 exp(lambda+(2-q)nu) S(r) dr_km. |
| 38 | RE10b tolerance | Relative 1e-10, predeclared. |
| 39 | RE10b independence | Yes; independent series primitive and high-precision quadrature; expected value does not call production integrator. |
| 40 | RE10b mutants killed | 17 distinct transformed-input cases, 9 integral families at q=6/8; one wrong-q/extra-lapse alias counted once; not source-code mutation builds. |
| 41 | Wrong-domain nonvacuity | S=1e-40(2+x^2)>0 outside support; wrong whole-domain integral rejected. |
| 42 | Innermost-first assertion | Explicit increasing partition validation and reversed-order rejection. |
| 43 | Reaction-response object | RotochemicalReactionResponse. |
| 44 | Xi implementation | eta_inf_MeV/(BoltzmannMeVPerK*Tinf_K). |
| 45 | K_B authority | Zaki::Physics::K_BOLTZ_EV; MeV/K scaling and derived erg/K via governed conversion. |
| 46 | Reaction-rate units | count/s; Ltilde*Tinf^(q-1)*H/kB_erg. |
| 47 | Reaction sign | H odd with xi H>=0; rates share eta sign for positive active normalization. |
| 48 | Eta R positivity | Positive/negative eta tests passed. |
| 49 | Equilibrium Lnu | Ltilde*Tinf^q, channel-separated response result. |
| 50 | Nonequilibrium Lnu | Ltilde*F(xi)*Tinf^q, channel-separated response result. |
| 51 | Delta Lnu | Ltilde*(F-1)*Tinf^q; stable pure increment polynomial. |
| 52 | Chemical heating | Internal sum eta_inf*R [MeV/s], nonnegative in active tests. |
| 53 | MeV-to-erg boundary | RotochemicalThermalPower::From, once using governed conversion. |
| 54 | Same-Ltilde RE9 | Response/ledger-level identity passed; actual coupled cooling replacement proof pending. |
| 55 | Placeholder double counting prevented | Controlled replacement was draft only; no controlled coupled run. Legacy thermal code unchanged. |
| 56 | Chemical RHS | Ratified -Z R+2 W Omega OmegaDot; draft archived, production coupling pending. |
| 57 | Z source | Attempted fresh ChemicalImbalanceResponse from GlobalChemicalNumberResponse; never baseline/runtime copied values. |
| 58 | W source | Unchanged RotochemicalSpinDrive::Compute refused fresh radial10000 numerical budget. No W result consumed. |
| 59 | Spin-only oracle | Pending; no actual coupled/evolution oracle completed. |
| 60 | Reaction-only oracle | Pending actual cross-Z ODE; response signs/rates tested. |
| 61 | Lyapunov oracle | Pending active/dead cross-Z ODE coverage. |
| 62 | Cross-Z tested | Existing governed coefficient tests retained; new secular cross-Z tests pending. |
| 63 | Thermal RHS | Intended (all net power)/(Tinf*Cstar); archived draft only. |
| 64 | Photon cooling changed | NO. |
| 65 | Heat capacity changed | NO; scratch fixed-background entropy adapter created but not evolved. |
| 66 | ODE solver | Existing RKF45 retained unchanged; draft component-tolerance adapter archived. |
| 67 | Relative tolerance | Intended 1e-7; refinement 1e-9; not executed. |
| 68 | Per-state absolute tolerances | Intended (1e-12,1e-18,1e-18), refinement /100; not executed. |
| 69 | Stiffness diagnostic | Not reached; no evidence RKF45 unsuitable. |
| 70 | ODE convergence | Not reached. |
| 71 | Ltilde convergence | Manufactured oracle passed; physical-fixture quadrature/radial refinement not reached. |
| 72 | Coupled benchmark completed | NO; stopped before any trajectory. |
| 73 | Final Tinf | Not generated. |
| 74 | Final eta_npe | Not generated. |
| 75 | Final eta_npmu | Not generated. |
| 76 | Final xi_e | Not generated. |
| 77 | Final xi_mu | Not generated. |
| 78 | Heating activated | Not assessed in a trajectory; positive response heating tested. |
| 79 | Enhanced neutrinos activated | Not assessed in a trajectory; positive response increments tested. |
| 80 | Incremental beta thermal sign | Four source roots verified; time-history sign behavior not generated. |
| 81 | Quasi-steady/tracking | Not reached. |
| 82 | Expected 1/7 scaling | Frozen fit criterion retained, no fit/run. |
| 83 | Transient behavior | Not generated. |
| 84 | Mutation count | 17 distinct transformed-input integrator cases / 9 families; other coupling families pending. |
| 85 | Mutation failures detected | All 17 distinct cases were distinguished beyond 1e-10 by independent oracle; support/domain/currentness negative controls also passed, not added as inflated mutation count. |
| 86 | Retained-caveat closure | Twelve-row table in implementation record; response closures do not close full candidate. |
| 87 | PHASE5D1_RESPONSE_SHA | 14b74cba1a0c7b6faec591aedf6f5c9a76fa42d6. |
| 88 | PHASE5D1_EVOLUTION_SHA | NOT CREATED. |
| 89 | Implementation record | docs/validation/PHASE5D1_CONTROLLED_ROTOCHEMICAL_EVOLUTION_IMPLEMENTATION.md. |
| 90 | Focused analytic suite | 2/2 PASS, raw rc=0; final focused run 3.06 seconds; 17 distinct mutation cases. |
| 91 | Focused coupled suite | No completed/registered coupled suite; not passed or counted as skipped. |
| 92 | Phase-5B regression | Entry separately passed raw rc=0, also included in final full verification. |
| 93 | Phase-5C regression | Entry separately passed raw rc=0, also included in final full verification. |
| 94 | Data-free suite | 47/47 PASS in disjoint serial selections: 25 self-contained + 22 complementary non-external tests; each raw rc=0. |
| 95 | Full authenticated suite | 70/70 PASS, raw rc=0, 1730.8 seconds. |
| 96 | Raw return codes | Final build and all completed validation suites raw rc=0. W assembly refusal raw rc=1 retained explicitly; it is not a passing CTest. |
| 97 | Failures/skips | Zero CTest failures/skips. Four diagnostic W-assembly attempts returned rc=1 for the same explicit refusal; the unused debugger attempt was terminated (rc=143). No coupled run exists; incomplete work is not hidden as a test skip. |
| 98 | Ten baselines unchanged | Yes, byte hashes match entry; no new baseline. |
| 99 | Phase-5B hash | 7588f0e9cd62f5b6be48bb725e4d0ba6e47b64f50d843117e70c17a78beeb5fa unchanged. |
| 100 | Phase-5C hash | 7027aa6179fe9111d76586d467d007c4223ec9d26ece0fee773327973c77bfe7 unchanged; candidate a6b430b2356987b54a7c82c87d0b00f3e1d9698f6f299c32c80efcc40c63833b unchanged. |
| 101 | EOS/data diff | None; scratch generated fixture outside governed data. |
| 102 | Literature diff | None; no downloads, all 22 manifest entries authenticated. |
| 103 | INV-11b | Partial typed response-state implementation; integrated evolved ownership pending. |
| 104 | INV-11c | Response sign/order analytic tests pass; full secular coupling pending. |
| 105 | INV-11d | Response thermal-ledger implementation/test only; coupled no-double-counting pending. |
| 106 | INV-11e | Global coefficient currency implemented; frozen Z/W run semantics pending. |
| 107 | INV-11f | UNRESOLVED / IMPLEMENTATION + VALIDATION PENDING; blocked before coupling. |
| 108 | Global INV-11 | UNRESOLVED / awaiting completion, independent review and owner ratification. |
| 109 | Realistic A18 implemented | NO. |
| 110 | Realistic FR2005 normalization claimed | NO. |
| 111 | BNV begun | NO. |
| 112 | PHASE5D1_CANDIDATE_SHA | NOT CREATED; terminal evidence commit is a blocker record, not a complete candidate. |
| 113 | Branch local/upstream/live | Published blocker-record tip (not a candidate SHA); exact local/upstream/live equality is recorded in the final delivery artifact and response. |
| 114 | Merged | NO; canonical master unchanged. |
| 115 | Final disposition | B — RESPONSE MACHINERY VALIDATED; COUPLED SECULAR EVOLUTION BLOCKED ON NUMERICAL ARCHITECTURE (upstream frozen W-accuracy gate, not a demonstrated ODE limitation). |
| 116 | Exact remaining caveats | Radial10000 W budget fails; background uncertainty not freshly qualified; archived coupling needs six architecture guards and full RE9–14/16–18 integration, mutation, convergence and trajectory evidence. No complete candidate, no source-realism claim, no adjudicated ADR contradiction. |
| 117 | Exactly one next action | Review the pre-result blocker and authorize a revised predeclaration only after freshly qualifying background/structural resolution against unchanged governed Phase-5C goals, before resuming coupled evolution. Not begun automatically. |
