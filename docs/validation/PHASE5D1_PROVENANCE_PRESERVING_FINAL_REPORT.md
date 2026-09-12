# Phase-5D-1C2 final disposition and 130-field report

Candidate evidence only. The implementation-team specialist reviews are distinct from the required future independent scientific review. No next action below has begun. Post-commit remote identities are resolved in the final user response.

| # | Required field | Result |
|---:|---|---|
| 1 | Canonical master SHA | d019ae390be4f5e3daba05039903485cb497e397 |
| 2 | Branch / worktree | physics/phase5d-controlled-rotochemical-evolution / /Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase5d-controlled-evolution |
| 3 | Entry branch SHA | f7116c1408c06f976527f86d4397ad6d4540dedf |
| 4 | Sole writer | Root / integrator; exactly one repository writer |
| 5 | Read-only specialists | Provenance, thermal and numerical specialists, Astra HIGH; no repository writes; implementation-team review only |
| 6 | Protected governed path count | 33 |
| 7 | All protected entry hashes recorded? | YES, before implementation; all 33 included in candidate protected_entry |
| 8 | Archived draft audited? | YES; all 112 manifest hashes and tracked patch verified; no whole patch applied |
| 9 | Archived files reused | None restored verbatim; quadrature/fixture/exponential/validator ideas reaudited and adapted |
| 10 | Archived components reimplemented | 22; exact per-path mapping in candidate archived_draft_audit.decisions |
| 11 | Archived components rejected | 6: both governed Analysis source edits and four prior draft implementation/status/evidence files; exact paths in candidate mapping |
| 12 | Frozen context path / API | CompactStar/Physics/Rotochemical/FrozenRotochemicalRunContext.hpp; constructor consumes semantic authorities, RequireFullCurrent, RequireCheapCurrent, RequireOwners, Evaluate, typed Z/W/I/Ltilde access |
| 13 | Owns semantic Z? | YES, shared const ChemicalImbalanceResponse and its G authority |
| 14 | Owns semantic W? | YES, const RotochemicalSpinDrive computed once from exact retained Z/fixed authorities |
| 15 | Owns fixed-baryon authority? | YES, private const copy with retained metadata and G lifetime coverage |
| 16 | Owns Ltilde authorities? | YES, shared GlobalUrcaChannelCoefficient; all controlled terms use it |
| 17 | Owns thermal source bytes? | YES, four retained paths/bytes/SHA256 plus const parsed table |
| 18 | Owns shared spin owner? | YES, one shared_ptr<const PrescribedSpinHistory> |
| 19 | Lifetime strategy | Strong semantic/G/provider/revision/star/source/thermal/spin chain; private typed snapshots; caller-handle destruction passes |
| 20 | Currentness strategy | Full factory/source validation at construction and before runs; cheap revision/profile/structural/owner/spin checks before and after every RHS evaluation |
| 21 | Runtime baseline JSON? | NO |
| 22 | Analysis source buffering? | NO |
| 23 | Radial resolution | 80000 |
| 24 | EOS resolution | 8192 |
| 25 | rho_c | 1.10e15 g cm^-3 |
| 26 | SMe | 1e-51 erg cm^-3 s^-1 K^-8 |
| 27 | SMmu | 2e-51 erg cm^-3 s^-1 K^-8 |
| 28 | Ltilde_Me | 1.812909813491381e-32 erg s^-1 K^-8 |
| 29 | Ltilde_Mmu | 9.8504968582611409e-34 erg s^-1 K^-8 |
| 30 | Physical Ltilde convergence | PASS: GL refinement relative 4.8309713681580895e-14 / 6.9461243999055514e-15 <=1e-6; radial comparison 2.3628271923805636e-4 / 4.8732697187797193e-8 <=5e-3 |
| 31 | State vector | (ln(Tinf/1e8 K), eta_npe_inf [MeV], eta_npmu_inf [MeV]); typed semantic channels |
| 32 | Component-tolerance adapter | GSL scaled_new(1,rtol,1,0,component_atol,3): atol_i + rtol*abs(y_i); actual GSL error levels tested |
| 33 | Baseline rtol | 1e-7 |
| 34 | Baseline atols | (1e-12,1e-18,1e-18) |
| 35 | Refined rtol | 1e-9 |
| 36 | Refined atols | (1e-14,1e-20,1e-20) |
| 37 | Callback-currentness guard | PASS; executable positive/refusal controls, with governed regression gates separately passed |
| 38 | Shared-spin guard | PASS; executable positive/refusal controls, with governed regression gates separately passed |
| 39 | Semantic-Z guard | PASS; executable positive/refusal controls, with governed regression gates separately passed |
| 40 | Empty-channel guard | PASS; executable positive/refusal controls, with governed regression gates separately passed |
| 41 | Dependency-before-RHS guard | PASS; executable positive/refusal controls, with governed regression gates separately passed |
| 42 | Radial/provenance guard | PASS; executable positive/refusal controls, with governed regression gates separately passed |
| 43 | Thermal-source-byte guard | PASS; executable positive/refusal controls, with governed regression gates separately passed |
| 44 | Governed-source-immutability guard | PASS; executable positive/refusal controls, with governed regression gates separately passed |
| 45 | Controlled equilibrium cooling owner | Same context-owned GlobalUrcaChannelCoefficient / private FrozenReactionEvaluator |
| 46 | Historical placeholder handling | Disabled in this controlled context; actual placeholder injection detected; global legacy behavior unchanged |
| 47 | Same-Ltilde coupled RE9 | PASS, exact zero-imbalance reduction including accumulated derivative |
| 48 | Chemical RHS | dot eta = -Z R + 2 W Omega OmegaDot, both cross terms, no dot(Z) |
| 49 | Thermal RHS | dot ln(Tinf/1e8) = [LH-DeltaLnu-Lnu_eq_controlled-Lgamma-Lnu_other]/(Cstar*Tinf) |
| 50 | MeV-to-erg boundary | Exactly once through RotochemicalThermalPower::From |
| 51 | Spin-only oracle | PASS; both signs and nonzero initial eta |
| 52 | Spin-only maximum relative error | 4.1907178240110105e-16 |
| 53 | Reaction-only oracle | PASS; positive/negative eta, two-temperature small-xi matrix exponential |
| 54 | Reaction-only maximum relative error | 1.0690936789558023e-09 |
| 55 | Cross-Z result | PASS, both rows retained; separate row omissions detected |
| 56 | Active Lyapunov | PASS |
| 57 | Dead-channel Lyapunov | PASS, cross-Z motion retained |
| 58 | Coupled process/support controls | PASS: disabled De/Dmu despite open triangle; disconnected support does not activate closed interior; currentness failure precedes write |
| 59 | Unit mutations | All 7 named unit cases detected; double-kB interpretations share one conservative family |
| 60 | Sign / structural mutations | All 11 named transformed sign/index/cross-Z cases detected, including full-map transpose and true state-slot swap; sign aliases counted once; stale/foreign owners separately refused |
| 61 | Thermal-ledger mutations | All 4 detected: full instead of increment, omitted heating, double equilibrium, actual legacy placeholder; full/double-equilibrium alias counted once |
| 62 | Pretrajectory protected hashes | PASS, 33/33 equal entry |
| 63 | Pretrajectory Phase-5B regression | PASS, 1/1 |
| 64 | Pretrajectory Phase-5B raw rc | 0 |
| 65 | Pretrajectory Phase-5C regression | PASS, 1/1 |
| 66 | Pretrajectory Phase-5C raw rc | 0 |
| 67 | First trajectory completed? | YES |
| 68 | Start time | 0 s (0 yr) |
| 69 | End time | 3.15576e17 s (1e10 yr) |
| 70 | Final Tinf | 1030441.4596289885 K |
| 71 | Final Tsurf_inf | 75082.805533038438 K |
| 72 | Final eta_npe | 0.021510445901899329 MeV |
| 73 | Final eta_npmu | 0.024299132258749929 MeV |
| 74 | Final xi_e | 242.24409542012097  |
| 75 | Final xi_mu | 273.64943248317354  |
| 76 | Final Omega | 2347.5475293446275 rad/s |
| 77 | Final OmegaDot | -3.2002469209335231e-15 rad/s^2 |
| 78 | Final R_e | 1.752286823396162e+36 count/s |
| 79 | Final R_mu | 2.2321577268888389e+35 count/s |
| 80 | Final Lnu_eq_controlled | 24296377871460988 erg/s |
| 81 | Final DeltaLnu | 2.5954427679404978e+28 erg/s |
| 82 | Final LH | 6.9080121470151791e+28 erg/s |
| 83 | Final DeltaP_beta | 4.3125693790746813e+28 erg/s |
| 84 | Accepted steps | 8828 |
| 85 | Rejected steps | 1618 |
| 86 | RHS evaluations | 71505 |
| 87 | Minimum step [s] | 1 |
| 88 | Maximum step [s] | 65697735351472 |
| 89 | Stiffness conclusion | RKF45 adequate for this controlled benchmark; max steps/output 362 baseline / 364 refined; no late-time collapse; baseline rejection fraction 0.15489182462186482 |
| 90 | Refined trajectory completed? | YES; 10415 accepted, 2183 rejected, 86004 RHS evaluations |
| 91 | ODE convergence metrics | Tinf_K=1.7408584682611815e-05; eta_e_MeV=1.5866651419012378e-07; eta_mu_MeV=1.9479546821964215e-07; xi_e=1.7387445455634339e-05; xi_mu=1.7448183109581896e-05; LH=1.9340242731514387e-05; DeltaLnu=2.5373918741539931e-05; DeltaPbeta=0.00013316625368904651; frozen T/eta <=2e-4 PASS |
| 92 | Incremental sign crossing? | YES, both channels and summed beta power |
| 93 | Incremental crossing time / xi | Electron estimated 459453.1567 yr; muon 287057.1372 yr at xi=4.909710028924132. Summed beta power bracket 334965.4392..354813.3892 yr; separate from channel roots. All estimates are checkpoint-bracket/log-interpolation diagnostics. |
| 94 | Full-root crossing behavior | Both cross xi=5.633717467648343 later: electron estimated 486585.5972 yr, muon 305933.1006 yr; exact brackets in candidate. Summed full-power bracket: 354813.3892..375837.4043 yr. |
| 95 | Quasi-steady reached? | YES, 35 eligible checkpoints in fixed 1e9..1e10 yr window |
| 96 | 1/7 result | PASS |
| 97 | 1/7 metric | Slopes 0.1413043725241868 / 0.14166915596806406; asymptote errors 0.0047678802070558746 / 0.003784293384081039 |
| 98 | PHASE5D1_EVOLUTION_SHA | d3670f6d4e021def0483909b6d2fdeed1c6973a4 |
| 99 | Candidate artifact | docs/validation/phase5d1_controlled_evolution_candidate.json |
| 100 | Candidate artifact SHA256 | 47da15fa9e32be095a78d1b78ed17c3ffa14e5ef9be5528db3780a76afd0c079 |
| 101 | Retained response suite | PASS, 2/2 |
| 102 | Coupled suite | PASS, 3/3 registered tests plus recorded physical-run validation |
| 103 | Trajectory validation | PASS, raw rc=0; six complete runs and fixed-checkpoint comparisons |
| 104 | Final Phase-5B regression | PASS, 1/1, raw rc=0 |
| 105 | Final Phase-5C regression | PASS, 1/1, raw rc=0 |
| 106 | Complete data-free suite | PASS, all 50/50 actually executed within the complete authenticated invocation |
| 107 | Complete authenticated suite | PASS, 73/73; final test-only expansion rerun separately with raw rc=0 |
| 108 | Raw return codes | All required build, response/coupled/trajectory, pretrajectory, complete-suite and final governed gates: 0. Shared suite projections refer to the same actual complete invocation; the final expanded coupled test was additionally rerun. |
| 109 | Failures / skips | Zero in final required validation. Earlier development compilation errors were repaired; one development assembly was cancelled before any trajectory and is not counted as a test pass. |
| 110 | Final protected hashes unchanged? | YES, exact 33/33 |
| 111 | Ten governed baselines unchanged? | YES |
| 112 | Phase-5C candidate unchanged? | YES, a6b430b2356987b54a7c82c87d0b00f3e1d9698f6f299c32c80efcc40c63833b |
| 113 | EOS/data diff | NONE; authenticated bytes unchanged |
| 114 | Literature diff | NONE; authenticated bytes unchanged |
| 115 | INV-11b | IMPLEMENTED / CANDIDATE-VALIDATED |
| 116 | INV-11c | IMPLEMENTED / CANDIDATE-VALIDATED |
| 117 | INV-11d | IMPLEMENTED / CANDIDATE-VALIDATED |
| 118 | INV-11e | IMPLEMENTED / CANDIDATE-VALIDATED for frozen-v1 semantics |
| 119 | INV-11f | IMPLEMENTED / CANDIDATE-VALIDATED for controlled v1 |
| 120 | Global INV-11 | UNRESOLVED / AWAITING INDEPENDENT REVIEW AND OWNER RATIFICATION |
| 121 | Realistic normalization claimed? | NO |
| 122 | A18 begun? | NO |
| 123 | BNV begun? | NO |
| 124 | PHASE5D1_CANDIDATE_SHA | The commit containing this record with exact message docs: validate controlled rotochemical evolution candidate; resolved SHA is reported after commit in the final response (avoids self-referential commit hash). |
| 125 | Branch local/upstream/live equality | Post-push fetch/live verification is reported in the final response; this committed record precedes that operation. |
| 126 | Canonical unchanged? | YES, d019ae390be4f5e3daba05039903485cb497e397; verified before candidate commit; post-push verification is reported in the final response |
| 127 | Merged? | NO |
| 128 | Final disposition | A — PHASE-5D CONTROLLED NON-SUPERFLUID ROTOCHEMICAL EVOLUTION IMPLEMENTED AND CANDIDATE-VALIDATED — FIRST PROVENANCE-PRESERVING COUPLED T/ETA TRAJECTORY COMPLETE — READY FOR INDEPENDENT SCIENTIFIC REVIEW |
| 129 | Exact remaining caveats | Candidate evidence only; no independent review, owner ratification or canonical integration. Global INV-11 remains unresolved. Controlled mathematical/architecture normalizations; realistic FR2005/A18 remains source-limited and blocked. Frozen v1 coefficients and fixed background; no dot(Z), DU, superfluid or BNV evolution. Crossing times are checkpoint-bracket/log-interpolation estimates, not event-localized roots. Net thermal residual is a cancellation: its maximum relative baseline/refined difference is about 1.02 percent, while the predeclared T/eta convergence gates pass. Mutation evidence uses transformed production evaluations; aliases do not receive independent credits. Quadrature and radial diagnostics measure stability of the declared representation, not a continuum error certificate. Reproduction requires the authenticated frozen qualification inputs and entry manifest; missing inputs fail rather than skip or trigger requalification. Thermal entropy uses the predeclared fixed-background free-gas adapter; its extension outside degeneracy is mathematical, not finite-temperature EOS authority. Existing heat-capacity cache interpolation is retained. |
| 130 | Exactly one recommended next action | Run an independent scientific and numerical review of the complete provenance-preserving Phase-5D controlled-evolution candidate before owner ratification or canonical integration. The reviewer must independently verify frozen-context lifetime/currentness, the 33-path governed-source immutability gate, semantic Z/W consumption, same-Ltilde thermal bookkeeping, spin/reaction signs, component-scaled RKF45 tolerances, spin-only/reaction-only/Lyapunov oracles, first trajectory energy bookkeeping, stiffness and convergence, modified-Urca thermal sign crossings, quasi-steady 1/7 scaling, mutation coverage, and protected upstream regressions. Do not begin realistic A18 closure or BNV until the standard controlled rotochemical candidate passes that review. |
