# Final targeted C-CUT disposition

**C-CUT remains incomplete; restricted P2 is not ready for evolution
authorization.** The actual fixed-B response has now been derived and evaluated
from retained data. Its first-order size is extremely small, but the required
finite-cut remainder has not been established. No missing allowance was assigned
zero and no qualification row was promoted. The next action is the single local
surface-transfer witness in [NEXT_CALCULATION.md](NEXT_CALCULATION.md), reusing
the existing stellar interiors. New global stellar solves are not justified.

## What is now established

[DERIVATION.md](DERIVATION.md) supplies every observable's terminal response,
the implicit central adjustment `Hc_s=-B_s/B_Hc`, and the fixed-B response
`Q_s|B=Q_s|Hc-Q_Hc*B_s/B_Hc`. It includes the thermal nonlinearities, redshift,
chemical inverse, exact identities, historical scales, baseline normalization,
frozen fit, and correlated path-tangent response. Fixed-central zero shell terms
for neutron/muon counts are explicitly not fixed-B zero responses.

At the anchors, Hc_s is about 5.7126e-15 per unit ln(Ps). The corresponding
first-order central adjustment at 2.1e-12 is about 1.20e-26. The actual original
pressure law is replayed: a=8.383928153884364 and maximum log-pressure residual
2.0676793610618915e-12. Source: `evidence/response.json.gz:1`, derived from the
original profile terminal pressures and retained achieved coordinates.

The following are maxima over the three retained targets and original component
scales. The second column is |partial Q/partial ln(Ps)|_B / S; tangents have S=1.
The last column is its product with 2.1e-12, **not a completed error allowance**.
Full signed, dimensional and component records are in `evidence/response.json.gz`;
the complete machine summary is `evidence/summary.json`.

| Quantity | Sampled normalized sensitivity | First-order product only |
|---|---:|---:|
| t_n (path) | 2.47296e-14 | 5.19321e-26 |
| t_e (path) | 2.66636e-14 | 5.59935e-26 |
| t_mu (path) | 1.93408e-15 | 4.06158e-27 |
| N_n | 1.39819e-14 | 2.93620e-26 |
| N_e | 2.20680e-12 | 4.63428e-24 |
| N_mu | 7.77892e-13 | 1.63357e-24 |
| Ltilde_Me | 3.29399e-14 | 6.91738e-26 |
| Ltilde_Mmu | 3.62499e-13 | 7.61248e-25 |
| mu_B | 1.55532e-15 | 3.26617e-27 |
| P0_average | 1.55532e-15 | 3.26617e-27 |
| P1_average | 2.60821e-15 | 5.47724e-27 |
| radius | 6.17179e-5 | 1.29608e-16 |
| Z_row_npe | 6.60715e-9 | 1.38750e-20 |
| Z_row_npmu | 3.13014e-11 | 6.57330e-23 |
| Cstar | 7.89229e-8 | 1.65738e-19 |
| envelope | 3.21596e-5 | 6.75351e-17 |
| Tsurface_inf | 2.69567e-5 | 5.66091e-17 |
| surface_gravity | 1.28639e-4 | 2.70141e-16 |
| mu_n_profile (equilibrium scope) | 1.55532e-15 | 3.26617e-27 |

P2 remains exactly zero. Species support remains conditional on accepted
topology. Metric rows retain their consumer-only limitation.

## Evidence classification and precise stop

- **Analytic identities:** terminal TOV/redshift/integral derivatives, central
  fixed-B correction, chemical inverse derivative, and the alternative on-law
  expression for the path-tangent sensitivity. Tested independently; these are
  identities, not error bounds on retained stellar states.
- **Convergence-supported sampled scalar sensitivities:** all 16 non-tangent
  groups pass the three-width second-order and three-grid fourth-order
  diagnostics at all three targets. This supports the sampled derivatives;
  it does not establish a finite transverse remainder over the whole window.
- **Inconclusive mixed tangent estimates:** all nine tangent grid-order checks
  fail S's signed/order admissibility test; the three neutron width checks also
  fail it. Raw increments and orders are retained. This does not invalidate the
  previously accepted C-TAN estimates: these are new pressure sensitivities of
  tangents, involving an additional mixed derivative.
- **Accepted conditionals retained:** C-REF, C-TAN, C-CONT, C-INPUT, C-SUPPORT,
  C-MUN, ADR-0022 and the metric restriction. Their domains are unchanged.
- **Unsupported:** a uniform finite transfer U_Q(B) and U_T(B) across
  |delta ln(Ps)|<=2.1e-12 with fixed-B matching, baseline and endpoint coverage.
  Thus full C-CUT terms, their normalization and endpoint dependencies stay null.

The exact missing quantity can be expressed as the finite map remainder after
`delta_s*D_Q`: its Taylor form needs `partial_s D_Q` in the cut strip (and the
corresponding tangent response). No retained pressure-direction levels establish
it. The existing full Ps-to-zero shell bounds concern a different displacement
and fixed-central geometry. The recorded small pressure residual is an input
bound, not an output bound. No multiplier or threshold headroom substitutes for
the missing response.

There is a narrower route than new interior mixed derivatives: the on-law
identity in DERIVATION.md expresses D_T using T_x and derivatives of the terminal
baryon response. Accepted C-CONT count curvature bounds T_x on the law. This
removes the need to rescue the noisy N_hh-t*B_hh Richardson estimate; a validated
finite terminal map and fixed-B remapping are still required. The assessment
does not assert that C-CONT bounds an arbitrary off-law Hessian.

## Qualification and execution

The historical ledger is unchanged and was independently re-audited:
399 approximate rows remain unqualified, 21 exact P2 rows, 21 conditionally
qualified support rows, and 21 metric consumer-only rows. There is no new
qualification matrix because no supported finite C-CUT allowance was obtained.
The authoritative current matrix remains
`../static-completion/evidence/matrix.tsv.gz`, with the complete uncertainty
ledger at `../static-completion/evidence/ledger.json.gz`.

All 42 original C-REF failures remain historical failures and still satisfy only
the amended necessary partial C-REF test. All **353 utilization-change expectation
violations** remain explicitly recorded in the earlier ledger. No threshold,
fit, historical certificate, production tangent or scientific policy changed.
Source: `evidence/independent-ledger-audit.json:1` and the unchanged ledger.

EKU Slurm job **19154**, phyn004, completed exit 0 at execution commit
`c7a406803f57603b7c623404acfcd0eba7cc7115`. It ran seven independent C-CUT tests,
nine existing static refusal tests, 5,697 independent ledger checks and five
successful in-memory mutation refusals. Postprocessing at 70 and 90 digits
agrees to relative 2.58e-58; this does not add precision to retained stellar
outputs. Receipt: 6.951918153 CPU seconds, 7.261055 wall seconds, 241744 KiB peak
child RSS. One CPU was requested and affinity was one CPU; Slurm reports two CPU
allocation entries. Accounting storage is disabled; controller completion and
process receipts are retained instead. No solver or trajectory ran.

Earlier attempts are retained, not hidden: 19150 exposed the literal rounded-h0
convention in the new identity tests; it was corrected by differentiating the
literal EOS. 19151 completed the preliminary analysis. 19152 refused the
pressure-law replay because a shell-witness reconstructed pressure had been
extracted in place of the historical profile pressure; the final package keeps
both, uses the historical one, and reproduces the accepted residual exactly.
19153 refused its source preflight before numerical work; its combined guard
did not record which predicate failed, so the cause is not asserted. The final
runner records both commit IDs and clean status; 19154 passed that guard.

The delivery adds only assessment source/tests, derivation, input/provenance
snapshots and generated evidence. Change classification and governing authority
are in PREDECLARATION.md. All earlier tracked files are byte-identical to the
starting commit. The final delivery commit follows the execution commit with
documentation/evidence only; its exact ID and clean status accompany the handoff.
No Mac worktree or substantial Mac file was created. No merge or promotion was
performed. Scientific status remains **UNVERIFIED SCIENTIFIC CANDIDATE**.
