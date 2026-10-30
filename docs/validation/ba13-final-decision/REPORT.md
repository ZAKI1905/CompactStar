# BA13 final evidence and qualification decision

**Disposition E: the available evidence does not complete full BA13 qualification.**
The scientific decision is ready for review; the certificate is **not ready for
independent acceptance as a qualified BA13 certificate**. No measured result
establishes disposition C (historical production accuracy cannot meet the original
requirement). Option A remains numerically plausible. This is an UNVERIFIED
SCIENTIFIC CANDIDATE, with no owner ratification, promotion, merge or trajectory.

The refusal is not 441 unrelated problems or a demand for rigorous interval
arithmetic. Four common missing models account for the approximate rows: joint
reference bias, joint tangent remainder, physical continuation, and nonlinear
input/normalization propagation. Full metric_structure also lacks nodewise
reference data. Their necessary accuracy is quantified below. More central-star
refinement would not supply those models. The mandated evidence stop was reached;
no further stellar calculation or observer was run.

Evidence: `evidence/summary.json:1`, `evidence/tangent-reassessment.json:1`,
`evidence/conditional-intervals.json.gz`, and the 462-row matrix. Derivations,
quantity scope and exact missing components are in `METHODS.md:1`.

## Authenticated basis and protected records

Entry is **84e1062db78b91e4fb0988074941b1aea6019917**, on
`candidate/linux-ba13-evidence-20261026`, clean at entry. Work proceeds only in
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/ba13-final-decision`,
on `candidate/linux-ba13-final-decision-20261030`. The first commit d5b40b3 fixes
the claims, norms, error model, tests and 600 CPU-second ceiling before new
numerical assessment. Executable assessment was frozen in 5957154; 328ccf4 adds
methods and strengthens the missing-error mutation. The final delivery commit is
the head named in the verified incremental bundle and handoff receipt.

All **1,993 manifest entries** checked (including overlap between manifests) pass
SHA-256 authentication; the original certificate is
`4212cad05bb4244ab75e80b8dfea023e6d0fcee69086c9993bd48d26a40e41fb` and fitter is
`9b3ef1ce17df7e1f8d2d7c40e111ee0844abfddc2cc54b8f93cb3f2ce1322370`.
All other local branch tips at entry are ancestors. This is content/ancestry
authentication; it does not assert cryptographically signed authorship.
The latest available independent review is **2026-10-20, reviewing 093d121**;
there is no located independent review of 84e1062 or this package. Replaying that
review's numerical script does not constitute a new independent review.
Source: `evidence/authentication.json:1`.

Original X1 failures, X2 INCONCLUSIVE reference qualification, double/Newton
N512 failure, historical BA12/BA12R failures and every older BA13 outcome are
unchanged. The historical producer, fitter, certificate, production physics,
accepted dependency locks and canonical master are unchanged. Change classes:
validation numerical-method, engineering, documentation and generated evidence.
Authority: task authorization, GOVERNANCE §§1–5, INV-02/03/04/06/14/15 and accepted
preflight §11. ADR-0021 remains PROPOSED. Production numbers do not change.

## Actual new work and what it establishes

No new stellar evaluations, solver development, nodewise observer, full BNV
trajectory or thermal regression. On one pinned Slurm core, jobs **19138–19140**
performed only the predeclared static assessment and controls: **8.363868 measured
CPU-seconds**, including test subprocesses, versus a 600-second ceiling.
The initial assessment benchmark used 2.657819 CPU s; expansion was permitted
only after it passed. Scheduler accounting storage is disabled; allocation and
process CPU receipts are retained. The initial sandboxed scheduler submission
could not resolve ekucn1 and was interrupted before submission; the permitted
scheduler retry succeeded. No scientific run failed in this task. Intentional
mutation failures and the original control version are retained.
Source: `evidence/resources.json:1`, `evidence/receipt-benchmark.json:1`,
`evidence/receipt-validate.json:1`, `evidence/receipt-controls-v2.json:1`.

- **1,042 assessment checks pass:** all historical row identities, 21-target
  production/reference comparisons, raw tangent ratio/order checks, frozen
  holdout criteria and 20-gap coverage records.
- **55 tests pass:** all 43 inherited tests plus 12 new tests. The strengthened
  12-test suite was rerun after its mutation correction; this is not 67 unique tests.
- **Seven source mutations are detected:** missing error replaced by zero, wrong
  evidence authority, blanket correlation cancellation, wrong reference,
  sampled-curvature promotion, omitted denominator propagation, missing dependency.
- The independent review's 50-digit chemical inverse calculation replays
  successfully, including both endpoint Neumann remainders.

These checks establish arithmetic, provenance and refusal behavior. They do not
create scientific authority or complete missing physical uncertainty estimates.

## Scientific findings and required accuracy

The existing twelve offset-H calculations and N512 endpoints independently
reproduce all nine retained tangent ratios. Observed orders span 3.9483–4.0031;
largest finest-grid tail estimates are 1.57273e-14 (t_n), 1.26247e-14 (t_e),
2.68173e-15 (t_mu), all in absolute count/count units. Those estimates cover the
correlated GL4 discretization component. They do not encompass an arbitrary common
reference error's derivative, higher stencil remainder or moving-cut response.
The exact missing joint residual is |delta N_span − t delta B_span|, together
with a separated denominator and stencil/surface correction, not independent
endpoint tails or automatic cancellation.

Both prospective holdouts at targets 5 and 15 remain successful for all 15
non-tangent groups with their own achieved B and pressure cuts. Their largest
criterion ratios remain 0.117275 and 0.116597. Every achieved-B gap has an explicit
conditional remainder. Neither two holdouts nor three sampled curvatures supply
continuous authority. The attained endpoint misses -1e-6 by approximately
9.94e-16 in x, far below the existing coordinate uncertainty; it belongs in the
same continuation/input model, not a new high-precision root campaign.
Sources: `docs/validation/ba13-evidence/REPORT.md:45`, new tangent and interval
artifacts, and `evidence/summary.json:1`.

The matrix retains **21 exact-P2 PASS, 0 FAIL, 441 UNQUALIFIED**. The largest
finite added ledger consumes **3.07495%** of historical headroom (Ltilde_Me,
mostly the measured production representation discrepancy). N_mu retains the
highest historical utilization, **0.54075348**. Finite estimates are not complete
physical budgets. Complete added uncertainty and proposed full utilization remain
null for all approximate/support rows.

The following are *single-axis requirements*, not estimates. They are sufficient
caps on an **additional uniform pointwise error** if the existing finite terms
are valid and other additional gap/fit errors are zero. All unknown contributions
must share the full joint inequality in METHODS; these caps cannot be allocated
independently and then summed.

| Group | Historical headroom | Remaining pointwise-error cap | Units |
|---|---:|---:|---|
| t_n | 9.62563e-5 | 1.33972e-5 | absolute dN_n/dB |
| t_e | 2.89156e-6 | 4.02519e-7 | absolute dN_e/dB |
| t_mu | 2.10670e-7 | 2.93318e-8 | absolute dN_mu/dB |
| N_mu | 4.59247e-5 | 6.39337e-6 | historical scaled norm |
| Z_row_npmu | 6.49078e-5 | 9.03708e-6 | historical row norm |
| Ltilde_Me | 9.87733e-5 | 1.33289e-5 | historical scaled norm |

Thus the evidence gap is **not a need for 1e-15 accuracy**. A much coarser,
well-supported quantitative joint estimate could suffice. Conversely, unused
headroom is not evidence that an unknown contribution fits inside it. The full
22-group summary and every row's joint accuracy region are in `MATRIX.md:1`.

## Few common blockers, precise missing evidence and importance

| ID | Missing evidence | Expected importance and smallest route to resolution |
|---|---|---|
| C-REF | Joint output-level X2 error model covering stage/residual propagation, EOS/inversion, centre, arithmetic/summation and common bias, in each quantity's actual scope | Existing primary-output ablations/cross-method differences are tiny compared with headroom; likely small under a justified model, but no quantitative bound follows from that comparison alone. Supply a validated residual/sensitivity and shared-bias argument or an explicit supported conditional model. No new stellar solver is required by the proposed standard. |
| C-TAN | Bound/validated estimate for the joint span error, stencil remainder and fixed-cut versus moving-cut derivative, with achieved-B alignment | GL4 tails and same-method geometry residuals are many orders below the absolute caps. The unresolved part is their transfer to the complete derivative, especially cancellation. A validated remainder argument can suffice; another general solver campaign is not warranted. |
| C-CONT | Quantity-specific continuous curvature envelope, positive B_H, fixed species topology and surface-cut derivative/mixed-response model, including endpoints | Sampled topology margins, stencil stability and both holdouts strongly support smooth continuation. They cannot quantify arbitrary unseen behavior. Owner acceptance must name the numerical-continuation assumptions; support/positive lepton tangents remain separate exact/inequality requirements. |
| C-INPUT | Complete achieved-coordinate/reference normalization error and nonlinear fit design/weight/input remainder, including pressure conversion | Fixed-projection propagation is algebraically settled and all 21 rows are used. Missing input ranges/response conditioning prevent completion. They may be bounded jointly within the reported accuracy region; tests or covariance cannot supply the inputs. |
| NODEWISE | Independent r,m,nu,lambda values at the original 2049 support coordinates, moving surface included, with interior interpolation validation | Magnitude is not known. P−R about 1.27e-6 is only representation error. The passive observer is gated by C-REF/C-CONT evidence, not just missing owner signature, and was not run. |

The chemical inverse majorant is verified but alternative to direct X2 comparison;
it is not double-counted. Cstar's absolute cache bias up to 6.55536e-4 is kept
separate from its exact common-factor cancellation at equal thermal knots.
mu_n_profile uses the equilibrium redshifted identity only within its actual
eta=0 reference-field definition. All species numbers/tangents, both Z rows and
Urca coefficients, all thermal knots, gravity/envelope/temperature, mu_B/P0/P1,
radius, full metric structure and support remain in the inventory.
Sources and equations: `METHODS.md:106` onward; inherited scope in
`docs/validation/ba13-uncertainty/SCIENTIFIC_REVIEW.md:90`.

## Owner decision and next action

`OWNER_DECISION.md:1` is the exact proposed ratification record for **Option A,
interpretation (2)**: A full-interval coverage; B analytical bounds plus independently
validated numerical estimates and explicit conditional acceptance; C physical
equilibrium response including historical production error; D authenticated
sidecar; E independently checked achieved-B reference transfer. None is marked
ACCEPTED. A/C clarify the old scientific contract; B is evidence policy, D record
format, E validation method. Investigation category shares do not become new
historical physical thresholds.

A is not supported because the construction is incomplete. B alone would imply a
specific already-defined narrow calculation closes the gap; the joint models are
not yet specified, so that claim is not justified. C is not supported by the
measured discrepancies. Restricted sampled consistency and exact identities are
defensible, but not a completed D uncertainty certificate. **E is therefore the
full-claim disposition, with Option A still viable subject to the named evidence.**

The smallest next action is **one independent scientific review of this decision
package and the quantitative C-REF/C-TAN/C-CONT/C-INPUT models, alongside the exact
owner selections**. That review should either supply/accept supported numerical
models within the explicit joint accuracy regions, or select a clearly narrower
claim. It should not request another general refinement campaign. If the models
are supported, the remaining bounded calculation is the predeclared passive
nodewise observer with an independent interior-interpolation benchmark. Only a
subsequent complete sidecar, independent acceptance review and explicit owner
approval could enable the already scoped certificate-dependent controlled-BNV
pilots. This package authorizes no physical-rate selection or long evolution.

Reproduction and office-Mac bundle import are in `REPRODUCE.md:1`. Large raw matrix
and logs remain external under the SHA-256 manifest. Changes are limited to this
new validation package/tool/tests and an append-only ADR proposal; no historical
scientific record is edited or retroactively promoted.
