# Surface-transfer decision

**The finite interior C-CUT transfer is supported, conditional on the accepted
input bands. Full C-CUT qualification remains blocked by the fixed-B endpoint
extension. Restricted P2 is not yet ready for evolution authorization.**

The authorized local calculation has been completed. No global stellar
integration, 21-star refinement, thermal trajectory, certificate promotion or
master merge was performed. This is an unratified scientific candidate under
GOVERNANCE §§1–5; it does not introduce or approve a scientific condition.

## Result and exact remaining uncertainty

The complete finite terminal map and its first variation have been evaluated
using retained centers and offsets at targets 0, 10 and 20, for both signs of
|delta ln Ps|<=2.1e-12. The calculation transports structural, baryon, thermal
and susceptibility integrals, reconstructs the nonlinear temperature and
chemical inverse, and includes the central adjustment at fixed B. Tangents
use the accepted on-law count curvature identity. See [METHOD.md](METHOD.md).

The uniform finite bounds require an on-law starting coordinate x' with

    |x'-x| <= eta = 3.021576760982882e-26,
    |terminal log-pressure length| <= 2.100000000000254e-12.

These displayed limits are rounded upward; exact values are in
[bounds.json](evidence/assessment/bounds.json). Existing accepted C-CONT bands
cover [-1e-6,0]. One pressure sign at each endpoint needs the family beyond
that interval. No existing supported remainder establishes the requisite
continuation. **Smallness of eta is not evidence for that continuation.**

The retained offsets were reused. All six endpoint grid curvature diagnostics
are inconclusive under S; four width diagnostics are order-compatible and two
are inconclusive. These are sampled fixed-cut curvatures, not a uniform
on-law endpoint enclosure. Their failure is preserved in
[endpoint-diagnostics.json](evidence/assessment/endpoint-diagnostics.json).
It does not invalidate the accepted on-law identity inside its domain.

The missing quantity is a uniform bound on the on-law first derivatives of
the normalized observables and on the second derivatives of the species
counts, together with the existing pointwise bands, positive B_Hc and topology,
on the two remapping endcaps. The terminal pressure vector field and its
finite remainder are no longer missing. The baseline endpoint also enters
every relative normalization; absolute tangent rows still use the baseline
in their paired difference and fit.

## Observable-specific finite bounds

U below bounds the finite **fixed-B** change in the original normalization
when both the initial and final family coordinates are in the accepted
interval. Values are rounded upward for this table; the complete terminal
and remapping components are in the generated ledger. Tangents use absolute
units. Exact terminal zeros for Nn/Nmu are support identities, not missing
contributions: their nonzero remapping terms are included.

| Quantity | Interior U |
|---|---:|
| t_n | 1.618e-24 |
| t_e | 1.585e-24 |
| t_mu | 4.381e-26 |
| N_n | 2.948e-26 |
| N_e | 4.946e-24 |
| N_mu | 1.838e-24 |
| Ltilde_Me | 6.920e-26 |
| Ltilde_Mmu | 9.046e-25 |
| mu_B, P0_average | 9.103e-27 |
| P1_average | 1.170e-26 |
| radius | 1.297e-16 |
| Z_row_npe | 1.409e-20 |
| Z_row_npmu | 7.119e-23 |
| Cstar | 1.658e-19 |
| envelope | 6.756e-17 |
| Tsurface_inf | 7.849e-17 |
| surface_gravity | 2.703e-16 |
| mu_n_profile, maximum over original nodes | 1.580e-26 |

Pressure remains positive, D>0.8557, and the terminal region remains pe-only
(enthalpy margin to neutron onset >8.199e-4). Fixed-B remapping denominators
exceed 0.999999999999879; the chemical inverse denominator exceeds
0.99999999999999999998. These domain checks cover both pressure signs.
The local variation excess is bounded by 8.316e-13. Finite exponential,
rank-one inverse, structural remainder and tangent bounds are included.

## Qualification matrix and original budgets

The updated complete [462-row matrix](evidence/assessment/matrix.tsv) and
[uncertainty ledger](evidence/assessment/ledger.json.gz) retain all prior records
verbatim, including their numeric JSON types. The full-cut and full-utilization
fields remain null wherever the endpoint continuation is required.

| Disposition | Rows |
|---|---:|
| UNQUALIFIED_ENDPOINT_EXTENSION | 399 |
| EXACT_IDENTITY: P2 | 21 |
| CONDITIONALLY_QUALIFIED: support | 21 |
| CONSUMER_SCOPE_ONLY: metric | 21 |

The **42 original C-REF failures remain recorded**; their satisfaction of the
amended necessary partial test is unchanged, with no full-PASS promotion.
The **353 utilization-change expectation violations remain recorded**, with
unchanged row identities. Original thresholds, fits, the full ten-times C-REF
allowance, worst-case 2e increment, baseline normalization and all correlated
C-TAN/C-INPUT remainders are preserved. The existing thermal total-reference
fit bounds under ADR-0022 remain unchanged. There is no new B-dependent
shared-error assumption.

**The maximum full C-CUT effect on the original budgets is still unsupported.**
It is not the small interior number in the table. In the explicitly
unauthorized scenario extending the same bounds to the two endcaps, the
largest added utilization is **1.940607143781908e-11**, at surface_gravity row
20, and the largest total utilization is **0.7824646077302390** (N_mu).
This calculation includes the paired increment, all-row fit, final sliver
and nonlinear baseline denominator. It demonstrates the consequence of a
named condition; it neither supplies evidence for that condition nor authorizes
qualification. See [summary.json](evidence/assessment/summary.json).

## Smallest defensible next decision

For a narrower **conditional pilot claim**, the exact additional owner
condition would be:

> For RF-BNV-P2-20261108 only, assume that the already accepted reference
> pointwise bands, observable slope/count-curvature bounds, positive B_Hc and
> fixed topology continue along the same surface law to
> [-1e-6-eta, eta], with eta equal to the bound in this delivery's bounds.json.
> Retain the full finite C-CUT allowances and every previous condition. This
> extension is a scientific assumption, not a numerically validated endpoint
> result or a rigorous enclosure.

That is the smallest change to the scientific claim that the present evidence
can support discussing. It does not shorten the required depletion interval,
assume constant shared errors, change thresholds or launch an evolution.
It requires explicit scientific owner approval; this delivery does not select it.

If a numerical rather than conditional endpoint claim is required, the single
remaining evidence object is a **two-ended equilibrium-family variational
enclosure** on those endcaps, reusing the retained endpoint interiors. It must
bound Q_x and N_species,xx with a finite continuation remainder and separated
B_Hc/topology, then reproduce the transfer/ledger checks. Another terminal
pressure-flow calculation cannot supply that interior-family information.
No global re-shoot or 21-star campaign has been shown indispensable, and none
is proposed as an automatic next step.

The accepted metric consumer restriction remains in force. Native identity
bindings and Bdot reconciliation are not reopened. Source/control snapshot
preparation and a separate owner evolution authorization still follow static
qualification; actual-path R18/R20, nonlinear accuracy and signal resolution
remain future evolution checks, not results of this calculation.

## Provenance, validation and resources

Starting commit: `5498d6036b7369b06021bb7f98988961356291fd`.
Final numerical source: `30e24f488a6b0dbd5f1b7c24129256c67db8e184`, clean at
launch. Branch: `candidate/results-first-bnv-pilot-20261108`, worktree
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot`.
New files only: scoped local-transfer implementation/runner/audit, six tests,
this documentation and generated evidence. Change classes: validation
numerical method, engineering, documentation and generated artifact. Historical
sources, scientific policy, runtime behavior and all existing evidence are unchanged.

EKU phyn004 job **19157** completed successfully: **53.209438844 CPU seconds**,
**53.860937409 wall seconds**, **290192 KiB** maximum child RSS. One CPU was
requested and the process was pinned to one logical CPU; Slurm allocated two
CPU entries. Limits: 1 GiB, five minutes wall, 120 aggregate CPU seconds.
Controller/process receipts are retained; Slurm accounting storage is disabled.

- All **22 tests** pass: seven existing cut identities, six local map/interval/
  nonlinear inverse/refusal tests, nine existing static tests.
- **18 local maps** at each of 90 and 120 digits, each with 1/2/4 RK4 steps;
  six propagated variations. Maximum relative precision disagreement in local
  observable changes: **7.503e-65**. Of 324 observable convergence diagnostics,
  288 are order-compatible; 36 are inconclusive for the exactly unchanged
  neutron/muon terminal core counts. No numerical spread is used in place of
  the analytic finite majorant.
- Historical ledger audit: **5697 checks**, five mutation refusals.
- Independent Decimal successor audit: **5321 checks**, five mutation refusals,
  including baseline/endcap omission and unauthorized qualification.
- Earlier job **19155 failed** in interval/derivative operand dispatch before
  generating a map (0.733831753 CPU seconds). The fix was verified by job
  **19156**, which passed in 44.469975127 CPU seconds. Its preliminary successor
  serialized old numeric fields as strings; the final source corrects that and
  adds the independent preservation audit. Both attempts and their logs remain
  in `attempts/`; only job 19157 is the final evidence.

Total CPU across these three attempts is less than 99 seconds. Substantial
outputs remain on cluster storage; no Mac-side files were created. See
[REPRODUCE.md](REPRODUCE.md), input hashes and SHA256SUMS. The delivery commit
contains this report and the already generated final evidence; its exact ID
and clean-worktree check are reported at handoff.
