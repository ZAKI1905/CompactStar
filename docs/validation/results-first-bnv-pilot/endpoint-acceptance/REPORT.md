# Conditional endpoint acceptance — restricted P2 disposition

**C-CUT and the restricted static matrix now pass conditionally. Not every
restricted-pilot readiness gate passes, so no paired-evolution authorization
request is issued.** The native source/control assembly, non-evolving snapshots
and restricted runner are still incomplete; nonlinear time-evolution accuracy
has not been established. See the complete [readiness matrix](READINESS.md).

## Owner decision and scientific scope

The approval is recorded verbatim in [OWNER_APPROVAL.md](OWNER_APPROVAL.md:1)
and the accepted endpoint amendment of ADR-0022. Authority was committed as
`b7f1209` before the assessor consumed it. Starting source was
`afd4b06074a5f597ef4cbe38ce3c9cb04386f556`, in the clean EKU worktree
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot`,
branch `candidate/results-first-bnv-pilot-20261108`. Every other local branch tip
was an ancestor; other worktrees and master were unchanged.

For this experiment only, C-ENDPOINT extends the previously accepted bands,
curvatures, positive baryon derivative and fixed topology over the required
fixed-B remapping endcaps. The retained eta is about 3.021576760982882e-26,
below the approved cap 3.022e-26. The physical interval stays [-1e-6,0]. This is
**an owner-accepted conditional scientific assumption, not a newly established
rigorous numerical bound**. It applies to no other star, regime or larger
interval. All S/C-REF/C-TAN/C-CONT/C-INPUT/C-SUPPORT/C-CUT/C-MUN conditions,
ADR-0022's thermal-fit amendment and metric-consumer limitation remain explicit
in every new row (`evidence/assessment/ledger.json.gz:1`).

The transition audit checked 30 retained endpoint/offset states and 12 positive
baryon secants. Every sampled state retained the same ordering of terminal,
neutron, envelope-base, muon and central enthalpies. Minimum sampled terminal-neutron
margin is 8.199277440638066e-4; minimum central-muon margin is about 0.01970708.
The validated terminal tube retains its positive topology margin and separated
pressure, fixed-B, horizon and chemical denominators. Retained production
selection indices are constant. **No known endpoint transition contradicts the
condition in the inspected evidence.** Sampled checks are not a theorem about
unsampled continuation. Continuous selection remains unproven and all six
endpoint grid-curvature checks remain inconclusive
(`evidence/assessment/endpoint-transition-audit.json:1`,
`../surface-transfer/evidence/assessment/endpoint-diagnostics.json:1`).

## Completed static result

| Disposition | Rows |
|---|---:|
| CONDITIONALLY_QUALIFIED | 420 |
| EXACT_IDENTITY, P2 | 21 |
| CONSUMER_SCOPE_ONLY, metric | 21 |

All 399 formerly endpoint-blocked approximate rows now have a complete
conditional C-CUT total. The other 21 conditional rows are species support.
The 21 metric rows remain restricted to the accepted consumer functionals;
there is no nodewise metric or full BA13 promotion.

Maximum full utilization is **0.782464607730239**, for N_mu. The maximum
C-CUT-added utilization is **1.940607143781908e-11**, for surface_gravity at
index 20. Both are evaluated against the unchanged historical thresholds.
Envelope and Tsurface_inf maxima are approximately 0.5202796022275094 and
0.5186840987342625. Complete results:
[matrix.tsv](evidence/assessment/matrix.tsv),
[ledger.json.gz](evidence/assessment/ledger.json.gz),
[summary.json](evidence/assessment/summary.json).

The assessor reuses the validated local witness without new stellar solves.
It propagates the finite fixed-B U through the paired increment, unchanged
all-row fit, baseline error, final attained-coordinate sliver and nonlinear
baseline denominator. For numerator N, baseline b and fit norm f it evaluates
(N+(f+2)U+endpoint)/(1-b-U), using fU at the origin and omitting the added
baseline U only for absolute tangents. Endpoint=2U*sliver/h at index 20.
No unsupported quantity is set to zero; actual-path uncertainty remains null.
The inherited ledger is recursively byte-value/type-equivalent, with only a
new successor entry attached (`tools/linux/results_first/endpoint_acceptance.py:110`, independent audit).

The **42 historical C-REF failures**, **355 original and 353 amended
utilization-change expectation violations**, six inconclusive endpoint grid
checks, and BA12/BA12R failures are preserved. The old unauthorized scenario
remains unauthorized **as a historical record**; the new entry derives its
conditional authority from the separately committed owner decision. Historical
FAIL/UNQUALIFIED rows are not rewritten. Runtime certificates, fits, thresholds,
production code and guard behavior remain unchanged.

## Verification and resources

Final numerical source: `b0a080db9f334d6825412e1fb3b102f5bd0dba7b`.
EKU job **19162**, phyn004, completed exit 0 with a clean exact-source preflight.
It used **34.027985 CPU seconds**, **37.948161 wall seconds**, and **269484 KiB**
maximum child RSS. Requested limits: one CPU, 1 GiB, five minutes wall,
120 CPU seconds. Affinity was one logical CPU; Slurm lists two allocation CPU
entries. Receipts: `evidence/receipt.json:1`, `evidence/controller-final.txt:1`.

- **27 numerical/authority/refusal tests pass**: seven cut identities, six local
  surface tests, nine original static tests and five endpoint scope/transition
  tests. Historical authority tests run against a copied historical ADR; the
  current accepted amendment and old manifests are both preserved.
- Historical ledger audit: **5,697 checks and five mutation refusals**.
- Historical transfer audit: **5,321 checks and five mutation refusals**.
- New independent audit: **5,804 checks and ten mutation refusals**. Independent
  mpmath arithmetic at 130 digits checks the Decimal 160-digit ledger with an
  expanded denominator correction. Mutations test baseline/endcap omission,
  zero cut error, changed threshold, historical promotion, lost condition,
  rigorous-bound promotion, metric promotion, nonlinear-readiness promotion,
  and replacing unknown actual-path uncertainty with zero. Mutations are
  in-memory copies only (`evidence/assessment/independent-audit.json:1`).

EKU prerequisite job **19159** built and passed all **13** affected native-stack
synthetic/analytic tests. Build: 81.632563 wall seconds, 224.174694 CPU seconds,
318904 KiB maximum child RSS. Tests: 1.862179 wall seconds, 1.585077 CPU seconds,
31488 KiB maximum child RSS. Four CPUs/8 GiB/15 minutes were requested; controller
runtime was 1m32s. Source at launch was `3b3f0a4`; later changes concern endpoint
accounting/launcher files only, not any compiled scientific source or those tests.
See `evidence/prerequisites/prerequisites.log:1`, `build.log:1`,
`controller-final.txt:1` and `artifact-hashes.json:1` in that evidence directory.
These tests establish their synthetic contract scope, not native paired snapshots
or actual-path R18/R20/nonlinear accuracy.

Earlier attempts are retained in `attempts/`:

- **19158** failed while a broad Git archive requested an unrelated omitted
  promisor blob from this sparse checkout, before numerical tests. The runner
  now copies materialized inputs and retrieves only the pinned old ADR.
- **19160** passed all tests and old audits, generated the matrix, then failed
  the independent relative check for the tiny utilization increase. Decimal
  precision was increased from 110 to 160 to preserve guard digits during
  subtraction; the independent 1e-95 relative arithmetic tolerance and every
  scientific threshold were unchanged.
- **19161** refused the combined clean/source preflight before numerical work.
  That version did not record which predicate failed, so no cause is asserted.
  The final runner records actual/requested commits and full status before
  refusing; **19162** passed. No preflight guard was waived.

The initial sandboxed Slurm submission could not resolve the controller and was
interrupted without submitting a job; network-enabled submission then succeeded.
No automatic-approval rejection remains unresolved. Slurm accounting storage is
disabled; controller and process receipts are retained instead. All substantial
artifacts remain on EKU; no Mac-side files were created.

## Readiness outcome and remaining work

The static nonlinear surface-transfer requirement passes conditionally. The
nonlinear **time-evolution** requirements are still unexecuted and must remain
trajectory-time acceptance gates. Before requesting those trajectories, native
matched pair preparation itself must be completed. The existing builder uses
the old input identity/Bdot and invokes new central/sequence solves; using it
here would violate the current restriction (`READINESS.md:1`). The successful
native ordinary-context binding test did not qualify a paired builder.

The next implementation target is authenticated retained-data reconstruction of
the owning central/sequence/tangent dependencies, followed by the already
predeclared source/control snapshots and a restricted native runner with all
historical refusals intact. A generic import declaration or a copied scalar
coefficient is insufficient evidence. No new scientific failure demonstrates a
need for more stellar calculations, and none was performed. Full readiness and
therefore the success condition for a paired-evolution authorization request
have **not** been met.

Files changed: accepted ADR-0022 amendment; new endpoint-acceptance authority,
reports and evidence; endpoint assessor/auditor/runners; endpoint refusal tests.
Classes: scientific qualification semantics, validation numerical-method,
engineering, documentation and generated artifacts. Governing authority:
GOVERNANCE §§1–5; ADR-0015/0017/0022; recorded owner decisions. Scientific impact:
new conditional uncertainty accounting and qualification status only; no physical
solver output, equation, stellar model or numerical threshold changed. This
remains an AI-authored scientific candidate within the explicitly accepted
conditional scope. No master merge, historical certificate promotion, new star,
BNV trajectory or temperature curve was produced.
