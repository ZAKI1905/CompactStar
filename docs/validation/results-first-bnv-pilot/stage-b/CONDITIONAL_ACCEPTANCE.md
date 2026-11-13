# Stage B conditional acceptance — OWNER RATIFIED (restricted)

**Owner conditionally ratified on 2026-11-12.** The owner explicitly accepted the
independent **B — CONDITIONAL PASS** disposition for the fixed-input,
full-retention mathematical P2 benchmark at the exact reviewed commit below.
C1–C7 are retained verbatim as the review restrictions; the original draft
language in C7 records why only the explicit owner decision below ratifies it.
This is not a microscopic BNV rate, astrophysical constraint, generalized
stellar-background certificate, or longer-horizon accuracy claim.

**Object of review:** `49253200a801ad9b336789dffb4e93d535a8e903`, branch
`candidate/results-first-bnv-pilot-20261108`, descended from accepted Stage A
`c90b2205e553afcd102c7acd5f84b253e0a67f1b`. The exact original review is
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/stage-b-independent-review-20261111/REVIEW.md`;
its unchanged copy and selected evidence are bound by `review-evidence/IMPORT.json:1`.

**Result offered for acceptance:** source minus exact-zero control endpoint
ΔT∞ = +15,325,803.204209456 K, conditional HP enclosure 7.47e−9 K; original
24 configurations, 24 schedule repeats, 96 diagnostic grids, six restarts,
ten native references and 58/58 aggregate checks. Native endpoint bounds remain
separate from the HP enclosure. Gross heating +5.19 µK, net weak reaction
−2.59 µK, actual-potential storage −70.5 µK, combined chemistry **−73.1 µK**.
The chemical values are reviewer term-suppressed diagnostic integrations; only
the isolated gross/net reaction values also have candidate analytic intervals.
Sources: `review-evidence/REVIEW.md:41`, `:281`, `:397`; `evidence/qualification.json:1`.

| Condition | Exact acceptance restriction |
|---|---|
| C1 — Fixed inputs | Cstar, Urca coefficients, tangent, photon/direct anchors and frozen background accepted as given. Depletion 2.048e−9; maximum frozen utilization 0.001108 (the review's ~1.1e−3). BA13/static uncertainty is not requalified. |
| C2 — Model | Prescribed full-retention P2 sink, not a microscopic BNV decay or astrophysical constraint. Numerical precision concerns the declared real RHS; omitted finite-temperature weighting has ~6.26 K scale. |
| C3 — Arithmetic | 110-digit evaluation and 1e−80 per-step/PI padding are premises, corroborated independently, not formally machine-proved. |
| C4 — Domain/horizon | Exactly 2,048 Julian years (64,629,964,800 s), x ∈ [−0.03,0.2], \|η_e\| < 4e−7 MeV, \|η_mu\| < 8e−7 MeV; exp(LH)<2. Muon imbalance uses ~66% of its box. Regenerate constants for any extension; no automatic longer-horizon acceptance. |
| C5 — Events | Only two transversal upward source crossings, knots 94 and 95. No downward, near-tangential or control crossings. Native knot offsets and old-cell stage extensions are bounded as described in PIECEWISE.md. |
| C6 — Chemistry | Microkelvin results are fixed-input mathematical diagnostics. Isolated analytic reaction bounds exclude storage. R20 cannot test the chemical ledger; pointwise identities and independent model agreement support it. |
| C7 — Owner | Only the owner can ratify under GOVERNANCE §5. This draft, the independent AI review and Git commits do not ratify Stage B. |

Authority: `review-evidence/REVIEW.md:417`; numerical utilization details:
`evidence/final-summary.json:1`. Mathematical bounds conditional on premises,
independent numerical estimates, and fixed-input/model-form limitations retain
the distinct classifications in `review-evidence/REVIEW.md:397`.

D1–D4 are addressed in [REPORT.md](REPORT.md) and [PIECEWISE.md](PIECEWISE.md):
complete chemistry table; native knots/extension strips; R20 decomposition and
Stage A/B distinction; three separate uncertainty classes next to ΔT. The R20
source/control totals remain 4.793855e−5 / 4.036164e−5 N against 2e−4 N.
No curvature-bound replacement, new numerical campaign or threshold change is made.

The reviewer re-executed assessment code, but did not regenerate native traces,
rerun preflight contracts/unit tests, or rerun event-neighbour guards. Native
references share production RHS/anchors; Decimal and reviewer mpmath models are
independent implementations of the same declared model. These qualifications
remain part of acceptance (`review-evidence/REVIEW.md:406`, `:495`). Historical
BA12 and all BA12R failures remain failures (`review-evidence/REVIEW.md:141`).
No master merge, certificate promotion, new physics or next-study approval follows.

**Owner decision:** CONDITIONAL RATIFICATION of Stage B, retaining C1–C7.
**Ratifying authority/date:** project owner, explicit instruction of 2026-11-12.
**Exact reviewed commit:** `49253200a801ad9b336789dffb4e93d535a8e903`.
**Accepted endpoint:** ΔT∞(2048 yr) = **+15,325,803.204209456 K**.
The recording commit is discoverable with `git log --follow -- this-file`;
it is distinct from the reviewed commit. No numerical rerun was performed.
Consistency verification: [entry-authentication.json](../../../research/bnv-rotochemical-competition/evidence/entry-authentication.json).
The corrected report, its evidence, and the imported review were checked before
recording this decision; the reviewed production source, tests, tools and Stage B
numerical evidence are unchanged at the task entry `81b808f`.
