# Phase-6A-1 controlled BNV numerical failure evidence

Classification: **CONTROLLED MATHEMATICAL / ARCHITECTURE BNV FAILURE EVIDENCE;
NOT A PHYSICAL BNV MODEL; NOT A GOVERNED BASELINE; NOT AN OWNER-RATIFIED
NUMERICAL RESULT**.

Disposition: **C. CONTROLLED TRAJECTORY / NUMERICAL VALIDATION FAILED —
CANDIDATE NOT ACCEPTABLE**.

## Authority and immutable inputs

- canonical entry: `961dfa0de6f76df71df4cb98edc8e1b35a5c21b1`
- predeclaration: `6a3c3e8653f407ba9734968bb7a1802cb78a5549`
- pretrajectory gate: `460ca713a6bd2633b4afbf3926e959b751cea0c6`
- trajectory execution implementation: `6fc3b067365725064bc9debda9489997806b8d77`
- numerical validator implementation: `b6df489c7ab7c15a3b92493618ae0fc3246e9a57`

The four exact run cards, baseline/refined RKF45 tolerances, output grids, R20
budgets and frozen boundary are unchanged from the predeclaration. No physical BNV
rate, lifetime, mixing, matrix element, particle mass, decay channel or model was
introduced.

## Execution evidence

The authoritative raw execution completed all four declared source runs and all
matched controls at both tolerance levels:

| card | rows per source/control and tolerance | baseline accepted source steps | refined accepted source steps |
|---|---:|---:|---:|
| `RF-P0-TRANSIENT-v1` | 1025 | 1056 | 1056 |
| `CPL-P0-TRANSIENT-v1` | 2049 | 2090 | 2090 |
| `CPL-P1-TRANSIENT-v1` | 2049 | 2090 | 2090 |
| `CPL-P2-LINEAR-QSS-v1` | 8193 | 8318 | 8322 |

The scratch execution log SHA-256 is
`c4e565f905c47db7216624c958a765a538cc2b6b6fa3b0f2f6f9f560a56b88aa`.
The 32-file raw-evidence manifest has aggregate SHA-256
`81314ddceba77eaaf29f493928560e93f3a6b45d59f302c11fbec0e799da0945`.
These build-tree files are failure evidence, not a candidate artifact.

All serialized checkpoints remained inside the frozen certificate. The maximum
observed utilization was `0.27037946636901355`, limited by `N_mu`, at the P2 final
checkpoint. The maximum `abs(DeltaB)/B0` was `4.999999999602668e-7`; maximum
`abs(DeltaN_mu/N_mu)` was `2.7037908300962217e-5`. No frozen quantity failed, so
there is no first frozen-failure time.

## Mandatory BA12 stop

The accepted BA12 metric is

```text
abs(y_baseline-y_refined)
----------------------------------------------- <= 1.
atol_baseline + rtol_baseline max(abs(y_baseline),abs(y_refined))
```

For `CPL-P2-LINEAR-QSS-v1`, the maximum was
`1.7255917120989046` in `x_state` at
`t=23113476562.5 s = 732.421875 yr`, with baseline value
`0.051132258160285306` and refined value `0.05113226698535251`.
The independent compiled comparator also stopped on `x_state`. All other cards
and all matched controls passed their component comparison; the P2 matched control
was bit-identical across tolerances.

The primary validation log SHA-256 is
`d89a8df15092db9fee93144b9e2e5b7fca4a53bc29a43e0652fc86f4f781dc33`.
The independent BA12 log SHA-256 is
`f506f5fbf957e97ff83140a4bfd38fcfaf39a36a134ede33cd9607e7c0931e38`.

This is an immutable-threshold failure, not authority to alter the run card or
solver tolerances. R20 candidate acceptance, BA16 QSS classification and the
trajectory form of BA17 were not continued after the stop. No file named
`docs/validation/phase6a1_controlled_bnv_candidate.json` was generated.

## Gate disposition

- BA1-BA10b: PASS before trajectory, with M1-M21 detected under their accepted
  primary/supplementary assignments.
- BA11: raw coupled execution completed; full R20 candidate acceptance was not
  completed before the BA12 mandatory stop.
- BA12: **FAIL**, exact evidence above.
- BA13: static certificate PASS; runtime checkpoints remained valid through the
  failed numerical experiment.
- BA14: raw owner serializer/schema checks completed; no candidate artifact was
  eligible for final BA14 acceptance.
- BA15: repository protection/regression evidence is recorded separately in the
  final campaign report; it cannot cure BA12.
- BA16: NOT RUN after BA12 STOP; no QSS label.
- BA17: pretrajectory mathematical bound PASS; trajectory diagnostic NOT RUN
  after BA12 STOP.

No candidate numerical result is retained or interpreted. The exact recommended
next action is to return the BA12 convergence conflict to the owner for a new,
separately accepted numerical plan; do not retune this campaign in place.
