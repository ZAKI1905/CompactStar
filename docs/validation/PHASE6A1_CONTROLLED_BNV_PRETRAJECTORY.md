# Phase-6A-1I controlled abstract BNV pretrajectory gate

**Disposition:** PRETRAJECTORY PASS. The owner-accepted mathematical run cards
may be executed after this record is committed. At the time this record was
written, no BNV trajectory had been generated.

This is an architecture and validation record for a controlled mathematical
experiment. It is not a physical BNV model, rate, lifetime, mixing parameter,
matrix element, decay specialization, governed numerical baseline, or
owner-ratified numerical result.

## Authenticated state

- canonical entry and unchanged canonical master:
  `961dfa0de6f76df71df4cb98edc8e1b35a5c21b1`;
- predeclaration commit:
  `6a3c3e8653f407ba9734968bb7a1802cb78a5549`;
- production semantic commits:
  `ac4a858c956ed1f0b171418ccd698adc69d6431e` and
  `a96d366207b9c5a209afdf1fb816f2e5c5af8e34`;
- validation implementation commits:
  `fe0f4c03ffa72609c78338412b99657fb4886055` and
  `89edd398b4404ff69aab81c20297997f07c14c33`;
- predeclared temperature-grid completion:
  `95d162eb06e02644e6752bde3665689f22dc2355`.

The committed run-card drive values, durations, depletion targets, ODE
tolerances, checkpoint grids, target-B tolerance, R18/R20 definitions and all
acceptance thresholds are unchanged. The authorized cards remain:

| Run card | Fractional drive | Duration | Predicted depletion | Partition |
|---|---:|---:|---:|---|
| `RF-P0-TRANSIENT-v1` | `-1e-13 yr^-1` | `1e6 yr` | `1e-7` | P0 |
| `CPL-P0-TRANSIENT-v1` | `-1e-13 yr^-1` | `1e5 yr` | `1e-8` | P0 |
| `CPL-P1-TRANSIENT-v1` | `-1e-13 yr^-1` | `1e5 yr` | `1e-8` | P1 |
| `CPL-P2-LINEAR-QSS-v1` | `-1e-12 yr^-1` | `5e5 yr` | `5e-7` | P2 |

No physical rate or physical channel input was selected.

## Static frozen certificate

The final certificate uses 21 achieved target-B stars at the committed nominal
grid from zero through `-1e-6`. Every target residual and final bracket width is
at most `5e-11 B0`. The independent comparison uses achieved `B_solved` values,
the full outer-half indices 10 through 20, uncertainty-weighted centered linear
fits and the fixed criterion
`abs(residual) <= 3 u_res + 0.10 T_X`. No sample or fit window was removed.

Every target star recomputes `t`, `Z`, enabled `Ltilde`, direct-energy
quantities, particle numbers, source-support metric/potential profiles, 13 exact
thermal knots, surface gravity and envelope outputs through the named governed
owners. Each star's unclamped thermal owner consumes a fresh current-star
free-gas entropy adapter made by the unchanged Phase-5D formula. The candidate
run remains frozen on the authenticated B0 adapter.

- static certificate: PASS;
- maximum utilization: `0.54075347894390757`;
- limiting sample/quantity: index 20, `N_mu`;
- runtime just-inside ceiling: PASS at
  `DeltaB/B0=-9.9999997155751408e-7`;
- runtime just-outside ceiling: refused at
  `DeltaB/B0=-1.0000000283979203e-6`;
- stale and over-budget certificates: refused;
- M20: detected before RHS publication.

Certificate evidence hashes:

- JSON: `6643054d8cba36035af15a15615baf9920d9a56adf3f5211f13ea9f2745f4409`;
- production-loader TSV:
  `9217e278eb0a4b06b193a2a02f4d59f285d74e4033427b87ccd51181c3f7beab`.

Intermediate stopped validation-tool attempts are preserved in build scratch:
one initially sampled outside neutron source support; one successful sequential
star was interrupted solely to shard the expensive unchanged thermal-cache
construction; one correctly refused the B0-only thermal table on perturbed
stars. None contributed a row to the final certificate, and no refusal was
bypassed with clamping or extrapolation.

## BA1 through BA10b

| Gate | Result |
|---|---|
| BA1 source semantic atomicity/currentness/refusals | PASS |
| BA2 typed Phase-5B tangent, canonical `B_B=B_n+B_p`, closure/budget | PASS |
| BA3 sole moving-reference projection, baryon null and exact lift | PASS |
| BA4 neutron-sink signs and arithmetic | PASS |
| BA5 physical sliding null and forbidden raw-G/k routes | PASS |
| BA6 positive/negative reaction-free two-channel transient | PASS |
| BA7 P0/P1/P2, exact relativistic R10, NR limit, nonzero-eta R18 | PASS |
| BA8 open-system/beta ledger, R-a/R-b/R-c and unit boundaries | PASS |
| BA9 forbidden heat, unit, fate and double-count mutations | PASS |
| BA10a governed spin-on zero-BNV bit identity | PASS |
| BA10b spin-off matched control and transferred qualification | PASS |

Final source/projection log SHA-256:
`09ff8da22b250b2c384c68e1cef203ca09d9561f00550144b2d6d12a58b6d896`.
Final matched-control log SHA-256:
`52b718560edd236c3ecf51fe92ef47bb8286c124e0e4a12dd452ab372636d273`.

## M1 through M21

Every required pretrajectory mutant fired its accepted primary detector:

| Mutant | Detector result |
|---|---|
| M1 raw `G_y S_y` | DETECTED by sliding-null/source oracle |
| M2 diagnostic `k` substituted for `t` | DETECTED by tangent/sliding oracle |
| M3 wrong sign on `t Bdot` | DETECTED by baryon neutrality |
| M4 omitted `t Bdot` | DETECTED by sliding null |
| M5 omitted sigma channel | DETECTED by exact lift |
| M6 transposed/cross-Z fault | DETECTED by asymmetric two-channel oracle |
| M7 eta channel swap | DETECTED by named asymmetric channel oracle |
| M8 Fermi-hole double count | DETECTED by R-a/R-b/R-c oracle |
| M9 `Echem_dot` inserted as heat | DETECTED by reservoir-only ledger |
| M10 PdV/gravity inserted as heat | DETECTED by production-input/ledger oracle |
| M11 omitted `DeltaLnu` | DETECTED by beta ledger |
| M12 doubled equilibrium neutrinos | DETECTED by zero-source/beta ledger |
| M13 omitted MeV-to-erg | DETECTED at typed direct boundary |
| M14 doubled MeV-to-erg | DETECTED at typed direct boundary |
| M15 fluid/star escape confusion | DETECTED by fate-energy closure |
| M16 terminal fate double booking | DETECTED by unique fate/event identity |
| M17 stale tangent | DETECTED by currentness refusal |
| M18 source/tangent domain mismatch | DETECTED before projection |
| M19 `Bdot != b^T S` | DETECTED by atomic source closure |
| M20 ignored frozen failure | DETECTED by pre-RHS runtime refusal |
| M21 equilibrium neutron potential at nonzero eta | DETECTED by R18/O13 |

No R20 discrimination credit is claimed for M21.

## Protected regression and provenance disposition

- governed baselines: 11/11 entry hashes unchanged;
- protected Phase-5D paths: 33/33 entry hashes unchanged;
- Phase-5 authority/candidate artifacts: 14/14 unchanged;
- EOS/data inputs: 9/9 unchanged;
- authenticated literature manifest: 22/22 entries pass;
- Phase-5B structural response regression: PASS;
- Phase-5C chemical coefficient regression: PASS;
- focused direct and thermal ledgers: PASS;
- protected-manifest test: PASS;
- physical BNV rate selected: NO;
- physical BNV model selected: NO.

Focused CTest log SHA-256:
`f2218306f68e7ddd72484c679ea58f2703706bc4ccbdb1b3ebe2ad070727ebff`.
Hash-check log SHA-256:
`69f510db78442acd01b1265735925d4877cbd41c36c57f44e78a179baac134ec`.

## Authorization boundary after this commit

Only the four already committed controlled mathematical run cards may now be
executed. A run-card value, tolerance, output grid, certificate window or
interpretation may not be changed in response to its result. Any coupled-run,
ODE-refinement, R20, QSS, beta-bound or runtime frozen-validity failure stops
candidate acceptance; it does not authorize retuning. No governed baseline may
be created, and no result may be labelled physical or owner-ratified.
