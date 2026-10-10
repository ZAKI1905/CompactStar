# Campaign-readiness evidence

Status: incomplete nonlinear qualification. No full campaign was run.
Accepted implementation: 0c3f9269c884161c9d89163660a47b141a93ace6.
Predeclaration commit: d99ebb4. Production implementation bytes remain unchanged.

- `environment.json`: read-only source, baseline, input, package, compiler and
  standalone-test audit; unavailable actual hashes are null, never empty hashes.
- `recovery-*.json`: bounded local searches, including rejected/nonmatching copies.
- `recovered-inputs.json`: four copied and rehashed inputs, outside Git; no links
  were redirected. Frozen BNV certificate and coefficients are not recovered.
- `configure.log`: actual failed qualified-stack configuration. The required 13
  tests were not run; the historical 13/13 log remains in the accepted evidence.
- `linear-stability/`: 72 frozen-Jacobian diagnostic cases and high-precision
  reference assessment. This is not the nonlinear BNV model or its energy ledger.
- `linear-stability-dev01/`: retained failed input serialization attempt, before
  any numerical output. Fixed by converting a NumPy scalar to a Python float.
- `edge-probe.tsv`: 14 GSL-only cases, preserving accepted production behavior.
- `retained-summary.json`: reanalysis of authenticated archived windows; includes
  explicit BASELINE ledger and REFINED/ULTRA schedule coverage limitations.
- `verification.json`: fresh build/test commands, source/executable identities,
  schema checks and exact accepted-source comparison.
- `SHA256SUMS`: complete file hashes at publication (excluding this hash list).

The first provenance script invocation used the original Phase-5D generator and
correctly refused ChemicalResponse.cpp's old pin. The final audit records all 33
old pins alongside accepted-candidate hashes, including its three pre-existing
approved deltas; it never updates the old pins. The initial summary invocation
refused a zero division at stored endpoints. The final summary requires U=0 when
F=0 and excludes those exact endpoints from reconstruction-budget division;
nonzero uncertainty with zero allocation still fails. Neither repair changes a
scientific budget or any accepted implementation.

## Reproduction

From the readiness worktree, build the three dependency-independent diagnostics
with Apple clang and the existing MacPorts GSL. Exact invocations are in
`verification.json`; no optimization or fast-math flag is used. Run:

```
python3 -B tests/bnv/readiness_environment.py --output docs/validation/evidence/phase6-campaign-readiness
/Users/keeper/miniforge3/bin/python -B tests/bnv/readiness_stability.py --output /tmp/phase6-readiness-fresh-linear --probe build/readiness-standalone/linear_probe
/Users/keeper/miniforge3/bin/python -B tests/bnv/readiness_summary.py
```

The stability output directory must be fresh. Use the generated assessment and
TSV there for independent comparison; do not overwrite retained evidence. The
summary intentionally consumes the committed linear assessment and archived
nonlinear evidence. Python identities are recorded; the default framework Python
lacked mpmath, so the existing miniforge environment was used. No package installed.
The standalone diagnostics require no Zaki/CONFIND or frozen BNV inputs and do not
claim those missing dependencies were validated.

After Scriptorium recovery, authenticate with `cmake/verify_dependencies.py` and
all six input hashes before configuring. The exact regression selector is:

```
ctest --test-dir build/readiness-qualified --output-on-failure -R '^(chemical_production_contract|chemical_exact_oracles|chemical_curved_gc9|particle_number_analytic|passive_checkpoint_output_contract|phase6a1_direct_energy|phase6a1_thermal_ledger|heat_capacity_v1|cache_contract|cache_thermal_contract|phase5d_component_tolerances|phase6_breakpoint_contract|phase6_currentness_bytes)$'
```

Build those targets in the qualified Debug configuration first. A subsequent
bounded nonlinear test harness must implement the committed predeclaration before
it runs; the present linear diagnostic is not a substitute. No command here
starts the full campaign.

A final relink produced a different Mach-O executable hash. Its cause was not
assumed to be numerical change. The first successful run remains in
`linear-stability-build1/`; `linear-stability/` repeats the full experiment with
the final executable, and the validator requires every numerical/count field
(except wall time) to match. It records each executable's actual identity.

Diagnostic source files were uncommitted during the runs and were subsequently
committed unchanged as `dfc0bac`. The recorded source hashes, compile commands,
executable hashes and exact accepted-production comparison bind the executions;
the earlier `head` field alone does not identify these new diagnostic files.
