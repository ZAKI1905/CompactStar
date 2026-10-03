# Phase-6 breakpoint candidate evidence

Status: **BOUNDEDLY QUALIFIED CANDIDATE / NOT RATIFIED / CAMPAIGN STOPPED**.

`assessment.json` is the separately proposed successor assessment, with original
BA12R state/ledger/rejection verdicts retained alongside it. `supplement-assessment.json`
contains the complete valid run4 comparison and the equilibrium-tail stability probe.
`manifest.json` authenticates every compressed run/check payload, its uncompressed
bytes and its permitted use. Each run's execution record identifies the executable,
source hashes, exact command, immutable inputs, return code and publication status.
No executable or build tree is included.

The final physical implementation is fa01369654f93a3ae473e6f51fe2a480d0bdfc26.
The later crossing-only harness adds explicit reference uncertainty and roundoff
reporting, without changing integration. `final07` is the fresh full bounded run on
the final implementation; `crossings08` supplies the stricter 37-crossing assessment.
`supplement01` uses the identical numerical main/reference paths for complete run4
comparisons and the tail probe. The exact-zero diagnostic constraint introduced
later does not affect those integrations or references.

Failed and interrupted development runs are preserved and excluded as indicated
in the manifest. In particular, run04 failed the one-ULP event test; run05 failed
zero-control energy-uncertainty propagation. The narrow repairs and their positive
and negative tests are documented in the candidate report and amendment. These
failures have not been relabeled as passes.

The original independent review is separately preserved, byte-unchanged, in
`../obs32-independent-review/`. Its validity index controls historical use: never
use run2 split M-* results, run3's corrupted control observation file, or metadata's
incorrect control-prediction sign. Run4's complete observations were compared at
all 490 source and 1277 control times. No historical rows were repaired.

Reproduction (from the worktree, with the authenticated input paths available):

```sh
cmake --build build/obs32 --target phase6_breakpoint_qualification phase6_breakpoint_contract passive_checkpoint_output_contract
python3 tests/bnv/run_breakpoint_qualification.py --root /absolute/fresh/root --mode all
python3 tests/bnv/assess_breakpoint_qualification.py /absolute/fresh/root/output
python3 tests/bnv/run_breakpoint_qualification.py --root /absolute/fresh/supplement --mode supplement
python3 tests/bnv/assess_breakpoint_supplement.py /absolute/fresh/supplement/output
```

The checked runner requires fresh roots and available space, authenticates all
six immutable inputs and the imported bundle, snapshots the executable, checks
artifact schemas/completeness/time order, and publishes the output directory only
after successful writes and read-back. It does not equate calculation completion
with numerical qualification. Individual payloads are deterministic gzip streams;
unpack them into scratch space to inspect the TSV/JSON records.
