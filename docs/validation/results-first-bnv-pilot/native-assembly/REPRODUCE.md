# Reproduce the assembly-input stop and retained contracts

These commands do not build a star, evaluate an actual BNV RHS or evolve a thermal
state. They reproduce the input inventory and six existing synthetic/analytic
contracts. They are **not** commands for native paired assembly, which is blocked
as documented in `REPORT.md:1`. Use a fresh cluster output directory.

```bash
export RF_SOURCE=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot
export RF_OUTPUT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/native-assembly-reproduction
mkdir "$RF_OUTPUT"
cd "$RF_SOURCE"
python3 -B tools/linux/results_first/inventory_native_assembly.py \
  --source "$RF_SOURCE" --output "$RF_OUTPUT/input-inventory.json"
# Expected exit 2, BLOCKED_EMPTY_KNOWN_OWNER_DIRECTORIES, input_hashes_match=True.
sha256sum -c docs/validation/results-first-bnv-pilot/endpoint-acceptance/SHA256SUMS
```

Authenticate the retained executables before running them. These are exactly the
prior prerequisite binaries, not a newly compiled native experiment. Run from the
fresh output directory so no retained test logs are overwritten.

```bash
python3 -B - <<'PY'
import hashlib, json, os, subprocess
from pathlib import Path
source = Path(os.environ['RF_SOURCE'])
output = Path(os.environ['RF_OUTPUT'])
pins = json.loads((source / 'docs/validation/results-first-bnv-pilot/endpoint-acceptance/evidence/prerequisites/artifact-hashes.json').read_text())
names = {'passive_checkpoint_output_contract', 'phase6a1_direct_energy',
         'phase6a1_thermal_ledger', 'phase5d_component_tolerances',
         'phase6_breakpoint_contract', 'phase6_currentness_bytes'}
for path, expected in pins.items():
    executable = Path(path)
    if executable.name not in names:
        continue
    if hashlib.sha256(executable.read_bytes()).hexdigest() != expected:
        raise RuntimeError('Executable authentication failed: ' + path)
    with (output / (executable.name + '.log')).open('x') as log:
        subprocess.run([path], cwd=output, stdout=log,
                       stderr=subprocess.STDOUT, timeout=30, check=True)
PY
```

The existing file and history inventory can be reproduced without lazy fetching:

```bash
rg --files --hidden --no-ignore /mnt/sdd/zaki/CompactStar > "$RF_OUTPUT/cluster-file-list.txt"
GIT_NO_LAZY_FETCH=1 git log --all --no-renames --format= --name-only \
  > "$RF_OUTPUT/git-historical-paths.txt"
GIT_NO_LAZY_FETCH=1 git rev-list --objects --all --missing=print \
  > "$RF_OUTPUT/git-object-paths.txt"
tar -tzf /mnt/sdd/zaki/CompactStar/linux-transition-20261010/ba13-recovery/native-static-inputs.tar.gz \
  > "$RF_OUTPUT/native-static-members.txt"
tar -tzf /mnt/sdd/zaki/CompactStar/linux-transition-20261010/ba13-recovery/native-ba13-candidate-20261012.tar.gz \
  > "$RF_OUTPUT/native-candidate-members.txt"
```

The expanded investigation receipts, complete archive/member lists, old Phase-5B
source mismatch, central-density coverage and composition-inverse witness are
retained externally with hashes in `evidence/external-manifest.json:1`. Their
compact results are committed alongside this file. The additional arithmetic
checks read retained data only; they are not new numerical qualifications.

The proposed future trajectory matrix is `NONLINEAR_PLAN.md:1`. Do not invoke the
historical native-binding fixture, Phase-5B producer, paired builder or campaign
to reproduce this report: those paths perform new stellar solves or trajectories.
