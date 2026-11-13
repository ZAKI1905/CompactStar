# First BNV–rotochemical competition study

**Completed S1 static map and S2 short simultaneous-driving tests, 2026-11-12.**
New scientific results are independently cross-checked candidates awaiting
scientific review. The owner approved restricted S1/S2 and separately ratified
Stage B; no longer-horizon or microscopic calculation is included.

Start with the [physics report](REPORT.md). At P=5 ms, Pdot=1e-20 on Structure-1,
a neutron-loss rate around 1e-12 yr^-1 cancels much of the rotochemical drive.
Surface emission is suppressed before rising again after drive reversal.
Near electron cancellation, about 0.413 keV deposited per destroyed baryon
restores the RC-only temperature. The formal stationary minimum is outside
the existing temperature and depletion limits.

| Deliverable | Location |
|---|---|
| Restricted Stage B owner ratification | [CONDITIONAL_ACCEPTANCE.md](../../validation/results-first-bnv-pilot/stage-b/CONDITIONAL_ACCEPTANCE.md) |
| Approved scope and equations | [ADR-0025](../../adr/ADR-0025-bnv-retained-energy-study.md) |
| Physics, assumptions, interpretation and next calculation | [REPORT.md](REPORT.md) |
| Execution, validation and retained failures | [EXECUTION.md](EXECUTION.md) |
| Standalone solver, checks and plotting | [tools/research/bnv_competition](../../../tools/research/bnv_competition) |
| Input values and original provenance | [inputs/manifest.json](inputs/manifest.json), [structure.json](inputs/structure.json) |
| All 5,580 combined stationary states, compact | [map.csv](data/map.csv) |
| Matching 5,580 BNV-only stationary states | [bnv-only.csv](data/bnv-only.csv) |
| 30 RC controls and zero-drive boundary | [controls.csv](data/controls.csv) |
| Full fiducial comparison and interactions | [fiducial.csv](data/fiducial.csv) |
| Channel-null detail and timing examples | [cancellation-zoom.csv](data/cancellation-zoom.csv), [representative-periods.csv](data/representative-periods.csv) |
| Four short production trajectories | [accepted.csv](evidence/short/accepted.csv), [endpoints.csv](evidence/short/endpoints.csv) |
| Independent static / evolution checks | [verification.json](evidence/verification.json), [independent-check.json](evidence/short/independent-check.json) |
| Large external data and hashes | [external-receipt.json](evidence/external-receipt.json), [external-SHA256SUMS](evidence/external-SHA256SUMS) |

Six candidate scientific figures have vector PDF and 200 dpi PNG exports:
[A: competition](figures/A_competition.pdf),
[B: surface temperature](figures/B_temperature.pdf),
[C: reaction powers and imbalance signs](figures/C_channels.pdf),
[D: applicability](figures/D_applicability.pdf),
[E: reachability](figures/E_reachability.pdf), and
[F: scaling and deposition](figures/F_scaling_deposition.pdf).
Each identifies the nonphysical benchmark and candidate status. The report
explains unresolved narrow cancellation cores and saturated plotting ranges.

## Reproduction on EKU

The static solver reads only committed inputs. It needs Python 3.11, NumPy
1.26.4, SciPy 1.14.1, Matplotlib 3.9.4 and mpmath 1.3.0; pinned requirements
are [alongside the scripts](../../../tools/research/bnv_competition/requirements.txt).
The supplied Slurm scripts use the existing cluster Python environment with
isolated SciPy/mpmath dependencies under the study r1 directory. They do not
modify the old environment. Update those paths to an equivalent environment
when reproducing elsewhere on the cluster.

From this worktree, use a **fresh** output directory under /mnt/sdd. For example:

    export COMP_SOURCE="$PWD"
    export COMP_OUTPUT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/bnv-competition-reproduction-new
    export COMP_COMMIT="$(git rev-parse HEAD)"
    mkdir -p "$COMP_OUTPUT"
    sbatch --partition=phyall --nodelist=phyn004 --cpus-per-task=12 --mem=8G \
      --export=ALL --output="$COMP_OUTPUT/static-%j.log" \
      tools/research/bnv_competition/static.sbatch

After that job succeeds, submit check.sbatch and supplement.sbatch with one CPU
and 4 GB each, using the same exports and separate output logs. They check
stationarity, root/stability independence and the cancellation region, and
render the figures. They may be ordered with Slurm afterok dependencies.
The example resources match the measured workload; no wall-time cap is imposed.
Choose a currently healthy node if cluster health changes. phyn001/phyn002 were
avoided because of previously observed NFS stalls.

S2 additionally needs the retained complete stellar archive and the qualified
toolchain/input paths specified in the existing source snapshot. The native
runner authenticates that archive and does not rebuild the sequence:

    sbatch --partition=phyall --nodelist=phyn004 --cpus-per-task=16 --mem=32G \
      --export=ALL --output="$COMP_OUTPUT/build-%j.log" \
      tools/research/bnv_competition/build.sbatch

After the build/contracts pass, submit short.sbatch with one CPU and 12 GB,
then short_check.sbatch with one CPU and 4 GB. Keep COMP_COMMIT set; compute
nodes lack git in their default PATH. The first script loads one retained
owner graph, constructs the matched contexts, and runs all four one-year
cases at two tolerances. The second independently integrates the declared
system and its thermal/chemical energy integrals. Neither invokes Stage B
or S3. The four-state native runner is deliberately not a fixture-free CTest.

package.py assembles the review bundle for this particular dated study and
writes its external manifest. It is for maintainers, not necessary to reproduce
the mathematics; it expects retained r1/r2/r3 siblings. Do not point it at a new
reproduction directory without adapting its evidence-selection policy.

## CSV interpretation

All values use infinity time/energy/temperature. Rates Gamma are per Julian
year; D and R are counts/s; eta and q_dep are MeV; powers are erg/s; times are
years. Columns ending keV report energy per destroyed baryon at infinity.
map.csv and bnv-only.csv align row by row; use the map's d and Gamma to identify
each four-corner comparison, matching RC controls by d. The common 00 boundary
has all thermal/power observables zero. Committed CSV line endings are normalized to LF without changing cells.
Full unprojected stationary-full.csv
and bnv-only.csv remain outside Git and are hashed.

- beta is undefined for spin-off cases (NaN); undefined per-event quantities
  at D=0 are also NaN. Infinite times at the undriven boundary are deliberate.
- difference_X = X_both-X_RC. interaction_X = X_both-X_RC-X_BNV+X_neither.
  Normalized interactions divide by abs(X_RC)+abs(X_BNV).
- stable means local frozen-model stability, not physical reachability.
- temperature_in_certificate tests only 1e6<=T_inf<=1e9 K.
  The legacy column name reachable_in_certificate means **only** temperature
  plus Gamma max(t_build,3 tau_slow)<=1e-6. It is not certificate admission.
- build_time_yr is a constant-drive energy-norm lower bound from eta=0.
  reach_time_estimate_yr is a heuristic local estimate, not a trajectory.
- rotation_e/mu, composition_eta_e/mu and composition_B_e/mu are distinct
  fractional-composition diagnostics, not new certified thresholds.
- thermal_residual, chemical_residual and subsystem_residual use the scales
  specified explicitly in model.py. They measure numerical equations, not
  uncertainty in stellar inputs or microscopic physics.
- The grid prescribes beta exactly at mathematical nulls. Gamma is its rounded
  display coordinate; recomputing cancellation by subtracting large double
  precision drives can leave a spurious tiny residual.
- k_B T and the illustrative smooth finite-T scale are diagnostics, not a
  demonstrated physical deposit or rigorous omission bound.

The reproducible audit.py script checks input/source hashes, untouched C1–C7,
local report links and all external-file hashes. Run it with Python 3.11 after
packaging; [packaging-check.json](evidence/packaging-check.json) records the pass.
[Visual inspection](evidence/visual-review.json) is recorded separately from
numerical verification.

Stage A/B evidence directories and all nine worktrees are preserved.
