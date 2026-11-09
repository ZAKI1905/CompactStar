# Native paired assembly: retained-data stop

**Assembly is blocked; no executable source/control snapshots were produced.**
The native runner is **not ready** for nonlinear actual-path qualification.
This delivery completes input/coverage investigation and a bounded run proposal,
not the requested paired assembly. `evidence/input-inventory.json:1` records null
source/control snapshots and zero stellar solves, actual BNV RHS evaluations and
thermal trajectories. No historical builder or scientific guard was changed.

## Authenticated starting point and scope

Repository: `/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot`.
Branch: `candidate/results-first-bnv-pilot-20261108`.
Starting HEAD: `afeb0418ab697829f3fec5722022f70db65f4b5b`, clean and explicitly
selected by the owner. All 14 other local branch tips are ancestors; nine
worktrees exist (`evidence/source-preflight.json:1`). Master was not changed.

Change classes: engineering evidence collection, documentation and compact
generated evidence. Authority: GOVERNANCE §§1–5; ADR-0015 §§7–9,14–16;
ADR-0016 §§1.1,3–5; ADR-0017 Amendment 1; owner instructions for this assembly
and subsequent archive/reconstruction investigation. Proposed ADR-0023 records
the ownership stop. No scientific/numerical source, baseline, coefficient,
certificate or historical failure was modified. Candidate status is unpromoted.

The conditional BA13 result remains 420 conditional, 21 exact and 21
metric-consumer-only rows, maximum utilization 0.782464607730239. It was read and
authenticated, not reassessed (`evidence/input-inventory.json:1`). This remains
the restricted P2 mathematical benchmark; no microscopic rate was introduced.

## What the retained material supplies

All 14 checked input files match their recorded SHA-256 identities: the eleven
original pinned files plus native tangent, target-star and thermal-surface tables
from the original fitter binding (`evidence/input-inventory.json:1`). The central
profile SHA is `62b8dc594fd9c30d413ccb1d906459e2599049966b0bec097284b43598a07a56`;
the chemical certificate is `337ed8ae8d6168d8a612f1344f204dc927ac9e0867a47542d5b943dbf27b8537`;
the coefficient table is `ec3942227c98532af4802db4bdd229603715556af9d44929082368b3afc7d831`.
The complete thermal/model/EOS/frozen-certificate hashes are in that inventory.
These authenticate inputs; they do not certify reconstruction of an owning star.

The common proposed initial state is t=0, Tinf=Tref=1e8 K, eta=(0,0) MeV,
spin=0. Native B0 is 7.616906518718875e56 and derived Bdot is
-2.413652026364133e37 count/s. The historical Mac literal is
-2.4136520263641375e37, produced from its different B0; the exact binary64 and
80-digit reconciliation remains `../preparation-evidence/bdot-witness.json:1`.
No literal or table was substituted. Source/control configuration remains as
predeclared in `../EXPERIMENT.md:11`; there are **no snapshot SHA-256s to report**.

Source inspection confirms the intended common owners: MU-only normalizations
SMe=1e-51 and SMmu=2e-51 erg cm^-3 s^-1 K^-8, coupled Z in MeV/count, thermal
Cstar in erg/K, eta in MeV and x=ln(Tinf/Tref). The photon owner selects the
default iron envelope at rho_b=1e10 g/cm^3, xi=0, full emitting surface. Its
`EnvelopePotekhin2003_Iron` implementation delegates to the existing 1997 fit:
Ts^4=1e24 g14 (1.81 Tb/1e8 K)^2.42; Tsurface_inf=Ts exp(nu_surface),
Lgamma_inf=4 pi (R_km*1e5)^2 sigma Ts^4 exp(2 nu_surface).
This is source inspection, not a new runtime verification
(`tests/rotochemical/coupled.hpp:13`;
`CompactStar/Physics/Rotochemical/FrozenRotochemicalRunContext.hpp:104`;
`CompactStar/Physics/Driver/Thermal/PhotonCooling.hpp:247`;
`CompactStar/Physics/Driver/Thermal/Boundary/src/EnvelopePotekhin2003.cpp:29`;
`CompactStar/Physics/Driver/Thermal/Boundary/src/EnvelopePotekhin1997.cpp:89`;
`CompactStar/Physics/Driver/Thermal/src/PhotonCooling_Details.cpp:311`).

## Completed coverage investigation

The external evidence directory is
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/results-first-native-assembly-20261109-r1`.
Its detailed file/member lists and source/interior analyses are pinned by
`evidence/external-manifest.json:1`. Small receipts are also committed here.

* The 15 target-0 sequence scratch directories in both run1 and repeat have no
  files. The prior native-current and native-verify central/sequence directories
  also have no files (`evidence/input-inventory.json:1`).
* The two native archives contain 249 and 284 regular files, respectively;
  their member lists contain chemical profiles/tables/certificates, not complete
  owner snapshots (`evidence/archive-coverage.json:1`). The broader supplied
  handoff/confirmation archives were also inspected: 15 archives were readable;
  the remaining historical-source.tar is a zero-byte prior-attempt artifact.
  The full scan receipt is in the external manifest.
* A complete filename scan included 51,620 cluster files, including ignored
  build scratch and retained worktrees. There are 319 large chemical-profile
  files with the expected schema and 105 distinct central-energy ratios.
  Their range is [0.999993822756231, 1.0000000417050623] relative to native B0's
  central energy. None supplies the eight nonzero stencil offsets
  ln(rho/rho0)=+/-0.00025, +/-0.0005, +/-0.001, +/-0.002.
  See external `interior-coverage.json:1` and `evidence/reconstruction-summary.json:1`.
* A Linux Debug Phase-5B artifact **was found** in toolchain build scratch. It
  has derivative/stencil/tail metadata, but no radial owner arrays. It uses EOS
  SHA `6505db3821cc83c29b6a4301809cb6b1bcbe439e57d5e7486f8901540620e64a`, whereas
  qualified native EOS is `f6b9c4068a430cc356b54e39c21a0ca17e17b8c2b35d17b836638298d1112d7f`.
  Its provider and table-generator source hashes differ, and B0 is
  7.616906518718873e56. It is ineligible for substitution
  (`evidence/phase5b-candidate-review.json:1`).
* `independent-review-20261011/experiments/hartle-captures/` contains four-row
  DS-CMF-1 surface regression summaries at radial10000, not free-gas radial80000
  interiors. The producer's provenance claims do not serialize those interiors
  (`tests/analysis/produce_particle_number_reference.py:79`).
* Git history enumeration completed: 5,745 distinct historical paths. The
  relevant stored artifacts are the Phase-5B summaries above and old DS-CMF
  `NStar_Profile` outputs, whose producers identify a different EOS
  (`main/Test/ns_build_main.cpp:63`, `main/Test/spin_therm_evol_2_main.cpp:116`).
  The local partial repository has 4,916 omitted objects; a content scan hit
  the unavailable `OMITTED_OBJECTS_NOT_SUPPLIED` promisor. A no-lazy-fetch tree/
  historical-path scan completed successfully, but **omitted historical content
  was not claimed inspected** (`evidence/reconstruction-summary.json:1`).

This establishes insufficient *available, authenticated* reconstruction material,
not a claim that no copy can exist elsewhere. No remote branch was fetched,
changed or published. No Mac Scriptorium mount is visible from this session.

## Reconstruction feasibility and smallest regeneration

The retained central profile supplies exact r, m, nu, nB, pressure, energy and
rounded nodal species densities. Its thermal geometry, chemical quadrature and
MU coefficients can in principle be rebuilt without a global stellar solve and
checked against the retained rows. That alone does not construct the existing
full context (`tests/bnv/frozen_sensitivity.cpp:134`).
The `nn`/`np`/`ne`/`nmu` columns are local provider evaluations; the distinct
`*_profile` columns retain the actual profile density products. Substituting the
former for the latter would change composition (`tests/bnv/frozen_sensitivity.cpp:111`).

Missing central data include original binary64 species fractions, nu-prime and
EOS d(eps)/dp; original first-order/monopole arrays are also absent. The nodal
product cannot uniquely restore the species fraction: at row 0 two adjacent
binary64 neutron fractions yield the same retained n_n. An explicit arithmetic
witness is `evidence/composition-reconstruction.json:1`. Re-evaluating EOS fields
and integrating Hartle responses on a loaded fixed background is conceivable,
but needs a new loader and validation of every rebuilt field, source ownership,
and bit/tolerance agreement; it is not demonstrated exact reconstruction here.
The ordinary constructor rebuilds nu and first-order rotation, while the monopole
requires EOS derivative data (`CompactStar/Core/src/NStar.cpp:303`;
`CompactStar/Core/src/RotationSolver.cpp:1381`).

More decisively, the **eight distinct noncentral stellar interiors are absent**.
Neither the narrow one-sided BA13 family nor the aggregate t row can reconstruct
the accepted two-sided three/five-point derivative and its uncertainty.
Integrating a homogeneous sensitivity equation on the central interior would be
a different structural method, explicitly an independent oracle rather than the
accepted production owner (`../../PHASE5B0_INV09_GLOBAL_RESPONSE_PREFLIGHT.md:350`).
No no-global-solve reconstruction preserving the current guards has been
established. A numerical t or sidecar cannot be attached in place of that owner.

If no complete matching archive is recovered, the smallest straightforward
regeneration using the **unchanged accepted sequence producer** is one
15-contributor sequence at rho0=1.10e15, h=(.001,.0005,.00025), k=-2..2,
radial80000, using the existing authenticated EOS bytes. That is **15 canonical
stellar solves**, not a 21-star BA13 rerun. Select its first k=0 contributor
(ordinal 2) as the nominal central star, after exact profile/B0 equality checks;
reuse that one sequence for both the fixed-B response and the tangent. Computing
the central A response on this owned star requires no additional TOV solve.
This assembly change must retain central-first chemical lifetime ownership and
every sequence contributor. The recipe
contains only nine distinct densities (eight noncentral); reducing the number of
solves through reuse would itself require an explicitly validated loader/replay
change, so eight is a lower bound, not a ready executable plan
(`tests/rotochemical/fixture.hpp:44`; `tests/bnv/campaign_fixture.hpp:19`;
`CompactStar/Analysis/src/ParticleNumberResponse.cpp:369`).

This calculation is **NOT AUTHORIZED AND WAS NOT RUN**. Before requesting it,
implement a complete capture path that exports all profile columns, EOS derivative,
Hartle arrays, achieved states, tails, A/B/K values/errors and stencil provenance;
retain a manifest and validate a fresh-process reconstruction. Require exact
native central-profile/B0 and retained Z/Ltilde/t/error agreement with the existing
accepted budgets; mismatch stops, with no pin or budget update. This would be a
bounded data-regeneration request, not a new stellar-solver investigation or a
reopening of BA13. It cannot proceed under the present no-new-solves instruction.

## Verification and remaining execution gates

Six authenticated retained Linux executables pass: passive checkpoint contract,
direct energy, thermal ledger, component tolerances, breakpoint contract and
byte-currentness. Their SHA-256s match the prior build receipts; dynamic-library
resolution reports no missing bindings. No C++ was changed or rebuilt. These are
synthetic/analytic tests, **not native-pair tests**
(`evidence/retained-contract-results.json:1`; external per-test logs/ldd records).

The inventory exits 2 with `BLOCKED_EMPTY_KNOWN_OWNER_DIRECTORIES`; this is a
recorded assembly stop, not a failed BA13 test. No pair was loaded, no native-pair
initial RHS/R18 was evaluated, and no pair serialization, isolation or restart
result is claimed. Finite-interval R20, nonlinear stability, accumulated state
error, event/reconstruction error, rejection fraction and resolved matched
temperature signal all still require actual trajectories
(`NONLINEAR_PLAN.md:1`). Passive replay is not full solver-state restart.

No Slurm job was submitted in this task. Proposed build ceiling is 4 CPUs/8 GiB/
15 minutes; each future trajectory stage is 1 CPU/4 GiB/30 minutes, single-thread
libraries on phyall/phyn004. These are resource ceilings, not runtime forecasts.
The exact proposed 1-year and 2048-year run matrices, norms and stops are in
`NONLINEAR_PLAN.md:1`. No 14,648-year or 500 kyr trajectory is launched or requested
while snapshots are absent. The next necessary decision is recovery of complete
owner data or a separately reviewed bounded regeneration, not another static assay.
