# Phase-5D-1 structural/background resolution qualification plan

**Status:** PRE-RESULT / IMMUTABLE QUALIFICATION PROTOCOL. This record is committed before any new Phase-5D resolution-qualification output is generated.

## 1. Entry and blocker

- Canonical entry: `d019ae390be4f5e3daba05039903485cb497e397`.
- Implementation branch: `physics/phase5d-controlled-rotochemical-evolution`.
- Branch blocker SHA before qualification: `30a15b299ddd7eedfb2f5ab824179b48a20ca3ea`.
- Preserved original predeclaration SHA: `697ec14610d63640866f3c0aed0a51d567cb037d`.
- Preserved response implementation SHA: `14b74cba1a0c7b6faec591aedf6f5c9a76fa42d6`.
- Worktree: `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase5d-controlled-evolution`.

The retained original declaration is
`docs/validation/PHASE5D1_CONTROLLED_EVOLUTION_PREDECLARATION.md`; the blocker and
full final report are retained in
`docs/validation/PHASE5D1_CONTROLLED_ROTOCHEMICAL_EVOLUTION_IMPLEMENTATION.md` and
`docs/validation/PHASE5D1_FINAL_REPORT.md`. The repository copies of the exact
refusal are `docs/validation/phase5d1_w_assembly_refusal.json` and
`docs/validation/phase5d1_refusal_excerpt.txt`; authenticated scratch remains
under `build/phase5d-audit/`.

At the originally predeclared radial resolution 10000, fresh structural W
assembly refused before the ODE. The positive structural contributions
`sum_j |Z_ij| E_I_j` were `8.918515424724095e-12 MeV s^2` for Npe and
`3.263850588763438e-11 MeV s^2` for NpMu. They alone exceeded the unchanged
Phase-5C numerical goals `1e-12 MeV s^2` and `4e-12 MeV s^2` by factors
`8.918515424724095` and `8.159626471908595`. The attempted 10000-point
certificate also reused background/provider/anchor characterization from
40000/80000 work while recomputing only part of the new-resolution evidence.
That mixed-resolution certificate is scientifically inadmissible.

This blocker does not show that ADR-0014 is inconsistent, that the W formula is
wrong, that governed Phase-5C regressed, that RKF45 is unsuitable, or that the
controlled benchmark is impossible. It shows only that (a) the predeclared
10000-point structural construction cannot satisfy the unchanged W numerical
budget and (b) the attempted 10000-point background certificate was not
independently qualified at that resolution. The original predeclaration
correctly stopped before any result.

## 2. Candidate fixed before output

The sole replacement production candidate is fixed now at
`radial_resolution = 80000`. This is the governed Phase-5B/Phase-5C
structural-fixture resolution underlying the accepted coefficient machinery;
it is not selected by inspecting a Phase-5D trajectory or new qualification
output. The EOS resolution remains 8192 intervals and the central mass-energy
density remains `rho_c = 1.10e15 g cm^-3`. The physical fixture remains the
ordinary NStar Track-R Structure-1 midpoint, whole star, with the existing
finite-cut/tail contract and unchanged source/provider authority.

There is no post-result resolution selection. If 80000 fails any gate, this
qualification stops. It will not try 120000, 160000, or any other radial
resolution, and it will not select 20000 or 40000 as a Phase-5D production
resolution after observing output.

## 3. Governing numerical/error method

The calculation uses the accepted Phase-5C construction without redefining an
error. Controlling records are ADR-0013,
`PHASE5C1_NUMERICAL_BUDGET_AND_IMPLEMENTATION_PLAN.md`,
`PHASE5C1R_STRUCTURAL_UNCERTAINTY_SEMANTICS_RATIFICATION.md`,
`PHASE5C2_PRODUCTION_ACCEPTANCE_PREDECLARATION.md`,
`PHASE5C2_CORRECTED_CHEMICAL_COEFFICIENT_IMPLEMENTATION.md`, and
`PHASE5C2_CORRECTED_CHEMICAL_COEFFICIENT_INTEGRATION.md`, together with the
governed Phase-5B and Phase-5C producers/regressions.

`numerical_error`, `validation_envelope`, and any certified analytic bound are
distinct. The qualification will preserve the accepted meaning and accounting
of each. It will not reinterpret a numerical error as a certified bound or a
validation envelope, and will not invent a new W-error definition.

The accepted chemical-background characterization is regenerated in fresh
scratch with its governed ladder:

- EOS 4096 / radial 40000: lower EOS characterization only;
- EOS 8192 / radial 40000: lower radial/background characterization only;
- EOS 8192 / radial 80000: the sole production background and finest characterization.

The accepted Phase-5B structural radial characterization additionally uses
radial 20000 and 40000 evidence where the governed PB6/producer method calls for
the 20000/40000/80000 ladder. These are characterization runs only. The
production `A/B/K/I` construction and W factory use radial 80000.

## 4. Fresh dependency reconstruction

All resolution-dependent quantities required by the accepted production
certificate will be recomputed from governed source inputs in new, separate
scratch directories. No stale object, runtime baseline JSON, reviewed candidate
JSON, copied certificate, or copied `G`, `Q`, `Z`, `I`, or `W` value may seed the
calculation. In particular the qualification freshly reconstructs:

- the EOS table/profile/background for every governed characterization member;
- provider nodal comparisons and equilibrium-anchor characterization;
- onset/refusal windows, center contribution, surface/tail continuation and
  tail enclosure;
- the central and neighboring sequence stars, structural `A`, `B`, `K`,
  `K_error`, `I_phys`, and `E_I` quantities;
- global `G`, every accepted G error component, support/rank/conditioning and
  partition diagnostics;
- derived `Q`, `E_Q`, `Z`, `E_Z`, solve/arithmetic contributions;
- `W`, each numerical-error contribution, total `E_W`, and the separately
  classified validation envelope; and
- compiler, architecture, build configuration, EOS/source hashes, provider and
  profile identities/versions, chemical domain, quadrature, tail policy,
  node/partition counts, currentness and lifetime/refusal evidence.

Resolution-independent authorities may be reused only as immutable contracts,
not numerical seeds: physical constants and unit definitions, the governed EOS
source generator/model identity, source-code bytes, fixed absolute goals, the
predeclared empirical `V_I_validation` envelope, and the accepted compiler-only
portable-provenance comparison rule. Their source authority and classification
will be recorded field by field. No radial-dependent cached evidence is
authorized for reuse.

The final 80000 certificate must identify the production resolution, or the
resolution-independence authority, for every field. A radial-dependent
contribution characterized only at another production resolution is a failure.
There will be no copied certificate with a changed label and no mixed-resolution
certificate.

## 5. Execution isolation

Characterization members may run concurrently only when they have distinct
scratch/output directories and do not share mutable star/profile/cache state or
result files. They use the same compiled executables/toolchain and governed
EOS/source inputs, with deterministic resolution-labelled run identities.
Production assembly is a separate fresh process after its required
characterization completes. Git writes, staging, commits, and pushes remain
serial under exactly one repository writer.

No coupled-evolution executable or trajectory will run. The archived coupling
patch under `build/phase5d-audit/uncommitted-coupling/` will not be restored.

## 6. Frozen pass/fail criteria

The unchanged production factories must all accept:
`GlobalChemicalNumberResponse`, `ChemicalImbalanceResponse`, and
`RotochemicalSpinDrive`. No bypass or test-only construction of refused W is
permitted.

The complete accepted W numerical-error sum will be reported term by term,
including `|Z| E_I`, `E_Z |I|`, `E_Z E_I`, arithmetic error, and every other
governed nonnegative contribution. Required goals remain:

| Gate | Npe | NpMu |
|---|---:|---:|
| `E_W_numerical` | `<= 1e-12 MeV s^2` | `<= 4e-12 MeV s^2` |
| `V_W_validation` | `<= 1e-10 MeV s^2` | `<= 3e-10 MeV s^2` |

The unchanged absolute G/Q/Z goals from the Phase-5C predeclaration remain in
force, as do support, rank, conditioning, refusal, tail, currentness, and
provenance requirements. The governed result comparison is validation only and
occurs after fresh production; neither the Phase-5C baseline nor reviewed
candidate may seed production. Central values with unexplained drift from the
governed fixture fail closed.

Fresh `phase5b_structural_response_regression` and
`phase5c_chemical_coefficient_regression` must each pass with raw rc=0. The
governed baselines and reviewed candidate must remain byte-identical. The
previous focused response 2/2, data-free 47/47, and authenticated 70/70 evidence
is authenticated but not rerun unless executable/test logic changes.

Any factory refusal, W numerical or validation-envelope exceedance, G/Q/Z gate
failure, mixed-resolution need, unavailable fresh radial-dependent field,
unexplained central-value drift, upstream regression failure, need for another
radial resolution, or need to modify another benchmark parameter stops the
task. Goals and tolerances will not be relaxed.

## 7. Downstream boundary

Only if every fresh 80000 gate passes may a new addendum supersede the original
10000 radial choice. Such an addendum may change only the production radial
resolution to 80000. SMe, SMmu, enabled/disabled processes, T0, eta0, spin
history, run interval, RKF45 intent and component tolerances, refinement policy,
quasi-steady criteria, source-function roots, RE10b tolerance, and every other
benchmark-defining value remain frozen.

Even after a pass, this task will not enter an ODE, generate a secular
trajectory, restore coupling, create an evolution/candidate SHA, implement A18,
begin BNV, alter Z/W physics, change benchmark normalization, relax a Phase-5C
goal, or merge. Qualification failure leaves the original fail-closed history
and coupling block in force.
