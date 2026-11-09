# First nonlinear P2 failure: chemical-energy uncertainty floor

**Reproduced: Echem is the only failing diagnostic.** At checkpoint 1,
**t = 131,490 s**, its estimated uncertainty exceeds the unchanged allocation by
**8.077630769832355 times**. O1 and O2 chemical energies are identical. The refusal
is driven by propagation of the prescribed chemical-state uncertainty floor,
not by measured O1/O2 disagreement. No tolerance, formula, solver or acceptance
criterion was changed (`evidence/checkpoint-1/diagnostics.tsv:1`,
`evidence/checkpoint-1/failure.txt:1`).

## Exact failure

| Quantity | Value | Units |
|---|---:|---|
| Echem, both O1 and O2 | 2.8395321373652751e28 | MeV |
| U | 4.5873384329419162e17 | MeV |
| F | 2.8395321373652749e17 | MeV |
| Allowed 0.20 F | 5.6790642747305504e16 | MeV |
| U/F | 1.6155261539664711 | dimensionless |
| Utilization U/(0.20 F) | **8.077630769832355** | dimensionless |
| Absolute O1/O2 disagreement d | **0** | MeV |
| 2 max(d, 64 ulp) | 5.62949953421312e14 | MeV |
| Propagated state contribution | **4.581708933407703e17** | MeV |

Echem is the **whole-star chemical free energy**, not temperature, chemical
imbalance, luminosity or a first-law residual. It is stored in MeV and converted
to erg at the R20 boundary. G = 3.6549335738821012e34 erg/s is the gross-power scale
for power diagnostics; Echem does **not** use that dimensionally different scale.
The diagnostic uses F = 1e-11 |Echem| here, hence its allowed relative uncertainty
is 2e-12. The measured estimate U/|Echem| is about 1.6155e-11
(`evidence/checkpoint-1/diagnostics.tsv:1`, repository
`CompactStar/Physics/BNV/src/PassiveCheckpointOutput.cpp:623`, `:117`).

All three reconstructed state components pass their state criteria. Their
utilizations are x: 0.094046, eta_e: 0.091510, eta_mu: 0.082131. Of 18 diagnostic
channels, 17 pass and only Echem fails. There is no new R18 or R20 failure inferred
from this result; the original R20 computation was never reached
(`evidence/checkpoint-1/state.tsv:1`, `evidence/checkpoint-1/diagnostics.tsv:1`).

## Precise estimator construction

For a diagnostic Q, the unchanged implementation defines

```text
d_Q = |Q(O2) - Q(O1)|
M_Q = max(|Q(O1)|, |Q(O2)|)
floor_Q = max(1e-11 max(1, M_Q), 64 ulp(M_Q))
F_Q = max(1e-11 G, floor_Q) for the enumerated power channels;
      floor_Q otherwise (including Echem).
P_Q = 2 max over eight corners |Q(O2_state +/- u_state) - Q(O2_state)|
U_Q = 2 max(d_Q, 64 ulp(M_Q)) + P_Q
Require U_Q <= 0.20 F_Q.
```

These are observed implementation definitions, not a newly proposed estimator.
The existing `Lout` branch is outside the enumerated gross-power cases; it is
zero and passes here, and was not changed or investigated further. The chemical
energy box uses all three state components, although frozen Echem depends only
on the two chemical imbalances (`PassiveCheckpointOutput.cpp:605` and `:623`).

The state radii passed to that box are

```text
D_O2,i = atol_O2,i + rtol_O2 max(|O1_i|, |O2_i|)
F_O,i  = max(D_O2,i, 64 ulp(M_O,i))
U_O,i  = 2 max(|O2_i-O1_i|, F_O,i)
u_state,i = U_O,i + U_event,i
```

The accepted owner-certified exact-zero-control exception is inapplicable to the
nonzero driven source. Here the primary state disagreements are all zero, and
F_O = D_O2, so the **requested tolerance floor** determines U_O. The fixed O2
rtol is 1e-13, chemical atols are 1e-24 MeV; O1 remains 1e-12 and 1e-23 MeV.
Independent local witnesses use rtol 1e-14, chemical atols 1e-25 MeV, with
half-span rk8pd and quarter-span RKF45 caps. No witness was retuned
(`evidence/checkpoint-1/checkpoint.json:1`, `evidence/assessment.json:1`, repository
`CompactStar/Physics/BNV/src/PassiveCheckpointOutput.cpp:373`, `:436`, `:519`).

| State component | U_O | U_event | Propagated radius |
|---|---:|---:|---:|
| x, dimensionless | 2.0617071852343281e-18 | 3.4940109074239889e-21 | 2.0652011961417521e-18 |
| eta_e, MeV | 2.0907862781724682e-24 | 3.2311742677852644e-27 | 2.0940174524402534e-24 |
| eta_mu, MeV | 2.2302905957642751e-24 | 1.2924697071141057e-26 | 2.2432152928354161e-24 |

## Why the energy estimate fails

The unchanged frozen diagnostic is

```text
Echem = [Z11 eta_e^2 - (Z01+Z10) eta_e eta_mu + Z00 eta_mu^2]
        / [2 (Z00 Z11 - Z01 Z10)]   MeV.
```

It is quadratic and tiny relative to its later-time value, while the chemical
absolute-tolerance floor remains fixed. Propagating that finite box into a strict
relative energy budget produces the refusal. The propagated box accounts for
**99.8773% of total U**. Linearized decomposition shows the absolute-atol floor
contributes about 4.345434117186622e17 MeV, or **7.65167 allocations by itself**;
the rtol floor contributes 2.27162570989222e16 MeV, or **0.40000 allocations**.
The atol contribution explains about 94.84% of the propagated box; local-witness/
event and floating-point effects account for the remaining small terms
(`evidence/assessment.json:1`, repository
`CompactStar/Physics/BNV/src/FrozenControlledBnvRunContext.cpp:161`).

For both chemical components, the recorded U_event is exactly its 64-ulp floor:
the independent local witness discrepancy term is zero. The x witness has a
nonzero discrepancy, but frozen Echem does not depend on x
(`evidence/assessment.json:1`).

This estimate must be distinguished from measured disagreement:

- O1/O2 chemical-energy disagreement is exactly zero.
- Maximum O2 disagreement with any of the five retained original references is
  **3.0786325577728e13 MeV**, only **0.0005421** of the allowed energy budget.
- The original main's entire accepted-state history reproduces byte-for-byte.
  Its enclosing step is ordinal 9, [85,212.386075900038, 417,541.69177191023] s.
  Full states, local RK error, suggested step and tolerances are persisted.

References are corroboration, not a certified accumulated-roundoff bound. This
replay does not establish the full Stage A accuracy gates
(`evidence/assessment.json:1`, `evidence/main.accepted_states.tsv:1`,
`evidence/checkpoint-1/checkpoint.json:1`).

An offline 100-decimal-digit evaluation of the **same algebra at the already
recorded binary64 states and coefficients**, with no new RHS or trajectory,
independently reproduces the energy box: 4.581672737754628e17 MeV versus recorded
4.581708933407703e17 MeV. Center algebra rounding is about 1.8812e12 MeV; maximum
rounding among recorded center/corner evaluations is 6.7319e12 MeV. The quadratic
numerator cancellation condition is 1.04052 and determinant condition 1.00114.
Thus catastrophic cancellation in this energy formula is not the dominant cause.
This calculation checks local algebra only; it does not bound accumulated ODE
roundoff (`evidence/assessment.json:1`, `tools/linux/native_qualification/assess_checkpoint_failure.py:1`).

## Classification and smallest next action

**E: an incompatibility between a conservative chemical-state tolerance floor and
the stricter relative chemical-energy allocation near initial equilibrium.**

- **A, genuine nonlinear integration error:** not demonstrated; observed main,
  O1/O2 and independent-reference agreement is much better than the failed budget.
- **B, reference-trajectory uncertainty:** not the term triggering this refusal;
  the diagnostic criterion does not use the five global references.
- **C, accumulated roundoff/cancellation:** a global bound remains outstanding,
  but recorded local rounding terms are far too small to explain this factor-eight
  failure; high-precision algebra corroborates the propagated-box mechanism.
- **D, diagnostic implementation/normalization defect:** no incorrect formula,
  unit conversion or acceptance normalization is established. The lost failure
  report was an observability defect and is corrected in the diagnostic collector.

No solver-algorithm repair is indicated by these measurements. The smallest
prospective numerical change is **a tighter O2 chemical reconstruction absolute
tolerance derived from the existing Echem budget**, leaving the main solver,
physical inputs, F and 0.20 criterion unchanged. This is a numerical-method change
requiring justification and separate authorization; it was **not implemented**.
Replacing the floor with zero merely because witnesses agree is unjustified.
Using the power scale G to normalize Echem would mix energy and power units.

The linearized necessary common O2 chemical-atol bound at this failed checkpoint
is about **7.84e-26 MeV**, before reserving witness/event/roundoff costs. Merely
fixing that point is insufficient: the first nonzero checkpoint on the existing
1920-interval grid is t=16,436.25 s, where the same algebra using the retained
reference requires a bound below approximately **9.80e-27 MeV**, versus current
1e-24 MeV. These are allocation estimates, not accepted replacement settings or
a guarantee of PASS (`evidence/assessment.json:1`).

Next bounded work should predeclare a chemical-atol choice with margin below the
finest-grid budget, retain this failure reporter, and validate the first nonzero
checkpoints and independent error/roundoff bounds before another matched Stage A
attempt. State accuracy, reference interpolation, R18/R20, restart, accumulated
roundoff and source/control criteria all remain required. This diagnosis alone
does not authorize or qualify a new 1-year pair, Stage B, or a full campaign.

## Execution, observability and preservation

Original snapshot identities and all 509 bound files are unchanged. The separate
collector extracts `native_main` verbatim from authenticated historical source,
links unchanged retained numerical archives, and uses the original source initial
state/algorithm/settings. It replays only one original-t0 source main, then stops
at the first failing passive checkpoint; the five retained references are read,
not integrated again. It evaluates no later checkpoints or control trajectory.
All original authentication and ownership/scientific guards remain active
(`PREDECLARATION.md:1`, `evidence/build-manifest.json:1`, `evidence/invocation.json:1`).

The reusable reporter writes all state and diagnostic channel budgets, not just
the first failed name, to TSV; writes time/index/bracket/tolerances to JSON; and
flushes a concise `failure.txt` before termination. Its delegating context records
the original local RHS evaluations and diagnostic corners without adding calls.
The two primary contexts had 15 RHS calls each; refinement and independent
contexts had 29 each. These are RHS-stage counts, not accepted-step counts.
The historical collector remains intact for original snapshot reproducibility;
future qualification collectors should use the new reporter
(`tests/diagnostics/checkpoint_evidence.hpp:1`, `tests/diagnostics/checkpoint_replay.cpp:1`,
`evidence/assessment.json:1`).

Regression tests deliberately produce a failed Echem diagnostic and failed state
components. They verify exact persisted U/F/threshold, units, checkpoint/time,
and that reporting changes neither numerical values nor refusal status.
Tests passed under Slurm job **19175**. Build/test: 27.145 s, 504,332 KiB maximum
child RSS. The sole authorized diagnostic replay was job **19176**: **113.739 s**,
**435,676 KiB** maximum child RSS, explicit affinity **[0]**, 4 GiB allocation and
30-minute ceiling. Slurm reserves two SMT threads of a core; the replay process
was restricted to one logical CPU. Exit 1 is the preserved scientific refusal,
not a replay PASS (`evidence/build.log:1`, `evidence/regression/diagnostic/failure.txt:1`,
`evidence/replay.log:1`, `evidence/cpu-affinity.txt:1`, `evidence/replay-final.txt:1`).

Changes are limited to new diagnostic collectors/reporting tests, build/replay/
assessment scripts and documentation/evidence. No production numerical file or
snapshot-bound historical collector was edited. Candidate status, all BA12/BA12R
failures and conditional BA13 assumptions remain unchanged. No stellar solve,
BA13 investigation, certificate promotion, master merge or unauthorized trajectory
occurred. Substantial artifacts are under
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/native-p2-checkpoint-diagnosis-20261109-r1`;
compact receipts and SHA-256 manifests are Git-tracked. See `REPRODUCE.md:1`.
