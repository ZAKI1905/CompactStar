# Phase-5D-1 structural/background resolution qualification

**Disposition A: PHASE-5D STRUCTURAL/BACKGROUND RESOLUTION FRESHLY QUALIFIED AT RADIAL 80000 — REVISED PREDECLARATION AUTHORIZED — READY TO RESUME COUPLED EVOLUTION.**

This disposition qualifies only the controlled fixture's structural/background
construction and authorizes a revised predeclaration. It does not resume
coupling, enter an ODE, generate a trajectory, restore archived code, implement
A18, begin BNV, or resolve global INV-11.

## 1. Identity, history, and fixed owner decision

| Item | Authenticated value |
|---|---|
| Canonical entry | `d019ae390be4f5e3daba05039903485cb497e397` |
| Branch | `physics/phase5d-controlled-rotochemical-evolution` |
| Worktree | `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase5d-controlled-evolution` |
| Branch HEAD before qualification | `30a15b299ddd7eedfb2f5ab824179b48a20ca3ea` |
| Original predeclaration | `697ec14610d63640866f3c0aed0a51d567cb037d` |
| Response implementation | `14b74cba1a0c7b6faec591aedf6f5c9a76fa42d6` |
| Blocker record | `30a15b299ddd7eedfb2f5ab824179b48a20ca3ea` |
| Initial plan-history commit | `36cfe1e2f7307f366175001a7d94f0b22036f1fc` |
| Controlling corrected qualification plan | `9a74b245014ea1418b24c483dc066aefb7aca14d` |
| Production EOS resolution | 8192 intervals |
| Production radial resolution | **80000**, fixed before admissible output |
| Central mass-energy density | `1.10e15 g cm^-3` |
| Fixture | ordinary NStar Track-R Structure-1 midpoint, whole star |

The initial plan correctly fixed 80000 but initially classified the historical
`V_I_validation` as reusable. Source inspection before admissible factory
assembly found that it is resolution dependent. A second plan commit corrected
the protocol, excluded the early default-producer/PB6 probes, and required a
fresh M1/M2/PB reconstruction. Every result used below was generated afterward
in new scratch. This preserved the fail-closed history and prevented the
unmodified producer's stored envelope from qualifying freshness.

The 80000 selection was an owner decision, not an output fit: it is the governed
Phase-5B/Phase-5C structural-fixture resolution. No 80000 trajectory existed.
No 120000/160000 trial or post-result choice among lower rungs occurred.

## 2. Original 10000 refusal and blocker classification

The retained 10000-point production W factory refused before the ODE with
`AccuracyGoalUnmet` at `CompactStar/Analysis/src/ChemicalResponse.cpp:804`, raw
rc=1. Its positive structural contribution `sum_j |Z_ij| E_I_j` was:

| Channel | 10000 structural term (MeV s^2) | Frozen goal (MeV s^2) | Term / goal |
|---|---:|---:|---:|
| Npe | `8.918515424724095e-12` | `1e-12` | `8.918515424724095` |
| NpMu | `3.263850588763438e-11` | `4e-12` | `8.159626471908595` |

That term alone made the predeclared 10000 calculation impossible under the
unchanged numerical gates. Separately,
`build/phase5d-audit/prepare.py` combined freshly recomputed 10000 tail evidence
with 40000/80000 background/provider/anchor characterization. The resulting
mixed-resolution certificate is scientifically inadmissible. Repository copies
of the exact refusal are `phase5d1_w_assembly_refusal.json` and
`phase5d1_refusal_excerpt.txt`; the authoritative 117-field prior report is
`PHASE5D1_FINAL_REPORT.md`.

This is not evidence that ADR-0014 is inconsistent, that W is formulated
incorrectly, that governed Phase-5C regressed, that RKF45 is unsuitable, or that
the controlled benchmark is impossible. It establishes only failure of the
10000 structural numerical budget and failure to independently characterize its
background. The original predeclaration correctly stopped before a result.

## 3. Fresh execution and field inventory

All admissible work is retained under
`build/phase5d-audit/qualification-fresh-80000/`. The fresh structural envelope,
PB6, PB7, PB9--11, PB12, PB13, and chemical characterization ran as independent
processes with separate directories and result/log files. They used the same
compiled Debug/assertions-enabled toolchain and governed sources. No process
shared a mutable StarProfile or cache. Git writes remained serial under one
writer. The raw prerequisite return codes are all zero.

| Field family | Fresh source resolution / authority | Reused numerical input? |
|---|---|---|
| EOS/profile background | EOS4096/r40000 and EOS8192/r40000 characterization; EOS8192/r80000 production | No |
| Provider and anchor | fresh EOS8192/r80000 nodal/provider/anchor evaluation | No |
| Center/onsets/refusals/tail | fresh EOS8192/r80000 profile and analytic comparison certificates | No |
| Structural A/B/K/I | complete fresh r20000/r40000/r80000 M1 construction; r80000 production | No |
| Structural validation envelope | fresh M1/M2 + PB6/PB7/PB9--11/PB12/PB13 | No |
| G/Q/Z | fresh production factories at r80000 | No |
| W and both uncertainty tracks | fresh r80000 factory with fresh `I`, `E_I`, and `V_I` | No |
| Absolute goals | accepted Phase-5C owner contract; resolution independent | Contract only |
| Constants, units, model/source code | immutable governed authority | Contract only |
| Governed baseline/candidate | read only after the primary factory passed, by regression comparison | Validation only |

No radial-dependent cached number was reused. The scratch assembly read only the
`goals` field from `phase5c2_preproduction_evidence.json`; it did not access its
stored `I`, `E_I`, `V_I`, G/Q/Z/W, or certificate numbers. The final production
certificate contains no mixed-resolution contribution. The lower-resolution
calculations are explicitly characterization evidence, not alternative Phase-5D
production candidates.

## 4. Provenance

| Item | Fresh value |
|---|---|
| Compiler | Apple LLVM 21.0.0 (`clang-2100.1.1.101`) |
| Language/configuration | C++17 (`201703`), assertions enabled |
| Platform/architecture | Darwin arm64 |
| GSL | 2.7.1 |
| Python stack | Python 3.12.10, NumPy 2.3.1, SciPy 1.16.0, mpmath 1.4.1 |
| EOS identity | FR2005 equilibrium free gas; canonical Structure-1 table generator |
| Local provider | `track-r-fernandez-reisenegger-2005-free-gas-local`; R2006 corrected charge-neutral interpretation |
| EOS8192 table SHA-256 | `7cd44c92e1e7206e0e68e3fed7e3f0ca68e79ab4517d02b96ff78b9be23d3f1a` |
| r80000 profile SHA-256 | `e9cd03b0b8449806f6c9883d75de1d3dff0cf1481675efc45a56519655d40890` |
| Chemical quadrature | onset-aware segment Gauss-Legendre orders 8/16/32 plus bisected 32 characterization; production default 16 |
| Accumulation | ordered Neumaier segment/node sums |
| Production nodes | 324176 |
| Validation nodes | 2431320 |
| Profile versions | central plus 15 sequence profiles, all version 2 |
| Chemical domain | whole positive profile with explicit pe/npe/npemu onset/refusal partition |
| Tail policy | positive monotone pe comparison; finite cut plus analytic enclosure |

The local-provider maximum relative comparison spread is
`4.829745737565275e-14`; the maximum unscaled local H condition is
`9463.508699362937`; the maximum backward absolute residual is
`5.18960522309443e-16`. The fresh independent background comparison and
provider/anchor contributions are included in `E_G`, not used as a replacement
central value. Existing compiler-portability semantics exclude only the
compiler version string from equality; no new provenance exception exists.

## 5. Fresh structural response

Species order is n, p, e, mu. A and B have units count km^2; K has units count
km^2; I has units count s^2.

| Species | A | B | K | K numerical error | I_phys | E_I numerical |
|---|---:|---:|---:|---:|---:|---:|
| n | `3.9647718619373394e59` | `1.6255270466118125e59` | `1.1793897334337809e58` | `5.760213721014721e52` | `1.3122480530141585e47` | `6.409102119567852e41` |
| p | `1.8404826326625264e57` | `5.761376707019385e57` | `-1.179389733433783e58` | `2.5487178899023356e51` | `-1.3122480530141607e47` | `2.835831103064694e40` |
| e | `1.8292534582866583e57` | `5.192854858061696e57` | `-1.0459711462551772e58` | `1.871991303390107e51` | `-1.1637998545112904e47` | `2.082871228648889e40` |
| mu | `1.1229174375868009e55` | `5.685218489556692e56` | `-1.3341858717812756e57` | `8.681980762319852e50` | `-1.484481985023382e46` | `9.660006381851622e39` |

Common reduction values are
`A_B=3.9831766882639647e59`, `B_B=1.6831408136820063e59`, and
`d epsilon_c/dq=-2.3665142309456835`. The requested structural errors are
`K_error,e=1.871991303390107e51` and
`K_error,mu=8.681980762319852e50` count km^2. The fresh empirical envelopes are
`V_K,e=4.386661132887441e53` and `V_K,mu=4.927777951659812e52`, producing
`V_I=(4.880818755395441e42,5.482892414134075e41)` count s^2. PB13 is already
inside K numerical error and was not double counted in the validation envelope.

## 6. Fresh G, Q, and Z

All three unchanged production factories accepted. Matrices use source order
n/e/mu for G and channel order Npe/NpMu for Q and Z.

G (count/MeV):

```text
[[ 2.904621518418346e55,  0,                      0                    ],
 [ 0,                     2.201380364481809e53, -1.0354115533873709e51],
 [ 0,                    -1.0354115533873709e51,  9.746440586307364e51 ]]
```

E_G (count/MeV):

```text
[[1.0219404103274857e50, 1.8349389199423292e37, 1.128861506525617e36 ],
 [1.3796979262051092e37, 1.5824359678370865e46, 2.9270347720524113e44],
 [5.4936893113781254e35, 2.927034771968926e44,  2.745206319837477e45 ]]
```

Q and E_Q (count/MeV):

```text
Q   = [[ 2.184981541996094e53,  -1.1006095941901028e51],
       [-1.1006095941901028e51,  9.743848458410244e51 ]]
E_Q = [[2.17914539749605e46,  5.478831452061045e44 ],
       [5.478831451947084e44, 2.7560651427082615e45]]
```

Z and E_Z (MeV/count):

```text
Z   = [[4.5793031807026964e-54, 5.172519910278805e-55 ],
       [5.1725199102788054e-55, 1.0268727975139168e-52]]
E_Z = [[4.499198710915918e-61, 4.285424700647402e-61],
       [4.285424702045283e-61, 2.908175030741565e-59]]
```

Every entrywise error is within the frozen Phase-5C absolute goals:

```text
G goals = [[2e50,1e40,1e40],[1e40,4e46,7e44],[1e40,7e44,6e45]]
Q goals = [[5e46,2e45],[2e45,8e45]]
Z goals = [[2e-60,2e-60],[2e-60,9e-59]]
```

Support is Neutron/Electron/Muon, rank 3. G eigenvalues are
`(9.741345082771019e51,2.2014313195171725e53,2.904621518418346e55)`
count/MeV; the infinity-norm condition estimate is `2995.7010369732106`; the
inverse perturbation rho is `3.5183255579973513e-6`. Support, rank,
conditioning, currentness, and refusal gates accepted.

The fresh C++ background-ladder quantity
`2*(|G_8192,r40000-G_4096,r40000|+|G_8192,r80000-G_8192,r40000|)` is enclosed
componentwise by the fresh production `E_background`; no lower-resolution value
is substituted into the r80000 result.

## 7. Fresh W numerical budget and validation envelope

The r80000 production W factory accepted without bypass. Channel order is
Npe/NpMu and units are MeV s^2.

| Quantity | Npe | NpMu |
|---|---:|---:|
| W | `-5.406177501704724e-7` | `-1.5845719480103649e-6` |
| structural `sum |Z| E_I` | `1.0037764595781469e-13` | `1.0027334706341596e-12` |
| Z-error `sum E_Z |I|` | `5.87233038180986e-14` | `4.815871106906613e-13` |
| cross `sum E_Z E_I` | `1.3510974542660785e-20` | `2.898558813794819e-19` |
| arithmetic | `3.8413201520659786e-21` | `1.1259060869479959e-20` |
| total E_W numerical | `1.59100967128208e-13` | `1.4843208824397632e-12` |
| frozen E_W goal | `1e-12` | `4e-12` |
| E_W / goal | `0.159100967128208` | `0.3710802206099408` |
| goal / E_W margin | `6.28531691573045` | `2.6948350907960283` |
| V_W validation | `2.2634352552795975e-11` | `5.882694393677142e-11` |
| frozen V_W goal | `1e-10` | `3e-10` |
| V_W / goal | `0.22634352552795975` | `0.19608981312257143` |

The independently summed W numerical terms reproduce the factory total to
floating summation roundoff (one ulp in Npe, exact in NpMu). Numerical error and
validation envelope remain separately classified; neither is called a
certified continuum bound.

## 8. Tail and certificate acceptance

The fresh pe tail certificate has `R_cut=12.766174760730493 km`,
`M_cut=0.9209327364386665 km`, `R_upper=12.768154903424522 km`, and
`M_total_upper=0.9209327364387035 km`. Its `G_ee` enclosure is
`3.681810845442163e45 count/MeV`; fresh direct same-cut and independent
high-precision shell values are `2.4255976742539016e45` and
`2.42559767427363e45 count/MeV`, both enclosed. Neutron and muon refusal
windows were freshly reconstructed from source guards and entire containing
profile cells. Center, refusal, tail, and background components were all
transported separately. Tail acceptance passed.

The final certificate JSON SHA-256 is
`a57401f867a45f9eb4a52064133ba352a6b201abd4aee38f98424ac1fe824d16`;
transport SHA-256 is
`7fc892b50bfbcdd2c5d963e0ff866777f3628f40fc282e1ce5a0b0364200a453`.
The consolidated fresh qualification evidence SHA-256 is
`fea0015c6e5bec54edddbc9802c1ad79940c0b247e2c1fb4bf28da592eca4550`.

## 9. Governed comparison and regressions

Only after the independent fresh-V_I production factory passed, the existing
governed machinery regenerated and compared the fixture. The equality-bearing
scientific result matches the governed Phase-5C baseline; the sole permitted
portable field is the compiler version. The governed baseline was never a
production input.

| Gate | Result | Raw rc |
|---|---|---:|
| `phase5b_structural_response_regression` | 1/1 PASS | 0 |
| `phase5c_chemical_coefficient_regression` | 1/1 PASS | 0 |

The protected artifacts remain unchanged:

- Phase-5B baseline SHA-256
  `7588f0e9cd62f5b6be48bb725e4d0ba6e47b64f50d843117e70c17a78beeb5fa`;
- Phase-5C baseline SHA-256
  `7027aa6179fe9111d76586d467d007c4223ec9d26ece0fee773327973c77bfe7`;
- reviewed Phase-5C candidate SHA-256
  `a6b430b2356987b54a7c82c87d0b00f3e1d9698f6f299c32c80efcc40c63833b`.

Because this task changes documentation and scratch evidence only, the prior
focused response 2/2, data-free 47/47, and authenticated 70/70 evidence remains
applicable and was not needlessly rerun. EOS/data and literature have no task
diff. Production response source/tests have no task diff.

## 10. Scope and authorization

All frozen goals passed without modification. The qualification changes no
Z/W physics, numerical-error definition, benchmark normalization, process set,
thermal input, spin law, run interval, ODE tolerance, refinement rule,
quasi-steady criterion, source-function root, or RE10b tolerance. No coupled
ODE was entered; no trajectory was generated; no archived coupling was
restored; no production response code was changed by this task; no evolution or
candidate SHA exists; there is no merge. A18 and BNV were not begun. Global
INV-11 remains **UNRESOLVED** pending the later controlled-evolution validation
and independent review/owner ratification.

All fresh r80000 structural/background, G/Q/Z, W numerical, W validation,
support/rank/conditioning, refusal/tail, provenance/currentness, governed
comparison, and upstream regression gates pass. Therefore amendment of the
Phase-5D controlled-evolution predeclaration is authorized, limited solely to
the replacement production radial resolution 80000.
