# Phase-6A-0 Cowling baryon-direction diagnostic

**DERIVED / NUMERICALLY SUPPORTED / NOT RATIFIED AS PRODUCTION PHYSICS.**
**Date:** 2026-09-13; R3 precision/scope clarification 2026-09-14. **Scope:** non-governed Phase-6 scientific note; no coefficient replacement.
Canonical master `0a7418aecb7314cfa472a78f1faf477be8456a94`; draft entry
`5a6bf7cb9455d684ddb6fccb22ad2b9fec940b3a`.
Companions: [preflight with R3 corrections](PHASE6A0_BNV_THERMAL_FIRST_LAW_PREFLIGHT.md) and
[PROPOSED ADR-0015](../adr/ADR-0015-bnv-open-system-thermal-ledger.md).

## 1. What this note establishes

Phase-5C is correct and governed for its declared R2006/Cowling contract. The new Phase-6
question is how total baryon loss moves the equilibrium star. A fixed-metric chemical response
cannot determine that hydrostatic baryon direction. This does not invalidate the governed
fixed-baryon beta machinery or authorize rewriting Phase-5C. R2006 Appendix/footnote 4's
Cowling justification concerns baryon-conserving perturbations. The baryon-direction mismatch
is PHYSICAL-MODEL uncertainty/scope, not `numerical_error`; PROPOSED ADR-0015 section 7 explicitly
narrows only ADR-0013 Q1 / ADR-0014 section 3.17's forward-looking BNV seam.

Evidence owners: ADR-0013 defines G_y/Z; ADR-0011 defines structural response; R1 independently
extracts t and negative oracles from unchanged governed baselines. The supplied independent
report `/Users/keeper/Downloads/PHASE6A0R0_CONSOLIDATED_REPORT.md` (R0), SHA-256
`c1df790a336a0b46aeb37154818a09ff44c5bcdb366f61697b3b6d5f00d41e92`, Q21-Q25/Q54 G13,
reports off-beta-equilibrium TOV results. Its underlying full G_true matrix and scratch scripts
were not supplied, so those numerical results are **attributed R0 evidence**, not an R1
regeneration. The equations below are independently checked. No new baseline is installed.

## 2. Physical tangent versus fixed-metric direction

Use y=(N_n,N_e,N_mu), b=(1,1,1), L=[[-1,-1],[1,0],[0,1]], P selecting e/mu.
The physical tangent is `t=(partial N_eq/partial B)_Omega`, b^Tt=1. At spin OFF,
`t_i=B_i/B_B` from independent neighboring-star structural derivatives. The Cowling direction
is `k=G_y b/(b^T G_y b)`; it changes redshifted baryon potential at fixed metric.

ON THE STRUCTURE-1 FREE-GAS FIXTURE, the long digits below are arithmetic-reproducibility
oracles on governed bytes, NOT physical precision. Drive coefficients are eta_dot/abs(Bdot)
in MeV/count; multiplying by abs(Bdot) in count/s gives MeV/s. Ratios are fixture-specific.

| Quantity | R1 governed-byte arithmetic oracle |
|---|---|
| t | (0.9657700849496014,0.030852171225661786,0.0033777438247248118) |
| k | (0.9922178920029234,0.007484539490838682,0.0002975685062378332) |
| t_e/k_e; t_mu/k_mu | approximately 4.12;11.35 |
| Physical neutron-sink eta_dot/abs(Bdot) | (-1.4302859054e-55,-3.6280967205e-55) MeV/count |
| REJECTED raw-G neutron-sink eta_dot/abs(Bdot) | (-3.442789339881122e-56,-3.442789339881122e-56) MeV/count |
| REJECTED raw-G physical-slide eta_dot/abs(Bdot) | (+1.0860069714e-55,+3.2838177865e-55) MeV/count |

The physical slide is S=t Bdot, giving sigma=P(S-t Bdot)=0. For baryon-changing S, the
rejected raw route `-L^T G_y^-1 S=-Z(S_l-k_l Bdot)` agrees with the physical route only if
t=k. Baryon-conserving sources agree regardless. Thus the old number is explicitly a
**COWLING k-ROUTE NEGATIVE ORACLE**, never a physical BNV drive.

Structural provenance: `CompactStar/Analysis/src/ParticleNumberResponse.cpp:365` and `:450`;
ADR-0011 sections 3-4; Phase-5B baseline hash
`7588f0e9cd62f5b6be48bb725e4d0ba6e47b64f50d843117e70c17a78beeb5fa`.
G/Z owner: `CompactStar/Analysis/src/ChemicalResponse.cpp:705`; Phase-5C baseline hash
`7027aa6179fe9111d76586d467d007c4223ec9d26ece0fee773327973c77bfe7`.
Conservative ratio numerical budgets are(1.443e-7,4.353e-9,2.179e-9), not certified intervals;
relative t budgets approximately (1.5e-7,1.4e-7,6.5e-7), supporting about 6–7 significant
digits. Drive budgets must also propagate governed Z errors, with Cowling model uncertainty
separate. The fixture's d ln N_mu^eq/d ln B approximately 54 (R2 E-9) requires an explicit
muon population/t_mu/Z_npmu sensitivity budget for depletion; it is not a universal EOS value.
Sum-rule residuals are approximately 1.2e-14. See preflight section 6 for the propagation formula.

## 3. Projection identity, valid also for an indefinite global response

For symmetric invertible G, a_G=b^TGb nonzero and invertible Z=L^TG^-1L, define
H=G-Gbb^TG/a_G. Hb=0, so H=L C L^T with C=P H P^T. Multiplication gives
`ZC=L^T[I-bb^TG/a_G]P^T=I`. Hence

```text
Z^-1=P[G-Gbb^TG/(b^TGb)]P^T.                                 (D1)
```

R1 verifies D1 on governed G to relative residual approximately 2e-16. It is an algebraic
identity, not a proof that any G supplies the physical equilibrium direction. D1 is the
ADR-0013 section 3.3 Q identity; symmetry is not needed for D1 but is needed for D2 below.
For a true equilibrium response, g_eq=mu_B b gives
`G_true b=(dB/dmu_B)t`. Therefore

```text
Z_true^-1=G_true^(e,mu)-(dB/dmu_B)t_l t_l^T.                  (D2)
```

The projection removes the total-baryon direction before the inverse beta response is formed.
A true equilibrium-energy Hessian need not be positive in the total-baryon direction even
when the fixed-baryon beta curvature is positive. Symmetry is evidence of integrability/
first-law consistency, not by itself proof of a global thermodynamic theorem.

## 4. New R0 off-equilibrium TOV evidence, not production authority

R0 reports independently solved free-gas TOV configurations with three constant neutral
redshifted potentials and differentiation of their global species counts:

| Reported diagnostic | R0 result / provenance |
|---|---|
| G_true symmetry | about 1e-6 in one calculation;1.4e-5 in another rerun (R0 sections 3,6) |
| Eigenvalues | approximately(-4.1e54,+9.9e51,+2.5e53) count/MeV |
| b^T G_true b | approximately-4.66e54 count/MeV |
| dmu_B/dB | approximately-2.147e-55 MeV/count |
| G_true b/(b^T G_true b) | approximately(0.965778,0.030851,0.003370), near structural t, not within-budget agreement |
| Fixed-metric contribution to t | reported -6.28k, with metric/volume remainder +7.28 in the baryon sum |
| Cowling Z diagonal error | approximately 1.9% npe,0.9% npmu |
| Off-diagonal | Z_true approximately 2.34e-55 vs Cowling5.17e-55 MeV/count, ratio about 0.45 |
| Impact on the reported BNV drive/W | only a few percent, despite much larger relative off-diagonal error |
| Projection explanation | about 93% of off-diagonal discrepancy attributed to total-baryon projection; hybrid true projection approximately reproduces Z_true |

The displayed muon direction differs from structural t_mu by about 0.24%, much larger than
the Phase-5B propagated budget; "near" must not be read as precision agreement.
These reported values have limited displayed precision and no newly established R1 convergence
certificate. A final reviewer should obtain the full R0 scratch/matrix or independently solve
the off-equilibrium TOV problem before accepting quantitative error claims. They are not
production tolerances or replacement coefficients. The t/k physical mismatch and raw-map refusal
already follow independently from the two authenticated governed input sets.

## 5. Energy consequences and future oracles

Actual moving-reference individual potentials follow the first-order identity
`delta g=-(I-bt^T)P^Teta`, in particular delta mu_n=t_l^Teta. Cowling reconstruction would
instead use k_l and leave a spurious O(eta Bdot) residual. The physical controlled ledger has
no such residual to count as heat. A future check may deliberately install the k reconstruction
to detect this error. It may not install that reconstruction as a physical thermal model.

REJECTED predecessor objects: fixed-reference alpha/a/alpha^2/(2a), E_2 and a free F_ref
with Delta(delta g^T F) as a whole-star thermal-flow addend. A fixed-background algebra fixture
can have those mathematical state/work terms, but it is not the physical moving-sequence
contract. The true sequence curvature lives in E_eq(B). Retain only
`E_chem=eta^TZ^-1eta/2` for fixed-current-B disequilibrium, with the R1 variable-Z chain-rule
correction in the preflight section 8.

Proposed future oracles: t sum rules/currentness; physical sliding null; raw-k negative values;
D1/D2 projection; independently converged off-equilibrium star derivatives; derivative symmetry;
fixed-B positive beta curvature versus unrestricted baryon curvature; actual-potential ledger
closure. None is implemented here. Phase-5B/C/D code, governed baselines, historical candidates,
EOS and all literature bytes remain unchanged. This note is a Phase-6 scientific finding for
review, not Phase-5 replacement authority.
