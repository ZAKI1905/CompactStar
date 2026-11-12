# Proposed study: when BNV chemistry competes with deposited power

**PROPOSED — owner consideration, 2026-11-12. No implementation or run authorized
by this document.** Study the competition as a function of baryon-loss drive,
retained energy and thermal/chemical state, first as a conservative phenomenological
experiment. A microscopic decay prediction is a separate milestone. The proposed
minimum is one energy-partition provider at the existing interface, with unchanged
composition and thermal-ledger owners; alternatives are in
[ADR-0025](../../../adr/ADR-0025-bnv-retained-energy-study.md).

## Governing equations and the four powers

Use infinity-coordinate time and redshifted T, η and potentials; energies per
baryon are MeV and global powers erg/s. Let c be MeVToErg, D=−Bdot>0, b=(1,1,1),
P select the e/μ components, t=dN_eq/dB, and Z the symmetric positive chemical
matrix. Spin is off and Z,t,μ_B are frozen. The governed equations are

```text
sigma = P(S_y − t Bdot),                 etadot = −Z(R(T,eta) + sigma)
S_y = (−D,0,0)  =>  sigma = D t_l,       mu_n,actual = mu_B + t_l·eta
Echem = (c/2) eta^T Z^−1 eta,            LH = c eta·R >= 0
A = −c eta·sigma,                       Edot_chem = A − LH
P0 = c mu_B D − Lout,fluid,              Pactual = P0 − A
Cstar Tdot = Pactual + LH − DeltaLnu − Lnu,eq − Lgamma
           = P0 − Edot_chem − DeltaLnu − Lnu,eq − Lgamma.             (1)
```

Here **A is source work into the chemical reservoir**, not its net derivative;
P0 denotes equilibrium-potential net deposition after fluid exit (not the named
P0 fixture). Storage is positive Edot_chem; release is negative. Extra neutrino
cooling is DeltaLnu=Lnu,full−Lnu,eq at the **same state**, not the source/control
neutrino difference caused by different temperatures. Thus the four requested
powers are Pactual, LH, DeltaLnu and Edot_chem, with the transfer A linking them.
Do not add −Edot_chem to a thermal equation already containing Pactual+LH.
Authority: `docs/adr/ADR-0015-bnv-open-system-thermal-ledger.md:144`, `:211`, `:233`;
implementation: `CompactStar/Physics/BNV/src/FrozenControlledBnvRunContext.cpp:150`.

For prompt escape, Lout,fluid=Lesc,star and no hidden reservoir is required:

```text
Edot_eq = −c mu_B D
( Eeq + U + Echem )dot = −Lesc,star − Lnu,full − Lgamma.             (2)
```

If products remain in a nonthermal X reservoir, use Lout,fluid=Lesc,star+J_X and
Edot_X=J_X (minus separately owned releases/losses). Ignoring X while treating
its transfer as star escape would change the physics (ADR-0015:93).

## Dimensionless comparison (derived here from equation 1)

Normalize to P*=c μ_B D, which stays useful even at zero retained power. Define

```text
p=P0/P*, d=Pactual/P*, a=A/P*, h=LH/P*, n=DeltaLnu/P*, s=Edot_chem/P*
a=−t_l·eta/mu_B,        d=p−a,           s=a−h
Cstar Tdot/P* = d+h−n−ell = p−s−n−ell,
ell=(Lnu,eq+Lgamma)/P*.                                           (3)
```

The direct comparison is (h/d, n/d, s/d), retaining the sign of s and reporting
a/d as the source-transfer diagnostic. Gross reaction heating competes when
h~d; net weak reactions compete when |h−n|~d. The **combined** chemical response
relative to P0 competes when |s+n|~|p| and heats only if s<−n. At chemical
stationarity s=0, it is −DeltaLnu, even if LH−DeltaLnu>0. Ratios with d=0 or
p=0 are undefined: report the P*-normalized vector instead. With D=0 use raw
powers or normalize to passive cooling; an initial stored imbalance can still
relax and heat. For observational relevance also compare all powers to
Lnu,eq+Lgamma: parity with a tiny BNV deposition need not matter thermally.

A useful per-destroyed-baryon form is q_dep=Pactual/(cD),
a_E=A/(cD)=−t_l·eta, and r=DeltaLnu/LH when LH>0:

```text
LH/Pactual = (eta·R/D)/q_dep
DeltaLnu/Pactual = r LH/Pactual
Edot_chem/Pactual = a_E/q_dep − LH/Pactual.                        (4)
```

At a reached stationary state R=−D t_l: LH/(cD)=a_E and the gross threshold is
q_dep~a_E; net weak parity is q_dep~|1−r| a_E. These are instantaneous power
criteria, not integrated ΔT fractions or proof of reachability.

## Drive, temperature, imbalance and equilibration

The existing nonsuperfluid modified-Urca model has R_j=K_j T^7 H_M(xi_j),
DeltaLnu_j=Ltilde_j T^8[F_M(xi_j)−1], xi_j=eta_j/(kB T), with consistent
normalization K_j=Ltilde_j/(c kB). Its odd H_M contains xi,xi^3,xi^5,xi^7;
F_M−1 contains xi^2,xi^4,xi^6,xi^8. The accepted coefficient authority is
ADR-0015:268 and ADR-0014; the exact independent implementation used in the review
is `scripts/model.py:67` at the original review root (`review-evidence/IMPORT.json:1`).
Published context is [Fernández & Reisenegger (2005)](https://arxiv.org/abs/astro-ph/0502116),
*Rotochemical Heating in Millisecond Pulsars: Formalism and Nonsuperfluid Case*;
the BNV moving-reference extension is governed by ADR-0015, not inferred from that paper.

| Regime | Derived scaling and implication |
|---|---|
| Small \|xi\|, fixed T | R≈Λ(T)η with Λ diagonal and positive, Λ∝T^6. ηdot=−ZΛη−ZD t_l. For initial zero η, η(t)=[I−exp(−ZΛt)]η_* with η_*=−Λ^−1 D t_l. Relaxation times are inverse eigenvalues of ZΛ. |
| Early drive-dominated growth | η≈−ZD t_l t. A∝D²t; LH∝D²T^6 t². At fixed q_dep their ratios to direct power scale as D t/q_dep and D T^6 t²/q_dep respectively. |
| Reached linear stationary state | LH*=cD² t_l^TΛ^−1t_l ∝D²T^−6; LH*/Pactual∝D/(T^6 q_dep). Lower T raises stationary imbalance, but also increases relaxation time and can prevent reaching it. |
| Small \|xi\| neutrino balance | r→22020/14680=3/2, so net weak reactions cool: LH−DeltaLnu≈−LH/2. |
| Nonlinear, large \|xi\| | R∝eta^7, LH∝eta^8 and r→9/24=3/8. At fixed coefficients stationary \|eta_*\|∝D^(1/7), LH*∝D^(8/7), LH*/Pactual∝D^(1/7)/q_dep. The small-xi T^−6 scaling no longer applies. |

The per-channel net-heating threshold is |xi|≈4.909710 for this modified-Urca
model; mixed channels require summing powers, not testing only the largest xi.
Use tau_relax=1/eigenvalue[Z (∂R/∂eta)] in the applicable regime and compare to
both T/|Tdot| and D/|Ddot| (infinite for constant D), as well as elapsed time.
Temperature feedback requires coupled evolution. No superfluid or direct-Urca
threshold is inferred (ADR-0015:268).

The retained review gives tau=13,196–16,170 yr at T_end, linear stationary
η*= (−1.73,−3.50)e−6 MeV and LH*/P*=7.5222e−11. Thus a_E*≈6.53e−8 MeV
(0.0653 eV) and full-retention chemistry is negligible. At 2,048 yr a=9.597e−12
and h=1.234e−12. These are existing fixed-temperature diagnostics, not new
long-horizon predictions (`review-evidence/outputs/longterm.json:1`,
`review-evidence/outputs/trajectory-check.json:1`). The review's ~25 GeV imbalance
for parity with full retention lies far outside the accepted small-departure
model. Its negative conclusion applies to the driven P2 benchmark, not arbitrary
initial reservoirs, stronger sources or microscopic processes.

## Retention semantics and the minimum extension

**Recommended option for owner consideration:** interpolate the existing actual-
potential P0/P2 energy partitions with a constant 0≤f≤1, preserving D and sigma:

```text
Eesc,fluid = Eesc,star = (1−f) mu_n,actual     [MeV per event]
EX = Ein = 0
q_dep = f mu_n,actual
Pactual/P* = f(1−a),       P0/P* = f+(1−f)a.                       (5)
```

This is a fraction of the **actual chemical-energy budget**, not a particle
branching ratio, not a fraction of product kinetic energy, and not an independent
extra heat source. At f=1 it recovers Stage B; at f=0 it recovers the existing
cold P0 cancellation with zero direct residual, while chemistry still evolves.
Pactual−P0=c eta·sigma and equation (2) hold for every f. Gross weak parity near
stationarity occurs at f≈a/(1−a), ~7.5e−11 at the review's fixed T; net weak parity
in its linear regime is ~3.8e−11. These numbers only seed an eventual predictor:
reducing retention changes T, eta, r and the achieved state.

The tempting alternative **P0=f μ_B cD** instead requires
Eesc=(1−f)μ_B, q_dep=fμ_B+t_l·eta, and d=f−a. It is a different, also algebraically
conservative experiment; in a filling state f<a gives negative direct residual.
It cannot label f=0 an ordinary nonnegative cold deposition channel without a
separate energy-extraction interpretation. Equation (5) avoids that issue and
connects existing fixtures. Neither option may multiply only Pactual in the
thermal RHS while leaving Lout and P0 unchanged. Owner selection is pending.

The existing `IDirectBnvEnergyPartition::Evaluate` already receives actual
potentials (`CompactStar/Physics/BNV/DirectEnergyLedger.hpp:33`); the ledger owns
q=−nu·mu_actual−Eesc+Ein and the sole conversion to erg
(`CompactStar/Physics/BNV/src/DirectEnergyLedger.cpp:73`). P0/P2 providers already
set Eesc to mu_n,actual/zero (`tests/bnv/energy_partition_fixtures.hpp:28`).
The smallest extension is therefore a parameterized partition plus explicit
run-card/provenance and prompt-escape fate binding, not a new ODE or a new
chemical-work term. Keep the original full-retention run card and evidence immutable.
A future implementation must evaluate small residuals stably: for f~1e−10,
subtracting two ~868 MeV numbers loses precision. An algebraically factored
residual in the energy-ledger owner, checked against independent high-precision
q/P0/R18 calculations, may be necessary; never hide cancellation by relaxing gates.

## What makes a microscopic study different

A physical process must jointly supply stoichiometry, occupied-state event rates,
product spectra, medium/Pauli effects and escape/absorption or accumulation.
Particle and energy integrals use the same events and their respective lapse
weights (counts e^Phi, powers e^(2Phi)). Its D, sigma and q_dep generally cannot
be varied independently. Fermi-hole energy is already inside mu_actual−Eesc;
adding it again double-counts work (ADR-0015:93, :295).
[McKeen, Pospelov & Raj (2021)](https://arxiv.org/abs/2105.09951) provides a concrete
neutron-disappearance heating example involving refilled Fermi-sea vacancies;
it does not validate equation (5) as a microscopic decay model. An invisible
product is not generally zero heating, and retained massive products may require
an X reservoir and feedback rather than prompt-escape bookkeeping.

The finite-temperature issue becomes decisive at small q_dep: the ~0.065 eV
stationary competition scale above is below even the illustrative smooth-source
1.2–2.4 eV/event floors at 10^8 K, and far below a surface-weighted O(kBT) term.
Stage B's ~6.26 K omitted scale is specific to its history, not an error bar to
reuse after varying f. The new partition must declare its weighting assumption
and omitted-power scale independently of f; it must not multiply the omission
by f merely because the cold residual is small. Below that floor the study maps
**cold-model competition**, not a physical thermal sign. A microscopic prediction
needs finite-T source/state derivatives or an explicit demonstration that they
are negligible (ADR-0015:60). The catalog and README at the AGENTS.md macOS
literature path were unavailable on EKU; accepted local contracts and the linked
primary papers are used, with no claim to have inspected those local PDFs.

## Smallest proposed study after a separate approval

1. At fixed admissible T and D, solve R(T,η*)=−D t_l and evaluate relaxation,
   energy ratios and the omitted-term floor. Select f around gross/net parity
   and choose horizons from relaxation and frozen-validity budgets, not an
   arbitrary long run. Start with zero initial imbalance; reservoir-release
   experiments need a separately declared initial energy and driving history.
2. Only then propose a small coupled matrix varying D, initial T and f, including
   f=1 regression, f=0 mathematical control, and exact-zero source. Report the
   four powers, A, all closure identities, passive-cooling comparison and ΔT;
   keep channel diagnostics distinct from physical counterfactuals.
3. Regenerate domain boxes, Cstar event coverage and numerical constants for
   each approved scope. The old source exp(LH)<2 bound reaches its limit around
   3,202 yr; ~50,000 yr stationary tests are not qualified by Stage B. Lower
   retention can cause downward crossings, not covered by Stage B's two upward
   events. Retain frozen depletion/population/response guards; stop or propose
   structural evolution if they fail. The current thermal certificate spans
   10^6–10^9 K (`CompactStar/Physics/BNV/src/FrozenControlledBnvRunContext.cpp:105`).

No trajectory, predictor solve or qualification campaign was run for this proposal.
Resources for future approved work should follow measured needs and EKU scheduling,
with independent calculations parallelized where useful and no arbitrary CPU cap.
