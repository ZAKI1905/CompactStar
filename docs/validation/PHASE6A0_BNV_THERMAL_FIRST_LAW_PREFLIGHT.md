# Phase-6A-0 BNV thermal first-law and chemical-response preflight

**Status: PROPOSED SCIENTIFIC PREFLIGHT; NOT HUMAN-RATIFIED.**
**Date:** 2026-09-13. **Change class:** documentation containing proposed
scientific-semantic and architecture contracts; no executed scientific behavior changes.
**Canonical entry:** `0a7418aecb7314cfa472a78f1faf477be8456a94`
(`docs: close controlled rotochemical evolution`).
**Branch:** `analysis/phase6a0-bnv-thermal-first-law-preflight`.
**Worktree:** `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase6a0-bnv-thermal-preflight`.
**Decision companion:** [ADR-0015](../adr/ADR-0015-bnv-open-system-thermal-ledger.md), PROPOSED, NOT ACCEPTED.

## 1. Entry, authority, and scope

Entry authentication found the canonical checkout at
`/Users/keeper/Documents/CompactStar/repo/CompactStar`, clean on master, and
local master = origin/master = live `refs/heads/master` = the entry SHA above.
Existing worktrees and branches were inspected first; neither requested Phase-6 path nor
branch existed. The new worktree was created directly from that exact commit. No other
branch contained either new document. The eleven governed baselines were hashed at entry.
The Phase-5D baseline `tests/baselines/phase5d1_controlled_evolution.json` remains
SHA-256 `2606916915b2da5c051b1a06637c0a371a74751c6e77b2775bc629a63bc9f6dd`.

Canonical Phase-5D controlled frozen-v1 is implemented, validated, independently reviewed,
human-ratified, fresh-context reproducible, governed, canonically integrated, and closed for
its declared controlled scope. Some upstream headers still describe the pre-fast-forward
branch state: their dated history is preserved, while the authenticated entry commit and
owner's entry instruction establish current canonical inclusion. This does not close global
INV-11, realistic FR2005/A18, or authorize production BNV
(`docs/validation/PHASE5D1_CONTROLLED_ROTOCHEMICAL_EVOLUTION_INTEGRATION.md:20`,
`docs/validation/PHASE5D1_CONTROLLED_ROTOCHEMICAL_EVOLUTION_INTEGRATION.md:210`,
`docs/SCIENTIFIC_INVARIANTS.md:1004`).

The eight requested standard-authority documents were read completely before derivation:
ADR-0013, ADR-0014, Phase-5C2 integration, Phase-5D1 ratification and integration,
SCIENTIFIC_INVARIANTS, MODERNIZATION_ROADMAP, and CURRENT_ARCHITECTURE. GOVERNANCE and
AGENTS were also read. Relevant Analysis/Rotochemical APIs were inspected. The source and
claim hierarchy remains `GOVERNANCE.md:19`; proposal is not accepted authority.

### Specialist arrangement and actual runtime limitation

E (root/lead) is the **sole repository writer**, including all Git actions and both documents.
A, B, C were launched at the outset as read-only Astra HIGH specialists: A thermodynamics/GR,
B entire Zakeri corpus, C unreduced response/storage. The runtime permits four total concurrent
slots including E and subsequently rejected a fourth specialist thread even after A completed
(`agent thread limit reached`). A therefore performed the separately labeled D energy-ledger
red-team pass while B/C continued. This is **three specialist identities covering four roles**,
not four independent concurrent reviewers. A and D must not be counted as independent identities.
E independently recomputed the matrix/sign/storage checks. This internal preflight is not the
requested future independent review. No Fable escalation was warranted by a source-backed
unresolved disagreement; missing sources are not a model-voting question.

### No implementation

No production BNV process, rate, lifetime, cross section, operator coefficient, dark coupling,
EOS, A18 reconstruction, trajectory, numerical scan, test, or baseline is introduced. No
Phase-5B/5C/5D machinery or governed document is amended. Legacy BNV extension files are not
activated or promoted to this new framework's authority
(`docs/architecture/CURRENT_ARCHITECTURE.md:341`). All proposed formulas below are derivations
for review, not claims that an existing API executes them.

## 2. Authenticated source manifest

Shared read-only root: `/Users/keeper/Documents/CompactStar/literature`.
The library has 20 catalogued PDFs. Its checksum list contains 22 entries, including README
and catalog. All **22/22 passed** `shasum -a 256 -c literature/SHA256SUMS.txt` from the parent
CompactStar directory. Individual source hashes were checked before use. The checksum list's
own hash is a recorded authentication fingerprint, not a self-signature or independent trust
certificate. Catalog roles and supersession remain controlling (`literature/README.md:13`).

| Manifest | SHA-256 |
|---|---|
| `catalog.tsv` | `285040a9931bd3be20eedcc69900a0afd520ee702844c2d7771e140de12223f4` |
| `SHA256SUMS.txt` | `567f21e661bc16dd4ca0dc188ce3188d2e1f8167077d02b2cb7b3a6327181d9d` |

All source paths below are relative to that root. Page numbers are one-based PDF pages unless
explicitly labeled printed/journal. Source text and equations are distinguished from new
algebraic consequences. No acquired source was installed; bounded web metadata discovery is
not catalog promotion.

| ID | Exact source PDF | SHA-256 | Role |
|---|---|---|---|
| R95 | `rotochemical/1995-Reisenegger-Rotochemical-Heating.pdf` | `9af85e37c7a52fd5b704c0ba07cc0ad89741d23b049df31cb6867d501d91d0ff` | Original standard framework |
| R97 | `rotochemical/1997-Reisenegger-Constraining-Dense-Matter-Superfluidity-through-Thermal-Emission-from-Millisecond-Pulsars.pdf` | `19a10133511aefc05ece33d7454cba60def03de85d9980c55fd7229058c2d08b` | Supporting context |
| FR05 | `rotochemical/2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf` | `f184d7d1d7030b61a021eb5c7ac14b1f1b30c7ea69e9d53473d153cfb069ea88` | Primary non-superfluid formalism |
| R06 | `rotochemical/2006-Reisenegger-Rotochemical-Heating-of-Neutron-Stars-Rigorous-Formalism-with-Electrostatic-Potential-Perturbations.pdf` | `a286f15e083e52becd95b3000cbb5ec3ed97148681cf10a43f1a1cc5c4d23ae8` | Corrected electrostatic/neutral response; supersedes naive FR05 inverse |
| JRF06 | `rotochemical/2006-Jofre-Reisenegger-Fernandez-Gravitochemical-Heating-Model.pdf` | `2dd5444d19cebae12509fe4ecb7dac31957d332e2131813894a10cceac403109` | Fixed-baryon external-forcing precedent only |
| Y20 | `rotochemical/2020-Yanagi-NS-Therm-Thesis.pdf` | `69590539c275fa679a5521a9c5abedd9fdc58718b1554827786c8f707bc618cc` | Supporting crosscheck, not convention origin |
| H67 | `rotation/1967-Hartle-I.pdf` | `ed263946e9bc13842399b5c9e9c2eae31823e7323bc81b456fb5174697cefc35` | Cold first integral and baryon measure |
| H70 | `rotation/1970-Hartle-IV.pdf` | `2836d50173580f8923fb3b21c8e5005dfadfbc8bf0f6b22daeb1c4a92dde4cc7` | Restricted rotating-star variational identity |
| B22 | `bnv/zakeri/2022-Berryman.pdf` | `4446c3c28c09205f267d71fd67ce7d375499c997096f458c9a30b9ec6f66162d` | User-authored BNV corpus; 62 pages |
| A24 | `bnv/zakeri/2024-Allahverdi.pdf` | `f66470df4d217892f8efdc8fe7dd6cc1ea8300e3fe9d0924d219b21c05adae49` | User-authored BNV corpus; 21 pages |
| G24 | `bnv/zakeri/2024-Gardner.pdf` | `9611f56e18708db559eb6e868587edf0b7e578a20d38338b2248497f738c72ee` | User-authored BNV corpus; 33 pages |
| Z24 | `bnv/zakeri/2024-Zakeri.pdf` | `ecd364655f04d5585c895183fa38259520bdb9e80a8de2e9f68680f0d24f1a55` | User-authored BNV corpus; 29 pages |

B audited all four BNV PDFs (145 pages total), including references, using complete text reads
and targeted rendered equation checks. This is the exact relevant Zakeri corpus in the
catalog, not a claim to cover every BNV publication or untracked owner manuscript.

## 3. Standard contract and notation

Use the repository metric `(-,+,+,+)`,
`ds^2=-exp(2nu)c^2 dt^2+exp(2Lambda)dr^2+r^2 dOmega^2`, with lapse normalized to
infinity and `exp(Lambda)=(1-2m/r)^(-1/2)`. This is INV-03/04
(`docs/SCIENTIFIC_INVARIANTS.md:227`, `CompactStar/Geometry.hpp:1`).
Set `c=1` for tensor contractions; use seconds, cm, erg, and Kelvin in the local ledger.
The symbol `s` is entropy density, not entropy per particle; species chemical potentials
include rest energy. Gamma is signed creation, whereas channel event density R_a is
nonnegative for a declared direction; reverse channels may be separate events. Net beta
rate DeltaGamma_l can have either sign. A dot denotes infinity-coordinate time unless D is
shown, and `dV` always means proper volume.

Chemical matrices use MeV and counts; thermal powers use erg/s. Every displayed `g S`,
`eta R`, or chemical energy in the global energy ledger is multiplied once by
`C_M=Units::MEV_FM3_TO_ERG_CM3/10^39` when crossing into erg. No new numerical unit constant
is defined (`CompactStar/Physics/Rotochemical/ChemicalImbalanceState.hpp:13`). For readability
some algebra uses a common energy unit and writes that conversion implicitly.

The standard signs, rate integral, heat, and neutrino increment are fixed by
`docs/adr/ADR-0014-secular-rotochemical-evolution-contract.md:110`, `:150`, `:191`, `:253`.
Cold response and named source basis are fixed by
`docs/adr/ADR-0013-corrected-rotochemical-chemical-coefficients.md:74`, `:113`, `:142`.
No FR05 naive full four-species inverse is revived; R06 eqs. (10)-(19), PDF pp.2-4,
remain corrected authority. No second electrostatic projection is allowed.

## 4. Independent local first-law derivation

Let `D=u^mu nabla_mu`, `theta=nabla_mu u^mu`, `u.u=-1`, and describe the locally
thermodynamic ordinary-matter subsystem by

```text
T_m^{mu nu} = (epsilon+p) u^mu u^nu + p g^{mu nu}
nabla_mu T_m^{mu nu} = Q_m^nu
nabla_mu(n_i u^mu) = Gamma_i                          (P1)
q_E = -u_nu Q_m^nu .                                 (P2)
```

Positive Gamma creates ordinary species. Positive q_E transfers **total energy into ordinary
matter**, including rest energy in the same convention as epsilon and mu. It is not defined
as heat with the chemical contribution already removed. Calling q_E "nonchemical deposited
heat" and then subtracting mu Gamma would be a semantic error.

Contract P1 with u_nu. The normalization derivative `u_nu D u^nu=0` gives
`u_nu nabla_mu T_m^{mu nu}=-D epsilon-(epsilon+p)theta`, so

```text
D epsilon+(epsilon+p)theta=q_E.                       (P3)
```

Apply the assumed local Gibbs differential and Euler extensivity identity:

```text
d epsilon=T ds+sum_i mu_i dn_i
 epsilon+p=T s+sum_i mu_i n_i.
```

Substitute `D n_i=Gamma_i-n_i theta` into P3. The compression terms cancel by Euler,
leaving, without any assumed sign of the residual,

```text
q_s := T [D s+s theta]
     = T nabla_mu(s u^mu)
     = q_E-sum_i mu_i Gamma_i.                        (P4)
```

This derivation was made independently by A and checked by E; the standard-beta limit below
is a separate falsifier. Gamma has units count cm^-3 proper-s^-1, mu erg/count, s erg K^-1
cm^-3, and q_E/q_s erg cm^-3 proper-s^-1. For multiple transport currents or a dissipative
stress tensor, the entropy current and constitutive terms must be extended explicitly.

P4 is an **open-matter entropy balance**. It need not be nonnegative: escaping products carry
energy and entropy outside this subsystem. It is not the total irreversible entropy production
of matter plus products. A claim about the second law requires all entropy reservoirs. It
also does not put an unresolved nonthermal reservoir into an equilibrium epsilon(T,n).

### 4.1 Entropy source versus temperature source

Within the proposed leading-degeneracy, cold-chemistry, frozen-thermal-coefficient fixture,
call P4 `q_direct` after prompt thermalization and use it in the controlled thermal ledger.
This is a declared approximation compatible with Phase-5D, not an exact finite-T EOS claim
(`docs/validation/PHASE5D1_CONTROLLED_ROTOCHEMICAL_EVOLUTION_RATIFICATION.md:262`).
At fixed local volume, an exact finite-temperature EOS instead gives

```text
c_V D T = q_E-sum_i [mu_i+T (partial s/partial n_i)_T] Gamma_i.       (P5)
```

Derivation: write `ds=(c_V/T)dT+sum_i s_,ni dn_i` in P4. For expansion add the
corresponding reversible entropy/density-volume terms; for global changing geometry add
profile, lapse and heat-capacity derivatives. P5 is not silently installed into Phase-5D.
For Sommerfeld `s=gamma(n)T`, `mu(T)=mu_0-gamma_,n T^2/2` and the bracket becomes
`mu_0+gamma_,n T^2/2`. Neglected composition-entropy terms are O(T^2 Gamma) and can dominate
when `mu-E_esc` nearly cancels. A zero direct entropy residual is an exact zero-direct-heat
oracle only **within the declared controlled thermal model**; it is not an unconditional
zero-temperature-change theorem.

## 5. Per-event particle and energy ledger

For directed channel a declare its identity, ordinary-matter stoichiometry nu_ia,
local proper event density R_a(r,t), and

```text
Gamma_i^a=nu_ia R_a;   Gamma_i=sum_a Gamma_i^a.         (P6)
```

The matter boundary must be declared: which equilibrated particles, fields, photons, and
nonthermal products are inside it, and which are tracked separately. In the simplest prompt
closure define external/input energy E_in,a, final truly escaping energy E_esc,a, and net
energy J_X,a transferred into an explicitly separate retained nonthermal reservoir. Then

```text
q_E^a=(E_in,a-E_esc,a-J_X,a)R_a
q_direct^a=[E_in,a-E_esc,a-J_X,a-sum_i mu_i nu_ia]R_a. (P7)
```

J_X is a net transfer, not necessarily a permanent loss. Later X thermalization returns
energy through q_E once; X escape returns through its own transport ledger. If X is retained
inside the declared matter T_m but is not equilibrated, P4 needs extra state variables rather
than using P7 as though X were heat. For the prompt controlled toy, set J_X=0 by declaration.
Do not identify stored **composition** energy with X: composition is already in mu_i and the
state function of section 10.

Event conservation relates intermediate deposition to initial/final particle energies. It
forbids independently choosing every quantity called "energy per event."
For destruction of one occupied particle of total energy E_occ, with no other matter
participant, and independently supplied E_in,

```text
Gamma_i=-R
q_E=(-E_occ+E_dep)R
E_occ+E_in=E_dep+E_esc+J_X
q_direct=(mu_i-E_occ+E_dep)R
        =(mu_i+E_in-E_esc-J_X)R.                     (P8)
```

Here E_dep is energy actually returned by products to equilibrated matter, not a name for
the whole first-law residual. General many-particle events require all participant energies
and stoichiometry. These are conditional event identities, not a microscopic BNV model.
E_occ, mu_i, escaping energy and deposited energy use the **same energy zero and rest-mass
convention**. A change of energy zero must transform q_E and mu Gamma together.

The primitive semantic inputs are process identity/stoichiometry, event measure, and an
energy-partition model that closes event conservation. One may parameterize E_occ/E_esc and
derive E_dep, or supply deposition and solve for escape; one cannot treat two conservation-
related residuals as independent heat inputs. External incident energy is independent only
when there is an actual external owner, such as an incident particle or controlled reservoir.

## 6. Fermi-hole double-count adjudication

For an occupied state below the Fermi surface, the microscopic rearrangement contribution
would be `(mu_i-E_occ)R` under P8's assumptions. Adding independently defined product deposition
`E_dep R=(E_occ-E_esc)R` gives `(mu_i-E_esc)R`. Therefore the hole contribution is a **partial
microscopic representation of the already derived open-system chemical residual**, not an
additional source on top of P7. The thermodynamic chemical term `+mu_i R` is not by itself the
hole energy: the occupied-state removal energy must be subtracted too.

| Term or proposed usage | Classification | Exact rule |
|---|---|---|
| Independently incident/deposited external energy | INDEPENDENT | Add only with distinct external event/boundary owner |
| True product deposition E_dep in P8 | MODEL-DEPENDENT partition input | Add to `mu-E_occ`, not again to `mu-E_esc` |
| Hole `mu-E_occ` compared with full P7 | PARTIAL OVERLAP | It is one component of that residual |
| Hole plus product deposition compared with P7 | SAME TERM / DIFFERENT REPRESENTATION | Equality P8; one representation only |
| "Deposited energy" defined as `mu-E_esc` | SAME TERM / DIFFERENT REPRESENTATION | It already is the complete prompt residual |
| An added hole term on top of P7 | Forbidden double count | Unless a different state/reservoir contribution is proved and subtracted from its previous owner |
| Named Fermi-hole formula attributed to authenticated Zakeri corpus | SOURCE-LIMITED | No such named term/formula located in the four audited PDFs; no invented equation attribution |
| Specific absent primary's claimed heating efficiency | SOURCE-LIMITED | Must acquire/authenticate/audit that exact source before equation-level verdict |

**Formal double-count rule is resolved by P8. Published equation-level Fermi-hole attribution
is not resolved by the authenticated corpus.** No audited paper is accused of adding the same
term twice without an actual equation demonstrating it. A generic proof cannot establish
what an absent author's formula means. Full per-paper evidence is in section 17.

## 7. Local-to-infinity derivation and spatial sources

For a static emitter, `d tau=e^nu dt`. A freely propagating escaping product has conserved
Killing energy `E_infinity=e^nu E_local`; this holds for massive as well as massless freely
escaping particles. Counting events uses time dilation alone; power uses both factors:

```text
S_y=dot N_y^BNV=int_D e^nu Gamma_y(r,t) dV                         (P9)
P_direct^infinity=int_D e^(2nu) q_direct(r,t) dV                   (P10)
L_esc,BNV^infinity=int_D e^(2nu) sum_a E_esc,a(r,t)R_a(r,t) dV.     (P11)
```

If R is count cm^-3 proper-s^-1 and the radial measure is km^3, multiply by 10^15.
If density is fm^-3, the count integration uses 10^54. These conversions belong to the
consumer, not to a second geometry owner. The thermal identity follows alternatively from
`T_infinity dot S=int e^(2nu) T div(su)dV` in a Tolman-isothermal static domain.
FR05 PDF p3 eqs.(1)-(6), p4 eq.(15), and R06 PDF p2 eq.(6) agree with the independently
derived factors; `G_y`'s inverse lapse is a different operation.

E_esc means energy of products that really leave the accounted system, evaluated at the
stated production/accounting frame, after retention and reabsorption are resolved. A massive
product gravitationally bound to the star is not escaping merely because it is weakly
interacting. If it scatters, radiates, converts or escapes later, transport and storage own
those stages. Rotating emitters require `-p_t`, directional transport, and angular-momentum
accounting rather than blindly using the static formula. Spatially separated production and
deposition must each use their own radius/lapse or a conserved-Killing-energy transfer map.

A whole-star S_y suffices for frozen diffusive **eta drive**. It generally does not determine
thermal power because local mu, E_esc, deposition fraction, and event depth vary. It suffices
for a complete ledger only when the relevant redshifted conjugates and energy-per-event
weights are supplied as consistently integrated quantities, with declared support and prompt
relaxation. Multiple diffusive domains or unresolved transport require separate states/fluxes.

## 8. Mandatory standard-beta limiting case

BNV=0. Positive net beta direction is neutron decay,
`Gamma_n=-DeltaGamma_l`, `Gamma_p=Gamma_l=+DeltaGamma_l`. Spectators cancel. Thus

```text
-sum_i mu_i Gamma_i^beta=(mu_n-mu_p-mu_l)DeltaGamma_l
                       =eta_l DeltaGamma_l.                         (P12)
```

No external input does **not** set total q_E=0 while neutrinos escape. The beta matter-energy
source includes `q_E,beta=-Q_nu,beta`, so the full entropy/thermal source is
`sum eta_l DeltaGamma_l-Q_nu,beta`. With uniform `eta_l^infinity=e^nu eta_l`,

```text
R_l=int_D e^nu sum_(a in l) DeltaGamma_a dV
L_H,beta^infinity=C_M sum_l eta_l^infinity R_l.                       (P13)
```

This is exactly the governed eta R identity, not an alternate sign convention. Subtract only
the same-temperature, same-Ltilde equilibrium term already owned by ordinary cooling:

```text
DeltaL_nu,beta=sum_a Ltilde_a [F_a(xi_l)-1] T_infinity^q_a
DeltaP_beta=L_H,beta-DeltaL_nu,beta.                                 (P14)
```

The implementation confirms `eta*rate` and `h-d` respectively
(`CompactStar/Physics/Rotochemical/RotochemicalReactionResponse.hpp:35`, `:40`, `:53`;
`CompactStar/Physics/Rotochemical/FrozenRotochemicalRunContext.hpp:136`).
Small-xi modified Urca gives a negative quadratic coefficient:
`xi H_M-(F_M-1)=(-7340/(11513 pi^2))xi^2+O(xi^4)`.
Thus the incremental beta term **can cool**. The controlled historical trajectory demonstrates
this before its sign crossing
(`docs/validation/PHASE5D1_CONTROLLED_ROTOCHEMICAL_EVOLUTION_RATIFICATION.md:156`).
No sign, state, normalization, or redshift contradiction with Phase-5 was found.

## 9. Unreduced chemical response and source map

Neutrality eliminates proton count as an independent coordinate:

```text
y=(n_n,n_e,n_mu)^T, N_y=(N_n,N_e,N_mu)^T, n_p=n_e+n_mu
b=(1,1,1)^T, N_B=b^T N_y
g_y=(mu_n,mu_p+mu_e,mu_p+mu_mu)^T
L=[[-1,-1],[1,0],[0,1]].                                             (P15)
```

These follow from differentiating the neutral energy and ADR-0013's x/y transform, rather
than copying the prompt. The actual code constructs L at
`CompactStar/Analysis/src/ChemicalResponse.cpp:705` and `Z=L^T G^-1 L` at `:715`.
Species and beta axes are `CompactStar/Analysis/ChemicalResponse.hpp:16` and `:22`.

For a fixed metric, `delta n_y=C_y delta g_y(r)`, `delta g_y(r)=e^-nu delta g_y^infinity`.
Integrating counts gives

```text
G_y=int_D e^-nu C_y dV,  z:=delta N_y=G_y delta g_y^infinity
eta^infinity=-L^T delta g_y^infinity.                                (P16)
```

In P16 use dV in fm^3 when C_y is in fm^-3/MeV, or multiply a km^3
measure by 10^54 (cm^3 by 10^39); G then has units count/MeV.
G is symmetric positive definite only on qualified supported physical directions; no padded
inverse, gauge pseudoinverse, or unsupported species source is allowed. The frozen npe-mu
oracle has full three-dimensional global support even though some radial branches have fewer
active species. Production integrates the inverse-lapse weight at
`CompactStar/Analysis/src/ChemicalResponse.cpp:579`.

For actual whole-star particle counts with no material boundary flux,
`dot N_y=L R+S_BNV`. Reversible structural readjustment alone does not create particles.
For departures from a moving equilibrium reference,

```text
dot z=L R+S_BNV+S_boundary-dot N_y,ref.
F_ref := -dot N_y,ref                     (reference forcing, not creation).  (P17)
```

At fixed G and fixed L, differentiating P16 yields

```text
D_BNV := -L^T G_y^-1                  [2 by 3, MeV/count]
dot eta|_S=D_BNV S
D_BNV L=-Z;  Z=L^T G_y^-1 L.                                      (P18)
```

For fixed-baryon spin forcing `dot N_ref=L I_Omega 2Omega OmegaDot`, P18 gives
`+2Z I_Omega Omega OmegaDot=+2W Omega OmegaDot`, recovering Phase-5D.
D's rows are `(Npe,NpMu)` and columns `(Neutron,Electron,Muon)`. It is rectangular and has no
matrix symmetry claim. It should be a derived, qualified solve/view retaining the complete
G owner, basis, branch support, domain, metric, provider bytes, equilibrium anchor, numerical
uncertainty, lifetime and currency. Do not independently store it as another authority. It
requires distinct output-channel and input-number typing; the Phase-5C symmetric Z concession
explicitly does not authorize an asymmetric expanded map
(`docs/adr/ADR-0013-corrected-rotochemical-chemical-coefficients.md:590`).

### 9.1 Proof Z alone is insufficient

`b^T L=0`, `rank(L)=2`, hence `range(L)=ker(b^T)`. Neutron loss `S=(-R,0,0)` has
`b^T S=-R`, so it cannot equal Lr. More strongly, let A=G^-1. Replacing A by
`A+b v^T+v b^T` leaves `L^T A L` unchanged but changes D by `-(L^T v)b^T`.
For sufficiently small perturbations positive definiteness survives. Thus identical Z can
correspond to different BNV drives. No reduction from a generic source to the beta plane may
precede its unreduced mapping.

### 9.2 Neutron-removal sign oracle

For ideal free gas the neutral Hessian is
`H_y=diag(d_n,0,0)+[[0,0,0],[0,d_p+d_e,d_p],[0,d_p,d_p+d_mu]]`,
where `d_i=p_F,i^2/(3mu_i n_i)>0`. Neutron response is block-decoupled locally and after
integration. Removing neutrons gives `dot eta_e=dot eta_mu=-R/G_nn<0`. Independently, removing
occupied neutrons lowers mu_n before weak readjustment while other number densities remain
fixed (`CompactStar/EOS/src/LocalThermodynamics.cpp:125`,
`CompactStar/EOS/src/TrackRFreeGasThermodynamics.cpp:566`).

The existing governed Phase-5C baseline was used only as a mathematical oracle, not as runtime
source data or a BNV model. Its SHA-256 is
`7027aa6179fe9111d76586d467d007c4223ec9d26ece0fee773327973c77bfe7`; G is at
`tests/baselines/phase5c_chemical_coefficients.json:142566`. Independent scratch calculations give

```text
G_nn=2.904621518418346e55 count/MeV
D=[[ 3.442789339881122e-56, -4.544875287303884e-54, -4.828240976290691e-55],
   [ 3.442789339881122e-56, -4.828240976290691e-55, -1.0265285185799286e-52]]
dot eta_e/R=dot eta_mu/R=-3.442789339881122e-56 MeV/count.
```

No rate magnitude was selected. These values characterize this frozen oracle only; interacting
EOS cross derivatives need not preserve every source-direction sign.

### 9.3 Charge-changing channels

A neutral-chart source implies `S_p=S_e+S_mu`. Pure proton disappearance with electrons
unchanged does not lie in the chart and is **not** repaired by simply assigning a three-vector.
Charge must be conserved by the full channel, outgoing charged currents/products, surface
charge, and electromagnetic fields. Rapid EM redistribution does not alone determine whether
electrons or muons are removed; lepton conversion and transport need specified processes.

An explicitly neutral paired removal `(nu_n,nu_p,nu_e,nu_mu)=(0,-1,-1,0)` maps to
`nu_y=(0,-1,0)`. It is a permissible controlled aggregate event with direct residual
`(mu_p+mu_e+E_in-E_esc-J_X)R`; both eta derivatives are positive in the free-gas oracle.
It is not a statement that all proton BNV has this source. Charged outgoing products or
positron annihilation must include their compensating electron source and energy exactly once.
Unclosed charge currents or an unknown post-EM source are refusal/source-limited cases.

## 10. Chemical free energy is a state function

Within the fixed-geometry cold linear response define A=G^-1 and z=N-N_ref. Integrating the
local Hessian quadratic form gives

```text
E_2^infinity=(1/2)int e^nu delta n_y^T H_y delta n_y dV
            =(1/2) z^T A z.                                         (P19)
```

In P19, delta n has units fm^-3, H_y MeV fm^3, and dV fm^3; convert a
km^3 measure by 10^54 or cm^3 by 10^39 before integration.
Units are MeV (convert once to erg). E_2 is **tangent-subtracted quadratic thermodynamic
storage**, not the entire cold energy and not solely beta disequilibrium when baryon number
changes. At a beta-equilibrium reference `g_0^infinity=mu_B0^infinity b`,
`E_cold(N)-E_cold(N_ref)=g_0^T z+E_2+O(z^3)`.
The linear rest/chemical term cannot be omitted from a total-energy claim for baryon removal.
This fixed-metric functional is not the varying-star ADM mass functional.

Let `a=b^T G b`, `k=G b/a`, and `alpha=b^T z`. Then uniquely
`z=k alpha+L ell`, with `eta=-Z ell` and vanishing cross term because `A k=b/a`:

```text
E_2=alpha^2/(2a)+(1/2)eta^T Z^-1 eta
E_beta=(1/2)eta^T Z^-1 eta=V/2.                                    (P20)
```

The first term is reversible baryon-direction curvature; E_beta is energy relaxable by beta
reactions at the same baryon count within this model. Eta alone loses alpha, so future
complete chemical work/storage needs integrated baryon departure or full z in addition to eta.
This is required even if D already supplies the immediate eta derivative.

For frozen coefficients and reference forcing F_ref,

```text
dot E_2 = -eta^T R+delta g^T S_BNV+delta g^T F_ref
          +delta g^T S_boundary
(dot E_beta)_beta=-eta^T R <=0
(dot E_beta)_S=eta^T Z^-1 D S
             =(delta g-alpha b/a)^T S
(dot E_beta)_spin=eta^T Z^-1(2W Omega OmegaDot)
                 =2eta^T I_Omega Omega OmegaDot.                    (P21)
```

Nonincrease is not unconditional strict decrease: an inactive dissipative direction can
remain frozen, as ADR-0014 `:132` already qualifies. **No stored energy or its source-driven
growth is added as an independent heat term.** Its beta decrease is the source of eta R;
adding both the decrease and eta R doubles the same transfer.

## 11. Changing background, dot(G), dot(Z), and valid frozen scope

Within the linear model, at fixed L but G(t),

```text
dot eta=-L^T A(dot N-dot N_ref)+L^T A dot G A z
 dot E_2=z^T A(dot N-dot N_ref)-(1/2)z^T A dot G A z.                 (P22)
```

If L itself changes, add `-dot L^T A z`; if reference imbalance is nonzero, also differentiate
that reference. These are derivatives of the stated model, not permission to evolve an
unqualified background. For the fixed-baryon chart `eta=-Z delta N_l`,

```text
dot eta=+dot Z Z^-1 eta-Z(dot N_l-dot N_l,eq)
dot E_beta=eta^T Z^-1 dot eta
          -(1/2)eta^T Z^-1 dot Z Z^-1 eta.                          (P23)
```

The plus sign exactly matches ADR-0014 `:380`. Generic baryon-changing P22 cannot be replaced
by P23 alone. In particular the real hydrostatic sequence tangent `partial N_eq/partial N_B`
is not generally k: k is equilibrium at fixed geometry, not an evolving stellar sequence.

A first controlled fixture may freeze background/G/Z/W/Ltilde/thermal model only over a
perturbative domain: `abs(Delta N_B)/N_B <<1`, small local composition shifts within the
qualified branch, `abs(eta)<<mu_eq`, fast diffusive/EM relaxation relative to forcing, and
small omitted coefficient effects. Small global depletion is necessary, **not sufficient**:
localized depletion or a threshold can invalidate it. Quantify, for example,

```text
||G_0^(-1/2) Delta G G_0^(-1/2)|| <<1
|Delta C_*|/C_* <<1; weighted changes of W, Ltilde, lapse, g_surface and I small. (P24)
```

Near cancellation use absolute power/error budgets, not relative error against zero. No
numerical tolerance is selected here. Before setting one, measure same-model structural
sequence derivatives, response/thermal coefficient variation, source-depth dependence,
active-support movement, and induced errors in eta/T/power against predeclared science goals.
Long secular depletion requires updated structure, G, Z/W/Ltilde where applicable, heat
capacity, surface gravity, support, rotation and their provenance; it cannot omit P22/P23.

## 12. Matched no-BNV control and terminology

Construct trajectories B and C with the same initial star, T_infinity, eta, baryon departure,
spin history, standard beta microphysics and enabled channels, same Ltilde, photon/envelope,
standard neutrino/thermal models, solver and numerical tolerances. C sets **only the BNV
particle and associated energy sources** to zero. No unrelated equilibrium star substitutes.
In a future state-coupled model background/spin differences caused by BNV are outcomes under
the same model, not an excuse to choose different initial controls.

Write `Delta X=X_B-X_C` at the same infinity time. Define
`Delta T_BNV=T_B-T_C`, `Delta Lgamma_BNV=Lgamma(T_B)-Lgamma(T_C)`, and

```text
P_direct,BNV = int e^(2nu)q_direct dV
P_beta-mediated,BNV = Delta[L_H,beta-DeltaL_nu,beta]
Delta P_total,BNV = P_direct,BNV+P_beta-mediated,BNV
                   -Delta Lnu_eq-Delta Lnu_other-Delta Lgamma
                   +Delta P_other,owned.                           (P25)
```

P25 is the difference of complete thermal RHS powers. The equilibrium-neutrino feedback at
different temperatures does not disappear merely because microphysics is matched. On each
trajectory DeltaL_nu,beta means enhancement over its **own same-temperature** equilibrium;
outer `Delta[...]` instead means difference between trajectories. These operations differ.

For fixed common heat-capacity function define `U(T)=int^T C_*(T')dT'`. Then
`d[U(T_B)-U(T_C)]/dt=Delta P_total,BNV`. It is not generally `C_* d(Delta T)/dt`.
Report both signed `int P_direct dt` and `int Delta P_total dt`, with endpoints, energy units,
and control identity. Neither positive direct power nor positive beta chemical dissipation
alone establishes a positive total thermal effect. Photon response is a loss/observable, not
another independent heating mechanism.

## 13. Required thermal ledger and ownership

All global energies/powers below are at infinity, converted to erg/erg s^-1 once. Local
quantities use cm^-3 proper-s^-1. "Derived" means not an extra RHS contribution. `P4`, etc.
refer to derivations here; Phase-5 authority is ADR-0014 §§3.3-3.9. No new row is implemented.

| Row | Definition / sign / units | Local-global and redshift | Owner / independent or derived | Double-count partner | Authority / implementation status |
|---|---|---|---|---|---|
| A particle source | Gamma_i=sum nu_ia R_a; positive creation; count cm^-3 s^-1 or S_i count/s | P9, one lapse | Channel stoichiometry + source integration; primitive particle input | beta source, boundary current or reference forcing mislabeled as BNV | P6/P15-P18; proposed, no BNV API |
| B escaping BNV energy | E_esc>=0; L_esc>=0 outward, subtract from matter; erg/s | P11, two lapses | Escape/transport model; partition input | C, J_X, later product escape, beta neutrinos | P7/P11; B22/A24/G24/Z24 qualified audit; proposed |
| C direct product transfer | q_E includes destruction and returned deposition; signed erg cm^-3 s^-1 | P2/P7/P8; two lapses | Total matter-energy transfer; derived from closed partition or independent external transfer | B, D, hole, deposition already defined as residual | P3/P8; no new implementation |
| D open-system chemical term | -sum mu_i Gamma_i; either sign; local erg cm^-3 s^-1 | two lapses; -C_M (g_y^infinity)^T S if uniform conjugates | Local thermodynamics + particle source; derived | hole representation; E state change | P4/P8; proposed BNV use |
| E stored chemical energy | E_2 or explicitly E_beta; nonnegative quadratic state, derivative either sign; erg | P19 single energy lapse, P20/P21 global | Chemical state diagnostic; derived state, NOT heat | F and source-driven storage counted twice | ADR-0013 plus P19-P21; proposed diagnostic |
| F beta dissipation | L_H=C_M sum eta_l R_l>=0; erg/s | P13; two local lapses | Existing standard reaction response; derived transfer from E_beta | E_beta decrease; H | ADR-0014/P13; governed implementation reused unchanged |
| G neutrino enhancement | DeltaLnu=sum Ltilde(F-1)T^q>=0; erg/s | two local lapses already in Ltilde | Existing same-Ltilde response; derived loss | full equilibrium+enhanced neutrino term; H; direct-event neutrinos | ADR-0014/P14; governed implementation unchanged |
| H incremental beta effect | F-G; either sign; erg/s | global; no new lapse | Thermal boundary; derived subtotal | F/G must not also be added if H is added | ADR-0014/P14; governed implementation unchanged |
| I photon response | DeltaLgamma=B-C; either sign; subtract in P25; erg/s | surface luminosity already at infinity | Existing photon model evaluated at each T; derived feedback | direct deposition / total power subtotal | ADR-0014 §3.9/P25; existing owner unchanged |
| J structural dissipation | Only an authenticated irreversible mechanism; absent in controlled scope | mechanism-dependent; two lapses for local dissipative density | Separate dissipative physics owner, if justified | reversible PdV, gravity, rotational work | No source-backed mechanism established here; excluded |
| K total matched thermal power | P25, either sign; erg/s | already global | Comparison of complete RHS, derived | all component/subtotal rows | P25; proposed matched comparison |
| L integrated heat/cooling | int K dt, signed erg; also report int P_direct dt separately | infinity coordinate dt; no extra lapse | Diagnostic integral, derived | endpoint Delta U; not a second energy transfer | P25/P27; proposed oracle |

## 14. Structural and rotational work; global stellar first-law audit

Quasistatic PdV and gravitational readjustment are reversible work/state changes unless a
specific dissipative mechanism is provided. Neither a change in M nor a change in binding
energy is automatically heat. The BNV corpus discusses quasistatic structure, deposition and
neutrinos; section 17 does not establish a separate universally dissipated gravitational
heating term. No such term enters row J.

BNV products can carry angular momentum. A24 Appendix D and Z24 §§2-3 concern mass/rotation
and orbital timing; they are not thermal luminosity formulas. In a Newtonian illustrative
state identity, `E_rot=J^2/(2I)` implies `dE_rot=Omega dJ-(Omega^2/2)dI`; neither term is
implicitly dissipated. For relativistic rotation a source-qualified stellar first law is
required. Prescribed-spin v1 treats spin forcing as externally controlled; it does not
calculate BNV torque, I evolution or angular momentum carried by products.

### 14.1 What authenticated global authority supports

H67 PDF p15 / journal p1019, eqs.(83)-(85), supplies isentropic injection chemical potential
and its first integral; H67's metric symbol nu is twice this repository's nu, so its
`e^(nu/2)` becomes `e^nu`. H67 PDF p18 / journal p1022 eqs.(109)-(110) supplies baryon number
and `dn/n=d epsilon/(epsilon+p)`. H70 PDF pp3-4 / journal pp113-114 eqs.(10),(15),(16) derives
`delta M=int Omega delta j` holding baryon number and entropy of each fluid ring fixed.
For uniform Omega this is rotational work `Omega delta J` under those constraints. It is
not an independent heat source or a general finite-entropy multi-species theorem.

### 14.2 Derived static cold check

A independently varied the authenticated cold count/mass integrals; E checked the logic.
In geometric units let `f=1-2m/r`, `mu=(epsilon+p)/n`, `d epsilon=mu dn`, `dp=n dmu`.
TOV implies `mu'=-mu nu'`; therefore `mu_infinity=mu e^nu` is constant. With a vacuum-matched
zero-pressure surface,

```text
M=int_0^R 4pi r^2 epsilon dr;  N=int_0^R 4pi r^2 n f^(-1/2) dr
 delta N=int_0^R 4pi r^2 delta epsilon K(r) dr
K(r)=1/[mu(r)sqrt(f(r))]+int_r^R 4pi s n(s) f(s)^(-3/2) ds.
```

Using TOV, `d[1/(mu sqrt f)]/dr=4pi r n f^(-3/2)`, hence K'=0 and
`K(R)=1/mu_infinity`. Restoring units yields the restricted derived corollary

```text
delta(M c^2)=mu_B^infinity delta N_B.                               (P26)
```

For a moving surface the leftover term is `-4pi R^2 p_R delta R`; it vanishes only at p_R=0.
The repository finite p_cut surface requires tail/boundary work for an exact whole-star claim.
P26 relates cold equilibrium baryon loss to redshifted chemical energy, consistent with P8;
it does not make all that energy heat or set E_esc equal to mu without a channel model.

### 14.3 Source-limited generalization

The broad schematic `d(Mc^2)=T_infinity dS+Omega dJ+sum mu_i^infinity dN_i` is **not installed
as an authenticated theorem**. It requires specification of thermal/diffusive equilibrium,
which species are independent conserved quantities, electrochemical potentials, rotation,
boundary conditions and allowed variations. Hartle & Sharp (1967), Bardeen (1970), and
Thorne (1977) are cited/discoverable primary leads but absent from the authenticated catalog.
Their general theorem must be audited before a realistic changing-star energy closure.
The restricted P26 is a derivation from catalogued authority, explicitly not a quotation of
those absent works. A18 is unnecessary for this abstract controlled check.

## 15. Symbolic finite-interval conservation oracle

Use a closed **controlled fixture**, not a closed ordinary-matter system: the fixture ledger
includes escaping products and a prescribed forcing reservoir. Freeze G, L, g_0 and the
thermal coefficient model; represent parameter forcing explicitly by F in
`dot z=L R+S+F`. This affine reference is a manufactured controlled model; it is not a fully
self-consistent changing hydrostatic sequence. Let `delta g=A z`, `g=g_0+delta g`, and
`U(T)=int C_*(T)dT`. No separate X reservoir here. All terms use the same energy units.

```text
dot U = P_in-L_esc-g^T S+L_H-Lnu_full-Lgamma-Lother
 dot E_2=delta g^T(S+F)-L_H
=> dot(U+E_2)=P_in-L_esc-g_0^T S+delta g^T F
               -Lnu_full-Lgamma-Lother.                            (P27)
```

P27 cancels beta chemical dissipation internally and counts escape once. The input available
relative to the cold tangent is `P_avail=P_in-g_0^T S`; for disappearance `-g_0^T S` is
positive reference chemical energy released, not externally injected energy. Reporting only
P_in as "BNV injected energy" would omit the changing cold reservoir. Equivalently include the linear reservoir `g_0^T(N_y-N_y,initial)` using actual
particle numbers in the state and move its derivative to the left. Using z instead would
also require the reference-motion work, because F acts on z rather than actual counts.

Subtract the matched control (initial differences zero) and integrate [0,t]:

```text
int [P_in-g_0^T S+Delta(delta g^T F)] dt
 = int L_esc,BNV dt + Delta E_2(t) + Delta U(t)
   +int {Delta[DeltaLnu_beta]+Delta Lnu_eq+Delta Lother+Delta Lgamma} dt. (P28)
```

Common F does not cancel its work because delta g differs between trajectories. For the
strict no-spin/no-parameter-forcing oracle set F=0; it then closes with no forcing-work term.
For a separately stored X add Delta E_X and its terminal escape/thermalization exactly once.
For changing g_0/reference add their explicit potential work; for changing G include P22.
The finite-temperature composition-entropy correction in P5 must have its own consistent
state-function terms if retained. These additions cannot be hidden in a positive efficiency.

Frozen geometry cannot close the full ADM mass, gravitational boundary work, or state-coupled
rotational budget. P28 is an exact **symbolic identity of the declared controlled model**,
not a certification of a real star's global energy conservation. No conservation baseline
was generated. A future implementation must test this identity against independent integrated
states and emitted energies, including near-cancellation absolute residuals.

## 16. Future controlled channels and predeclared analytic oracles

| Future channel | Neutral source / expected drive | Prompt thermal residual | Required oracle / double-count rule |
|---|---|---|---|
| 1 neutron disappearance | S=(-R,0,0); both eta derivatives negative in free gas | `(mu_n+E_in-E_esc)R` | Independent dmu_n/dn_n and full-G solve; no extra hole |
| 2 neutral proton-electron paired removal | S=(0,-R,0), with nu_p=-1; both positive in free gas | `(mu_p+mu_e+E_in-E_esc)R` | Charge conservation and lepton ownership; distinct beta reaction stoichiometry; no automatic proton-only projection |
| 3 no-beta-drive aggregate | S=k dot N_B, k=Gb/(b^TGb); D S=0 | `P_in-Lesc-C_M g^T S` globally | Exact null solve; fractional aggregate, not single integer microscopic reaction; baryon curvature can change with eta=0 |

For the governed free-gas mathematical oracle,
`k=(0.9922178920029234,0.007484539490838682,0.0002975685062378332)` and `b^T k=1`.
No actual BNV rate, event energy value or lifetime is selected for any row.

| Oracle | Predeclared expectation |
|---|---|
| Zero BNV source | Same initial state and same governed stack reproduce exactly the Phase-5D trajectory; no extra rounding path should activate |
| Nonzero S, zero prompt residual | Chemical state may change; direct controlled thermal source is zero; later beta response may have either sign |
| Escape equals chemical event release | With no external/X input, `E_esc=-sum mu nu` yields zero P7 residual; conditional physical realizability must be checked |
| No escape, complete retention | `q_direct=-sum mu nu R` (plus independent input); for disappearance `mu_i R`; a thermalized upper endpoint for fixed channel/input, not a universal heat bound |
| No weak reactions, frozen coefficients | `eta(t)=eta(0)+D int S dt+int f_spin dt` |
| Weak reactions restored, forcing stopped | `dot E_beta=-eta R<=0`; strictness only on active dissipative directions |
| Null source | Dk=0 even though baryon count changes; retain alpha and reversible storage |
| Rank/conservation | b^T L=0, rank L=2, DL=-Z; baryon-changing source cannot be Lr |
| Energy ledger | P27/P28 close with escape, cold tangent, beta neutrinos, thermal state, photons and forcing work counted once |
| Energy-zero change | Shift mu and q_E consistently; residual unchanged |
| Charge violation | Unclosed charged source refuses before neutral projection |
| Nonthermal retention | Reservoir growth delays heat; later release not counted at production and again at deposition |
| Variable coefficient | P22/P23 terms required; frozen result recovered when derivatives vanish |
| Redshift | One lapse for counts, two for power, inverse lapse for G, single lapse for stored energy |
| Finite-T cancellation | Zero entropy residual does not establish zero exact T drive; P5 or an omission bound required |

These are proposed analytic tests only. They do not create governed numerical baselines.

## 17. Published-corpus audit and author clarification

### Three evidence classes; no inferred author intent

**A. PUBLISHED SOURCE CLAIMS** means only what the authenticated papers explicitly derive or
state. The table below exhausts the four catalogued papers' relevant BNV thermal statements;
it does not define a complete BNV thermal theory.

**B. AUTHOR WORKING NOTES / HYPOTHESES** are non-authoritative scientific inputs, separately
audited below when available. They cannot override the first law or become published claims.
The author clarified during this preflight that the published Zakeri papers do not contain
the author's complete intended BNV thermal theory, and that additional unpublished work is
unfinished. Absence from the papers is therefore no evidence of intent to exclude a term.

**C. INDEPENDENT PREFLIGHT DERIVATION** comprises P1-P28 and their stated assumptions, derived
from the first law and governed Phase-5 definitions. The objective is a complete, consistently
owned ledger for the declared scope, not merely reproducing the published treatment. Its
source-limited extensions remain explicitly open.

The canonical repository's document/text/TeX files, `notes.txt`, relevant legacy BNV headers,
and shared literature/data/external filename inventory were checked for thermal working
notes, Fermi-hole terminology and unpublished manuscripts. The only identifiable local
planning note is the architecture sketch `notes.txt`; no complete author thermal working
manuscript was located. This is a bounded project search, not a claim about private files or
other worktrees. The following author-provided qualification is retained verbatim:

> The author reports additional unpublished BNV thermal working notes that were not available/audited in this preflight. Therefore this document may assess the published corpus exhaustively, but it must not claim to exhaust the author's unpublished intended theory.

**AUTHOR WORKING NOTES — NON-AUTHORITATIVE / NOT PUBLISHED**

| Available local idea / evidence | Classification against independent derivation | Limit |
|---|---|---|
| `notes.txt:40`: BNV source drives d eta/dt and “heating” | Chemical drive independently derived/confirmed by P18; positive thermal sign unresolved and cannot be inferred from this label | Architecture sketch, authorship of detailed physical claim not established; not the unavailable thermal manuscript |
| `notes.txt:44`: chemical Gamma eta term | Independently derived/confirmed as beta chemical dissipation P13/P21 | Must also retain enhanced neutrino cooling; not total BNV response |
| `notes.txt:36`: BNV torque if supplied by a process | Plausible but not yet derived for a specific channel | Separate rotational reservoir; no heat equivalence |
| Legacy `CompactStar/Physics/BNV.hpp:61`: heating bounded by spin-down power, “model assumptions apply” | Source-limited as a general thermal bound | Cannot identify rotational work with heat without a model; not promoted to Phase-6 authority |

No unavailable idea is labeled inconsistent or double-counted. If the author's additional
notes become available, audit **each** idea as independently derived/confirmed, plausible but
not yet derived, overlapping/double-counted, inconsistent with the first-law ledger,
source-limited, or unresolved. That separate audit should precede independent Opus review /
owner ratification when the notes are available; their absence did not delay this derivation.

### Exact corpus and per-paper coverage

B22 is Berryman, Gardner, Zakeri, *Neutron Stars with Baryon Number Violation, Probing Dark
Sectors*, Symmetry 14, 518 (2022), DOI 10.3390/sym14030518, arXiv 2201.02637.
A24 is Allahverdi, Thompson, Zakeri, *Insights from Binary Pulsars and Laboratories into Baryon
Number Violation: Implications for GeV Dark Matter*, arXiv 2409.08178v1, 12 September 2024.
G24 is Gardner, Zakeri, *Probing Dark Sectors with Neutron Stars*, Universe 10, 67 (2024), DOI
10.3390/universe10020067, arXiv 2311.13649v2.
Z24 is Zakeri, *Pulsar Timing Anomalies: A Window into Baryon Number Violation*, arXiv
2311.05586v2, 6 May 2024; JCAP 05 (2024) 052 as identified by A24 reference 10.
Z24 pages below distinguish printed pages from PDF pages; the others coincide.

The omission assessment refers only to the published calculation compared with the
independently derived ledger. “Not addressed” makes no inference about author intent.

| Paper / process / equation-page | Ordinary particles destroyed/created | Initial-state energy | Escaping products | Deposited products | Named thermal term / Fermi-hole? / chemical potential? | Beta disequilibrium / enhanced neutrinos | Stellar readjustment | Double-count risk and assessment |
|---|---|---|---|---|---|---|---|---|
| B22: generic BNV; nn to e- e+ pp10-11; n to chi gamma Eq42 p19; n to 3chi Eq43 p19 | Channel dependent; full cascade must include annihilated ambient electron and surviving products | Cold constrained energy Eq4 p6; no occupied-state-to-heat formula | Dark particles and neutrinos; opacity Eq21 p10; p31 separates three neutrino origins | Photons/charged products qualitatively thermalize | Approximate retained-energy heating pp30-31; no named hole formula; delta mu=mu_n-mu_p-mu_e pp8-9 | Explicit direction/timescales Eqs12-18; no coupled governed eta/T ledger or explicit FR05 enhancement calculation | Eqs22-28 pp11-12; mass/rotation Eqs29-40 | Deposition and beta response recognized. Full quantitative partition not addressed; retained-energy estimates approximation-level. No published duplicate thermal addend established |
| G24: n to chi gamma/phi; chi chi to phi_B phi_B; induced decay, pp11-17,23 | Production and depletion chain determine source; cannot count first vertex alone | Local occupied E(p), Eq8 p17; CM energy Eq9 distinct | Light products after efficient depletion, pp15,23 | Detailed sequence required, explicitly p11 | Kinetic/annihilation heating in review context; no BNV hole formula; chemical equilibrium regime assumed | p12 distinguishes chemical-equilibrium and slower-chemical-response regimes; no evolved eta/F/H thermal model | Eqs2,10, table1; work explicitly p11 | Supports separate neutrino/thermal/work owners. Need for detailed partition is explicitly identified; quantitative partition is not calculated here; do not infer excluded physics |
| Z24: generic baryon loss, escaping chi or conversion back to SM; Eq3.7 printed6/PDF7 | Per-baryon loss times local n; subsequent sequence must be specified | Mean emitted chi energy in Q_chi=Gamma n Ebar_chi | Escaping chi/neutrinos carry angular momentum | No complete local deposition calculation | Emissivities used in angular-momentum estimate, not a named heat/hole formula; no mu-based BNV heat | Adjacent prose explicitly mentions BNV-induced beta enhancement; Eq3.8 uses equilibrium MU T^8 approximation | Evolving I/R/spin and equilibrium sequences | Emission is an energy-loss owner, not heat. Thermal closure not addressed; equilibrium emissivity is approximation-level. No duplicate thermal terms established |
| A24: B to psi gamma (B=n,Lambda), then psi n to pi- K+; alternate n to psi pi0, pp5-9 | Complete principal chain removes two ordinary baryons; charged products require full cascade | Eq25 E_B*=sqrt(m*^2+p^2), E_B=E_B*+Sigma0; its symbol mu is a mass ratio | Heavy psi bound/depleted; final neutrinos described; Appendix E escape calculation | Cascade energy partition not completed | Explicit distinction from heating studies p7; no own hole/heat formula | Urca-fast regime assumed; no eta evolution/enhanced-neutrino thermal calculation | EqD5 p20 sequence mass/rotation | Particle rate independent of energy model. Missing thermal closure not addressed, equilibrium regime approximation-level. Printed escape conventions require clarification, below |

### Equation-level crosscheck

These are published inclusions and omissions relative to class C, with omission status
explicit. Ledger classifications do not accuse an unpublished theory of any omission.

| Published locus | What is actually included | Relation to derived ledger; omission assessment |
|---|---|---|
| B22 Eq4 p6 | Cold energy with baryon/charge constraints gives mu_i=B_i mu_n-Q_i mu_e | Chemical-work authority; SAME TERM / DIFFERENT REPRESENTATION only when used within the first law. Independent thermal conversion not addressed |
| B22 Eqs12-18 pp8-9 | Direction-specific Urca relaxation estimates; phase-space functions nonzero at zero imbalance | Forward-minus-reverse required to obtain governed signed net R. Approximation-level, not a contradiction with Phase-5; complete coupled eta/T not addressed |
| B22 Eq21 p10 | Neutrino mean free path, energy/density dependent | MODEL-DEPENDENT escape transport; neither heat nor universal escape efficiency; complete transport not addressed |
| B22 pp10-11, nn to e- e+ | Positron annihilation and photon deposition qualitatively | MODEL-DEPENDENT complete stoichiometry/product partition. Ambient electron ownership required; exact thermal residual not addressed |
| B22 Eqs22-28 pp11-12 | Quasistatic equilibrium-sequence changes after baryon loss | INDEPENDENT background response; no dissipative heat mechanism specified. Irreversible thermalization not addressed |
| B22 Eqs32-33 p16; G24 Eq10 p21; A24 EqD5 p20 | Sequence mass/rotational energy, e.g. Mdot_eff=(M_,Ec+Omega^2 I_,Ec/2) Bdot/B_,Ec | INDEPENDENT mechanical reservoir/work, not thermal power. Installing it as heat risks duplicating other energy owners unless a dissipation model proves otherwise; no such model derived there |
| B22 section5.4 pp30-31, unnumbered thermal estimates | Observational heating scales from an assumed large retained energy per destroyed neutron; high temperature if most/all energy is trapped | Controlled approximation / MODEL-DEPENDENT. Reproduced only with the same retained residual assumption; occupied-state, redshift and chemical partition not calculated. No numerical rate chosen here |
| B22 p31, n to 3chi | States that invisible disappearance evades quoted heating constraints | Invisibility alone does not prove mu_n-E_esc=0. MODEL-DEPENDENT / SOURCE-LIMITED pending absent primary [164]; exact residual not addressed here, not a demonstrated model error |
| B22 p31 three neutrino origins | Direct BNV, restoration of chemical equilibrium, hotter-star thermal emission | Distinct owners. PARTIAL OVERLAP risk if a cascade is counted again as beta neutrinos. A closed L_H-DeltaLnu trajectory is not given; beta response is explicitly recognized |
| G24 pp11-12 section2.3 | Neutrino emission, thermal energy and work; detailed reaction sequence needed; chemical/hydro time hierarchy | MODEL-DEPENDENT partition; need for detailed examination explicitly identified, quantitative partition not calculated here. Supports complete ledger; no intention to exclude beta/structure may be inferred |
| G24 Eqs8-9 p17 | Occupied local E(p)=sqrt(m*^2+p^2)+Sigma0 and distinct CM energy | Energy-convention translation. Neither general occupied E nor CM energy is automatically mu; heat calculation not addressed |
| Z24 Eqs3.2-3.6 printed4-5/PDF5-6 | Rotational-energy and magnetodipole relations | INDEPENDENT rotational work. A symbol E_BNV in this derivation is not a thermal source |
| Z24 Eq3.7 and text printed6/PDF7 | Q_chi=Gamma n Ebar_chi for angular-momentum transport | INDEPENDENT escaping-energy owner if truly escaping; contributes negative q_E. Chemical residual needed for thermal interpretation, not calculated there |
| Z24 Eq3.8 printed6/PDF7 | Equilibrium neutron-branch modified-Urca Q_nu proportional to T_9^8 for small angular-loss estimate | Approximation-level equilibrium emissivity, despite explicit nearby enhancement caveat; no F(eta/kT). Signed denominator in the smallness ratio requires magnitude interpretation, never positive heating efficiency |
| A24 Eqs24-25 pp8-9 | dn_B/dtau proper density-loss rate and kinematic mass-ratio variables | INDEPENDENT particle-source model, not thermal input; exact product-energy partition not addressed. mu=m_psi/m_B* is not chemical potential |
| A24 EqD3 p19 | Bdot=-int e^nu Gamma_nm n dV, Gamma_nm proper per-baryon rate | EXACTLY REPRODUCED one-lapse number integral P9 |
| A24 EqE1 p20 | Prints E_psi^NS=e^nu E_psi^nm > E_esc^NS=m_psi e^-nu | Printed frame discrepancy: infinity threshold is m_psi, local threshold m_psi e^-nu. Requires convention resolution; not evidence underlying numerical calculation is wrong |
| A24 EqE2 p20 | Boosted local two-body product energy using total local occupied E_B | INDEPENDENT local product-energy input, not heat |
| A24 EqE3 p20 | Angular escape condition prints actual E_psi^nm in numerator | Solved threshold requires escape-threshold energy; actual energy makes it self-referential. SOURCE NOTATION AMBIGUITY; clarify before reuse |
| A24 EqE4 p21 | Escape rate called per baryon; prefactor epsilon_Bpsi^2 (m_B*)^2 times dimensionless integrals | Natural units energy^4, a density rate as Eq24, without explicit 1/n_B. Printed normalization discrepancy requires clarification; no inference about numerical plots |
| A24 EqE5 p21 | Ratio of integrals e^nu Gamma_esc n dV / e^nu Gamma n dV | Correct NUMBER escape fraction if both Gamma truly per baryon; conditional on E4 normalization. Not an energy escape fraction/luminosity |
| A24 Eq28 and adjacent p12 estimate | Visible (m_psi-m_p)/2 used to approximate cosmological 21-cm bound | Approximation in a different system, not neutron-star thermal deposition |

The Appendix E findings were checked in rendered PDF equations, not extraction alone.
No rate, escape fraction or published numerical constraint is implemented or recalculated.
Other heating/cooling discussion (B22 sections6.3/6.5; G24 sections2.1/2.2) concerns review
context: dark-boson cooling, incoming dark matter kinetic energy and annihilation. It supplies
no further BNV hole entropy formula. Incoming energy would be a separate external owner.

The exhaustive four-paper search found no standalone formula named Fermi-hole,
Fermi-sea-rearrangement or degeneracy heating; literal “hole” occurrences concern black holes.
This establishes a published-corpus scope fact only. It **does not establish author intent**,
absence from unpublished work, or that a future complete theory should omit rearrangement.
No paper in this corpus is shown to add duplicate hole and chemical-residual heating.
The prospective no-double-count rule in section6 is independently resolved; attribution of a
named microscopic term to additional primary literature remains SOURCE-LIMITED.

## 18. Proposed architecture and candidate invariants

All layers below are **PROPOSED, NOT IMPLEMENTED**. Preserve the Phase-5D standard subsystem
and frozen-context provenance rather than adding process-specific conditions to generic
execution (`CompactStar/Physics/Rotochemical/FrozenRotochemicalRunContext.hpp:136`, ADR-0014).

| Layer | Proposed semantic owner / provenance |
|---|---|
| 1 process identity | Directed channel ID, full ordinary stoichiometry, charge/current closure, participant and intermediate identities |
| 2 local particle source | R_a(r,t), Gamma_i=nu_ia R_a, proper-time/volume units; no hidden heat efficiency |
| 3 energy partition | E_occ distribution, E_in, E_esc, actual product deposition, nonthermal reservoir and thermalization model; event conservation |
| 4 global integration | One-lapse particle source, two-lapse power, explicit domain/profile and support; independent energy-weighted integration |
| 5 unreduced chemical response | Authenticated G_y and D=-L^T G_y^-1 derived on demand; ordering/units/conditioning, same star/context/hash; retain baryon coordinate |
| 6 storage diagnostics | E_2, alpha, E_beta and derivatives; cold tangent reference and parameter work |
| 7 standard beta coupling | Reuse governed signed rates, eta R and same-Ltilde enhancement; no BNV-specific branches in generic beta driver |
| 8 thermal ledger | q_E and chemical residual once; finite-T scope, nonthermal delay, distinct beta and escape owners |
| 9 matched control | Same setup with only BNV sources disabled; store all difference terms at equal times |
| 10 diagnostics | Residual conservation, source/energy provenance, sign-map outputs and uncertainty bounds |

The same particle-source model must support alternative energy-partition hypotheses without
changing chemistry. D may be cached only as a derived view with exact G/context provenance,
never as an independently fitted or separately governed physical coefficient. The existing
Phase-5D producer/manifest/diagnostic pattern can be reused after review; no present API is
claimed to satisfy this proposed BNV contract.

Candidate invariants, proposed here and in ADR-0015, **not added to the governed invariant file**:

| ID | Proposed invariant |
|---|---|
| BNV-1 | Particle source and matter-energy/escape/deposition source are distinct semantic objects linked by event identity and conservation |
| BNV-2 | Generic sources enter the unreduced charge-neutral G_y basis before beta reduction; retain baryon-changing direction |
| BNV-3 | Stored chemical energy is a state function, never an independent thermal addend; separate baryon curvature from beta-relaxable storage |
| BNV-4 | Beta-mediated matched response is Delta[L_H-DeltaLnu], with either sign; total response also includes cooling feedback |
| BNV-5 | Escaping energy is subtracted once, with declared frame, normalization, transport boundary and correct lapse |
| BNV-6 | A named hole term may be installed only after translating its occupied-state and deposition conventions and proving no overlap with the existing residual |
| BNV-7 | Net thermal sign is computed with an error bound, never assumed from “heating” terminology |
| BNV-8 | Matched no-BNV control differs only by BNV source objects; shared-input differences require explicit physical necessity and ownership |
| BNV-9 | Charge-changing channels close ordinary/product/field currents before neutral-chart mapping; no hidden charged disappearance |
| BNV-10 | q_E is total matter-energy transfer; q_E-mu Gamma is entropy energy, and temperature-RHS identification declares composition-entropy approximation |
| BNV-11 | Reversible structural/rotational work is not heat without a source-backed dissipative mechanism |
| BNV-12 | Frozen-background validity requires measured coefficient/residual variation bounds; changing G/reference/Z terms cannot be silently omitted outside that scope |
| BNV-13 | Published claims, non-authoritative author hypotheses and independent derivations remain separately labeled; omission never establishes author intent |

## 19. Future thermal-sign map, without a numerical scan

Predeclare separate signs for prompt P_direct, matched beta-mediated response, and total
matched instantaneous power; also report Delta T and Delta Lgamma, which have memory and need
not track the instantaneous sign. Candidate axes are event energy partition (escape,
deposition, delayed retention), charge-consistent species direction, xi_e=eta_e/kT,
xi_mu=eta_mu/kT, temperature, source depth/radius, and weak-process/support regime.
Initial state, elapsed time, spin forcing and matched-control identity are conditional inputs.
No BNV lifetime/coupling/rate amplitude or numerical threshold is selected here.

Outputs: **NET COOLING**, **NET HEATING**, **NEAR ZERO** only when the absolute residual and
its uncertainty support that classification. Overlay **SIGN UNRESOLVED** if finite-T,
energy-partition, source or numerical uncertainty spans zero. A small fractional error on
large individually cancelling terms cannot establish the residual sign. Rates affect the
state history, whereas energy partition changes the prompt residual at fixed source;
the map must keep those dependencies inspectable. This is a design, not computed domains.

## 20. Missing authority and scientific disposition

| Missing primary / item | Published pointer or required use | Consequence |
|---|---|---|
| Berryman/Gardner/Zakeri, PRD109023021 (2024), arXiv2305.13377 | G24[37], A24[9], Z24[14]; microscopic/macroscopic translation | Not the supplied B22 paper; authenticate before extending its specific rate/energy claims |
| McKeen/Pospelov/Raj, PRL127061805 (2021), Neutron Star Internal Heating Constraints on Mirror Matter | B22[185] p20 | Most direct identified heating-specific primary for named-hole interpretation; absent |
| McKeen/Pospelov/Raj, PRD103115002 (2021), Cosmological and astrophysical probes of dark baryons | B22[184] | Additional dark-baryon thermal interpretation absent |
| Ema/McGehee/Pospelov/Ray arXiv2405.18472; Fox/Hostert/Menzo/Pospelov/Zupan arXiv2407.03450 | A24[56],[57], explicit heating references p7 | Metadata located, no authenticated equation audit; not authority here |
| Strumia arXiv2112.09111 / JHEP02(2022)067 | B22[164], G24[92] | Specific n to 3chi retention/heating assertion source-limited |
| Goldman/Nussinov JHEP08(2010)091, arXiv0907.1555; Goldman/Mohapatra/Nussinov PRD100123021(2019), arXiv1901.07077 | B22/G24/Z24/A24 timing and mirror-loss references | Detailed energy-to-timing interpretation beyond audited sources remains source-limited |
| Haensel A&A262131 (1992); Baym/Pethick/Sutherland opacity reference | B22[115],[118] | Underlying directional-rate/opacity approximation details absent; no substitution for governed beta model |
| General applicable relativistic stellar first-law primary | H67/H70 support restricted results, not complete finite-entropy multispecies rotating law | General global closure source-limited; static cold derivation P26 remains available |
| Author's additional unpublished thermal working notes | Author clarification during this task | Not available/audited; audit separately if supplied before review/ratification; does not block first-principles derivation |

No new source was acquired, installed in `_incoming`, or promoted to catalog authority.
Additional acquisition must preserve provenance in `_incoming` and pass catalog/governance
before being treated as authority. Metadata discovery is not that approval.

**Chosen disposition D:**

**PHASE-6 BNV ENERGY LEDGER SOURCE-LIMITED — ADDITIONAL PRIMARY FIRST-LAW / BNV AUTHORITY REQUIRED.**

This disposition concerns the requested complete source-backed preflight, not failure of the
controlled algebra. The entropy balance, redshift factors, beta limit, full-G map, neutron sign,
storage decomposition and prospective no-double-count rule are derived and internally checked.
There is no established contradiction with Phase-5 and no demonstrated published duplicate heat
term. The unavailable named-hole primary and general global first-law authority prevent claiming
complete source closure. The author's unpublished intentions are not judged by these gaps.

**Blocking for disposition A:** authenticate and audit the missing thermal/hole primary needed
for exact named-term adjudication; obtain applicable general stellar-first-law authority or
explicitly review a narrower global scope. Resolve any Appendix E energy/normalization
convention before adopting that escape model. These do not prevent the abstract local ledger.
**Material:** finite-T composition entropy near cancellation; nonthermal/dark retention and
cascade transport; charge closure; changing structure/G/reference/spin and boundary work;
source-specific energy partitions. **Nonblocking for the abstract frozen oracle:** realistic
A18/FR2005 remains blocked, global INV-11 unresolved, no realistic rate, no numerical sign map,
and D shares A's identity. No genuine source-backed disagreement requires Fable.

**Recommended next action (not executed):** acquire with provenance and authenticate the
identified primary BNV heating/hole and applicable stellar-first-law references; resolve the
printed Appendix E frame and rate-normalization ambiguities before reuse; complete the
source-specific crosscheck against this proposed ledger. Audit the author's unpublished notes
separately if they become available, before independent Opus review / owner ratification.
Then request independent scientific review of the completed proposed contract. Do not implement
a BNV rate, A18 reconstruction or production BNV trajectory before that review and authorization.

## 21. Scratch verification and protected standard stack

All mathematical checks and build products were kept outside the repository in temporary
scratch directories. They are not governed baselines. C independently derived the free-gas
signs and full-space storage decomposition, with exact rational identity checks; E separately
used the authenticated governed G_y numeric payload. A/D checked the thermodynamic and finite-
interval cancellations independently, including nine exact-rational manufactured cases after
an unavailable SymPy import was replaced with standard-library arithmetic.

E's double-precision checks found DL+Z exactly zero in its chosen evaluation, relative Z-to-
governed-payload difference 1.156e-17, normalized Dk residual 4.80e-20, and relative quadratic
storage-split residual 1.226e-16. C's independent reconstruction agreed to roundoff. These are
scratch algebra checks, not precision claims for new physics. One-lapse count versus two-lapse
power and inverse-lapse susceptibility were independently dimension-checked.

A fresh scratch Debug build used the existing standard targets. Initial selected CTest run
passed 4/5; `phase5d_independent_oracles` could not import mpmath under the initially selected
system Python. Reconfiguration to the already installed miniforge Python, with mpmath, required
no repository/dependency change. The final run passed **5/5**, exit code 0:

- phase5d_response
- phase5d_independent_oracles
- phase5d_component_tolerances
- phase5d_protected_manifest
- phase5d_harness_controls

The protected manifest covers the 33 governed source paths, ten predecessor baselines and
three special files; the eleventh Phase-5D baseline was independently hashed. Entry hashes of
all 2896 tracked paths were compared successfully: zero changed entry paths. No multi-hour
Phase-5D trajectory suite was rerun; no governed code changed. Final document checks passed the exact two-document allowlist,
valid repository path/line references and companion links, PROPOSED ADR status, unchanged entry
hashes, `git diff --check` and `git diff --cached --check` before commit. The branch is committed
with the requested message and pushed non-force only after these checks. Its resulting commit SHA and
local/upstream/live equality belong in the final execution report, avoiding a self-referential
SHA inside this commit. Canonical master remains the entry SHA; no merge is authorized.
