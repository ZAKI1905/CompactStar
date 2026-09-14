# Phase-6A-0 BNV thermal first-law contract preflight — R3 review reconciliation

**Status: FINAL-REVIEW CORRECTIONS APPLIED / READY FOR BOUNDED INDEPENDENT RE-REVIEW.**
**ADR-0015: PROPOSED / NOT ACCEPTED / NOT OWNER-RATIFIED / NOT CANONICALLY INTEGRATED.**
**Date:** 2026-09-14 (R3); R1 derivation record dated 2026-09-13. **Change class:** documentation of proposed scientific-semantic and
architecture contracts; no numerical behavior or implementation changes.
**Canonical master:** `0a7418aecb7314cfa472a78f1faf477be8456a94`.
**PHASE6A0_DRAFT_ENTRY_SHA:** `5a6bf7cb9455d684ddb6fccb22ad2b9fec940b3a`.
**PHASE6A0_R2_REVIEWED_SHA (R3 entry):** `7a31862f1e8a5cf046316882e05d76e4924e27d9`.
**Branch:** `analysis/phase6a0-bnv-thermal-first-law-preflight`.
**Worktree:** `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase6a0-bnv-thermal-preflight`.
**Companions:** [PROPOSED ADR-0015](../adr/ADR-0015-bnv-open-system-thermal-ledger.md) and
[non-governed Cowling diagnostic](PHASE6A0_COWLING_BARYON_DIRECTION_DIAGNOSTIC.md).

## 1. Entry, evidence hierarchy and rewrite boundary

At R1 entry, local master, origin/master and live master were authenticated at the canonical SHA above;
local writer HEAD, upstream and live writer branch matched the draft entry. The worktree was
clean. Existing worktrees and histories touching these paths were inspected; the earlier draft
is the only differing ancestor, not a competing authority. The first live check encountered
sandbox DNS restriction; the permitted network escalation succeeded. No master mutation occurs.

E/root was the sole R1 repository writer. Read-only A/B/C specialists checked thermodynamics,
source/author-note interpretation and response algebra; their work is internal verification,
not the final independent Opus review. The owner-supplied five-agent hardening is separate
external review evidence, not these agents' output. No Fable adjudication is needed.

R3 is a documentation-only sole-writer correction pass by root/E; no new specialist review is
claimed. The clean worktree, local HEAD/upstream/live branch matched PHASE6A0_R2_REVIEWED_SHA;
local/origin/live master matched the canonical SHA. The R2 final report (158 lines) and
Reviewer E (360 lines) were read completely before edits. Reviewer E's adjudicated E-1–E-27
list is the exact correction checklist; section 21 records every disposition and source hash.
The reviewed equations and numerical oracles are preserved; R3 reconciles contract scope,
notation and declared approximation/observable limits without implementation.

Project authority remains GOVERNANCE (`GOVERNANCE.md:14`), ADR-0011, ADR-0013, ADR-0014 and the
scoped invariants. Phase-5D controlled frozen-v1 is closed and canonically integrated; global
INV-11 remains unresolved for broader physics (`docs/SCIENTIFIC_INVARIANTS.md:1004`). Phase-5C
is correct and governed for its declared R2006/Cowling contract. This revision does not replace
its coefficients or reopen Phase 5. The previous preflight is **superseded as a draft**, not
an authority to defend. Its fixed-reference baryon-response and storage construction are retired.

Evidence classes:

| Class | Meaning and limit |
|---|---|
| Project authority | Accepted Phase-5 ADRs/invariants and authenticated implementation definitions |
| Catalogued published source | Role-qualified literature authority; individual hashes checked |
| Uncatalogued local published source | Read and byte-fingerprinted supporting evidence only; `_incoming` remains unpromoted |
| Author working notes | UNPUBLISHED / NON-AUTHORITATIVE; useful hypotheses checked against derivation |
| Earlier Opus consultation | Working derivation; superseded by consolidated hardening where they disagree |
| Consolidated hardening R0 | Independent scientific evidence supplied by owner, not a governed contract |
| R1 derivation/check | Explicit algebra or scratch result here; proposal, not ratification |

The published Zakeri papers do not contain the author's complete intended thermal theory.
Omissions do not establish author intent. The previously unavailable Thermo_BNV and whiteboards
are now audited separately in section 18; none is silently corrected in place.

### 1.1 Input manifests

Shared literature root: `/Users/keeper/Documents/CompactStar/literature` (read-only).
Catalog SHA-256 `285040a9931bd3be20eedcc69900a0afd520ee702844c2d7771e140de12223f4`.
Checksum-manifest SHA-256 `567f21e661bc16dd4ca0dc188ce3188d2e1f8167077d02b2cb7b3a6327181d9d`.
All 22 listed entries passed. All 43 current literature files, including `_incoming`, were
fingerprinted at R1 entry for immutability; a fingerprint is not catalog promotion.

R0 means `/Users/keeper/Downloads/PHASE6A0R0_CONSOLIDATED_REPORT.md`, SHA-256
`c1df790a336a0b46aeb37154818a09ff44c5bcdb366f61697b3b6d5f00d41e92`, read completely (252 lines).
Its Q53 correction list and Q54 G1-G20 are the explicit coverage checklist below. R0's statement
calling commit 5a6bf7c an "uncommitted draft" is editorial: Git authenticates it as committed.
Its scratch scripts are described in R0, not included in the supplied report; new G_true TOV
numbers below are **reported R0 evidence**, not falsely claimed to be regenerated in R1.
Earlier consultation `/Users/keeper/Downloads/BNV_thermal_consultation.md`, SHA-256
`accfe7b20a198518814a980fc4a4db9209518a1b91a63593bcc3c9ccfe745b8a`, was also read completely.
It supplies history, not the final spin sign, source mapping, or QSS conclusions.

Catalog source paths below are relative to the shared root. R1 retains the earlier complete
four-paper/145-page audit and rechecks the affected interpretations; it does not claim a new
full rereading of every catalogued PDF. Published equation locators appear in section 17.
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


| ID | Exact additional input | SHA-256 | Evidence status |
|---|---|---|---|
| Thorne | `/Users/keeper/Documents/CompactStar/literature/_incoming/Thorne-Lecture GR Stellar Struc and Dyn.pdf` | `b5d20da0f9dbacfad2f2d05820d506adeacd343947b80571346090be0f1a6e51` | Uncatalogued published lecture; supporting only |
| MPR21 | `/Users/keeper/Documents/CompactStar/literature/_incoming/McKeen_2105.09951v2.pdf` | `213b15826f8753832437fab24db438f7d4e4e300cc94eca70ec88e51cbff9338` | Uncatalogued published preprint; supporting only |
| G22 | `/Users/keeper/Documents/CompactStar/literature/_incoming/Goldman_2208.03771v1.pdf` | `30c3713f1779943d326a252f886df5de5e1c721176cce196110ead2c4ae8bc0a` | Uncatalogued published preprint; supporting only |
| GR10 | `/Users/keeper/Documents/CompactStar/literature/_incoming/2010-Internal Heating of Old Neutron Stars- Contrasting Different Mechanisms.pdf` | `6db48d156bd140977b9f672d33a57f643c082d99660e61da4ebee790ceaa1399` | Uncatalogued published source; supporting only |
| Thermo_BNV | `/Users/keeper/Downloads/Thermo_BNV.pdf` | `8750d3ba5c36c4bd95a5b2a04b6df8279976784f418ae3369dcf366b8119ec92` | Author working notes, unpublished |
| WB-Oct | `/Users/keeper/Documents/CompactStar/literature/_incoming/Notes 20231009_173916.jpg` | `2a1c9f160593aa9505177fda20dd1b7ada6808f80785e60d112d762ba578f23a` | Author whiteboard, unpublished |
| WB-Dec | `/Users/keeper/Documents/CompactStar/literature/_incoming/20231219_190420.jpg` | `0289482d28becc5ddb0c3723655c61b227d26b0dec5d7d8d072d1dd7f98a02d9` | Author whiteboard, unpublished |

No new acquisition or catalog mutation is part of R1. Missing primary papers are not represented
as read. Derived equations do not acquire missing-paper authority merely because R0 cites them.

## 2. Notation, domain and units

Metric signature (-,+,+,+), lapse exp(Phi)=exp(nu_repo),
`ds^2=-e^(2Phi)c^2 dt^2+e^(2Lambda)dr^2+r^2dOmega^2`,
`dV_proper=4pi r^2 e^Lambda dr`. Hartle's nu_H=2Phi. Use c=1 in contractions and restore Mc^2
for stellar energy (`docs/SCIENTIFIC_INVARIANTS.md:227`). D=u.grad is proper-time derivative;
dots are infinity-coordinate time. rho is total energy density, s entropy density, s_b=s/n_B.
All mu and particle energies include rest mass with one common energy zero. n and Gamma use
count cm^-3 and count cm^-3 proper-s^-1; q uses erg cm^-3 proper-s^-1.

Global neutral coordinates are `N_y=(N_n,N_e,N_mu)^T`, `N_p=N_e+N_mu`,
`g_y=(mu_n,mu_p+mu_e,mu_p+mu_mu)^T`, `b=(1,1,1)^T`, `B=b^T N_y`;
`L=[[-1,-1],[1,0],[0,1]]`, `P=[[0,1,0],[0,0,1]]`, `PL=I_2`, `b^TL=0`.
Eta=-L^T g_y^infinity, order npe then npmu, is measured in MeV; positive net R creates leptons
(neutron decay), R in count/s. Reserve R for this global beta-rate vector; local directed
BNV event-rate density is R_a in count cm^-3 proper-s^-1. These definitions are ADR-0013 sections 3-4 and
`CompactStar/Analysis/ChemicalResponse.hpp:16`, `CompactStar/Analysis/src/ChemicalResponse.cpp:705`.

G_y is count/MeV, Z is MeV/count, t and k are count/count, sigma count/s. Chemical energies
are MeV; multiply exactly once by C_(MeV->erg) for erg. Its implementation owner is
`MeVToErg` (`CompactStar/Physics/Rotochemical/ChemicalImbalanceState.hpp:14`); the exactly-once
thermal boundary is `CompactStar/Physics/Rotochemical/RotochemicalReactionResponse.hpp:48`
and `:53`. Local thermodynamic R1–R5 use mu/rho/q in consistent erg units; R7 and subsequent
MeV chemical calculations show conversion explicitly when combined with erg/s luminosities.
E_chem and the equilibrium function in section 8 are in MeV; thermal section 9 uses
E_eq=M_eq c^2 in erg and C_(MeV->erg)E_chem for its chemical reservoir. Global powers/luminosities
are at infinity even when their superscript is suppressed. T_infinity
is Kelvin and k_B is MeV/K in xi. I_Omega denotes the **two lepton components** of the
fixed-baryon spin response; its three-component lift is L I_Omega. This avoids confusing
three-axis target motion with two-axis W=Z I_Omega.

The first toy is a whole-star, diffusive, nonrotating, non-superfluid free-gas fixture. Local
thermalization and a single ordinary-fluid velocity are assumed. Full rotating/time-dependent
metrics, finite-temperature structural evolution and multi-fluid transport are future scope.

## 3. G1, G15: local first law, worldline entropy and finite temperature

With `T_m^{mu nu}=(rho+p)u^mu u^nu+p g^{mu nu}`, define

```text
nabla_mu(n_i u^mu)=Gamma_i                    (creation positive)
nabla_mu T_m^{mu nu}=Q^nu,  q_E=-u_nu Q^nu     (total energy INTO ordinary fluid)
D rho+(rho+p)theta=q_E,   theta=nabla_mu u^mu.                  (R1)
```

The contraction uses u.u=-1 and u.Du=0. Gibbs `d rho=T ds+sum mu_i dn_i`, Euler
`rho+p=Ts+sum mu_i n_i`, and `D n_i=Gamma_i-n_i theta` yield

```text
T(Ds+s theta)=T nabla_mu(su^mu)=q_E-sum mu_i Gamma_i.           (R2)
```

Explicit cancellation: `-theta sum mu_i n_i+(rho+p)theta=Ts theta`. Thus reversible PdV,
hydrostatic/gravitational adjustment and binding-energy change are changes of the equilibrium
state, not extra entropy/heat inputs. R2 is an open-subsystem balance and has no general
nonnegative sign theorem. Do not call q_E “nonchemical heat.” Changing the energy zero by
constants c_i requires `rho'=rho-sum c_i n_i`, `mu_i'=mu_i-c_i`,
`q_E'=q_E-sum c_i Gamma_i`; R2 is unchanged.

Since `D n_B+n_B theta=Gamma_B`, s=n_B s_b gives

```text
T n_B D s_b=q_E-sum mu_i Gamma_i-T s_b Gamma_B.                 (R3)
```

The last term is denominator bookkeeping for entropy per baryon, not a new energy flow or a
separately established entropy flux carried by destroyed particles. With baryon creation or
loss, enclosed current B is not a comoving shell label; use a worldline label or explicit
source-aware shell transport instead of differentiating at fixed enclosed B.

Let `s_i=(partial s/partial n_i)_T`, `c_V=T(partial s/partial T)_n`. Substituting
`Ds=(c_V/T)DT+sum s_i(Gamma_i-n_i theta)` into R2 proves

```text
c_V DT=q_E-sum_i(mu_i+T s_i)Gamma_i
       -T(s-sum_i n_i s_i)theta.                              (R4)
```

Thus zero entropy-energy residual is not zero exact temperature drive. The chain rule applied to Gibbs gives
`mu_i+Ts_i=(partial rho/partial n_i)_T`. For a nonrelativistic degenerate free gas,
`s=(pi^2/2)n k_B^2 T/E_F,kin`, so `Ts_,n=(pi^2/6)(k_BT)^2/E_F,kin`.
At 10^8 K this is 1.22 eV for E_F,kin=100 MeV (roughly 0.9-2 eV for the representative
interior range, larger in low-density layers); the compression coefficient is `(2/3)Ts`
per NR species. Ultrarelativistic species have a different coefficient (e.g. Ts/3).
These eV/event corrections are small against 10-30 MeV/event but matter at a sign floor.
They are omitted by controlled frozen-v1, not identically zero in a physical finite-T star.
For P0 Gamma_n=-R_a and q_E=-mu_n R_a in common erg units, fixed-volume R4 gives **+Ts_,n R_a**.
For a uniform occupied-sea sink, Sommerfeld expansion gives
`<mu-E>_T=(2/5)E_F,kin-(pi^2/3)(k_BT)^2/E_F,kin`.

The finite-T floor is CHANNEL-WEIGHTING DEPENDENT (R2 E-3/E-26). Smooth P0 and uniform-sea
P1 have scales (pi^2/6) and (pi^2/3) times (k_BT)^2/E_F,kin, about 1.2 and 2.4 eV/event
at 10^8 K and E_F,kin=100 MeV. A source concentrated within k_BT of the Fermi surface can
have an O(k_BT) residual of either sign, approximately 8.6 keV/event at 10^8 K.
Every physical channel must redeclare its own finite-T weighting floor; eV/event is NOT universal.

The frozen controlled thermal function Udot=C_* T_infinity_dot is evaluated at fixed B;
Eeq_dot=C_(MeV->erg)mu_B^infinity Bdot uses the cold equilibrium tangent. They explicitly
omit finite-T B-dependent state derivatives, including (partial U_th/partial B)_T Bdot
and the finite-T correction to mu_B. These are in the same small finite-T class for the
smooth toy (the NR U_th derivative example is approximately 0.6 eV/event at those parameters).
A future exact physical conservation oracle must own these terms along with the R4 terms.

Extensions (not installed): heat flux contributes `-nabla_mu q^mu-q^mu a_mu` to the energy
projection; viscous stress contributes `-pi^{mu nu}nabla_mu u_nu`, positive for dissipative
Navier-Stokes constitutive laws. Diffusion j_i gives `+sum mu_i nabla_mu j_i^mu` in R2 before
redefining the entropy current. It is NOT generally boundary-only: integration by parts leaves
chemical-gradient work; the full entropy production contains `-sum j_i.grad(mu_i/T)`.
Only equilibrium/no-flux restrictions remove that bulk term. Reaction bulk viscosity represents
the same chemical reaction dissipation already owned by eta R; enhanced neutrino radiation is
still its separate loss, not identical to mechanical viscosity. No second beta-viscous channel.
These qualifications correct overbroad shorthand in R0 Q8/Q15, without changing the single-fluid toy.

## 4. G2, G19: local-to-infinity factors and standard thermal limit

A static local event has `d tau=e^Phi dt`, `E_infinity=e^Phi E_local`. Thus

```text
S_i=dot N_i^src=int e^Phi Gamma_i dV
T_infinity dot S_entropy=int e^(2Phi)[q_E-sum mu_i Gamma_i]dV.   (R5)
```

The second line follows also by integrating the entropy-current divergence over a comoving
worldtube with no unowned boundary flux and using T_infinity=e^Phi T. The quasistatic static-
slice approximation is declared; an arbitrary rotating/dynamic geometry is not justified by
copying these lapse factors. Uniform redshifted neutral conjugates permit the chemical part
`-C_(MeV->erg)g_y^{infinity T}S_y` for MeV conjugates. The escaping luminosity at infinity uses
only products leaving the STAR: `L_esc,star^infinity=C_(MeV->erg)int e^(2Phi)sum E_esc,star,a R_a dV`
when event energies are MeV.

| Quantity | Local-to-global measure | Units / reason |
|---|---|---|
| Number source | e^Phi dV | count/s; proper time to coordinate time |
| Entropy-energy/power | e^(2Phi) dV | erg/s; clock plus energy redshift |
| Fixed-metric test/reservoir energy | e^Phi dV | erg; energy redshift only, not the self-gravitating ADM mass functional |
| Cowling susceptibility | e^-Phi C_y dV | count/MeV; delta g_local=e^-Phi delta g_infinity |
| Geometry conversion | km^3 to cm^3:10^15; km^3 to fm^3:10^54 | chemical C_y in fm^-3/MeV needs fm^3 measure |

For Gamma_B=0, Gamma_i=n_B DZ_i and R3 becomes
`T Ds_b+sum mu_i DZ_i=q_E/n_B`. For radial luminosity flux,
`q_E=-e^-Lambda e^-2Phi/(4pi r^2) d(Le^(2Phi))/dr` (plus neutrino/local sources).
This is the full-chemical-potential form of Thorne (3.11-5); Thorne's alternative rest-mass
conversion q in (3.10) cannot be added to a rest-mass-inclusive version a second time.
Supporting Thorne pages are specified in section 11; R5 agrees with FR05 (4),(5),(15), R06(6).

For a positive beta neutron-decay rate, Gamma=(-1,+1,+1)DeltaGamma in the **(n,p,l) species
basis**; in governed y use the corresponding L column times DeltaGamma_l. Then
`-sum mu_i Gamma_i=eta_l DeltaGamma` in common energy units (C_(MeV->erg) multiplies the
MeV eta expression for erg power). Beta neutrinos have q_E=-Q_nu, not zero merely because
there is no external injection. Since eta_l(r)=e^-Phi eta_l^infinity,

```text
R_l=int e^Phi DeltaGamma_l dV
L_H^infinity=C_(MeV->erg) sum eta_l^infinity R_l
Delta P_beta=L_H-Delta Lnu.                                  (R6)
```

This exactly recovers FR05 (38), ADR-0014 sections 3.5-3.7 and the existing reaction/thermal
implementation, with no new normalization or sign convention.

## 5. G3, G16, G20: direct events, hole identity and cold bracket

For directed event a, `Gamma_i^a=nu_ia R_a`, R_a>=0. E_esc,fluid is the energy leaving the
ORDINARY THERMAL FLUID, which may go to infinity or to a retained nonthermal X reservoir.
With no external input,

```text
q_E,a=-C_(MeV->erg)E_esc,fluid,a R_a
q_dir,a=C_(MeV->erg)[-sum_i nu_ia mu_i-E_esc,fluid,a]R_a.                  (R7)
```

Here mu and event energies are MeV and q is erg cm^-3 proper-s^-1. External incident energy,
if independently owned, adds C_(MeV->erg)E_in R_a. Incoming charges and angular momentum
also require the external-inflow owner (section 12). Split
`E_esc,fluid=E_esc,star+E_X` at production; J_X is the redshifted transfer into X, not star
luminosity. Later X release returns energy once; reservoir accumulation and later escape have
their own ledger. Do not subtract both inclusive fluid escape and E_X again.

Neutron disappearance satisfies

```text
Q_dir=mu_n-E_esc,fluid=(mu_n-E_n)+E_dep,
E_n=E_dep+E_esc,fluid.                                       (R8)
```

**FERMI-HOLE HEATING IS NOT AN ADDITIONAL INDEPENDENT TERM.** The microscopic hole is
(mu_n-E_n)R_a in MeV per volume/time; the full chemical term +mu_n R_a alone is not the hole until occupied-
state energy removal is included. MPR21 p 2/Eq 9 and preceding discussion decompose collision
nN -> n' N into converted-neutron hole `(mu_n-E_n^i)`, spectator hole `(mu_N-E_N^i)`, and kicked
spectator `(E_N^f-mu_N)`. Energy conservation `E_n^i+E_N^i=E_n'+E_N^f` makes their sum
`mu_n-E_n'`. Equal neutron/mirror masses permit their common rest-subtracted convention.
A bound n' can leave the thermal fluid without leaving the star; binding is an independent
kinematic conclusion (when the product Lorentz factor gamma_X<e^-Phi) or G22 p 2's statement, not MPR prose. This is supporting
uncatalogued source evidence plus independent R8 algebra, not a newly accepted rate.

Exactly ONE representation per process:

| Representation | Components permitted together |
|---|---|
| R-a | Chemical -sum mu_i Gamma_i plus total energy transfer -E_esc,fluid,a R_a (+ external input), converted together once |
| R-b | [sum_removed(mu_i-E_i)+sum_created(E_j-mu_j)+other independently deposited energy] R_a, with stoichiometric multiplicities and one common rest-inclusive zero; neutron specialization is [(mu_n-E_n)+E_dep]R_a |
| R-c | Complete [-sum_i nu_ia mu_i-E_esc,fluid,a]R_a; neutron specialization (mu_n-E_esc,fluid)R_a |

R-a/R-b/R-c use one event rate and energy convention; convert their MeV powers once for the
thermal boundary. Adding a separate hole or E_dep atop R-a/R-c, dot E_chem as heat, dot M_eff as heat, generic
fluid work, duplicate reaction-bulk-viscosity dissipation or bound-product energy twice is
forbidden. The H1-H10 refusal table in section 18 makes these process-level tests.

### Illustrative channel, not the contract foundation

For n -> chi gamma, using one canonical local energy convention, E_n=E_chi+E_gamma:
photon absorbed/chi leaves fluid gives `Q=mu_n-E_chi=(mu_n-E_n)+E_gamma`; both leave fluid gives
`Q=mu_n-E_n`; photon absorbed and fraction f_chi deposited gives `Q=mu_n-(1-f_chi)E_chi`.
A bound inert chi is an X-reservoir transfer, not L_esc,star. Reverse conversion/Pauli blocking/
chemical equilibration may require Regime II. RMF canonical vs kinetic energies and actual
channel rates remain unratified; Sigma0 cancels in mu_n-E_n but can affect vertex kinematics.

### Cold bound and uniform-sea diagnostic

For a cold ordinary neutron sink with no external energy and no hidden extraction of another
reservoir, occupied E_n<=mu_n, E_dep>=0, hence
`0<=Q_dir`, `Q_dir>=mu_n-E_n`. If outgoing freely escaping products have total rest energy
m_out c^2, every truly escaping massive product must satisfy `e^Phi E_local>=m c^2`.
The sharper pointwise bound is `Q_dir<=mu_n-e^-Phi m_out c^2` in common local units.
At infinity the upper power is `C_(MeV->erg)mu_n,actual^infinity |dot B|
-L_esc,rest,star^infinity-J_X^infinity`; nonnegative J_X only tightens the bound.
At equilibrium, or for the specified neutron-sink branch with delta mu_n^infinity<=0,

```text
0<=P_dir^infinity<=|dot E_eq|-L_esc,rest,star^infinity,  E_eq=M_eq c^2. (R9)
```

Here `L_esc,rest,star^infinity=sum_a int e^Phi R_a (sum_truly_escaped m_j c^2)dV`: the
minimum Killing energy of each truly escaping massive product at infinity times the
coordinate event rate. This has **one lapse**, not two lapses on local rest mass. A two-lapse
local-rest-energy lower bound would be weaker and is not this definition. Use rest energies
in erg in this luminosity integral, or convert MeV rest energies exactly once.
Do not extend this equilibrium-mass upper bound to arbitrary initial disequilibrium, charged
stoichiometry, external injection, or bound-X energy definitions without rederiving it.
“Invisible” does not generally mean “no heating”: a process sampling occupied states below
the Fermi surface retains hole energy. No universal strictly positive lower bound is claimed
for a threshold-tuned Fermi-surface process.

For uniform momentum occupation of a nonrelativistic sea, <p^2>=3p_F^2/5, giving
`<mu-E>=(2/5)E_F,kin`. For relativistic free particles,

```text
E_F=sqrt(p_F^2+m^2)
<E>=3/[8p_F^3] {p_F E_F(2p_F^2+m^2)-m^4 asinh(p_F/m)}
<mu-E>=E_F-<E>.                                              (R10)
```

R1 scratch obtains 54.0379 MeV at p_F=530 MeV,m=939.56542052 MeV, ratio 0.38827 to E_F,kin.
R0 reports approximately 22/33/46 MeV at n_B=0.16/0.30/0.50 fm^-3 and54 MeV centrally.
These are local free-gas diagnostics, not universal constants or a measured whole-star mean.
Rate-weighting and radial integration are required for that mean; a momentum-dependent rate
cannot borrow the uniform-sea average without checking it.

## 6. G4, G10, G11: moving equilibrium chemistry and structural t

Define the physical departure at the actual B and spin:

```text
delta N_y=N_y-N_y^eq(B,Omega)=L ell,   ell=P delta N_y
b^T delta N_y=0 identically,  eta=-Z ell
S_y=int e^Phi Gamma_y^BNV dV,  dot B=b^T S_y
 t=(partial N_y^eq/partial B)_Omega,  b^T t=1
 dot N_y=L R+S_y
 dot N_y^eq=t dot B+L I_Omega 2Omega dotOmega
 dot ell=R+sigma-I_Omega 2Omega dotOmega
 sigma=P(S_y-t dot B)=S_l-t_l dot B.                         (R11)
```

Since `b^T(S_y-t dot B)=0`, the source difference is uniquely L sigma. This is not an
arbitrary projection of a baryon-changing source: the physical moving-reference subtraction
has already made it baryon-neutral. Therefore Z **plus t** suffice for controlled chemistry:

```text
eta_dot=-Z(R+sigma)+2W Omega dotOmega+Zdot Z^-1 eta,
W=Z I_Omega; Zdot=0 for controlled frozen-v1.                 (R12)
```

This follows by differentiating eta=-Z ell; the positive Zdot sign is governed by
`docs/adr/ADR-0014-secular-rotochemical-evolution-contract.md:380`. Zero BNV recovers that ADR
exactly; this proposal neither alters its implementation nor activates its future Zdot term.
If eta_initial=0 and spin forcing=0, **S_y=t dot B implies eta=0 identically**, even along a
sliding coefficient sequence. This is the physical negative control.

For neutron-only removal S=(dot B,0,0), dot B<0, `sigma_l=t_l|dot B|>0` on the free-gas fixture;
`eta_dot=-Z t_l |dot B|` has both components negative. H_l is odd and has the sign of eta_l,
so R_l<0 (p+l -> n+neutrino capture), eta_l R_l>=0. Removal makes the star neutron-poor
relative to its new equilibrium target; simply holding the old metric/target gives the wrong
magnitude. Lepton-only S_l without the target would give zero drive and is also wrong.

### Construction and error provenance of t

At Omega=0, `t_i=B_i/B_B`, where `B_i=dN_i^eq/d epsilon_c` and `B_B=dB/d epsilon_c` are
from the existing Phase-5B structural owner. Ratios are unchanged if the common derivative
coordinate is ln epsilon_c. This is a new derived Phase-6 view of governed inputs, not an
already ratified BNV coefficient object. Authority: ADR-0011 sections 3-4,
`CompactStar/Analysis/src/ParticleNumberResponse.cpp:365` (independently solved neighbors,
five/three-point stencils and ladder) and `:450` (baryon denominator/error).
Retain star/provider bytes, domain/surface policy, stencil/step ladder, denominator condition,
error semantics, charge/baryon sum rules and input currency; never extract a lifetime-free ratio.

R1 independently re-extracted the governed baselines as scratch mathematical oracles
**ON THE STRUCTURE-1 FREE-GAS FIXTURE**. Long digit strings for t, k and drives are
arithmetic-reproducibility oracles on governed bytes, NOT physical precision. All displayed
drive coefficients mean eta_dot/|Bdot| in MeV/count; multiplying by |Bdot| in count/s gives MeV/s:

| Quantity | R1 arithmetic oracle / propagated budget |
|---|---|
| t in n/e/mu order | (0.9657700849496014, 0.030852171225661786, 0.0033777438247248118) |
| B_B; numerical error | 1.6831408136820063e59; 1.2652133015974696e52 in the baseline derivative units |
| Conservative ratio numerical budgets | (1.4428870e-7, 4.3530477e-9, 2.1785265e-9) absolute |
| b^Tt-1; t_p-t_e-t_mu | approximately -1.21e-14; +1.20e-14 |
| eta_dot_e / abs(dot B) | -1.4302859054377823e-55 MeV/count |
| eta_dot_mu / abs(dot B) | -3.628096720520591e-55 MeV/count |

The ratio budget is `(error_Bi+|t_i|error_BB)/(|B_B|-error_BB)`, not a certified interval.
R0's (1e-7,3e-9,2e-9) was rounded; the slightly larger explicit propagated values above
prevent invented precision: relative t budgets are about (1.5e-7,1.4e-7,6.5e-7), supporting
roughly 6–7 significant digits rather than the printed 16–17. Propagate both t and governed Z
errors for a physical drive budget. Quote fixture drives as (-1.4303e-55,-3.6281e-55) MeV/count;
these rounded values do not remove the separate Cowling physical-model uncertainty.
Phase-5B hash `7588f0e9cd62f5b6be48bb725e4d0ba6e47b64f50d843117e70c17a78beeb5fa`;
Phase-5C hash `7027aa6179fe9111d76586d467d007c4223ec9d26ece0fee773327973c77bfe7`.

### Future charts and domains

Use physical Omega in s^-1 and q_phys=Omega_phys^2. Define
`K_phys=K_repo/c_km_s^2=L I_Omega` (count s^2), distinct from the governed
geometric K_repo (count km^2). Then `N_eq(B,q_phys)=N_0(B)+q_phys K_phys(B)+O(q_phys^2)`
and `t_i(B,q_phys)=N_0,i'(B)+q_phys K_phys,i'(B)+O(q_phys^2)`.
With J=I(B,Omega)Omega, chain rule gives

```text
t^J=t^Omega-2Omega^2 K_phys (partial_B ln I)_Omega
                    /[1+Omega(partial_Omega ln I)_B]+O(Omega^4). (R13)
```

Fixed-J energy derivatives use Mc^2; fixed-Omega use Mc^2-Omega J. Spin-off makes these
charts coincide. No torque, rotating sequence or O(Omega^2) BNV extension is implemented.

For a fixed-isobar subdomain D, `b^T t^D=dB_D^eq/dB`, not necessarily1. If its required
baryon boundary flux is carried at a declared composition Y_c with b^TY_c=1, define
v=S^D-t^D dot B and boundary flux=-b^T v; then the baryon-neutral forcing is
`v-Y_c b^T v`, giving `sigma^D=P(I-Y_c b^T)(S^D-t^D dot B)`.
This conditional boundary model must match ADR-0011 PN8 and account for actual transport;
it is not a universal consequence of the whole-star algebra. First toy uses whole-star D.

## 7. G5, G13: retire raw-G chemistry; t versus k diagnostic

Within its declared Phase-5 scope, G_y remains the R2006/Cowling source used to build `Z=L^T G_y^-1 L`. Its fixed-metric
baryon direction `k=G_y b/(b^T G_y b)` is diagnostic only:
`k=(0.9922178920,0.00748453949,0.000297568506)`,
`t_e/k_e approximately 4.12`, `t_mu/k_mu approximately 11.35` **ON THE STRUCTURE-1 FREE-GAS
FIXTURE**; these ratios and drive magnitudes are not EOS-independent numbers.
Because b^Tk=1 and L^TG^-1k=0, decompose S=k dot B+L P(S-k dot B), proving

```text
D_old S=-L^T G_y^-1 S=-Z(S_l-k_l dot B).                     (R14, REJECTED physical route)
```

For baryon-changing sources this agrees with R12 only if t=k (invertible Z); for dot B=0
both routes already agree. Their difference is `-Z(t_l-k_l)dot B`. R1 reproduces
`(+1.0860069714e-55,+3.2838177865e-55)|dot B|` on the physical sliding control, where the
correct answer is zero. The old `-3.442789339881122e-56 MeV/count` in both channels is a
**COWLING k-ROUTE NEGATIVE ORACLE**, not physical neutron-sink response. Do not install D_old,
use k as the physical null, or reconstruct actual individual changing-B potentials with G_y^-1.
The forward-looking BNV-seam language in ADR-0013 Q1
(`docs/adr/ADR-0013-corrected-rotochemical-chemical-coefficients.md:212`) and ADR-0014 section 3.17
(`docs/adr/ADR-0014-secular-rotochemical-evolution-contract.md:472`–`:478`) is explicitly narrowed
by PROPOSED ADR-0015 section 7, subject to owner ratification. "Cannot be represented in the
fixed-baryon two-channel space" applies to RAW S_y. The moving-reference source
Sigma_y=S_y-t Bdot obeys b^T Sigma_y=0 and Sigma_y=L sigma identically for every charge-consistent
source. The physical seam is {Bdot,sigma}, sigma=P(S_y-t Bdot), using qualified t and governed Z;
the raw-G_y map is rejected. G_y retains its unreduced coefficient/Cowling ownership and k
diagnostic role. Phase-5 coefficient mathematics, Q/Z/W ownership, governed baselines and all
Phase-5 results remain intact. The future owner request must explicitly acknowledge this
narrowing; this document neither edits nor ratifies the accepted ADRs.

For invertible symmetric G, nonzero a_G=b^TGb and invertible Z, set
H_proj=G-Gbb^TG/a_G. It annihilates b, so H_proj=L C L^T with C=P H_proj P^T.
Then `ZC=L^TG^-1 H_proj P^T=I`, because L^Tb=0 and L^TP^T=I. Thus

```text
Z^-1=P[G-Gbb^TG/(b^TGb)]P^T.                                (R15)
```

Positive definiteness of the full G is not required for this identity. If a true global
response obeys `G_true b=(dB/dmu_B^infinity)t`, then
`Z_true^-1=G_true^(e,mu)-(dB/dmu_B^infinity)t_l t_l^T`.
The dedicated diagnostic records the R0 numerical evidence, its missing underlying scratch,
and R1's independent algebra/protected-input checks. G_true is NOT a new governed coefficient.

## 8. G6, G8: individual potentials and chemical energy reservoir

At fixed B, delta N=L ell and `E_chem=(1/2)ell^T Z ell=(1/2)eta^T Z^-1 eta>=0`.
The linear term vanishes because g_eq=mu_B b and b^TL=0. It is a state function/reservoir,
never a thermal-RHS addend. V=eta^T Z^-1 eta=2E_chem is the Phase-5 Lyapunov quantity.
E_eq(B)=M_eq(B)c^2 already contains baryon-direction curvature.
No additional moving-baryon quadratic store from the old fixed-reference draft survives.

For the spin-off first-order potential response, write
`E(N)=E_eq(B)+(1/2)ell^T Z(B)ell`, `ell=P[N-N_eq(B)]`.
To first order in ell, varying N gives
`dE=mu_B b^TdN-eta^TP(I-tb^T)dN+O(ell^2 dB)`; hence

```text
delta g^infinity=-(I-bt^T)P^T eta+O(eta^2),
delta mu_n^infinity=t_l^T eta+O(eta^2),
-delta g^T S=eta^T sigma+O(eta^2 dot B).                     (R16)
```

Equivalently, the equilibrium Hessian obeys H_true t=mu_B' b; symmetry and b^TL=0 imply
t^T delta g=0, and L^T delta g=-eta fixes R16 uniquely. R16 is exact in the declared frozen
quadratic model; it is not an exact all-orders formula for a changing star. If Z(B) variation
is retained through quadratic order, the gradient additionally contains
`+(1/2)b ell^T Z_,B ell`. Thermal B-dependence and rotation likewise require their derivatives.
No raw Cowling inverse supplies those individual physical potentials.

### Corrected variable-Z derivative: explicit R0 erratum

Direct chain rule and R12 give two equivalent forms:

```text
Echem_dot=eta^T Z^-1 eta_dot
         -(1/2)eta^T Z^-1 Zdot Z^-1 eta                      (R17a)
        =-eta^T(R+sigma)+2Omega dotOmega eta^T I_Omega
         +(1/2)eta^T Z^-1 Zdot Z^-1 eta.                     (R17b)
```

The **minus half belongs to the chain-rule form** R17a. R0 Q29/G8 and the owner's summarized
candidate formula kept that minus sign after substituting an eta_dot that already includes
+Zdot Z^-1eta. That is algebraically inconsistent; R17b has **plus half**. Independently,
`Echem_dot=-eta^T ell_dot+(1/2)ell^T Zdot ell` proves the same result. Scratch scalar example
ell=3,Z=2+t gives dE/dt=+4.5; the erroneous substituted form gives -4.5. A separate matrix
finite difference also agrees with R17a/R17b. This is a resolved correction, not a choice
between conventions or a Phase-5 contradiction. Zdot=0 makes the frozen results identical.

Thus beta contributes -eta^TR<=0, BNV contributes -eta^Tsigma, and spin contributes
+2Omega dotOmega eta^T I_Omega. The full filling derivative need not be positive for arbitrary
initial states/history. “While storage fills” is a condition, not a universal inequality.

## 9. G7: actual versus equilibrium direct power and closed controlled ledger

Use `L_out,fluid^infinity=L_esc,star^infinity+J_X^infinity`; in the prompt-escape first toy
J_X=0. For general charge-consistent S,

```text
P_dir^infinity(actual)=-C_(MeV->erg)g_actual^{infinity T}S-L_out,fluid^infinity
P_dir^infinity(eq)=-C_(MeV->erg)mu_B^infinity dot B-L_out,fluid^infinity
P_dir^infinity(actual)=P_dir^infinity(eq)+C_(MeV->erg)eta^T sigma                          (R18)
```

to R16's order, exactly in the frozen quadratic model. For a neutron sink,
`P_dir^infinity(actual)=C_(MeV->erg)mu_n,actual^infinity |dot B|-L_out,fluid^infinity`;
`mu_n,actual^infinity=mu_B^infinity+eta^Tt_l` and sigma=t_l|dot B| prove the **plus** sign
in R18. Eta<0 means eta^Tsigma<0: the actual hole is shallower. No extra heat is created by
renaming the same energy decomposition.

The controlled thermal equation uses actual potentials:

```text
Udot=C_* T_infinity_dot
    =P_dir^infinity(actual)+C_(MeV->erg)eta^T R-Lnu_eq-Delta Lnu-Lgamma-Lother
    =P_dir^infinity(eq)+C_(MeV->erg)eta^T(R+sigma)-Lnu_full-Lgamma-Lother.          (R19)
```

For spin OFF, frozen Z,t, a consistent equilibrium-energy tangent and the declared thermal
model, R17 cancels its chemical transfer and gives

```text
Eeq_dot+C_(MeV->erg)Echem_dot+Udot=-L_out,fluid^infinity-Lnu_full-Lgamma-Lother,
Eeq_dot=C_(MeV->erg)mu_B^infinity dot B.                                 (R20)
```

Here U and E_eq are erg, E_chem is MeV, and all powers on the thermal RHS are erg/s at infinity.
The strongest closure proof is R24 + R5 with actual potentials around the hydrostatic/diffusive-
equilibrium star, independent of the quadratic split; section 3 explicitly declares the finite-T
state derivatives omitted by this cold/frozen implementation proposal.

For retained X, add EX_dot=J_X^infinity-P_return^infinity-L_X,escape^infinity and include P_return once in R19; R20 then
closes with star-escaping terminal luminosities. A bound product cannot be called an infinity
luminosity at birth. Finite-interval identity for prompt escape is
`Delta(E_eq+C_(MeV->erg)E_chem+U)=-int(L_esc,star^infinity+Lnu_full+Lgamma+Lother)dt`, with initial/end states and
common energy zero explicitly specified; subtract the matched no-BNV identity for differences.
There is no freely selectable equilibrium-reference work input. Prescribed spin, variable Z,
variable heat capacity and structural evolution require the actual state derivatives/rotational
work, not an extrapolation of this spin-off identity. No physical conservation oracle can
silently omit those terms outside the controlled scope.

Relative to P_dir(eq), at spin OFF/frozen Z the correction is
`C_(MeV->erg)eta^T(R+sigma)-Delta Lnu=-C_(MeV->erg)Echem_dot-Delta Lnu`.
It is nonpositive **while Echem_dot>=0**; in QSS it is exactly -Delta Lnu. Relative to
P_dir(actual), the beta increment remains C_(MeV->erg)eta^TR-Delta Lnu, which can have either sign.
For P0, E_esc,fluid=actual mu_n is itself state-dependent and P_dir(actual)=0; it must not
be interpreted as an absolute fixed escape-energy model.

## 10. G9, G14: beta polynomials, bounds and reachable QSS

The roots, minima, B1, B2 and ratio bounds below are properties of the governed NON-SUPERFLUID
R1995/FR2005 polynomial model, not universal under superfluidity.
Let u=xi/pi, xi=eta/(k_B T_infinity), F-1=f, and h=xi H. The governed functions
(`CompactStar/Physics/Rotochemical/UrcaImbalanceFunctions.hpp:13`) give

```text
f_D=(1071u^2+315u^4+21u^6)/457
h_D=(714u^2+420u^4+42u^6)/457
f_M=(22020u^2+5670u^4+420u^6+9u^8)/11513
h_M=(14680u^2+7560u^4+840u^6+24u^8)/11513
p_D=h_D-f_D=(-357u^2+105u^4+21u^6)/457
p_M=h_M-f_M=(-7340u^2+1890u^4+420u^6+15u^8)/11513.             (R21)
```

Each process contributes Ltilde T^q p to Delta P_beta. The small-xi coefficients are
-357/(457pi^2), -7340/(11513pi^2). Positive incremental roots and minima, independently
recomputed in R1 at 60 decimal digits, are:

| Process | Positive nonzero root | Minimum p | xi at minimum |
|---|---|---|---|
| Modified | 4.9097100289241306 | -0.4676589500024851 | 3.6125406881312751 |
| Direct | 4.7870134733368997 | -0.5277746053621574 | 3.4972939247699681 |

Eta=0 gives zero, not strict cooling. For0<|xi|<root the increment cools, independent of sign
eta. Different channels can occupy different sign domains; report the total weighted sum.
Minimizing these explicit polynomials proves instantaneous **B1** bounds:
`-Delta P_beta<=0.467659 Lnu_eq,M+0.527775 Lnu_eq,D`, summing enabled process normalizations.
The ratio p/h is a positive-weight average of the monomial ratios: M=(-1/2,1/4,1/2,5/8),
D=(-1/2,1/4,1/2). Thus `-1/2<=Delta P_beta/L_H<=5/8` for mixtures where L_H>0;
pure direct upper limit1/2. At eta=0 use the limit rather than dividing by zero.

**B2 is a source-driven scale, not an unconditional bound for arbitrary stored initial eta.**
For a neutron-sink capture branch with both imbalances in the cooling interval and nondecreasing
chemical storage, R17 implies `L_H=C_(MeV->erg)eta^TR<=-C_(MeV->erg)eta^Tsigma` (spin OFF/frozen Z). Combining
`-Delta P_beta<=L_H/2` with sigma=t_l|dot B| gives
`-Delta P_beta/[C_(MeV->erg)|dot B|]<=(1/2)xi_root k_BT(t_e+t_mu)` in MeV/baryon.
At10^8 K the M value is **0.72411 keV/baryon**, D0.70601 keV. The assumptions matter:
a releasing pre-existing chemical reservoir, arbitrary forcing or vanishing source rate invalidates
an unconditional per-source instantaneous bound. Compared with 10-30 MeV/event, this cooling
scale is O(10^-4) or smaller. It does not prove every history heats at every time.

Rates have the existing semantic normalization
`R_l=C_(MeV->erg)^-1 sum_(a in l)(Ltilde_a/k_B)H_a(xi)T^(q_a-1)`.
For spin-off QSS with invertible Z, `R=-sigma` independently of Z; every driven channel must
have nonzero reaction normalization. In the large-|xi| modified
regime H_M approximately C_H xi^7, C_H=24/(11513pi^8); let
A_l=C_(MeV->erg)^-1 sum Ltilde_a over the applicable M channels. Then

```text
|eta_qs,l|=(sigma_l k_B^8/(A_l C_H))^(1/7), eta_qs,l<0
xi_qs,lin=-sigma_l k_B/[A_l H'_M(0) T^7]
J_reaction=Z diag(A_l H'_M(0) T^6/k_B^2).                    (R22)
```

Relaxation times are inverse eigenvalues of J_reaction; a single-channel reduction has
`tau_lin=k_B^2/[Z A H'(0)T^6]`. R0's extra symbol C in this expression must mean the energy
conversion if Ltilde is in erg; it is not heat capacity. The explicit A definition removes
that ambiguity. A drive-dominated estimate `t_build~|eta_qs|/|Zsigma|` scales as Gamma^(-6/7),
while |eta_qs| scales Gamma^(1/7). It is an estimate, not a replacement for reaching QSS.

R0's corrected **LARGE-|xi| MODIFIED-URCA QSS** build estimates (attributed evidence,
not R1 trajectories; R2 E-17 reproduced their asymptotic scaling):

| Abstract Gamma, yr^-1 | Large-xi build time, yr | Applicable interpretation |
|---|---|---|
| 1e-17 | approximately 3.4e12 | Large-xi asymptotic QSS not built in 1e10 yr; says nothing against hot linear QSS |
| 1e-14 | approximately 9e9 | Large-xi build estimate marginal on stellar age; check the actual regime |
| 1e-10 | approximately 3e6 | Large-xi scaling example, not selected physical rate |

At 10^8 K the corresponding asymptotic xi estimates are approximately 0.43, 1.15 and 4.27;
none is a deep large-|xi| example. R2 E-17 (Reviewer E rederivation 10) instead evaluates
fixture linear relaxation at approximately 3.5e4 yr (electron) and 2.9e4 yr (muon), with
tau_lin proportional to T^-6. These scales inherit the governed controlled normalizations.
Hot linear QSS may be rapidly reached; cooling makes tau_relax grow strongly, allowing
freeze-out. Later reaction-free/drive-dominated evolution can become applicable; large-|xi|
QSS is a distinct later/asymptotic possibility requiring its own reachability check.
EVERY QSS statement requires tau_relax(T,xi) << relevant evolution time IN THE APPLICABLE REGIME.

When reactions remain negligible (early t<<tau_lin(T), or after freeze-out with a verified
small reaction term), `eta(t)=eta(0)-Z int sigma dt` is the appropriate transient
oracle. Conditional photon-balanced asymptotics: reached large-xi M QSS with P0 gives net
power proportional to |dot B|^(8/7), hence T_s proportional to |dot B|^(2/7); fixed positive
MeV direct energy per event gives T_s proportional to |dot B|^(1/4). Neither is a claim about
an unreached QSS or a transient cooling trajectory.

## 11. G12: stellar first law, local authority and derivation scope

Thorne lecture printed 172/PDF 8 defines full chemical potentials at fixed entropy and volume;
printed 183/PDF 19 Eq 3.3 uses static Killing energy; printed 189-190/PDF 25-26 derives luminosity
redshift and proper time separately. Printed207/PDF 43 Eq 3.62 is injection with the **same local
composition and specific entropy**, `delta E=e^Phi(rho+p)delta B/n_B`. Euler decomposes this
particular injection into T_infinity delta S+sum g_infinity delta N. Thorne explicitly says
injection energy generally varies with radius; Eq 3.62 alone is not a proof of arbitrary
independent multispecies/entropy variations. These pages were visually checked in R1.

A stronger static result is independently derived here. Assume spherical hydrostatic equilibrium,
Tolman temperature and uniform neutral redshifted conjugates, local neutrality, consistent
handling/cancellation of electrostatic contributions, and a vacuum-matched p_R=0 surface.
The background may be beta-disequilibrated; beta equilibrium is needed only when specializing
to g_eq=mu_B b. Thus actual diffusive-equilibrium potentials are valid in R24 and the closure.
In geometric units, let f=1-2m/r, M=int4pi r^2 rho dr, and integrate S and N_y with f^-1/2.
Varying the density integrals, including the metric-volume variation induced by delta m, and
using local Gibbs gives

```text
T_infinity delta S+g_y^{infinity T}delta N_y
 =int_0^R4pi r^2 delta rho(r) K(r)dr
K(r)=e^Phi(r)/sqrt(f(r))
     +int_r^R4pi x e^Phi(x)[rho(x)+p(x)]f(x)^(-3/2)dx.        (R23)
```

TOV implies `d(e^Phi/sqrt f)/dr=4pi r e^Phi(rho+p)f^-3/2`, so K'=0 and vacuum matching
K(R)=1. Restoring a common energy convention proves

```text
delta(Mc^2)=T_infinity delta S+g_y^{infinity T}delta N_y.      (R24)
```

R23–R24 use g in the same energy convention as Mc^2 and T delta S; convert interface MeV
conjugates once by C_(MeV->erg) if Mc^2 is expressed in erg.

A moving finite-pressure surface leaves `-4pi R^2p_R delta R` in delta(Mc^2) relative to the
right side. The finite p_cut fixture must retain/qualify that boundary when making a complete
stellar claim. A future finite-p_cut dM/dB oracle must specify the convention for
-4pi R^2 p_R dR/dB rather than silently using a vacuum-surface formula. Cold catalyzed g_eq=mu_B b gives `d(M_eq c^2)/dB=mu_B^infinity` at p_R=0.
R0 reports independent fixture agreement to approximately 3e-10; that nonlinear TOV computation
was not rerun in R1. Hartle1967 Eq 83-85 and109-110 support the cold first integral/count measure;
R06 Eq 2 supplies the neutral/electrochemical convention.

Uniform-rotation extension: for a stationary family with K=partial_t+Omega partial_phi,
thermal/diffusive conjugates are the appropriate co-rotating Killing multipliers (mu/u^t and
T/u^t), not an unqualified static e^Phi. If the constrained stationary variational problem is
extremal at fixed species counts, entropy and J, stationarity of
`Mc^2-Omega J-T_infinity S-g_infinity^T N` yields, along the family,
`delta(Mc^2)=Omega delta J+T_infinity delta S+g_infinity^T delta N`.
This is a **conditional first-order variational extension**, not a quotation of an absent
Bardeen theorem or a proof that the static kernel already covers arbitrary rotation.
Hartle1970 Eq 10/15/16 and Thorne printed 209/PDF 45 support constrained rotation/adiabatic
variations, not the complete unconstrained multispecies theorem by themselves.
The source-limited general stationary proof/quotable theorem is future scope; the explicit
static neutral derivation closes the first spin-off contract. This distinction corrects R0's
stronger source-attribution wording without blocking the controlled abstract ledger.

Reversible structural energy stays in E_eq. BNV spin work and product angular momentum require
separate Omega dot J ownership when enabled. A generic "work on fluid" thermal row is absent.
GR10 discusses small crust-strain heating for its spin-down assumptions; no universal BNV
strain/dissipation bound follows without a strain history. No crust or viscous mechanism is
installed. Numerical ratios for arbitrary damping are not promoted from R0 estimates.

## 12. G17: generic input classes and product-fate regimes

The generic ledger does not depend on n -> chi gamma. Its information dependencies are:

| Input class | Determines | What remains process-dependent |
|---|---|---|
| A sequence and total dot B | E_eq(B), mu_B, t, radius, I, surface gravity, conditional mass/rotation observables | EOS, sequence chart and validity |
| B ordinary stoichiometry | S_y, sigma, eta drive and standard beta response | Full channel/cascade charge and species sources |
| C event energy partition | E_esc,fluid, E_dep, direct thermal value | Energy zero, occupied-state/radial weighting, thresholds |
| D product fate | Escape, local SM thermalization, retained X or interacting X ownership | Transport/decay/capture/retention histories |
| E hidden-sector interactions/state | Added thermal transport, radiation, chemical constraints, ordinary weak-rate/kinematics feedback, eventually EOS | Couplings, N_X(t), T_X, composition and structure |
| F external inflow | Incoming energy, conserved charges and angular momentum for capture-induced or externally driven BNV | Incident flux/charge/spin ledger, linked to event identity without double-counted deposition |

The generic source callback may depend on time and state (T, eta, structure, accumulated sector):
Sigma_y(t,state)=S_y(t,state)-t_sequence Bdot(t,state). This Sigma_y is the baryon-neutral
moving-reference source; keep raw S_y and sigma=P Sigma_y separately typed/owned. This retains
ADR-0014 section 3.17's state-dependence while explicitly narrowing its raw-source seam meaning.

Each future channel MUST declare product fate: **PROMPT_ESCAPE**, **SM_THERMALIZATION**,
**BOUND_INERT**, or **BOUND_INTERACTING**, with mixed outcomes represented by separate weighted
branches. PROMPT_ESCAPE requires massive terminal products to satisfy e^Phi E_local>=m c^2
at production, plus an actual escape/transport model; energy below threshold belongs to retained X.
The flag covers terminal products, not just the first intermediate. SM_THERMALIZATION includes
absorption elsewhere inside ordinary matter: conserve Killing energy in transit and count its
deposition once. A tiny N_X/B is insufficient for ordinary-star validity.

Regime I requires every relevant omitted effect below a **declared measured error budget**:

| Requirement | Proposed measurable conditions |
|---|---|
| Accumulation/Pauli | Per-species source/depletion ledger; no appreciable occupation blocking in the produced phase-space support; E_F,X small relative to available production energy where that criterion applies |
| Stress/EOS/chemistry | Energy-density/pressure backreaction small; X-induced ordinary chemical shifts below a declared absolute imbalance/thermal resolution; no rapid reverse-equilibration condition mu_n=mu_X |
| Heat/transport | C_X/C_* small; omitted X radiative/heat-transfer power below both cooling and relevant direct/residual error budgets; negligible X opacity/conduction |
| Ordinary weak feedback from X | Independently bound X-mediated effective-mass shifts, new spectator/catalyst channels and modified Urca thresholds/normalizations; may remain even as DeltaB/B tends to zero |
| Frozen coefficients | abs(Delta B)/B small AND changes of Z,W,Ltilde,C_*,t,radius,g_s,I within declared budgets; moving support and threshold changes checked |

For a smooth nonzero coefficient X, estimate `|Delta B/B| |d ln X/d ln B|`; for matrices,
zeros, near-zero eta or P0 use scaled norms and absolute tolerances, not singular fractional
conditions such as |delta mu|<<|eta| at eta=0 or L_X<<P_dir=0. Time-dependent source support and
accumulation require direct monitoring. X-mediated ordinary weak feedback is a separate check:
it is NOT measured by an estimator proportional to DeltaB/B.
ON THE STRUCTURE-1 FREE-GAS FIXTURE, R2 E-9 finds d ln B/d ln epsilon_c approximately 0.18
and d ln N_mu^eq/d ln B approximately 54 (near-threshold muons). The large muon drift is
fixture-specific and must enter the depletion budget; measure N_mu^eq, t_mu and Z_npmu
sensitivities separately, not assign all three the same derivative. No universal numerical
tolerance is selected here.
Failure of any frozen condition requires exit from this **controlled frozen Regime-I contract**.
Ordinary coefficient drift alone can lead to a future evolving-background ordinary Regime I;
it does not by itself prove a physical hidden-sector Regime II. Failure of X thermal,
mechanical or chemical conditions triggers the corresponding Regime II subclass.

Thermal Regime II: an accumulated sector affects C, transfer or radiation while mechanically
small. Mechanical/EOS Regime II: stress/gravity/pressure significant. Chemical-equilibrium
Regime II: rapid reversible conversion changes the equilibrium constraints. X-mediated changes
to ordinary weak kinematics/normalization require the corresponding chemical/thermal feedback
owner even without appreciable baryon depletion. State-space seam:
(B,B_X,J,S,S_X), adding multiple X counts/temperatures as necessary. Proposed progression:
ordinary-star model -> **MixedStar-lite** (test-fluid structure, N_X ledger, separate or slaved
T_X) -> full two-fluid MixedStar. Existing MixedStar is not promoted or activated.
Slaving T_X requires its relaxation time short enough that P_exchange approximately L_X,rad;
otherwise integrate X thermal storage. Nothing in this progression is implemented.

A LINEARIZED TOY accumulation model illustrates why abundance alone fails: N_X approximately f_X R_B t,
scattering transfer `L_X approximately c f f' n_e sigma_eX k_BT R_B(t/2)` while
`P_dir=Q R_B`. Thus

```text
L_X/P_dir approximately c f f' n_e sigma_eX k_BT(t/2)/Q.      (R25)
```

The rate cancels only with linear accumulation/transfer, fixed f' and fixed Q. Goldman's actual
accumulated degenerate mirror-fluid model is NOT rate-independent: f' depends on N_X through
E_F'(N_X), so the cancellation is spoiled (R2 E-21; Reviewer E, rederivation 11). In its stated
degenerate scaling, fixed-temperature E_F' proportional to N_X^(2/3) makes L_X proportional to
N_X^(1/3), rather than N_X. G22 supports the qualitative conclusion that tiny N_X/N need not
mean negligible thermal influence. Coupling, accumulation time and STATE-DEPENDENT transport
control that influence; quantitative normalizations remain source-noted in section 17.

With added hidden loss L_X=lambda T^p and standard cooling still present, persistent deposition
sets the floor by `P_dir=L_std(T_floor)+L_X(T_floor)`. For nonnegative standard loss and p>0,
`T_floor<=(P_dir/lambda)^(1/p)`; equality requires lambda T^p to be the total loss.
Relative to a Regime-I BNV prediction the new loss can cool;
relative to the passive matched control it initially cools only where its net new loss exceeds
deposition at the control state. For a finite persistent loss law vanishing at T=0, the floor
prevents an unconditional claim of cooling to zero. Evolving rates/loss laws need a new check.

## 13. G18: matched control, thermal sign, lag and integrated energy

The comparator is the **PASSIVE SAME-INITIAL-B0 NO-BNV STAR**. Same initial star, T, eta,
ordinary beta/photon/neutrino models, spin history, solver and
tolerances; disable only BNV sources/product sector in the control. For the first toy both
runs use the corresponding frozen B0 coefficients, with the target departure defined by R11.
A future sliding-background comparison must explicitly evolve the B-dependent inputs and state
energy. In interpreting DeltaT_s, separately report changes due to radius, g_s, the envelope
relation, C_* and other structure from genuine thermal changes relative to that same-initial-B0
control; do not silently compare unrelated equilibrium stars.

At equal coordinate time define Delta T, Delta T_s, Delta Lgamma and

```text
Delta P_total=P_dir(actual)+Delta[L_H-Delta Lnu]
              -Delta Lnu_eq-Delta Lgamma-Delta Lother
              +Delta P_additional,owned
 d[U(T_BNV)-U(T_control)]/dt=Delta P_total.                   (R26)
```

Outer Delta is the difference of the two trajectories; inner Delta Lnu is each trajectory's
disequilibrium enhancement at its own T. Neither equilibrium neutrino feedback nor photons
cancel merely because the model is shared. U(T)=int C_*(T)dT is common in the frozen toy;
Delta U is not generally C(T)Delta T. Report int P_dir dt, int Delta P_total dt, Delta U,
Delta T, Delta T_s and Delta Lgamma together.

If cooling locally behaves as L approximately T^a_cool, linearized thermal lag is
`tau_th=C_*T/(a_cool L)`. The quoted approximately 1.4e5-1e6 yr range is instead the
energy-content/luminosity timescale `tau_energy=C_*T/L` on the controlled trajectory (R2 E-25);
a_cool is about 2.4–2.8. Both diagnostics inherit free-gas C_* and placeholder controlled
Ltilde and are not physical predictions. Do not attach tau_energy's numbers to tau_th.
Temperature remembers past power; early integrated cooling can later be erased by heating.
A power sign is not automatically an instantaneous surface-temperature difference.
Source-driven sliding introduces O(Delta B/B) coefficient effects, which are within the frozen
error band only after sensitivity checks.

Every NET HEATING / NET COOLING / NEAR ZERO / SIGN UNRESOLVED classification must specify
its observable AND time/window: DeltaT_infinity(t), DeltaT_s(t), DeltaLgamma(t), DeltaU_th(t),
or integrated DeltaP on a stated [t0,t1]. Use SIGN UNRESOLVED when the error band crosses zero. R0's measured temperature resolution1.74e-5 motivates a provisional1e-4 relative
reporting floor for that fixture, NOT a universal tolerance. Re-measure the future solver floor,
combine it with channel-weighting-specific finite-T, energy-partition, t/Z/Cowling and coefficient-drift errors,
and use absolute budgets near cancellation. A generic positive efficiency is forbidden.

Precise scientific sign statements: restricted cold Regime-I direct disappearance power is
nonnegative; beta-mediated Delta P_beta can cool; ordinary10-30 MeV direct terms dominate the
roughly 0.7 keV source-driven cooling scale. Therefore **generic Regime-I BNV cooling is not
supported**. Cooling needs a tuned small-direct-residual corner and applicable history/bounds,
not merely invisible products. Arbitrary pre-existing imbalance, imposed spin and mixed sources
require their own comparison. Regime-II cooling must name its comparator and energy-loss owner.

Future sign-map axes: event partition, source species direction, eta_e/kT, eta_mu/kT, T,
source radius/support, weak regime, finite-T threshold/thermal-tail weighting, N_X(t), product
fate and coupling/age. A thermal-tail channel can cool without violating the cold theorem. Separate direct,
beta-mediated and total signs. No numerical sign scan was run.

## 14. Proposed first toy and G1-G20 oracle coverage

Propose only: governed Phase-5D free-gas fixture, spin OFF, whole-star diffusive domain, uniform
abstract neutron sink per local proper time (`Gamma_n=-gamma n_n`, S_n=int e^Phi Gamma_n dV).
Use t from the qualified Phase-5B owner and Z from Phase-5C. Coordinate total-loss fraction
|dot B|/B is not automatically the local gamma; normalize and record the lapse/particle fraction.

| Partition | Definition | Oracle / limit |
|---|---|---|
| P0 | E_esc,fluid=**actual current local mu_n** | Zero cold direct entropy-energy residual; artificial isentropic Fermi-surface removal, not n->chi gamma; finite-T R4 remains |
| P1 | E_esc,fluid=E_n, uniform occupied-sea weighting | Hole-only heat; R10 analytic average; requires radial integration |
| P2 | E_esc,fluid=0 | Maximal ordinary retention mu_n/event within chosen no-other-reservoir toy |
| P1gamma later | Actual two-body n->chi gamma partition | Needs source-qualified kinematics/transport; not ratified here |

R0's two-tier **mathematical test-design examples**, not chosen physical rates or constraints:
P1/P2 around 1e-17 yr^-1; P0 approximately 1e-13-1e-12 yr^-1 over 1e7-1e8 yr, or use source-only
transient oracles. R1 selects/implements no physical rate, lifetime, coupling or trajectory.
The faster P0 examples still require depletion/coefficient checks; mere small Gamma times age
does not waive them. **THE CURRENT P0 TIER IS PRIMARILY AN ETA-EVOLUTION / MOVING-REFERENCE /
TRANSIENT ORACLE.** At 10^8 K and these rates, reachable beta thermal power is approximately
171/17 times below the omitted +T s_,n R_a term (R2 E-18), and far below the current thermal
sign-resolution floor. It is NOT a resolved thermal-sign experiment. A future purely
mathematical P0 thermal-sign stress test must include finite-T terms and remeasure resolution;
no new physical rate is chosen. Apply source-only tests before relaxation or after validated
freeze-out; the large-|xi| build estimate cannot rule out hot linear QSS.
A charge-neutral paired proton/electron removal is a later source-direction test, not the first toy.

| G item | Derivation / independently checkable target | R1 status |
|---|---|---|
| G1 | R1-R3; heat flux, diffusion and viscous ownership | DERIVED; R0 diffusion/bulk shorthand corrected |
| G2 | R5, clock+energy factors; Thorne/FR05 reductions | DERIVED for declared static/quasistatic scope |
| G3 | R7-R8; MPR three-piece sum; actual-potential power | DERIVED; supporting preprint independently checked |
| G4 | R11-R12; b^T(S-tBdot)=0; S=0 standard limit | DERIVED; qualified Z+t seam |
| G5 | R14; physical slide must reject old k-route | DERIVED / reproduced numeric negative oracle |
| G6 | R16 first-order gradient and delta mu_n | DERIVED; quadratic changing-Z extension stated |
| G7 | R18-R20; actual/eq representations, no unowned work | DERIVED for spin-off controlled closure |
| G8 | R17a/R17b, state function not heat | DERIVED; substituted metric sign corrected from R0 |
| G9 | R21 minima/roots/ratios; B1 and conditional B2 | Independently recomputed; B2 scope restricted |
| G10 | Neutron drive/capture signs and numbers | Reproduced from governed t/Z inputs |
| G11 | B_i/B_B, error, sum rules; chart/boundary | DERIVED / reproduced; rotating/subdomain future conditional |
| G12 | R23-R24 cold/static neutral first law | DERIVED; rotating extension conditional, quotable theorem source-limited |
| G13 | R15 projection; true specialization | DERIVED; G_true TOV numbers attributed to R0, not regenerated |
| G14 | R22 QSS and transient/reachability scaling | DERIVED; build-time numbers R0 diagnostic evidence |
| G15 | R4, Sommerfeld and compression corrections | DERIVED; eV scratch check |
| G16 | R9-R10 actual-potential bracket and sea average | DERIVED under declared event conditions |
| G17 | R25 accumulation cancellation, criteria/floor | Conditional model derivation, no hidden-sector rate ratification |
| G18 | R26 matched-state energy/lag/sign floors | DERIVED; future solver floor must be measured |
| G19 | Section4 units/redshift table | Independently dimension-checked |
| G20 | R-a/R-b/R-c; H1-H10 refusal table | DERIVED and source/notes audited |

Mandatory future tests additionally include zero-source bit identity with Phase-5D; physical
slide S=t dot B; t baryon/charge sum rules; initial neutron eta<0,R<0,eta R>=0; no-reaction
eta=eta0-Z int sigma dt; reached QSS R=-sigma; identical direct power under all three allowed
representations; chemical storage excluded from thermal RHS; P0 current-potential rule; P1
analytic average; finite-interval R20/R26 residuals; variable-Z scalar counterexample; product
fate/charge closure refusals. No test code or governed baseline is created in R1.

## 15. Proposed architecture and invariant ownership

Future layers: process identity/charge-closed stoichiometry; state-dependent local particle
source; external inflow of energy/charges/angular momentum; separate
energy partition; mandatory product fate and X reservoir; global source/energy integration;
qualified structural t; governed Z with sigma coupling; chemical-state diagnostics; unchanged
standard beta response; thermal ledger; matched control; diagnostics/sign map. Reuse generic
Phase-5D machinery after review. No reaction-specific conditions belong in generic execution.
For the BNV response, G_y builds Z and optional k diagnostics; its Phase-5 scoped Cowling
ownership is intact. No raw-G baryon drive or individual potential
reconstruction. No new production class is claimed to exist.

Candidate BNV-14 through BNV-25 replace incompatible old draft BNV-2/BNV-3 architecture clauses;
none is added to SCIENTIFIC_INVARIANTS or marked accepted:

| Candidate | Proposed statement |
|---|---|
| BNV-14 | Cold direct bracket uses actual potentials and stated event/fate assumptions; equilibrium upper bound conditional; each physical channel redeclares its finite-T weighting floor (smooth-toy eV is not universal), including omitted state-derivative terms |
| BNV-15 | Moving-reference delta N is baryon-neutral; physical forcing S-t dot B; physical sliding null |
| BNV-16 | t derives from the qualified equilibrium-sequence owner on the SAME declared domain D as its paired G_y/Z/reaction quantities, with chart/error/currentness and sum rules; long digits are governed-byte arithmetic oracles, not physical precision |
| BNV-17 | Explicitly narrows ADR-0013 Q1 (`docs/adr/ADR-0013-corrected-rotochemical-chemical-coefficients.md:212`) and ADR-0014 section 3.17 (`docs/adr/ADR-0014-secular-rotochemical-evolution-contract.md:472`–`:478`) upon owner ratification as proposed in ADR-0015 section 7: raw S_y is not two-channel, but Sigma_y=S_y-t Bdot=L sigma identically; physical seam {Bdot,sigma}; G_y retains Phase-5 coefficient/Cowling ownership and k diagnostics, never the physical raw-BNV map; Phase-5 mathematics/baselines intact |
| BNV-18 | First-order individual-potential identity R16; controlled ledger has no free reference-work input; higher-order state derivatives explicit |
| BNV-19 | E_chem is state storage; R17a/R17b are equivalent; never a heat addend |
| BNV-20 | Exactly one direct-energy representation; no duplicated hole/deposit/storage/mechanical/reaction-viscous/product energy |
| BNV-21 | Product-fate flag and per-species mechanical/chemical/thermal validity checks mandatory |
| BNV-22 | Roots/minima/B1/conditional B2/ratio bounds belong to the governed NON-SUPERFLUID R1995/FR2005 polynomials; every QSS statement requires tau_relax(T,xi) << relevant evolution time IN THE APPLICABLE REGIME, with nonzero normalization in each driven channel for R=-sigma |
| BNV-23 | One common rest-mass-inclusive energy zero; E_esc,fluid=E_esc,star+E_X distinguishes fluid exit from star escape/retained storage; P0 uses actual current local mu_n; cold over-extraction requires explicit external/reservoir owner |
| BNV-24 | R2006 Appendix/footnote 4 justifies Cowling for baryon-conserving perturbations; the Phase-6 baryon-direction mismatch is PHYSICAL-MODEL uncertainty/scope, not `numerical_error`; Phase-5C coefficient mathematics, scoped Cowling contract and baseline unchanged; forward BNV seam narrowing is explicit in BNV-17 |
| BNV-25 | PASSIVE SAME-INITIAL-B0 NO-BNV STAR comparator; matched energy identity; every sign names observable/time-window and resolution budget; distinguish tau_energy=C_*T/L from linear-response tau_th=C_*T/(a_cool L) |

Retained principles from old candidates: distinct particle/energy sources; once-only escape;
computed net sign; charge closure; evidence-class separation. They do not rescue the retired
raw-source map or fixed-reference baryon storage.

## 16. Scientific diagnostic evidence and R0 corrections

The companion note derives R15 and separates the reported off-equilibrium TOV evidence from
fresh R1 baseline arithmetic. R0 reports G_true symmetric at approximately 1e-6 (one independent
rerun1.4e-5), baryon direction near t (not agreement within its propagated budget), b^TG_true b approximately-4.66e54 count/MeV,
mu_B' approximately-2.147e-55 MeV/count, and eigenvalues of both signs. Cowling's positive
baryon curvature cannot be used for true total-B response. The fixed-baryon Z discrepancy is
much smaller on diagonals: about 1.9%/0.9%; true off-diagonal approximately 2.34e-55 versus
Cowling5.17e-55 MeV/count; BNV drives/W differ only a few percent on this fixture. Reported93%
of the off-diagonal difference arises from the baryon projection. These statements have R0
provenance and await final independent review, not a new production coefficient.

R1 does not treat R0 as infallible. Resolved corrections/qualifications are explicit:

| R0 statement | R1 adjudication |
|---|---|
| Substituted Echem_dot retains negative half-Zdot | Incorrect algebra; R17a negative, R17b positive; independent finite differences agree |
| E13 exact for any variable quadratic model | First order; add (1/2)b ell^T Z_,B ell at quadratic order |
| Diffusion boundary-only; bulk viscosity equals L_H-DeltaLnu | Gradient entropy production can remain; effective reaction viscosity overlaps chemical dissipation, radiation separate |
| Two-lapse expression has no geometry qualification | Static/quasistatic Killing/slice assumptions required |
| Injection Eq 3.62 proves arbitrary global multispecies variations | Restricted same-composition injection; R23 supplies static independent derivation; rotating theorem conditional |
| Equilibrium-mass bound universal | Actual-potential bound primary; equilibrium replacement needs equilibrium or demonstrated sink sign |
| B2/source cooling scale universal instantaneous bound | Requires capture/filling/history assumptions; B1 is unconditional within the polynomial model |
| All MeV channels heat at all times | No such arbitrary-history theorem; generic direct-dominated controlled comparison retained |
| Goldman Eq 9 has an extra4pi/3 | Not corroborated: N_e,N_e' are counts and collision rate includes division by volume; source's rough density/count conversion is the ambiguity |
| Author-note title scientifically wrong | Retired as this contract's organizing headline; a title's emphasis alone is not a false equation |

These adjudications harden the proposal rather than preserve either predecessor's mistakes.
They introduce no unresolved contradiction in the spin-off controlled ledger.

## 17. Published-source audit and source notes

B22/G24/A24 pages are PDF=printed; Z24 printed pages are one less than PDF after cover.
The four catalogued papers' relevant coverage remains distinct from the author's later notes:

| Source | Published process/energy content | Comparison with complete ledger |
|---|---|---|
| B22 Eq 4 p 6; Eqs12-18 pp 8-9 | Constrained chemical equilibrium; direction-specific Urca estimates | Not signed net rates without forward-minus-reverse. Beta disequilibrium is explicitly recognized; no complete eta/T ledger calculated |
| B22 Eq 21 p 10; pp 10-11 | Neutrino opacity; charged/photon cascade including ambient-electron annihilation | Full cascade/transport ownership needed; no independent universal deposition fraction |
| B22 Eqs22-28 pp 11-12, Eqs32-33 p 16 | Quasistatic sequence and mass/rotation changes | Mechanical state energy, not an independent thermal addend |
| B22 Eqs42-43 p 19; pp 30-31 | n->chi gamma, n->3chi; retained-energy heating estimates; direct/weak-restoration/thermal neutrinos distinguished | Published approximations; inclusive energy partition and beta coupling must be specified |
| B22 p 31 invisible n->3chi statement | Says channel evades heating constraints | Invisibility alone does not establish evasion; sub-Fermi residual can retain heat. Whether the cited n->3chi model weakens or evades a particular constraint remains source-limited without its primary source and event weighting |
| G24 pp 11-12 section 2.3 | Neutrino emission, thermal energy, work; detailed depletion chain needed; chemical/hydro timescale regimes | Recognizes these owners; quantitative complete partition not calculated, no inferred intention to exclude it |
| G24 Eqs8-9 p 17; Eq 10 p 21; p 23 | Occupied canonical energy versus invariant CM energy; sequence mass/rotation; escaping terminal products | Energy conventions must be translated; Eq 10 is not heat |
| Z24 Eqs3.2-3.6 printed 4-5 | Rotational-energy/magnetodipole accounting | Not BNV thermal power |
| Z24 Eq 3.7/text and3.8 printed 6/PDF 7 | Q_chi=Gamma n Ebar for angular-momentum transport; equilibrium MU T^8 estimate; prose notes beta enhancement | Emission is a loss owner. Equilibrium approximation is not a completed enhanced-neutrino thermal model |
| A24 pp 5-9, Eqs24-25 | B->psi gamma then psi n->pi- K+ removes two baryons per complete chain; occupied energies, density rate, mu symbol is mass ratio | Particle chemistry and energy partition separate; charge closure and full chain required |
| A24 EqD3 p 19; D5 p 20 | One-lapse number rate; sequence mass/rotation response | Number convention recovered; mass derivative is not heat |
| A24 EqE1-E3 p 20 | Escape-frame/threshold notation and boosted energy | Correct static threshold is e^Phi E_local>=m c^2 (Thorne3.3). E1 prints local threshold as infinity threshold; E3 actual-energy vs threshold notation unresolved |
| A24 EqE4-E5 p 21 | Escape rate called per baryon but density-rate dimensions; number-weighted escape ratio | Clarify density normalization; E5 is number, not energy escape fraction. No verdict about underlying numerical plots |
| A24 Eq 28/adjacent p 12 | Cosmological21-cm visible-energy estimate | Different system; not neutron-star deposition authority |

Published omissions are classified only as explicit approximation, explicitly identified need
for detailed treatment, or not addressed. No absence establishes author intent. The four papers
contain no standalone named hole formula and no demonstrated additive hole double count.
The new MPR supporting evidence and Thermo_BNV microscopic formula supply that interpretation
separately; they are not retroactively attributed to a Zakeri publication.

### MPR21 supporting source note

PDFp 2 immediately preceding Eq 9 gives the three-piece identity in R8; Eq 9 rate-weights
`mu_tilde_n-p_n'^2/(2m)` under equal masses, cold occupations, no n' blocking and final-spectator
Pauli restrictions. Eq 10 p 3 lacks explicit GR lapse/proper-volume factors and must be translated.
The report's approximate n' speed/binding is not literally stated there; binding must be an
independent channel conclusion or sourced to G22, not misquoted as MPR prose.

Eq 8 p 2 visibly prints `1/(1.2e11 yr)`. A crude consistency check against its volume averages
and Eq 11 implies about 5.74e-7 MeV/event, inconsistent in scale with its stated tens-MeV Fermi
energies. Replacing11 by19 would yield approximately 57 MeV, but is **not an authenticated
repair**. The J2144 ceiling differs by roughly an order of magnitude (approximately 8–10)
from this proposed 10^19 reading; Cas A is consistent with that reading (R2 E-14).
Classification: SOURCE NOTE / POSSIBLE TYPO, normalization unresolved; no repaired rate,
constraint or prefactor is installed. The PRL version remains absent.

### Goldman22 supporting source note

PDFp 2 assigns approximately 30 MeV to ordinary Fermi-sea relaxation and describes bound mirror
products; pp 3-4 introduces accumulated e'/D' and transfer to invisible radiation. These are
successive energy owners, not a universal gravitational heat term. Eq 8 p 4 compares transfer
power to luminosity; blackbody capacity does not itself prove sufficient energy supply.
Literal cooling to zero is not established by a coupled solution with persistent deposition.

R1 corrects/extends the R0 normalization notes:

- Eq 9 p 4 `collision rate=c f f' N_e N_e' sigma/V`, V=4pi R_c^3/3, is correct for **counts**.
  R0's claimed extra volume denominator is not corroborated; preceding rough count-density
  conversion omits an order-one factor.
- Eq 6 p 3 prints an inverse `(m_e' c)^3` prefactor; substitution p=m_ec X into its own S2
  requires a positive cubic factor. The h versus hbar state-counting convention also needs
  resolution before quantitative density use.
- S8/S9, after the stated effective inertial-mass replacement, give `(16e^2/9)E_F T`;
  printed S10 drops e^2. This is an internal algebraic source note, not a repaired screening law.
- The main/supplement coupling scales have rounding/convention issues; the transport/screening
  model needs primary authority. R1 does not endorse R0's broad quantitative normalization
  closure or any cross-section/coupling constraint.

Qualitative accumulated-sector thermal transfer is supported; quantitative Regime-II
predictions remain source-limited. None of these source files was changed.

## 18. AUTHOR WORKING NOTES — UNPUBLISHED / NON-AUTHORITATIVE

Thermo_BNV has11 PDF pages, printed page=PDF-1. Both whiteboards and critical displayed
formulas were visually checked. Classifications concern supplied equations, never the author's
unprovided intentions or a claim that these notes exhaust the intended theory.

| Note locus | Classification | Evidence / required interpretation |
|---|---|---|
| (1.1), printed 1 | Independently confirmed | Product rule for entropy with s per baryon; not a universal external-heat identity |
| (1.2)-(1.3) | Correct but incomplete | Full thermodynamic chemical potentials and declared rest-energy convention required |
| (1.4) | Independently confirmed | Adiabatic/composition-fixed identity |
| (1.5) | RETIRE as a physical source | Euler numerator vanishes for consistent definitions; epsilon intensive. A vacuous0=0 is not categorically a false identity, but cannot create BNV heat |
| (1.6), printed 1 | Correct but incomplete / overlap hazard | Enthalpy-per-baryon source lacks explicit escape; finite-T Ts term is normalization bookkeeping; material-coordinate interpretation required |
| (2.1)-(2.6), printed 2-4 | Compatible under stated approximations | Standard redshift/diffusion/chain rule; governed ADR-0013 replaces any old singular full-matrix inversion |
| Procedure1(b)ii | Incomplete | Frozen-background validity needs integrated depletion and coefficient changes, not a rate alone |
| Procedure1(c) | Incomplete if actual sources omit BNV | “Various reactions” must include S_y; pure neutron lepton-reduction can hide omission; no author-intent inference |
| (2.7)-(2.8), (2.10) | Restricted / incomplete | Neutron-only specialization; t_l=Y_l+B dY_l/dB, not automatically Y_l. General zero drive is 2Omega dotOmega I_Omega,l=+sigma_l |
| (2.9) | Conditional source formula | Rotochemical asymptotic assumptions and reachability remain required |
| (2.11)-(2.12) | Restricted beta-mediated estimate | Not a generic total-BNV thermal or rate bound; direct partition and reached QSS needed |
| (3.1) first line, printed 5 | Independently confirmed under stated energy reading | Read E*_cm as the canonical CM energy m*sqrt(A), with A=1+sigma^2+2sigma x, so its dilation factor reproduces the authenticated A24 expression; then the hole factor E_F-E with the same event weighting plus deposited photons equals complete mu-E_chi |
| (3.1) second line | Inconsistent rate weight / typo-level | With A=1+sigma^2+2sigma x, printed x(1+sigma x)A must be x A^2 to agree with its own3.11 and A24 Eq 24; ratio A/(1+sigma x). No claim about plotted implementation |
| (3.2) | Independently confirmed | Isotropic phase-space integral definition |
| (3.3) | Effective-mass NOTATION defect (not dimensional) | Prints m_B^3; x=E*/m* requires (m_B*)^3 Jacobian, consistent own3.18 |
| (3.4)-(3.5) | Consistent definitions | Energy-weighted integrals, conditional on amplitude/normalization |
| (3.6) | Correctable frame/mass ambiguity | Specify CM energy and effective versus vacuum mass; no RMF rate ratification |
| (3.7)-(3.9), (3.11)-(3.18) | Structure independently supported | Boost/kinematic chain conditional on definitions;3.10 is empty label,3.14 collinear not general boost;3.17 polynomial verified at exact-rational sample points |
| (3.19), (3.21)-(3.22), printed 7 | Independently confirmed structure | Two-lapse proper-volume luminosity; prose “per baryon” needs separate normalization |
| (3.20) | Direct-only approximation | Extend with governed beta heating and full neutrinos; if Lnu already includes enhancement do not subtract it again |
| Photon formula immediately below 3.20 | Additional dimensional defect | Printed Lgamma=4pi sigma_SB T_s^4 lacks R^2 |
| (3.23) | Envelope source-limited attribution | Supplied reference placeholder unresolved; R0 identifies GPE iron form, not FR05 fully accreted49; no envelope swap |
| delta m>0, printed 8-9 | Confirmed with comparison defined | Stable near-equilibrium fixed-conserved-state chemical storage; not arbitrary mass-loss comparison |
| delta mu growth at Gamma mu_B | Scale hypothesis, incomplete | Actual source drive is -Zsigma, not an insertion-energy estimate |
| “work on fluid” / title | Retired contract framing | No generic extra work heat; title emphasis alone is not a false scientific equation or proof of intended exclusion |
| Oct whiteboard | Inconsistent source coefficient | Distinguishes rest mass mu_B from full mu_bar but uses rest-only mu_B d(delta A) in total first law |
| Dec whiteboard | Incomplete; sign inconsistency | n=delta a/delta V product rule valid, no escape/dissipation law; (1-Gamma delta t) destruction conflicts with Gamma=Bdot/B unless redeclared |

Historical notes remain untouched. In particular the internally consistent hole+photon
representation is not itself double counting. The error is combining it with an already
inclusive source or using inconsistent rates/energy conventions.

| Hazard | Proposed refusal / precise classification |
|---|---|
| H1 HYPOTHETICAL combination of note (1.6) full retention + section 3 complete emissivity | Refuse same-boundary full-retention mu R_a plus (mu-E_chi)R_a; it yields 2mu-E_chi at T=0. The historical notes do NOT demonstrate this combination; section 3 explicitly neglects those chemical/work terms |
| H2 extra hole over complete residual | Same energy/different representation; refuse duplicate |
| H3 dot E_chem as heat | State change, not another positive thermal addend |
| H4 dot M_eff as heat | State/rotational work, not intrinsically dissipation; possible overlap |
| H5 rest-mass q plus rest-inclusive energy | Refuse inconsistent energy-zero accounting |
| H6 unequal hole/photon event weights | Inconsistent partition, not necessarily an additive duplicate; common event identity/rate required |
| H7 generic work heat | No independent irreversible mechanism; microscopic hole “work” already included |
| H8 bound-product energy twice | Explicit outbound/inbound transfers are legitimate; duplicate positive heat attribution without reservoir is refused |
| H9 reaction bulk viscosity again | Refuse duplicate dissipation from the same resolved weak reactions; radiative increment remains separate |
| H10 retired fixed-reference storage/work | REJECTED as primary moving-star physics; old alpha^2/(2a), E_2, free F_ref and Delta(delta g^T F) survive only in predecessor history, not as physical flows |

## 19. Remaining source limitations and final review gate

Not read because absent as authenticated primary sources: Bardeen1970; Bardeen-Carter-Hawking1973;
Thorne1977 ApJ (distinct from supplied lecture); Misner-Sharp1964; Alford et al2018; Haensel1992;
MPR2020 and the final PRL MPR21 version; Strumia2022; Berezhiani2021; GMN2019; canonical-RMF
textbook authority; supporting transport/superfluid/two-fluid sources listed in R0 Q49.
Some review material exists in `_incoming`; it is not silently substituted for a missing
primary. No missing-paper quotation is used to support a core equation.

Load-bearing controlled equations R1-R12,R16-R21,R23-R26 are derived with stated assumptions;
core readiness is not blocked merely by lack of a quotable general theorem. Material future
limits: actual channel kinematics/escape and A24 AppendixE; MPR prefactor; Goldman transport
normalizations and screening; realistic-EOS t and Cowling error; superfluidity; global rotating
variational proof; changing-background second-order potentials, heat capacity and boundaries;
Regime-II dynamics; whole-star rate-weighted hole energy. R0's G_true full matrix/scratch is
not supplied, so its numerical diagnostic is attributed rather than regenerated/ratified.

**R3 disposition A:** PHASE-6A-0 FINAL-REVIEW CORRECTIONS APPLIED — MATERIAL GOVERNANCE
CONFLICT EXPLICITLY RECONCILED — READY FOR BOUNDED FRESH-CONTEXT RE-REVIEW.
ADR-0015 remains PROPOSED / NOT ACCEPTED / NOT OWNER-RATIFIED / NOT CANONICALLY INTEGRATED.
BNV implementation is NOT BEGUN. The R2 review passed the load-bearing physics; this correction
pass does not itself certify a clean independent re-review. No BNV rate, trajectory, n->chi gamma
implementation, Regime-II/MixedStar implementation, realistic A18 work or canonical merge.

Exact next action, **not executed** (the commit alias is reported with its resolved hash outside
this document to avoid a self-referential commit):

> Run a BOUNDED fresh-context Claude Opus 5 XHIGH re-review of ONLY the delta from
> 7a31862f1e8a5cf046316882e05d76e4924e27d9 to
> PHASE6A0_R3_CORRECTION_SHA, plus the immediately surrounding affected
> paragraphs.
>
> The reviewer must verify:
>
> - E-1 governance reconciliation is explicit and semantically consistent;
> - all E-2...E-27 corrections are applied without changing load-bearing physics;
> - ADR-0013/0014 coefficient mathematics and Phase-5 baselines remain
>   untouched;
> - no new contradiction or accidental ratification has been introduced.
>
> It does NOT need to rederive the entire BNV framework again unless the delta
> changed a load-bearing equation.
>
> If that bounded re-review passes with zero BLOCKING/MATERIAL findings, return
> to the owner for explicit ratification of PROPOSED ADR-0015, including explicit
> acknowledgment of the narrowing of ADR-0013 Q1 / ADR-0014 §3.17's
> forward-looking BNV seam.
>
> Do not begin the re-review or ratification automatically.

## 20. R1 validation record

Permanent allowlist: this preflight, PROPOSED ADR-0015, and the dedicated Cowling diagnostic.
No status-link changes were needed. Entry hashes cover 2898 tracked paths (including the two
old drafts), all 43 literature files and external/author inputs identified above. All11 governed
baselines,33 protected Phase-5 source paths, historical candidate artifacts and EOS/data remain
unchanged; final comparison excludes only the two authorized draft edits and new diagnostic.

R1 scratch checks: independent t/error/source-map/physical-null arithmetic; projection identity
(relative residual approximately 2e-16); direct actual/eq identity; scalar and matrix variable-Z
finite differences;60-digit beta roots/minima/B2; relativistic sea-average and eV correction.
One root scratch attempt named its script math.py and shadowed the standard module; renaming
the scratch script resolved the import, without repository or dependency edits. All final
scratch checks passed. No new governed baseline or production trajectory was generated.

Existing lightweight standard checks reran5/5 successfully: phase5d_response,
phase5d_independent_oracles, phase5d_component_tolerances, phase5d_protected_manifest,
phase5d_harness_controls. The existing scratch build uses the same unchanged production inputs
and miniforge Python with mpmath. No multi-hour trajectory suite was rerun.

Final checks passed: exact documentation-only diff, companion Markdown links and repository path/line
references, fenced-block integrity, G1-G20/BNV14-25 coverage, stale-concept labeling, immutable
entry hashes, git diff --check and staged diff --check. Commit is a new descendant with exact
message `docs: harden bnv thermal first-law preflight`, never amend/squash. Resulting SHA and
local/upstream/live equality are reported outside the commit to avoid self-reference. Push is
non-force; fetch and live master authentication follow. Canonical master remains the entry SHA.


## 21. R3 final-review correction record and validation

This is the exact E-1–E-27 reconciliation checklist from R2 final-report item 108 / Reviewer E
Part V, under the user's Phase-6A-0R3 instructions. Reviewer E is the final overlap adjudicator;
the review is evidence, not ratification. R2 reported 0 BLOCKING, 1 MATERIAL, 26 NONBLOCKING
and 21 NOTE findings. R3 records the writer's applied corrections, not a new reviewer verdict.

### Review-input provenance

Review root: `/private/tmp/claude-501/-Users-keeper-Documents-CompactStar/a0872564-d9f7-4ca8-80ef-b0756a592e40/scratchpad`.

| Review input relative to that root | SHA-256 | R3 use |
|---|---|---|
| PHASE6A0R2_FINAL_REVIEW_REPORT.md | 20996d48959c302f9f70f28635078004e6da65751d98191a6fc7dcbc5b6a0d14 | Read completely, 158 lines; exact consolidated checklist |
| reports/reviewer_E.md | 93a4c8167b021e1d50c7aa4d99d73e029ad0fb66006eb5ab10ccaaa89d0e673e | Read completely, 360 lines; final adjudication and precision evidence |
| reports/reviewer_A.md | e2d86abb0dd69bc4c4824871c051ff693f4dd49377b6f57e0ac7174faf33b87b | Fingerprinted; E's adjudication sufficed; no separate full-read claim |
| reports/reviewer_B.md | e866e6f151053a6d62da5ee987026e172dca266293ea95f84bdb50988ee21898 | Fingerprinted; E's adjudication sufficed; no separate full-read claim |
| reports/reviewer_C.md | 76512bf55e4f3cd1fff7cce620a3d67c4373bccbf5e93d652ddda25f4227b111 | Fingerprinted; E's adjudication sufficed; no separate full-read claim |
| reports/reviewer_D.md | e77341a067ba828865a3056aa6828a2d920f966a42b789f63d6effbea1f11bb7 | Fingerprinted; E's adjudication sufficed; no separate full-read claim |

The scratch paths are exact audit provenance, not a claim of permanent publication or catalog
promotion. R2 supplied the new qualification scales cited in R3; no new physical-rate model or
whole-star trajectory was computed. Published claims, author working notes and independent
preflight derivation remain separate evidence classes.

### Mandatory corrections — all applied

| ID | R2 class | Applied correction / inspection location |
|---|---|---|
| E-1 | MATERIAL | ADR sections 1/7/20, preflight section 7 and BNV-17 explicitly narrow ADR-0013 Q1:212 / ADR-0014 section 3.17:472–478 only for the forward BNV seam; raw S_y versus Sigma_y=L sigma; future owner acknowledgment mandatory |
| E-2 | NONBLOCKING | C_(MeV->erg), MeVToErg owner and once-only boundary explicit; A2/A6 and R7/R18–R20 thermal interfaces show conversion; chemical state remains MeV |
| E-3 | NONBLOCKING | ADR section 4, preflight section 3 and BNV-14: channel-weighting-specific finite-T scale; smooth P0/P1 versus Fermi-surface O(k_BT); NR-only compression coefficient |
| E-4 | NONBLOCKING | R reserved for global beta vector; local directed BNV density R_a, including P0 and microscopic/hazard formulas |
| E-5 | NONBLOCKING | Unified E_esc,fluid / E_esc,star / E_X and infinity L_out,fluid / L_esc,star / L_esc,rest,star; removed undefined aliases; infinity convention explicit |
| E-6 | NONBLOCKING | Preflight section 4 labels (-1,+1,+1) as (n,p,l) and gives the governed y-basis L column |
| E-7 | NONBLOCKING | Preflight section 3 identifies chain rule applied to Gibbs |
| E-8 | NONBLOCKING | Preflight section 6 and diagnostic section 2 label long digits arithmetic oracles, propagate t budget and require Z/model budgets; drive coefficients consistently MeV/count per abs(Bdot) |
| E-9 | NONBLOCKING | ADR section 7, preflight sections 6/7/12 and diagnostic: Structure-1 free-gas qualifier, fixture-specific muon drift approximately 54 and depletion/sensitivity requirement |
| E-10 | NONBLOCKING | Preflight section 5: binding belongs to independent kinematics/G22, not MPR prose |
| E-11 | NONBLOCKING | Cold-bound and PROMPT_ESCAPE definitions include e^Phi E_local>=mc^2; pointwise upper bound includes e^-Phi rest threshold |
| E-12 | NONBLOCKING | Thermo_BNV (3.3) classified effective-mass notation defect, not dimensional; photon formula dimensional defect remains separate |
| E-13 | NONBLOCKING | Thermo_BNV (3.1) first line requires canonical CM energy E*_cm=m*sqrt(A) to reproduce authenticated A24 weighting |
| E-14 | NONBLOCKING | MPR source note: J2144 differs by roughly an order of magnitude from possible 10^19 reading; Cas A consistent; no repaired prefactor installed |
| E-15 | NONBLOCKING | H1 explicitly HYPOTHETICAL combination; no claim the historical notes performed it |
| E-16 | NONBLOCKING | ADR section 5 / preflight R-b now multi-species removed/created/deposition expression with one rate and common energy zero |
| E-17 | NONBLOCKING | Large-xi modified-Urca build table qualified; 0.43/1.15/4.27 asymptotic xi at 10^8 K; linear scales 3.5e4/2.9e4 yr, T^-6; hot QSS/cooling/freeze-out/transient sequence and applicable-regime criterion |
| E-18 | NONBLOCKING | Current P0 tier is eta/moving-reference/transient oracle, not resolved thermal sign; finite-T contribution and resolution dominate its reachable beta signal |
| E-19 | NONBLOCKING | Class F external incoming energy/charges/angular momentum; state-dependent source callback with raw/moving/reduced outputs distinguished |
| E-20 | NONBLOCKING | X-mediated ordinary weak kinematic/normalization feedback separately budgeted, including DeltaB/B->0; corresponding Regime-II owner |
| E-21 | NONBLOCKING | R25 is linearized fixed-f' toy; Goldman's actual accumulated degenerate fluid is not rate-independent |
| E-22 | NONBLOCKING | Passive same-initial-B0 no-BNV comparator; every sign names observable/time-window; future surface-temperature interpretation separates structure |
| E-23 | NONBLOCKING | ADR section 10, preflight section 10 and BNV-22 explicitly NON-SUPERFLUID R1995/FR2005 polynomial roots/minima/B1/B2/ratio bounds |
| E-24 | NONBLOCKING | Floor solves P_dir=L_std+L_X; hidden-only power-law expression is an upper bound |
| E-25 | NONBLOCKING | 1.4e5–1e6 yr attached to C_*T/L; linear response C_*T/(a_cool L) separate; free-gas C_* / placeholder Ltilde dependency stated |
| E-26 | NONBLOCKING | Cold/frozen Udot and Eeq_dot omit finite-T B derivatives and mu_B correction; future exact conservation oracle must own them |
| E-27 | NONBLOCKING | BNV-16 same domain D; BNV-23 common rest-inclusive zero and fluid/star/X split plus actual-local-mu P0; BNV-24 R2006 baryon-conserving Cowling basis and physical-model error class |

### NOTE-class disposition without scope expansion

E-28/E-29: clarified beta-disequilibrated hydrostatic/diffusive R24 and its direct closure with
R5 and actual potentials. E-30: finite-p_cut future oracle must declare surface work.
E-31: affected typography cleaned. E-32/E-33: diagnostic names the Q identity and R2006
baryon-conserving Cowling basis; D2 symmetry requirement distinguished. E-34: "near t" is not
within-budget agreement. E-35: V=2E_chem explicit. E-36/E-37/E-38: passed null/charge/variable-Z
statements preserved. E-39: absorption elsewhere in ordinary matter is SM_THERMALIZATION.
E-40: finite-T thermal-tail weighting added to sign-map axes. E-41: nonnegative J_X explicitly
tightens the upper bound. E-42: extra Goldman S9 count-density slip deferred; no new published
source claim needed for this bounded correction. E-43: nonzero normalization per driven QSS
channel stated. E-44: ADR Status key and diagnostic status clarified. E-45: pre-existing Phase-5B
baseline status-string housekeeping deferred; no baseline metadata edit. E-46: revised R0 P0
tier retained with corrected purpose. E-47: status-banner reconciliation belongs to a future
explicitly authorized owner-ratification task; implementation remains NOT BEGUN. E-48: passed
physics preserved, not reratified.

### R3 verification boundary

R3 entry fingerprints cover 2899 tracked files and all 43 literature files. The only allowed
changes are these three Phase-6 documents. Every one of the other 2896 tracked files is
compared by SHA-256, including accepted ADR-0013/0014, production/test code, CMake, all 11
governed baselines, Phase-5B/C/D artifacts, historical candidates and tracked EOS/data.
The review inputs and literature bytes are also compared to their R3 entry fingerprints.
No scientific trajectory suite is required for this documentation-only correction.

Final R3 results: PASS — exact three-file documentation allowlist against both reviewed entry
and canonical master; all 2896 protected tracked paths and 43 literature files unchanged;
all six review inputs unchanged; 22/22 literature manifest entries authenticated; all 11 governed
baselines unchanged. The existing `tests/rotochemical/manifest.py:46` generator/checker passed
33 protected paths, ten predecessor baselines and three special artifacts; the eleventh
Phase-5D baseline separately matches `2606916915b2da5c051b1a06637c0a371a74751c6e77b2775bc629a63bc9f6dd`.

Documentation checks passed: all 27 checklist rows, companion links, repository path/line
references, Markdown fences/table columns and stale-alias searches. Of 34 numbered equation
blocks, 26 remain byte-identical. A2/A6/R6/R7/R9/R18/R19/R20 change only the declared conversion,
local-rate or infinity/fluid/star labels; manual comparison and an exact-rational thermal-boundary
cancellation check preserve their physical identities. Chemical matrices, signs, beta polynomials,
roots/minima and numerical fixture coefficients are unchanged. No trajectory suite or new
production test was run; this pass uses byte protection and documentation/scratch checks.

Two validation-command mistakes were corrected without touching protected data: the scratch
baseline counter initially counted only the three JSON files, omitting eight TSV baselines;
and the first manifest command ran inside literature although its paths are parent-relative.
The corrected all-format count and `shasum -a 256 -c literature/SHA256SUMS.txt` from the
CompactStar parent both passed. No hash mismatch or scientific test failure occurred.
`git diff --check` passed; the staged version is required before commit.

The correction commit uses
`docs: reconcile bnv thermal preflight review`, a new descendant of PHASE6A0_R2_REVIEWED_SHA.
Resulting PHASE6A0_R3_CORRECTION_SHA, push/fetch equality and canonical parity are reported
outside the commit; no amend, squash, force push or merge is authorized.
