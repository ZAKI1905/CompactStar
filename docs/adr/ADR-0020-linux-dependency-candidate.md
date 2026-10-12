# ADR-0020: Explicit native Linux dependency candidate

Status: PROPOSED FOR INDEPENDENT REVIEW; implementation authorized by the owner's
2026-10-10 Linux-transition task and 2026-10-11 correction task, not scientifically
ratified or integrated.

Dependency/build and structural changes extend ADR-0019's explicit external
package authority to native Linux. Darwin retains the byte-identical Mac lock.
Linux must select an absolute `COMPACTSTAR_LINUX_DEPENDENCY_LOCK`; absent locks,
wrong platform/mode/source/compiler/flags, changed package bytes or changed
GSL/zlib/runtime/toolchain closure are fatal before package CONFIG files execute.
Cross compilation and implicit package fallback are unsupported. Exact Zaki
2.0.1 and CONFIND 2.0.0 source identities remain those accepted by ADR-0019.

A lock is an integrity record for a candidate build, not a signature, numerical
acceptance, or campaign authorization. The versioned isolated prefixes are
reviewed with this candidate; changing a prefix or its content requires new
provenance and requalification. Absolute paths intentionally bind this cluster
installation. A relocatable design would need explicit reauthentication.

Alternatives: overwriting the Mac lock would destroy accepted authority; accepting
any package with a matching version would omit source and arithmetic provenance;
vendoring Linux binaries would blur external dependency ownership. Retain the
separate externally built packages and authenticate before loading CMake code.

`Tags.hpp` only qualifies existing enum names to honor C++17. Numerical equations,
constants, baselines and failure records are unchanged. Linux harness adapters
retain the historical Mac verification paths; they must not counterfeit Mac
platform evidence. Current proposed cross-platform scoring follows the coherent
qualification amendment below. `docs/validation/linux-transition/PREDECLARATION.md`
remains the historical starting plan, not a competing current scoring contract.
A combined independent
scientific/numerical review precedes any canonical integration or long campaign.

Candidate evidence tracked under `docs/validation/linux-transition` is restricted
to provenance JSON, inventories, text reports and compressed textual test logs.
Native archives, executables, scientific input trees and disposable builds remain
outside Git in the isolated workspace. Original failure logs are preserved.

## Proposed qualification amendment — 2026-10-11

Status remains PROPOSED, not owner-ratified. The owner explicitly authorizes this
implementation pass, static certificate construction and short qualification
pilots. No baseline promotion, Mac policy change or full BNV campaign is authorized.
This section supersedes numerical scoring suggestions in the prior Linux
PREDECLARATION and independent review section 9; it preserves their failure evidence.
Governance §§2–5, INV-14/15, ADR-0011/0012/0013/0014/0017/0019 remain the physical
and ownership authorities. Changes are numerical-method (PE-V6), engineering
(U11, source reads and compiler selection), dependency/build and documentation.

Five distinct claims must be reported:

A. Same platform/build/runtime: exact payload reproducibility across repeats;
CONFIND also exact worker-count invariance. Debug and Release are separate
platform baselines unless independently shown equal. Hexadecimal float spelling
is not a payload; decode IEEE binary64, retaining signed zero. Require all integer,
identity, ordering and topology fields exactly. A display normalization may never
change a bit-pattern field.

B. Cross-platform agreement: scalar functions use independent arbitrary precision
and operation/conditioning bounds. For the neutron onset retain the existing
3e-16 relative accuracy requirement at governed constants. Its near-neighbour
inverse-density diagnostics use 16 binary64 epsilon times |n|, from subtraction
conditioning (mu is well separated from mass), square root, cubing, division,
and rounded hbar-c/pi inputs; this is not an onset-gate relaxation. Libm direct
pow/log10 probes use one ULP absolute error at exact binary inputs (including
extra held-out inputs), never a promise of correct rounding. Cancellation-prone
quantities require absolute or scale-based errors. No blanket 10*rtol rule.
Structures and derived coefficients require resolution/refinement evidence,
validated error estimates and propagated uncertainty. Missing evidence is
UNQUALIFIED, even when drift is tiny. Historical same-grid bounds do not silently
become cross-grid tolerances. Coefficient agreement is assessed componentwise
against the sum of independently validated numerical error envelopes; this is
conditional on the existing refinement tests passing for the current producer.
Discrete layer indices and unambiguous topology remain exact. At near-ties,
inspect raw scalar accuracy and geometry; a different point count is a documented
conditioning result, not an automatic numerical pass or permission to redesign.

C. Independent accuracy: PE-V6 uses >=80 decimal digits from exact binary input
constants, with a second precision run, and the algebraic two-body momentum
formula. The reviewed 2026-09 literal and native long double are not sole oracles.
Additional neighbouring mass/mu inputs and held-out libm arguments guard against
fitting to observed differences. U11 constructs below/equal/above in the actual
geometric comparison units, using nextafter; no round-trip identity assumption.

D. Provenance: hashes of every protected source, input, profile, coefficient,
certificate and producer remain exact. The historical 33-path manifest remains
unchanged. A separately versioned upstream successor records each transition and
must keep all 33 paths, ten baselines and three special artifacts, with negative
mutation controls. Its passing result is explicitly candidate-only. This is the
separate upstream extension required by the frozen-context plan §8, not a waiver.
A native certificate must bind to its own native star; no tolerance bypass of the
Mac exact-profile guard. Record loaded libm/libc and dependency hashes, CPU flags,
compiler/effective FP flags, and float/double/long-double size and mantissa widths.
Test representative compute nodes, including actual arithmetic results.

E. Acceptance: Linux native goldens are candidates until owner acceptance after
independent review. Existing Mac-only exact tests and their original Linux failures
remain visible. A proposed N/A classification requires replacement scientific
coverage; absence of that coverage is UNQUALIFIED, never a pass-count gain.
Actual AppleClang verification is pending separately. Preserve the known Mac
implicit contraction policy; cluster double-width simulations are not Mac runs.

PE-V6 will use explicit fma in mu_e and binary64 inverse-density arithmetic,
preserving equations and constants. Other extended-precision sites are audited
per quantity, not mechanically rewritten. Record changed onset metadata and EOS
hashes from fresh producers. U11 changes fixtures only; FindTbIndex's outermost
rho >= threshold rule remains exact. CONFIND contour algorithms remain unchanged.

The legacy T2 lifetime sensitivity (2.2–4.3%, within-Mac 0.59–0.66%) remains
scientifically unresolved. Establish call-path independence before excluding it
from controlled Phase-6 scope; no whole-library equivalence claim is permitted.

Execution limits: Slurm for builds/solves; enumerate tests and inspect horizons.
Run the 13 existing readiness tests and static diagnostics. Short source/control/
event pilots require an authenticated native fixture; never attach Mac sensitivities
to another star. Candidate linear stability ceilings are not nonlinear validation.
1 Myr, 10 Gyr and 500 kyr campaigns are excluded. Record wall/CPU/memory costs.
