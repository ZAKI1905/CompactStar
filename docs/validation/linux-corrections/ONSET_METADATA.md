# Derived metadata transition

Constants remain the exact accepted Zaki binary64 inputs; the independently
computed neutron threshold is 7.35672890373283169196...e-9 fm^-3.
The corrected production result is its correctly rounded binary64 value:
`0x1.f98cf18cad246p-28` = 7.3567289037328318e-9 (17-digit rendering).
The prior Linux result was `0x1.f98cf18cad24dp-28`, seven ULP higher.
The error against the exact-input arbitrary-precision result is 8.43e-18 relative;
PE-V6's 3e-16 bound is unchanged. This is numerical-method, not new physics.

The inverse-density correction also changes the Linux muon threshold from
`0x1.d3f3d3280fa06p-2` = 0.45698480541241826 to
`0x1.d3f3d3280fa07p-2` = 0.45698480541241832 fm^-3 (one ULP).
It is still different from the historical Mac 0.45698480541241987.
No constant, threshold selection rule or expected Mac table is substituted.

Source: imported review §§4.3/5.1 and linux-production-probe.json. Provider domain
metadata is assembled from its live thresholds in
CompactStar/EOS/src/TrackRFreeGasThermodynamics.cpp:78, and the characterization
producer writes fresh model.txt/freegas.tsv/profile.tsv. The coefficient evidence
producer hashes those actual files (tests/analysis/chemical_production_evidence.py:177).
The static qualification job regenerates these artifacts and keeps their native
source/provenance strings; historical Mac tables are untouched. Actual regenerated
hashes and refinement status must be read from that job's receipt, never inferred
from this scalar agreement or overwritten with the Mac strings.

The table writer also contained hardcoded source_commit/provider_commit fields
from f4ae22d/933494d. They are now explicitly historical_generator_origin_commit
and historical_provider_origin_commit, with the candidate arithmetic contract
and actual source-authority location disclosed. The protected successor records
this seventh source transition. No historical provenance string is forged to
make a newly generated table appear to come from an old producer.
