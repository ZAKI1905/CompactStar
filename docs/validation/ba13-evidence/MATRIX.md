# BA13 complete matrix

Full 21-by-22 data: [matrix.tsv.gz](evidence/matrix.tsv.gz) and [matrix.json.gz](evidence/matrix.json.gz).
Every row retains historical drift, threshold, utilization and assessment status. Unknown total justified error and physical headroom are null.
The table below shows the row with largest finite-added/headroom fraction for each quantity. Finite sums include unaccepted numerical and conditional estimates; the remaining column is partial headroom after those terms, not a complete physical enclosure.

| Quantity | Sample | Historical utilization | Finite added terms | Justified total added error | Partial remaining headroom | Status |
|---|---:|---:|---:|---:|---:|---|
| t_n | 20 | 0.00332093366 | 5.95022e-10 | UNRESOLVED | 9.62557e-05 | UNQUALIFIED |
| t_e | 20 | 0.062762656 | 4.7709e-10 | UNRESOLVED | 2.89109e-06 | UNQUALIFIED |
| t_mu | 20 | 0.376275898 | 3.46282e-10 | UNRESOLVED | 2.10324e-07 | UNQUALIFIED |
| N_n | 20 | 0.0097197128 | 4.01026e-12 | UNRESOLVED | 9.9028e-05 | UNQUALIFIED |
| N_e | 20 | 0.0488365517 | 8.72779e-11 | UNRESOLVED | 9.51163e-05 | UNQUALIFIED |
| N_mu | 20 | 0.54075348 | 4.8313e-09 | UNRESOLVED | 4.59198e-05 | UNQUALIFIED |
| Ltilde_Me | 20 | 0.0122667955 | 3.03723e-06 | UNRESOLVED | 9.57361e-05 | UNQUALIFIED |
| Ltilde_Mmu | 20 | 0.25681093 | 1.26309e-09 | UNRESOLVED | 7.43176e-05 | UNQUALIFIED |
| mu_B | 20 | 0.00188263921 | 2.39759e-12 | UNRESOLVED | 9.98117e-05 | UNQUALIFIED |
| P0_average | 20 | 0.00188263921 | 2.39759e-12 | UNRESOLVED | 9.98117e-05 | UNQUALIFIED |
| P1_average | 20 | 0.00261457001 | 3.47973e-12 | UNRESOLVED | 9.97385e-05 | UNQUALIFIED |
| P2_average | 0 | 0 | 0 | 0 | 0.0001 | PASS |
| radius | 20 | 0.0128306445 | 1.47617e-11 | UNRESOLVED | 9.87169e-05 | UNQUALIFIED |
| Z_row_npe | 20 | 0.0216301769 | 1.00279e-10 | UNRESOLVED | 9.78369e-05 | UNQUALIFIED |
| Z_row_npmu | 20 | 0.350922112 | 4.22268e-09 | UNRESOLVED | 6.49036e-05 | UNQUALIFIED |
| Cstar | 20 | 0.0145900852 | 1.1599e-08 | UNRESOLVED | 9.85294e-05 | UNQUALIFIED |
| envelope | 20 | 0.00975368197 | 6.48223e-08 | UNRESOLVED | 9.89598e-05 | UNQUALIFIED |
| Tsurface_inf | 20 | 0.00787060676 | 6.48195e-08 | UNRESOLVED | 9.91481e-05 | UNQUALIFIED |
| surface_gravity | 20 | 0.0370511514 | 6.86913e-11 | UNRESOLVED | 9.62948e-05 | UNQUALIFIED |
| mu_n_profile | 20 | 0.00188265152 | 5.80869e-11 | UNRESOLVED | 9.98117e-05 | UNQUALIFIED |
| metric_structure | 20 | 0.0330124185 | 1.26865e-06 | UNRESOLVED | 9.54301e-05 | UNQUALIFIED |
| species_support | 0 | 0 | 0 | UNRESOLVED | 1 | UNQUALIFIED |

For every individual sample, the full files list each finite term, evidence class, source hashes, accepted-authority field, missing remainder and conditional dependency.
All approximate X2 rows retain unresolved reference stage/arithmetic, unseen bias, finite-sensitivity remainders, moving-cut input and nonlinear fit/input terms. Interpolation, positive-B alignment, stencil/normalization and continuous curvature are explicit conditionals. The registry accepts none.
mu_n_profile additionally records the equilibrium-definition condition; metric_structure lacks R-Q nodewise evidence; species_support lacks continuous no-crossing acceptance. Only P2 is an exact closed identity.
