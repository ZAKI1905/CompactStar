# BA13 final decision matrix

Updated inventory: **462 rows (21 stars × 22 groups)**. Full [JSON](evidence/matrix.json.gz) and [TSV](evidence/matrix.tsv.gz).
Historical values, thresholds, norms, terms, classes and provenance are retained; all missing full uncertainties stay null.
This successor adds common-blocker identifiers and a joint required-accuracy region. It does not rescore old outcomes.

| Group | Historical max utilization | Max finite addition / headroom | Minimum partial remaining | Required pointwise-error cap* |
|---|---:|---:|---:|---:|
| Cstar | 0.0145900852 | 0.0117708% | 9.8529393e-05 | 1.3719089e-05 |
| Ltilde_Me | 0.0122667955 | 3.07495% | 9.5736086e-05 | 1.3328918e-05 |
| Ltilde_Mmu | 0.25681093 | 0.00169956% | 7.4317644e-05 | 1.034788e-05 |
| N_e | 0.0488365517 | 9.17591e-05% | 9.5116258e-05 | 1.3243849e-05 |
| N_mu | 0.54075348 | 0.0105201% | 4.5919821e-05 | 6.3933745e-06 |
| N_n | 0.0097197128 | 4.04962e-06% | 9.9028025e-05 | 1.3788469e-05 |
| P0_average | 0.00188263921 | 2.40212e-06% | 9.9811734e-05 | 1.3897643e-05 |
| P1_average | 0.00261457001 | 3.48886e-06% | 9.973854e-05 | 1.3887451e-05 |
| P2_average | 0 | 0% | 0.0001 | — |
| Tsurface_inf | 0.00787060676 | 0.0653337% | 9.914812e-05 | 1.3805208e-05 |
| Z_row_npe | 0.0216301769 | 0.000102496% | 9.7836882e-05 | 1.3622667e-05 |
| Z_row_npmu | 0.350922112 | 0.00650567% | 6.4903566e-05 | 9.0370771e-06 |
| envelope | 0.00975368197 | 0.0654608% | 9.8959809e-05 | 1.3778988e-05 |
| metric_structure | 0.0330124185 | 1.31196% | 9.5430112e-05 | — |
| mu_B | 0.00188263921 | 2.40211e-06% | 9.9811734e-05 | 1.3897643e-05 |
| mu_n_profile | 0.00188265152 | 5.81965e-05% | 9.9811677e-05 | 1.3897635e-05 |
| radius | 0.0128306445 | 1.49536e-05% | 9.8716921e-05 | 1.3745203e-05 |
| species_support | 0 | 0% | 1 | — |
| surface_gravity | 0.0370511514 | 7.13343e-05% | 9.6294816e-05 | 1.3407952e-05 |
| t_e | 0.062762656 | 0.0164994% | 2.8910855e-06 | 4.0251918e-07 |
| t_mu | 0.376275898 | 0.164372% | 2.1032381e-07 | 2.9331769e-08 |
| t_n | 0.00332093366 | 0.000618164% | 9.625572e-05 | 1.3397208e-05 |

*Caps are single-axis feasibility limits, not estimates: additional gap/fit errors are set to zero for this display, and inherited finite terms are assumed valid. The full joint inequality is in every applicable JSON row and METHODS.md. Tangents are absolute dN_i/dB; other numeric caps use the historical scaled component norm. Metric structure has only P−R, not an independent reference discrepancy. Species support is categorical.

Disposition: **21 exact-P2 PASS, 0 FAIL, 441 UNQUALIFIED**. Approximate full utilization remains null. No finite sum is a complete qualification. Source: `evidence/summary.json:1`; formulas: `METHODS.md`, final section.
