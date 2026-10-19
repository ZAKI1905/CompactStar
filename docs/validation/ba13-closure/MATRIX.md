# Complete frozen-threshold BA13 matrix

All 462 rows are in [matrix.json.gz](evidence/matrix.json.gz) and [matrix.tsv.gz](evidence/matrix.tsv.gz).
The JSON carries terms, evidence classes, hypotheses, units, scales, source pointers, correlations and refusal reasons.
Component values remain in the external component-data.json.gz, authenticated by the manifest.
Finite partial utilization is a diagnostic; it is never a qualified total.

| Quantity | PASS | FAIL | UNQUALIFIED | Minimum headroom | Largest partial utilization |
|---|---:|---:|---:|---:|---:|
| t_n | 0 | 0 | 21 | 9.62563146e-05 | 0.0282863307 |
| t_e | 0 | 0 | 21 | 2.89156255e-06 | 0.0850172911 |
| t_mu | 0 | 0 | 21 | 2.10670097e-07 | 0.394434096 |
| Z_row_npe | 0 | 0 | 21 | 9.78369823e-05 | 0.0216550257 |
| Z_row_npmu | 0 | 0 | 21 | 6.49077888e-05 | 0.351966968 |
| Cstar | 0 | 0 | 21 | 9.85409916e-05 | 0.0147498278 |
| Ltilde_Me | 0 | 0 | 21 | 9.87733204e-05 | 0.0427041026 |
| Ltilde_Mmu | 0 | 0 | 21 | 7.4318907e-05 | 0.257659528 |
| mu_B | 0 | 0 | 21 | 9.98117361e-05 | 0.00188313618 |
| mu_n_profile | 0 | 0 | 21 | 9.98117348e-05 | 0.00188320377 |
| P0_average | 0 | 0 | 21 | 9.98117361e-05 | 0.00188313618 |
| P1_average | 0 | 0 | 21 | 9.9738543e-05 | 0.00261692036 |
| P2_average | 21 | 0 | 0 | 0.0001 | 0 |
| N_n | 0 | 0 | 21 | 9.90280287e-05 | 0.00972058324 |
| N_e | 0 | 0 | 21 | 9.51163448e-05 | 0.0489388629 |
| N_mu | 0 | 0 | 21 | 4.5924652e-05 | 0.542381252 |
| species_support | 0 | 0 | 21 | 1 | 0 |
| metric_structure | 0 | 0 | 21 | 9.66987581e-05 | 0.0456988819 |
| radius | 0 | 0 | 21 | 9.87169356e-05 | 0.0128341342 |
| surface_gravity | 0 | 0 | 21 | 9.6294885e-05 | 0.0370615575 |
| envelope | 0 | 0 | 21 | 9.90246318e-05 | 0.0104046912 |
| Tsurface_inf | 0 | 0 | 21 | 9.92129394e-05 | 0.0085210798 |

Total: **21 PASS / 0 FAIL / 441 UNQUALIFIED**. No approximate group has complete evidence.
Muon-sector numerical rows: **84 UNQUALIFIED**; exact support/domain rows: **21 UNQUALIFIED**.
Disposition **E**; [scientific report](REPORT.md) and [owner decisions](OWNER_DECISIONS.md).
