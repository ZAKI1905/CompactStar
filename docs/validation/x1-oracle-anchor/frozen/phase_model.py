"""P7 audit: a-priori inverse-interpolation curvature model for the linear-in-nu support edge,
coefficient-free Hermite edge, and attribution of the residual to common-B rho_c phase noise.

Error model (derived, not fitted):  r_lin - r_true = -(nu''/2nu') theta(1-theta) h^2 + O(h^3),
nu''/(2nu') = -(r-M)/(r(r-2M)) when P, eps are negligible at the onset (exterior form).
Post-hoc diagnostic on retained pilot data; NOT a validation of any estimator."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
from common import PHASES, profile, star, analysis, edges, write

DLNL_DLNB, DLNR_DLNB = -1.809209e-6 / 1e-6, -1.257426e-6 / 1e-6   # independent oracle paired drifts at index 20


def main():
    out = {'status': 'POST_HOC_DIAGNOSTIC_NOT_VALIDATION', 'cases': {}, 'summary': {}}
    for i in (0, 10, 20):
        for c in PHASES + ['r640000']:
            d = profile(c, i)
            e = edges(d)
            a = analysis(c, i)
            s = a['edge_sensitivity']                       # dlnL/dr_e (1/km)
            e.update(L_lin=a['L_nu'], L_herm=a['L_nu'] * (1 - s * (e['r_lin'] - e['r_herm'])),
                     L_pred_corrected=a['L_nu'] * (1 - s * e['predicted_lin_minus_true']),
                     edge_sensitivity=s, rho_c=star(c, i)['rho_c'])
            out['cases'][f'{i}:{c}'] = e
    # a-priori maxima of the curvature term at nominal h (theta(1-theta) <= 1/4)
    for i in (0, 10, 20):
        e0 = out['cases'][f'{i}:nominal']
        rows = [out['cases'][f'{i}:{c}'] for c in PHASES]
        rng = lambda key: (max(x[key] for x in rows) - min(x[key] for x in rows)) / abs(rows[0][key])
        out['summary'][str(i)] = dict(
            coef_nu2_over_2nu1=-e0['coef'], a_priori_max_edge_range_rel=e0['coef'] * e0['h'] ** 2 / 4 / e0['r_herm'],
            a_priori_max_L_range_rel=e0['edge_sensitivity'] * e0['coef'] * e0['h'] ** 2 / 4,
            observed_edge_range_rel=rng('r_lin'), observed_L_range_rel=rng('L_lin'),
            hermite_edge_range_rel=rng('r_herm'), hermite_L_range_rel=rng('L_herm'),
            max_abs_lin_minus_herm_vs_prediction=max(abs((x['r_lin'] - x['r_herm']) / x['predicted_lin_minus_true'] - 1)
                                                     for x in rows if x['predicted_lin_minus_true'] > 1e-12))
    # rho_c attribution of the Hermite residual (sequence sensitivity, not fitted)
    rho0, rho20 = star('nominal', 0)['rho_c'], star('nominal', 20)['rho_c']
    dlnB_dlnrho = (-1e-6) / (rho20 / rho0 - 1)
    sL = DLNL_DLNB * dlnB_dlnrho
    out['sensitivities'] = dict(dlnB_dlnrho_c=dlnB_dlnrho, dlnL_dlnrho_c=sL, dlnr_e_dlnrho_c=DLNR_DLNB * dlnB_dlnrho)
    for i in (0, 10, 20):
        rows = [out['cases'][f'{i}:{c}'] for c in PHASES]
        res = [(x['L_herm'] / rows[0]['L_herm'] - 1) - sL * (x['rho_c'] / rows[0]['rho_c'] - 1) for x in rows]
        out['summary'][str(i)].update(hermite_L_residual_after_rho_c=res, rho_c_jitter=[x['rho_c'] / rows[0]['rho_c'] - 1 for x in rows])
    # paired: observed vs a-priori curvature-difference prediction (correct sign: L_lin = L_true (1 + s err))
    for i in (10, 20):
        obs, pred, herm = [], [], []
        for c in PHASES:
            a0, ai = out['cases'][f'0:{c}'], out['cases'][f'{i}:{c}']
            obs.append(ai['L_lin'] / a0['L_lin'] - 1)
            pred.append(a0['edge_sensitivity'] * (ai['predicted_lin_minus_true'] - a0['predicted_lin_minus_true']))
            herm.append(ai['L_herm'] / a0['L_herm'] - 1)
        out['summary'][f'paired_{i}'] = dict(
            observed_paired_L_range=max(obs) - min(obs), predicted_curvature_paired_range=max(pred) - min(pred),
            observed_minus_prediction_range=max(o - p for o, p in zip(obs, pred)) - min(o - p for o, p in zip(obs, pred)),
            hermite_paired_L_values=herm, hermite_paired_L_range=max(herm) - min(herm))
    write('phase_model.json', out)
    for k, v in out['summary'].items():
        print(k, {kk: (f'{vv:.4e}' if isinstance(vv, float) else ([f'{x:.2e}' for x in vv] if isinstance(vv, list) else vv)) for kk, vv in v.items()})
    print('sensitivities', out['sensitivities'])


if __name__ == '__main__':
    main()
