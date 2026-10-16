"""Absolute-accuracy closure check and frozen-candidate predictions for an oracle-anchor test.

Production-side analytic-EOS limit, from existing data only:
  X_lim = X_640 + (X_640 - X_320)/3          Richardson with p = 2 fixed a priori
          (L and r_e use the coefficient-free Hermite edge, removing the curvature term)
  X_inf = X_lim * (1 + EOS offset)           offset = X(65536c)/X(8192) - 1 at r80000
Compared post hoc with the EXISTING independent oracle (tolerance 2e-12); the comparison is
a consistency check, not blind.  The X_inf values are the predictions a future, converged
oracle must meet (see the report, section 9)."""
import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import PILOT, profile, star, analysis, edges, write


def hermite_L(case, i):
    e = edges(profile(case, i))
    a = analysis(case, i)
    return e['r_herm'], a['L_nu'] * (1 - a['edge_sensitivity'] * (e['r_lin'] - e['r_herm'])), a['C_nu']


def main():
    metrics = json.loads((PILOT / 'metrics.json').read_text())
    out = {'status': 'POST_HOC_CONSISTENCY_AND_FROZEN_CANDIDATE_PREDICTIONS', 'targets': {}}
    for i in (0, 10, 20):
        r3, L3, C3 = hermite_L('r320000', i)
        r6, L6, C6 = hermite_L('r640000', i)
        s3, s6 = star('r320000', i), star('r640000', i)
        x320 = dict(L=L3, C=C3, r_e=r3, R=s3['R'], M=s3['M'])
        x640 = dict(L=L6, C=C6, r_e=r6, R=s6['R'], M=s6['M'])
        eos = {k: analysis('e65536c', i)[f] / analysis('nominal', i)[f] - 1 for k, f in (('L', 'L_nu'), ('C', 'C_nu'), ('r_e', 'r_nu'))}
        eos.update({k: star('e65536c', i)[k] / star('nominal', i)[k] - 1 for k in ('R', 'M')})
        last_eos = {k: analysis('e65536c', i)[f] / analysis('e32768c', i)[f] - 1 for k, f in (('L', 'L_nu'), ('C', 'C_nu'), ('r_e', 'r_nu'))}
        last_eos.update({k: star('e65536c', i)[k] / star('e32768c', i)[k] - 1 for k in ('R', 'M')})
        o = metrics['oracle'][f'{i}:2e-12']
        oc = metrics['oracle'][f'{i}:2e-10']
        oracle = dict(L=o['L_nu'], C=o['C_nu'], r_e=o['r_nu'], R=o['R'], M=o['M'])
        oracle_coarse = dict(L=oc['L_nu'], C=oc['C_nu'], r_e=oc['r_nu'], R=oc['R'], M=oc['M'])
        t = {}
        for k in x640:
            lim = x640[k] + (x640[k] - x320[k]) / 3
            pred = lim * (1 + eos[k])
            t[k] = dict(x640=x640[k], richardson_tail_rel=lim / x640[k] - 1, eos_offset_rel=eos[k],
                        last_eos_increment_rel=last_eos[k], predicted_analytic_limit=pred,
                        existing_oracle_2e12=oracle[k], oracle_minus_prediction_rel=oracle[k] / pred - 1,
                        oracle_tolerance_change_rel=oracle[k] / oracle_coarse[k] - 1)
        out['targets'][str(i)] = t
    write('anchor.json', out)
    for i, t in out['targets'].items():
        print('target', i, '  '.join(f"{k}: oracle-pred {v['oracle_minus_prediction_rel']:+.2e} (tol-change {v['oracle_tolerance_change_rel']:+.2e})" for k, v in t.items()))


if __name__ == '__main__' and len(sys.argv) == 1:
    main()


# d ln X / d ln B along the equilibrium sequence (independent oracle paired drifts over dB/B = -1e-6)
DLN_X_DLN_B = dict(L=-1.809209, C=-1.470015, r_e=-1.257426, R=-1.283062, M=0.950700)
HP_ONSET_KM = 0.045      # pressure scale height at the neutron onset, P/|dP/dr| (km)
KINK_DLNP = 10 * 1e-10   # allowance: 10 x the production rk8pd relative tolerance for the onset crossing


def bands():
    """Model-derived acceptance bands for the proposed oracle-anchor experiment (relative units).
    Components are a-priori or independently measured; none uses the oracle-vs-prediction differences.
      q-term : |q_obs - 1/4| / (1 - q_obs) * |X640 - X320|, q_obs from the 160k/320k/640k ladder
      eos    : |X(65536c)/X(32768c) - 1|, the last measured EOS increment (retained at face value)
      rho    : (1/8) * (nominal jitter rho_c range) * |dlnX/dlnrho_c|   (O(h) phase term at 640k)
      kink   : 10 * rtol(1e-10) * H_P(onset) mapped to X (outer-layer shift allowance)"""
    out = json.loads((Path(__file__).resolve().parent / 'anchor.json').read_text())
    pm = json.loads((Path(__file__).resolve().parent / 'phase_model.json').read_text())
    sens = pm['sensitivities']
    for i in (0, 10, 20):
        r1, L1, C1 = hermite_L('r160000', i)
        s1 = star('r160000', i)
        x160 = dict(L=L1, C=C1, r_e=r1, R=s1['R'], M=s1['M'])
        t = out['targets'][str(i)]
        rho_range = max(pm['summary'][str(i)]['rho_c_jitter']) - min(pm['summary'][str(i)]['rho_c_jitter'])
        s_edge = analysis('r640000', i)['edge_sensitivity']          # d ln L / d r_e  (1/km)
        for k, v in t.items():
            x640 = v['x640']
            x320 = x640 - 3 * v['richardson_tail_rel'] * x640         # inverts lim = x640 + (x640 - x320)/3
            q = (x640 - x320) / (x320 - x160[k])
            qterm = abs(q - 0.25) / (1 - q) * abs(x640 - x320) / abs(x640)
            rho = rho_range / 8 * abs(DLN_X_DLN_B[k]) * sens['dlnB_dlnrho_c']
            shift = KINK_DLNP * HP_ONSET_KM                           # outer-layer radial shift allowance, km
            kink = {'L': s_edge * shift, 'r_e': shift / x640, 'R': shift / x640}.get(k, 0.0)
            band = qterm + abs(v['last_eos_increment_rel']) + rho + kink
            v.update(q_obs_640=q, band_qterm=qterm, band_eos=abs(v['last_eos_increment_rel']), band_rho=rho,
                     band_kink=kink, band_total_excluding_oracle=band,
                     existing_oracle_within_band=abs(v['oracle_minus_prediction_rel']) <= band)
    write('anchor.json', out)
    for i, t in out['targets'].items():
        print('target', i, '  '.join(f"{k}: band {v['band_total_excluding_oracle']:.2e} (q {v['band_qterm']:.1e}, eos {v['band_eos']:.1e}, rho {v['band_rho']:.1e}, kink {v['band_kink']:.1e}) existing |diff| {abs(v['oracle_minus_prediction_rel']):.2e}" for k, v in t.items()))


if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'bands':
    bands()
