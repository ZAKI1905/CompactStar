"""P5 audit: separate EOS-table refinement from shooting, structure and onset effects.

(1) EOS-ladder increments (8192 -> 16384c -> 32768c -> 65536c) of structure (R, M, N_mu,
    rho_c, nu_surface) and of the estimator quantities (L, C, r_e, mu_inf) per target.
(2) Common-B residual and bracket of every new star (shooting error).
(3) Node-level ln P difference between 32768c and 65536c at the identical radial grid around
    the neutron onset: where the perturbation enters and the outward offset it leaves."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
from common import PILOT, TON, profile, star, analysis, write

LADDER = ['e2048', 'e4096', 'nominal', 'e16384c', 'e32768c', 'e65536c']


def main():
    out = {'status': 'POST_HOC_DIAGNOSTIC', 'increments': {}, 'shooting': {}, 'onset_localization': {}}
    for i in (0, 10, 20):
        inc = {}
        for key in ('R', 'M', 'N_mu', 'rho_c', 'nu_surface'):
            v = [star(c, i)[key] for c in LADDER]
            inc[key] = [(v[j + 1] - v[j]) / abs(v[-1]) for j in range(len(v) - 1)]
        for key in ('L_nu', 'C_nu', 'r_nu', 'mu_inf'):
            v = [analysis(c, i)[key] for c in LADDER]
            inc[key] = [(v[j + 1] - v[j]) / abs(v[-1]) for j in range(len(v) - 1)]
        out['increments'][str(i)] = inc
        a = profile('e32768c', i, (0, 3, 4))
        b = profile('e65536c', i, (0, 3, 4))
        n = min(len(a), len(b)) - 1
        assert np.array_equal(a[:n, 0], b[:n, 0])
        r = a[:n, 0]
        dlp = np.log(b[:n, 2] / a[:n, 2])
        k = int(np.searchsorted(-a[:, 1], -TON))
        inner = dlp[(r > 12.0) & (r < a[k - 1, 0] - 0.05)]
        outer = dlp[(r > a[k, 0]) & (r < a[k, 0] + 0.02)]
        out['onset_localization'][str(i)] = dict(
            first_pe_node=k, node_offsets=list(range(-6, 6)), dlnP_near_onset=[float(dlp[j]) for j in range(k - 6, k + 6)],
            mean_dlnP_interior=float(inner.mean()), mean_dlnP_outside_onset=float(outer.mean()),
            outward_offset_jump=float(outer.mean() - inner.mean()),
            R_change_rel=float(b[-1, 0] / a[-1, 0] - 1) if len(a) == len(b) else None)
    for p in sorted((PILOT / 'runs').glob('*.receipt.json')):
        name = p.name.replace('.receipt.json', '')
        case, i = name.rsplit('-', 1)
        st = star(case, int(i))
        out['shooting'][name] = dict(residual_rel=st['residual'] / st['B_target'], bracket_rel=st['bracket'] / st['B_target'],
                                     evaluations=st['evaluations'])
    out['max_abs_shooting_residual_rel'] = max(abs(v['residual_rel']) for v in out['shooting'].values())
    write('eos_decomposition.json', out)
    for i in ('0', '10', '20'):
        print('target', i, {k: ' '.join(f'{x:+.2e}' for x in v) for k, v in out['increments'][i].items()})
        o = out['onset_localization'][i]
        print('   onset localization: interior', f"{o['mean_dlnP_interior']:+.2e}", 'outside', f"{o['mean_dlnP_outside_onset']:+.2e}", 'R change', o['R_change_rel'])
    print('max |shooting residual|', out['max_abs_shooting_residual_rel'])


if __name__ == '__main__':
    main()
