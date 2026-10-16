"""Common-B shooting through an approximate baryon-number functional.

The production star is shot to B_prod(star) = B_target, where B_prod is a nodal quadrature.
If B_prod = B_true (1 + e_B), the solved star has the wrong true B and central density:
d ln rho_c = - e_B / (d ln B / d ln rho_c).  We estimate e_B by a phase-free equilibrium-
coordinate quadrature B_nu and test whether it predicts (i) the radial ladder of rho_c and
(ii) the jitter-phase variation of rho_c.  Post-hoc diagnostic; not a certified error."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import PHASES, profile, star, B_nu, write

LADDER = ['r20000', 'r40000', 'nominal', 'r160000', 'r320000', 'r640000']


def main():
    rho0, rho20 = star('nominal', 0)['rho_c'], star('nominal', 20)['rho_c']
    sens = (-1e-6) / (rho20 / rho0 - 1)                    # d ln B / d ln rho_c along the sequence
    out = {'status': 'POST_HOC_DIAGNOSTIC', 'dlnB_dlnrho_c': sens, 'ladder': {}, 'jitter': {}}
    for i in (0, 20):
        rows = []
        for c in LADDER:
            st = star(c, i)
            e = (st['B_solved'] - B_nu(profile(c, i, (0, 1, 2, 3)))) / st['B_solved']
            rows.append(dict(case=c, e_B=e, rho_c=st['rho_c']))
        for j in range(1, len(rows)):
            rows[j]['rho_c_increment'] = rows[j]['rho_c'] / rows[j - 1]['rho_c'] - 1
            rows[j]['predicted_increment'] = (rows[j - 1]['e_B'] - rows[j]['e_B']) / sens
        out['ladder'][str(i)] = rows
        rows = []
        for c in PHASES:
            st = star(c, i)
            rows.append(dict(case=c, e_B=(st['B_solved'] - B_nu(profile(c, i, (0, 1, 2, 3)))) / st['B_solved'], rho_c=st['rho_c']))
        for x in rows:
            x['rho_c_vs_nominal'] = x['rho_c'] / rows[0]['rho_c'] - 1
            x['predicted_rho_c_vs_nominal'] = -(x['e_B'] - rows[0]['e_B']) / sens
        out['jitter'][str(i)] = rows
    write('b_functional.json', out)
    for i in ('0', '20'):
        print('ladder', i, ' '.join(f"{x['case']}:{x.get('rho_c_increment', 0):+.3e}/{x.get('predicted_increment', 0):+.3e}" for x in out['ladder'][i]))
        print('jitter', i, ' '.join(f"{x['case']}:{x['rho_c_vs_nominal']:+.2e}/{x['predicted_rho_c_vs_nominal']:+.2e}" for x in out['jitter'][i]))


if __name__ == '__main__':
    main()
