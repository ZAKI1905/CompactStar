"""Radial ladders 20k..640k of the support edge and Ltilde_Me with the linear-in-nu edge (pilot)
and the coefficient-free Hermite edge.  Post-hoc diagnostic on discovery + pilot profiles."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import profile, analysis, edges, write

LADDER = ['r20000', 'r40000', 'nominal', 'r160000', 'r320000', 'r640000']


def main():
    out = {'status': 'POST_HOC_DIAGNOSTIC'}
    for i in (0, 10, 20):
        rows = []
        for c in LADDER:
            e = edges(profile(c, i))
            a = analysis(c, i)
            rows.append(dict(case=c, theta=e['theta'], r_lin=e['r_lin'], r_herm=e['r_herm'], L_lin=a['L_nu'],
                             L_herm=a['L_nu'] * (1 - a['edge_sensitivity'] * (e['r_lin'] - e['r_herm']))))
        res = {'rows': rows}
        for key in ('r_lin', 'r_herm', 'L_lin', 'L_herm'):
            v = [x[key] for x in rows]
            inc = [(v[j + 1] - v[j]) / v[-1] for j in range(len(v) - 1)]
            res[key] = dict(increments=inc, ratios=[inc[j + 1] / inc[j] for j in range(len(inc) - 1)])
        out[str(i)] = res
        print(f'target {i}:', '  '.join(f"{k} ratios " + ' '.join(f'{x:.4f}' for x in res[k]['ratios']) for k in ('r_lin', 'r_herm', 'L_lin', 'L_herm')))
    write('hermite_ladder.json', out)


if __name__ == '__main__':
    main()
