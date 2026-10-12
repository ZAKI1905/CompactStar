import sys
def load(p):
    rows=[l.rstrip('\n').split('\t') for l in open(p) if not l.startswith('#') and l.strip()]
    return rows[0], [[float(x) for x in r] for r in rows[1:]]
def drift(a, b, label):
    h, A = load(a); h2, B = load(b); assert h == h2 and len(A) == len(B)
    worst = {}
    for ra, rb in zip(A, B):
        for name, x, y in zip(h, ra, rb):
            r = abs(x - y) / max(abs(x), abs(y), 1e-300) if x != y else 0.0
            worst[name] = max(worst.get(name, 0.0), r)
    top = sorted(worst.items(), key=lambda t: -t[1])[:6]
    print('%-34s identical_cols=%d/%d  worst: %s' % (label, sum(1 for v in worst.values() if v == 0), len(worst), ', '.join('%s=%.1e' % (k, v) for k, v in top)))
S = sys.argv[1]
macD = 'sources/CompactStar/tests/baselines/hartle_monopole_dscmf1_debug.tsv'
macR = 'sources/CompactStar/tests/dependency_migration/fixtures/Release-hartle_monopole_regression.tsv'
linD = S + '/hartle/Debug/hartle_monopole_dscmf1_debug.tsv'; linR = S + '/hartle/Release/hartle_monopole_dscmf1_debug.tsv'
drift(linD, macD, 'Linux Debug vs Mac Debug'); drift(linR, macR, 'Linux Release vs Mac Release')
drift(macD, macR, 'CONTROL Mac Debug vs Mac Release'); drift(linD, linR, 'CONTROL Linux Debug vs Linux Release')
