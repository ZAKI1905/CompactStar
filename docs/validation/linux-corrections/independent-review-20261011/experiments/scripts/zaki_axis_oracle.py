# Independent oracle: for every Zaki log-axis node in the historical fixed table, reconstruct the
# platform-independent input x_i = log10(min) + i*(log10(max)-log10(min))/res (double arithmetic),
# then compare Mac value (Apple libm pow), Linux value (glibc pow) and the correctly-rounded 10^x.
import ctypes, math, mpmath as mp, random
mp.mp.dps = 50
libm = ctypes.CDLL('/mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain/sysroot/usr/lib64/libm.so.6')
for f in ('pow',): getattr(libm, f).restype = ctypes.c_double; getattr(libm, f).argtypes = [ctypes.c_double, ctypes.c_double]
libm.log10.restype = ctypes.c_double; libm.log10.argtypes = [ctypes.c_double]
def cr(v):
    f = float(v); c = [f, math.nextafter(f, -math.inf), math.nextafter(f, math.inf)]
    return min(c, key=lambda t: abs(mp.mpf(t) - v))
def cr_log10(x): return cr(mp.log10(mp.mpf(x)))
def cr_pow10(x): return cr(mp.mpf(10) ** mp.mpf(x))
H = '/mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain/historical/'
def load(path, section):
    out = []
    with open(path) as fh:
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if p[0] == section: out.append(float.fromhex(p[2]))
    return out
# Rebuild the axis_historical inputs: 12 axes from tests/fp_preservation/characterize.cpp
cases = [((1e15,4e15),99,'Linear'),((3e13,5e15),99,'Linear'),((1e15,4e15),49,'Linear'),((3e13,5e15),49,'Linear'),
         ((8e14,5e15),99,'Log'),((5e13,5e16),150,'Log'),((1e15,4e15),255,'Linear'),((3e13,5e15),255,'Linear'),
         ((1e15,4e15),95,'Linear'),((3e13,1.5e15),95,'Linear'),((1e15,2.3e15),159,'Linear'),((3e13,1.4e15),159,'Linear')]
class Inputs:
    def __init__(s, st): s.state = st
    def next(s):
        s.state = (s.state*6364136223846793005 + 1442695040888963407) & ((1<<64)-1)
        return float(s.state >> 11) * 2.0**-53
rng = Inputs(0x5eeda11ce)
rand_axes = []
for k in range(300):
    lo = -3 + 19*rng.next()
    mn = libm.pow(10.0, lo); mx = libm.pow(10.0, lo + 0.5 + 3.5*rng.next())
    n = 10 + int(290*rng.next()); rand_axes.append(((mn, mx), n, 'Log'))
def analyse(section, axes):
    mac, lin = load(H+'reference-fixed.tsv', section), load(H+'Zaki-Debug-fixed.tsv', section)
    assert len(mac) == len(lin), (len(mac), len(lin))
    k = 0; stats = dict(n=0, log=0, diff=0, mac_cr=0, lin_cr=0, both_cr=0, neither=0, input_dep=0, lin_match=0)
    for (mn, mx), res, scale in axes:
        for i in range(res+1):
            m, l = mac[k], lin[k]; k += 1; stats['n'] += 1
            if scale != 'Log': continue
            stats['log'] += 1
            # input with glibc log10 and with correctly-rounded log10
            a, b = libm.log10(mn), libm.log10(mx); x_g = a + i*(b - a)/res
            a2, b2 = cr_log10(mn), cr_log10(mx); x_c = a2 + i*(b2 - a2)/res
            if x_g != x_c: stats['input_dep'] += 1
            if libm.pow(10.0, x_g) == l: stats['lin_match'] += 1
            if m == l: continue
            stats['diff'] += 1
            truth = cr_pow10(x_c)
            mc, lc = (m == truth), (l == truth)
            stats['mac_cr'] += mc; stats['lin_cr'] += lc; stats['both_cr'] += (mc and lc); stats['neither'] += (not mc and not lc)
    assert k == len(mac), (k, len(mac))
    print(section, stats)
analyse('axis_historical', cases)
analyse('axis_random', rand_axes)
# Survey: glibc pow(10,x) vs correctly rounded over random x in the stellar log range [log10(8e14), log10(5e16)]
random.seed(1); bad = 0; N = 20000
for _ in range(N):
    x = random.uniform(math.log10(8e14), math.log10(5e16))
    if libm.pow(10.0, x) != cr_pow10(x): bad += 1
print("glibc pow(10,x) survey in stellar range: %d/%d not correctly rounded" % (bad, N))
