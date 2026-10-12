import ctypes, math, mpmath as mp
mp.mp.dps = 50
libm = ctypes.CDLL('/mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain/sysroot/usr/lib64/libm.so.6')
libm.pow.restype = ctypes.c_double; libm.pow.argtypes = [ctypes.c_double, ctypes.c_double]
libm.log10.restype = ctypes.c_double; libm.log10.argtypes = [ctypes.c_double]
def cr(v):
    f = float(v); c = [f, math.nextafter(f, -math.inf), math.nextafter(f, math.inf)]
    return min(c, key=lambda t: abs(mp.mpf(t) - v))
def cr_log10(x): return cr(mp.log10(mp.mpf(x)))
def cr_pow10(x): return cr(mp.mpf(10) ** mp.mpf(x))
H = '/mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain/historical/'
def load(path, section):
    return [float.fromhex(l.split('\t')[2]) for l in open(path) if l.split('\t')[0] == section]
class Inputs:
    def __init__(s, st): s.state = st
    def next(s):
        s.state = (s.state*6364136223846793005 + 1442695040888963407) & ((1<<64)-1)
        return float(s.state >> 11) * 2.0**-53
rng = Inputs(0x5eeda11ce); axes = []
for k in range(300):
    lo = -3 + 19*rng.next(); mn = libm.pow(10.0, lo); mx = libm.pow(10.0, lo + 0.5 + 3.5*rng.next())
    n = 10 + int(290*rng.next()); axes.append(((mn, mx), n))
mac, lin = load(H+'reference-fixed.tsv', 'axis_random'), load(H+'Zaki-Debug-fixed.tsv', 'axis_random')
k = 0; cat = {}
def bump(key): cat[key] = cat.get(key, 0) + 1
# NOTE: the random axis endpoints mn,mx themselves come from pow(10,lo) in the harness, so on Mac they may differ too.
for (mn, mx), res in axes:
    for i in range(res+1):
        m, l = mac[k], lin[k]; k += 1
        if m == l: continue
        a, b = libm.log10(mn), libm.log10(mx); x_g = a + i*(b - a)/res
        a2, b2 = cr_log10(mn), cr_log10(mx); x_c = a2 + i*(b2 - a2)/res
        lin_ok = (l == cr_pow10(x_g))              # glibc pow correctly rounded at glibc's own input?
        if x_g == x_c:
            mac_ok = (m == cr_pow10(x_c))
            bump('same input: glibc pow CR=%s, Mac pow CR=%s' % (lin_ok, mac_ok))
        else:
            mac_ok_c = (m == cr_pow10(x_c))        # Mac consistent with a CR log10 input and CR pow
            bump('glibc log10 input not CR: glibc pow CR at own input=%s; Mac == CR(10^CR-input)=%s' % (lin_ok, mac_ok_c))
for key, n in sorted(cat.items(), key=lambda t: -t[1]): print('%4d  %s' % (n, key))
print('total differing:', sum(cat.values()))
