import ctypes, math, mpmath as mp
mp.mp.dps = 60
libm = ctypes.CDLL('/mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain/sysroot/usr/lib64/libm.so.6')
libm.pow.restype = ctypes.c_double; libm.pow.argtypes = [ctypes.c_double, ctypes.c_double]
libm.exp10.restype = ctypes.c_double; libm.exp10.argtypes = [ctypes.c_double]
def cr10(x):
    v = mp.mpf(10) ** mp.mpf(x)
    f = float(v); c = [f, math.nextafter(f, -math.inf), math.nextafter(f, math.inf)]
    best = min(c, key=lambda t: abs(mp.mpf(t) - v))
    # distance of true value from the midpoint between best and its neighbour toward v, in ULP
    return best, v
L = float.fromhex('0x1.22fad6cb53501p-10')
MAC140, MAC13F = float.fromhex('0x1.00a7b7380b14p+0'), float.fromhex('0x1.00a7b7380b13fp+0')
mid = (mp.mpf(MAC140) + mp.mpf(MAC13F)) / 2
ulp = mp.mpf(MAC140) - mp.mpf(MAC13F)
inputs = {'L (exact level)': L, "L' = 1-(1-L) [center->corner edge]": 1.0 - (1.0 - L), "2-(2-L)": 2.0 - (2.0 - L), "(L/2)*2": (L / 2) * 2}
print("Mac table: RAW0/2 x=%a (pow at L); RAW4/5 x=%a; Linux RAW4/5 x=%a" % (MAC140, MAC140, MAC13F))
for name, x in inputs.items():
    g = libm.pow(10.0, x); e = libm.exp10(x); best, v = cr10(x)
    print("%-38s x=%a\n   glibc pow=%a  glibc exp10=%a  correctly-rounded=%a  (10^x - midpoint)/ULP=%+.4f" % (name, x, g, e, best, float((v - mid) / ulp)))
