import mpmath as mp, struct
mp.mp.dps = 60
MN, MP, ME, K = 939.56542052, 938.27208816, 0.51099895000, 5.067730716156395e-3  # exact binary doubles
mn, mpv, me, k = (mp.mpf(x) for x in (MN, MP, ME, K))
mu_e = ((mn - mpv) * (mn + mpv) + me * me) / (2 * mn)
p = mp.sqrt((mu_e - me) * (mu_e + me))
hbarc = 1 / k
truth = p**3 / (3 * mp.pi**2 * hbarc**3)
print("truth (exact constants, 60 digits):", mp.nstr(truth, 30))
reviewed = mp.mpf('7.3567289037328326656352e-9')
print("reviewed literal rel dev from truth: %.3e" % float((reviewed - truth) / truth))
def ulp(x):
    import math
    return math.ulp(x)
for name, v in [("Linux production", 7.356728903732837544237787e-09), ("Linux independent(80-bit)", mp.mpf('7.356728903732831692581188e-09'))]:
    d = (mp.mpf(v) - truth) / truth
    print("%-28s rel dev from truth = %+.3e  (= %+.2f double ULPs)" % (name, float(d), float((mp.mpf(v) - truth) / ulp(float(v)))))
# correctly rounded double of truth
t = float(truth)
import math
cands = [t, math.nextafter(t, 0), math.nextafter(t, 1)]
best = min(cands, key=lambda c: abs(mp.mpf(c) - truth))
print("correctly rounded double of truth: %.25e  (%a)" % (best, best))
