# Simulate the arm64 Mac path: long double == double, AppleClang default -ffp-contract=on.
import mpmath as mp, math
mp.mp.dps = 60
MN, MP, ME, K = 939.56542052, 938.27208816, 0.51099895000, 5.067730716156395e-3
def rd(x):  # round mpf to nearest double
    return float(x)
def fma(a, b, c):  # single rounding of a*b+c
    return rd(mp.mpf(a) * mp.mpf(b) + mp.mpf(c))
reviewed = mp.mpf('7.3567289037328326656352e-9')
mac_logged_indep = 7.356728903732833408335e-09   # PHASE5A5_TRACKR_PE_INDEPENDENT_REVIEW.md:643
def independent(contract):
    hbar_c = 1.0 / K
    d, s = MN - MP, MN + MP
    mu_e = (fma(d, s, ME * ME) if contract else (d * s + ME * ME)) / (2.0 * MN)
    mom = math.sqrt((mu_e - ME) * (mu_e + ME))
    pi = math.acos(-1.0)
    return mom * mom * mom / (3.0 * pi * pi * hbar_c * hbar_c * hbar_c)
def production(contract):
    d, s = MN - MP, MN + MP
    mu_e = (fma(d, s, ME * ME) if contract else (d * s + ME * ME)) / (2.0 * MN)
    hbar_c = 1.0 / K
    p_F = math.sqrt((mu_e - ME) * (mu_e + ME)); pi = math.acos(-1.0)
    ratio = p_F / hbar_c
    return ratio * ratio * ratio / (3.0 * pi * pi)
for contract in (False, True):
    i, p = independent(contract), production(contract)
    print("contract=%s independent=%.22e production=%.22e" % (contract, i, p))
    print("   matches Mac logged independent? %s" % (i == mac_logged_indep))
    print("   |indep-rev|/rev=%.3e (<2e-16?)  |prod/indep-1|=%.3e (<3e-16?)" % (abs((mp.mpf(i)-reviewed)/reviewed), abs(p/i-1)))
