"""Shared read-only helpers for the independent onset-pilot review calculations.

All inputs are read from the authenticated evidence trees; nothing is written there.
Outputs (JSON) are written only to the directory containing this file.
Run every script with:  <toolchain python> -I -B <script>   from this directory.
"""
import csv
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

ROOT = Path('/mnt/sdd/zaki/CompactStar/linux-transition-20261010')
PILOT = ROOT / 'ba13-onset-pilot'           # onset pilot evidence (read only)
STUDY = ROOT / 'ba13-uncertainty' / 'runs'   # previous study runs (read only)
HERE = Path(__file__).resolve().parent

# Exported binary constants of the Track-R free-gas model (model.txt of every profile).
MN, MP, ME, MM = 939.56542051999998, 938.27208815999995, 0.51099894999999995, 105.65837550000001
HC = 197.32698045930249
TON = 7.3567289037328318e-9                     # production neutron-onset density, fm^-3
MUM = float(np.hypot(MP, np.sqrt(MM * MM - ME * ME)) + MM)   # mu_B at muon onset
PHASES = ['nominal'] + [f'r8000{k}' for k in range(1, 6)]
OLD_CASES = {'nominal', 'r20000', 'r40000', 'r160000', 'r320000', 'e2048', 'e4096', 'e16384'}


def run_dir(case, i):
    return (STUDY if case in OLD_CASES else PILOT / 'runs') / f'{case}-{i}'


def profile(case, i, cols=(0, 1, 2, 3, 4)):
    """r[km], m[km], nu, nB[fm^-3], P_geom[km^-2] (columns of the retained profile)."""
    return np.loadtxt(run_dir(case, i) / 'profile' / 'profile.tsv', skiprows=1, usecols=cols)


def star(case, i):
    row = next(csv.DictReader((run_dir(case, i) / 'star.tsv').open(), delimiter='\t'))
    return {k: float(v) for k, v in row.items()}


def analysis(case, i):
    """Pilot estimator output for a case (discovery.json holds the pre-existing profiles)."""
    if case in OLD_CASES:
        return json.loads((PILOT / 'discovery.json').read_text())[f'{i}:{case}']
    return json.loads((PILOT / 'analysis' / f'{case}-{i}.json').read_text())


def mu_pe(n):
    p = HC * np.cbrt(3 * np.pi ** 2 * n)
    return np.hypot(MP, p) + np.hypot(ME, p)


def baryon_density(mu):
    """n_B of free npe(mu) matter at baryon chemical potential mu (MeV), fm^-3."""
    mu = np.asarray(mu, float)
    out = np.zeros_like(mu)
    a = mu < MUM
    v = mu[a]
    p2 = ((v - MP - ME) * (v + MP + ME) * (v - MP + ME) * (v + MP - ME)) / (4 * v * v)
    pp = np.sqrt(np.maximum(p2, 0))
    pn = np.sqrt(np.maximum(v * v - MN * MN, 0))
    out[a] = (pn ** 3 + pp ** 3) / (3 * np.pi ** 2 * HC ** 3)
    if np.any(~a):
        v = mu[~a]
        lo = np.full(v.shape, MM)
        hi = v - MP
        for _ in range(60):
            mid = (lo + hi) / 2
            pe = np.sqrt(np.maximum((mid - ME) * (mid + ME), 0))
            pm = np.sqrt(np.maximum((mid - MM) * (mid + MM), 0))
            pos = np.hypot(MP, np.cbrt(pe ** 3 + pm ** 3)) + mid > v
            hi = np.where(pos, mid, hi)
            lo = np.where(pos, lo, mid)
        lam = (lo + hi) / 2
        pe = np.sqrt((lam - ME) * (lam + ME))
        pm = np.sqrt(np.maximum((lam - MM) * (lam + MM), 0))
        pn = np.sqrt(v * v - MN * MN)
        out[~a] = (pn ** 3 + pe ** 3 + pm ** 3) / (3 * np.pi ** 2 * HC ** 3)
    return out


def anchor_mu_inf(d):
    r, m, nu, nb = d[:, :4].T
    k = int(np.searchsorted(-nb, -TON))      # first pe node outside the neutron onset
    return k, float(np.exp(nu[k]) * mu_pe(nb[k]))


def edges(d):
    """Linear-in-nu edge (pilot estimator), coefficient-free Hermite edge using the exact
    TOV derivative nu'=(m+4 pi r^3 P)/(r(r-2m)) at the two bracketing nodes, and the a-priori
    leading curvature prediction for the linear edge error."""
    r, m, nu, nb, P = d.T
    k, mui = anchor_mu_inf(d)
    a = k - 1
    non = np.log(mui / MN)
    rlin = float(np.interp(non, nu, r))
    dn = (m + 4 * np.pi * r ** 3 * P) / (r * (r - 2 * m))
    h = r[k] - r[a]

    def herm(s):
        return ((2 * s ** 3 - 3 * s ** 2 + 1) * nu[a] + (s ** 3 - 2 * s ** 2 + s) * h * dn[a]
                + (-2 * s ** 3 + 3 * s ** 2) * nu[k] + (s ** 3 - s ** 2) * h * dn[k])

    def dherm(s):
        return ((6 * s ** 2 - 6 * s) * nu[a] + (3 * s ** 2 - 4 * s + 1) * h * dn[a]
                + (-6 * s ** 2 + 6 * s) * nu[k] + (3 * s ** 2 - 2 * s) * h * dn[k])

    s = (rlin - r[a]) / h
    for _ in range(60):
        s -= (herm(s) - non) / dherm(s)
    rherm = float(r[a] + s * h)
    theta = (rherm - r[a]) / h
    M = float(np.interp(rherm, r, m))
    coef = (rherm - M) / (rherm * (rherm - 2 * M))          # = -nu''/(2 nu') with P, eps -> 0
    predicted = coef * theta * (1 - theta) * h * h           # r_linear - r_true, km (always >= 0)
    return dict(k=k, h=float(h), theta=float(theta), M_edge=M, coef=float(coef), mu_inf=mui,
                r_lin=rlin, r_herm=rherm, predicted_lin_minus_true=float(predicted))


def B_nu(d, order=24):
    """Phase-free baryon number: composition from mu_B = mu_inf exp(-nu), onset-split
    Gauss-Legendre with t^2 substitution next to each onset; linear nu, m between nodes."""
    r, m, nu, nb = d[:, :4].T
    k, mui = anchor_mu_inf(d)
    ons = [float(np.interp(np.log(mui / v), nu, r)) for v in (MN, MUM)]
    cuts = np.unique(np.r_[0., r, ons])
    x, w = leggauss(order)
    t = (x + 1) / 2
    w = w / 2
    tot = np.longdouble(0)
    for s in range(0, len(cuts) - 1, 2048):
        A = cuts[s:min(s + 2048, len(cuts) - 1), None]
        B = cuts[s + 1:min(s + 2048, len(cuts) - 1) + 1, None]
        rr = A + (B - A) * t
        ww = np.broadcast_to((B - A) * w, rr.shape).copy()
        for o in ons:
            L = B[:, 0] == o
            R = A[:, 0] == o
            rr[L] = B[L] - (B[L] - A[L]) * t * t
            ww[L] = (B[L] - A[L]) * 2 * t * w
            rr[R] = A[R] + (B[R] - A[R]) * t * t
            ww[R] = (B[R] - A[R]) * 2 * t * w
        mm = np.interp(rr, r, m)
        mm = np.where(rr < r[0], m[0] * (rr / r[0]) ** 3, mm)
        nv = np.interp(rr, r, nu)
        tot += np.sum(ww * 4 * np.pi * rr ** 2 / np.sqrt(1 - 2 * mm / rr) * baryon_density(mui * np.exp(-nv)),
                      dtype=np.longdouble)
    tot += 4 * np.pi / 3 * r[0] ** 3 * baryon_density(np.array([mui * np.exp(-nu[0])]))[0]
    return float(tot) * 1e54


def write(name, obj):
    (HERE / name).write_text(json.dumps(obj, indent=1, default=float) + '\n')
