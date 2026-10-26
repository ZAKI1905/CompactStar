"""Reference assembly; reconstruction R is never an independent reference."""
import csv
import gzip
import hashlib
import json
import pathlib
import sys
import numpy as np

S = pathlib.Path(__file__).resolve().parents[3]
R = S.parents[1]
O = R/'ba13-evidence-20261026'
OLD = R/'ba13-closure-20261019'
F = OLD/'assessment-final'
RAW = R/'ba13-recovery/run1/raw'
LD = np.longdouble
sys.path.insert(0, str(S/'tools/linux/ba13_closure'))
from assess import fields


def read(p):
    return json.loads(p.read_text())


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def arr(v):
    return np.atleast_1d(np.asarray(v, dtype=LD))


def load_reference():
    with gzip.open(F/'component-data.json.gz', 'rt') as f:
        comp = json.load(f)
    m1, m3, m4 = [read(F/(n+'.json')) for n in ('m1', 'm3', 'm4')]
    revpath = R/'ba13-closure-independent-review-20261020/outputs/reference-diagnostics.json'
    diag = read(revpath)
    x = np.array(m4['achieved_x'], dtype=LD)
    refs = {n: np.stack([arr(comp['X'][str(i)][n]) for i in (0, 10, 20)])
            for n in comp['X']['0']}
    # Tangent reference changes from N512 to measured N2048, still at the same
    # frozen narrow width. This does not rescore any historical row.
    tg = read(O/'tangents.json')
    for n in ('t_n','t_e','t_mu'):
        refs[n] = np.array([[LD(tg['records'][f'{i}:{n}']['t']['2048'])]
                            for i in (0,10,20)])
    K = {n: max(m3['curvatures'][n]['K_hypothesis'],
                m3['curvatures'][n]['wider_stencil_K'],
                10*diag['B_curvature_and_interpolated_reference'][n]['target_curvature_max'])
         for n in refs}
    return comp, m1, m3, m4, diag, x, refs, K


def interpolation(x, knots, values, *, reference_kind='X2'):
    if reference_kind != 'X2':
        raise ValueError('X2 self-error cannot qualify a reconstruction reference')
    order = np.argsort(knots)
    knots = np.asarray(knots, dtype=LD)[order]
    values = np.asarray(values, dtype=LD)[order]
    x = np.atleast_1d(np.asarray(x, dtype=LD))
    if np.any(x < knots[0]) or np.any(x > knots[-1]):
        raise ValueError('reference extrapolation forbidden')
    j = np.clip(np.searchsorted(knots, x, side='right')-1, 0, len(knots)-2)
    w = (x-knots[j])/(knots[j+1]-knots[j])
    return values[j]+w[:,None]*(values[j+1]-values[j]), j, w


def normalized_fields(data, comp):
    values = fields(data)
    temperatures = [LD(v['Tinf_K']) for v in csv.DictReader(
        (RAW/'thermal_surface.tsv').open(), delimiter='\t') if int(v['index'])==0]
    c0 = LD(read(R/'ba13-headroom/profile-v2-0.json')['C_decache'])
    out = {}
    for n,v in values.items():
        scales = arr(comp['original_units'][n]['scales'])
        if n == 'Cstar':
            out[n] = np.full(13,LD(v)/c0,dtype=LD)
        elif n in ('envelope','Tsurface_inf'):
            out[n] = LD(v)*np.array(temperatures)**LD('.605')/scales
        else:
            out[n] = arr(v)/scales
    return out


def paired_comparison(production, reference, root_band, *, reference_kind='X2'):
    if reference_kind != 'X2':
        raise ValueError('wrong independent reference')
    if root_band is None or not np.isfinite(root_band) or root_band < 0:
        raise ValueError('missing common-B alignment')
    error = np.asarray(production,dtype=LD)-np.asarray(reference,dtype=LD)
    paired = np.maximum.accumulate(np.max(abs(error-error[0]),axis=1))
    return error, paired, float(root_band)
