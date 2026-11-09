#!/usr/bin/env python3
"""Retained-state C-CUT derivatives; no stellar solver. Scientific candidate.

Outputs are sensitivities, not automatically finite-error bounds. See the
predeclaration and derivation in docs/validation/results-first-bnv-pilot/cut-assessment.
"""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import mpmath as mp
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
DOC = ROOT / 'docs/validation/results-first-bnv-pilot/cut-assessment'
WIDTHS = ('0.0001', '0.00005', '0.000025')
TARGETS = (0, 10, 20)


def read(path):
    with (gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)) as f:
        return json.load(f, parse_float=str)


def wide(value):
    # Match X2's physical constants: binary64 constants widened to Real.
    return mp.mpf(float(value))


def eos(H):
    """Exact zero-temperature pe branch, independent closed Fermi integrals."""
    proton, electron = wide('938.27208816'), wide('0.51099895000')
    hc = wide(1 / float('5.067730716156395e-3'))
    h0 = wide(float(proton) + float(electron))
    h = h0 * mp.exp(H)
    dh = h0 * mp.expm1(H)
    k2 = dh*(h+h0)*(h-proton+electron)*(h+proton-electron)/(4*h*h)
    k = mp.sqrt(k2)
    n = k**3 / (3*mp.pi**2*hc**3)
    pressure, energy, chis, mus = mp.mpf(0), mp.mpf(0), [], []
    for mass in (proton, electron):
        mu = mp.sqrt(mass**2+k2)
        mus.append(mu)
        energy += (k*mu*(2*k2+mass**2)-mass**4*mp.asinh(k/mass))/(8*mp.pi**2*hc**3)
        pressure += (k*mu*(2*k2-3*mass**2)+3*mass**4*mp.asinh(k/mass))/(24*mp.pi**2*hc**3)
        chis.append(mu*k/(mp.pi**2*hc**3))
    chi = chis[0]*chis[1]/sum(chis)
    # X2 widens the rounded binary64 sum h0 separately from the two masses.
    # Differentiate that literal convention; do not erase its tiny EOS defect.
    k_H = (h*h-h0*h0*(proton-electron)**2/(h*h))/(4*k)
    p_H = sum(k**4/(3*mp.pi**2*hc**3*mu) for mu in mus)*k_H
    n_H = k*k*k_H/(mp.pi**2*hc**3)
    return dict(p=pressure, e=energy, n=n, cv=mp.pi**2*sum(chis)/3,
                chi=chi, h=h, hc=hc, k=k, masses=(proton, electron),
                p_H=p_H,n_H=n_H,chemical_enthalpy=sum(mus))


def terminal(H, R, M, K):
    q = eos(H)
    geo = wide(1.602176634e33*6.673e-8/float(2.99792458e10)**4*1e10)
    D = 1-2*M/R
    V = M/R**3 + 4*mp.pi*q['p']*geo
    A = D/V
    alpha = q['p']/q['p_H']
    W = 4*mp.pi*R*A/mp.sqrt(D)
    R_s = -A/R*alpha
    M_s = 4*mp.pi*R**2*q['e']*geo*R_s
    B_s = -W*alpha*q['n']*mp.mpf('1e54')
    K_s = 4*mp.pi*(q['e']+q['p'])*geo/V*alpha
    g_log_s = M_s/M-2*R_s/R-K_s+alpha
    cvfactor = mp.mpf('1e15')*wide('8.61733326214518e-11')**2*wide('1.602176634e33')
    return dict(R=R_s, M=M_s, B=B_s, K=K_s, alpha=alpha,
        g_log=g_log_s, C_shell=-W*alpha*mp.exp(H-K)*q['cv']*cvfactor,
        G_shell=-W*alpha*mp.exp(H-K)*q['chi']*mp.mpf('1e54'),
        thermodynamic_residual=(q['e']+q['p'])/(q['n']*q['h'])-1,
        literal_EOS_derivative_factor=(q['e']+q['p'])/q['p_H'],
        K_identity_residual=K_s-((-M_s/R+M*R_s/R**2)/D+alpha),
        pressure_cgs=q['p']*wide('1.602176634e33'))


def fields(raw):
    p = list(map(mp.mpf, raw['primary']))
    g = list(map(mp.mpf, raw['G_packed']))
    K, Hs, Hb = (mp.mpf(raw[n]) for n in ('K', 'Hs', 'Hb'))
    det = g[1]*g[3]-g[2]**2
    inv = mp.matrix([[g[3]/det, -g[2]/det], [-g[2]/det, g[1]/det]])
    z = inv + mp.ones(2)/g[0]
    gravity = p[2]/p[1]**2/mp.sqrt(1-2*p[2]/p[1])*mp.mpf('2.99792458e10')**2/mp.mpf('1e19')
    env = (mp.mpf('1e24')*gravity*(mp.mpf('1.81e-8')*mp.exp(Hb-K))**mp.mpf('2.42'))**mp.mpf('.25')
    values = dict(B=p[0], mass=p[2], radius=p[1], Cstar=p[5], Ltilde_Me=p[4],
        Ltilde_Mmu=mp.mpf(raw['Lmu']), mu_B=p[7], P0_average=mp.mpf(raw['P0']),
        P1_average=mp.mpf(raw['P1']), N_n=mp.mpf(raw['N'][0]),
        N_e=mp.mpf(raw['N'][2]), N_mu=mp.mpf(raw['N'][3]), surface_gravity=gravity,
        envelope=env, Tsurface_inf=env*mp.exp(K-Hs),
        Z_row_npe_0=z[0, 0], Z_row_npe_1=z[0, 1],
        Z_row_npmu_0=z[1, 0], Z_row_npmu_1=z[1, 1])
    surf = terminal(Hs, p[1], p[2], K)
    fixedH = dict(B=surf['B'], mass=surf['M'], radius=surf['R'],
        Cstar=-surf['K']*p[5]+surf['C_shell'],
        Ltilde_Me=-6*surf['K']*p[4], Ltilde_Mmu=-6*surf['K']*values['Ltilde_Mmu'],
        N_n=mp.mpf(0), N_mu=mp.mpf(0), N_e=surf['B'],
        surface_gravity=gravity*surf['g_log'],
        envelope=env*(mp.mpf('.25')*surf['g_log']-mp.mpf('.605')*surf['K']),
        Tsurface_inf=values['Tsurface_inf']*(mp.mpf('.25')*surf['g_log']+mp.mpf('.395')*surf['K']-surf['alpha']))
    for n in ('mu_B', 'P0_average', 'P1_average'):
        fixedH[n] = surf['K']*values[n]
    for name, i in (('Z_row_npe', 0), ('Z_row_npmu', 1)):
        for j in range(2):
            fixedH[f'{name}_{j}'] = surf['K']*z[i, j]-surf['G_shell']*inv[i, 0]*inv[0, j]
    return values, fixedH, surf


def stencil(center, lo, hi):
    """Differentiate the quadratic at the actual retained central H."""
    hm = mp.mpf(center['primary'][8])-mp.mpf(lo['primary'][8])
    hp = mp.mpf(hi['primary'][8])-mp.mpf(center['primary'][8])
    if hm <= 0 or hp <= 0:
        raise ValueError('central H outside stencil')
    return (lambda q0, qm, qp: (hm**2*(qp-q0)+hp**2*(q0-qm))/(hp*hm*(hp+hm)),
            lambda q0, qm, qp: 2*((qp-q0)/hp-(q0-qm)/hm)/(hp+hm))


def mixed_response(BH, Bs, BsH, Bss, NH, NHH, BHH, shell, B0, slope):
    t = NH/BH
    tH = (NHH-t*BHH)/BH
    ts = (shell-t)*BsH/BH-Bs*tH/BH
    path = t+slope/B0*(shell-t)*Bs
    path_s = ts+slope/B0*((shell-t)*(Bss-Bs*BsH/BH)-ts*Bs)
    return dict(fixed_cut=t, path=path, fixed_cut_s=ts, path_s=path_s,
                H_derivative=tH, mixed_count=NHH-t*BHH)


def path_response_from_onlaw(T, Tx, Bs, Bsx, B0, slope, shell):
    den=1-slope*Bs/B0
    return -Tx*Bs/(B0*den)+(shell-T)*Bsx/(B0*den**2)


def response(center, lo, hi, B0, slope):
    v, fh, surf = fields(center)
    vl, _, _ = fields(lo)
    vh, _, _ = fields(hi)
    d1, d2 = stencil(center, lo, hi)
    H = {n: d1(v[n], vl[n], vh[n]) for n in v}
    HH = {n: d2(v[n], vl[n], vh[n]) for n in v}
    if H['B'] <= 0:
        raise ValueError('nonpositive B_H')
    fixedB = {n: fh[n]-H[n]*surf['B']/H['B'] for n in v}
    Hs, R, M, K = mp.mpf(center['Hs']), v['radius'], v['mass'], mp.mpf(center['K'])
    # Derivatives of the closed terminal baryon response, not count differences.
    f = lambda h, r, m: terminal(h, r, m, K)['B']
    dH, dR, dM = (mp.diff(f, (Hs, R, M), order) for order in ((1,0,0), (0,1,0), (0,0,1)))
    B_sH = dR*H['radius']+dM*H['mass']
    B_ssH = dH*surf['alpha']+dR*surf['R']+dM*surf['M']
    Hc_x=(B0-slope*surf['B'])/H['B']
    Bsx=dR*(H['radius']*Hc_x+slope*surf['R'])+dM*(H['mass']*Hc_x+slope*surf['M'])+dH*surf['alpha']*slope
    tangents = {}
    for n, species, shell in (('t_n','N_n',0), ('t_e','N_e',1), ('t_mu','N_mu',0)):
        tangents[n] = mixed_response(H['B'],surf['B'],B_sH,B_ssH,H[species],
                                    HH[species],HH['B'],shell,B0,slope)
        t=tangents[n]['fixed_cut']
        TH=tangents[n]['H_derivative']*(1-slope*surf['B']/B0)+slope/B0*(shell-t)*B_sH
        TsH=(shell-t)*B_sH/H['B']*(1-slope*surf['B']/B0)+slope/B0*(shell-t)*B_ssH
        Tx=TH*Hc_x+TsH*slope
        alt=path_response_from_onlaw(tangents[n]['path'],Tx,surf['B'],Bsx,B0,slope,shell)
        tangents[n]['onlaw_identity_residual']=alt-tangents[n]['path_s']
        if abs(alt-tangents[n]['path_s'])>mp.mpf('1e-55'):
            raise ValueError('on-law mixed response identity')
    return dict(values=v, fixed_H_s=fh, H_derivative=H, fixed_B_s=fixedB,
        Hc_s=-surf['B']/H['B'], terminal=surf, tangents=tangents,
        baryon_identity_residual=(fixedB['N_n']+fixedB['N_e']+fixedB['N_mu']))


def diagnostic(values, formal):
    a, b, c = values
    u, v = b-a, c-b
    order = mp.log(abs(u/v), 2) if u and v else None
    ok = u*v > 0 and order is not None and abs(order-formal) <= mp.mpf('.5')
    return dict(increments=[u, v], order=order, formal_order=formal,
                status='ORDER_COMPATIBLE' if ok else 'INCONCLUSIVE',
                formal_tail=abs(v)/(2**formal-1) if ok else None)


def serialize(value):
    if isinstance(value, mp.mpf):
        return mp.nstr(value, 50)
    raise TypeError(type(value).__name__)


def calculate(inputs, digits):
    with mp.workdps(digits):
        B0 = mp.mpf(inputs['cuts'][0]['B'])
        # Original reviewer E_surface: binary64 np.polyfit of log geometric P
        # at the retained achieved coordinates. No new surface law is fitted.
        retained_x=read(ROOT/'docs/validation/ba13-closure/evidence/m4.json')['achieved_x']
        xf=np.array(retained_x,dtype=float)
        pf=np.array([s['Pgeom'] for s in inputs['cuts']],dtype=float)
        yf=np.log(pf)
        coeff=np.polyfit(xf,yf,1)
        slope=wide(coeff[0])
        residual=wide(np.max(np.abs(yf-np.polyval(coeff,xf))))
        expected=mp.mpf(read(ROOT/'docs/validation/ba13-evidence/evidence/continuation.json')['moving_surface']['lnPs_linear_fit_residual_max'])
        if abs(residual-expected)>mp.mpf('1e-26'):
            raise ValueError('original pressure-law replay differs')
        x=[wide(v) for v in xf]
        y=[mp.log(wide(v)) for v in pf]
        xm, ym = sum(x)/len(x), sum(y)/len(y)
        hp_slope=sum((a-xm)*(b-ym) for a,b in zip(x,y))/sum((a-xm)**2 for a in x)
        if residual > mp.mpf('2.1e-12'):
            raise ValueError('retained pressure law exceeds C-CUT')
        records, convergence, normalized = {}, {}, {}
        comp = read(ROOT/'docs/validation/ba13-evidence/evidence/comparison.json.gz')
        for target in TARGETS:
            ti = str(target)
            c = inputs['centers'][ti]
            records[ti] = {}
            for width in WIDTHS:
                endpoints = inputs['widths'][ti][width]
                records[ti]['w'+width] = response(c['512'], endpoints['-1'], endpoints['1'], B0, slope)
            for grid in ('512','1024','2048'):
                endpoints = inputs['grids'][ti][grid]
                records[ti]['n'+grid] = response(c[grid], endpoints['-1'], endpoints['1'], B0, slope)
            finest = records[ti]['n2048']
            for n, detail in comp.items():
                keys = ([n+'_0',n+'_1'] if n.startswith('Z_row_') else
                        ['mu_B'] if n == 'mu_n_profile' else [n])
                conv = []
                norm = []
                for j, key in enumerate(keys):
                    def get(rec):
                        return rec['tangents'][key]['path_s'] if key.startswith('t_') else rec['fixed_B_s'][key]
                    conv.append(dict(component=j,
                        width=diagnostic([get(records[ti]['w'+w]) for w in WIDTHS], 2),
                        grid=diagnostic([get(records[ti]['n'+g]) for g in ('512','1024','2048')], 4)))
                    if key.startswith('t_'):
                        norm.append(abs(get(finest)))
                    else:
                        # Reference values already carry the exact historical component scales,
                        # including each thermal knot and equilibrium profile node.
                        ref = detail['reference_values'][TARGETS.index(target)]
                        refscale = max(abs(mp.mpf(v)) for v in ref) if len(keys)==1 else abs(mp.mpf(ref[j]))
                        norm.append(abs(get(finest)/finest['values'][key])*refscale)
                convergence[f'{target}:{n}'] = conv
                normalized[f'{target}:{n}'] = dict(sensitivity=max(norm),
                    linearized_product=max(norm)*mp.mpf('2.1e-12'), finite_allowance=None)
        return dict(digits=digits, slope=slope, law_intercept_geom=wide(coeff[1]),
                    high_precision_slope_diagnostic_only=hp_slope,
                    law_residual=residual, records=records,
                    convergence=convergence, normalized=normalized)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args=parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('EKU Slurm required')
    manifest=read(DOC/'inputs-manifest.json')
    for name, record in manifest.items():
        path=ROOT/name
        if hashlib.sha256(path.read_bytes()).hexdigest()!=record['sha256']:
            raise ValueError('changed input '+name)
    inputs=read(DOC/'retained-inputs.json.gz')
    result=calculate(inputs,70)
    higher=calculate(inputs,90)
    max_rel=mp.mpf(0)
    for key,value in result['normalized'].items():
        other=higher['normalized'][key]['sensitivity']
        max_rel=max(max_rel,abs(value['sensitivity']-other)/max(abs(other),mp.mpf('1e-100')))
    if max_rel > mp.mpf('1e-45'):
        raise ValueError('postprocessing precision disagreement')
    args.output.mkdir()
    (args.output/'response.json').write_text(json.dumps(result,default=serialize,indent=2)+'\n')
    groups={}
    for name in read(ROOT/'docs/validation/ba13-evidence/evidence/comparison.json.gz'):
        rec=[result['normalized'][f'{t}:{name}'] for t in TARGETS]
        failures=[f'{t}:{j}:{axis}' for t in TARGETS
                  for j,r in enumerate(result['convergence'][f'{t}:{name}']) for axis in ('width','grid')
                  if r[axis]['status']!='ORDER_COMPATIBLE']
        groups[name]=dict(max_sampled_normalized_sensitivity=max(r['sensitivity'] for r in rec),
                         max_linearized_product=max(r['linearized_product'] for r in rec),
                         inconclusive_order_checks=failures, finite_cut_allowance=None)
    summary=dict(status='SENSITIVITY_ASSESSMENT_ONLY',cut_closed=False,
        precision_relative_disagreement=max_rel, groups=groups,
        missing='validated finite-cut transfer/remainder, including correlated mixed tangent response',
        stellar_calculations=0,trajectories=0)
    (args.output/'summary.json').write_text(json.dumps(summary,default=serialize,indent=2)+'\n')
    print(json.dumps(summary,default=serialize))


if __name__=='__main__':
    main()
