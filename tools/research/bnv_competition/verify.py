#!/usr/bin/env python3
"""Independent polynomial equations, Newton roots and Jacobians for S1.

This is a computational cross-check, not an independent human scientific review.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import math
import mpmath as mp
import numpy as np
from scipy.optimize import root, brentq
from model import Model


def main():
    p = argparse.ArgumentParser(); p.add_argument('out', type=Path); a = p.parse_args()
    assert os.environ.get('SLURM_JOB_ID')
    m = Model(); mp.mp.dps = 55
    table = list(csv.DictReader((a.out/'stationary-full.csv').open()))
    fid = list(csv.DictReader((a.out/'fiducial.csv').open()))
    checks = {}; records = []
    null_mu = m.I[1]*m.t[0]/(m.I[0]*m.t[1])
    selected = []
    for drive in [1e-16, 4*np.pi**2*1e-20/.005**3, 1e-9]:
        for beta in [.1, .3, 1., null_mu, 2., 3., 100.]:
            selected.append(min(table, key=lambda r: abs(math.log(float(r['drive_rad2_s3'])/drive))+abs(float(r['beta_e'])-beta)))
    # Independent powers in dimensional eta,T, expanded rather than H(xi).
    def independent(T, e, u):
        rates, nu = [], mp.mpf(0)
        for L, eta in zip(m.L, [e,u]):
            rate = mp.mpf(L)/mp.mpf(m.kberg)*sum(mp.mpf(c)*eta**k*T**(7-k)/(11513*mp.pi**(k+1)*mp.mpf(m.kb)**k)
                     for c,k in [(14680,1),(7560,3),(840,5),(24,7)])
            rates.append(rate)
            nu += mp.mpf(L)*(T**8+sum(mp.mpf(c)*eta**k*T**(8-k)/(11513*mp.pi**k*mp.mpf(m.kb)**k)
                      for c,k in [(22020,2),(5670,4),(420,6),(9,8)]))
        heat = mp.mpf(m.erg)*(e*rates[0]+u*rates[1])
        photon = mp.mpf(m.Lg0)*(T/mp.mpf(1e8))**mp.mpf(m.power_exponent)
        return rates, heat, nu, photon

    for row in selected:
        drive, gamma, beta = [float(row[k]) for k in ['drive_rad2_s3','Gamma_yr','beta_e']]
        # Recover exact designated coordinates from stored gamma/drive metadata.
        if abs(beta-1) < 1e-14: beta = 1.
        if abs(beta-null_mu) < 1e-14: beta = null_mu
        _, spin, bn, target = m.drivers(drive,gamma,beta)
        scale = max(abs(spin).max(), abs(bn).max())
        T = float(row['Tinf_K']); eta = np.array([float(row['eta_e_MeV']),float(row['eta_mu_MeV'])])
        power_scale = float(row['LH_erg_s'])
        def equations(z, e, u):
            rates, heat, nu, photon = independent(mp.exp(z)*T,e*.02,u*.02)
            return ((rates[0]-mp.mpf(target[0]))/scale,
                    (rates[1]-mp.mpf(target[1]))/scale,
                    (heat-nu-photon)/power_scale)
        initial = (mp.mpf(0),mp.mpf(eta[0])/.02,mp.mpf(eta[1])/.02)
        sol = mp.findroot(equations, initial, tol=mp.mpf('1e-43'), maxsteps=60)
        ref = np.array([float(mp.exp(sol[0])*T),float(sol[1]*.02),float(sol[2]*.02)])
        diff = max(abs(ref-np.r_[T,eta])/np.array([T,.02,.02]))
        tightT,tighte = m.solve(drive,gamma,xtol=2e-14,beta_e=beta)
        convergence = max(abs(np.r_[tightT,tighte]-np.r_[T,eta])/np.array([T,.02,.02]))
        # Independent mpmath differentiation of the original vector field.
        def vector(z,e,u):
            TT=mp.exp(z)*T; ee=e*.02; uu=u*.02
            rates, heat, nu, photon=independent(TT,ee,uu)
            rr=mp.matrix([rates[0]-target[0],rates[1]-target[1]])
            ed=-mp.matrix(m.Z.tolist())*rr
            return ((heat-nu-photon)/(m.capacity(T)*TT),ed[0]/.02,ed[1]/.02)
        Ji=mp.matrix(3)
        args=(mp.mpf(0),mp.mpf(eta[0])/.02,mp.mpf(eta[1])/.02)
        for i in range(3):
            for j in range(3):
                Ji[i,j]=mp.diff(lambda v: vector(*[v if k==j else args[k] for k in range(3)])[i],args[j])
        S=np.array([T,.02,.02]);Ja=m.jacobian(T,eta)*S[None,:]/S[:,None]
        jdiff=float(max(abs(Ji[i,j]-Ja[i,j]) for i in range(3) for j in range(3))/max(abs(v) for v in Ji))
        vals=mp.eig(Ji,left=False,right=False);realmax=max(mp.re(v) for v in vals)
        tau=float(-1/realmax/m.year)
        tauerr=abs(tau/float(row['tau_slow_yr'])-1)
        # A different binary64 solver, perturbed seed, at nonsingular states.
        scipy_error = None
        if min(abs(eta)) > 1e-6:
            def fun(v): return np.array([float(q) for q in equations(*v)])
            answer=root(fun, np.array([.03,eta[0]/.02*1.05,eta[1]/.02*.95]),method='hybr',tol=1e-10)
            assert answer.success and max(abs(fun(answer.x)))<1e-8
            scipy_error=max(abs(answer.x-np.array([float(s) for s in sol])))
        records.append(dict(drive=drive,beta_e=beta,root_error_scaled=diff,tightening_error_scaled=convergence,
                            jacobian_relative_error=jdiff,stable=bool(realmax<0),
                            independent_tau_yr=tau,tau_relative_difference=tauerr,hybr_error=scipy_error))
    checks['independent_55_digit_roots']=max(r['root_error_scaled'] for r in records)<1e-9
    checks['solver_tightening']=max(r['tightening_error_scaled'] for r in records)<1e-9
    checks['independent_jacobian']=max(r['jacobian_relative_error'] for r in records)<1e-9
    checks['independent_stability']=all(r['stable'] for r in records)
    checks['slow_eigenvalue_accuracy']=max(r['tau_relative_difference'] for r in records)<1e-6
    checks['alternative_hybr']=max(r['hybr_error'] or 0 for r in records)<1e-7
    # Exact odd/even and channel-null limits at fixed T.
    et=np.array([.02,-.03]);ra,ha,na,qa=m.reactions(2e6,et);rb,hb,nb,qb=m.reactions(2e6,-et)
    checks['rate_odd_power_even']=np.array_equal(ra,-rb) and ha==hb and na==nb
    checks['zero_drive_boundary_only']=m.solve(0.,0.)[0]==0
    checks['electron_exact_null']=all(float(r['eta_e_MeV'])==0 for r in table if abs(float(r['beta_e'])-1)<1e-14)
    checks['muon_exact_null']=all(float(r['eta_mu_MeV'])==0 for r in table if abs(float(r['beta_e'])-null_mu)<1e-14)
    checks['mixed_sign_window']=all(float(r['eta_e_MeV'])<0<float(r['eta_mu_MeV']) for r in table if 1.00001<float(r['beta_e'])<null_mu-1e-5)
    # Conventional RC power scaling over the photon dominated, large-xi range.
    controls=list(csv.DictReader((a.out/'controls.csv').open()))[1:]
    fit=[r for r in controls if 0<float(r['drive_rad2_s3'])<=1e-11]
    x=np.log([float(r['drive_rad2_s3']) for r in fit])
    slopes={k:float(np.polyfit(x,np.log([float(r[k]) for r in fit]),1)[0]) for k in ['Lgamma_erg_s','Tsurface_inf_K','eta_e_MeV','tau_slow_yr']}
    checks['RC_asymptotic_exponents']=abs(slopes['Lgamma_erg_s']-8/7)<.01 and abs(slopes['Tsurface_inf_K']-2/7)<.003
    # Thresholds defined by 10% surface suppression and return to RC luminosity.
    drive=4*np.pi**2*1e-20/.005**3
    critical=-2*m.I[0]*drive/m.t[0]/m.B*m.year
    rc=m.diagnostics(drive,0)
    ratio=lambda beta:m.surface(m.solve(drive,critical*beta,beta_e=beta)[0])/rc['Tsurface_inf_K']
    thresholds={'surface_10_percent_suppression_beta':brentq(lambda b:ratio(b)-.9,.01,.9),
                'return_to_RC_beta':brentq(lambda b:ratio(b)-1,1.2,3)}
    checks['dimensions']=abs(m.kberg/(m.kb*m.erg)-1)<1e-15
    checks={k:bool(v) for k,v in checks.items()}
    result=dict(status='PASS' if all(checks.values()) else 'FAIL',checks=checks,representatives=records,
                RC_scaling_slopes=slopes,thresholds=thresholds,job=os.environ['SLURM_JOB_ID'],
                scope='Independent equations/solver/differentiation; no independent human review or certificate extension.')
    (a.out/'verification.json').write_text(json.dumps(result,indent=2,default=lambda x:x.item())+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='representatives'},indent=2))
    assert all(checks.values()), checks


if __name__=='__main__': main()
