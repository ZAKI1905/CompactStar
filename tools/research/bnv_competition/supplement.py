#!/usr/bin/env python3
"""Resolve cancellation locality, representative periods and applicability."""
import argparse,json,os
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
from model import Model
from scan import write_csv

p=argparse.ArgumentParser();p.add_argument('out',type=Path);a=p.parse_args()
assert os.environ.get('SLURM_JOB_ID');m=Model()
drive=4*np.pi**2*1e-20/.005**3;critical=-2*m.I[0]*drive/m.t[0]/m.B*m.year
mu=m.I[1]*m.t[0]/(m.I[0]*m.t[1]);rc=m.diagnostics(drive,0)
rat=lambda b:m.surface(m.solve(drive,critical*b,beta_e=b)[0])/rc['Tsurface_inf_K']
minimum=minimize_scalar(rat,bounds=(1,mu),method='bounded',options={'xatol':1e-13})
betas=np.unique(np.r_[1.,mu,minimum.x,*[v+s*np.logspace(-12,-1,45) for v in [1.,mu] for s in [-1,1]]])
zoom=[]
for b in betas:
 r=m.diagnostics(drive,b*critical,beta_e=b);r['surface_ratio']=r['Tsurface_inf_K']/rc['Tsurface_inf_K']
 # Exact deposition required to restore the RC temperature, allowing the
 # stationary chemical imbalances to adjust at that temperature.
 target=m.drivers(drive,b*critical,b)[-1];eta=m.eta_at(rc['Tinf_K'],target)
 _,h,n,q=m.reactions(rc['Tinf_K'],eta)
 r['q_restore_RC_keV']=(m.photon(rc['Tinf_K'])+n+q-h)/(m.erg*r['D_count_s'])*1000
 zoom.append(r)
write_csv(a.out/'cancellation-zoom.csv',zoom)
representatives=[]
for period in [.0025,.005,.01,.02]:
 d=4*np.pi**2*1e-20/period**3;gc=-2*m.I[0]*d/m.t[0]/m.B*m.year
 base=m.diagnostics(d,0,P_s=period)
 for b in [0.,.1,.3,1.,mu,2.,3.,10.]:
  r=m.diagnostics(d,gc*b,P_s=period,beta_e=b);r['surface_ratio']=r['Tsurface_inf_K']/base['Tsurface_inf_K']
  r['rotation_10percent_diagnostic']=max(r['rotation_e'],r['rotation_mu'])<.1
  r['local_relaxation_under_1Gyr']=r['tau_slow_yr']<1e9
  r['driving_quasistatic_10percent_diagnostic']=r['relax_over_drive']<.1
  representatives.append(r)
write_csv(a.out/'representative-periods.csv',representatives)
# Along the QSS chemical manifold sample the thermal balance, including all
# changes of sign, as a check for an overlooked positive-temperature root.
thermal_checks=[]
for d in [1e-16,drive,1e-9]:
 gc=-2*m.I[0]*d/m.t[0]/m.B*m.year
 for b in [.1,1.,mu,3.,100.]:
  target=m.drivers(d,b*gc,b)[-1];signs=[]
  for T in np.logspace(0,10,301):
   eta=m.eta_at(T,target);_,h,n,q=m.reactions(T,eta)
   signs.append(np.sign(h-n-q-m.photon(T)))
  changes=int(sum(x!=y for x,y in zip(signs,signs[1:])));thermal_checks.append(dict(drive=float(d),beta_e=float(b),positive_T_sign_changes=changes))
assert all(r['positive_T_sign_changes']==1 for r in thermal_checks)
assert all(r['stable'] for r in zoom+representatives)
summary=dict(job=os.environ['SLURM_JOB_ID'],minimum_beta_e=float(minimum.x),minimum_surface_ratio=float(minimum.fun),
             zoom_points=len(zoom),representative_points=len(representatives),thermal_root_bracket_scan=thermal_checks,
             reachability='energy-norm build lower bound and local small-amplitude eigenmodes; no S3 or actual attainment claim',
             warning='Exact-null slow modes are much narrower than main map cells. Three local e-folds are a diagnostic, not a rigorous settling time or depletion bound.',
             stable=True)
(a.out/'supplement.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
