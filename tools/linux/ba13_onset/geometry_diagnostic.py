"""POST-HOC explanation of failed phase limits; never a repaired prospective PASS."""
import json,resource
from pathlib import Path
import numpy as np
S=Path(__file__).resolve().parents[3];R=S.parents[1];out=R/'ba13-onset-pilot';hp=json.loads((out/'constants.json').read_text());dis=json.loads((out/'discovery.json').read_text());result={}
for i in [0,10,20]:
 for c in ['nominal']+[f'r{n}' for n in range(80001,80006)]+['r640000']:
  o=dis[f'{i}:{c}'] if c=='nominal' else json.loads((out/'analysis'/f'{c}-{i}.json').read_text())
  root=R/('ba13-uncertainty' if c=='nominal' else 'ba13-onset-pilot');d=np.loadtxt(root/'runs'/f'{c}-{i}/profile/profile.tsv',skiprows=1,usecols=(0,1,2,4))
  k=o['k'];re=o['r_nu'];h=d[k,0]-d[k-1,0];theta=(re-d[k-1,0])/h;m=float(np.interp(re,d[:,0],d[:,1]));p=hp['P_on']*1.3235754613110692e-6;e=hp['eps_on']*1.3235754613110692e-6
  den=re*(re-2*m);A=m+4*np.pi*re**3*p;nup=A/den;mp=4*np.pi*re**2*e;pp=-(e+p)*nup
  nusecond=((mp+12*np.pi*re**2*p+4*np.pi*re**3*pp)*den-A*(2*re-2*m-2*re*mp))/den**2
  err=-.5*nusecond/nup*theta*(1-theta)*h*h
  result[f'{i}:{c}']=dict(theta=theta,h_km=h,nu_prime=nup,nu_second=nusecond,predicted_linear_inverse_error_km=err,original_edge_km=re,geometry_subtracted_edge_km=re-err,L_nu=o['L_nu'],geometry_subtracted_L=o['L_nu']*(1-o['edge_sensitivity']*err),status='POST_HOC_DIAGNOSTIC_NOT_VALIDATED_CORRECTION')
summary={}
for i in [0,10,20]:
 rows=[result[f'{i}:{c}'] for c in ['nominal']+[f'r{n}' for n in range(80001,80006)]]
 summary[str(i)]={k:float(np.ptp([r[k] for r in rows])/abs(rows[0][k])) for k in ['original_edge_km','geometry_subtracted_edge_km','L_nu','geometry_subtracted_L']}
print(json.dumps(summary,indent=2));(out/'geometry-diagnostic.json').write_text(json.dumps(dict(summary=summary,cases=result,status='POST_HOC_MECHANISM_DIAGNOSIS_FROZEN_P7_REMAINS_FAIL'),indent=2)+'\n')
u=resource.getrusage(resource.RUSAGE_SELF);(out/'geometry-resource.json').write_text(json.dumps(dict(cpu_seconds=u.ru_utime+u.ru_stime,max_rss_kib=u.ru_maxrss),indent=2)+'\n')
