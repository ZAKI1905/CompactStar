#!/usr/bin/env python3
"""Independent free-gas Cv/positive-sum/cache interpolation decomposition."""
import csv,json,math,sys
from pathlib import Path
import numpy as np
import mpmath as mp
from scipy.integrate import simpson
S=Path(__file__).resolve().parents[3];R=Path(__import__('os').environ.get('BA13_STUDY_OUT',str(S.parents[1]/'ba13-uncertainty')));OLD=S.parents[1]/'ba13-recovery'
KB=8.61733326214518e-5*1e-6;CONV=1.602176634e33;EPS=sys.float_info.epsilon
assert '8.61733326214518e-5' in (S.parents[1]/'sources/Zaki/Zaki/Physics/src/Constants.cpp').read_text()
Tgrid=np.array([math.exp((1-i/159)*math.log(1e-5)+(i/159)*math.log(1.)) for i in range(160)])
def factor(T):
 x=T*KB;k=int(np.searchsorted(Tgrid,x)-1);assert 0<=k<159
 w=(math.log(x)-math.log(Tgrid[k]))/(math.log(Tgrid[k+1])-math.log(Tgrid[k]))
 return ((1-w)*Tgrid[k]+w*Tgrid[k+1])/x
records=[];precision=[]
roots={'native':OLD/'run1/raw','mac':OLD/'historical-package/BA13-original-evidence-20261012/original/phase6a1-frozen-certificate-final'}
for platform,root in roots.items():
 thermal=list(csv.DictReader((root/'thermal_surface.tsv').open(),delimiter='\t'))
 for i in range(21):
  p=root/'chemical-profiles'/f'target-{i}';d=np.loadtxt(p/'profile.tsv',skiprows=1)
  meta=(p/'model.txt').read_text().splitlines();mass=np.array([float(l.split()[0]) for l in meta[:4]]);hc=float(meta[0].split()[1])
  r,m,nu=d[:,:3].T;ns=d[:,6:10]
  pf=hc*np.cbrt(3*np.pi**2*ns);density_slope=np.sum(pf*np.hypot(pf,mass),axis=1)/(3*hc**3)
  integrand=density_slope*np.exp(-nu)*4*np.pi*r*r/np.sqrt(1-2*m/r)*1e15*KB*KB*CONV
  terms=(integrand[:-1]+integrand[1:])*.5*np.diff(r)
  slope=math.fsum(terms);simpson_slope=simpson(integrand,x=r)
  gamma=len(terms)*EPS/(1-len(terms)*EPS)
  for k,row in enumerate(thermal[i*13:(i+1)*13]):
   T=float(row['Tinf_K']);actual=float(row['Cstar_erg_K']);f=factor(T);exact_linear=slope*T;predicted=exact_linear*f
   records.append(dict(platform=platform,index=i,knot=k,T=T,Cstar=actual,nodes=len(r),analytic_nodal_trapezoid=exact_linear,
    cache_factor=f,cache_relative_bias=f-1,predicted_cache=predicted,producer_minus_prediction_relative=(actual-predicted)/actual,
    sequential_sum_gamma=gamma,arithmetic_bound_fixed_nodal_terms=(gamma+128*EPS)*abs(predicted),
    original_uncertainty=float(row['u_Cstar_erg_K']),simpson_minus_trapezoid_relative=(simpson_slope-slope)/slope,
    absolute_accuracy_status='UNBOUNDED_STRUCTURE_EOS_INTERPOLATION',
    correlation='Cache factor shared by all stars at each T; cancels in exact linear-source relative depletion ratios. Node/arithmetic/background errors need separate treatment.'))
  if platform=='native' and i in (0,10,20):
   for j in np.linspace(0,len(r)-1,9,dtype=int):
    answers=[]
    for prec in (80,112):
     with mp.workdps(prec):
      ps=[mp.mpf(hc)*(3*mp.pi**2*mp.mpf(float(n)))**(mp.mpf(1)/3) if n else mp.mpf(0) for n in ns[j]]
      ans=mp.fsum(pj*mp.sqrt(pj*pj+mp.mpf(float(ma))**2) for pj,ma in zip(ps,mass))/(3*mp.mpf(hc)**3)
      answers.append(float(ans))
    assert answers[0]==answers[1]
    diff=abs(density_slope[j]-answers[0]);bound=128*EPS*answers[0];assert diff<=bound
    precision.append(dict(index=i,node=int(j),relative_difference=diff/answers[0],bound_relative=128*EPS,status='PASS_FIXED_INPUT_ARITHMETIC'))
report=dict(status='MEASURED_NOT_FULL_ACCURACY_BOUND',records=records,high_precision=precision,
 constants=dict(kB_MeV_K=KB,MeV_fm3_erg_cm3=CONV),
 derivation='Cv/Tinf = kB^2*MeVfm3_to_cgs * integral sum(pF*sqrt(m^2+pF^2))/(3 hc^3) exp(-nu) dV; entropy adapter linear in T',
 caveat='gamma bounds positive summation on fixed terms; 128 epsilon covers tested term arithmetic, not a rigorous libm or structural bound')
(R/'thermal-diagnostic.json').write_text(json.dumps(report,indent=2)+'\n')
print('max cache bias',max(x['cache_relative_bias'] for x in records),'max prediction discrepancy',max(abs(x['producer_minus_prediction_relative']) for x in records))
print('fixed input arithmetic heldout checks',len(precision))
