#!/usr/bin/env python3
"""Common-q mapping diagnostics and target alignment. Spreads are not bounds."""
import csv,json,sys
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
S=Path(__file__).resolve().parents[3];R=Path(__import__('os').environ.get('BA13_STUDY_OUT',str(S.parents[1]/'ba13-uncertainty')))
sys.path.insert(0,str(S/'tests/analysis'));from chemical_trackr_budget import Model
q=np.linspace(0,1,2049);q2=np.linspace(0,1,4097)
def mapped(p):
 d=np.loadtxt(p/'profile/profile.tsv',skiprows=1);model=Model(p/'profile');r,m,nu,nb=d[:,:4].T
 lam=-.5*np.log1p(-2*m/r);pf=model.hc*np.cbrt(3*np.pi**2*d[:,6]);mu=np.exp(nu)*np.hypot(model.mass[0],pf)
 names=['r','m','nu','lambda','nB','nn','np','ne','nmu','mu_n_inf']
 fields=[r,m,nu,lam,nb,*d[:,6:10].T,mu];x=r/r[-1]
 output={};interp={};doubled={};raw_holdout={}
 for n,f in zip(names,fields):
  # Regular-center values, r=0,m=0,lambda=0; central nu/n remain finite.
  a=np.r_[0,x];v=np.r_[0 if n in ('r','m','lambda') else f[0],f]
  lin=np.interp(q,a,v);pch=PchipInterpolator(a,v)(q);scale=max(abs(v))
  output[n]=lin;interp[n]=float(np.max(abs(lin-pch))/scale)
  doubled[n]=float(np.max(abs(np.interp(q2,a,v)-PchipInterpolator(a,v)(q2)))/scale)
  keep=np.arange(0,len(a),2);keep=np.unique(np.r_[keep,len(a)-1]);omit=np.setdiff1d(np.arange(len(a)),keep)
  raw_holdout[n]=float(np.max(abs(np.interp(a[omit],a[keep],v[keep])-v[omit]))/scale) if len(omit) else 0.
 dims=d[:,5].astype(int);top=[int(dims[0])]+[int(dims[j]) for j in range(1,len(dims)) if dims[j]!=dims[j-1]]
 edges=np.interp(model.onsets,nb[::-1],r[::-1])/r[-1]
 return output,dict(interpolation_spread=interp,doubled_spread=doubled,withheld_nodes_spread=raw_holdout,topology=top,onset_q=edges.tolist(),nodes=len(r))
records=[];base={}
for i in (0,5,10,15,20):
 if (R/'runs'/f'nominal-{i}'/'COMPLETE.json').exists():base[i]=mapped(R/'runs'/f'nominal-{i}')
for p in sorted((R/'runs').glob('*')):
 if not (p/'COMPLETE.json').exists():continue
 i=int(p.name.rsplit('-',1)[1]);v,meta=mapped(p);nom,bmeta=base[i];st=next(csv.DictReader((p/'star.tsv').open(),delimiter='\t'));st={k:float(v) for k,v in st.items()}
 drift={k:float(np.max(abs(v[k]-nom[k]))/max(abs(nom[k]))) for k in v}
 alignment={k:abs(st[k+'_hi']-st[k+'_lo']) for k in ('R','M')}
 records.append(dict(case=p.name,index=i,common_q_normalized_drift=drift,**meta,
  target_residual_fraction=abs(st['residual'])/st['B_target'],bracket_fraction=st['bracket']/st['B_target'],rho_c=st['rho_c'],
  endpoint_alignment_spread=alignment,topology_matches_nominal=meta['topology']==bmeta['topology'],
  support_edge_shift=(np.array(meta['onset_q'])-bmeta['onset_q']).tolist(),status='DIAGNOSTIC_NO_GLOBAL_BOUND'))
(R/'profile-diagnostic.json').write_text(json.dumps(records,indent=2)+'\n');print('profiles',len(records),'topology matches',all(v['topology_matches_nominal'] for v in records))
