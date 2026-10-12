#!/usr/bin/env python3
"""Independent BA13 table parser and 80/112-digit normal-equation comparator.

Does not import the recovered fitter. Arithmetic checks do not validate the
physical provenance or sufficiency of its assigned numerical uncertainties.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import sys
import mpmath as mp

NAMES=('t_n','t_e','t_mu','Z_row_npe','Z_row_npmu','Cstar','Ltilde_Me',
       'Ltilde_Mmu','mu_B','mu_n_profile','P0_average','P1_average','P2_average',
       'N_n','N_e','N_mu','species_support','metric_structure','radius',
       'surface_gravity','envelope','Tsurface_inf')
EPS=sys.float_info.epsilon

def table(p):
    with p.open() as f:return [{k:float(v) for k,v in r.items()} for r in csv.DictReader(f,delimiter='\t')]

def load(raw):
    t=table(raw/'target_stars.tsv');d=table(raw/'tangent.tsv');c=table(raw/'coefficients.tsv')
    h=table(raw/'thermal_surface.tsv');p=table(raw/'profile_metrics.tsv')
    assert [len(v) for v in (t,d,c,h,p)]==[21,21,21,273,43029]
    for rows in (t,d,c):assert [r['index'] for r in rows]==list(range(21))
    b=t[0]['B_solved'];x=[(r['B_solved']-b)/b for r in t]
    for i in range(21):
        assert d[i]['B_solved']==c[i]['B_solved']==t[i]['B_solved']
        assert abs(t[i]['residual'])<=5e-11*b and 0<=t[i]['bracket_width']<=5e-11*b
    hg=[[r for r in h if r['index']==i] for i in range(21)]
    pg=[[r for r in p if r['index']==i] for i in range(21)]
    assert all(len(g)==13 for g in hg) and all(len(g)==2049 for g in pg)
    result={}
    def add(name,values,errors,kind='relative',scale=None):
        result[name]=(values,errors,scale or [max(abs(v),1e-300) for v in values[0]],kind)
    def scalar(name,rows,key,error=None,kind='relative',scale=None):
        add(name,[[r[key]] for r in rows],[[r[error] if error else 0.] for r in rows],kind,scale)
    for name in ('t_n','t_e','t_mu'):scalar(name,d,name,'u_'+name,'t',[1.])
    for name in ('N_n','N_e','N_mu'):scalar(name,t,name,'u_'+name,'number')
    for name in ('Ltilde_Me','Ltilde_Mmu'):scalar(name,c,name,'u'+name)
    for name,key,err in [('mu_B','mu_B_inf','u_mu_B'),('P0_average','P0_inf','u_P0'),('P1_average','P1_inf','u_P1'),('P2_average','P2_inf','u_P2')]:
        scalar(name,c,key,err,scale=[1.] if name=='P2_average' else None)
    scalar('radius',t,'radius_km')
    for j,name in enumerate(('Z_row_npe','Z_row_npmu')):
        v=[[r[f'Z{j}{k}'] for k in range(2)] for r in c]
        add(name,v,[[r[f'uZ{j}{k}'] for k in range(2)] for r in c],'z',[max(map(abs,v[0]))]*2)
    for name,key,err in [('Cstar','Cstar_erg_K','u_Cstar_erg_K'),('envelope','Tsurface_local_K',None),('Tsurface_inf','Tsurface_inf_K','u_Tsurface_inf_K')]:
        add(name,[[r[key] for r in g] for g in hg],[[r[err] if err else 0. for r in g] for g in hg])
    scalar('surface_gravity',[g[0] for g in hg],'surface_g14','u_surface_g14')
    add('mu_n_profile',[[r['mu_n_inf_MeV'] for r in g] for g in pg],[[0.]*2049 for _ in pg])
    keys=('r_km','m_km','nu','lambda')
    add('metric_structure',[[r[k] for k in keys for r in g] for g in pg],[[0.]*8196 for _ in pg],
        scale=[max(max(abs(r[k]) for r in pg[0]),1e-300) for k in keys for _ in range(2049)])
    def topology(g):
        v=[int(r['active_dimension']) for r in g];assert min(v)>0
        return [v[i] for i in range(len(v)) if i==0 or v[i]!=v[i-1]]
    add('species_support',[[0. if topology(g)==topology(pg[0]) else 2.] for g in pg],[[0.] for _ in pg],'topology',[1.])
    assert set(result)==set(NAMES)
    return x,result

def limit(kind,value,u):
    if kind=='t':return max(5*u,1e-4*abs(value))
    if kind=='z':return max(5*u,1e-4)
    if kind=='number':return 1e-4
    if kind=='topology':return 1.
    return 1e-4+5*u

def solve(x,y,u,dps):
    with mp.workdps(dps):
        x=list(map(mp.mpf,x));y=list(map(mp.mpf,y));u=list(map(mp.mpf,u))
        w=[1/(z*z) for z in u]
        s=mp.fsum(w);sx=mp.fsum(a*b for a,b in zip(w,x));sxx=mp.fsum(a*b*b for a,b in zip(w,x))
        sy=mp.fsum(a*b for a,b in zip(w,y));sxy=mp.fsum(a*b*c for a,b,c in zip(w,x,y))
        det=s*sxx-sx*sx;assert det>0
        a=(sxx*sy-sx*sxy)/det;b=(s*sxy-sx*sy)/det
        cov=[[sxx/det,-sx/det],[-sx/det,s/det]]
        residual=[v-a-b*z for z,v in zip(x,y)]
        ur=[mp.sqrt(e*e+cov[0][0]+2*z*cov[0][1]+z*z*cov[1][1]) for z,e in zip(x,u)]
        return [float(a),float(b)],[[float(v) for v in r] for r in cov],list(map(float,residual)),list(map(float,ur))

def main():
    if sys.flags.optimize:raise RuntimeError('independent validation requires active assertions')
    p=argparse.ArgumentParser();p.add_argument('--raw',type=Path,required=True);p.add_argument('--certificate',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();x,series=load(a.raw);cert=json.loads(a.certificate.read_text())
    assert cert['outer_half_indices']==list(range(10,21))
    assert [s['achieved_abscissa'] for s in cert['samples']]==x
    records=[];quantities={};max_arithmetic=0.;failures=[]
    for name in NAMES:
        values,errors,scales,kind=series[name]
        drift=[max(abs(v-v0)/s for v,v0,s in zip(row,values[0],scales)) for row in values]
        uncertainty=[max(math.hypot(e,e0)/s for e,e0,s in zip(row,errors[0],scales)) for row in errors]
        thresholds=[limit(kind,values[i][0],uncertainty[i]) for i in range(21)]
        nonlinear=0.; historical_nonlinear=0.;zero_errors=sum(e==0 for row in errors for e in row)
        for j,scale in enumerate(scales):
            xx=x[10:];yy=[row[j] for row in values[10:]]
            uu=[max(row[j],64*EPS*max(1.,abs(y))) for row,y in zip(errors[10:],yy)]
            first=solve(xx,yy,uu,80);second=solve(xx,yy,uu,112)
            assert first==second,'precision instability'
            co,cv,res,ur=first;f=cert['fits'][name][j]
            # Conservative forward arithmetic allowance for eleven-term binary64
            # sums, weighted centering, and prediction cancellation. It is used
            # only to compare algorithms, never added to a scientific threshold.
            sy=max(max(map(abs,yy)),1e-300);spread=max(xx)-min(xx)
            allowance=1024*EPS*sy
            for actual,ref,bound in zip((f['intercept'],f['slope']),co,(allowance,allowance/spread)):
                ratio=abs(actual-ref)/bound;max_arithmetic=max(max_arithmetic,ratio)
                if ratio>1:failures.append([name,j,'coefficient arithmetic',ratio])
            for k in range(2):
                for l in range(2):
                    if abs(f['covariance'][k][l]-cv[k][l])>1024*EPS*max(abs(cv[k][l]),1e-300):failures.append([name,j,'covariance'])
            for k,(rr,uu_res,check) in enumerate(zip(res,ur,f['checks'])):
                i=k+10;rn=abs(rr)/scale;un=uu_res/scale
                if rn>3*un+.1*thresholds[i]:failures.append([name,j,'independent scientific residual',i])
                if abs(check['residual_normalized']-rn)>allowance/scale:failures.append([name,j,'residual arithmetic',i])
                if abs(check['uncertainty_normalized']-un)>4096*EPS*max(un,1e-300):failures.append([name,j,'residual covariance',i])
                if check['index']!=i or check['allowed']!=3*check['uncertainty_normalized']+.1*thresholds[i]:failures.append([name,j,'historical criterion',i])
                nonlinear=max(nonlinear,rn);historical_nonlinear=max(historical_nonlinear,check['residual_normalized'])
            records.append(dict(quantity=name,component=j,coefficients=co,covariance=cv,residual=res,u_res=ur))
        envelope=0.;maxutil=0.
        for i in range(21):
            envelope=max(envelope,drift[i]);bound=envelope+historical_nonlinear
            sample=cert['samples'][i]
            assert sample['drift_bound'][name]==bound and sample['threshold'][name]==thresholds[i]
            assert sample['utilization'][name]==bound/thresholds[i]
            independent_util=(envelope+nonlinear)/thresholds[i];maxutil=max(maxutil,independent_util)
            if independent_util>1:failures.append([name,'independent envelope',i])
        quantities[name]=dict(components=len(scales),zero_error_entries=zero_errors,error_entries=21*len(scales),max_independent_utilization=maxutil,
            largest_numerical_uncertainty=max(uncertainty),nonlinearity=nonlinear)
    report=dict(status='PASS_INDEPENDENT_HISTORICAL_SEMANTICS_REFIT' if not failures else 'FAIL',
        precision_digits=[80,112],precision_repeat='exact binary64 outputs',fit_count=len(records),failures=failures,
        method='separate raw parser; uncentered 2x2 high-precision normal equations',
        arithmetic_allowance='1024 epsilon times value scale; slope additionally divided by abscissa span; comparison only',
        max_coefficient_arithmetic_allowance_utilization=max_arithmetic,quantities=quantities,
        certificate_sha256=hashlib.sha256(a.certificate.read_bytes()).hexdigest(),
        scientific_qualification=False,limitation='Reproduces historical zero-error assignments and envelope formula; does not validate omitted errors')
    with a.output.open('x') as f:json.dump(report,f,indent=2,allow_nan=False);f.write('\n')
    with gzip.open(str(a.output)+'.fits.json.gz','wt') as f:json.dump(records,f,allow_nan=False)
    print(report['status'],'fits',len(records),'failures',len(failures))
    if failures:raise SystemExit(1)

if __name__=='__main__':main()
