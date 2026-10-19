#!/usr/bin/env python3
"""Candidate closure of frozen BA13; no historical fitter or solver mutation.

All spreads are estimates. No registry of owner acceptance is supplied. Numerical
algebra and measured smallness therefore cannot silently promote a conditional row.
"""
import csv
import gzip
import hashlib
import json
import math
import os
import pathlib
import resource
import subprocess
import sys
import time
from collections import Counter
import numpy as np

S=pathlib.Path(__file__).resolve().parents[3]
R=S.parents[1]
O=pathlib.Path(os.environ.get('BA13_CLOSURE_OUTPUT_ROOT',R/'ba13-closure-20261019'))
RESULT=O/'assessment-final'
HR=R/'ba13-headroom'
RAW=R/'ba13-recovery/run1/raw'
sys.path.insert(0,str(S/'tools/linux/ba13_headroom'))
from enclosure import FRACTIONS,L,budget,inverse_majorant,qualify,root_error,projection
sys.path.insert(0,str(S/'tools/linux'))
from ba13_refit import load,NAMES
LD=np.longdouble
TARGETS=(0,10,20)
WIDTHS=('0.0001','0.00005','0.000025')
MUON=('t_mu','N_mu','Z_row_npmu','Ltilde_Mmu','species_support')

def read(p,precise=False):
    return json.loads(p.read_text(),parse_float=LD if precise else float)

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def serial(v):
    if isinstance(v,np.ndarray):return v.tolist()
    if isinstance(v,np.floating):return float(v)
    if isinstance(v,np.integer):return int(v)
    raise TypeError(type(v).__name__)

def write(name,value):
    with (RESULT/name).open('x') as f:
        json.dump(value,f,indent=2,allow_nan=False,default=serial);f.write('\n')

def arr(v):return np.atleast_1d(np.asarray(v,dtype=LD))
def maximum(v):return float(np.max(np.abs(v)))

def fields(v):
    """Actual physical quantities in historical units (thermal values below)."""
    p=np.array(v['primary'],dtype=LD);g=np.array(v['G_packed'],dtype=LD)
    det=g[1]*g[3]-g[2]*g[2]
    z=np.array([[1/g[0]+g[3]/det,1/g[0]-g[2]/det],
                [1/g[0]-g[2]/det,1/g[0]+g[1]/det]],dtype=LD)
    k=LD(v['K']);gs=p[2]/p[1]**2/np.sqrt(1-2*p[2]/p[1])*LD('2.99792458e10')**2/LD('1e19')
    return dict(N_n=v['N'][0],N_e=v['N'][2],N_mu=v['N'][3],
                Ltilde_Me=p[4],Ltilde_Mmu=v['Lmu'],Cstar=p[5],radius=p[1],
                surface_gravity=gs,envelope=(LD('1e24')*gs*(LD('1.81e-8')*np.exp(LD(v['Hb'])-k))**LD('2.42'))**LD('.25'),
                Tsurface_inf=(LD('1e24')*gs*(LD('1.81e-8')*np.exp(LD(v['Hb'])-k))**LD('2.42'))**LD('.25')*np.exp(k-LD(v['Hs'])),
                mu_B=p[7],P0_average=v['P0'],P1_average=v['P1'],Z_row_npe=z[0],Z_row_npmu=z[1])

def band(a,b,c):
    """Prespecified order diagnostic, no floor for zero/rounding-limited spreads."""
    a,b,c=map(arr,(a,b,c));d1=b-a;d2=c-b
    p=[];valid=True
    for u,v in zip(d1,d2):
        if u==0 or v==0 or u*v<=0:
            p.append(None);valid=False
        else:
            order=float(np.log2(abs(u/v)));p.append(order)
            valid=valid and 3.5<=order<=4.5
    return dict(order=p,status='ORDER_COMPATIBLE' if valid else 'INCONCLUSIVE',
                coarse_tail=maximum(d1)*16/15,fine_tail=maximum(d2)/15,
                increments=[d1,d2],evidence_type='conditional_estimate',
                condition_id='gl4-asymptotic-remainder',authority=None)

def residual(xx,yy,uu):
    """Historical outer-half floor, all-sample residual, no covariance charge."""
    xx=np.asarray(xx,float);yy=np.asarray(yy,float);uu=np.asarray(uu,float)
    sig=np.maximum(uu[10:],64*np.finfo(float).eps*np.maximum(1.,abs(yy[10:])))
    w=(sig.min()/sig)**2;xm=np.sum(w*xx[10:])/w.sum();ym=np.sum(w*yy[10:])/w.sum()
    slope=np.sum(w*(xx[10:]-xm)*(yy[10:]-ym))/np.sum(w*(xx[10:]-xm)**2)
    return yy-(ym+slope*(xx-xm))

def main():
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Slurm required')
    os.sched_setaffinity(0,{min(os.sched_getaffinity(0))})
    RESULT.mkdir(exist_ok=True)
    start=time.process_time();wall=time.monotonic()
    x,series=load(RAW);x=np.array(x,dtype=LD)
    cert=read(R/'ba13-recovery/run1/fit/certificate.json')
    assert sha(S/'tests/bnv/compare_candidate.py')=='9b3ef1ce17df7e1f8d2d7c40e111ee0844abfddc2cc54b8f93cb3f2ce1322370'
    assert sha(R/'ba13-recovery/run1/fit/certificate.json')=='4212cad05bb4244ab75e80b8dfea023e6d0fcee69086c9993bd48d26a40e41fb'
    assert read(O/'receipt-m1.json')['status']=='PASS'
    profiles={i:read(HR/f'profile-v2-{i}.json') for i in range(21)}
    supp={v['index']:v for v in read(HR/'supplement.json')['records']}
    th=list(csv.DictReader((RAW/'thermal_surface.tsv').open(),delimiter='\t'))
    temp=np.array([LD(v['Tinf_K']) for v in th if int(v['index'])==0])
    Bp=np.array([LD(p['star']['B_solved']) for p in profiles.values()]);B0=Bp[0]
    Br=np.array([LD(p['N_nu'][0])+LD(p['N_nu'][1]) for p in profiles.values()])
    berr=(Br-Bp)/B0
    design_x=Br/Br[0]-1
    dx_design=design_x-x
    def norm(n,v):
        if n=='Cstar':return np.full(13,LD(v)/LD(profiles[0]['C_decache']),dtype=LD)
        if n in ('envelope','Tsurface_inf'):
            return LD(v)*temp**LD('.605')/np.array(series[n][2],dtype=LD)
        return arr(v)/np.array(series[n][2],dtype=LD)
    P={n:np.asarray(series[n][0],dtype=LD)/np.array(series[n][2],dtype=LD) for n in NAMES}
    recon={}
    for i,p in profiles.items():
        d=dict(N_n=p['N_nu'][0],N_e=p['N_nu'][2],N_mu=p['N_nu'][3],
               Cstar=p['C_nu'],Ltilde_Me=p['L_nu'],Ltilde_Mmu=p['Lmu_nu'],
               mu_B=p['mu_inf'],P0_average=p['mu_inf'],P1_average=supp[i]['integrals']['16']['P1'],
               Z_row_npe=p['Z_nu'][0],Z_row_npmu=p['Z_nu'][1])
        recon[i]={n:norm(n,v) for n,v in d.items()}
        for n in ('radius','surface_gravity'):recon[i][n]=P[n][i].copy()
        for n in ('envelope','Tsurface_inf'):
            # Ts_node/Ts_physical = exp(.605*(nu_physical-nu_node)).
            recon[i][n]=P[n][i]/(1+LD(p['Ts_offset']))
        recon[i]['mu_n_profile']=LD(p['mu_inf'])/np.array(series['mu_n_profile'][2],dtype=LD)
        recon[i]['metric_structure']=np.array(supp[i]['mapped_metric'],dtype=LD)/np.array(series['metric_structure'][2],dtype=LD)
    X={};raws={};m1={};derivatives={};roots={};x2_B={};reference={}
    for i in TARGETS:
        raws[i]={n:read((HR if n==512 else O)/f't{i}-n{n}.json',True) for n in (512,1024,2048)}
        X[i]={n:{k:norm(k,v) for k,v in fields(d).items()} for n,d in raws[i].items()}
        reference[i]={}
        derivatives[i]={}
        hs=LD(raws[i][512]['primary'][8])
        stencils={}
        for w in WIDTHS:
            parent=O if w==WIDTHS[-1] else HR
            lo=read(parent/f't{i}-h{w}--1.json',True);hi=read(parent/f't{i}-h{w}-1.json',True)
            a=hs*LD(w);bhi=LD(hi['primary'][0]);blo=LD(lo['primary'][0])
            t=(np.array(hi['N'],dtype=LD)-np.array(lo['N'],dtype=LD))/(bhi-blo)
            fhi=fields(hi);flo=fields(lo)
            stencils[w]=dict(lo=lo,hi=hi,a=a,B1=(bhi-blo)/(2*a),
                             B2=(bhi-2*LD(raws[i][512]['primary'][0])+blo)/a**2,t=t,
                             q1={n:(norm(n,fhi[n])-norm(n,flo[n]))/(2*a) for n in fhi},
                             q2={n:(norm(n,fhi[n])-2*X[i][512][n]+norm(n,flo[n]))/a**2 for n in fhi})
        narrow=stencils[WIDTHS[-1]]
        for j,n in ((0,'t_n'),(2,'t_e'),(3,'t_mu')):
            for resolution in X[i]:X[i][resolution][n]=arr(narrow['t'][j])
        x2_B[i]=LD(raws[i][2048]['primary'][0])
        for n in X[i][2048]:
            if n.startswith('t_'):
                j={'t_n':0,'t_e':2,'t_mu':3}[n]
                tv=[s['t'][j] for s in stencils.values()]
                delta1=tv[1]-tv[0];delta2=tv[2]-tv[1]
                ratio=float(delta1/delta2) if delta2 else None
                stencil_status='ORDER_COMPATIBLE' if ratio is not None and delta1*delta2>0 and 3<=ratio<=5 else 'INCONCLUSIVE'
                nb=band(*[raws[i][z]['N'][j] for z in (512,1024,2048)])
                bb=band(*[raws[i][z]['primary'][0] for z in (512,1024,2048)])
                bspan=abs(LD(narrow['hi']['primary'][0])-LD(narrow['lo']['primary'][0]))
                integration=2*(nb['coarse_tail']+abs(float(tv[-1]))*bb['coarse_tail'])/float(bspan)
                stencil=float(abs(delta2)/3)
                eoscenter={}
                for v in ('eos','center'):
                    vl=read(HR/f't{i}-derivative-{v}--1.json',True);vh=read(HR/f't{i}-derivative-{v}-1.json',True)
                    t=(LD(vh['N'][j])-LD(vl['N'][j]))/(LD(vh['primary'][0])-LD(vl['primary'][0]))
                    eoscenter[v]=float(abs(t-stencils[WIDTHS[1]]['t'][j]))
                reference[i][n]=dict(integration=integration,eos_equilibrium=eoscenter['eos'],center=eoscenter['center'],
                                     stencil=stencil,order=nb['order'],order_status=nb['status'],
                                     stencil_ratio=ratio,stencil_status=stencil_status,stencil_values=tv,
                                     band_complete=False,missing=['endpoint GL4-error transfer to offset H','implicit-stage and arithmetic propagation','EOS/center unseen remainders'])
            else:
                r=band(*[X[i][z][n] for z in (512,1024,2048)])
                ab={v:norm(n,fields(read(HR/f't{i}-{v}.json',True))[n]) for v in ('eos192','center-small','center-large')}
                reference[i][n]=dict(integration=r['fine_tail'],eos_equilibrium=maximum(ab['eos192']-X[i][512][n]),
                                     center=max(maximum(ab[v]-X[i][512][n]) for v in ('center-small','center-large')),
                                     stencil=0.,order=r['order'],order_status=r['status'],
                                     band_complete=False,missing=['implicit-stage and arithmetic propagation','EOS/center unseen remainders'])
            derivs=[]
            for w in WIDTHS:
                st=stencils[w];b1=st['B1'];b2=st['B2']
                if not n.startswith('t_'):
                    q1=st['q1'][n];q2=st['q2'][n]
                    derivs.append(dict(width=w,slope=q1*B0/b1,curvature=B0**2*(q2*b1-q1*b2)/b1**3))
            if n.startswith('t_'):
                j={'t_n':0,'t_e':2,'t_mu':3}[n]
                for wa,wb in zip(WIDTHS,WIDTHS[1:]):
                    a,b=stencils[wa],stencils[wb]
                    n1a=(LD(a['hi']['N'][j])-LD(a['lo']['N'][j]))/(2*a['a'])
                    n1b=(LD(b['hi']['N'][j])-LD(b['lo']['N'][j]))/(2*b['a'])
                    n1=(4*n1b-n1a)/3;n2=(LD(b['hi']['N'][j])-2*LD(raws[i][512]['N'][j])+LD(b['lo']['N'][j]))/b['a']**2
                    n3=6*(n1a-n1b)/(a['a']**2-b['a']**2)
                    b1=(4*b['B1']-a['B1'])/3;b2=b['B2'];b3=6*(a['B1']-b['B1'])/(a['a']**2-b['a']**2)
                    t1=n2/b1-n1*b2/b1**2
                    t2=n3/b1-n1*b3/b1**2-2*n2*b2/b1**2+2*n1*b2*b2/b1**3
                    derivs.append(dict(width_pair=[wa,wb],slope=arr(t1*B0/b1),curvature=arr(B0**2*(t2*b1-t1*b2)/b1**3)))
            derivatives[i][n]=derivs
        b1=stencils[WIDTHS[1]]['B1'];bfine=band(*[raws[i][z]['primary'][0] for z in (512,1024,2048)])
        residualB=abs(float(x2_B[i]-Bp[i]))
        roots[i]=dict(computed_B_residual=residualB,functional_GL4_tail_estimate=bfine['fine_tail'],
                      slope_secants=[float(stencils[w]['B1']) for w in WIDTHS],
                      local_slope_hypothesis=float(b1*(1-LD('.001'))),
                      H_error_estimate=root_error(residualB,bfine['fine_tail'],float(b1*(1-LD('.001')))),
                      evidence_type='conditional_estimate',condition_id='positive-B-slope-and-reference-band',authority=None,
                      missing='local derivative enclosure and EOS/center/stage/arithmetic B remainder')
    for n in X[0][2048]:
        m1[n]={}
        for i in (10,20):
            discrepancy=maximum((P[n][i]-P[n][0])-(X[i][2048][n]-X[0][2048][n]))
            bx=sum(reference[t][n][k] for t in (0,i) for k in ('integration','eos_equilibrium','center','stencil'))
            share=.23*(cert['samples'][i]['threshold'][n]-cert['samples'][i]['drift_bound'][n])
            orders_ok=all(reference[t][n]['order_status']=='ORDER_COMPATIBLE' and reference[t][n].get('stencil_status','ORDER_COMPATIBLE')=='ORDER_COMPATIBLE' for t in (0,i))
            gate='FAIL' if discrepancy>share else 'INCONCLUSIVE' if bx>=discrepancy or discrepancy+bx>share or not orders_ok else 'PASS_ESTIMATOR_ONLY'
            m1[n][i]=dict(paired_production_minus_X2=discrepancy,reference_partial_estimate=bx,
                         R_estimate=discrepancy+bx,representation_share=share,
                         estimator_gate=gate,reference_band_complete=False,qualification='UNQUALIFIED',
                         order_gate_satisfied=orders_ok,reason='reference remainders/evidence acceptance unresolved')
    # M2's full-grid prerequisite is checked from actual retained raw profiles.
    grids=[];grid0=None
    for i in range(21):
        p=RAW/f'chemical-profiles/target-{i}/profile.tsv'
        grid=np.loadtxt(p,skiprows=1,usecols=(0,))
        if grid0 is None:grid0=grid
        same=grid.shape==grid0.shape and np.array_equal(grid,grid0)
        match=grid.shape==grid0.shape
        first=int(np.flatnonzero(grid!=grid0)[0]) if match and not same else None
        grids.append(dict(index=i,nodes=len(grid),identical_to_0=same,first_different_node=first,
                          interior_identical_to_0=match and np.array_equal(grid[:-1],grid0[:-1]),
                          max_radial_difference_km=maximum(grid-grid0) if match else None,
                          r0=float(grid[0]),surface=float(grid[-1]),profile_sha256=sha(p),
                          grid_sha256=hashlib.sha256(grid.astype('<f8').tobytes()).hexdigest()))
    grid_ok=all(g['identical_to_0'] for g in grids)
    interior_grid_ok=all(g['interior_identical_to_0'] for g in grids)
    selection={key:dict(values=[profiles[i][key] for i in range(21)],
                        constant=all(profiles[i][key]==profiles[0][key] for i in range(21)))
               for key in ('k','neutron_cell','envelope_node')}
    edge=np.array([profiles[i]['r_nu'] for i in range(21)])
    a,b=profiles[0]['neutron_cell'];edge_margin=min(min(v-a,b-v) for v in edge)
    en=np.array([np.log1p(profiles[i]['Tb_offset']) for i in range(21)])
    ew=np.array([profiles[i]['envelope_cell_logT_jump'] for i in range(21)])
    node_margin=float(min(np.minimum(en,ew-en)))
    window=dict(physical_thresholds={i:{k:float(raws[i][2048][k]) for k in ('Hs','Hb','Hn','Hmu')} for i in TARGETS},
                Hc_values=[float(raws[i][2048]['primary'][8]) for i in TARGETS],
                neutron_cell_min_sample_margin_km=float(edge_margin),neutron_edge_sample_motion_km=float(np.ptp(edge)),
                neutron_edge_sample_monotone=bool(np.all(np.diff(edge)>0)),
                envelope_node_min_sample_margin=float(node_margin),envelope_sample_motion=float(np.ptp(en)),
                envelope_sample_monotone=bool(np.all(np.diff(en)>0) or np.all(np.diff(en)<0)),
                selection=selection,continuous_selection_proven=False,
                physical_response='no production-node jump charged to smooth physical response',
                production_response='sample monotonicity and margins do not bound unsampled adaptive-controller changes')
    window['sufficient_no_switch_conditions']=dict(
        neutron_edge_Lipschitz_limit_km_per_x=2*float(edge_margin)/float(max(abs(np.diff(x)))),
        neutron_edge_observed_secant_km_per_x=float(max(abs(np.diff(edge)/np.diff(x)))),
        envelope_Lipschitz_limit_per_x=2*node_margin/float(max(abs(np.diff(x)))),
        envelope_observed_secant_per_x=float(max(abs(np.diff(en)/np.diff(x)))),
        rationale='nearest achieved sample is at most half the largest gap away; L*gap/2 < margin suffices without monotonicity',
        derivative_bounds_established=False)
    # M5: independent endpoint error sets, full structure, both nonlinear tails.
    inv={};zscale={n:float(series[n][2][0]) for n in ('Z_row_npe','Z_row_npmu')}
    base=profiles[0]
    for i in range(21):
        endpoints=[]
        for j in (0,i):
            p=profiles[j];g=np.array(p['G_fixed_profile']);e=sum(np.array(p[k]) for k in ('background_inherited','center_inherited','tail_error'))
            aa,total,rem,rho=inverse_majorant(g,e)
            endpoints.append(dict(index=j,G=g,error_set=e,jacobian=aa,rho=rho,
                                  total=abs(L).T@total@abs(L),nonlinear=abs(L).T@rem@abs(L)))
        total=endpoints[0]['total']+endpoints[1]['total'] if i else np.zeros((2,2))
        nonlinear=endpoints[0]['nonlinear']+endpoints[1]['nonlinear'] if i else np.zeros((2,2))
        inv[i]=dict(endpoints=endpoints,total_Z_matrix=total,nonlinear_Z_matrix=nonlinear,
                    row_majorants={n:float(max(total[k]))/zscale[n] for k,n in enumerate(zscale)},
                    evidence_type='conditional_estimate',condition_id='per-target-chemical-background-magnitude',authority=None,
                    correlation='arbitrary independent endpoint error sets; identical point pairs exactly cancel',
                    additive_in_primary_route=False,reason='alternative structural/reference route; do not double count with direct production-X2 comparison')
    # Curvature hypothesis is evaluated on every achieved-B gap, including
    # uncertainty in coordinates estimated from the reconstruction and X2 B.
    curvature={};common={};fit={};transfer={}
    for n in NAMES:
        if n not in X[0][2048]:continue
        curv=[maximum(derivatives[i][n][-1]['curvature']) for i in TARGETS]
        spread=max(curv)/min(curv) if min(curv)>0 else None
        K=10*max(curv)
        sens=np.stack([derivatives[i][n][-1]['slope'] for i in TARGETS])
        si=np.stack([np.array([np.interp(j,TARGETS,np.asarray(sens[:,c],float)) for c in range(sens.shape[1])],dtype=LD) for j in range(21)])
        coord=si*berr[:,None]
        common[n]=dict(sensitivities=si,coordinate_source=berr,paired=np.max(abs(coord-coord[0]),axis=1),
                       x2_B_minus_production={i:float((x2_B[i]-Bp[i])/B0) for i in TARGETS},
                       evidence_type='conditional_estimate',condition_id='B-functional-and-sensitivity-remainder',authority=None)
        gaps=[]
        for j in range(20):
            width=float(abs(x[j+1]-x[j]))
            u=float(abs(dx_design[j])+abs(dx_design[j+1]))
            gaps.append(dict(left=j,right=j+1,achieved_spacing=width,coordinate_estimate=u,
                             tested_spacing=width+u,remainder_estimate=K*(width+u)**2/8,
                             evidence_type='conditional_estimate',condition_id='global-curvature-continuation',authority=None))
        wide_K=10*max(maximum(derivatives[i][n][0]['curvature']) for i in TARGETS)
        curvature[n]=dict(K_hypothesis=K,factor=10,three_target_curvatures=curv,
                          wider_stencil_K=wide_K,narrower_width_prediction_satisfied=all(v<=wide_K for v in curv),
                          cross_target_ratio=spread,estimator_gate='INCONCLUSIVE' if spread is None or spread>2 else 'PASS_ESTIMATOR_ONLY',
                          width_checks={i:derivatives[i][n] for i in TARGETS},intervals=gaps,
                          validated_supremum=False,missing='global derivative/continuation enclosure; fixed-cut H derivatives versus moving surface-cut response')
        if n.startswith('t_'):
            rep=np.zeros_like(P[n]);str_samples={i:P[n][i]-X[i][2048][n]-coord[i] for i in TARGETS}
        else:
            rep=np.stack([P[n][i]-recon[i][n] for i in range(21)])
            str_samples={i:recon[i][n]-X[i][2048][n]-coord[i] for i in TARGETS}
        mag=max(maximum(v) for v in str_samples.values())
        transfer[n]=dict(representation_signed=rep,structural_samples=str_samples,
                         independent_endpoint_magnitude=mag,grid_identity=grid_ok,
                         transfer_validated=False,evidence_type='conditional_estimate',condition_id='all-target-structural-transfer',authority=None,
                         interior_grid_identity=interior_grid_ok,
                         reason=('interior grids agree; only moving surface node differs; full-grid transfer precondition fails and structural-error supremum remains unvalidated'
                                 if interior_grid_ok and not grid_ok else 'complete grids differ; magnitude transfer unvalidated' if not grid_ok else 'grid identity alone does not establish an error-magnitude bound'),
                         common_B_removed_before_structural_estimate=True)
        vals,errs,scales,_=series[n];vals=np.asarray(vals,float);errs=np.asarray(errs,float);scales=np.asarray(scales,float)
        hfit=0.;design_change=0.;weight_change=0.;inner_excess=0.
        for c in range(vals.shape[1]):
            y=vals[:,c];u=errs[:,c];sc=scales[c]
            correction=np.array(-rep[:,c]*sc,float)
            old=residual(x,y,u);new=residual(x,y+correction,u)
            moved=residual(design_x,y+correction,u)
            hh=projection(np.asarray(x[10:],float),u[10:],y[10:])
            hfit=max(hfit,maximum((np.eye(11)-hh)@correction[10:])/sc)
            design_change=max(design_change,maximum(moved-new)/sc)
            # The difference includes the value-dependent sigma floor. The
            # fixed-projection value contribution is separately observable.
            weight_change=max(weight_change,maximum(new[10:]-old[10:]-(np.eye(11)-hh)@correction[10:])/sc)
            co=cert['fits'][n][c];r=abs(y-(co['intercept']+co['slope']*np.asarray(x,float)))/sc
            inner_excess=max(inner_excess,max(0.,float(max(r))-float(max(r[10:]))))
        fit[n]=dict(deterministic_projection=hfit,design_perturbation=design_change,
                    weight_floor_perturbation=weight_change,inner_half_excess=inner_excess,
                    estimate=hfit+design_change+weight_change+inner_excess,
                    shared_B0_coordinate_error=dx_design,fit_covariance_added=False,
                    evidence_type='conditional_estimate',condition_id='uncertain-design-and-source-remainder',authority=None,
                    missing='structural/reference source perturbation of residuals, design uncertainty remainder, and continuous fit remainder')
    rows=[]
    for n in NAMES:
        vals,errs,scales,kind=series[n]
        raw=np.max(abs(P[n]-P[n][0]),axis=1)
        for i in range(21):
            history={k:cert['samples'][i][k][n] for k in ('drift_bound','threshold','utilization')}
            hr=history['threshold']-history['drift_bound'];alloc=budget(hr)
            term={k:dict(value=None,evidence_type='unresolved',source=None,missing_estimate='quantity-specific evidence absent') for k in FRACTIONS}
            def put(k,value,source,condition=None,missing='',charged_to=None):
                term[k]=dict(value=value,evidence_type='unresolved' if value is None else 'conditional_estimate' if condition else 'analytic_bound',source=source,
                             missing_estimate=missing,condition_id=condition,authority=None,charged_to=charged_to)
            breakdown={}
            if n in transfer:
                tr=transfer[n];rep=tr['representation_signed']
                representation=maximum(rep[:i+1]-rep[0])
                structural=2*tr['independent_endpoint_magnitude'] if i else 0.
                put('production_representation',representation+structural,'P-R deterministic paired correction plus independent aligned R-X endpoints',
                    'all-target-structural-transfer',tr['reason'])
                breakdown=dict(deterministic_representation=representation,structural_independent_endpoints=structural,
                               physical_reference='X2 fixed-Hc N=2048; t uses third-width N=512 number stencil',
                               reconstruction_reference='none; production tangent compared directly to X2 number stencil' if n.startswith('t_') else 'fixed production profile; never independent stellar truth',
                               reference_band_complete=False,common_B_subtracted_from_R_minus_X=True)
                for k in ('integration','eos_equilibrium','center'):
                    val=max(reference[t][n][k] for t in TARGETS)*(2 if i else 0)
                    put(k,val,'X-Q reference contribution only; conditional all-target endpoint transfer',
                        'X2-reference-'+k, 'unseen remainder and quantity-dependent all-target transfer unvalidated')
                put('implicit_stages',None,'X2 reference stage error',None,'no propagated stage-residual/nonlinear remainder for this quantity')
                put('floating_point',None,'production binary64 reconstruction and X2 long-double reference kept separate',None,
                    'no enclosed binary64 elementary-function/reconstruction error or full reference arithmetic propagation; no 2^-64 production allowance')
                put('surface_events',None,'fixed finite-pressure surface and physical support boundaries',None,
                    'moving-cut derivative and event/reference remainder not established; finite-cut surface itself is the governed quantity')
                cb=maximum(common[n]['paired'][:i+1])
                put('common_B',cb,'B-functional discrepancy times component sensitivity; same signed source counted once',
                    'B-functional-and-sensitivity-remainder','unbounded functional and sensitivity remainder')
                if n.startswith('t_'):
                    put('normalization',0.,'absolute tangent component, S=1')
                    put('differentiation_inverse',2*max(reference[t][n]['stencil'] for t in TARGETS),
                        'three centered widths and actual attained B denominator','tangent-stencil-remainder',
                        'third-width h^2 test is numerical; no uniform paired stencil remainder')
                else:
                    u0=maximum(P[n][0]-X[0][2048][n])+sum(reference[0][n][k] for k in ('integration','eos_equilibrium','center'))
                    added=sum(t['value'] for t in term.values() if t['value'] is not None)
                    nv=(float(max(raw[:i+1]))+added)*u0/(1-u0) if u0<1 else None
                    put('normalization',nv,'shared B0 component denominator, original historical scale','B0-reference-normalization',
                        'reference absolute error and transferred endpoint magnitude remain conditional')
                    put('differentiation_inverse',0.,'complete nonlinear quantity is evaluated before constructing the reference band; M5 is an alternative, not an additive duplicate',charged_to='integration/production_representation')
                put('nonlinear_fit',fit[n]['estimate'],'all-21 residuals, historical projection, perturbed achieved-B design and weight floor',
                    'uncertain-design-and-source-remainder',fit[n]['missing'])
                bs=max((g['remainder_estimate'] for g in curvature[n]['intervals'][:i]),default=0.)
                put('between_samples',bs,'K=10 max three curvatures; all 20 achieved-B gaps evaluated',
                    'global-curvature-continuation',curvature[n]['missing'])
            elif n in ('metric_structure','mu_n_profile'):
                d=np.stack([recon[j][n]-P[n][j] for j in range(21)])
                put('production_representation',maximum(d[:i+1]-d[0]),'all-target deterministic support remapping/chemical identity reconstruction',
                    'nodewise-structural-reference','no independent X2 nodewise structural evidence at original support coordinates')
            if n in ('P2_average','species_support'):
                assert np.all(P[n]==0)
                if n=='P2_average':assert np.all(np.asarray(errs)==0)
                for k in FRACTIONS:put(k,0.,'exact original P2 identity' if n=='P2_average' else 'exact sampled species topology identity')
                if n=='species_support':
                    put('between_samples',None,'analytic threshold ordering at retained Hc',None,
                        'containing Hc interval/positive B slope and no-crossing continuation require explicit acceptance; sampled monotonicity is insufficient')
            q=qualify(history,term,exact=n in ('P2_average','species_support'),frozen_historical={k:cert['samples'][i][k][n] for k in history})
            unit='absolute dN_i/dB component (dimensionless count/count)' if n.startswith('t_') else 'historical scaled maximum-component drift'
            rows.append(dict(index=i,quantity=n,historical=history,headroom=hr,norm_kind=kind,units=unit,
                             component_count=len(scales),component_scales=dict(source='unchanged historical table parser',min=min(scales),max=max(scales)),
                             values_source=dict(historical=f'component-data.json.gz#/original_units/{n}/values/{i}',
                                                normalized_production=f'component-data.json.gz#/P/{n}/{i}',
                                                reconstruction=f'component-data.json.gz#/R/{i}/{n}' if n in recon[i] else None,
                                                independent_reference=f'component-data.json.gz#/X/{i}/{n}' if i in TARGETS and n in X[0][2048] else None),
                             allocated_budget=alloc,contributions=term,reference_decomposition=breakdown,
                             correlation='deterministic P-R signed pairs; aligned R-X independent endpoints; shared B0 coordinate counted once; no fit covariance duplication',
                             normalization='unchanged historical scales; Cstar same-knot cache factor cancels exactly; t has S=1',
                             alternative_inverse_majorant=inv[i]['row_majorants'].get(n),alternative_inverse_added=False,
                             evidence_standard_accepted=False,continuous_interval=[float(x[i]),0.],**q))
    assert len(rows)==462 and len({(r['index'],r['quantity']) for r in rows})==462
    summary=dict(schema_id='compactstar.ba13.closure.v1',entry_commit='286aafb63002093b3052e23a3972903011f36035',
                 predeclaration_commit='fdffaaf',candidate=True,owner_ratified=False,disposition='E',
                 counts={v:sum(r['status']==v for r in rows) for v in ('PASS','FAIL','UNQUALIFIED')},full_BA13_closure=False,
                 grid_identity=grid_ok,interior_grid_identity=interior_grid_ok,
                 muon_counts={v:sum(r['status']==v for r in rows if r['quantity'] in MUON) for v in ('PASS','FAIL','UNQUALIFIED')},
                 numerical_expansion_stopped=True,stop_reason='M2 transfer unproved; M1 reference and M3 global-continuation evidence incomplete',
                 remaining=['all-target structural-error transfer','X2 EOS/center/stage/arithmetic and tangent endpoint remainders',
                            'global curvature or accepted numerical continuation standard','uncertain-design and source-residual completion',
                            'metric and chemical-potential profile nodewise reference','five owner decisions'],
                 certificate_sha256=sha(R/'ba13-recovery/run1/fit/certificate.json'),fitter_sha256=sha(S/'tests/bnv/compare_candidate.py'))
    write('m1.json',dict(reference=reference,paired=m1,root_alignment=roots))
    write('m2.json',dict(grids=grids,grid_identity=grid_ok,interior_grid_identity=interior_grid_ok,transfer=transfer,qualification='UNQUALIFIED'))
    write('m3.json',dict(window=window,curvatures=curvature))
    write('m4.json',dict(common_B=common,fit=fit,reconstructed_B=Br,production_B=Bp,achieved_x=x,design_x=design_x))
    write('m5.json',inv)
    with gzip.open(RESULT/'component-data.json.gz','xt') as f:
        json.dump(dict(P=P,R=recon,X={i:X[i][2048] for i in TARGETS},
                       original_units={n:dict(scales=series[n][2],values=series[n][0]) for n in NAMES}),f,default=serial,allow_nan=False)
    write('matrix.json',dict(**summary,rows=rows));write('summary.json',summary)
    with (RESULT/'matrix.tsv').open('x') as f:
        cols=['index','quantity','historical_drift','threshold','headroom','units']+list(FRACTIONS)+['partial_utilization','proposed_total_utilization','status','missing']
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t');w.writeheader()
        for r in rows:
            w.writerow(dict(index=r['index'],quantity=r['quantity'],historical_drift=r['historical']['drift_bound'],threshold=r['historical']['threshold'],
                            headroom=r['headroom'],units=r['units'],**{k:t['value'] if t['value'] is not None else 'UNRESOLVED' for k,t in r['contributions'].items()},
                            partial_utilization=r['partial_utilization'],proposed_total_utilization=r['proposed_total_utilization'],status=r['status'],
                            missing='; '.join(k+': '+v for k,v in r['refusal_reasons'].items())))
    write('receipt-static.json',dict(cpu_seconds=time.process_time()-start,numerical_cpu_seconds=time.process_time()-start,build_cpu_seconds=0.,
          wall_seconds=time.monotonic()-wall,peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,job=os.environ['SLURM_JOB_ID'],
          affinity=sorted(os.sched_getaffinity(0)),source_sha256=sha(pathlib.Path(__file__))))
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
