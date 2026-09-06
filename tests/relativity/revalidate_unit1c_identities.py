#!/usr/bin/env python3
"""U3/U13: physical-CGS identities, mass measures and representation controls.
Usage: python revalidate_unit1c_identities.py TRACKR_EXPORT CMF_EXPORT OUT_JSON
Acceptance: existing 8e-15 algebra bound; Track-R's predeclared 1e-7/3e-7
accumulated bounds at 40000. For tabulated CMF, bound the mass integral by
monotone endpoint epsilon measures, not a tolerance inferred from its residual.
The CMF crust partition does not imply monotone convergence of interpolants.
"""
from pathlib import Path
import sys,json
sys.dont_write_bytecode = True
import numpy as np
import revalidate_background as b


def stats(lhs,rhs,r,mask):
    ii=np.flatnonzero(mask);err=np.abs(lhs-rhs);rel=err/np.maximum(abs(rhs),np.finfo(float).tiny)
    i=ii[np.argmax(err[ii])];j=ii[np.argmax(rel[ii])]
    return dict(max_absolute=float(err[i]),absolute_at_r_km=float(r[i]),max_relative=float(rel[j]),relative_at_r_km=float(r[j]),nodes=len(ii),roundoff_absolute_bound=float(36*np.finfo(float).eps*np.max(abs(rhs[ii]))),roundoff_relative_bound=float(36*np.finfo(float).eps))


def pointwise(a,track):
    r,m,e,p,nb,nu,nup=a[:,:7].T
    # CGS TOV right-hand sides evaluated independently of RelativityUnits.
    rc,grams,rho,pc=a[:,18:22].T
    G,c=b.G,b.c
    npcg=G/(rc*rc*c*c)*(grams+4*np.pi*rc**3*pc/c**2)/(1-2*G*grams/(rc*c*c))
    mg=G*grams/c**2/1e5;pg=G*pc/c**4*1e10
    dm=4*np.pi*rc*rc*rho*G/c**2*(m/mg)
    dp=-(rho*c*c+pc)*npcg*G/c**4*1e15*(p/pg)
    F=(m+4*np.pi*r**3*p)/(r*(r-2*m))
    pairs={'metric':(np.exp(-2*a[:,17]),1-2*m/r),'mass':(dm,4*np.pi*r*r*e),'hydrostatic':(dp,-(e+p)*F),'nuprime':(nup,F),'pressure_nuprime':(dp,-(e+p)*nup)}
    bands={'all':r>0,'center_side':r<.01,'bulk_core':(r>.01)&(rho>1e14),'near_surface':r>.999*r[-1]}
    if track:
        bands.update(muon_onset=abs(nb-.456984805412419)<.02,neutron_onset=(nb>3e-9)&(nb<2e-8),outer_pe=nb<7.35672890373299e-9)
    else:
        # Locate composition onset from the exported species, never transplant
        # Track-R thresholds into the crust EOS. Absent transitions are explicit.
        for name,k in [('muon_onset',16),('neutron_onset',13)]:
            jumps=np.flatnonzero((a[:-1,k]>0)!=(a[1:,k]>0))
            if len(jumps):
                mask=np.zeros(len(a),bool)
                for j in jumps:mask[max(0,j-2):min(len(a),j+4)]=True
                bands[name]=mask
        bands['outer_pe']=(a[:,13]==0)&(a[:,16]==0)&(a[:,14]>0)&(a[:,15]>0)
    out={n:{key:stats(lh,rh,r,mask) for key,(lh,rh) in pairs.items()} for n,mask in bands.items() if np.any(mask)}
    out['regimes_absent']=[n for n in ['muon_onset','neutron_onset','outer_pe'] if n not in out]
    return out


def integral(a,order):
    x,w=np.polynomial.legendre.leggauss(order);z=a[:-1,None,:]+(x+1)[None,:,None]/2*np.diff(a,axis=0)[:,None,:]
    wt=np.diff(a[:,0])[:,None]*w/2;r,m,e,p=z[:,:,:4].transpose(2,0,1)
    return np.cumsum(np.sum(wt*4*np.pi*r*r*e,axis=1)),np.cumsum(np.sum(wt*(e+p)*(m+4*np.pi*r**3*p)/(r*(r-2*m)),axis=1))


def measures(a):
    mass,pres=integral(a,8);m16,p16=integral(a,16)
    r=a[:,0];out={}
    for name,diff,scale,floor in [('mass',a[1:,1]-a[0,1]-mass,a[-1,1],max(abs(m16-mass))),('pressure',a[1:,3]-a[0,3]+pres,a[0,3],max(abs(p16-pres)))]:
        i=np.argmax(abs(diff));out[name]=dict(max_absolute=float(abs(diff[i])),relative_to_star_scale=float(abs(diff[i])/scale),at_r_km=float(r[i+1]),gauss_8_16_absolute_floor=float(floor),scale=float(scale))
    # Monotonic density on each physical segment bounds its integral exactly.
    # Center mass is the integrator's carried regular series, not extrapolation.
    assert np.all(np.diff(a[:,2])<=0),'nonmonotone epsilon invalidates endpoint bound'
    vol=4*np.pi/3*np.diff(r**3);lo=np.cumsum(vol*a[1:,2]);hi=np.cumsum(vol*a[:-1,2]);actual=a[1:,1]-a[0,1]
    numerical_floor=1e-8*a[-1,1] # 100 times canonical 1e-10 relative RK tolerance.
    out['endpoint_mass_measure']=dict(max_lower_violation=float(max(0,np.max(lo-actual))),max_upper_violation=float(max(0,np.max(actual-hi))),numerical_floor_absolute=float(numerical_floor),max_envelope_relative=float(np.max(hi-lo)/a[-1,1]))
    assert np.all(actual>=lo-numerical_floor) and np.all(actual<=hi+numerical_floor),'U3 mass measure failure'
    # Monotone p, epsilon and m give endpoint bounds on every positive TOV
    # source factor. This is independent of the linear interpolant quadrature.
    dr=np.diff(r);r0,r1=r[:-1],r[1:];m0,m1=a[:-1,1],a[1:,1]
    e0,e1=a[:-1,2],a[1:,2];p0,p1=a[:-1,3],a[1:,3]
    assert np.all(np.diff(a[:,1])>=0) and np.all(np.diff(a[:,3])<=0)
    denlo=r0*(r0-2*m1);denhi=r1*(r1-2*m0)
    assert np.all(denlo>0),"endpoint pressure bound requires positive denominators"
    plo=np.cumsum(dr*(e1+p1)*(m0+4*np.pi*r0**3*p1)/denhi)
    phi=np.cumsum(dr*(e0+p0)*(m1+4*np.pi*r1**3*p0)/denlo)
    pactual=a[0,3]-a[1:,3];pfloor=1e-8*a[0,3]
    out['endpoint_pressure_measure']=dict(max_lower_violation=float(max(0,np.max(plo-pactual))),max_upper_violation=float(max(0,np.max(pactual-phi))),numerical_floor_absolute=float(pfloor),max_envelope_relative=float(np.max(phi-plo)/a[0,3]))
    assert np.all(pactual>=plo-pfloor) and np.all(pactual<=phi+pfloor),"U3 pressure measure failure"
    return out


def main():
    td,cd,out=map(Path,sys.argv[1:]);result={'status':'IN PROGRESS','stars':{}}
    cases=[(td,n,True) for n in ['g8192r20000','g8192r40000']]+[(cd,'cmf-r'+str(n),False) for n in [10000,20000,40000,80000]]
    for d,n,t in cases:
        a=np.loadtxt(d/(n+'.tsv'),skiprows=1);z={'pointwise':pointwise(a,t),'integral':measures(a),'negative':{}}
        for mode in ['mass_only','energy_only','pressure_only','former_production']:
            neg=a.copy()
            if mode in ['mass_only','former_production']:neg[:,1]*=1.4766250380501249/(b.G*b.Ms/b.c**2/1e5)
            if mode in ['energy_only','former_production']:neg[:,2]*=6.67430e-8/b.G
            if mode in ['pressure_only','former_production']:neg[:,3]*=6.67430e-8/b.G
            neg[:,17]=-.5*np.log(1-2*neg[:,1]/neg[:,0])
            zz=pointwise(neg,t)['all'];z['negative'][mode]=zz
            assert max(v['max_relative'] for v in zz.values())>1e-6,'negative control did not discriminate'
        result['stars'][n]=z;out.write_text(json.dumps(result,indent=2))
        assert max(v['max_relative'] for v in z['pointwise']['all'].values())<8e-15,'U3 algebraic failure'
        if n=='g8192r40000':assert z['integral']['mass']['relative_to_star_scale']<1e-7 and z['integral']['pressure']['relative_to_star_scale']<3e-7
        print(n,'PASS',flush=True)
    result['status']='U3/U13 identity and mass-measure checks PASS; CMF pressure also passes independent monotone endpoint source bounds'
    out.write_text(json.dumps(result,indent=2))
if __name__=='__main__':main()
