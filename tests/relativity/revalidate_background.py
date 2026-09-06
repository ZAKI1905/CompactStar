#!/usr/bin/env python3
"""ADR-0012 U3/U12/U13 independent test-side revalidation.

Usage: python revalidate_background.py BACKGROUND_DIRECTORY
Run relativistic_unit_background DIR --campaign first. Requires NumPy/SciPy.
No production conversion imports, response API or baseline writes. Homogeneous
ODE and proper-volume quadrature follow the accepted Unit-0 methodology.
"""
from pathlib import Path
import sys,json
import numpy as np
from scipy.integrate import solve_ivp
D=Path(sys.argv[1]) if __name__ == "__main__" else Path("."); G=6.673e-8; c=29979245800.; Ms=1.98892e33
x,w=np.polynomial.legendre.leggauss(8)

def load(name): return np.loadtxt(D/(name+'.tsv'),skiprows=1)
def mixed(a):
    a=a.copy()
    # Deliberate rejected representation, test-side only. Nominal IAU GM and
    # modern G are independent negative inputs; never positive conversion owners.
    beta=(1.3271244e26/c**2/1e5)/(G*Ms/c**2/1e5)
    alpha=6.67430e-8/G
    old_surface=.5*np.log(1-2*beta*a[-1,1]/a[-1,0])
    a[:,1]*=beta; a[:,2:4]*=alpha
    a[:,5]+=old_surface-a[-1,5]
    a[:,17]=-.5*np.log(1-2*a[:,1]/a[:,0])
    return a

def quad(a):
    z=a[:-1,None,:]+(x+1)[None,:,None]/2*np.diff(a,axis=0)[:,None,:]
    return z,np.diff(a[:,0])[:,None]*w/2

def count(a):
    z,wt=quad(a); r=z[:,:,0];m=z[:,:,1]
    W=4*np.pi*r*r/np.sqrt(1-2*m/r)
    return np.einsum('ij,ijk->k',wt*W,z[:,:,13:17])+4*np.pi*a[0,0]**3/3*a[0,13:17]

def homogeneous(a, rtol=2e-10, atol=1e-13, partitions=20000):
    r=a[:,0];ec,pc,dc=a[0,[2,3,7]]
    # Regular homogeneous response normalized to d ln(epsilon_c)=1.
    init=[4*np.pi/3*ec*r[0]**3,ec/(ec+pc)/dc]
    def rhs(t,y):
        k=np.clip(np.searchsorted(r,t)-1,0,len(r)-2)
        b=a[k]+(t-r[k])/(r[k+1]-r[k])*(a[k+1]-a[k])
        m,e,p,nup=b[[1,2,3,6]];es=(a[k+1,2]-a[k,2])/(r[k+1]-r[k]);mh,ph=y
        return [-4*np.pi*t*t*ph/nup*es,-mh*(1+8*np.pi*t*t*p)/(t-2*m)**2-4*np.pi*(e+p)*t*t*ph/(t-2*m)]
    sol=solve_ivp(rhs,(r[0],r[-1]),init,method='DOP853',rtol=rtol,atol=atol,t_eval=r,max_step=(r[-1]-r[0])/partitions)
    assert sol.success
    h=sol.y.T;z,wt=quad(a)
    hz=h[:-1,None,:]+(x+1)[None,:,None]/2*np.diff(h,axis=0)[:,None,:]
    rr=z[:,:,0];m=z[:,:,1];np_=z[:,:,6];ns=z[:,:,13:17]
    W=4*np.pi*rr*rr/np.sqrt(1-2*m/rr);xi=hz[:,:,1]/np_
    dn=np.diff(a[:,13:17],axis=0)/np.diff(r)[:,None]
    density=np.einsum('ij,ik->k',-wt*W*xi,dn)
    metric=np.einsum('ij,ijk->k',wt*W*hz[:,:,0]/(rr-2*m),ns)
    shell=4*np.pi*r[-1]**2/np.sqrt(1-2*a[-1,1]/r[-1])*h[-1,1]/a[-1,6]*a[-1,13:17]
    return density+metric+shell

def identities(a):
    r,m,e,p,nb,nu,nup=a[:,:7].T;rr,grams,rho,pc=a[:,18:22].T
    mg=G*grams/c**2/1e5;eg=G*rho/c**2*1e10;pg=G*pc/c**4*1e10
    F=(m+4*np.pi*r**3*p)/(r*(r-2*m));fg=(mg+4*np.pi*r**3*pg)/(r*(r-2*mg))
    dm=4*np.pi*r*r*eg*(m/mg);dp=-(eg+pg)*fg*(p/pg)
    pairs={'mass':(dm,4*np.pi*r*r*e),'hydrostatic':(dp,-(e+p)*F),'nuprime':(nup,F),'pressure_nuprime':(dp,-(e+p)*nup),'metric':(np.exp(-2*a[:,17]),1-2*m/r)}
    bands={'all':np.ones(len(r),bool),'center':r<.01,'bulk':(r>=.01)&(nb>.46),'muon':abs(nb-.456984805412419)<.02,'neutron':(nb>3e-9)&(nb<2e-8),'outer_pe':nb<7.35672890373299e-9,'surface':r>.999*r[-1]}
    out={b:{key:float(np.max(abs((lhs-rhs)[mask]/rhs[mask]))) for key,(lhs,rhs) in pairs.items()} for b,mask in bands.items() if np.any(mask)}
    z,wt=quad(a);r,m,e,p=z[:,:,:4].transpose(2,0,1)
    mass=np.sum(wt*4*np.pi*r*r*e,axis=1)
    hydro=np.sum(wt*(e+p)*(m+4*np.pi*r**3*p)/(r*(r-2*m)),axis=1)
    # Epsilon is linear per profile segment. Gauss quadrature is exact for the
    # polynomial r^2*epsilon source; hydrostatic quadrature is independently refined.
    out['integral_mass']=float(np.max(abs(a[1:,1]-a[0,1]-np.cumsum(mass)))/a[-1,1])
    out['integral_pressure']=float(np.max(abs(a[1:,3]-a[0,3]+np.cumsum(hydro)))/a[0,3])
    return out

if __name__ == "__main__":
    result={'species':['n','p','e','mu'],'status':'IN PROGRESS'}
    a=load('g8192r40000');result['U3']=identities(a);result['U13_mixed']=identities(mixed(a))
    assert max(result['U3']['all'].values())<8e-15
    # Predeclared from Unit-0 quadrature/refinement: <1e-7 mass, <3e-7 pressure.
    assert result['U3']['integral_mass']<1e-7 and result['U3']['integral_pressure']<3e-7
    assert result['U13_mixed']['all']['mass']>2e-4 and result['U13_mixed']['all']['nuprime']>5e-5
    result['sequence']={}
    for g in [4096,8192]:
        counts={j:count(load(f'seq{g}_{j}')) for j in [-4,-2,-1,1,2,4]}
        B=(8*(counts[1]-counts[-1])-(counts[2]-counts[-2]))/.006
        B2=(8*(counts[2]-counts[-2])-(counts[4]-counts[-4]))/.012
        result['sequence'][str(g)]={'Bx_h0005':B.tolist(),'Bx_h001':B2.tolist(),'step_relative':(B2/B-1).tolist()}
        assert np.max(abs(B2/B-1))<2e-4
        if g==8192:
            for mode in ['production_A1','negative_mixed']:
                aa=a if mode=='production_A1' else mixed(a)
                BB=B if mode=='production_A1' else (8*(count(mixed(load('seq8192_1')))-count(mixed(load('seq8192_-1'))))-(count(mixed(load('seq8192_2')))-count(mixed(load('seq8192_-2')))))/.006
                H=homogeneous(aa);err=H/BB-1
                result[mode]={'Bx_sequence':BB.tolist(),'Bx_homogeneous':H.tolist(),'relative_discrepancy':err.tolist(),'worst':float(max(abs(err)))}
                (D/'revalidation.json').write_text(json.dumps(result,indent=2))
                print(mode,result[mode],flush=True)
                if mode=='production_A1':assert max(abs(err))<2e-4
                else:assert abs(err[0])>8e-4
    result['table_relative']=(np.array(result['sequence']['4096']['Bx_h0005'])/np.array(result['sequence']['8192']['Bx_h0005'])-1).tolist()
    assert max(abs(np.array(result['table_relative'])))<2e-4
    result['status']='ADR-0011 Q4 UNIT-BOUNDARY PREREQUISITE VALIDATED'
    result['claim_boundary']='Candidate unit prerequisite only. Production Phase-5B B_i does not exist; PB7 itself is not passed.'
    (D/'revalidation.json').write_text(json.dumps(result,indent=2));print(result['status'])
