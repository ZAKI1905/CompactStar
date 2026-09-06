#!/usr/bin/env python3
"""U10 independent count routes and propagated quadrature bounds.
Production integrates linearly interpolated nodal wV*nB. The Gauss route
integrates nB and m interpolants instead. Their difference is representation
error; monotone endpoint bounds contain both without fitted tolerances.
"""
from pathlib import Path
import sys,json
sys.dont_write_bytecode = True
import numpy as np
import revalidate_background as b
out=Path(sys.argv[3]);result={'status':'IN PROGRESS','stars':{}}
for folder,names in [(sys.argv[1],['g8192r20000','g8192r40000']),(sys.argv[2],['cmf-r10000','cmf-r20000','cmf-r40000','cmf-r80000'])]:
    meta={s.split()[0]:s.split() for s in (Path(folder)/'meta.tsv').read_text().splitlines()}
    for name in names:
        a=np.loadtxt(Path(folder)/(name+'.tsv'),skiprows=1);r,m,nb=a[:,0],a[:,1],a[:,4];B=float(meta[name][12]);F=4*np.pi*r*r*nb/np.sqrt(1-2*m/r)
        nodal=float(np.trapezoid(F,r)*1e54)
        z,w=b.quad(a);rr,mm,nn=z[:,:,0],z[:,:,1],z[:,:,4]
        gauss=float(np.sum(w*4*np.pi*rr*rr*nn/np.sqrt(1-2*mm/rr))*1e54)
        # Independent physical CGS integral, recovering the physical TOV mass.
        rc=a[:,18];cg=4*np.pi*rc*rc*nb*1e39/np.sqrt(1-2*b.G*a[:,19]/(rc*b.c*b.c))
        cgs=float(np.trapezoid(cg,rc));floor=16*len(a)*np.finfo(float).eps
        assert abs(cgs/B-1)<floor and abs(nodal/B-1)<floor,'U10 count/physical-volume inconsistency'
        assert np.all(np.diff(nb)<=0) and np.all(np.diff(m)>=0)
        r0,r1=r[:-1],r[1:];m0,m1=m[:-1],m[1:];n0,n1=nb[:-1],nb[1:]
        assert np.all(1-2*m1/r0>0)
        lo=float(np.sum(np.diff(r)*4*np.pi*r0*r0*n1/np.sqrt(1-2*m0/r1))*1e54)
        hi=float(np.sum(np.diff(r)*4*np.pi*r1*r1*n0/np.sqrt(1-2*m1/r0))*1e54)
        assert lo<=min(B,gauss,cgs)<=max(B,gauss,cgs)<=hi,'U10 independent quadrature enclosure failure'
        old=b.mixed(a);oldN=float(np.trapezoid(4*np.pi*r*r*nb/np.sqrt(1-2*old[:,1]/r),r)*1e54)
        center=4*np.pi*r[0]**3*nb[0]/3*1e54
        row=dict(B=B,nodal_B=nodal,physical_cgs_B=cgs,gauss_B=gauss,nodal_relative=nodal/B-1,physical_relative=cgs/B-1,roundoff_relative_bound=floor,gauss_relative=gauss/B-1,endpoint_bounds=[lo,hi],old_B_nodal=oldN,historical_shift=nodal/oldN-1,omitted_center_estimate=center,center_fraction=center/B,surface_r_km=r[-1],surface_nb=nb[-1],cutoff_policy='Canonical finite p_cut surface; no unmodeled vacuum extension; production starts at r0. Center omission quantified separately.')
        if name.startswith('g'):
            row['n_plus_p_max_relative']=float(max(abs((a[:,13]+a[:,14])/nb-1)))
            # Existing Track-R active-branch table interpolation preserves baryon
            # algebra to ~600 eps; use summation+division conditioning allowance.
            assert row['n_plus_p_max_relative']<1024*np.finfo(float).eps,'whole-baryon species algebra'
        else:row['whole_baryon_semantics']='Use EOS n_B, including bound baryons in crust; free n+p alone is not whole-baryon density.'
        result['stars'][name]=row;out.write_text(json.dumps(result,indent=2));print(name,'PASS',flush=True)
# Analytic incompressible metric and constant number density, independent integral.
M,R,n=1.5,12.,.1;k=2*M/R**3;t=np.sqrt(k)*R
exact=2*np.pi*n/k**1.5*(np.arcsin(t)-t*np.sqrt(1-t*t))
x,w=np.polynomial.legendre.leggauss(64);r=(x+1)*R/2
q=np.sum(w*R/2*4*np.pi*n*r*r/np.sqrt(1-k*r*r))
assert abs(q/exact-1)<64*np.finfo(float).eps,'analytic proper-volume count'
result['analytic_relative']=q/exact-1
result['status']='U10 PASS; whole-baryon proper volume independently validated'
out.write_text(json.dumps(result,indent=2));print(result['status'])
