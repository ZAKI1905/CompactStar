#!/usr/bin/env python3
"""TEST-ONLY Track-R numerical-budget experiment; freshly produces every star.

No cached coefficients or review percentages are inputs. Dependencies: numpy,
scipy, mpmath (explicit failure if unavailable). Evidence persists outside source.
Independent momentum/common-potential mathematics supplies the reference route;
the old route integrates four intrinsic diagonal susceptibilities, with NO
neutral projector. Do not refactor old/new response ownership into one helper.
"""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import simpson, quad_vec, solve_ivp

L = np.array([[-1.,-1.],[1.,0.],[0.,1.]])
U = np.array([[1.,1.,1.],[0.,1.,0.],[0.,0.,1.]])

# PREDECLARED, before measurement: empirical enclosure = twice the sum of
# adjacent refinement differences plus independent route spread; no fitted
# coefficient tolerance. Such estimates are NUMERICALLY-BOUNDED, not rigorous.
# Correction acceptance: EACH named separation exceeds the sum of both route
# uncertainty estimates. A noncontracting quadrature ladder requires refinement
# or refusal unless both differences lie inside accumulated roundoff.
ORDERS = (8,16,32)
REFINEMENTS = ((4096,40000),(8192,40000),(8192,80000))
UROUND = np.finfo(float).eps/2


def norm(a):
    return float(np.linalg.norm(a,2))


def matrix(v):
    out = np.zeros(v.shape[:-1]+(3,3))
    out[...,0,0],out[...,1,1],out[...,1,2],out[...,2,2] = np.moveaxis(v,-1,0)
    out[...,2,1] = out[...,1,2]
    return out


def pack(a):
    return a[..., [0,1,1,2],[0,1,2,2]]


def zsolve(g):
    return L.T@np.linalg.solve(g,L)


def zerror(g, e):
    # Entrywise |Delta G|<=e; |invG| e bounds arbitrary signed perturbations.
    gi = np.linalg.solve(g,np.eye(3))
    b = np.abs(gi)@e
    rho = norm(b)
    if rho >= 1:
        raise AssertionError('uncertainty prevents reliable inversion')
    # Neumann majorant preserves tiny neutron-only perturbations through L.
    major = np.linalg.solve(np.eye(3)-b,b@np.abs(gi))
    if np.any(major<0):
        raise AssertionError('negative uncertainty majorant')
    return np.abs(L).T@major@np.abs(L),rho


def qerror(g,e):
    gx=U@g@U.T;ex=np.abs(U)@e@np.abs(U).T
    a=gx[0,0];h=gx[1:,0];ea=ex[0,0];eh=ex[1:,0]
    if not a>ea:
        raise AssertionError('uncertain baryon denominator')
    return (ex[1:,1:]+(np.outer(np.abs(h),eh)+np.outer(eh,np.abs(h))+np.outer(eh,eh))/(a-ea)
            +np.abs(np.outer(h,h))*ea/(a*(a-ea)))


class Model:
    def __init__(self, directory):
        lines=(directory/'model.txt').read_text().splitlines()
        self.mass=np.array([float(line.split()[0]) for line in lines[:4]])
        self.hc=float(lines[0].split()[1])
        self.onsets=np.array(list(map(float,lines[4].split())))
        self.geom=float(lines[5])
        self.identity=lines[6:]
        # High-precision independent thresholds from exported binary constants.
        with mp.workdps(80):
            mn,mp_,me,mm=map(mp.mpf,self.mass)
            hc=mp.mpf(self.hc)
            lam=((mn-mp_)*(mn+mp_)+me*me)/(2*mn)
            self.lambda0=float(lam)
            self.true_n0=float(((lam-me)*(lam+me))**mp.mpf('1.5')/(3*mp.pi**2*hc**3))
        self.p0=self.density(self.lambda0,self.mass[2])
        self.kp0=self.kinetic(self.p0,self.mass[1])

    def density(self,mu,mass):
        return np.maximum((mu-mass)*(mu+mass),0)**1.5/(3*np.pi**2*self.hc**3)

    def kinetic(self,n,mass):
        p=self.hc*np.cbrt(3*np.pi**2*n)
        return p*p/(np.hypot(mass,p)+mass)

    def state(self,n):
        n=np.asarray(n)
        low=np.full(n.shape,self.lambda0)
        high=np.full(n.shape,150.)
        def at(lam):
            e=self.density(lam,self.mass[2]); mu=self.density(lam,self.mass[3]); p=e+mu
            # Stable neutron kinetic energy relative to the analytic onset.
            dn=(lam-self.lambda0)+(self.kinetic(p,self.mass[1])-self.kp0)
            nn=np.maximum(dn*(2*self.mass[0]+dn),0)**1.5/(3*np.pi**2*self.hc**3)
            return np.stack([nn,p,e,mu],axis=-1)
        for _ in range(58):
            mid=(low+high)/2; ns=at(mid)
            below=ns[...,0]+ns[...,1]<n
            low=np.where(below,mid,low); high=np.where(below,high,mid)
        ns=at((low+high)/2)
        pe=n<=self.true_n0
        return np.where(pe[...,None],np.stack([n*0,n,n,n*0],axis=-1),ns)

    def susceptibilities(self,n):
        ns=self.state(n)
        p=self.hc*np.cbrt(3*np.pi**2*ns)
        return np.hypot(self.mass,p)*p/(np.pi**2*self.hc**3)

    def response(self,n):
        s=self.susceptibilities(n)
        den=s[...,1]+s[...,2]+s[...,3]
        safe=np.where(den>0,den,1.)
        return np.stack([s[...,0],s[...,2]*(s[...,1]+s[...,3])/safe,
                         -s[...,2]*s[...,3]/safe,s[...,3]*(s[...,1]+s[...,2])/safe],axis=-1)

    def pe_at_h(self,h):
        # Analytic inversion of sqrt(mp²+p²)+sqrt(me²+p²)=S.
        p_mass,e_mass=self.mass[1:3]
        excess=(p_mass+e_mass)*np.expm1(max(h,0.))
        s=p_mass+e_mass+excess
        p2=excess*(2*(p_mass+e_mass)+excess)*(s*s-(p_mass-e_mass)**2)/(4*s*s)
        p=np.sqrt(max(p2,0.))
        n=(p/self.hc)**3/(3*np.pi**2)
        # Momentum integral, independent of the production antiderivative/series.
        x,w=np.polynomial.legendre.leggauss(16)
        t=(x+1)/2
        eps=sum(3*n*np.dot(w/2,t*t*np.hypot(m,p*t)) for m in (p_mass,e_mass))
        pressure=sum(n*p*p*np.dot(w/2,t**4/np.hypot(m,p*t)) for m in (p_mass,e_mass))
        chi=np.hypot([p_mass,e_mass],p)*p/(np.pi**2*self.hc**3)
        c=0. if p==0 else chi.prod()/chi.sum()
        return n,eps*self.geom,pressure*self.geom,c

    def mp_state(self,n):
        """Independent 70-digit equilibrium and chi for onset measure bounds."""
        with mp.workdps(70):
            mn,p,e,mu=map(mp.mpf,self.mass);hc=mp.mpf(self.hc);n=mp.mpf(n)
            def density(u,ma):
                return max((u-ma)*(u+ma),mp.mpf(0))**mp.mpf('1.5')/(3*mp.pi**2*hc**3)
            def chemical(v,ma):
                return mp.sqrt(ma*ma+hc*hc*(3*mp.pi**2*v)**(mp.mpf(2)/3))
            zero=((mn-p)*(mn+p)+e*e)/(2*mn)
            if n<=density(zero,e):
                ns=[mp.mpf(0),n,n,mp.mpf(0)];bmu=chemical(n,p)+chemical(n,e)
            else:
                a=zero;b=mp.mpf(150)
                for _ in range(240):
                    lam=(a+b)/2;ne=density(lam,e);nm=density(lam,mu);np_=ne+nm
                    bnmu=chemical(np_,p)+lam;nn=density(bnmu,mn)
                    if nn+np_<n:a=lam
                    else:b=lam
                lam=(a+b)/2;ne=density(lam,e);nm=density(lam,mu);np_=ne+nm
                bmu=chemical(np_,p)+lam;ns=[density(bmu,mn),np_,ne,nm]
            chi=[chemical(v,ma)*hc*(3*mp.pi**2*v)**(mp.mpf(1)/3)/(mp.pi**2*hc**3)
                 if v else mp.mpf(0) for v,ma in zip(ns,[mn,p,e,mu])]
            return +mp.log(bmu),[+v for v in chi]

    def neutron_source_upper(self):
        # The returned Npe root has |F|<=5e-11 (source), in addition to 64 eps
        # subtraction guards. Bound its composition uncertainty by residual/Dmin.
        # Dp+De decreases with density; evaluate at 1.001*n0 for a conservative
        # lower derivative across this tiny neighbourhood, verified below.
        with mp.workdps(70):
            mn,p,e,_=map(mp.mpf,self.mass);hc=mp.mpf(self.hc)
            ulp=mp.mpf(float(np.spacing(self.onsets[0])))
            nn_guard=mp.mpf(2)**30*ulp
            probe=mp.mpf(self.onsets[0])*mp.mpf('1.001')
            pf=hc*(3*mp.pi**2*probe)**(mp.mpf(1)/3)
            dmin=sum(mp.pi**2*hc**3/(mp.sqrt(ma*ma+pf*pf)*pf) for ma in [p,e])
            residual=mp.mpf('5e-11')+128*mp.mpf(np.finfo(float).eps)*(2*mn+e)
            nn=nn_guard+residual/dmin+2*ulp
            bnmu=mp.sqrt(mn*mn+hc*hc*(3*mp.pi**2*nn)**(mp.mpf(2)/3))
            lam=((bnmu-p)*(bnmu+p)+e*e)/(2*bnmu)
            ne=((lam-e)*(lam+e))**mp.mpf('1.5')/(3*mp.pi**2*hc**3)
            assert nn+ne<probe
            return float(np.nextafter(float(nn+ne),np.inf)),float(nn_guard),float(residual/dmin)

    def muon_source_interval(self):
        with mp.workdps(70):
            mn,p,e,mu=map(mp.mpf,self.mass);hc=mp.mpf(self.hc)
            def density(u,ma):
                return ((u-ma)*(u+ma))**mp.mpf('1.5')/(3*mp.pi**2*hc**3)
            def chemical(v,ma):
                return mp.sqrt(ma*ma+hc*hc*(3*mp.pi**2*v)**(mp.mpf(2)/3))
            ne=density(mu,e);centre=chemical(ne,p)+mu
            probe=mp.mpf(self.onsets[1])*mp.mpf('1.001')
            # Encloses both source 64*eps guards and elementary subtraction noise.
            margin=128*mp.mpf(np.finfo(float).eps)*(2*chemical(probe,mn)+chemical(probe,e)+mu)
            return (float(np.nextafter(float(density(centre-margin,mn)+ne),-np.inf)),
                    float(np.nextafter(float(density(centre+margin,mn)+ne),np.inf)),float(margin))


def weight(r,m,nu):
    return 4*np.pi*r*r*np.exp(-nu)/np.sqrt(1-2*m/r)


def integrate(d, model, order, onset_aware=True, intrinsic=False, partition_level=0):
    r,m,nu,nb=d[:,:4].T
    ons=np.interp(model.onsets,nb[::-1],r[::-1])
    cuts=np.unique(np.r_[r,ons]) if onset_aware else r
    for _ in range(partition_level):
        cuts=np.sort(np.r_[cuts,(cuts[:-1]+cuts[1:])/2])
    x,w=np.polynomial.legendre.leggauss(order)
    t=(x+1)/2; w=w/2
    total=np.zeros(4)
    for start in range(0,len(cuts)-1,1024):
        a=cuts[start:min(start+1024,len(cuts)-1),None]
        b=cuts[start+1:min(start+1024,len(cuts)-1)+1,None]
        rr=a+(b-a)*t; ww=np.broadcast_to((b-a)*w,rr.shape).copy()
        if onset_aware:
            for onset in ons:
                left=(b[:,0]==onset); right=(a[:,0]==onset)
                rr[left]=b[left]-(b[left]-a[left])*t*t
                rr[right]=a[right]+(b[right]-a[right])*t*t
                ww[left]=(b[left]-a[left])*2*t*w
                ww[right]=(b[right]-a[right])*2*t*w
        mm=np.interp(rr,r,m); nv=np.interp(rr,r,nu); bn=np.interp(rr,r,nb)
        c=model.susceptibilities(bn) if intrinsic else model.response(bn)
        total+=np.sum((ww*weight(rr,mm,nv))[...,None]*c,axis=(0,1))
    c=model.susceptibilities(nb[:1])[0] if intrinsic else model.response(nb[:1])[0]
    # Central omission bracket supplied separately; this leading regular term is explicit.
    total+=4*np.pi/3*r[0]**3*np.exp(-nu[0])*c
    return total*1e54,len(cuts),order*(len(cuts)-1)


def old_z(intrinsic_integrals):
    # F2005 unprojected susceptibility/inversion; no corrected response helper.
    bn,bp,be,bmu=intrinsic_integrals
    shared=1/bn+1/bp
    return np.array([[shared+1/be,shared],[shared,shared+1/bmu]])


def tail(d,model,g):
    r,m,nu,nb=d[-1,:4]
    pf=model.hc*np.cbrt(3*np.pi**2*nb)
    h=np.log1p(sum(pf*pf/(np.hypot(ma,pf)+ma) for ma in model.mass[1:3])/sum(model.mass[1:3]))
    _,eps,pressure,cp=model.pe_at_h(h)
    assert h>0 and pressure>0 and nb<model.onsets[0]
    ru=2*m/(1-(1-2*m/r)*np.exp(2*h))
    volume=4*np.pi/3*(ru-r)*(ru*ru+ru*r+r*r)
    dm=volume*eps; mu=m+dm
    bound=volume*np.exp(-nu)*cp/np.sqrt(1-2*mu/r)*1e54
    def ode(q,y):
        rr,delta_mass,gee=y
        _,en,p,c=model.pe_at_h(q)
        mm=m+delta_mass
        dr=-rr*(rr-2*mm)/(mm+4*np.pi*rr**3*p)
        return [dr,4*np.pi*rr*rr*en*dr,
                weight(rr,mm,nu+h-q)*c*dr]  # units fm^-3/MeV * km³
    continuations=[]
    for tol in (1e-9,1e-11):
        sol=solve_ivp(ode,(h,0.),[r,0.,0.],method='DOP853',rtol=tol,
                      atol=[tol*1e-3,tol*1e-20,tol*1e-12],max_step=h/16)
        assert sol.success and sol.t[-1]==0
        continuations.append(sol.y[:,-1])
    rr,delta_mass,raw=continuations[-1]
    actual=raw*1e54
    direct_error=abs(continuations[0][2]-raw)*1e54*2
    # A second direct continuation integral: high-precision exact constant-M
    # Schwarzschild exterior r(h). Its neglected self-gravity is separately
    # bounded by the positive shell mass. No ODE or pe_at_h helper is called.
    with mp.workdps(60):
        rm,mm,hm,num=map(mp.mpf,[r,m,h,nu]);pm,em=map(mp.mpf,model.mass[1:3]);hc=mp.mpf(model.hc)
        def shell_integrand(q):
            s=(pm+em)*mp.exp(q)
            excess=(pm+em)*mp.expm1(q)
            p2=excess*(2*(pm+em)+excess)*(4*pm*em+2*(pm+em)*excess+excess*excess)/(4*s*s)
            pf=mp.sqrt(p2)
            chip=mp.sqrt(pm*pm+p2)*pf/(mp.pi**2*hc**3)
            chie=mp.sqrt(em*em+p2)*pf/(mp.pi**2*hc**3)
            cp=chip*chie/(chip+chie) if pf else 0
            radius=2*mm/(1-(1-2*mm/rm)*mp.exp(2*(hm-q)))
            dr=radius*(radius-2*mm)/mm
            return 4*mp.pi*radius**2*mp.exp(-num-hm+q)*cp/mp.sqrt(1-2*mm/radius)*dr
        shell_reference=float(mp.quad(shell_integrand,[0,hm])*mp.mpf(10)**54)
    self_gravity_allowance=bound*16*dm/(m*(1-2*mu/r))
    assert abs(shell_reference-actual)<=direct_error+self_gravity_allowance+4096*UROUND*bound
    # Floating rounding in mass comparisons is distinct from the positive mass bound.
    assert rr<=ru+64*UROUND*ru and 0<delta_mass<=dm*(1+1024*UROUND)
    assert 0<actual and actual+direct_error<bound
    cs=[model.pe_at_h(q)[3] for q in np.linspace(0,h,65)]
    assert np.all(np.diff(cs)>=0)
    lapse_error=dm/(r*(1-2*mu/r))
    e=np.zeros((3,3));e[1,1]=bound
    ez,_=zerror(g,e)
    return {'finite_cut_pressure_geom':pressure,'h_cut':h,'R_upper_km':ru,
        'M_upper_increment_km':dm,'direct_R_vacuum_km':rr,'direct_mass_increment_km':delta_mass,
        'direct_Gee_count_per_MeV':actual,'direct_numerical_error':direct_error,
        'independent_60_digit_exterior_shell_integral':shell_reference,
        'exterior_self_gravity_allowance':self_gravity_allowance,
        'bound_Gee_count_per_MeV':bound,'enclosed':True,'rank_one_axis':'ee',
        'old_intrinsic_tail_bounds':(bound*model.susceptibilities(np.array([nb]))[0]/cp).tolist(),
        'lapse_normalization_bound':lapse_error,'Z_relative_bound':norm(ez)/norm(zsolve(g)),
        'background_error_controlled_by_shell_mass':False},e,lapse_error


def refusal(d,model,g,directory):
    r,m,nu,nb=d[:,:4].T
    windows=np.loadtxt(directory/'windows.tsv',skiprows=1)
    results=[];total=np.zeros((3,3))
    for onset,side,last,first in windows:
        # Include rounding of the declared boundary, and inflate measured extent
        # before applying monotone model upper bounds to EVERY response component.
        width=2*max(abs(first-onset),abs(last-onset))+64*abs(np.spacing(onset))
        source=None
        if onset==model.onsets[0]:
            upper,nn_guard,root_error=model.neutron_source_upper()
            assert upper>max(last,first)
            width=upper-onset
            source={'nn_guard_fm3':nn_guard,'root_density_error_upper':root_error,
                    'source_derived_upper_nB':upper}
        else:
            lower,upper,margin=model.muon_source_interval()
            assert lower<min(last,first) and upper>max(last,first)
            width=(onset-lower) if side<0 else (upper-onset)
            source={'source_guard_potential_margin_MeV':margin,
                    'source_interval_lower_nB':lower,'source_interval_upper_nB':upper}
        lo,hi=sorted([onset,onset+side*width])
        radial=np.interp([hi,lo],nb[::-1],r[::-1])
        j=np.searchsorted(r,radial[0])-1
        assert j>=0 and radial[1]<=r[-1]
        rhi=r[min(j+2,len(r)-1)]
        mhi=np.max(m[max(j,0):min(j+3,len(m))])
        nulo=np.min(nu[max(j,0):min(j+3,len(nu))])
        wmax=4*np.pi*rhi*rhi*np.exp(-nulo)/np.sqrt(1-2*mhi/radial[0])
        # C <= diag(chi_n,chi_e,chi_mu) as PSD matrices. Entrywise offdiagonal
        # needs its own bound; |C_e,mu|<=min(chi_e,chi_mu).
        hlo,_=model.mp_state(lo);hhi,smp=model.mp_state(hi)
        with mp.workdps(70):
            dh=float(np.nextafter(float(hhi-hlo),np.inf))
        s=np.nextafter(np.array(list(map(float,smp))),np.inf)
        # Keep source-proven absent species exactly absent, not a positive floor.
        s[np.array([v==0 for v in smp])]=0.
        cmax=np.diag(s[[0,2,3]])
        cmax[1,2]=cmax[2,1]=min(s[2],s[3])
        # Independent differential TOV map: |dr/dh| <= r(r-2m_min)/m_min
        # for positive P. Use whole containing-cell extrema, not a point slope.
        mmin=float(np.min(m[max(j,0):min(j+3,len(m))]))
        physical_width=dh*rhi*(rhi-2*mmin)/mmin
        enclosed_width=max(float(radial[1]-radial[0]),physical_width)
        e=cmax*wmax*enclosed_width*1e54
        total+=e
        ez,rho=zerror(g,e)
        enn=np.zeros((3,3));enn[0,0]=e[0,0]
        eznn,_=zerror(g,enn)
        results.append({'onset':onset,'side':int(side),'provider_last_refused':last,
            'provider_first_available':first,'enclosing_density_width':width,
            'linear_profile_radial_width_km':float(radial[1]-radial[0]),
            'differential_TOV_radial_width_bound_km':physical_width,
            'source_guard_bound':source,'C_entrywise_upper':cmax.tolist(),
            'old_intrinsic_E_B':(s*wmax*enclosed_width*1e54).tolist(),
            'E_G':e.tolist(),'relative_G_nn_bound':e[0,0]/g[0,0],
            'E_Q':qerror(g,e).tolist(),
            'relative_Z_nn_only_bound':norm(eznn)/norm(zsolve(g)),
            'relative_Z_all_missing_response_bound':norm(ez)/norm(zsolve(g)), 'rho':rho,
            'coverage':'source guard and monotone model bounds; both radial maps enclosed'})
    return results,total


def main():
    exe=Path(sys.argv[1]).resolve();output=Path(sys.argv[2]).resolve();output.mkdir(parents=True,exist_ok=True)
    root=Path(tempfile.mkdtemp(prefix='run-',dir=output))
    (root/'predeclared.json').write_text(json.dumps({'orders':ORDERS,'refinements':REFINEMENTS,
        'estimator':'2 * sum adjacent differences + independent spread',
        'separation_acceptance':'each absolute separation > sum of both propagated uncertainties'},indent=2))
    datasets=[];gs=[];old_integrals=[];models=[]
    for table,radial in REFINEMENTS:
        directory=root/f't{table}-r{radial}'
        with (root/f'producer-{table}-{radial}.log').open('w') as log:
            subprocess.run([str(exe),str(directory),str(table),str(radial)],stdout=log,stderr=subprocess.STDOUT,check=True)
        d=np.loadtxt(directory/'profile.tsv',skiprows=1)
        model=Model(directory)
        datasets.append(d);models.append(model)
        gs.append(matrix(integrate(d,model,16)[0]))
        old_integrals.append(integrate(d,model,16,intrinsic=True)[0])
        print(f'fresh fixture {table}/{radial}: {len(d)} nodes',flush=True)
    d=datasets[-1];model=models[-1];r,m,nu,nb=d[:,:4].T
    ladder_results=[integrate(d,model,n) for n in ORDERS]
    ladder=[matrix(v[0]) for v in ladder_results]
    g=ladder[-1];z=zsolve(g)
    differences=[np.abs(ladder[i+1]-ladder[i]) for i in range(2)]
    bisected=matrix(integrate(d,model,16,partition_level=1)[0])
    partition_error=np.abs(bisected-g)
    roundoff=4096*UROUND*len(r)*np.abs(g)
    assert norm(differences[1])<=norm(differences[0])+norm(roundoff)
    quadrature=2*sum(differences)+2*partition_error+roundoff
    c=model.response(nb)
    weighted=weight(r,m,nu)[:,None]*c
    trap=matrix(np.trapezoid(weighted,x=r,axis=0)*1e54)
    simp=matrix(simpson(weighted,x=r,axis=0)*1e54)
    plain=matrix(integrate(d,model,16,onset_aware=False)[0])
    # Independent adaptive integration of a DECLARED DIFFERENT representation:
    # linear nodal response. Its interpolation spread must not masquerade as
    # quadrature error on the analytic-provider route.
    def interpolated(rr):
        cc=np.array([np.interp(rr,r,c[:,i]) for i in range(4)])
        return weight(rr,np.interp(rr,r,m),np.interp(rr,r,nu))*cc
    adaptive,adaptive_error=quad_vec(interpolated,r[0],r[-1],points=r[1:-1],
                                    epsabs=1e-16,epsrel=1e-9,limit=2*len(r))
    adaptive=matrix(adaptive*1e54)
    # Local provider independent comparison; only actual active H is solved.
    embeddings={1:np.array([[0.],[1.],[0.]]),2:np.array([[1.,-1.],[0.,1.],[0.,0.]]),
                3:np.array([[1.,-1.,-1.],[0.,1.,0.],[0.,0.,1.]])}
    own=np.zeros((len(r),3,3)); local_max=0.;condition_max=0.;residual_max=0.
    for dim,t in embeddings.items():
        mask=d[:,5]==dim;h=d[mask,12:21].reshape(-1,3,3)[:,:dim,:dim]
        inv=np.linalg.solve(h,np.broadcast_to(np.eye(dim),h.shape))
        own[mask]=t@inv@t.T
        condition_max=max(condition_max,float(np.linalg.cond(h).max()))
        residual_max=max(residual_max,float(np.max(np.abs(h@inv-np.eye(dim)))))
    available=d[:,5]>0
    if not np.all(available):
        raise AssertionError('nodal provider refusal requires interval adapter; no skipped node')
    delta=np.abs(own-matrix(c))
    local_max=float(np.max(np.linalg.norm(delta,axis=(1,2))/np.linalg.norm(matrix(c),axis=(1,2))))
    provider=np.trapezoid(weight(r,m,nu)[:,None,None]*delta,x=r,axis=0)*1e54
    # Match the actual profile composition, energy and pressure to the model.
    profile_ns=d[:,22:26]
    assert np.all(profile_ns>=0)
    pfs=model.hc*np.cbrt(3*np.pi**2*profile_ns)
    chis=np.hypot(model.mass,pfs)*pfs/(np.pi**2*model.hc**3)
    den=chis[:,1]+chis[:,2]+chis[:,3]
    profile_c=matrix(np.stack([chis[:,0],chis[:,2]*(chis[:,1]+chis[:,3])/den,
        -chis[:,2]*chis[:,3]/den,chis[:,3]*(chis[:,1]+chis[:,2])/den],axis=-1))
    anchor=2*np.trapezoid(weight(r,m,nu)[:,None,None]*np.abs(profile_c-matrix(c)),x=r,axis=0)*1e54
    anchor_report={'energy_max_relative_mismatch':float(np.max(np.abs(d[:,21]-d[:,10]*model.geom)/d[:,21])),
        'pressure_max_relative_mismatch':float(np.max(np.abs(d[:,4]-d[:,11]*model.geom)/d[:,4])),
        'species_max_absolute_density_mismatch':np.max(np.abs(profile_ns-d[:,6:10]),axis=0).tolist(),
        'species_max_density_scaled_mismatch':np.max(np.abs(profile_ns-d[:,6:10])/nb[:,None],axis=0).tolist(),
        'interpretation':'same model/bytes/constants; interpolation mismatch, not a universal tolerance'}
    # Grid spread is empirical background characterization, not shell-mass bound.
    background=2*(np.abs(gs[1]-gs[0])+np.abs(gs[2]-gs[1]))
    independent=np.loadtxt(directory/'independent-star.tsv',skiprows=1)
    independent_g=matrix(independent[:,2:6])
    independent_error=2*np.abs(independent_g[1]-independent_g[0])
    # Independent full-background spread includes its own vacuum tail. Retain it
    # conservatively without subtracting a guessed tail or cancelling errors.
    background+=np.abs(independent_g[-1]-g)+independent_error
    tail_report,tail_error,lapse_error=tail(d,model,g)
    refusal_report,refusal_error=refusal(d,model,g,directory)
    centre=np.abs(matrix(model.response(nb[:1])[0]))*4*np.pi/3*r[0]**3*np.exp(-nu[0])*1e54
    eg=quadrature+background+anchor+2*provider+tail_error+refusal_error+lapse_error*np.abs(g)+centre
    ez,rho=zerror(g,eg)
    intrinsic=integrate(d,model,32,intrinsic=True)[0]
    old=old_z(intrinsic)
    ei=2*(np.abs(old_integrals[1]-old_integrals[0])+np.abs(old_integrals[2]-old_integrals[1])+
          np.abs(intrinsic-old_integrals[-1]))+4096*UROUND*len(r)*intrinsic
    ei+=np.abs(independent[-1,6:10]-intrinsic)+2*np.abs(independent[1,6:10]-independent[0,6:10])
    # OLD uncertainty belongs to its FOUR intrinsic integrals, not to a normwise
    # reallocation of corrected G uncertainty. Bound each term in those axes.
    old_anchor=2*np.trapezoid(weight(r,m,nu)[:,None]*np.abs(chis-model.susceptibilities(nb)),x=r,axis=0)*1e54
    ei+=old_anchor+np.array(tail_report['old_intrinsic_tail_bounds'])
    ei+=sum(np.array(window['old_intrinsic_E_B']) for window in refusal_report)
    ei+=lapse_error*intrinsic
    assert np.all(intrinsic>ei)
    old_error=old_z(intrinsic-ei)-old_z(intrinsic)
    separation=np.abs(z-old);uncertainty=ez+old_error
    assert all(separation[i,j]>uncertainty[i,j] for i,j in [(0,0),(0,1),(1,1)])
    assert norm(z-old)>norm(uncertainty)
    assert np.linalg.eigvalsh(g)[0]>norm(eg)
    assert abs(g[0,1])<=eg[0,1] and abs(g[0,2])<=eg[0,2]
    qx=U@g@U.T
    q=qx[1:,1:]-np.outer(qx[1:,0],qx[0,1:])/qx[0,0]
    eq=qerror(g,eg)
    assert np.linalg.eigvalsh(q)[0]>norm(eq)
    assert norm(np.linalg.solve(q,np.eye(2))-z)<=4096*UROUND*norm(z)
    report={'classification':'TEST-ONLY numerical characterization, not production acceptance',
        'evidence_directory':str(root),'model_identity':model.identity,
        'runtime':{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'mpmath':mp.__version__},
        'files':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest()
                 for p in root.glob('*/profile.tsv')}|
                {str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest()
                 for p in root.glob('*/freegas.tsv')},
        'quadrature':{'owner':'onset-aware segment Gauss-Legendre with quadratic endpoint map',
            'orders':ORDERS,'G_relative_adjacent':[norm(v)/norm(g) for v in differences],
            'bisected_GL16_relative_G_difference':norm(partition_error)/norm(g),
            'candidate_Z_relative_spread':{name:norm(zsolve(v)-z)/norm(z) for name,v in
                [('A_raw_trapezoid',trap),('B_nonuniform_Simpson',simp),
                 ('C_segment_GL16',plain),('D_onset_GL16',ladder[1]),('E_adaptive_linear_response',adaptive)]},
            'adaptive_representation_quadrature_error':float(adaptive_error)*1e54,
            'node_count_order32':ladder_results[-1][2]},
        'G':g.tolist(),'Q':q.tolist(),'Z':z.tolist(),'old_Z':old.tolist(),
        'local':{'max_relative_provider_C_spread':local_max,'max_kappa_H_unscaled':condition_max,
                 'max_backward_absolute_residual':residual_max},
        'equilibrium_anchor':anchor_report,
        'independent_background':{'controls':[[1e-10,1e-8],[1e-12,1e-10]],
            'R_km':independent[:,0].tolist(),'M_km':independent[:,1].tolist(),
            'G':independent_g[-1].tolist(),
            'G_relative_refinement':norm(independent_error)/norm(g),
            'Z_relative_profile_difference':norm(zsolve(independent_g[-1])-z)/norm(z),
            'cross_background_Gee_difference':float(independent_g[-1,1,1]-g[1,1]),
            'cross_background_vs_trapezoid_Gee_difference':float(independent_g[-1,1,1]-trap[1,1]),
            'cross_background_difference_is_tail_measure':False},
        'background_relative_G_spread':[norm(gs[i+1]-gs[i])/norm(g) for i in range(2)],
        'tail':tail_report,'refusal_windows':refusal_report,
        'error_G_components':{name:value.tolist() for name,value in [('quadrature',quadrature),
            ('background_characterization',background),('provider_nodal_comparison',2*provider),
            ('equilibrium_anchor',anchor),('tail',tail_error),('refusal',refusal_error),('centre',centre)]},
        'E_G':eg.tolist(),'E_Q':eq.tolist(),'E_Z':ez.tolist(),'rho':rho,
        'old_intrinsic_error':ei.tolist(),
        'correction_gate':{'named_order':['Z_npe','Z_np','Z_npmu'],
            'signed_relative_to_old':[(z[i,j]-old[i,j])/old[i,j] for i,j in [(0,0),(0,1),(1,1)]],
            'frobenius_relative_to_old':float(np.linalg.norm(z-old)/np.linalg.norm(old)),
            'frobenius_relative_to_corrected':float(np.linalg.norm(z-old)/np.linalg.norm(z)),
            'frobenius_separation_over_uncertainty':float(np.linalg.norm(z-old)/np.linalg.norm(uncertainty)),
            'absolute_separations':[separation[i,j] for i,j in [(0,0),(0,1),(1,1)]],
            'combined_uncertainties':[uncertainty[i,j] for i,j in [(0,0),(0,1),(1,1)]],
            'minimum_separation_over_uncertainty':float(np.min(separation/uncertainty)),
            'PASS':True,'A18_substitute':False},
        'status':'PRECURSOR PASS; enclosure limitations must be assessed in planning record'}
    (root/'result.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(report,indent=2,allow_nan=False))


if __name__=='__main__':
    main()
