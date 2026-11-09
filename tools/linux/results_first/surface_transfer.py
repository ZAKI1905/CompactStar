#!/usr/bin/env python3
"""Bounded pe-terminal transfer; scientific candidate, no global stellar solve.

The finite interior majorant is distinct from unsupported family endcaps.
All prior source/ledger records stay immutable. See surface-transfer/PREDECLARATION.md.
"""
import argparse
import copy
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import mpmath as mp

ROOT=Path(__file__).resolve().parents[3]
DOC=ROOT/'docs/validation/results-first-bnv-pilot/surface-transfer'
spec=importlib.util.spec_from_file_location('cut',ROOT/'tools/linux/results_first/cut_response.py')
cut=importlib.util.module_from_spec(spec);spec.loader.exec_module(cut)
iv=mp.iv
EPS=mp.mpf('2.1e-12')


def I(a,b=None):
    return iv.mpf([str(a),str(a if b is None else b)])


def low(x): return mp.mpf(x.a._mpi_[0])
def high(x): return mp.mpf(x.b._mpi_[1])
def mag(x): return max(abs(low(x)),abs(high(x)))
def hull(x): return [mp.nstr(low(x),60),mp.nstr(high(x),60)]


class Jet:
    """Three interval first derivatives, only for this terminal (H,R,M) map."""
    def __init__(self,v,d=None):
        self.v=v if hasattr(v,'_mpi_') else I(v)
        self.d=[I(0) for _ in range(3)] if d is None else d
    def __add__(self,b):
        b=b if isinstance(b,Jet) else Jet(b)
        return Jet(self.v+b.v,[a+c for a,c in zip(self.d,b.d)])
    __radd__=__add__
    def __neg__(self):return Jet(-self.v,[-d for d in self.d])
    def __sub__(self,b):return self+-b
    def __rsub__(self,b):return -self+b
    def __mul__(self,b):
        b=b if isinstance(b,Jet) else Jet(b)
        return Jet(self.v*b.v,[a*b.v+self.v*c for a,c in zip(self.d,b.d)])
    __rmul__=__mul__
    def __pow__(self,n):
        n=mp.mpf(n);return Jet(self.v**I(n),[I(n)*self.v**I(n-1)*d for d in self.d])
    def __truediv__(self,b):return self*(b**-1 if isinstance(b,Jet) else Jet(b)**-1)
    def __rtruediv__(self,b):return self**-1*b
    def exp(self):
        v=iv.exp(self.v);return Jet(v,[v*d for d in self.d])


def series(z,power,offset,n=32):
    """Integral binomial series with absolute geometric value/derivative tails."""
    c=mp.mpf(1);total=Jet(0)
    for j in range(n):
        total+=c*z**j/(2*j+offset)
        c*= (mp.mpf(power)-j)/(j+1)
    zm=high(z.v)
    if not 0 < zm < 1:
        raise ValueError('pe series domain')
    # Subsequent |binomial| ratios are <=1 for the +/-1/2 series.
    c=abs(c)/(2*n+offset)
    tail=c*zm**n/(1-zm)
    dtail=c*zm**(n-1)*(n/(1-zm)+zm/(1-zm)**2)
    return Jet(total.v+I(-tail,tail),[d+I(-dtail,dtail)*dz for d,dz in zip(total.d,z.d)])


def eos_jet(H):
    p,e=map(cut.wide,('938.27208816','0.51099895000'))
    hc=cut.wide(1/float('5.067730716156395e-3'))
    h0=cut.wide(float(p)+float(e))
    h=h0*H.exp();dh=h-h0
    k2=dh*(h+h0)*(h-p+e)*(h+p-e)/(4*h*h)
    k=k2**mp.mpf('.5');norm=iv.pi**2*I(hc)**3
    energy,pressure=Jet(0),Jet(0);chi=[];mus=[]
    for mass in (p,e):
        z=k2/(mass*mass);mu=(mass*mass+k2)**mp.mpf('.5');mus.append(mu)
        energy+=k**3*mass/norm*series(z,mp.mpf('.5'),3)
        pressure+=k**5/(3*norm*mass)*series(z,mp.mpf('-.5'),5)
        chi.append(mu*k/norm)
    kH=(h*h-h0*h0*(p-e)**2/(h*h))/(4*k)
    pH=sum(k**4/(3*norm*mu) for mu in mus)*kH
    return dict(p=pressure,e=energy,n=k**3/(3*norm),
        cv=sum(chi)*(iv.pi**2/3),chi=chi[0]*chi[1]/sum(chi),alpha=pressure/pH)


def structural(box,scales,B0):
    H,R,M=[Jet(v,[I(scales[j]) if i==j else I(0) for j in range(3)]) for i,v in enumerate(box)]
    q=eos_jet(H);geo=cut.wide(1.602176634e33*6.673e-8/float(2.99792458e10)**4*1e10)
    D=1-2*M/R;V=M/R**3+q['p']*(4*iv.pi*geo)
    A=D/V;W=R*A/(D**mp.mpf('.5'))*(4*iv.pi)
    r=-A/R*q['alpha'];m=R**2*q['e']*r*(4*iv.pi*geo)
    b=-W*q['alpha']*q['n']*mp.mpf('1e54')/B0
    k=(q['e']+q['p'])/V*q['alpha']*(4*iv.pi*geo)
    gl=m/M-2*r/R-k+q['alpha']
    return dict(H=q['alpha'],R=r,M=m,b=b,k=k,gl=gl,W=W,eos=q,D=D,V=V)


def state(raw):
    p=list(map(mp.mpf,raw['primary']))
    return [mp.mpf(raw['Hs']),p[1],p[2],mp.mpf(raw['K']),p[0],p[5],
            *map(mp.mpf,raw['G_packed']),mp.mpf(raw['N'][2])]


def rhs(y):
    H,R,M,K,B,C,*rest=y
    t=cut.terminal(H,R,M,K)
    g=rest[:4]
    return [t['alpha'],t['R'],t['M'],t['K'],t['B'],-t['K']*C+t['C_shell'],
            *[-t['K']*v+(t['G_shell'] if j==1 else 0) for j,v in enumerate(g)],t['B']]


def flow(y,delta,n):
    y=list(y);h=delta/n
    for _ in range(n):
        k1=rhs(y);k2=rhs([v+h*k/2 for v,k in zip(y,k1)])
        k3=rhs([v+h*k/2 for v,k in zip(y,k2)]);k4=rhs([v+h*k for v,k in zip(y,k3)])
        y=[v+h*(a+2*b+2*c+d)/6 for v,a,b,c,d in zip(y,k1,k2,k3,k4)]
    return y


def output(raw,y):
    updated=copy.deepcopy(raw);p=list(updated['primary']);H,R,M,K,B,C=y[:6]
    dk=K-mp.mpf(raw['K']);p[0:3]=[B,R,M];p[4]=mp.mpf(p[4])*mp.exp(-6*dk)
    p[5]=C;p[7]=mp.mpf(p[7])*mp.exp(dk)
    updated.update(primary=p,Hs=H,K=K,G_packed=y[6:10],
        N=[raw['N'][0],mp.mpf(raw['N'][1])+B-mp.mpf(raw['primary'][0]),y[10],raw['N'][3]],
        Lmu=mp.mpf(raw['Lmu'])*mp.exp(-6*dk),P0=mp.mpf(raw['P0'])*mp.exp(dk),
        P1=mp.mpf(raw['P1'])*mp.exp(dk))
    return cut.fields(updated)[0]


def point_models(inputs,comparison,m1,ledger,m4):
    """Original normalized bands and their original physical-unit conversions."""
    first=inputs['centers']['0']['2048'];fields=cut.fields(first)[0]
    L=mp.mpf('1e-6');span=abs(mp.mpf(m4['achieved_x'][-1]));sliver=L-span
    models={}
    for name,d in comparison.items():
        base='mu_B' if name=='mu_n_profile' else name
        rows=[r['static_completion'] for r in ledger['rows'] if r['quantity']==name]
        extras=sum(mp.mpf(t['pointwise_allowance']) for t in rows[0]['terms'] if
            'pointwise_allowance' in t and t['name'] not in ('unseen_reference_bias_increment',))
        e=10*max(abs(mp.mpf(v)) for row in d['signed_error'] for v in row)
        if name=='Ltilde_Me':e=mp.mpf('9.1e-11')
        if name.startswith('t_'):
            # Existing ratio-level point contributions, not count errors / a span.
            tail=max(mp.mpf(t['pointwise_allowance']) for t in rows[0]['terms'] if t['name']=='reference_stage_remainder')*2/15
        else:
            tail=max(mp.mpf(m1['reference'][str(t)][base]['integration']) for t in (0,10,20))
        params=max(sum(mp.mpf(m1['reference'][str(t)][base][k]) for k in ('eos_equilibrium','center','stencil')) for t in (0,10,20))
        factor=mp.mpf(1)
        if name=='mu_n_profile':
            factor=max(abs(mp.mpf(v)) for row in d['reference_values'] for v in row)/max(abs(mp.mpf(v)) for row in comparison['mu_B']['reference_values'] for v in row)
        band=e+extras+(tail+params)*factor+max(mp.mpf(v) for v in d['root_each'].values())
        K=mp.mpf(d['K']);nc=len(d['reference_values'][0]);comps=[]
        for j in range(nc):
            vals=[mp.mpf(row[j]) for row in d['reference_values']]
            deriv=(abs(vals[2]-vals[0])+2*band)/span+K*span/2+K*sliver
            bound=band+K*L**2/8+deriv*sliver
            if name.startswith('t_'):scale=I(1)
            else:
                key=name+'_'+str(j) if name.startswith('Z_row_') else 'mu_B' if name=='mu_n_profile' else name
                raw=fields[key]
                # Enclose binary64 normalized serialization and long-double raw
                # serialization. These are conversion roundoff, not new physics.
                rr=mp.mpf(2)**-63;rn=mp.mpf(2)**-52
                scale=I(raw*(1-rr),raw*(1+rr))/I(vals[0]*(1-rn),vals[0]*(1+rn))
            comps.append(dict(scale=scale,lo=min(vals)-bound,hi=max(vals)+bound,derivative=deriv))
        models[name]=dict(components=comps,band=band,K=K)
    return models


def physical(model):
    c=model['components'][0]
    return I(c['lo'],c['hi'])*c['scale']


def pressure_root(p,guess):
    # Monotone EOS root bracket; interval endpoint pressures verify containment.
    lo,hi=guess/2,guess*2
    for _ in range(140):
        mid=(lo+hi)/2;q=eos_jet(Jet(mid))['p'].v
        if high(q)<p:lo=mid
        elif low(q)>p:hi=mid
        else:break
    return I(lo,hi)


def bounds(inputs,models,source_response):
    B0=mp.mpf(inputs['cuts'][0]['B']);a=mp.mpf(source_response['slope']);intercept=mp.mpf(source_response['law_intercept_geom'])
    geo=cut.wide(1.602176634e33*6.673e-8/float(2.99792458e10)**4*1e10)
    guess=mp.mpf(inputs['centers']['0']['2048']['Hs'])
    Hlo=pressure_root(mp.exp(intercept-a*mp.mpf('1e-6'))/geo,guess)
    Hhi=pressure_root(mp.exp(intercept)/geo,guess)
    H=I(low(Hlo),high(Hhi));R=physical(models['radius']);g=physical(models['surface_gravity'])
    y=g*R*I('1e19')/I(mp.mpf('2.99792458e10')**2)
    u=y*(iv.sqrt(1+y*y)-y);M=R*u
    # Include the actual retained offsets for map corroboration, without claiming
    # that those finite samples establish endcap curvature between samples.
    rawstates=[]
    for t in ('0','10','20'):
        rawstates+=list(inputs['centers'][t].values())
        for ends in inputs['widths'][t].values():rawstates+=list(ends.values())
        for ends in inputs['grids'][t].values():rawstates+=list(ends.values())
    initial=[H,R,M]
    for j in range(3):
        samples=[state(r)[j] for r in rawstates]
        initial[j]=I(min(low(initial[j]),*samples),max(high(initial[j]),*samples))
    scales=[high(v) for v in initial]
    # Proposed bootstrap tube. The proof below must verify strict Picard
    # containment; the padding itself is never used as the reported allowance.
    tube=[v+I(-2*EPS*s,2*EPS*s) for v,s in zip(initial,scales)]
    f=structural(tube,scales,B0)
    if not (low(f['D'].v)>0 and low(f['V'].v)>0 and low(f['eos']['p'].v)>0):raise ValueError('invalid terminal domain')
    hn=min(mp.mpf(r['Hn']) for r in rawstates)
    if not high(tube[0])<hn:raise ValueError('terminal topology')
    A=mag(f['b'].v);tau=EPS/(1-abs(a)*A);eta=A*tau
    if not abs(a)*A<1:raise ValueError('unseparated pressure remap')
    ff=[f[n] for n in ('H','R','M')]
    rates=[mag(v.v)/s for v,s in zip(ff,scales)]
    if not all(tau*rate<2*EPS for rate in rates):raise ValueError('Picard containment failed')
    J=[[mag(d)/scales[i] for d in v.d] for i,v in enumerate(ff)]
    L=max(sum(row) for row in J)
    if not tau*L<1:raise ValueError('Picard contraction failed')
    euler=[tau**2/2*sum(J[i][j]*rates[j] for j in range(3))*scales[i] for i in range(3)]
    # Initial on-law family derivative bounds, established only on [-1e-6,0].
    Rx=max(c['derivative']*high(c['scale']) for c in models['radius']['components'])
    gx=max(c['derivative']*high(c['scale']) for c in models['surface_gravity']['components'])
    rj=Jet(R,[I(1),I(0),I(0)]);gj=Jet(g,[I(0),I(1),I(0)])
    yy=gj*rj*mp.mpf('1e19')/mp.mpf('2.99792458e10')**2
    mass=rj*yy*((1+yy*yy)**mp.mpf('.5')-yy)
    Mx=mag(mass.d[0])*Rx+mag(mass.d[1])*gx
    v0=max(abs(a)*mag(f['H'].v)/scales[0],Rx/scales[1],Mx/scales[2])
    wx=tau*sum(mag(d) for d in f['b'].d)*mp.exp(L*tau)*v0
    remap_den=1-abs(a)*A-wx
    if not remap_den>0:raise ValueError('unseparated fixed-B map derivative')
    K=iv.log(f['D'].v)/2+tube[0];kr=mag(f['k'].v);gr=mag(f['gl'].v);alpha=mag(f['H'].v)
    cvfactor=mp.mpf('1e15')*cut.wide('8.61733326214518e-11')**2*cut.wide('1.602176634e33')
    shell=abs(f['W'].v*f['H'].v)*iv.exp(tube[0]-K)
    Cshell=mag(shell*f['eos']['cv'].v)*cvfactor
    Gshell=mag(shell*f['eos']['chi'].v)*mp.mpf('1e54')
    z00=high(physical(models['Z_row_npe']))
    def zmax(n,j):
        c=models[n]['components'][j];return high(I(c['lo'],c['hi'])*c['scale'])
    Z=[[z00,zmax('Z_row_npe',1)],[zmax('Z_row_npmu',0),zmax('Z_row_npmu',1)]]
    rank=tau*Gshell*mp.exp(kr*tau);chem_den=1-rank*z00
    if not chem_den>0:raise ValueError('unseparated chemical inverse')
    allowances={}
    for name,model in models.items():
        if name.startswith('t_'):
            species={'t_n':'N_n','t_e':'N_e','t_mu':'N_mu'}[name];c=1 if name=='t_e' else 0
            T=model['components'][0];Tmax=max(abs(c-T['lo']),abs(c-T['hi']))
            Nx=models[species];Tx=max(high(v['scale']) for v in Nx['components'])*Nx['K']/B0
            U=Tx*eta+Tmax*wx/remap_den
            allowances[name]=dict(interior_point_bound=U,terminal_part=Tmax*wx/remap_den,remap_part=Tx*eta)
            continue
        terminal=[];remap=[]
        for j,c in enumerate(model['components']):
            qmax=max(abs(c['lo']),abs(c['hi']));scale=low(c['scale'])
            if name in ('N_n','N_mu'):u=mp.mpf(0) # exact shell support identity
            elif name=='N_e':u=B0*eta/scale
            elif name=='radius':u=tau*mag(f['R'].v)/scale
            elif name.startswith('Ltilde_'):u=qmax*mp.expm1(6*kr*tau)
            elif name in ('mu_B','P0_average','P1_average','mu_n_profile'):u=qmax*mp.expm1(kr*tau)
            elif name=='Cstar':u=qmax*mp.expm1(kr*tau)+tau*Cshell*mp.exp(kr*tau)/scale
            elif name=='surface_gravity':u=qmax*mp.expm1(gr*tau)
            elif name=='envelope':u=qmax*mp.expm1((gr/4+mp.mpf('.605')*kr)*tau)
            elif name=='Tsurface_inf':u=qmax*mp.expm1((gr/4+mp.mpf('.395')*kr+alpha)*tau)
            elif name.startswith('Z_row_'):
                i=0 if name=='Z_row_npe' else 1
                u=qmax*mp.expm1(kr*tau)+mp.exp(kr*tau)*rank/chem_den*Z[i][0]*Z[0][j]/scale
            else:raise ValueError(name)
            terminal.append(u);remap.append(c['derivative']*eta)
        allowances[name]=dict(interior_point_bound=max(u+v for u,v in zip(terminal,remap)),
                              terminal_part=max(terminal),remap_part=max(remap))
    return dict(allowances=allowances,initial_box=[hull(v) for v in initial],tube=[hull(v) for v in tube],
        scales=scales,tau=tau,eta=eta,pressure_denominator=1-abs(a)*A,fixed_B_derivative_denominator=remap_den,
        chemical_denominator=chem_den,horizon_denominator=low(f['D'].v),pressure_lower=low(f['eos']['p'].v),
        topology_margin=hn-high(tube[0]),Jacobian_majorant=J,Jacobian_norm=L,
        euler_remainder=euler,variation_excess=mp.expm1(L*tau),baryon_flow_bound=A,
        family_derivative_norm=v0,baryon_flow_x_bound=wx,
        supported_final_x_interval=[-mp.mpf('1e-6')+eta,-eta],endpoint_extension_supported=False,
        missing='interior-family continuation on both remapping endcaps; local terminal flow cannot establish it')


def numerical_witness(inputs,bound,digits):
    results={}
    with mp.workdps(digits):
        for t in ('0','10','20'):
            cases={'center':inputs['centers'][t]['2048'],**inputs['grids'][t]['2048']}
            for label,raw in cases.items():
                y=state(raw);base=cut.fields(raw)[0]
                for sign in (-1,1):
                    delta=sign*mp.mpf('2.1e-12');ys=[flow(y,delta,n) for n in (1,2,4)]
                    qs=[output(raw,v) for v in ys];checks={}
                    for name in base:
                        if name=='B':continue
                        checks[name]=cut.diagnostic([q[name] for q in qs],4)
                    rate=rhs(y);residual=[ys[-1][i]-y[i]-delta*rate[i] for i in range(3)]
                    if any(abs(v)>u for v,u in zip(residual,bound['euler_remainder'])):
                        raise ValueError('finite-flow escaped analytic remainder')
                    results[f'{t}:{label}:{sign}']=dict(final=ys[-1],delta={n:qs[-1][n]-base[n] for n in base},
                        structure_linearization_residual=residual,orders=checks,
                        B_Ne_identity=(ys[-1][4]-y[4])-(ys[-1][10]-y[10]))
    return results


def variational_witness(inputs,bound,response):
    """Propagate the local map's first variation, never re-shoot the interior."""
    B0=mp.mpf(inputs['cuts'][0]['B']);a=mp.mpf(response['slope']);out={}
    for t in ('0','10','20'):
        raw=inputs['centers'][t]['2048'];y=state(raw)
        rec=response['records'][t]['n2048'];Hder=rec['H_derivative']
        BH=mp.mpf(Hder['B']);Bs=mp.mpf(rec['terminal']['B'])
        hc_x=(B0-a*Bs)/BH
        v=[mp.mpf(rec['terminal']['alpha'])*a,
           mp.mpf(Hder['radius'])*hc_x+a*mp.mpf(rec['terminal']['R']),
           mp.mpf(Hder['mass'])*hc_x+a*mp.mpf(rec['terminal']['M']),mp.mpf(0)]
        def f(z):
            H,R,M,w=z[:4]
            def vector(H,R,M):
                q=cut.terminal(H,R,M,mp.mpf(raw['K']))
                return [q['alpha'],q['R'],q['M'],q['B']/B0]
            value=vector(H,R,M)
            jac=[[mp.diff(lambda h,r,m:vector(h,r,m)[i],(H,R,M),tuple(1 if k==j else 0 for k in range(3)))
                  for j in range(3)] for i in range(4)]
            return value+[sum(jac[i][j]*z[4+j] for j in range(3)) for i in range(4)]
        for sign in (-1,1):
            z=y[:3]+[mp.mpf(0)]+v
            step=sign*EPS/4
            for _ in range(4):
                k1=f(z);k2=f([q+step*k/2 for q,k in zip(z,k1)])
                k3=f([q+step*k/2 for q,k in zip(z,k2)]);k4=f([q+step*k for q,k in zip(z,k3)])
                z=[q+step*(aa+2*b+2*c+d)/6 for q,aa,b,c,d in zip(z,k1,k2,k3,k4)]
            if abs(z[3])>bound['eta'] or abs(z[7])>bound['baryon_flow_x_bound']:
                raise ValueError('variational map escaped majorant')
            out[f'{t}:{sign}']=dict(baryon_change_over_B0=z[3],baryon_flow_x=z[7],
                mapped_coordinate_first_order=-z[3],central_adjustment_first_order=-Bs/BH*sign*EPS)
    return out


def successor(ledger,bound,comparison,m4):
    """A complete disposition matrix, with unsupported endcaps kept null."""
    result=copy.deepcopy(ledger);scenario_max=mp.mpf(0);impact_max=mp.mpf(0)
    for row in result['rows']:
        name,i=row['quantity'],row['index'];old=row['static_completion']
        if name not in bound['allowances']:
            row['surface_transfer']=dict(status=old['status'],full_utilization=old.get('full_utilization'))
            continue
        U=bound['allowances'][name]['interior_point_bound'];f=mp.mpf(comparison[name]['fit_operator_norm'])
        extra=(f+(0 if i==0 else 2))*U
        endpoint=mp.mpf(0)
        if i==20:
            h=abs(mp.mpf(m4['achieved_x'][20])-mp.mpf(m4['achieved_x'][10]))
            sliver=max(mp.mpf(0),mp.mpf('1e-6')+mp.mpf(m4['achieved_x'][20]))
            endpoint=2*U/h*sliver
        N=mp.mpf(old['known_numerator']);b=mp.mpf(old['known_baseline'])
        absolute=name.startswith('t_');den=1-b-(0 if absolute else U)
        if den<=0:raise ValueError('baseline scenario not separated')
        scenario=(N+extra+endpoint)/den
        threshold=mp.mpf(row['historical']['threshold']);util=scenario/threshold
        delta=(scenario-N/(1-b))/threshold
        scenario_max=max(scenario_max,util);impact_max=max(impact_max,delta)
        row['surface_transfer']=dict(status='UNQUALIFIED_ENDPOINT_EXTENSION',full_budget=None,full_utilization=None,
            interior_point_bound=U,full_interval_point_bound=None,baseline_cut_bound=None if not absolute else 'NOT_APPLICABLE',
            missing=['fixed-B remapping family extension at x=0 and x=-1e-6'],
            endpoint_extension_scenario=dict(authorized=False,qualification=False,
                additional_numerator=extra,additional_endpoint=endpoint,additional_baseline=0 if absolute else U,
                full_budget=scenario,utilization=util,utilization_increase=delta,
                assumption='same existing reference bands and C-CONT slope/curvature/topology bounds on remapping endcaps'))
    result['surface_transfer_summary']=dict(qualified=False,original_ledger_preserved=True,
        full_C_CUT_budget_effect=None,endpoint_extension_scenario_max_utilization=scenario_max,
        endpoint_extension_scenario_max_utilization_increase=impact_max,
        historical_CREF_failures=42,expectation_violations=353)
    return result


def serial(v):
    if isinstance(v,mp.mpf):return mp.nstr(v,65)
    raise TypeError(type(v).__name__)


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('EKU Slurm required')
    mp.mp.dps=90;iv.dps=90
    global EPS
    EPS=mp.mpf('2.1e-12')
    manifest=cut.read(DOC/'inputs.json')
    for name,sha in manifest.items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=sha:raise ValueError('changed input '+name)
    inputs=cut.read(cut.DOC/'retained-inputs.json.gz')
    response=cut.read(cut.DOC/'evidence/response.json.gz')
    comp=cut.read(ROOT/'docs/validation/ba13-evidence/evidence/comparison.json.gz')
    m1=cut.read(ROOT/'docs/validation/ba13-closure/evidence/m1.json')
    m4=cut.read(ROOT/'docs/validation/ba13-closure/evidence/m4.json')
    ledger=cut.read(ROOT/'docs/validation/results-first-bnv-pilot/static-completion/evidence/ledger.json.gz')
    models=point_models(inputs,comp,m1,ledger,m4)
    b=bounds(inputs,models,response)
    witnesses=[numerical_witness(inputs,b,d) for d in (90,120)]
    variation=variational_witness(inputs,b,response)
    precision=max(abs(witnesses[0][key]['delta'][name]-witnesses[1][key]['delta'][name])/
                  max(abs(witnesses[1][key]['delta'][name]),mp.mpf('1e-100'))
                  for key in witnesses[0] for name in witnesses[0][key]['delta'])
    if precision>mp.mpf('1e-45'):raise ValueError('local map precision check failed')
    new=successor(ledger,b,comp,m4)
    for old,row in zip(ledger['rows'],new['rows']):
        assert {k:v for k,v in row.items() if k!='surface_transfer'}==old
    args.output.mkdir()
    def write(name,value):
        data=(json.dumps(value,default=serial,indent=2)+'\n').encode()
        (args.output/name).write_bytes(gzip.compress(data,mtime=0) if name.endswith('.gz') else data)
    write('bounds.json',b);write('local-maps.json.gz',witnesses[1]);write('variational-map.json',variation);write('ledger.json.gz',new)
    summary=dict(status='FINITE_INTERIOR_TRANSFER_SUPPORTED_ENDPOINT_EXTENSION_UNSUPPORTED',
        cut_closed=False,precision_relative_disagreement=precision,
        **new['surface_transfer_summary'],eta=b['eta'],tau=b['tau'],allowances=b['allowances'],
        stellar_global_integrations=0,local_maps=len(witnesses[1]),trajectories=0)
    write('summary.json',summary)
    with (args.output/'matrix.tsv').open('w') as f:
        f.write('index\tquantity\tstatus\tfull_utilization\tinterior_cut_bound\tendcap_scenario_utilization\n')
        for row in new['rows']:
            s=row['surface_transfer'];scenario=s.get('endpoint_extension_scenario',{})
            f.write('\t'.join(map(str,(row['index'],row['quantity'],s['status'],s['full_utilization'],
                mp.nstr(s['interior_point_bound'],18) if 'interior_point_bound' in s else 'NA',
                mp.nstr(scenario['utilization'],18) if scenario else 'NA')))+'\n')
    print(json.dumps(summary,default=serial))


if __name__=='__main__':main()
