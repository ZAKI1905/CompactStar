#!/usr/bin/env python3
"""Comparison bounds for isolated LH and net-beta temperature contributions.

No new trajectory or production equation is evaluated. See CHEMICAL_RESPONSE.md.
Bounds use the actual coupled reference and a Lipschitz comparison inequality.
"""
import bisect,json,pathlib,sys
from decimal import Decimal as D,localcontext
from piecewise import MODEL,Coupled,domain,SCALE,END,PI
from run import save

def assess(out):
 with localcontext() as ctx:
  ctx.prec=110;m=Coupled(MODEL,'source');box=domain(m)
  hp=json.loads((out/'hp-source-1000000.json').read_text())['records']
  pre=json.loads((out/'preflight.json').read_text());forcing=D(pre['roles']['source']['normalization_rhs_bound'])
  conversion=m.erg*m.kb/m.kberg
  beta_coeff=[D(a)-conversion*D(b) for a,b in zip((22020,5670,420,9),(14680,7560,840,24))]
  umax=max(box['Y'])/(m.kb*D('1e8')*box['lo'].exp()*PI)
  assert beta_coeff[0]>sum(abs(beta_coeff[k])*umax**(2*k) for k in range(1,4)), 'net-beta sign not established'
  integrals={k:[D(0),D(0)] for k in ('LH','net_beta_loss')};weighted={k:D(0) for k in integrals}
  Jxx=max(c['J'][0][0] for c in box['cells'].values())
  cross=sum(max(c['J'][0][j] for c in box['cells'].values())*SCALE[j] for j in (1,2))
  # f0_x <= 0: C' > 0, positive net direct-minus-cooling power, cooling grows with x.
  assert all(m.slopes[c]>0 for c in box['cells'])
  _,_,_,_,nu_max,photon_max,_=m.powers((box['hi'],*box['Y']))
  pdir_min=m.pdir-m.erg*sum(box['Y'][j]*abs(m.sigma[j]) for j in range(2))
  assert pdir_min-nu_max-photon_max>0, 'thermal damping sign unproved'
  for left,right in zip(hp,hp[1:]):
   t=D(left['t_s']);h=D(right['t_s'])-t
   u=D(left['U'])+forcing*t*(box['L']*t).exp()
   reach=u+(box['M']+forcing)*h
   x=D(left['x']);xl=x-reach;xh=x+reach
   assert box['lo']<xl<xh<box['hi']
   tl=D('1e8')*xl.exp();th=D('1e8')*xh.exp()
   # Include all intersected capacity-cell endpoints and the interval endpoints.
   points=[xl,xh]+[k for k in m.knots if xl<k<xh]
   caps=[]
   for xx in points:
    c=bisect.bisect_right(m.knots,xx)-1
    caps.append(m.capacity[c]+m.slopes[c]*(xx-m.knots[c]))
   cmin=min(caps);cmax=max(caps)
   el=[max(D(0),abs(D(left[k]))-reach*SCALE[j+1]) for j,k in enumerate(('e','mu'))]
   eh=[abs(D(left[k]))+reach*SCALE[j+1] for j,k in enumerate(('e','mu'))]
   assert all(eh[j]+D('1e-10')*SCALE[j+1]<box['Y'][j] for j in range(2)), 'comparison chemical domain exhausted'
   assert box['lo']<xl-D('1e-10') and xh+D('1e-10')<box['hi'], 'comparison thermal domain exhausted'
   def heat(T,eta):return m.erg*sum(eta[j]*m.luminosity[j]/m.kberg*T**7*m.reaction(T,eta[j]) for j in range(2))
   hlo=heat(tl,el);hhi=heat(th,eh)
   net=[D(0),D(0)]
   for j in range(2):
    for k,b in enumerate(beta_coeff):
     a=m.luminosity[j]*b/(11513*(m.kb*PI)**(2*k+2))
     low=el[j]**(2*k+2)*tl**(6-2*k);high=eh[j]**(2*k+2)*th**(6-2*k)
     net[0]+=a*(low if a>=0 else high);net[1]+=a*(high if a>=0 else low)
   net[0]=max(D(0),net[0])
   for key,(lo,hi) in [('LH',(hlo,hhi)),('net_beta_loss',net)]:
    # Outward padding, far above 110-digit arithmetic/PI error.
    a=lo/(th*cmax)*(1-D('1e-75'));b=hi/(tl*cmin)*(1+D('1e-75'))+D('1e-90')
    integrals[key][0]+=h*a;integrals[key][1]+=h*b
    weighted[key]+=a*((-Jxx*(END-t-h)).exp()-(-Jxx*(END-t)).exp())/Jxx
  growth=(box['L']*END).exp();assert growth<2
  last=hp[-1];u=D(last['U'])+forcing*END*growth
  T=D('1e8')*D(last['x']).exp();Tlo=T*(-u).exp();Thi=T*u.exp()
  response={}
  for key,(lo,hi) in integrals.items():
   coupling=cross*END*growth*hi
   dmin=max((2-growth)*lo,weighted[key]-coupling);dmax=min(growth*hi,hi+coupling)
   assert dmax<D('1e-10'),'comparison tube not negligible; recheck domain'
   if key=='LH':bounds=[Tlo*(1-(-dmin).exp()),Thi*(1-(-dmax).exp())]
   else:bounds=[-Thi*(dmax.exp()-1),-Tlo*(dmin.exp()-1)]
   response[key]=dict(integral_forcing_bounds=integrals[key],delta_x_magnitude_bounds=[dmin,dmax],temperature_effect_K=bounds,relative_half_width=abs(bounds[1]-bounds[0])/abs(bounds[1]+bounds[0]),signal_gate_10_percent=abs(bounds[1]-bounds[0])<=D('.1')*abs(bounds[1]+bounds[0]),cross_chemical_feedback_bound=coupling)
  result=dict(status='CONDITIONAL_COMPARISON_BOUND_CANDIDATE',responses=response,Jxx=Jxx,thermal_chemical_cross_majorant=cross,L=box['L'],growth=growth,domain=box,source_reference_temperature_uncertainty_K=Thi-T,net_beta_sign='negative throughout the checked real domain',scope='Isolated mathematical response defined by omitting LH, or omitting LH-DeltaLnu, from the thermal RHS only, retaining coupled chemical kinetics and all other terms. No counterfactual trajectory was run. Not a new physical BNV model or a native Stage B qualification.',assumptions=['Fixed accepted stellar inputs and retained Decimal arithmetic pads.','Accepted interpolation is continuous; piecewise Jacobian majorants bound both original and term-suppressed vector fields.','The comparison tube stays in the declared domain; maximum extra displacement is explicitly checked below 1e-10 scaled units.'])
  save(out/'chemical-response.json',result)
if __name__=='__main__':assess(pathlib.Path(sys.argv[1]))
