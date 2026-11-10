#!/usr/bin/env python3
"""Actual event coverage against a piecewise reference; no extrapolated times."""
import csv,json,pathlib,sys
from decimal import Decimal as D,localcontext
from piecewise import MODEL,Coupled,domain,END,YEAR,rows,lift
from run import save
out=pathlib.Path(sys.argv[1]);result={}
with localcontext() as ctx:
 ctx.prec=110
 for role in ('source','control'):
  hp=json.loads((out/f'hp-{role}-1000000.json').read_text());m=Coupled(MODEL,role);box=domain(m)
  pre=json.loads((out/'preflight.json').read_text());forcing=D(pre['roles'][role]['normalization_rhs_bound'])
  crossings=[]
  for r in hp['crossings']:
   k=r['knot'];left=D(r['left_s']);right=D(r['right_s']);u=D(r['U'])+forcing*right*(box['L']*right).exp()
   assert D(r['left_x'])+u<m.knots[k]<D(r['right_x'])-u,'crossing sign interval not certified'
   # Global source lower thermal velocity from positive direct power and maximum cooling.
   _,_,_,_,nu,photon,_=m.powers((box['hi'],*box['Y']))
   direct=m.pdir-m.erg*sum(box['Y'][j]*abs(m.sigma[j]) for j in range(2))
   cmax=max(m.capacity[k] for k in range(min(box['cells']),max(box['cells'])+2))
   speed=(direct-nu-photon)/(D('1e8')*box['hi'].exp()*cmax)-forcing
   assert speed>0
   crossings.append(dict(knot=k,time_bracket_s=[left,right],time_bracket_year=[left/YEAR,right/YEAR],primary_state_bound=u,thermal_velocity_lower_per_s=speed,knot_x=m.knots[k],T_K=D('1e8')*m.knots[k].exp()))
  native=[]
  for path in sorted((out/'native'/role).glob('tier-*/events.tsv')):
   for r in rows(path):
    if int(r['knot'])<0:continue
    h=next(c for c in crossings if c['knot']==int(r['knot']));t=lift(r['time']);lo,hi=h['time_bracket_s']
    native.append(dict(configuration=path.parent.name,event=r,reference_time_error_bound_s=max(abs(t-lo),abs(t-hi))))
  result[role]=dict(HP_crossings=crossings,native_events=native,status='CANDIDATE_ACTUAL_EVENT_COMPARISON',scope='Reference crossing sign and transversality bounds. Native event-time comparison includes integration error; it is not the isolated local event uncertainty. Full original event allocations remain separate.')
 save(out/'events-assessment.json',result)
