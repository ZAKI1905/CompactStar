#!/usr/bin/env python3
"""Conditional continuous thermal enclosure using analytic RK4/Hermite remainders.

See preliminary-pair/ENCLOSURE.md. This adds a bound calculation around existing
offline RK4; it does not change any production numerical method or threshold.
"""
import collections,concurrent.futures,json,pathlib,sys,time
from decimal import Decimal as D,localcontext
from high_precision import Coupled,step,rows,lift,save
from roundoff import jacobian_bound
from verify_chemical_reference import PI,corner_uncertainty

def constants(model):
 R=D('.1');scale=[D(1),D('1e-6'),D('1e-6')]
 X=D('.0002')+R;Y=D('8e-10')+R*scale[1]
 Tlo=D('1e8')*(-X).exp();Thi=D('1e8')*X.exp()
 import bisect
 cell=bisect.bisect_right(model.knots,D(0))-1
 C0=model.capacity[cell]-model.slopes[cell]*model.knots[cell]
 Clo=C0-abs(model.slopes[cell])*X;assert Clo>0
 rates=[];nu=D(0);delta_nu=D(0)
 for j in range(2):
  a=[model.luminosity[j]*v/(model.kberg*11513*model.kb**(2*k+1)*PI**(2*k+2)) for k,v in enumerate((14680,7560,840,24))]
  rates.append(sum(a[k]*Y**(2*k+1)*Thi**(6-2*k) for k in range(4)))
  nu+=model.luminosity[j]*Thi**8
  extra=sum(model.luminosity[j]*v/(11513*(model.kb*PI)**(2*k+2))*Y**(2*k+2)*Thi**(6-2*k) for k,v in enumerate((22020,5670,420,9)))
  delta_nu+=extra;nu+=extra
 direct=abs(model.pdir)+model.erg*Y*sum(map(abs,model.sigma))
 heat=model.erg*Y*sum(rates);photon=model.photon*(lift(2.42)*X).exp()
 M=max([(direct+heat+nu+photon)/(Tlo*Clo)]+[(abs(model.source[i])+sum(abs(model.z[i][j])*rates[j] for j in range(2)))/scale[i+1] for i in range(2)])
 M=M*(1+D('1e-75'))+D('1e-90')
 radius=R/(2*M)*(1-D('1e-75'))
 J,box=jacobian_bound(model)
 L=max(sum(J[i][j]*scale[j]/scale[i] for j in range(3)) for i in range(3))
 powers=dict(P_dir_eq=D(0),P_dir_actual=direct,LH=heat,DeltaLnu=delta_nu,DeltaPbeta=heat+delta_nu,Lnu_eq=sum(model.luminosity)*Thi**8,Lnu_full=nu,Lgamma=photon,Lother=D(0),Lout=D(0),Pnet=direct+heat+nu+photon)
 return dict(R=R,M=M,analytic_time_radius_s=radius,L=L,scale=scale,J=J,real_box=box,cell=cell,complex_C_lower=Clo,complex_outgoing_power_bound=nu+photon,complex_gross_power_bound=direct+heat+nu+photon,complex_channel_bounds=powers)

def point(args):
 out,index=args;start=time.monotonic()
 with localcontext() as ctx:
  ctx.prec=110;model=Coupled(out,'source');c=constants(model)
  directory=next(p for p in (out/'sample'/f'point-{index}',out/'full'/f'point-{index}') if p.exists())
  cp=json.loads((directory/'packet/checkpoint.json').read_text())
  if cp['method']=='main-endpoint':return dict(point=index,status='EXACT_ENDPOINT')
  witnesses=rows(directory/'witnesses.tsv');chosen=next(r for r in witnesses if r['method']=='0' and r['subdivisions']=='8')
  trace=[r for r in rows(directory/'local-rhs.tsv') if r['context']==chosen['context']]
  lookup=collections.defaultdict(set)
  for r in trace:lookup[(lift(r['t_s']),*(lift(r[k]) for k in ('x','e','mu')))].add(lift(r['dx']))
  left=lift(cp['main_left_time_s']);y=tuple(lift(v) for v in cp['main_left_state']);nodes=[(left,y)]
  for r in rows(directory/'witness-steps.tsv'):
   if r['method']=='0' and r['subdivisions']=='8':nodes.append((lift(r['right']),tuple(lift(r[k]) for k in ('x','e','mu'))))
  errors=[];slopes=[];u=D(0);bounds=[];last=left;maximum_hp_remainder=D(0)
  for t,native in nodes:
   if t>last:
    h=(t-last)/64;q=h/c['analytic_time_radius_s'];assert q<D('.01')
    trunc=c['R']*q**5/(1-q)
    for _ in range(64):
     reach=u+c['M']*h
     assert c['real_box']['x_min']<y[0]-reach and y[0]+reach<c['real_box']['x_max']
     assert all(abs(y[i+1])+reach*c['scale'][i+1]<c['real_box']['eta_abs'][i] for i in range(2))
     y=step(model.coupled_rhs,y,h)
     u=(c['L']*h).exp()*u+trunc+D('1e-80')
    maximum_hp_remainder=max(maximum_hp_remainder,u)
   assert c['real_box']['x_min']<y[0]-u and y[0]+u<c['real_box']['x_max']
   assert all(abs(y[i+1])+u*c['scale'][i+1]<c['real_box']['eta_abs'][i] for i in range(2))
   native_dx=lookup[(t,*native)];assert len(native_dx)==1
   errors.append(abs(native[0]-y[0])+u)
   slopes.append(abs(next(iter(native_dx))-model.coupled_rhs(y)[0])+sum(c['J'][0][j]*u*c['scale'][j] for j in range(3)))
   last=t
  fourth=24*c['R']/c['analytic_time_radius_s']**4
  for i,(a,b) in enumerate(zip(nodes,nodes[1:])):
   h=b[0]-a[0]
   bounds.append(max(errors[i:i+2])+h/4*max(slopes[i:i+2])+fourth*h**4/384)
  sr={r['channel']:r for r in rows(directory/'packet/state.tsv')};dr={r['channel']:r for r in rows(directory/'packet/diagnostics.tsv')}
  ratios=[(abs(lift(cp['O2_state'][i])-y[i])+u*c['scale'][i])/lift(sr[k]['propagated_state_radius']) for i,k in enumerate(('x','eta_e','eta_mu'))]
  er=model.energy(y[1:]);ue,_=corner_uncertainty(model.energy,y[1:],tuple(u*c['scale'][i] for i in (1,2)))
  energy_ratio=(abs(lift(dr['Echem']['O2'])-er)+ue)/lift(dr['Echem']['U'])
  uniform=max(bounds)
  result=dict(point=index,status='CONDITIONAL_PASS' if uniform<=D('1e-12') and max(ratios+[energy_ratio])<=1 else 'FAIL',uniform_thermal_enclosure=uniform,thermal_band_limit=D('1e-12'),coupled_state_ratios=ratios,chemical_energy_ratio=energy_ratio,reference_final=y,reference_scaled_error_bound=u,maximum_RK4_scaled_remainder=maximum_hp_remainder,maximum_node_x_error=max(errors),maximum_node_dx_error=max(slopes),constants=c,elapsed_s=time.monotonic()-start,scope='Conditional on the exported frozen-anchor RHS; analytic Taylor tails, not measured agreement, bound RK4 and continuous Hermite interpolation.')
  save(out/'thermal-enclosure'/f'point-{index}.json',result)
  return {k:result[k] for k in ('point','status','uniform_thermal_enclosure','coupled_state_ratios','chemical_energy_ratio','reference_scaled_error_bound','elapsed_s')}

if __name__=='__main__':
 out=pathlib.Path(sys.argv[1]);workers=int(sys.argv[2]);(out/'thermal-enclosure').mkdir(exist_ok=True)
 indices=range(1921) if len(sys.argv)<4 else [int(x) for x in (out/'sample-indices.txt').read_text().split()]
 start=time.monotonic()
 with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as pool:result=list(pool.map(point,[(out,i) for i in indices],chunksize=1))
 save(out/('thermal-enclosure-all.json' if len(sys.argv)<4 else 'thermal-enclosure-sample.json'),dict(status='FAIL' if any(r['status']=='FAIL' for r in result) else 'CONDITIONAL_PASS',results=result,workers=workers,elapsed_s=time.monotonic()-start))
