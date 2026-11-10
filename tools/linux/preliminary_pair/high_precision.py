#!/usr/bin/env python3
"""Independent Decimal calculations; numerical estimates are explicitly conditional.

No production solver is imported. Reuses the reviewed chemical model, aligned
RK4 and thermal sensitivity calculation. Frozen binary64 inputs are lifted exactly.
"""
import argparse,bisect,collections,concurrent.futures,csv,json,math,os,pathlib,sys,time
from decimal import Decimal as D,localcontext
from run import ROOT,BASE,digest,save
sys.path.insert(0,str(ROOT/'tools/linux/native_qualification'))
from verify_chemical_reference import Chemistry,lift,rows,Hermite,convergence,corner_uncertainty,fidelity
from resolve_chemical_reference import aligned_rk4,sensitivity_bound,KS
COEFF=BASE/'ba13-recovery/run1/raw/coefficients.tsv'
TANGENT=COEFF.parent/'tangent.tsv'
SOURCE=BASE/'stage-a-parallel-20261110-r1/source-stage/run/trajectory/source'

class Coupled(Chemistry):
 def __init__(self,out,role):
  super().__init__(rows(COEFF)[0],rows(TANGENT)[0])
  self.role=role
  anchors={r['quantity']:lift(r['value']) for r in rows(out/'thermal-model/anchors.tsv')}
  cells=rows(out/'thermal-model/capacity.tsv')
  self.knots=[(lift(r['T_MeV'])/(D('1e8')*self.kb)).ln() for r in cells]
  self.capacity=[lift(r['C_erg_K']) for r in cells]
  self.slopes=[(b-a)/(y-x) for a,b,x,y in zip(self.capacity,self.capacity[1:],self.knots,self.knots[1:])]
  self.photon=anchors['Lgamma'];self.pdir=anchors['P_dir_eq'] if role=='source' else D(0)
  if role=='control':self.sigma=[D(0),D(0)];self.source=[D(0),D(0)]
 def powers(self,y):
  x,e,mu=y;T=D('1e8')*x.exp();eta=(e,mu)
  cell=bisect.bisect_right(self.knots,x)-1
  assert 0<=cell<len(self.slopes),'outside cache domain'
  C=self.capacity[cell]+self.slopes[cell]*(x-self.knots[cell])
  rates=[self.luminosity[j]/self.kberg*T**7*self.reaction(T,eta[j]) for j in range(2)]
  LH=self.erg*sum(eta[j]*rates[j] for j in range(2))
  nu=D(0)
  for j in range(2):
   from verify_chemical_reference import PI
   u2=(eta[j]/(self.kb*T*PI))**2
   nu+=self.luminosity[j]*T**8*(1+u2*(22020+u2*(5670+u2*(420+9*u2)))/11513)
  photon=self.photon*(lift(2.42)*x).exp()
  direct=self.pdir+self.erg*sum(eta[j]*self.sigma[j] for j in range(2))
  return T,C,direct,LH,nu,photon,rates
 def coupled_rhs(self,y):
  T,C,direct,LH,nu,photon,rates=self.powers(y)
  return ((direct+LH-nu-photon)/(T*C),*[self.source[i]-sum(self.z[i][j]*rates[j] for j in range(2)) for i in range(2)])

def step(rhs,y,h):
 a=rhs(y);b=rhs(tuple(v+h*w/2 for v,w in zip(y,a)))
 c=rhs(tuple(v+h*w/2 for v,w in zip(y,b)));d=rhs(tuple(v+h*w for v,w in zip(y,c)))
 return tuple(v+h*(a1+2*b1+2*c1+d1)/6 for v,a1,b1,c1,d1 in zip(y,a,b,c,d))

def global_reference(out,role,n,precision):
 start=time.monotonic();directory=out/'high-precision';directory.mkdir(exist_ok=True)
 with localcontext() as ctx:
  ctx.prec=precision;model=Coupled(out,role)
  main_times={lift(r['t_after_s']) for r in rows(SOURCE/'tier-0-cap-0/main.accepted_states.tsv')}
  times=sorted(main_times|{D(31557600)*i/1920 for i in range(1921)})
  y=(D(0),D(0),D(0));t=D(0);records=[]
  for target in times:
   if target>t:
    h=(target-t)/n
    for _ in range(n):y=step(model.coupled_rhs,y,h)
   records.append(dict(t_s=target,x=y[0],e=y[1],mu=y[2]))
   t=target
  save(directory/f'{role}-n{n}-p{precision}.json',dict(role=role,n=n,precision=precision,records=records,elapsed_s=time.monotonic()-start,job=os.environ.get('SLURM_JOB_ID'),assumptions=['Exact retained binary64 Cstar cell values, coefficient values, evaluated photon and equilibrium direct-power anchors.','Static normalization error and continuous RHS enclosure not discharged by this calculation.','Twice refinement differences estimate discretization only conditional on continuing contraction.']))

def local_point(args):
 out,point=args;start=time.monotonic()
 with localcontext() as ctx:
  ctx.prec=110;model=Chemistry(rows(COEFF)[0],rows(TANGENT)[0])
  directory=next(p for p in (out/'sample'/f'point-{point}',out/'full'/f'point-{point}') if p.exists())
  cp=json.loads((directory/'packet/checkpoint.json').read_text())
  if cp['method']=='main-endpoint':return dict(point=point,status='EXACT_ENDPOINT',elapsed_s=time.monotonic()-start)
  trace=rows(directory/'local-rhs.tsv');steps=rows(directory/'witness-steps.tsv');witnesses=rows(directory/'witnesses.tsv')
  fid=fidelity(model,trace)
  chosen=next(r for r in witnesses if r['method']=='0' and r['subdivisions']=='8')
  context=[r for r in trace if r['context']==chosen['context']]
  derivatives=collections.defaultdict(set)
  for r in context:derivatives[(lift(r['t_s']),*(lift(r[k]) for k in ('x','e','mu')))].add(lift(r['dx']))
  left=lift(cp['main_left_time_s']);initial=tuple(lift(v) for v in cp['main_left_state'])
  def node(t,y):
   ds=derivatives[(t,*y)];assert len(ds)==1,(point,t,ds)
   return t-left,y[0],ds.pop()
  nodes=[node(left,initial)]
  previous=left
  for r in steps:
   if r['method']=='0' and r['subdivisions']=='8':
    assert lift(r['left'])==previous
    previous=lift(r['right']);nodes.append(node(previous,tuple(lift(r[k]) for k in ('x','e','mu'))))
  assert previous==lift(cp['time_s'])
  path=Hermite(nodes)
  rho=max([abs(path(lift(r['right'])-left)-lift(r['x'])) for r in steps if r['subdivisions']=='16']+[abs(path(nodes[-1][0])-lift(cp['O2_state'][0]))])
  vals={8*k:aligned_rk4(model.rhs,nodes,initial[1:],k) for k in KS}
  conv=convergence(vals);ref=vals[512]
  discretization=tuple(2*abs(ref[i]-vals[256][i]) for i in range(2))
  sensitivity=sensitivity_bound(model,nodes,initial[1:])
  unc=tuple(v+sensitivity['uniform_thermal_bound_MeV'] for v in discretization)
  sr={r['channel']:r for r in rows(directory/'packet/state.tsv')}
  dr={r['channel']:r for r in rows(directory/'packet/diagnostics.tsv')}
  errors=[abs(lift(cp['O2_state'][i+1])-ref[i]) for i in range(2)]
  ratios=[(errors[i]+unc[i])/lift(sr[k]['propagated_state_radius']) for i,k in enumerate(('eta_e','eta_mu'))]
  er=model.energy(ref);ue,_=corner_uncertainty(model.energy,ref,unc)
  ratio=(abs(lift(dr['Echem']['O2'])-er)+ue)/lift(dr['Echem']['U'])
  result=dict(point=point,status='CONDITIONAL_PASS' if max(ratios+[ratio])<=1 and rho<=D('1e-13') and all(r['status']=='PASS' for r in conv) and fid['status']=='PASS' and sensitivity['status']=='PASS' else 'FAIL',rhs_fidelity=fid,convergence=conv,thermal_rho=rho,thermal_sensitivity=sensitivity,chemical_reference=ref,discretization_estimate=discretization,chemical_absolute_error=errors,chemical_state_ratios=ratios,chemical_energy_ratio=ratio,chemical_energy_error=abs(lift(dr['Echem']['O2'])-er),chemical_energy_uncertainty=ue,elapsed_s=time.monotonic()-start,input_hashes={str(p):digest(p) for p in directory.glob('*.tsv')},conditional_assumptions=['Uniform |x_true-x_Hermite|<=1e-12 throughout local interval.','Continued RK4 refinement contraction.','Exact retained left state; main-prefix uncertainty assessed separately.'])
  save(out/'local-hp'/f'point-{point}.json',result)
  return {k:result[k] for k in ('point','status','thermal_rho','chemical_state_ratios','chemical_energy_ratio','elapsed_s')}

def local(out,mode,workers):
 (out/'local-hp').mkdir(exist_ok=True)
 indices=[int(x) for x in (out/(mode+'-indices.txt')).read_text().split()]
 start=time.monotonic()
 with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as pool:
  results=list(pool.map(local_point,[(out,i) for i in indices],chunksize=1))
 save(out/f'local-hp-{mode}.json',dict(results=results,workers=workers,elapsed_s=time.monotonic()-start,status='FAIL' if any(r['status']=='FAIL' for r in results) else 'CONDITIONAL_PASS'))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('mode',choices=['global','sample','full']);p.add_argument('out',type=pathlib.Path);p.add_argument('--role',choices=['source','control']);p.add_argument('--n',type=int);p.add_argument('--precision',type=int,default=110);p.add_argument('--workers',type=int,default=24);a=p.parse_args()
 if a.mode=='global':global_reference(a.out,a.role,a.n,a.precision)
 else:local(a.out,a.mode,a.workers)
