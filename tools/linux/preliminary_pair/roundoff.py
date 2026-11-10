#!/usr/bin/env python3
"""A posteriori RK arithmetic residual enclosure relative to a declared real RHS.

Includes each retained endpoint/replay discrepancy. Does not bound RK truncation,
the static normalization, or floating-point time selection as a separate quantity.
All bounds are conditional on the exported binary64-anchor RHS model.
"""
import collections,functools,json,pathlib,re,sys,time
from decimal import Decimal as D,localcontext
from high_precision import Coupled,rows,lift,save,BASE
from verify_chemical_reference import PI

@functools.lru_cache(maxsize=2)
def tableau(method):
 source=(BASE/'toolchain/source/gsl-2.7.1/ode-initval2'/(method+'.c')).read_text()
 def rational(s):
  v=s.strip().split('/');return D(v[0].strip()) if len(v)==1 else D(v[0].strip())/D(v[1].strip())
 arrays={name:[rational(x) for x in value.split(',')] for name,value in re.findall(r'static const double (\w+)\[\]\s*=\s*\{([^}]+)\}',source)}
 scalars={name:rational(value) for name,value in re.findall(r'static const double (\w+)\s*=\s*([^;]+);',source)}
 if method=='rk8pd':return [[],[scalars['b21']]]+[arrays['b'+str(i)] for i in range(3,14)],arrays['Abar']
 return [[],[arrays['ah'][0]]]+[arrays['b'+str(i)] for i in range(3,7)],[scalars['c1'],D(0),*[scalars['c'+str(i)] for i in range(3,7)]]

def jacobian_bound(model):
 """Positive analytic majorants on a fixed box; coefficient rounding padded.

R_j=sum a_jk eta_j^(2k+1) T^(6-2k), Lnu_j=sum positive
even powers. Cstar is affine in x inside the single occupied cache cell.
"""
 xmin,xmax=D('-0.0001'),D('0.0002');Y=[D('3e-10'),D('8e-10')]
 tmin=D('1e8')*xmin.exp();tmax=D('1e8')*xmax.exp()
 import bisect
 cell=bisect.bisect_right(model.knots,xmin)-1
 assert cell==bisect.bisect_right(model.knots,xmax)-1
 slope=abs(model.slopes[cell]);cmin=min(model.capacity[cell]+model.slopes[cell]*(x-model.knots[cell]) for x in (xmin,xmax))
 rate=[];deta=[];dx=[];nu=[];neta=[];nx=[]
 for j in range(2):
  a=[model.luminosity[j]*v/(model.kberg*11513*model.kb**(2*k+1)*PI**(2*k+2)) for k,v in enumerate((14680,7560,840,24))]
  rate.append(sum(a[k]*Y[j]**(2*k+1)*tmax**(6-2*k) for k in range(4)))
  deta.append(sum((2*k+1)*a[k]*Y[j]**(2*k)*tmax**(6-2*k) for k in range(4)))
  dx.append(sum((6-2*k)*a[k]*Y[j]**(2*k+1)*tmax**(6-2*k) for k in range(4)))
  b=[model.luminosity[j]*v/(11513*(model.kb*PI)**(2*k+2)) for k,v in enumerate((22020,5670,420,9))]
  nu.append(model.luminosity[j]*tmax**8+sum(b[k]*Y[j]**(2*k+2)*tmax**(6-2*k) for k in range(4)))
  neta.append(sum((2*k+2)*b[k]*Y[j]**(2*k+1)*tmax**(6-2*k) for k in range(4)))
  nx.append(8*model.luminosity[j]*tmax**8+sum((6-2*k)*b[k]*Y[j]**(2*k+2)*tmax**(6-2*k) for k in range(4)))
 photon=model.photon*(lift(2.42)*xmax).exp();heat=model.erg*sum(Y[j]*rate[j] for j in range(2))
 direct=abs(model.pdir)+model.erg*sum(Y[j]*abs(model.sigma[j]) for j in range(2))
 denominator=tmin*cmin;power=direct+heat+sum(nu)+photon
 J=[[D(0)]*3 for _ in range(3)]
 J[0][0]=(model.erg*sum(Y[j]*dx[j] for j in range(2))+sum(nx)+lift(2.42)*photon+power*(1+slope/cmin))/denominator
 for j in range(2):J[0][j+1]=(model.erg*(abs(model.sigma[j])+rate[j]+Y[j]*deta[j])+neta[j])/denominator
 for i in range(2):
  J[i+1][0]=sum(abs(model.z[i][j])*dx[j] for j in range(2))
  for j in range(2):J[i+1][j+1]=abs(model.z[i][j])*deta[j]
 # 110-digit evaluation error is far below this explicit 1e-80 relative padding.
 J=[[v*(1+D('1e-80'))+D('1e-100') for v in row] for row in J]
 return J,dict(x_min=xmin,x_max=xmax,eta_abs=Y,cell=cell,T_min=tmin,T_max=tmax,C_min=cmin,decimal_padding_relative='1e-80')

@functools.lru_cache(maxsize=2)
def cached_model(out,role):return Coupled(out,role)

def audit(out,role,name,directory_override=None,initial_bound=None,method_override=None,write=True):
 start=time.monotonic()
 with localcontext() as ctx:
  ctx.prec=110;model=cached_model(out,role);J,box=jacobian_bound(model)
  directory=directory_override or out/f'{role}-history-traces'/name
  grouped=collections.defaultdict(list)
  for row in rows(directory/'local-rhs.tsv'):grouped[int(row['context'])].append(row)
  method=method_override or ('rk8pd' if name in ('reference-0','reference-1','reference-2') else 'rkf45')
  A,W=tableau(method);n=len(W);bound=list(initial_bound) if initial_bound is not None else [D(0)]*3;records=[];rhsmax=[D(0)]*3;replaymax=[D(0)]*3
  for mapping in rows(directory/'maps.tsv'):
   trace=grouped[int(mapping['context'])];segments=[]
   for r in trace:
    if r['cell']=='-1':segments.append([])
    assert segments,'missing cell-selection record'
    segments[-1].append(r)
   previous=None
   for segment in segments:
    if len(segment)==n+2:stages=segment[1:-1]
    elif len(segment)==n+1 and previous is not None:stages=[previous]+segment[1:-1]
    else:raise ValueError((name,mapping['context'],'unsupported segment',len(segment),n))
    initial=segment[0];endpoint=segment[-1];h=lift(endpoint['t_s'])-lift(initial['t_s']);assert h>0
    y0=[lift(initial[k]) for k in ('x','e','mu')];fs=[[lift(r[k]) for k in ('dx','de','dmu')] for r in stages]
    defects=[]
    for j,r in enumerate(stages):
     yy=tuple(lift(r[k]) for k in ('x','e','mu'))
     assert box['x_min']<yy[0]<box['x_max'] and all(abs(yy[i+1])+bound[i+1]<box['eta_abs'][i] for i in range(2))
     real=model.coupled_rhs(yy);defects.append([abs(fs[j][i]-real[i])+D('1e-90') for i in range(3)])
    errors=[]
    for j in range(n):
     residual=[abs(lift(stages[j][key])-y0[i]-h*sum(A[j][k]*fs[k][i] for k in range(j)))+D('1e-90') for i,key in enumerate(('x','e','mu'))]
     errors.append([bound[i]+residual[i]+h*sum(abs(A[j][k])*(defects[k][i]+sum(J[i][v]*errors[k][v] for v in range(3))) for k in range(j)) for i in range(3)])
     yy=[lift(stages[j][key]) for key in ('x','e','mu')]
     assert box['x_min']<yy[0]-errors[-1][0] and yy[0]+errors[-1][0]<box['x_max'],'stage thermal enclosure leaves real box'
     assert all(abs(yy[i+1])+errors[-1][i+1]<box['eta_abs'][i] for i in range(2)),'stage chemical enclosure leaves real box'
     for i in range(3):rhsmax[i]=max(rhsmax[i],defects[j][i])
    residual=[abs(lift(endpoint[key])-y0[i]-h*sum(W[j]*fs[j][i] for j in range(n)))+D('1e-90') for i,key in enumerate(('x','e','mu'))]
    bound=[bound[i]+residual[i]+h*sum(abs(W[j])*(defects[j][i]+sum(J[i][v]*errors[j][v] for v in range(3))) for j in range(n)) for i in range(3)]
    assert box['x_min']<lift(endpoint['x'])-bound[0] and lift(endpoint['x'])+bound[0]<box['x_max'],'endpoint thermal enclosure leaves box'
    assert all(abs(lift(endpoint[key]))+bound[i+1]<box['eta_abs'][i] for i,key in enumerate(('e','mu'))),'endpoint chemical enclosure leaves box'
    previous=endpoint
   mismatch=[abs(lift(mapping[key])-lift(mapping['replay_'+key])) for key in ('x','e','mu')]
   bound=[v+w for v,w in zip(bound,mismatch)]
   for i in range(3):replaymax[i]=max(replaymax[i],mismatch[i])
   assert bound[0]<D('1e-6') and max(bound[1:])<D('1e-12'),'enclosure leaves box'
   records.append(dict(time_s=mapping['right'],accumulated_arithmetic_enclosure=bound.copy(),retained_replay_difference=mismatch))
  result=dict(status='CONDITIONAL_ARITHMETIC_ENCLOSURE',role=role,history=name,method=method,J=J,box=box,records=records,maximum_rhs_defect=rhsmax,maximum_replay_difference=replaymax,final_enclosure=bound,elapsed_s=time.monotonic()-start,assumptions=['Exported binary64 Cstar cells and evaluated photon/direct-power anchors define the real RHS. Their static normalization uncertainty is separate.','GSL rational tableau defines the exact RK comparison; observed stage/endpoint residuals include binary64 coefficient arithmetic and reconstructed time-span effects.','Excludes RK truncation and reference observation interpolation.','Decimal elementary functions and arithmetic errors are covered by the stated conservative padding.'],inputs={str(p):__import__('run').digest(p) for p in (directory/'local-rhs.tsv',directory/'maps.tsv')} )
  if write:
   target=out/'roundoff';target.mkdir(exist_ok=True);save(target/f'{role}-{name}.json',result)
  return result

if __name__=='__main__':audit(pathlib.Path(sys.argv[1]),sys.argv[2],sys.argv[3])
