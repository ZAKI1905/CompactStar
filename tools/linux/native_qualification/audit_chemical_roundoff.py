#!/usr/bin/env python3
"""100-digit residual audit of recorded RK arithmetic, not a new trajectory.

The bound is conditional on the recorded thermal stage inputs. It must not be
advertised as a complete coupled-IVP or original-t0 accumulated-roundoff bound.
"""
import collections,csv,json,pathlib,re,sys
from decimal import Decimal,localcontext

def rows(p):
 with pathlib.Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def audit(root,coefficients,tangent,gsl):
 with localcontext() as ctx:
  ctx.prec=100
  D=lambda x:Decimal.from_float(float(x))
  def rational(s):
   v=s.strip().split('/')
   return Decimal(v[0].strip()) if len(v)==1 else Decimal(v[0].strip())/Decimal(v[1].strip())
  def tableau(method):
   source=(pathlib.Path(gsl)/(method+'.c')).read_text()
   arrays={name:[rational(x) for x in value.split(',')] for name,value in re.findall(r'static const double (\w+)\[\]\s*=\s*\{([^}]+)\}',source)}
   scalars={name:rational(value) for name,value in re.findall(r'static const double (\w+)\s*=\s*([^;]+);',source)}
   if method=='rk8pd':return [[],[scalars['b21']]]+[arrays['b'+str(i)] for i in range(3,14)],arrays['Abar']
   return [[],[arrays['ah'][0]]]+[arrays['b'+str(i)] for i in range(3,7)],[scalars['c1'],Decimal(0),*[scalars['c'+str(i)] for i in range(3,7)]]
  tables={m:tableau(m) for m in ['rk8pd','rkf45']}
  c=rows(coefficients)[0];tangent=rows(tangent)[0]
  Z=[[D(c['Z'+str(i)+str(j)]) for j in range(2)] for i in range(2)]
  sigma=[-D(tangent[k])*D(-2.413652026364133e37) for k in ['t_e','t_mu']]
  source=[-sum(Z[i][j]*sigma[j] for j in range(2)) for i in range(2)]
  L=[D(c[k]) for k in ['Ltilde_Me','Ltilde_Mmu']]
  kb=D(8.61733326214518e-5*1e-6);erg=D(1.602176634e33/1e39);kberg=D(float(kb)*float(erg))
  pi=Decimal('3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068')
  det=Z[0][0]*Z[1][1]-Z[0][1]*Z[1][0]
  def energy(y):
   e,m=y;return (Z[1][1]*e*e-(Z[0][1]+Z[1][0])*e*m+Z[0][0]*m*m)/(2*det)
  def reaction(T,eta):
   xi=eta/(kb*T);u2=(xi/pi)**2
   return xi/pi**2*(14680+u2*(7560+u2*(840+24*u2)))/11513
  def rhs(row):
   T=Decimal('1e8')*D(row['x']).exp();eta=[D(row[k]) for k in ['e','mu']]
   rates=[L[j]/kberg*T**7*reaction(T,eta[j]) for j in range(2)]
   return [source[i]-sum(Z[i][j]*rates[j] for j in range(2)) for i in range(2)]
  grouped=collections.defaultdict(list)
  for row in rows(pathlib.Path(root)/'local-rhs.tsv'):grouped[int(row['context'])].append(row)
  records=[]
  for context,trace in sorted(grouped.items()):
   offset=(context-1)%42+1
   method='rkf45' if (offset<=36 and (offset-1)%6==3) or offset>=40 else 'rk8pd'
   A,W=tables[method];n=len(W)
   # Chemical Jacobian bound on a box enlarged by 1%; temperature by 1 K.
   Tmax=max(Decimal('1e8')*D(r['x']).exp() for r in trace)+1
   etamax=[max(abs(D(r[k])) for r in trace)*Decimal('1.01') for k in ['e','mu']]
   slopes=[]
   for j in range(2):
    u2=(etamax[j]/(kb*Tmax*pi))**2
    slopes.append(L[j]/kberg*Tmax**6/(kb*pi**2)*(14680+3*7560*u2+5*840*u2**2+7*24*u2**3)/11513)
   J=[[abs(Z[i][j])*slopes[j] for j in range(2)] for i in range(2)]
   segments=[]
   for r in trace:
    if r['cell']=='-1':segments.append([])
    assert segments,'trace does not start with cell selection'
    segments[-1].append(r)
   bound=[Decimal(0),Decimal(0)];max_rhs=[Decimal(0),Decimal(0)];max_stage=[Decimal(0),Decimal(0)]
   supported=True;previous=None
   for segment in segments:
    if len(segment)==n+2:stages=segment[1:-1]
    elif len(segment)==n+1 and previous is not None:stages=[previous]+segment[1:-1]
    else:supported=False;break # rejected/truncated attempts cannot silently enter the audit
    assert len(stages)==n
    initial=segment[0];endpoint=segment[-1]
    h=D(endpoint['t_s'])-D(initial['t_s'])
    assert h>0
    y0=[D(initial[k]) for k in ['e','mu']]
    fs=[[D(r[k]) for k in ['de','dmu']] for r in stages]
    defects=[[abs(fs[j][i]-rhs(stages[j])[i]) for i in range(2)] for j in range(n)]
    errors=[]
    for j in range(n):
     residual=[abs(D(stages[j][key])-y0[i]-h*sum(A[j][k]*fs[k][i] for k in range(j))) for i,key in enumerate(['e','mu'])]
     error=[bound[i]+residual[i]+h*sum(abs(A[j][k])*(defects[k][i]+sum(J[i][v]*errors[k][v] for v in range(2))) for k in range(j)) for i in range(2)]
     errors.append(error)
     for i in range(2):max_rhs[i]=max(max_rhs[i],defects[j][i]);max_stage[i]=max(max_stage[i],residual[i])
    residual=[abs(D(endpoint[key])-y0[i]-h*sum(W[j]*fs[j][i] for j in range(n))) for i,key in enumerate(['e','mu'])]
    bound=[bound[i]+residual[i]+h*sum(abs(W[j])*(defects[j][i]+sum(J[i][v]*errors[j][v] for v in range(2))) for j in range(n)) for i in range(2)]
    assert all(bound[i]<etamax[i]/100 for i in range(2)), 'roundoff enclosure leaves Jacobian box'
    previous=endpoint
   record={'context':context,'point_position':(context-1)//42,'offset':offset,'method':method,'supported_trace':supported,'accepted_steps':len(segments),'conditional_on_recorded_thermal_stages':True}
   if supported:
    y=[D(trace[-1][k]) for k in ['e','mu']];E=energy(y)
    energy_bound=max(abs(energy([y[0]+se*bound[0],y[1]+sm*bound[1]])-E) for se in [-1,1] for sm in [-1,1])
    record.update(time_s=float(trace[-1]['t_s']),chemical_roundoff_enclosure_MeV=[str(x) for x in bound],maximum_rhs_defect_MeV_s=[str(x) for x in max_rhs],maximum_stage_residual_MeV=[str(x) for x in max_stage],Echem_roundoff_enclosure_MeV=str(energy_bound),Echem_relative_roundoff_enclosure=str(energy_bound/abs(E)))
   records.append(record)
  return {'schema':'conditional-chemical-roundoff-audit-v1','precision_digits':100,'records':records,
    'scope':'Chemical RK arithmetic and chemical RHS at recorded thermal stages; exact frozen binary64 inputs, exact rational RK coefficients; not an independent coupled or original-t0 error bound.',
    'excluded':['thermal-stage accumulated roundoff','main-prefix error','truncation error','reference interpolation error'],
    'stage_A_gate':False}

if __name__=='__main__':
 result=audit(*sys.argv[1:5]);pathlib.Path(sys.argv[5]).write_text(json.dumps(result,indent=2)+'\n')
 print('Audited',len(result['records']),'contexts; complete coupled accumulated-roundoff bound remains unestablished')
