#!/usr/bin/env python3
"""Piecewise R20 enclosure on unchanged scientific-grid trapezoids.

Missing grids remain incomplete. No native residual is corrected or replaced.
Rectangle bounds tolerate Cstar knots; thermal energy uses the exact piecewise
primitive only to bound quadrature error of the original stored trapezoid.
"""
import bisect,concurrent.futures,json,math,pathlib,sys
from decimal import Decimal as D,localcontext
from piecewise import MODEL,Coupled,domain,SCALE,END,PI
from run import save,digest
from high_precision import rows,lift
from verify_chemical_reference import corner_uncertainty

def capacity(m,x):
 c=bisect.bisect_right(m.knots,x)-1
 assert 0<=c<len(m.slopes)
 return m.capacity[c]+m.slopes[c]*(x-m.knots[c])

def primitive(m,x0,x1):
 if x0==x1:return D(0)
 sign=D(1) if x1>x0 else D(-1);lo=min(x0,x1);hi=max(x0,x1)
 points=[lo]+[x for x in m.knots if lo<x<hi]+[hi];total=D(0)
 for a,b in zip(points,points[1:]):
  cell=bisect.bisect_right(m.knots,(a+b)/2)-1;slope=m.slopes[cell]
  total+=D('1e8')*(b.exp()*(capacity(m,b)-slope)-a.exp()*(capacity(m,a)-slope))
 return sign*total

def outgoing(m,x,eta):
 T=D('1e8')*x.exp();nu=D(0)
 for j in range(2):
  u2=(eta[j]/(m.kb*T*PI))**2
  nu+=m.luminosity[j]*T**8*(1+u2*(22020+u2*(5670+u2*(420+9*u2)))/11513)
 return nu+m.photon*(D.from_float(2.42)*x).exp()

def outgoing_range(m,y,u):
 # Positive monomials in T and |eta|, with powers T^(8), T^(6),...,T^0.
 xl=y[0]-u;xh=y[0]+u
 el=[max(D(0),abs(y[j+1])-u*SCALE[j+1]) for j in range(2)]
 eh=[abs(y[j+1])+u*SCALE[j+1] for j in range(2)]
 return outgoing(m,xl,el)*(1-D('1e-75')),outgoing(m,xh,eh)*(1+D('1e-75'))

def assess_one(task):
 out,role,directory=task
 with localcontext() as ctx:
  ctx.prec=110;m=Coupled(MODEL,role);box=domain(m);r20=rows(directory/'R20-grids.tsv')
  if len(r20)!=4:return dict(role=role,configuration=directory.name,status='INCOMPLETE',grids=len(r20))
  hp=json.loads((out/f'hp-{role}-1000000.json').read_text())['records'];n=1920;dt=END/n
  data=rows(directory/'diagnostics-1920.tsv');unc=rows(directory/'uncertainties-1920.tsv');assert len(data)==len(unc)==1921
  pre=json.loads((out/'preflight.json').read_text());forcing=D(pre['roles'][role]['normalization_rhs_bound'])
  norm=json.loads((MODEL/'normalization.json').read_text())['roles'][role];dph=D(norm['photon_normalization_error'])
  ymax=[];output_error=[];quad_low=D(0);quad_high=D(0);hp_trap=D(0);node_power=[];max_r18=D(0);event_power=[]
  Cmax=max(m.capacity[k] for k in range(min(box['cells']),max(box['cells'])+2))
  for i,(r,h) in enumerate(zip(data,hp)):
   t=lift(r['t_s']);assert t==D(h['t_s']);y=tuple(D(h[k]) for k in ('x','e','mu'));u=D(h['U'])+forcing*t*(box['L']*t).exp()
   lo,hi=outgoing_range(m,y,u);normal=dph*(D.from_float(2.42)*(y[0]+u)).exp();lo-=normal;hi+=normal
   value=outgoing(m,y[0],y[1:]);node_power.append(value)
   native=sum(lift(r[k]) for k in ('Lnu_full_erg_s','Lgamma_erg_s','Lother_erg_s','L_out_fluid_inf_erg_s'))
   output_error.append(max(abs(native-lo),abs(native-hi)))
   G=max(D(1),sum(abs(lift(r[k])) for k in ('P_dir_actual_erg_s','LH_erg_s','DeltaLnu_erg_s','Lnu_eq_erg_s','Lgamma_erg_s','Lother_erg_s','L_out_fluid_inf_erg_s')))
   max_r18=max(max_r18,abs(lift(r['R18_residual_erg_s']))/G)
   eu=max(lift(unc[i][k])/SCALE[j] for j,k in enumerate(('Ex','Ee','Emu')))
   ey=tuple(lift(unc[i][k]) for k in ('x','e','mu'))
   el,eh=outgoing_range(m,ey,eu);ep=outgoing(m,ey[0],ey[1:]);event_power.append(max(abs(ep-el),abs(eh-ep)))
   if i<n:
    reach=u+(box['M']+forcing)*dt
    assert box['lo']<y[0]-reach and y[0]+reach<box['hi']
    assert all(abs(y[j+1])+reach*SCALE[j+1]<box['Y'][j] for j in range(2))
    lo,hi=outgoing_range(m,y,reach);normal=dph*(D.from_float(2.42)*(y[0]+reach)).exp()
    quad_low+=dt*(lo-normal);quad_high+=dt*(hi+normal)
  hp_trap=sum(dt*(a+b)/2 for a,b in zip(node_power,node_power[1:]))
  outgoing_quad=max(abs(hp_trap-quad_low),abs(hp_trap-quad_high))
  Umain=sum(dt*(a+b)/2 for a,b in zip(output_error,output_error[1:]))
  end=hp[-1];ye=tuple(D(end[k]) for k in ('x','e','mu'));ue=D(end['U'])+forcing*END*(box['L']*END).exp()
  native_T=lift(data[-1]['Tinf_K']);true_T=D('1e8')*ye[0].exp();temp_error=abs(native_T-true_T)+true_T*(ue.exp()-1)
  Umain+=Cmax*temp_error
  chemical_u=D(0)
  if role=='source':chemical_u=corner_uncertainty(m.energy,ye[1:],tuple(ue*SCALE[j] for j in (1,2)))[0]*m.erg
  Umain+=abs(lift(data[-1]['Echem_MeV'])*m.erg-m.energy(ye[1:])*m.erg)+chemical_u
  temps=[lift(r['Tinf_K']) for r in data];xs=[(T/D('1e8')).ln() for T in temps]
  native_thermal=sum((lift(a['Cstar_erg_K'])+lift(b['Cstar_erg_K']))*(tb-ta)/2 for a,b,ta,tb in zip(data,data[1:],temps,temps[1:]))
  ideal_thermal=sum((capacity(m,xa)+capacity(m,xb))*(tb-ta)/2 for xa,xb,ta,tb in zip(xs,xs[1:],temps,temps[1:]))
  thermal_quad=abs(ideal_thermal-primitive(m,xs[0],xs[-1]));capacity_eval=abs(native_thermal-ideal_thermal)
  # Independently recompute every original nested-grid operation in Decimal.
  rounding=[]
  for grid in r20:
   ds=rows(directory/f"diagnostics-{grid['intervals']}.tsv");h=END/D(grid['intervals'])
   du=sum((lift(a['Cstar_erg_K'])+lift(b['Cstar_erg_K']))*(lift(b['Tinf_K'])-lift(a['Tinf_K']))/2 for a,b in zip(ds,ds[1:]))
   outgoing_sum=sum(h/2*sum(lift(a[k])+lift(b[k]) for k in ('Lnu_full_erg_s','Lgamma_erg_s','Lother_erg_s','L_out_fluid_inf_erg_s')) for a,b in zip(ds,ds[1:]))
   eq=m.erg*lift(ds[0]['mu_B_inf_MeV'])*(lift(ds[-1]['B_count'])-lift(ds[0]['B_count']))
   chem=m.erg*(lift(ds[-1]['Echem_MeV'])-lift(ds[0]['Echem_MeV']))
   rounding.append(abs(lift(grid['residual'])-(du+outgoing_sum+eq+chem)))
  last=r20[-1];N=lift(last['N']);diffs=[abs(lift(a['residual'])-lift(b['residual'])) for a,b in zip(r20,r20[1:])]
  count_bound=D(0)
  if role=='source':
   # A full ulp bound for the difference of two independently rounded counts.
   count_bound=lift(math.ulp(float(data[0]['B_count'])))*abs(lift(data[0]['mu_B_inf_MeV']))*m.erg
  # Explicit additional event bound. Stored reconstruction already includes events;
  # this double-counting is disclosed and conservative.
  Uevent=sum(dt*(a+b)/2 for a,b in zip(event_power,event_power[1:]))
  Uevent+=sum(2*Cmax*lift(r['Tinf_K'])*(lift(u['Ex']).exp()-1) for r,u in zip(data,unc))
  if role=='source':
   for idx in (0,-1):
    eta=tuple(lift(data[idx][k]) for k in ('eta_e_MeV','eta_mu_MeV'))
    err=tuple(lift(unc[idx][k]) for k in ('Ee','Emu'))
    Uevent+=corner_uncertainty(m.energy,eta,err)[0]*m.erg
  Urecon=lift(last['reconstruction_uncertainty'])
  Uquad=2*max(diffs)+outgoing_quad+thermal_quad+capacity_eval+rounding[-1]
  total=(abs(lift(last['residual']))+count_bound+Umain+Urecon+Uevent+Uquad)/N
  contraction=diffs[-1]<=diffs[-2]/2 or all(diffs[i]<=rounding[i]+rounding[i+1] for i in range(3))
  result=dict(role=role,configuration=directory.name,status='CONDITIONAL_LEDGER_ENCLOSURE_CANDIDATE',R20_original=lift(last['residual']),N=N,U_main=Umain,U_reconstruction_including_stored_events=Urecon,U_event_extra=Uevent,U_quad=Uquad,U_count=count_bound,total_over_N=total,limit=D('2e-4'),reconstruction_event_allocation_over_N=(Urecon+Uevent)/N,allocation_limit=D('5e-6'),max_R18_over_G=max_r18,nested_differences=diffs,independent_arithmetic_bounds=rounding,contraction=contraction,outgoing_rectangle_quadrature_bound=outgoing_quad,thermal_piecewise_primitive_quadrature_bound=thermal_quad,capacity_evaluation_bound=capacity_eval,maximum_depletion=max(abs(lift(r['DeltaB_over_B0'])) for r in data),maximum_frozen_utilization=max(lift(r['max_frozen_utilization']) for r in data),finite_T_omitted_scale_K=max(lift(r['finite_T_omitted_floor_erg_s']) for r in data)*END/min(c['C_min'] for c in box['cells'].values()))
  result['gates']=dict(R20=total<=D('2e-4'),allocation=(Urecon+Uevent)/N<=D('5e-6'),R18=max_r18<=D('1e-10'),contraction=contraction)
  return result

if __name__=='__main__':
 out=pathlib.Path(sys.argv[1]);workers=int(sys.argv[2]);tasks=[(out,role,p) for role in ('source','control') for p in sorted((out/'native'/role).glob('tier-*')) if (p/'R20-grids.tsv').exists()]
 with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as pool:result=list(pool.map(assess_one,tasks))
 save(out/'ledger-assessment.json',dict(status='CANDIDATE_PARTIAL_UNTIL_ALL_MATRIX_GATES',workers=workers,configurations=result))
