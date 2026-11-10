#!/usr/bin/env python3
"""Piecewise extension of the accepted Cauchy/positive-Jacobian verification.

Frozen inputs, same real RHS, 110-digit Decimal with retained padding assumptions.
Crossing bridges use a Lipschitz local remainder, never a one-cell Cauchy tail.
This is verification code, not a replacement production integrator.
"""
import bisect,json,pathlib,sys
from decimal import Decimal as D,localcontext
ROOT=pathlib.Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'tools/linux/preliminary_pair'))
from high_precision import Coupled,step,lift,rows,BASE
from verify_chemical_reference import PI
SCALE=[D(1),D('1e-5'),D('1e-5')]
YEAR=D(31557600);END=2048*YEAR
MODEL=BASE/'preliminary-pair-20261110-r1'

def constants(model,xmin,xmax,Y,cell,R=D('.1'),scale=SCALE):
 """Same polynomial majorants as Stage A, on an explicit cell/domain."""
 tmin=D('1e8')*xmin.exp();tmax=D('1e8')*xmax.exp()
 slope=abs(model.slopes[cell]);cmin=min(model.capacity[cell]+model.slopes[cell]*(x-model.knots[cell]) for x in (xmin,xmax))
 assert cmin>0
 def majorants(t,Y):
  rate=[];de=[];dx=[];nu=[];ne=[];nx=[];dn=[]
  for j in range(2):
   a=[model.luminosity[j]*v/(model.kberg*11513*model.kb**(2*k+1)*PI**(2*k+2)) for k,v in enumerate((14680,7560,840,24))]
   rate.append(sum(a[k]*Y[j]**(2*k+1)*t**(6-2*k) for k in range(4)))
   de.append(sum((2*k+1)*a[k]*Y[j]**(2*k)*t**(6-2*k) for k in range(4)))
   dx.append(sum((6-2*k)*a[k]*Y[j]**(2*k+1)*t**(6-2*k) for k in range(4)))
   b=[model.luminosity[j]*v/(11513*(model.kb*PI)**(2*k+2)) for k,v in enumerate((22020,5670,420,9))]
   dn.append(sum(b[k]*Y[j]**(2*k+2)*t**(6-2*k) for k in range(4)))
   nu.append(model.luminosity[j]*t**8+dn[-1])
   ne.append(sum((2*k+2)*b[k]*Y[j]**(2*k+1)*t**(6-2*k) for k in range(4)))
   nx.append(8*model.luminosity[j]*t**8+sum((6-2*k)*b[k]*Y[j]**(2*k+2)*t**(6-2*k) for k in range(4)))
  return rate,de,dx,nu,ne,nx,dn
 rate,de,dx,nu,ne,nx,dn=majorants(tmax,Y)
 photon=model.photon*(lift(2.42)*xmax).exp();heat=model.erg*sum(Y[j]*rate[j] for j in range(2))
 direct=abs(model.pdir)+model.erg*sum(Y[j]*abs(model.sigma[j]) for j in range(2))
 denominator=tmin*cmin;power=direct+heat+sum(nu)+photon
 J=[[D(0)]*3 for _ in range(3)]
 J[0][0]=(model.erg*sum(Y[j]*dx[j] for j in range(2))+sum(nx)+lift(2.42)*photon+power*(1+slope/cmin))/denominator
 for j in range(2):J[0][j+1]=(model.erg*(abs(model.sigma[j])+rate[j]+Y[j]*de[j])+ne[j])/denominator
 for i in range(2):
  J[i+1][0]=sum(abs(model.z[i][j])*dx[j] for j in range(2))
  for j in range(2):J[i+1][j+1]=abs(model.z[i][j])*de[j]
 J=[[v*(1+D('1e-80'))+D('1e-100') for v in row] for row in J]
 L=max(sum(J[i][j]*scale[j]/scale[i] for j in range(3)) for i in range(3))
 Mreal=max([power/denominator]+[(abs(model.source[i])+sum(abs(model.z[i][j])*rate[j] for j in range(2)))/scale[i+1] for i in range(2)])*(1+D('1e-75'))+D('1e-90')
 Clo=cmin-slope*R;assert Clo>0
 Tc=D('1e8')*(xmax+R).exp();Yc=[Y[i]+R*scale[i+1] for i in range(2)]
 rates,_,_,nus,_,_,_=majorants(Tc,Yc)
 pc=abs(model.pdir)+model.erg*sum(Yc[j]*abs(model.sigma[j]) for j in range(2))+model.erg*sum(Yc[j]*rates[j] for j in range(2))+sum(nus)+model.photon*(lift(2.42)*(xmax+R)).exp()
 M=max([pc/(D('1e8')*(xmin-R).exp()*Clo)]+[(abs(model.source[i])+sum(abs(model.z[i][j])*rates[j] for j in range(2)))/scale[i+1] for i in range(2)])*(1+D('1e-75'))+D('1e-90')
 return dict(cell=cell,x_min=xmin,x_max=xmax,Y=Y,C_min=cmin,complex_C_min=Clo,J=J,L=L,M=M,Mreal=Mreal,R=R,r=R/(2*M)*(1-D('1e-75')),scale=scale)

def domain(model):
 lo,hi=(D('-.03'),D('.2'));Y=[D('4e-7'),D('8e-7')]
 # Positive majorants are valid also for the exact-zero control.
 cells=range(bisect.bisect_right(model.knots,lo)-1,bisect.bisect_right(model.knots,hi))
 cs=[constants(model,max(lo,model.knots[c]),min(hi,model.knots[c+1]),Y,c) for c in cells]
 return dict(lo=lo,hi=hi,Y=Y,cells={c['cell']:c for c in cs},L=max(c['L'] for c in cs),M=max(c['Mreal'] for c in cs))

def norm(y):return max(abs(y[i])/SCALE[i] for i in range(3))

def integrate(model,targets,max_step=D('1000000')):
 """Original-t0 offline reference; a checked tube and explicit crossing bridges.

Crossing remainder <= L M h^2 exp(2 L h) by exact/RK4 Euler deviations.
Single-cell remainder is R(h/r)^5/(1-h/r). Prefixes propagate with exp(L h).
"""
 box=domain(model);y=(D(0),)*3;t=D(0);u=D(0);records=[];crossings=[];steps=0;bridges=0
 L=box['L'];M=box['M'];lastcell=bisect.bisect_right(model.knots,y[0])-1
 for target in targets:
  while t<target:
   cell=bisect.bisect_right(model.knots,y[0])-1;c=box['cells'][cell]
   distance=min(y[0]-model.knots[cell],model.knots[cell+1]-y[0])
   h=min(max_step,target-t)
   # Keeping the entire true/reference tube within a cell justifies its analytic tail.
   safe=(distance-u)/(2*M)
   bridge=safe<D('.01')
   if bridge:h=min(h,D('.01'))
   else:h=min(h,safe)
   assert h>0 and t+h>t
   reach=u+M*h
   assert box['lo']<y[0]-reach and y[0]+reach<box['hi'],'thermal domain exhausted'
   assert all(abs(y[i+1])+reach*SCALE[i+1]<box['Y'][i] for i in range(2)),'chemical domain exhausted'
   q=h/c['r'];assert q<1
   local=L*M*h*h*(2*L*h).exp() if bridge else c['R']*q**5/(1-q)
   previous=y;y=step(model.coupled_rhs,y,h)
   u=(L*h).exp()*u+local+D('1e-80');t+=h;steps+=1;bridges+=int(bridge)
   newcell=bisect.bisect_right(model.knots,y[0])-1
   if newcell!=cell:
    assert bridge,'unbridged crossing'
    knot=max(cell,newcell)
    # A conservative event interval from positive transversal speed is supplied later.
    crossings.append(dict(knot=knot,left_s=t-h,right_s=t,left_x=previous[0],right_x=y[0],U=u))
   assert u<D('1e-10'),'enclosure exhausted'
  records.append(dict(t_s=t,x=y[0],e=y[1],mu=y[2],U=u))
 return dict(records=records,crossings=crossings,steps=steps,bridges=bridges,box=box)
