#!/usr/bin/env python3
"""Independent short DOP853 integration and separately scaled chemical ledger."""
import argparse,csv,json,os
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp,quad
from model import Model

p=argparse.ArgumentParser();p.add_argument('out',type=Path);a=p.parse_args()
assert os.environ.get('SLURM_JOB_ID');m=Model()
directory=a.out/'short'
native={r['quantity']:float(r['value']) for r in csv.DictReader((directory/'native-inputs.tsv').open(),delimiter='\t')}
endpoints=list(csv.DictReader((directory/'endpoints.csv').open()))
accepted=list(csv.DictReader((directory/'accepted.csv').open()))
m.I=np.array([native['I0'],native['I1']])
checks={};results=[]
for role in range(4):
 spin=role in [1,3];bnv=role in [2,3]
 D=native['Gamma_yr']*m.B/m.year if bnv else 0
 # Independent expanded-dimensional polynomial implementation.
 def powers(time,x,eta):
  T=1e8*np.exp(x);R=[];ex=0.
  for L,e in zip(m.L,eta):
   R.append(L/m.kberg*sum(c*e**k*T**(7-k)/(11513*np.pi**(k+1)*m.kb**k) for c,k in [(14680,1),(7560,3),(840,5),(24,7)]))
   ex+=L*sum(c*e**k*T**(8-k)/(11513*np.pi**k*m.kb**k) for c,k in [(22020,2),(5670,4),(420,6),(9,8)])
  R=np.array(R);h=m.erg*np.dot(eta,R);eq=sum(m.L)*T**8;lg=m.photon(T)
  period=np.sqrt(native['P0_s']**2+2*native['P0_s']*native['Pdot0']*time)
  omega=2*np.pi/period;omegadot=-2*np.pi*native['P0_s']*native['Pdot0']/period**3
  ref=-2*m.I*omega*omegadot if spin else np.zeros(2)
  sigma=D*m.t;AB=-m.erg*np.dot(eta,sigma);AO=-m.erg*np.dot(eta,ref)
  return T,R,h,ex,eq,lg,ref,sigma,AB,AO
 # State x, eta/(1e-9 MeV), integrated Pnet/1e41 and chemical works/1e30.
 def rhs(time,y):
  T,R,h,ex,eq,lg,ref,sigma,AB,AO=powers(time,y[0],y[1:3]*1e-9)
  net=h-ex-eq-lg
  return np.r_[net/(T*m.capacity(T)),-m.Z@(R+sigma+ref)/1e-9,
               net/1e41,AB/1e30,AO/1e30,h/1e30]
 refs=[]
 for tol in [2e-10,2e-12]:
  sol=solve_ivp(rhs,(0,m.year),np.zeros(7),method='DOP853',rtol=tol,atol=tol*1e-5,dense_output=True)
  assert sol.success;refs.append(sol)
 y=refs[-1].y[:,-1];eta=y[1:3]*1e-9
 n=next(r for r in endpoints if int(r['case'])==role and int(r['tier'])==1)
 native_y=np.array([float(n['x']),float(n['eta_e_MeV']),float(n['eta_mu_MeV'])])
 delta=np.r_[y[0],eta]-native_y
 energy=.5*m.erg*eta@m.Zinv@eta;integral=(y[4]+y[5]-y[6])*1e30
 chem_norm=max(abs(y[4]),abs(y[5]),abs(y[6]),1e-30)*1e30
 chemical_error=abs(energy-integral)/chem_norm
 T=1e8*np.exp(y[0]);dU=quad(lambda x:m.capacity(1e8*np.exp(x))*1e8*np.exp(x),0,y[0],epsabs=1e20,epsrel=1e-11)[0]
 thermal_error=abs(dU-y[3]*1e41)/max(abs(dU),1.)
 # Independent pointwise checks at production accepted endpoints.
 max_rhs=0.;max_chain=0.;max_R18=0.;max_escape=0.;max_power=0.;near_equilibrium_net_sign=True
 for r in accepted:
  if int(r['case'])!=role:continue
  e=np.array([float(r['eta_e_MeV']),float(r['eta_mu_MeV'])]);xx=float(r['x']);t=float(r['t_s'])
  TT,R,h,ex,eq,lg,ref,sigma,AB,AO=powers(t,xx,e)
  # Compare each small power on its own scale: total cooling cannot resolve it.
  for key,expected in [('LH',h),('DeltaLnu',ex),('Lnu_eq',eq),('Lgamma',lg),('AB',AB),('Aspin',AO),('Cstar',m.capacity(TT))]:
   max_power=max(max_power,abs(float(r[key])-expected)/max(abs(expected),1e-60))
  if np.any(e):near_equilibrium_net_sign=near_equilibrium_net_sign and h-ex<0
  ede=-m.Z@(R+ref+sigma);actual=np.array([float(r['eta_dot_e']),float(r['eta_dot_mu'])])
  max_rhs=max(max_rhs,float(max(abs(ede-actual))/max(max(abs(ede)),1e-30)))
  cs=max(abs(AB)+abs(AO)+abs(h),1.)
  max_chain=max(max_chain,abs(float(r['chain_residual']))/cs)
  max_R18=max(max_R18,abs(float(r['R18']))/cs)
  expected_escape=m.erg*D*(m.mu+m.t@e)
  max_escape=max(max_escape,abs(float(r['escape'])-expected_escape)/max(abs(expected_escape),1.))
 results.append(dict(case=role,native_minus_independent_state=(-delta).tolist(),
                     independent_tolerance_difference=(refs[0].y[:,-1]-refs[1].y[:,-1]).tolist(),
                     chemical_integrated_residual_relative=chemical_error,
                     thermal_integrated_residual_relative=thermal_error,
                     max_independent_chemical_rhs_relative=max_rhs,
                     max_chain_relative=max_chain,max_R18_relative=max_R18,max_escape_relative=max_escape,
                     max_independent_power_relative=max_power,near_equilibrium_net_weak_negative=bool(near_equilibrium_net_sign),
                     chemical_energy_erg=energy,integrated_AB_erg=y[4]*1e30,integrated_Aspin_erg=y[5]*1e30,integrated_LH_erg=y[6]*1e30))
 checks[f'case_{role}']=bool(abs(delta[0])<1e-12 and max(abs(delta[1:]))<1e-19 and chemical_error<1e-8 and thermal_error<1e-8 and max_rhs<1e-10 and max_chain<1e-10 and max_R18<1e-10 and max_escape<1e-13 and max_power<1e-10 and near_equilibrium_net_sign)
checks['all_short_validity_guards']=all(float(r['frozen_utilization'])<=1 and 1e6<=float(r['Tinf_K'])<=1e9 and float(r['depletion'])<=1e-6 for r in accepted)
checks['zero_direct_everywhere']=all(float(r['Pdir'])==0 for r in accepted)
result=dict(status='PASS' if all(checks.values()) else 'FAIL',checks=checks,results=results,job=os.environ['SLURM_JOB_ID'],
            scope='Four one-year cases, independent DOP853 and dimensional polynomials. Chemical and thermal identities tested separately; no global rotating-star closure.')
(directory/'independent-check.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2));assert all(checks.values())
