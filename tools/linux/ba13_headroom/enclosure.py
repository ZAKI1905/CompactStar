"""Candidate deterministic propagation; no fitter import, solver or threshold repair.

Algebraic bounds are conditional on supplied error sets. None of the functions
turn a quadrature spread or a sampled maximum into an established input bound.
"""
import math
import numpy as np

FRACTIONS=dict(eos_equilibrium=.08,center=.04,surface_events=.06,integration=.10,
 implicit_stages=.02,floating_point=.02,common_B=.10,production_representation=.23,
 differentiation_inverse=.10,normalization=.05,nonlinear_fit=.10,between_samples=.10)
TYPES={'analytic_bound','validated_numerical_estimate','conditional_estimate','unresolved'}
L=np.array([[-1.,-1.],[1.,0.],[0.,1.]])

def budget(headroom):
 if not math.isfinite(headroom) or headroom<0:raise ValueError('negative/nonfinite headroom')
 assert abs(sum(FRACTIONS.values())-1)<1e-14
 return {k:headroom*v for k,v in FRACTIONS.items()}

def paired_affine(coeff_i,coeff_0,independent_i=0.,independent_0=0.,same_point=False):
 """One bounded theta per shared source; independent means arbitrary, not random."""
 if same_point:return 0.
 if min(independent_i,independent_0)<0:raise ValueError('negative uncertainty')
 return float(np.sum(np.abs(np.asarray(coeff_i)-np.asarray(coeff_0))))+independent_i+independent_0

def normalized_drift(observed,paired,scale,reference_error):
 if not all(math.isfinite(v) for v in (observed,paired,scale,reference_error)):raise ValueError('nonfinite')
 if min(observed,paired,reference_error)<0 or scale<=reference_error:raise ValueError('uncertain normalization')
 return (observed+paired)/(scale-reference_error)

def root_error(residual,functional_error,min_slope):
 if not all(math.isfinite(v) for v in (residual,functional_error,min_slope)) or functional_error<0 or min_slope<=0:raise ValueError('invalid root enclosure')
 return (abs(residual)+functional_error)/min_slope

def inverse_majorant(g,e):
 g=np.asarray(g,float);e=np.asarray(e,float)
 if not np.all(np.isfinite(g)) or not np.all(np.isfinite(e)) or np.any(e<0):raise ValueError('invalid matrix')
 a=np.linalg.inv(g);k=np.abs(a)@e
 rho=float(max(abs(np.linalg.eigvals(k))))
 if rho>=1:raise ValueError('inverse not enclosed')
 linear=k@np.abs(a);total=np.linalg.solve(np.eye(len(g))-k,linear)
 if np.any(total < -1e-14*np.max(np.abs(total))):raise ValueError('negative majorant')
 remainder=np.maximum(total-linear,0.)
 return a,total,remainder,rho

def paired_inverse(g0,gi,e_shared,e0=None,ei=None):
 """Shared entrywise DeltaG, separate changing Jacobians, full Neumann tails."""
 es=np.asarray(e_shared,float);e0=np.zeros_like(es) if e0 is None else np.asarray(e0,float);ei=np.zeros_like(es) if ei is None else np.asarray(ei,float)
 a0,total0,rem0,rho0=inverse_majorant(g0,es+e0);ai,totali,remi,rhoi=inverse_majorant(gi,es+ei)
 linear=np.zeros((2,2));p0=L.T@a0;pi=L.T@ai;q0=a0@L;qi=ai@L
 for j in range(3):
  for k in range(3):linear+=np.abs(np.outer(pi[:,j],qi[k])-np.outer(p0[:,j],q0[k]))*es[j,k]
 linear+=np.abs(L).T@(np.abs(a0)@e0@np.abs(a0)+np.abs(ai)@ei@np.abs(ai))@np.abs(L)
 remainder=np.abs(L).T@(rem0+remi)@np.abs(L)
 return dict(linear=linear,remainder=remainder,total=linear+remainder,rho=max(rho0,rhoi))

def projection(x,errors,values):
 """Exact historical weight floor; scaled coordinates only improve arithmetic."""
 x=np.asarray(x,float);y=np.asarray(values,float);u=np.asarray(errors,float)
 if x.shape!=(11,) or y.shape!=(11,) or u.shape!=(11,):raise ValueError('outer half must contain 11 points')
 sigma=np.maximum(u,64*np.finfo(float).eps*np.maximum(1.,abs(y)))
 w=(min(sigma)/sigma)**2;z=(x-x.mean())/np.ptp(x);a=np.column_stack([np.ones(11),z])
 return a@np.linalg.solve(a.T@(w[:,None]*a),a.T*w)

def residual_error(h,shared,independent):
 r=np.eye(len(h))-h
 return np.sum(abs(r@np.asarray(shared)),axis=1)+abs(r)@np.asarray(independent)

def between_samples(second_derivative,spacing,jump=0.):
 if min(second_derivative,spacing,jump)<0 or not all(map(math.isfinite,(second_derivative,spacing,jump))):raise ValueError('invalid continuity data')
 return second_derivative*spacing**2/8+jump

def qualify(historical,terms,*,exact=False,accepted_authorities=None,frozen_historical=None):
 """No zero fill; no covariance accepted as physical source; no threshold update."""
 accepted_authorities=accepted_authorities or {}
 if frozen_historical is not None and historical!=frozen_historical:raise ValueError('historical threshold/drift identity changed')
 if not all(math.isfinite(historical[k]) for k in ('drift_bound','threshold')) or historical['drift_bound']<0 or historical['threshold']<=0:raise ValueError('invalid historical decision inputs')
 if set(terms)!=set(FRACTIONS):raise ValueError('incomplete contribution schema')
 for term in terms.values():
  if term['evidence_type'] not in TYPES:raise ValueError('unknown evidence type')
  if term.get('source')=='fit_covariance':raise ValueError('covariance is not physical uncertainty')
  v=term['value']
  if v is not None and (not math.isfinite(v) or v<0):raise ValueError('invalid uncertainty')
 reasons={}
 for k,v in terms.items():
  if v['value'] is None or v['evidence_type']=='unresolved' or v.get('missing_estimate'):
   reasons[k]=v.get('missing_estimate') or 'missing/unqualified evidence'
  elif v['evidence_type'] in ('conditional_estimate','validated_numerical_estimate'):
   field='condition_id' if v['evidence_type']=='conditional_estimate' else 'standard_id'
   name=v.get(field);authority=v.get('authority')
   # The caller supplies an explicit registry of cited, externally accepted
   # decisions. A term cannot create its own authority by carrying a string.
   if not name or not authority or accepted_authorities.get(name)!=authority:
    reasons[k]='explicit acceptance required for '+field+': '+str(name)
 missing=list(reasons)
 partial=historical['drift_bound']+sum(t['value'] or 0 for t in terms.values())
 conditional=[k for k,t in terms.items() if t['evidence_type']=='conditional_estimate']
 if exact and (historical['drift_bound']!=0 or any(t['value'] not in (0,None) for t in terms.values())):
  return dict(status='FAIL',proposed_total_utilization=None,partial_utilization=partial/historical['threshold'],missing=missing,conditional=conditional,refusal_reasons={'exact':'exact identity violated'})
 if missing:return dict(status='UNQUALIFIED',proposed_total_utilization=None,partial_utilization=partial/historical['threshold'],missing=missing,conditional=conditional,refusal_reasons=reasons)
 if exact and any(t['evidence_type']!='analytic_bound' for t in terms.values()):
  return dict(status='UNQUALIFIED',proposed_total_utilization=None,partial_utilization=partial/historical['threshold'],missing=['exact'],conditional=conditional,refusal_reasons={'exact':'numerical evidence cannot establish exact identity'})
 return dict(status='PASS' if partial<=historical['threshold'] else 'FAIL',proposed_total_utilization=partial/historical['threshold'],partial_utilization=partial/historical['threshold'],missing=[],conditional=conditional,refusal_reasons={})
