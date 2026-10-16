#!/usr/bin/env python3
"""Separate frozen-relative qualification. Retains all historical fields verbatim."""
import csv,hashlib,json,math,pathlib,resource,sys,time
import numpy as np
from enclosure import *
S=pathlib.Path(__file__).resolve().parents[3];R=S.parents[1];O=R/'ba13-headroom';D=S/'docs/validation/ba13-headroom'
sys.path.insert(0,str(S/'tools/linux'))
from ba13_refit import load,NAMES
start=time.process_time()
def read(p):return json.loads(p.read_text())
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(name,v):(O/name).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
x,series=load(R/'ba13-recovery/run1/raw');cert=read(R/'ba13-recovery/run1/fit/certificate.json')
headroom=read(S/'docs/validation/x2-fixed-step/evidence/headroom.json')
hist_summary=read(S/'docs/validation/ba13-recovery/evidence/fits/run1/certificate-summary.json')
assert digest(S/'docs/validation/ba13-recovery/evidence/fits/run1/certificate-summary.json')==headroom['source_sha256']
profiles={i:read(O/f'profile-v2-{i}.json') for i in range(21)}
supp={r['index']:r for r in read(O/'supplement.json')['records']}
for d in profiles.values():
 for p,h in d['hash_inputs'].items():assert digest(pathlib.Path(p))==h
pilots={};physical={};estimates={};sens={}
def matrix4(v):return np.array([[v[0],0,0],[0,v[1],v[2]],[0,v[2],v[3]]])
def fields(v):
 p=v['primary'];z=L.T@np.linalg.solve(matrix4(v['G_packed']),L)
 K=v['K'];g=p[2]/p[1]**2/math.sqrt(1-2*p[2]/p[1]);tb=math.exp(v['Hb']-K)
 return dict(N_n=v['N'][0],N_e=v['N'][2],N_mu=v['N'][3],Ltilde_Me=p[4],Ltilde_Mmu=v['Lmu'],Cstar=p[5],radius=p[1],surface_gravity=g,envelope=(g*tb**2.42)**.25,Tsurface_inf=(g*tb**2.42)**.25*math.sqrt(1-2*p[2]/p[1]),mu_B=p[7],P0_average=v['P0'],P1_average=v['P1'],Z_row_npe=z[0],Z_row_npmu=z[1])
for i in [0,10,20]:
 a=read(O/f't{i}-n128.json');b=read(O/f't{i}-n512.json');f=fields(b);fa=fields(a)
 variants={k:fields(read(O/f't{i}-{k}.json')) for k in ['center-small','center-large','eos192']}
 widths={};slopes={}
 for h in ['0.0001','0.00005']:
  lo=read(O/f't{i}-h{h}--1.json');hi=read(O/f't{i}-h{h}-1.json');db=hi['primary'][0]-lo['primary'][0]
  widths[h]=(np.array(hi['N'])-lo['N'])/db
  slopes[h]=(hi['primary'][0]-lo['primary'][0])/(hi['primary'][8]-lo['primary'][8])
  if h=='0.00005':
   fl,fh=fields(lo),fields(hi);sens[i]={k:(np.asarray(fh[k])-fl[k])/db for k in f}
   bspan=db
 for j,n in [(0,'t_n'),(2,'t_e'),(3,'t_mu')]:f[n]=widths['0.00005'][j]
 physical[i]=f;estimates[i]={}
 derivative_changes={}
 for variant in ['center','eos']:
  lo=read(O/f't{i}-derivative-{variant}--1.json');hi=read(O/f't{i}-derivative-{variant}-1.json')
  derivative_changes[variant]=(np.array(hi['N'])-lo['N'])/(hi['primary'][0]-lo['primary'][0])-widths['0.00005']

 for n in f:
  if n.startswith('t_'):
   j={'t_n':0,'t_e':2,'t_mu':3}[n]
   # Conditional h^2 truncation plus independent deterministic endpoint error.
   num=2*abs(b['N'][j]-a['N'][j])/255
   estimates[i][n]=dict(integration=2*num/abs(bspan),stencil=abs(widths['0.00005'][j]-widths['0.0001'][j])*4/3,
    eos_equilibrium=2*abs(derivative_changes['eos'][j]),center=2*abs(derivative_changes['center'][j]))
  else:
   scale=np.maximum(abs(np.asarray(f[n])),1e-300)
   estimates[i][n]=dict(integration=float(np.max(2*abs(np.asarray(f[n])-fa[n])/255/scale)),
    eos_equilibrium=float(np.max(2*abs(np.asarray(f[n])-variants['eos192'][n])/scale)),
    center=float(max(np.max(2*abs(np.asarray(f[n])-variants[k][n])/scale) for k in ['center-small','center-large'])))
 defect=read(O/f't{i}-defect.json');defect512=read(O/f't{i}-defect512.json');fine=read(R/f'x2-fixed-step/t{i}-n4096.json')['values'];p=defect['primary']
 du,dv=defect['defect_transport_uv'];uv=[p[1]**2,p[2]/p[1]**3]
 pilots[str(i)]=dict(defect512=defect512,derivative_ablation_changes={k:v.tolist() for k,v in derivative_changes.items()},dB_dH=slopes,dB_dH_relative_spread=abs(slopes['0.00005']/slopes['0.0001']-1),
  stencil={k:v.tolist() for k,v in widths.items()},quantities={k:np.asarray(v).tolist() for k,v in f.items()},
  errors=estimates[i],root_computed_residual_relative=b['primary'][0]/cert['samples'][i]['B_solved_count']-1,
  B_functional_discretization_estimate=2*abs(a['primary'][0]-b['primary'][0])/255,
  defect_transport_R_relative=du/(2*uv[0]),defect_transport_M_relative=dv/uv[1]+1.5*du/uv[0],
  fine_minus_128_R_relative=fine[1]/p[1]-1,fine_minus_128_M_relative=fine[2]/p[2]-1,
  defect_absolute_R_estimate=defect['defect_absolute_transport_uv'][0]/(2*uv[0]),
  defect_absolute_M_estimate=defect['defect_absolute_transport_uv'][1]/uv[1]+1.5*defect['defect_absolute_transport_uv'][0]/uv[0],
  stage_transport_R_estimate=defect['stage_transport_uv'][0]/(2*uv[0]),
  defect_quadrature_spread=defect['defect_quadrature_spread'],
  defect_status='CONDITIONAL: numerical Jacobians and quadrature; no derivative-tube or intra-step transport remainder enclosure')
# Reconstruction errors are signed, shared functions of the nominal profile.
# The finite differences below assess differences of those errors, never hypot.
recon={};absolute={};norm={};prod_delta={};Berr={}
for i,d in profiles.items():
 Berr[i]=sum(d['N_nu'][:2])-d['star']['B_solved']
 recon[i]=dict(Cstar=d['C_nu'],Ltilde_Me=d['L_nu'],Ltilde_Mmu=d['Lmu_nu'],N_n=d['N_nu'][0],N_e=d['N_nu'][2],N_mu=d['N_nu'][3],mu_B=d['mu_inf'],P0_average=d['mu_inf'],Z_row_npe=np.array(d['Z_nu'])[0],Z_row_npmu=np.array(d['Z_nu'])[1],P1_average=supp[i]['integrals']['16']['P1'])
 # Real nominal output, with only the mathematically common cache factor removed.
 prod_delta[i]={}
 for n,v in recon[i].items():
  prod=d['C_decache'] if n=='Cstar' else np.asarray(series[n][0][i])
  if np.size(prod)==1:prod=float(np.ravel(prod)[0])
  prod_delta[i][n]=np.asarray(v)-prod
for i in range(21):
 for n in ['radius','surface_gravity']:
  references=[physical[t][n] for t in [0,10,20]]
  if n=='surface_gravity':
   # Convert geometric acceleration to historical g14 using the shared c constant.
   references=[v*(2.99792458e10)**2/1e5/1e14 for v in references]
  recon[i][n]=float(np.interp(i,[0,10,20],references))
  prod_delta[i][n]=np.array(recon[i][n])-float(series[n][0][i][0])
# Fixed-point error relative to reconstructed quantities and paired error difference.
for n in NAMES:
 values,errs,scales,kind=series[n];s=np.asarray(scales);v=np.asarray(values)
 norm[n]=s
 absolute[n]=np.zeros(21)
 if n in recon[0]:
  s0=profiles[0]['C_decache'] if n=='Cstar' else float(series[n][2][0])
  absolute[n]=np.array([float(np.max(abs(prod_delta[i][n])))/s0 for i in range(21)])
 # Preserve *signed* endpoint corrections, including the changing scale.
 if n in recon[0]:
  paired=np.array([float(np.max(abs(prod_delta[i][n]-prod_delta[0][n])))/s0 for i in range(21)])
 else:paired=None
 estimates[n]=dict(absolute=absolute[n].tolist(),paired=paired.tolist() if paired is not None else None)
rows=[];inner={};missing_counts={}
for n in NAMES:
 values,errors,scales,kind=series[n];values=np.array(values);errors=np.array(errors);scales=np.array(scales)
 raw=np.max(abs(values-values[0])/scales,axis=1)
 fits=cert['fits'][n];co=np.array([[f['intercept'],f['slope']] for f in fits])
 residual=np.max(abs(values-(co[:,0][None,:]+np.array(x)[:,None]*co[:,1][None,:]))/scales,axis=1)
 outer=max(residual[10:]);allres=max(residual)
 inner[n]=dict(outer_max=float(outer),inner_max=float(max(residual[:10])),all_sample_max=float(allres),residuals=residual.tolist())
 # Actual sampled curvature is a diagnostic and cannot be used as an upper bound.
 slopes=np.diff(values,axis=0)/np.diff(x)[:,None]
 curvature=float(np.max(abs(np.diff(slopes,axis=0)/(np.diff(x)[:-1,None]+np.diff(x)[1:,None])*2)/scales))
 for i in range(21):
  sample=cert['samples'][i];hr=next(r for r in headroom['rows'] if r['target']==i and r['quantity']==n)
  assert hr['threshold']==sample['threshold'][n] and hr['drift_bound']==sample['drift_bound'][n]
  history={k:sample[k][n] for k in ['drift_bound','threshold','utilization']};h=hr['headroom'];alloc=budget(h)
  term={k:dict(value=None,evidence_type='unresolved',source=None,assumptions=[],missing_estimate='') for k in FRACTIONS}
  def setterm(k,value,source,assumptions=(),typ='conditional_estimate'):
   term[k]=dict(value=float(value) if value is not None else None,evidence_type=typ,source=source,assumptions=list(assumptions),missing_estimate='')
  def gap(k,why):term[k]['missing_estimate']=why
  # Do not extrapolate a formal bound from three points. These are labelled
  # transfer estimates, useful for headroom selection but subject to review.
  if n in physical[0]:
   base=estimates[0][n];pilotvals=[estimates[t][n] for t in [0,10,20]]
   for k in ['integration','center','eos_equilibrium']:
    vals=[d[k] for d in pilotvals if d.get(k) is not None]
    if vals:
     measured=2*max(vals)
     setterm(k,measured if measured else None,f'X2 three-target {k}',('observed refinement/ablation represents remaining error; transfer over depletion interval',))
     term[k]['measured_sensitivity']=measured
     if k=='eos_equilibrium':gap(k,'finite 96/192-node comparison is not a uniform EOS quadrature and equilibrium inversion remainder estimate')
     if k=='center':gap(k,'finite startup ablation is not a center-series and initial-integral remainder estimate')
   if n.startswith('t_'):
    setterm('differentiation_inverse',2*max(v['stencil'] for v in pilotvals),'two centered fixed-H stencils; attained B denominator',('h^2 expansion continues; numerical endpoint error charged in integration',))
    gap('center','measured derivative-of-center ablation lacks an unseen series-remainder estimate')
    gap('eos_equilibrium','measured derivative-of-EOS ablation lacks a uniform residual/inversion remainder estimate')
   else:setterm('differentiation_inverse',0,'no differentiation/inverse in this scalar',typ='analytic_bound')
   # Arithmetic of positive functional accumulation only, not input accuracy.
   gamma=(200*512*2**-64)/(1-200*512*2**-64)
   setterm('floating_point',2*gamma,'operation-count allowance for GL4 arithmetic',('200 rounded operations per step equivalent; elementary-function and EOS input errors excluded and remain separate',))
   if n=='radius':
    setterm('implicit_stages',2*max(p['stage_transport_R_estimate'] for p in pilots.values()),'passive stage-residual transport',('numerical stage tangent map; nonlinear remainder not enclosed',))
   else:
    gap('implicit_stages','quantity-specific functional response to stage residual; R proxy deliberately not used')
   gap('integration','conditional refinement estimate is checked, but defect transport has no enclosed intra-step/nonlinear remainder')
   gap('floating_point','operation-count allowance covers accumulation only; validated elementary-function/input roundoff propagation missing')
  else:
   for k in ['integration','center','eos_equilibrium','implicit_stages','floating_point']:
    gap(k,'profile/support or direct-average component needs its own transport from the structure/EOS error set; no scalar proxy adopted')
  # Independent fixed-profile representation differences are measured estimates.
  if estimates[n]['paired'] is not None:
   setterm('production_representation',max(estimates[n]['paired'][:i+1]),'signed all-target nominal versus onset-aware reconstructed quantities',('reconstruction background remains uncertain; paired corrections are not independent random noise',),'conditional_estimate')
   u0=estimates[n]['absolute'][0]
   if u0<1:setterm('normalization',(max(raw[:i+1])+max(estimates[n]['paired'][:i+1]))*u0/(1-u0),'exact ratio denominator propagation of observed reference correction',('reference correction uncertainty remains conditional on reconstruction accuracy',))
  elif n.startswith('t_'):
   # Comparison at validation targets, not substitution of oracle into nominal fit.
   differences=[abs(float(physical[t][n])-float(values[t,0])) for t in [0,10,20]]
   setterm('production_representation',2*max(differences),'nominal tangent versus independent GL4 number stencils',('three-target transfer; remaining derivative-background error not enclosed',))
   setterm('normalization',0,'absolute tangent norm; S=1',typ='analytic_bound')
  elif n=='metric_structure':
   delta=np.array([np.array(supp[j]['mapped_metric'])-values[j] for j in range(21)])
   setterm('production_representation',float(np.max(abs(delta[:i+1]-delta[0])/scales)),'all-target metric remapping from pinned edge to reconstructed neutron edge',('fixed-profile interpolation only; common background error remains',),'conditional_estimate')
   u0=float(np.max(abs(delta[0])/scales));setterm('normalization',max(raw[:i+1])*u0/(1-u0),'scaled-component denominator propagation',('reconstruction normalization error remains',))
  elif n=='mu_n_profile':
   delta=np.array([np.full(len(scales),profiles[j]['mu_inf'])-values[j] for j in range(21)])
   setterm('production_representation',float(np.max(abs(delta[:i+1]-delta[0])/scales)),'all nominal profile nodes versus constant redshifted chemical potential',('constant-potential reconstruction must share the bounded physical background',),'conditional_estimate')
   u0=float(np.max(abs(delta[0])/scales));setterm('normalization',max(raw[:i+1])*u0/(1-u0),'componentwise denominator propagation',('reconstruction reference uncertainty remains',))
  elif n in ['envelope','Tsurface_inf']:
   shifts=np.array([profiles[j]['Ts_offset'] for j in range(21)])
   setterm('production_representation',float(max(abs(shifts[:i+1]-shifts[0]))),'all-target physical boundary versus selected production node',('boundary interpolation and metric background estimates remain',),'conditional_estimate')
   setterm('normalization',max(raw[:i+1])*abs(shifts[0])/(1-abs(shifts[0])),'reference boundary offset denominator',('remaining background error',))
  else:
   gap('production_representation',{'metric_structure':'full scaled-node metric error, including physical-edge remapping and interpolation remainder','P1_average':'direct P1 fixed-background quadrature validated at three targets; remaining per-target event-measure/background error','radius':'production finite-cut surface error versus X2 at three targets; transfer/error enclosure missing','surface_gravity':'joint M,R surface propagation and nominal fixed-B alignment','species_support':'continuous exact support topology/no-crossing enclosure'}.get(n,'quantity-specific production discrepancy'))
   gap('normalization','uncertainty of the original B0 component scales under the same joint background perturbation')
  if n in sens[0]:
   correction=[]
   for j in range(21):
    q=np.array([np.interp(j,[0,10,20],[np.ravel(sens[t][n])[k] for t in [0,10,20]]) for k in range(np.size(sens[0][n]))])
    sc=float(np.max(abs(np.asarray(physical[0][n]))));correction.append(q*Berr[j]/sc)
   paired=max(float(np.max(abs(correction[j]-correction[0]))) for j in range(i+1))
   setterm('common_B',paired,'independent reconstructed B functional, not computed root residual',('local dQ/dB sensitivity interpolated from three targets; functional reconstruction uncertainty and derivative enclosure remain',))
  else:gap('common_B','quantity-specific dQ/dB or d(t)/dB with functional-error and achieved-coordinate uncertainty')
  if n in sens[0]:gap('common_B','uncertainty of reconstructed B functional and local sensitivity enclosure beyond measured secants')
  if n.startswith('Z_'):
   row=0 if n=='Z_row_npe' else 1;zs=max(abs(np.array(profiles[0]['Z_fixed_profile'])[row]));pairmax=0
   for j in range(i+1):
    # Shared inherited B0 background modeled ONCE; separate current centers/tails.
    e0=np.array(profiles[0]['center_inherited'])+profiles[0]['tail_error'];ei=np.array(profiles[j]['center_inherited'])+profiles[j]['tail_error']
    bound=paired_inverse(profiles[0]['G_fixed_profile'],profiles[j]['G_fixed_profile'],profiles[0]['background_inherited'],e0,ei)
    pairmax=max(pairmax,float(max(bound['total'][row]))/zs if j else 0.)
   setterm('differentiation_inverse',pairmax,'paired entrywise Neumann inverse majorant, changing Jacobian plus both remainders',('inherited background is exactly shared; per-target variation of that background is unresolved',))
   gap('differentiation_inverse','depleted-target change of the inherited background error; measured G differences do not bound this variation')
   setterm('surface_events',0,'tail enters inverse majorant only; not double counted',typ='analytic_bound')
  elif n!='species_support':
   setterm('surface_events',0,'finite-pressure surface is the defined quantity, not vacuum',('event location uncertainty belongs to integration/production representation; no independent missing event error may be inferred zero',))
  # Apply the fixed historical projection to the SAME signed paired error.
  residual_input=0.
  if n in recon[0]:
   correction=np.array([np.atleast_1d(prod_delta[j][n]) for j in range(21)])
   for component in range(values.shape[1]):
    col=0 if correction.shape[1]==1 else component
    correction_col=correction[:,col].copy()
    if n=='Cstar':correction_col*=values[0,component]/profiles[0]['C_decache']
    H=projection(np.array(x)[10:],errors[10:,component],values[10:,component])
    residual_input=max(residual_input,float(max(abs((np.eye(11)-H)@correction_col[10:])))/scales[component])
  setterm('nonlinear_fit',max(0.,allres-outer)+residual_input,'all 21 residuals plus (I-H)e for the same signed correction; no covariance added',('fixed design; uncertain achieved-B design and nonlinear response remain',),'conditional_estimate')
  if n!='P2_average':gap('nonlinear_fit','residual perturbation from uncertain achieved-B design, residual-input remainder, and continuous nonlinear fit remainder')
  gap('between_samples','no justified modulus for the paired discrepancy between achieved-B samples; production support/envelope node jumps require separate coverage')
  if n in ['envelope','Tsurface_inf','Ltilde_Me','metric_structure']:
   term['between_samples']['observed_cell_jump_scale']=max(profiles[j]['support_cell_relative'] if n=='Ltilde_Me' else .605*profiles[j]['envelope_cell_logT_jump'] for j in range(21)) if n!='metric_structure' else None
  if n=='species_support':
   assert np.all(values==0)
   for k in term:setterm(k,0,'exact topology agrees on all 21 original samples',typ='analytic_bound')
   term['between_samples']=dict(value=None,evidence_type='unresolved',source='exact support check',assumptions=[],missing_estimate='continuous support/no-crossing enclosure; 21 equal topologies alone do not prove no crossing between them')
  if n not in ['P2_average','species_support'] and not n.startswith('t_'):
   gap('normalization','B0 denominator uncertainty from the remaining shared structure/EOS/common-B error, beyond measured reconstruction correction')
  if n=='P2_average':
   assert np.all(values==0) and np.all(errors==0)
   for k in term:setterm(k,0,'exact P2=0 identity in original source and all samples',typ='analytic_bound')
  if n not in ['P2_average','species_support'] and term['production_representation']['value'] is not None:
   gap('production_representation','signed observed reconstruction discrepancy has no enclosed common-background/representation remainder')
  for k,v in term.items():
   if v['value'] is None and not v['missing_estimate']:v['missing_estimate']='quantity-specific '+k+' estimate not established'
  decision=qualify(history,term)
  row=dict(index=i,quantity=n,historical=history,historical_raw_drift=float(raw[i]),headroom=h,norm_kind=kind,component_scales=scales.tolist() if len(scales)<3 else dict(count=len(scales),min=float(min(scales)),max=float(max(scales))),
   allocated_budget=alloc,contributions=term,correlation='one signed correction per shared nominal/reconstruction error; no hypot; inverse uses distinct Jacobians',
   normalization='original B0 scales; additive denominator correction; t uses absolute S=1; C cache cancels only same knot',
   evidence_types=sorted(set(t['evidence_type'] for t in term.values())),
   continuous_coverage=dict(interval=[float(x[i]),0.],inner_half_evaluated=True,sampled_curvature=curvature,sufficient_curvature_if_no_jumps=8*alloc['between_samples']/max(abs(np.diff(x)))**2,assumption_validated=False),
   **decision)
  rows.append(row)
assert len(rows)==462 and len({(r['index'],r['quantity']) for r in rows})==462
from collections import Counter
summary=dict(schema_id='compactstar.ba13.headroom-uncertainty-assessment.v1',candidate=True,owner_ratified=False,entry_commit='89b629dd9961c29e127badbd723455d00dc38553',recommendation='B',counts=dict(Counter(r['status'] for r in rows)),
 fitter_sha256=digest(S/'tests/bnv/compare_candidate.py'),historical_certificate_sha256=digest(R/'ba13-recovery/run1/fit/certificate.json'),
 required_interpretation='ADR-0021 Option A and interpretation (2), proposed, right thresholds frozen',
 gaps=['continuous paired-error modulus and discrete production-node jump coverage','quantity-dependent nonlinear propagation remainder for numerical defect transport','remainder of measured tangent EOS/center perturbations and attained-B coordinate errors','variation of inherited chemical background at depleted targets'],
 all_profiles_complete=len(profiles),new_whole_sequence_solves=0,headroom_rows=462,
 cpu_seconds=time.process_time()-start,peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,source_sha256=digest(pathlib.Path(__file__)))
write('pilots.json',pilots);write('paired-discrepancies.json',estimates);write('inner-half.json',inner);write('matrix.json',dict(**summary,rows=rows));write('summary.json',summary)
with (O/'matrix.tsv').open('w') as f:
 cols=['index','quantity','norm_kind','historical_drift','threshold','headroom']+list(FRACTIONS)+['partial_utilization','proposed_total_utilization','status','missing']
 w=csv.DictWriter(f,fieldnames=cols,delimiter='\t');w.writeheader()
 for r in rows:w.writerow(dict(index=r['index'],quantity=r['quantity'],norm_kind=r['norm_kind'],historical_drift=r['historical']['drift_bound'],threshold=r['historical']['threshold'],headroom=r['headroom'],**{k:v['value'] if v['value'] is not None else 'UNRESOLVED' for k,v in r['contributions'].items()},partial_utilization=r['partial_utilization'],proposed_total_utilization=r['proposed_total_utilization'],status=r['status'],missing='; '.join(k+': '+r['contributions'][k]['missing_estimate'] for k in r['missing'])))
print(json.dumps(summary,indent=2))
