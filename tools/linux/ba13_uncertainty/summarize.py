#!/usr/bin/env python3
"""Assemble separate uncertainty evidence; never rewrite historical certificates."""
import csv,json,sys,math,re,hashlib
from pathlib import Path
import numpy as np
S=Path(__file__).resolve().parents[3];R=Path(__import__('os').environ.get('BA13_STUDY_OUT',str(S.parents[1]/'ba13-uncertainty')));OLD=S.parents[1]/'ba13-recovery'
sys.path.insert(0,str(S/'tools/linux'));from ba13_refit import load,NAMES
_,series=load(OLD/'run1/raw');cert=json.loads((OLD/'run1/fit/certificate.json').read_text())
def read(n):return json.loads((R/n).read_text())
def table(p):return list(csv.DictReader(p.open(),delimiter='\t'))
th=read('thermal-diagnostic.json');native=[v for v in th['records'] if v['platform']=='native'];mac=[v for v in th['records'] if v['platform']=='mac']
thermal=dict(exceedances=sum(abs(x['Cstar']-y['Cstar'])>x['original_uncertainty']+y['original_uncertainty'] for x,y in zip(native,mac)),
 max_platform_relative=max(abs(x['Cstar']-y['Cstar'])/x['Cstar'] for x,y in zip(native,mac)),
 max_independent_nodal_platform_relative=max(abs(x['analytic_nodal_trapezoid']-y['analytic_nodal_trapezoid'])/x['analytic_nodal_trapezoid'] for x,y in zip(native,mac)),
 max_cache_relative_bias=max(v['cache_relative_bias'] for v in native),max_reconstruction_relative=max(abs(v['producer_minus_prediction_relative']) for v in th['records']),
 fixed_term_summation_gamma=max(v['sequential_sum_gamma'] for v in native),high_precision_checks=len(th['high_precision']))
hold=read('holdout-validation.json');integrals=read('integral-diagnostic.json');profiles=read('profile-diagnostic.json');chem=read('chemical-diagnostic.json')['records'];units=read('fit-units.json')['quantities']
integ=[]
for i in (0,10,20):
 rows=[v for v in integrals if v['index']==i];base=next(v for v in rows if v['case']=='nominal' and v['order']==32 and not v['bisect']);coef=table(OLD/'run1/raw/coefficients.tsv')[i]
 nominal=table(R/'runs'/f'nominal-{i}'/'thermal.tsv')
 for key in ('Cstar_slope','nodal_density_Cstar_slope','mu_B','P0','P1','Ltilde_Me','Ltilde_Mmu'):
  vals=[v[key] for v in rows if v['case']=='nominal'];fine=next(v[key] for v in rows if v['case']=='r320000' and v['order']==32)
  production_key={'mu_B':'mu_B_inf','P0':'P0_inf','P1':'P1_inf'}.get(key,key)
  prod=float(coef[production_key]) if production_key in coef else None
  if key=='Cstar_slope':prod=float(nominal[0]['Cstar'])/float(nominal[0]['Tinf_K'])/native[i*13]['cache_factor']
  integ.append(dict(index=i,quantity=key,independent_fixed_profile=base[key],nominal_order_spread_relative=(max(vals)-min(vals))/abs(base[key]),
   radial_nominal_to_320k_relative=(fine-base[key])/base[key],original_to_independent_relative=(prod-base[key])/base[key] if prod is not None else None))
resources=[]
for p in sorted((R/'runs').glob('*.receipt.json')):
 v=json.loads(p.read_text());resources.append(dict(job=v['job'],case=v['case'],index=v['index'],cpu_seconds=v['cpu_seconds'],wall_seconds=v['wall_seconds'],max_rss_kib=v['max_rss_kib'],status=v['status'],kind='numerical'))
for name in ('build.log','eos-audit-build.log','thermal_diagnostic.log','fit_diagnostic.log','chemical-diagnostic.log','integral-diagnostic.log','profile-diagnostic.log','profile-diagnostic-first.log','eos-ordering.log','p1-validation.log'):
 p=R/name
 if not p.exists():continue
 for line in p.read_text().splitlines():
  if line.startswith('RESOURCE_USAGE_JSON '):
   v=json.loads(line.split(' ',1)[1]);resources.append(dict(log=name,job={'build.log':'14587','eos-audit-build.log':'14651','thermal_diagnostic.log':'14591','fit_diagnostic.log':'14591','chemical-diagnostic.log':'14596','integral-diagnostic.log':'14641','profile-diagnostic-first.log':'14651','profile-diagnostic.log':'14654','eos-ordering.log':'14651','p1-validation.log':'14654'}[name],cpu_seconds=v['user_cpu_seconds']+v['system_cpu_seconds'],wall_seconds=v['elapsed_seconds'],max_rss_kib=v['max_rss_kib'],status=v['exit_code'],kind='build' if 'build' in name else 'numerical'))
rows=[]
chemical_names={'t_n','t_e','t_mu','Z_row_npe','Z_row_npmu'}
for i in range(21):
 for name in NAMES:
  v,e,scales,kind=series[name];original=max(abs(a)/s for a,s in zip(e[i],scales));partial=[]
  mapkey={'radius':'R','N_n':'N_n','N_e':'N_e','N_mu':'N_mu','surface_gravity':'g14','Cstar':'Cstar','envelope':'Tsurface','Tsurface_inf':'Tsurface_inf'}.get(name)
  for h in hold:
   if h['index']==i and h['quantity']==mapkey and h['heldout_status']=='PASS_CONDITIONAL_PREDICTION':partial.append(dict(axis=h['axis'],knot=h['knot'],conditional_nominal_error=h['nominal_error_prediction'],normalized=h['nominal_error_prediction']/abs(h['values'][2]),status='EMPIRICAL_AXIS_ONLY'))
  missing=['global structure/EOS enclosure','common-B sensitivity enclosure','correlated reference/normalization propagation']
  source='structure, EOS, center/surface, interpolation, arithmetic'
  evidence=['profile-diagnostic.json','holdout-validation.json'] if i in (0,5,10,15,20) else ['authenticated nominal profile; no new structural refinement for this target']
  if name in chemical_names:
   source+='; tangent stencil/truncation' if name.startswith('t_') else '; response quadrature/inverse, inherited background'
   missing+=['new per-target tangent/background refinement']
   evidence+=['chemical-diagnostic.json']
  if name in ('Cstar','envelope','Tsurface_inf','surface_gravity'):
   source+='; thermal source/cache/envelope propagation';evidence+=['thermal-diagnostic.json'];missing+=['radial/onset thermal interpolation bound']
  if name in ('Ltilde_Me','Ltilde_Mmu','mu_B','P0_average','P1_average','mu_n_profile'):
   source+='; source support, fixed-profile quadrature, direct averages';evidence+=['integral-diagnostic.json' if i in (0,10,20) else 'original raw coefficients'];missing+=['independent whole-profile/event remainder']
  status='INCONCLUSIVE' if i in (0,5,10,15,20) else 'UNBOUNDED'
  proposed=None;impact=None
  if name=='species_support':missing+=['quantitative support-edge enclosure under all admissible errors'];evidence+=['original exact topology check; tested profile topology 3,2,1'];status='INCONCLUSIVE'
  if name=='P2_average':
   assert all(z==0 for z in v[i]+e[i]);status='PASS';proposed=0.;impact=cert['samples'][i]['utilization'][name];missing=[];source='exact analytic P2=0';evidence=['exact original P2 fields; recovered source formula'];partial=[]
  rows.append(dict(index=i,quantity=name,numerical_uncertainty_sources=source,independent_evidence=evidence,
   conservative_bound=proposed,missing_bounds=missing,correlation_treatment='shared signed structure/EOS perturbations; no independent-noise assumption; common cache factor cancels only in proved ratios',
   original_reported_uncertainty_normalized_max=original,original_reported_uncertainty_absolute_min=min(e[i]),original_reported_uncertainty_absolute_max=max(e[i]),components=len(scales),
   proposed_qualified_uncertainty=proposed,partial_evidence=partial,original_utilization=cert['samples'][i]['utilization'][name],qualified_utilization=impact,status=status))
from collections import Counter
summary=dict(disposition='B: specific additional bounded evidence required',qualified=False,thermal=thermal,
 holdout_counts=dict(Counter(v['heldout_status'] for v in hold)),matrix_counts=dict(Counter(v['status'] for v in rows)),
 completed_profile_cases=len(profiles),topology_checks_pass=all(v['topology_matches_nominal'] for v in profiles),
 max_target_residual_fraction=max(v['target_residual_fraction'] for v in profiles),max_bracket_fraction=max(v['bracket_fraction'] for v in profiles),
 center_estimate_ratio_range=[min(v['max_center_ratio'] for v in chem),max(v['max_center_ratio'] for v in chem)],
 integral_comparisons=integ,fit_units={n:{k:v for k,v in u.items() if k!='measurements'} for n,u in units.items()},
 cpu_seconds_by_kind={k:sum(v['cpu_seconds'] for v in resources if v['kind']==k) for k in ('numerical','build')},
 max_measured_rss_kib=max(v['max_rss_kib'] for v in resources),max_single_case_wall_seconds=max(v['wall_seconds'] for v in resources if 'case' in v))
for name,data in [('summary.json',summary),('disposition-matrix.json',rows),('resources.json',resources)]:
 tmp=R/(name+'.partial');tmp.write_text(json.dumps(data,indent=2)+'\n');tmp.replace(R/name)
with (R/'disposition-matrix.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['index','quantity','status','original_uncertainty_normalized_max','proposed_qualified_uncertainty','original_utilization','qualified_utilization','missing_bounds'])
 for v in rows:w.writerow([v['index'],v['quantity'],v['status'],v['original_reported_uncertainty_normalized_max'],v['proposed_qualified_uncertainty'],v['original_utilization'],v['qualified_utilization'],'; '.join(v['missing_bounds'])])
print(json.dumps({k:v for k,v in summary.items() if k not in ('integral_comparisons','fit_units')},indent=2))
