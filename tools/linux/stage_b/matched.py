#!/usr/bin/env python3
"""Matched budgets use sums, with no source/control uncertainty cancellation."""
import json,pathlib,sys
from decimal import Decimal as D,localcontext
from piecewise import ROOT  # Establish the retained Stage A helper import path.
from run import save
out=pathlib.Path(sys.argv[1]);accuracy=json.loads((out/'accuracy.json').read_text());ledger=json.loads((out/'ledger-assessment.json').read_text())
a={(r['role'],r['configuration']):r for r in accuracy['configurations']};l={(r['role'],r['configuration']):r for r in ledger['configurations']};results=[]
with localcontext() as ctx:
 ctx.prec=110
 for tier in range(3):
  for cap in range(4):
   name=f'tier-{tier}-cap-{cap}';source=a['source',name];control=a['control',name];points={}
   for kind in ('states','powers'):
    records=[]
    for s,c in zip(source['endpoint_budget'][kind],control['endpoint_budget'][kind]):
     assert s['channel']==c['channel'];value=D(s['value'])-D(c['value']);reference=D(s['reference'])-D(c['reference']);u=D(s['uncertainty'])+D(c['uncertainty']);budget=D(s['budget'])+D(c['budget'])
     records.append(dict(channel=s['channel'],matched_native=value,matched_reference=reference,summed_uncertainty=u,summed_budget=budget,budget_ratio=(abs(value-reference)+u)/budget,signal='RESOLVED' if reference and u<=D('.1')*abs(reference) else 'SIGN_UNRESOLVED'))
    points[kind]=records
   ts=D(source['endpoint_budget']['T_native']);tc=D(control['endpoint_budget']['T_native']);rs=D(source['endpoint_budget']['T_reference']);rc=D(control['endpoint_budget']['T_reference'])
   ut=D(0)
   for r in (source,control):
    e=r['endpoint_budget'];x=e['states'][0];T=D('1e8')*D(x['reference']).exp();eval_error=abs(T-D(e['T_reference']))
    ut+=abs(D(e['T_native'])-D(e['T_reference']))+T*(D(x['uncertainty']).exp()-1)+eval_error
   ls=l['source',name];lc=l['control',name];N=D(ls['N'])+D(lc['N'])
   uncertainty=sum(D(r[k]) for r in (ls,lc) for k in ('U_main','U_reconstruction_including_stored_events','U_event_extra','U_quad','U_count'))
   residual=D(ls['R20_original'])-D(lc['R20_original']);ratio=(abs(residual)+uncertainty)/N
   qualified=all(all(r['gates'].values()) and not r['missing'] for r in (source,control)) and all(all(r['gates'].values()) for r in (ls,lc))
   results.append(dict(configuration=name,status='CONDITIONAL_MATCHED_BUDGET_CANDIDATE' if qualified else 'PRELIMINARY_INCOMPLETE_OR_FAILED_COMPONENT',endpoint=points,temperature=dict(native_delta_K=ts-tc,reference_delta_K=rs-rc,total_error_bound_about_native_K=ut,signal='RESOLVED' if ut<=D('.1')*abs(rs-rc) else 'SIGN_UNRESOLVED'),R20=dict(residual_difference_erg=residual,summed_uncertainty_erg=uncertainty,summed_N_erg=N,ratio=ratio,limit=D('2e-4'))))
 save(out/'matched.json',dict(status='CANDIDATE',configurations=results,all_times_budget_argument='At every common scientific time, including the final endpoint, the matched discrepancy is bounded by the sum of the two individual discrepancies. Adding both positive uncertainties and both budgets proves the matched state/power inequality whenever both individual inequalities pass. The same triangle inequality applies to R18/R20 with summed gross normalizations. Each role separately verifies its own noncoincident adaptive endpoints. No uncertainty cancellation is used. Endpoint resolved-signal tests are explicit; no all-time chemical-power sign-resolution claim is made.'))
