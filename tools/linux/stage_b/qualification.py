#!/usr/bin/env python3
"""Collect completed gates without promoting verification to owner acceptance."""
import json,pathlib,sys
from decimal import Decimal as D
from piecewise import ROOT  # Establish the retained Stage A helper import path.
from run import save,digest

out=pathlib.Path(sys.argv[1]);inputs={};checks={}
def read(name):
 p=out/name
 inputs[name]=dict(path=str(p),sha256=digest(p))
 return json.loads(p.read_text())
expected={(role,f'tier-{tier}-cap-{cap}') for role in ('source','control') for tier in range(3) for cap in range(4)}
preflight=read('preflight.json')
checks['preflight_input_and_GSL_tableau_bindings']=preflight['status']=='PASS' and all(digest(pathlib.Path(p))==h for p,h in preflight['files'].items())
def matrix(name,records):
 keys=[(r['role'],r['configuration']) for r in records]
 checks[name+'_complete']=len(keys)==24 and set(keys)==expected

accuracy=read('accuracy.json');matrix('accuracy',accuracy['configurations'])
required={'power_refinement_contraction','accepted_powers','accepted_reference_power_allocation','scientific_states','reference_state_allocation','powers','reference_power_allocation','accepted_restart_states','accepted_reference_allocation','state_refinement_contraction'}
for r in accuracy['configurations']:
 key=r['role']+'/'+r['configuration']
 checks[key+'/accuracy']=set(r['gates'])==required and all(r['gates'].values()) and not r['missing'] and r['complete_arithmetic']
ledger=read('ledger-assessment.json');matrix('ledger',ledger['configurations'])
for r in ledger['configurations']:
 checks[r['role']+'/'+r['configuration']+'/ledger']=set(r.get('gates',{}))=={'R20','allocation','R18','contraction'} and all(r['gates'].values())
local=read('local-assessment.json');matrix('local_reconstruction',local['configurations'])
checks['all_local_reconstructions']=all(r['status']=='CONDITIONAL_PASS' and r['points']==1921 for r in local['configurations'])
coverage=read('coverage.json')
keys=[(r['role'],f"tier-{r['tier']}-cap-{r['cap']}") for r in coverage['configurations']]
checks['matrix_schedule_grid_restart_event_guards']=coverage['status']=='PASS_MATRIX_AND_NATIVE_GUARDS' and len(keys)==24 and set(keys)==expected and all(r['status']=='PASS' for r in coverage['configurations'])
roundoff=read('roundoff.json');hist=roundoff['histories']
checks['all_ten_reference_arithmetic_histories']=len(hist)==10 and {(r['role'],r['reference']) for r in hist}=={(role,k) for role in ('source','control') for k in range(5)} and all(r['status']=='CONDITIONAL_PIECEWISE_ARITHMETIC_ENCLOSURE' for r in hist)
interpolation=read('interpolation-roundoff-completion.json')
checks['observation_arithmetic']=interpolation['status']=='CONDITIONAL_PASS' and interpolation['observations']>=19210
events=read('events-assessment.json')
checks['actual_piecewise_crossings']=[r['knot'] for r in events['source']['HP_crossings']]==[94,95] and not events['control']['HP_crossings'] and len(events['source']['native_events'])==24 and not events['control']['native_events']
matched=read('matched.json');pairs=matched['configurations']
checks['all_twelve_matched_budgets']=len(pairs)==12 and {r['configuration'] for r in pairs}=={name for _,name in expected} and all(r['status']=='CONDITIONAL_MATCHED_BUDGET_CANDIDATE' and D(r['R20']['ratio'])<=D(r['R20']['limit']) and all(D(p['budget_ratio'])<=1 for rows in r['endpoint'].values() for p in rows) for r in pairs)
chemical=read('chemical-response.json')
responses={k:dict(temperature_effect_K=r['temperature_effect_K'],resolved=r['signal_gate_10_percent']) for k,r in chemical['responses'].items()}
save(out/'qualification.json',dict(status='CONDITIONAL_NUMERICAL_VERIFICATION_CANDIDATE' if all(checks.values()) else 'PRELIMINARY_FAILED_OR_INCOMPLETE_GATES',checks=checks,failed=[k for k,v in checks.items() if not v],inputs=inputs,chemical_comparison=responses,scope='Fixed accepted inputs, mathematical P2 and exact-zero control, 2048 Julian years. New verification remains an UNVERIFIED SCIENTIFIC CANDIDATE under GOVERNANCE section 5. This aggregation does not confer owner acceptance, absolute stellar accuracy, microscopic-rate validity, or qualification beyond this horizon. assessment.json is the earlier physical-curve generation receipt; this file records the subsequently completed numerical gates.'))
