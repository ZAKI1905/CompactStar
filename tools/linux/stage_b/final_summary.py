#!/usr/bin/env python3
"""Summarize saved verification results without changing scientific gates."""
import json,pathlib,sys
from decimal import Decimal as D
from piecewise import ROOT
from run import save,digest

out=pathlib.Path(sys.argv[1]);bindings={}
def read(name):
 p=out/(name+'.json');bindings[p.name]=digest(p)
 return json.loads(p.read_text())
qualification=read('qualification');accuracy=read('accuracy')['configurations']
ledger=read('ledger-assessment')['configurations'];local=read('local-assessment')['configurations']
matched=read('matched')['configurations'];events=read('events-assessment');coverage=read('coverage')
physical=read('assessment');chemical=read('chemical-response')
summary={}
for role in ('source','control'):
 a=[r for r in accuracy if r['role']==role];l=[r for r in ledger if r['role']==role]
 loc=[r for r in local if r['role']==role]
 details=[read('accuracy/'+role+'-'+r['configuration']) for r in a]
 accepted=[v for r in details for v in r['accepted']+(r['restart']['records'] if r['restart'] else [])]
 ap=[v for r in details for v in r['accepted_power']]
 summary[role]=dict(
  configurations=len(a),accuracy_maxima={k:max(D(r['maxima'][k]) for r in a) for k in a[0]['maxima']},
  accepted_restart_state_ratio=max(D(v) for r in accepted for v in r['prescribed_state_ratios']),
  accepted_restart_reference_allocation=max(D(v) for r in accepted for v in r['reference_allocations']),
  accepted_power_ratio=max(D(r['budget_ratio']) for r in ap),
  accepted_power_reference_allocation=max(D(r['reference_allocation_ratio']) for r in ap),
  ledger_maxima={k:max(D(r[k]) for r in l) for k in ('total_over_N','reconstruction_event_allocation_over_N','max_R18_over_G','maximum_depletion','maximum_frozen_utilization','finite_T_omitted_scale_K')},
  local_state_radius_ratio=max(D(v) for r in loc for v in r['worst_state']['state_ratios']),
  local_chemical_energy_radius_ratio=max(D(r['worst_Echem']['Echem_ratio']) for r in loc),
  local_points=sum(r['points'] for r in loc),
  endpoint=physical['roles'][role]['endpoint'],
  descriptive_integrated_powers_erg=physical['roles'][role]['integrated_powers_trapezoid_erg'])
table=[]
for r in matched:
 name=r['configuration'];la=[v for v in ledger if v['configuration']==name]
 table.append(dict(configuration=name,temperature=r['temperature'],R20=r['R20'],
  unresolved_endpoint_power_channels=[v['channel'] for v in r['endpoint']['powers'] if v['signal']!='RESOLVED'],
  source_R20_over_N=next(v['total_over_N'] for v in la if v['role']=='source'),
  control_R20_over_N=next(v['total_over_N'] for v in la if v['role']=='control')))
save(out/'final-summary.json',dict(status=qualification['status'],failed_checks=qualification['failed'],
 aggregate_checks=len(qualification['checks']),roles=summary,paired_HP_endpoint=physical['paired_endpoint'],
 matched_configurations=table,chemical_temperature_comparison=qualification['chemical_comparison'],
 source_crossings=events['source']['HP_crossings'],
 maximum_native_event_time_error_bound_s=max(D(r['reference_time_error_bound_s']) for r in events['source']['native_events']),
 coverage=coverage['status'],inputs=bindings,
 scope='Display summary of completed candidate verification. No new trajectory, modified threshold, owner acceptance, or absolute stellar accuracy. Individual integrated powers are descriptive trapezoids, not separately qualified integrals.'))
