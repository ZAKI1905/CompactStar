#!/usr/bin/env python3
"""Verify full matrix, endpoint schedule invariance, nested grids, and guards."""
import json,pathlib,sys
from decimal import Decimal as D
from piecewise import END,rows,lift
from run import digest,save
out=pathlib.Path(sys.argv[1]);records=[]
for role in ('source','control'):
 for tier in range(3):
  for cap in range(4):
   p=out/'native'/role/f'tier-{tier}-cap-{cap}';hist=rows(p/'main.accepted_states.tsv');repeat=rows(p/'endpoint-schedule-repeat/main.accepted_states.tsv')
   assert hist==repeat,'observation schedule changed accepted history'
   assert lift(hist[-1]['t_after_s'])==END
   assert int(hist[-1]['cumulative_rejected'])<=D('.1')*(len(hist)+int(hist[-1]['cumulative_rejected']))
   last=rows(p/'diagnostics-1920.tsv');assert len(last)==1921
   all_keys=list(last[0]);common=0
   for grid in (240,480,960,1920):
    data=rows(p/f'diagnostics-{grid}.tsv');unc=rows(p/f'uncertainties-{grid}.tsv');assert len(data)==len(unc)==grid+1
    for i,r in enumerate(data):
     reference=last[i*(1920//grid)]
     assert all(r[k]==reference[k] for k in all_keys),'nested-grid diagnostic mismatch'
     assert lift(r['t_s'])==END*i/grid
     if role=='control':assert lift(r['eta_e_MeV'])==lift(r['eta_mu_MeV'])==0
     assert D('1e6')<=lift(r['Tinf_K'])<=D('1e9')
     assert abs(lift(r['DeltaB_over_B0']))<=D('1e-6') and lift(r['max_frozen_utilization'])<=1
     common+=1
   ledger=rows(p/'R20-grids.tsv');assert len(ledger)==4
   restart=None
   if cap==3:
    tail=rows(p/'restart-tail/main.accepted_states.tsv');start=lift(tail[0]['t_before_s'])
    expected=next(r for r in hist if lift(r['t_after_s'])>=END/2)
    assert start==lift(expected['t_after_s']) and lift(tail[-1]['t_after_s'])==END
    restart=dict(start_s=str(start),accepted_steps=len(tail),tail_sha256=digest(p/'restart-tail/main.accepted_states.tsv'))
   events=[r for r in rows(p/'events.tsv') if int(r['knot'])>=0]
   assert [int(r['knot']) for r in events]==([94,95] if role=='source' else [])
   records.append(dict(role=role,tier=tier,cap=cap,status='PASS',accepted_steps=len(hist),events=len(events),nested_samples=common,history_sha256=digest(p/'main.accepted_states.tsv'),restart=restart))
save(out/'coverage.json',dict(status='PASS_MATRIX_AND_NATIVE_GUARDS',configurations=records,scope='All 24 configurations, 24 endpoint schedule repeats, 96 diagnostic grids, six restarts and actual event counts. Native completion additionally enforces event endpoint bypass and both nextafter neighbors. Numerical qualification still requires all independent error budgets.'))
