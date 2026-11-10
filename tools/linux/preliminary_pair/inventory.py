#!/usr/bin/env python3
"""Complete control coverage, native acceptance and schedule-invariance audit."""
import concurrent.futures,json,pathlib,sys,time
from high_precision import rows,lift,save,digest
out=pathlib.Path(sys.argv[1]);base=out/'control';start=time.monotonic()
ignored={'checkpoint_index','main_rtol','main_atol','main_maximum_step_s'}
canonical={};maximum_state=(0,None);maximum_diagnostic=(0,None)
for i in range(1921):
 p=base/f'tier-0-cap-0/checkpoints-1920/point-{i}'
 meta=json.loads((p/'checkpoint.json').read_text())
 assert meta['state'][1:]==[0,0] and meta['O2_state'][1:]==[0,0]
 for kind in ('state','diagnostics'):
  for r in rows(p/(kind+'.tsv')):
   assert r['pass']=='1',(p,kind,r)
   item=(float(r['utilization']),dict(point=i,channel=r['channel']))
   if kind=='state' and item[0]>maximum_state[0]:maximum_state=item
   if kind=='diagnostics' and item[0]>maximum_diagnostic[0]:maximum_diagnostic=item
 canonical[i]=({k:digest(p/k) for k in ('state.tsv','diagnostics.tsv')},{k:v for k,v in meta.items() if k not in ignored})

def audit(p):
 count=0;grids=0;dmax=0;umax=0
 accepted=rows(p/'main.accepted_states.tsv');assert len(accepted)==12
 assert lift(accepted[-1]['t_after_s'])==31557600 and all(int(r['cumulative_rejected'])==0 for r in accepted)
 assert digest(p/'main.accepted_states.tsv')==digest(p/'endpoint-schedule-repeat/main.accepted_states.tsv')
 assert all(int(r['knot'])==-1 for r in rows(p/'events.tsv'))
 for n in (240,480,960,1920):
  ds=rows(p/f'diagnostics-{n}.tsv');assert len(ds)==n+1
  for r in ds:
   assert r['valid_through_sample']=='1'
   assert 1e6<=float(r['Tinf_K'])<=1e9
   assert all(float(r[k])==0 for k in ('eta_e_MeV','eta_mu_MeV','LH_erg_s','DeltaLnu_erg_s','Echem_MeV','P_dir_actual_erg_s','Bdot_count_s'))
   dmax=max(dmax,abs(float(r['DeltaB_over_B0'])));umax=max(umax,float(r['max_frozen_utilization']))
  for i in range(n+1):
   q=p/f'checkpoints-{n}/point-{i}';hs,meta=canonical[i*(1920//n)]
   assert all(digest(q/k)==h for k,h in hs.items())
   actual=json.loads((q/'checkpoint.json').read_text());assert {k:v for k,v in actual.items() if k not in ignored}==meta
   assert 'status=QUALIFIED' in (q/'failure.txt').read_text()
   count+=1
  grids+=1
 assert dmax<=1e-6 and umax<=1
 return dict(tier=p.name,packets=count,grids=grids,accepted_steps=len(accepted),rejected=0,events=0,schedule_repeat_identical=True,main_sha256=digest(p/'main.accepted_states.tsv'),maximum_depletion=dmax,maximum_frozen_utilization=umax)

with concurrent.futures.ThreadPoolExecutor(max_workers=24) as pool:result=list(pool.map(audit,sorted(base.glob('tier-*'))))
assert len(result)==12 and sum(r['packets'] for r in result)==43248
assert len({r['main_sha256'] for r in result})==1
save(out/'control-inventory.json',dict(status='PASS_NATIVE_CHECKS_AND_COVERAGE; NOT_FORMAL_QUALIFICATION',packets=43248,grids=48,main_histories=12,endpoint_repeats=12,restart_tails=len(list(base.glob('tier-*/restart-comparison.tsv'))),unique_grid_points=1921,records=result,maximum_state_utilization=maximum_state,maximum_diagnostic_utilization=maximum_diagnostic,elapsed_s=time.monotonic()-start,workers=24))
print('PASS 43248 control packets, 12 histories, 48 grids, exact-zero chemistry')
