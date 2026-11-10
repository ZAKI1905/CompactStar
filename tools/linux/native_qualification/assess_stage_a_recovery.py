#!/usr/bin/env python3
"""Read-only actual-output inventory and first-failure assessment; no solver.

Never equate collection completion or raw reference agreement with accuracy.
Missing main-prefix, interpolation and roundoff bounds remain unqualified.
"""
import argparse
import csv
import hashlib
import json
from decimal import Decimal as D, localcontext
from pathlib import Path


def rows(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def lift(x):
    return D.from_float(float(x))


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def assess(run):
    source=run/'run/trajectory/source'
    log=(run/'stage.log').read_text() if (run/'stage.log').exists() else ''
    failures=[]; count=0; passed=0; max_state=None; max_energy=None
    for path in sorted(source.glob('tier-*/checkpoints-*/point-*/checkpoint.json')):
        meta=json.loads(path.read_text());count+=1
        state=rows(path.parent/'state.tsv');diagnostics=rows(path.parent/'diagnostics.tsv')
        record=dict(path=str(path.parent),time_s=meta['time_s'],index=meta['checkpoint_index'],
                    main_step=meta['main_step_ordinal'],O2_atol=meta['O2_atol'])
        bad_state=[r for r in state if r['pass']!='1']
        bad_diagnostics=[r for r in diagnostics if r['pass']!='1']
        if bad_state or bad_diagnostics:
            failures.append(dict(**record,state_failures=bad_state,diagnostic_failures=bad_diagnostics,
                checkpoint=meta,all_state_channels=state,all_diagnostic_channels=diagnostics,
                durable_failure=(path.parent/'failure.txt').read_text()))
        else:passed+=1
        for r in state:
            value=lift(r['utilization'])
            if max_state is None or value>max_state['ratio']:
                max_state=dict(ratio=value,channel=r['channel'],**record)
        for r in diagnostics:
            if r['channel']=='Echem':
                value=lift(r['utilization'])
                if max_energy is None or value>max_energy['ratio']:
                    max_energy=dict(ratio=value,**record,recorded=r)
    references=[]
    for path in sorted(source.glob('reference-*-observations.tsv')):
        rr=rows(path)
        references.append(dict(path=str(path),rows=len(rr),complete=len(rr)==1921 and float(rr[-1]['t'])==31557600,
                               initial=rr[0] if rr else None,final=rr[-1] if rr else None))
    spread={}
    if len(references)==5 and all(r['complete'] for r in references):
        rr=[rows(Path(r['path'])) for r in references]
        assert all([r['t'] for r in a]==[r['t'] for r in rr[0]] for a in rr)
        for field in rr[0][0]:
            if field=='t':continue
            spread[field]=max(max(lift(a[i][field]) for a in rr)-min(lift(a[i][field]) for a in rr) for i in range(1921))
    mains=[]
    for path in sorted(source.glob('tier-*/main.accepted_states.tsv')):
        rr=rows(path);repeat=path.parent/'endpoint-schedule-repeat/main.accepted_states.tsv'
        events=rows(path.parent/'events.tsv')
        r20=rows(path.parent/'R20-grids.tsv') if (path.parent/'R20-grids.tsv').exists() else []
        mains.append(dict(path=str(path),accepted_steps=len(rr),final=rr[-1] if rr else None,
            complete=bool(rr) and float(rr[-1]['t_after_s'])==31557600,
            schedule_repeat_identical=repeat.exists() and path.read_bytes()==repeat.read_bytes(),
            actual_events=sum(int(r['knot'])>=0 for r in events),R20=r20,
            restart_completed=(path.parent/'restart-comparison.tsv').exists()))
    r18=[]
    for path in sorted(source.glob('tier-*/diagnostics-*.tsv')):
        for r in rows(path):
            gross=max(D(1),sum(abs(lift(r[k])) for k in ('P_dir_actual_erg_s','LH_erg_s','DeltaLnu_erg_s',
                'Lnu_eq_erg_s','Lgamma_erg_s','Lother_erg_s','L_out_fluid_inf_erg_s')))
            r18.append(dict(path=str(path),time_s=r['t_s'],ratio=abs(lift(r['R18_residual_erg_s']))/gross,
                Tinf_K=r['Tinf_K'],valid_through_sample=r['valid_through_sample'],
                frozen_utilization=r['max_frozen_utilization'],depletion=r['DeltaB_over_B0']))
    status='FAIL_NATIVE_CHECKPOINT' if failures else ('SOURCE_COLLECTION_COMPLETE_UNQUALIFIED'
        if 'ROLE COLLECTION COMPLETE source' in log else 'INCOMPLETE_SOURCE_COLLECTION')
    result=dict(status=status,candidate=True,source_qualified=False,control_executed=(run/'run/trajectory/control').exists(),
        actual_checkpoint_packets=count,passed_native_packets=passed,failures=failures,
        maximum_native_state_utilization=max_state,maximum_native_Echem_utilization=max_energy,
        source_references=references,raw_reference_spread_not_a_bound=spread,mains=mains,
        published_R18_samples=len(r18),maximum_R18_ratio=max(r18,key=lambda x:x['ratio']) if r18 else None,
        remaining_accuracy_obligations={'chemical_high_precision_actual_checkpoints':'NOT_RUN_AFTER_NATIVE_FAILURE' if failures else 'PENDING',
          'main_prefix_and_accumulated_roundoff':'UNQUALIFIED', 'thermal_RHS_roundoff':'UNQUALIFIED',
          'global_reference_and_interpolation':'UNQUALIFIED; raw spreads do not bound accumulated error',
          'R20_total_uncertainty':'UNQUALIFIED; report original grid residuals separately',
          'restart_consistency':'NOT_RUN' if not any(x['restart_completed'] for x in mains) else 'PENDING',
          'paired_signal':'NOT_AVAILABLE until source qualifies and control passes'},
        assumptions=['Uniform thermal enclosure remains conditional; sampled agreement is insufficient.',
                     'Chemical reference allowance assumes continuation of refinement contraction.',
                     'The 13-point D1 study is not universal checkpoint qualification.'],
        stage_B_executed=False,stellar_regeneration=False,automatic_retry=False,
        source_log_sha256=sha(run/'stage.log') if (run/'stage.log').exists() else None)
    if 'RESOURCE_USAGE_JSON ' in log:
        result['resource_usage']=json.loads(log.split('RESOURCE_USAGE_JSON ',1)[1].splitlines()[0])
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('run',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    with localcontext() as ctx:
        ctx.prec=110
        result=assess(a.run)
    (a.output/'assessment.json').write_text(json.dumps(result,indent=2,default=str)+'\n')
    print(result['status'],'actual packets',result['actual_checkpoint_packets'],'failures',len(result['failures']))
