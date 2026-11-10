#!/usr/bin/env python3
"""Prove the prescribed restart is unavailable from an actual accepted history.

No integration. An inactive-cap induction is used, explicitly distinguished
from execution of the not-yet-run quarter-cap member of the fixed matrix.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[3]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def assess(history):
    with history.open() as f:
        rows=list(csv.DictReader(f,delimiter='\t'))
    repeat=history.parent/'endpoint-schedule-repeat/main.accepted_states.tsv'
    assert history.read_bytes()==repeat.read_bytes(),'schedule repeat differs'
    with (history.parent/'events.tsv').open() as f:
        events=list(csv.DictReader(f,delimiter='\t'))
    assert all(int(r['knot'])==-1 for r in events),'eventful path needs separate proof'
    end=31557600.;half=end/2
    assert float(rows[0]['t_before_s'])==0 and float(rows[-1]['t_after_s'])==end
    assert all(int(r['cumulative_rejected'])==0 for r in rows),'rejected trials not captured'
    caps=(38723807757.,19361903878.5,9680951939.25)
    checks=[]
    incoming=1.
    previous=0.
    for index,r in enumerate(rows,1):
        t=float(r['t_before_s'])
        assert int(r['sequence'])==index and t==previous,'history not contiguous'
        checks.append(dict(sequence=int(r['sequence']),left_s=t,right_s=float(r['t_after_s']),
            incoming_proposal_s=incoming,
            caps=[dict(cap_s=cap,proposal_unchanged=min(incoming,cap)==incoming,
                       gsl_target_unchanged=min(end,max(t+cap,math.nextafter(t,float('inf'))))==end)
                  for cap in caps]))
        incoming=float(r['suggested_next_h_s'])
        previous=float(r['t_after_s'])
    assert all(c['proposal_unchanged'] and c['gsl_target_unchanged'] for r in checks for c in r['caps'])
    stepper=ROOT/'CompactStar/Physics/BNV/BreakpointStepper.hpp'
    stage=ROOT/'tests/native/stage.cpp'
    text=stepper.read_text();collector=stage.read_text()
    assert 'std::max(std::min(h,max_step_),std::nextafter(t,INFINITY)-t)' in text
    assert 'Apply(t,y,h,std::min(target,std::max(t+max_step_,std::nextafter(t,INFINITY))));' in text
    assert 'return s.t_right_s>=end/2;' in collector
    assert 'it!=main.accepted_steps.end()&&it->t_right_s<end,"no interior restart endpoint"' in collector
    selected=next(r for r in rows if float(r['t_after_s'])>=half)
    valid=float(selected['t_after_s'])<end
    return dict(status='PASS_RESTART_AVAILABLE' if valid else 'FAIL_REQUIRED_RESTART_COVERAGE',
        role='source',tier='BASELINE',observed_policy='uncapped',required_restart_policy='quarter-cap',
        quarter_cap_executed=False,method='Exact inactive-cap induction on actual uncapped steps; not a new trajectory.',
        end_time_s=end,half_time_s=half,selected_endpoint=selected,
        last_endpoint_before_half_s=max(float(r['t_after_s']) for r in rows if float(r['t_after_s'])<half),
        restart_tail_duration_s=end-float(selected['t_after_s']),
        maximum_proposed_step_s=max([1.]+[float(r['suggested_next_h_s']) for r in rows]),
        minimum_planned_cap_s=min(caps),cap_checks=checks,
        proof=['The only numerical uses of max_step_ in Advance are min(h,max_step_) and the GSL target clamp.',
               'For every actual step, every planned finite cap leaves both operands unchanged; there were no rejected trials or events.',
               'Starting from the same exact state, stepper and deterministic RHS, induction gives the same GSL calls, accepted states and next-step proposals for quarter-cap.',
               'Therefore the required first accepted endpoint >=half-time is the final endpoint for that policy as well; the existing strict-interior restart guard must fail.',
               'No cap, restart time, numerical threshold or accepted-step history is altered to create coverage.'],
        source_locations={'cap_semantics':str(stepper),'restart_rule':str(stage)},
        hashes={str(p):sha(p) for p in (history,repeat,history.parent/'events.tsv',stepper,stage)},
        scientific_interpretation='Missing positive-length restart coverage, not a demonstrated temperature or chemical-state accuracy failure.',
        required_disposition='Stop the current retry and do not launch control; preserve partial output.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('history',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();result=assess(a.history)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'], 'restart tail seconds:',result['restart_tail_duration_s'])
