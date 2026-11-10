#!/usr/bin/env python3
"""Descriptive retained-output context after the Stage A stop; no solver calls."""
import argparse
import csv
import hashlib
import json
from decimal import Decimal as D, localcontext
from pathlib import Path

ROOT=Path(__file__).resolve().parents[3]


def rows(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))


def generate(run,output):
    source=run/'run/trajectory/source'
    a=json.loads((run/'assessment/assessment.json').read_text())
    old=ROOT/'docs/validation/results-first-bnv-pilot/native-runtime-retry/evidence/main.accepted_states.tsv'
    new=source/'tier-0-cap-0/main.accepted_states.tsv'
    assert new.read_bytes()==old.read_bytes()
    originalmanifest=ROOT/'docs/validation/results-first-bnv-pilot/native-runtime-retry/evidence/external-SHA256SUMS'
    oldhash={n:h for h,n in (l.split(maxsplit=1) for l in originalmanifest.read_text().splitlines())}
    ref_matches=[]
    oldrun=ROOT.parents[1]/'native-p2-qualification-1yr-r2/run/trajectory/source'
    for i in range(5):
        for suffix in ('history','observations'):
            prev=oldrun/f'reference-{i}-{suffix}.tsv';current=source/prev.name
            assert hashlib.sha256(prev.read_bytes()).hexdigest()==oldhash[str(prev)]
            ref_matches.append(dict(file=prev.name,identical=prev.read_bytes()==current.read_bytes()))
    assert all(r['identical'] for r in ref_matches)
    ref=rows(source/'reference-4-observations.tsv');main=rows(new)
    def lift(x):return D.from_float(float(x))
    H=lift(31557600)
    result={'reference4_end':ref[-1],
        'source_temperature_increase_K':lift(ref[-1]['Tinf'])-lift(ref[0]['Tinf']),
        'main_endpoint_difference_from_reference4':{k:abs(lift(main[-1][col])-lift(ref[-1][k]))
            for k,col in [('x','x_state'),('e','eta_e_MeV'),('mu','eta_mu_MeV')]},
        'restart_last_before_half_years':lift(a['restart_gate']['last_endpoint_before_half_s'])/H,
        'smallest_cap_years':lift(a['restart_gate']['minimum_planned_cap_s'])/H,
        'chemical_heating_over_direct_deposition':lift(ref[-1]['LH'])/lift(ref[-1]['P_dir_actual']),
        'raw_reference_spreads_are_not_error_bounds':True}
    result['reference4_1920_grid_direct_deposition_erg']=sum(
        (lift(b['t'])-lift(a['t']))*(lift(a['P_dir_actual'])+lift(b['P_dir_actual']))/2 for a,b in zip(ref,ref[1:]))
    result['reference4_1920_grid_outgoing_erg']=sum(
        (lift(b['t'])-lift(a['t']))*sum(lift(a[k])+lift(b[k]) for k in ('Lnu_full','Lgamma','Lother'))/2
        for a,b in zip(ref,ref[1:]))
    result['reference4_1920_grid_LH_erg']=sum(
        (lift(b['t'])-lift(a['t']))*(lift(a['LH'])+lift(b['LH']))/2 for a,b in zip(ref,ref[1:]))
    table=[]
    for p in sorted((source/'tier-0-cap-0/checkpoints-240').glob('point-*/checkpoint.json'),
                    key=lambda p:int(p.parent.name.split('-')[-1])):
        m=json.loads(p.read_text());s=rows(p.parent/'state.tsv');ds=rows(p.parent/'diagnostics.tsv')
        e=next((v for v in ds if v['channel']=='Echem'),None)
        table.append(dict(index=m['checkpoint_index'],time_s=m['time_s'],O2_x=m['O2_state'][0],
            O2_e=m['O2_state'][1],O2_mu=m['O2_state'][2],radius_e_MeV=s[1]['propagated_state_radius'],
            radius_mu_MeV=s[2]['propagated_state_radius'],Echem_MeV=e['O2'] if e else None,
            U_Echem_MeV=e['U'] if e else None,Echem_threshold_MeV=e['threshold'] if e else None,
            Echem_utilization=e['utilization'] if e else None,
            native_all_channels_pass=all(v['pass']=='1' for v in s+ds),high_precision_actual_output_status='NOT_REACHED'))
    assert len(table)==129 and all(v['native_all_channels_pass'] for v in table)
    with (output/'checkpoints.tsv').open('w') as f:
        wr=csv.DictWriter(f,fieldnames=list(table[0]),delimiter='\t',lineterminator='\n')
        wr.writeheader();wr.writerows(table)
    result['label']='DIAGNOSTIC ONLY; source not qualified, control not run'
    result['reference_and_main_regression']={'main_identical':True,'references':ref_matches}
    (output/'scientific-context.json').write_text(json.dumps(result,indent=2,default=str)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('run',type=Path);p.add_argument('output',type=Path)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    with localcontext() as ctx:
        ctx.prec=110
        generate(args.run,args.output)
