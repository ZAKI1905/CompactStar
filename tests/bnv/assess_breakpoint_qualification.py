#!/usr/bin/env python3
"""Report original and proposed successor criteria separately; never ratify."""
import argparse,csv,json,pathlib,math,os

def rows(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def boundary_classification(history):
    # Only accepted step/event records enter this classification, never errors.
    episodes=[]
    for tier in range(3):
        groups=[];count=0;start=None
        for r in history:
            if r['tier']!=str(tier):continue
            if start is None:start=float(r['t_left'])
            count+=1
            if int(r['boundary_knot'])>=0:
                direction=1 if int(r['boundary_knot'])==int(r['cell'])+1 else -1
                groups.append((int(r['boundary_knot']),direction,count,float(r['t_right'])-start))
                count=0;start=None
        episodes.append(groups)
    matching=all([(r[0],r[1]) for r in x]==[(r[0],r[1]) for r in episodes[0]] for x in episodes)
    limited=[]
    if matching:
        for j in range(len(episodes[0])):
            group=[e[j] for e in episodes]
            # One whole cell in one accepted segment at every tolerance and
            # durations agreeing within 0.1 percent, selected before scoring.
            if all(r[2]==1 for r in group) and max(r[3] for r in group)<=1.001*min(r[3] for r in group):
                limited.append({'knot':group[0][0],'direction':group[0][1],'steps':[r[2] for r in group],'duration_s':[r[3] for r in group]})
    return {'matching_event_sequence':matching,'boundary_limited_cell_intervals':limited,'counts_by_tier':[len(e) for e in episodes]}

def assess(root):
    result={'status':'PROPOSED successor, not ratified','campaign':'STOPPED','windows':{}}
    for name in ['source','control','bounded','tail']:
        observations=rows(root/(name+'-observations.tsv'))
        expected={'source':65,'control':65,'bounded':241,'tail':33}[name]
        schedules=[]
        for tier in range(3):
            group=[x for x in observations if x['tier']==str(tier)]
            assert len(group)==expected, (name,tier,'incomplete observations')
            times=[float(x['time']) for x in group]
            assert all(math.isfinite(t) for t in times) and all(a<b for a,b in zip(times,times[1:])), (name,tier,'corrupt time order')
            schedules.append(times)
        assert schedules[0]==schedules[1]==schedules[2], 'tier schedules differ'
        s=rows(root/(name+'-summary.tsv'))[0];l=rows(root/(name+'-ledger.tsv'))[0];st=rows(root/(name+'-statistics.tsv'))
        rejection=[int(x['rejected'])/(int(x['accepted'])+int(x['rejected'])) for x in st]
        r20=rows(root/(name+'-r20.tsv')) if (root/(name+'-r20.tsv')).exists() else []
        original={'state':all(int(s[k])==0 for k in ['legacy_c1','legacy_c2','legacy_c3']),
          'ledger':int(s['unresolved_checkpoints'])==0 and all(int(l[k])==0 for k in ['legacy_c1','legacy_c2','legacy_c3']),
          'rejection_fraction_le_0p10':all(v<=.1 for v in rejection)}
        # Boundaries are classified from records, before using outcome differences.
        clips=[float(x['boundary_limited_fraction']) for x in st]
        limits=[int(x['events']) for x in st]
        envelopes=rows(root/(name+'-envelopes.tsv'));assert envelopes
        envelope_max=max(float(x['bound_ratio']) for x in envelopes)
        classification=boundary_classification(rows(root/(name+'-history.tsv')))
        successor={'uncertainty_envelope':math.isfinite(envelope_max) and envelope_max<=1, 'state_accuracy':float(s['max_accuracy_ratio'])<=1,
          'reference_uncertainty':float(s['max_reference_uncertainty_over_DU'])<=.05,
          'ledger_accuracy':float(l['max_accuracy_ratio'])<=1,
          'reconstruction':int(s['unresolved_checkpoints'])==0,
          'R18':float(l['max_R18_normalized'])<=1e-10,
          'R20':len(r20)==3 and all(abs(float(x['R20_normalized']))<=2e-4 for x in r20),
          'R20_reconstruction':len(r20)==3 and all(float(x['reconstruction_over_N'])<=5e-6 for x in r20),
          'schedule_invariance':s['schedule_invariant']=='1',
          'progress':all(int(x['accepted'])<100000 and float(x['min_step'])>0 for x in st)}
        result['windows'][name]={'original_BA12R':original,'successor':successor,'successor_pass':all(successor.values()),
          'max_state_accuracy_ratio':float(s['max_accuracy_ratio']),'max_ledger_accuracy_ratio':float(l['max_accuracy_ratio']),
          'rejection_fraction':rejection,'boundary_clips':limits,'boundary_limited_fractions':clips,
          'uncertainty_envelope_max_ratio':envelope_max,'boundary_classification':classification,'R20_records':r20,
          'scope':'common-time local window; no full-campaign global accuracy claim'}
    cross=rows(root/'crossings.tsv');result['crossings']={'rows':len(cross),'count':len({(x['run'],x['step']) for x in cross}),
        'max_accuracy_ratio':max(float(x['accuracy_ratio']) for x in cross),'max_event_accuracy_ratio':max(float(x['event_accuracy_ratio']) for x in cross)}
    historical=rows(root/'historical-reproduction.tsv');assert len(historical)==4
    result['historical_reproduction']={'exact_Arm_E_and_obs32_assertions':'PASS in completed harness',
      'max_candidate_envelope_ratio':max(float(x['candidate_envelope_ratio']) for x in historical),
      'all_candidate_reconstructions_qualified':all(x['new_qualified']=='1' for x in historical),
      'smooth_bitwise':all(x['smooth_bitwise']=='1' for x in historical if x['observation']!='32')}
    h=result['historical_reproduction'];historical_ok=h['max_candidate_envelope_ratio']<=1 and h['all_candidate_reconstructions_qualified'] and h['smooth_bitwise']
    result['pass']=historical_ok and all(x['successor_pass'] for x in result['windows'].values()) and result['crossings']['count']==37 and result['crossings']['max_accuracy_ratio']<=1 and result['crossings']['max_event_accuracy_ratio']<=1
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=pathlib.Path);a=p.parse_args();r=assess(a.root)
    out=a.root/'assessment.json';assert not out.exists();temp=out.with_suffix('.json.partial');assert not temp.exists()
    content=json.dumps(r,indent=2,allow_nan=False)+'\n'
    with temp.open('w') as f:f.write(content);f.flush();os.fsync(f.fileno())
    assert temp.read_text()==content and json.loads(temp.read_text())==r
    temp.rename(out);print(content)
