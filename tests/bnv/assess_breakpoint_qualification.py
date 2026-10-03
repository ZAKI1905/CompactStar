#!/usr/bin/env python3
"""Report original and proposed successor criteria separately; never ratify."""
import argparse,csv,json,pathlib,math

def rows(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def assess(root):
    result={'status':'PROPOSED successor, not ratified','campaign':'STOPPED','windows':{}}
    for name in ['source','control','bounded','tail']:
        s=rows(root/(name+'-summary.tsv'))[0];l=rows(root/(name+'-ledger.tsv'))[0];st=rows(root/(name+'-statistics.tsv'))
        rejection=[int(x['rejected'])/(int(x['accepted'])+int(x['rejected'])) for x in st]
        r20=rows(root/(name+'-r20.tsv')) if (root/(name+'-r20.tsv')).exists() else []
        original={'state':all(int(s[k])==0 for k in ['legacy_c1','legacy_c2','legacy_c3']),
          'ledger':all(int(l[k])==0 for k in ['legacy_c1','legacy_c2','legacy_c3']),
          'rejection_fraction_le_0p10':all(v<=.1 for v in rejection)}
        # Boundaries are classified from records, before using outcome differences.
        clips=[float(x['boundary_limited_fraction']) for x in st]
        limits=[int(x['events']) for x in st]
        successor={'state_accuracy':float(s['max_accuracy_ratio'])<=1,
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
          'scope':'common-time local window; no full-campaign global accuracy claim'}
    cross=rows(root/'crossings.tsv');result['crossings']={'rows':len(cross),'count':len({(x['run'],x['step']) for x in cross}),
        'max_accuracy_ratio':max(float(x['accuracy_ratio']) for x in cross)}
    result['pass']=all(x['successor_pass'] for x in result['windows'].values()) and result['crossings']['count']==37 and result['crossings']['max_accuracy_ratio']<=1
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=pathlib.Path);a=p.parse_args();r=assess(a.root)
    out=a.root/'assessment.json';assert not out.exists();out.write_text(json.dumps(r,indent=2)+'\n');print(out.read_text())
