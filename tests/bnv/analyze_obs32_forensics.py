#!/usr/bin/env python3
"""Read-only reduction of local diagnostic traces; no ODE execution or acceptance promotion."""
import argparse,csv,ctypes,gzip,io,json,math
from pathlib import Path
fma=ctypes.CDLL(None).fma;fma.argtypes=[ctypes.c_double]*3;fma.restype=ctypes.c_double

def read(path):
    return path.read_text() if path.exists() else gzip.decompress(Path(str(path)+'.gz').read_bytes()).decode()
def rows(path):return list(csv.DictReader(io.StringIO(read(path)),delimiter='\t'))
def vec(r):return [float(r[k]) for k in ('x','eta_e','eta_mu')]
def delta(a,b):return [abs(x-y) for x,y in zip(a,b)]
def budgets(a,b,c,level=1):
    aa=[(1e-17,1e-23,1e-23),(1e-18,1e-24,1e-24),(1e-19,1e-25,1e-25)][level-1]
    ab=[(1e-18,1e-24,1e-24),(1e-19,1e-25,1e-25),(1e-20,1e-26,1e-26)][level-1]
    ra=[1e-12,1e-13,1e-14][level-1];rb=[1e-13,1e-14,1e-15][level-1];out=[]
    for i,k in enumerate(('x','eta_e','eta_mu')):
        mo=max(abs(a[i]),abs(b[i]));mi=max(mo,abs(float(c[k+'_L'])),abs(float(c[k+'_R'])));d=abs(a[i]-b[i]);D=fma(ra,mo,aa[i]);F=max(fma(1e-11,mi,(1e-16,1e-22,1e-22)[i]),64*math.ulp(mi));FO=max(fma(rb,mo,ab[i]),64*math.ulp(mo));U=2*max(d,FO)
        out.append({'component':k,'d':d,'D':D,'F_i':F,'F_O':FO,'U':U,'d_over_D':d/D,'U_over_point2F':U/(.2*F),'pass':d<=D and U<=.2*F})
    return out

def analyze(root,evidence,supplement=None):
    cases={r['index']:r for r in rows(evidence/'cases.tsv')};c=cases['32'];s={r['name']:r for r in rows(root/'solutions.tsv')};event=rows(root/'event.tsv')[0];cache=rows(root/'cstar-cache.tsv');ki=int(event['knot_index']);knot=cache[ki];Tk=float(knot['T_K']);xk=math.log(Tk/1e8);C=float(knot['Cstar'])
    slope_lo=(C-float(cache[ki-1]['Cstar']))/(math.log(float(knot['T_MeV']))-math.log(float(cache[ki-1]['T_MeV'])))/Tk
    slope_hi=(float(cache[ki+1]['Cstar'])-C)/(math.log(float(cache[ki+1]['T_MeV']))-math.log(float(knot['T_MeV'])))/Tk
    result={'scope':'BOUNDED LOCAL DIAGNOSTICS; HISTORICAL CAMPAIGN REMAINS FAIL; NO CANDIDATE','case':c,'event':event,'cache_knot':knot,'temperatures':{'T_L_K':1e8*math.exp(float(c['x_L'])),'T_R_K':1e8*math.exp(float(c['x_R']))},'Cstar':{'slope_below_erg_K2':slope_lo,'slope_above_erg_K2':slope_hi,'slope_jump_erg_K2':slope_hi-slope_lo,'relative_slope_jump':(slope_hi-slope_lo)/slope_lo},'solutions':{},'crossing_traces':{},'smooth':{}}
    ref=vec(s['split-right-O3'])
    for name,row in s.items():
        if name.startswith('root-'):continue
        result['solutions'][name]={**row,'difference_from_split_rk8pd_O3':delta(vec(row),ref) if float(row['t_end'])==float(c['t_obs']) else None}
    result['governed_unsplit_budgets']=budgets(vec(s['unsplit-O1']),vec(s['unsplit-O2']),c)
    # Independent arithmetic must reproduce every archived governed budget exactly.
    witness=json.loads(read(evidence/'saved-witness.json'))
    for b in result['governed_unsplit_budgets']:
        k=b['component']
        for field,archived in [('d','d_'),('D','D_O1_'),('F_i','F_'),('U','U_')]:
            assert b[field]==float(witness[archived+k]), (field,k)
    for l in (1,2):
        for k in ('x','eta_e','eta_mu'):
            assert float(s[f'unsplit-O{l}'][k])==float(witness[f'{k}_O{l}'])
    result['archived_endpoints_and_budgets_exact']=True
    result['split_O1_O2_budgets']=budgets(vec(s['split-right-O1']),vec(s['split-right-O2']),c)
    result['split_O2_O3_difference']=delta(vec(s['split-right-O2']),ref)
    result['diagnostic_unsplit_O3_O4_budgets']=budgets(vec(s['unsplit-O3']),vec(s['unsplit-O4']),c,3)
    for l in (1,2,3,4):
        name=f'unsplit-O{l}';attempts=rows(root/(name+'.steps.tsv'));rhs=rows(root/(name+'.rhs.tsv'))
        crossings=[r for r in attempts if float(r['x_before'])<xk<=float(r['x_after'])]
        accepted=[r for r in crossings if r['accepted']=='1'];rejected=[r for r in crossings if r['accepted']=='0'];assert accepted
        cell_counts={}
        for r in rhs:cell_counts[r['cell_after']]=cell_counts.get(r['cell_after'],0)+1
        # Every logged trial uses a valid containing cache cell, including rejected trials.
        for r in rhs:
            j=int(r['cell_after']);T=float(r['T_MeV'])
            assert float(cache[j]['T_MeV'])<=T<=float(cache[j+1]['T_MeV'])
        result.setdefault('all_trial_cells_contain_temperature',{})[name]=True
        first=accepted[0]; near=[r for r in attempts if abs(float(r['t'])-float(event['t_chosen']))<1e6 or r in crossings]
        result['crossing_traces'][name]={'first_rejected_crossing':rejected[0] if rejected else None,'first_accepted_crossing':first,'actual_rhs_cell_counts':cell_counts,'near_knot_attempts':near}
        assert sum(r['accepted']=='1' for r in attempts)==int(s[name]['accepted'])
        assert sum(r['accepted']=='0' for r in attempts)==int(s[name]['rejected'])
    for key,case in cases.items():
        if key=='32':continue
        assert not any(float(case['x_L'])<math.log(float(k['T_K'])/1e8)<float(case['x_R']) for k in cache[1:-1])
        a,b,d=[vec(s[f'smooth-{key}-O{l}']) for l in (1,2,3)];result['smooth'][key]={'case':case,'O1_O2':budgets(a,b,case),'O2_O3_difference':delta(b,d),'all_three_equal':a==b==d,'local_duration_ratio_to_failure':(float(case['t_obs'])-float(case['t_L']))/(float(c['t_obs'])-float(c['t_L']))}
    cells=rows(root/'cell-ownership.tsv');result['cell_ownership']=cells;at=[r for r in cells if r['point']=='0'];assert len(at)==2 and at[0]['Cstar']==at[1]['Cstar']
    rhs=rows(root/'rhs-continuity.tsv');lookup={float(r['offset_x']):r for r in rhs};zero=lookup[0.];f0=float(zero['xdot']);expected=-f0*(slope_hi-slope_lo)*Tk/C
    result['RHS']={'at_knot_nearest_packed_state':zero,'one_sided_derivatives_by_h':{},'predicted_jump_d_xdot_dx':expected,'predicted_jump_d_xdot_dt':expected*f0}
    fields=('xdot','eta_dot_e','eta_dot_mu','Cstar','Pnet','Pdir','LH','DeltaLnu','Lnu_eq','Lnu_full','Lgamma','sigma_e','sigma_mu','R_e','R_mu')
    for h in (1e-6,1e-7,1e-8,1e-9):
        values={}
        for k in fields:
            left=(float(zero[k])-float(lookup[-h][k]))/h;right=(float(lookup[h][k])-float(zero[k]))/h;values[k]={'left':left,'right':right,'jump':right-left}
        result['RHS']['one_sided_derivatives_by_h'][str(h)]=values
    routes=rows(root/'context-route.tsv');assert all(r['main_route']==r['replay_cold']==r['replay_warm'] for r in routes);result['context_route_values_exact']=len(routes)
    result['context_note']='Fresh disposable state objects share the immutable physics owners and mutable cache hint; replay_cold is a fresh-state label, not a cold thermal cache. Cold cache payload was separately rebuilt and matched all 160 entries.'
    if supplement:
        result['adjacent_RHS']=rows(supplement/'adjacent-rhs.tsv');result['main_right_comparison']=rows(supplement/'main-step-comparison.tsv')
        for p in ('below','at','above'):
            a,b=[r for r in result['adjacent_RHS'] if r['point']==p]
            assert all(a[k]==b[k] for k in fields)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--evidence',type=Path,required=True);p.add_argument('--supplement',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=analyze(a.root,a.evidence,a.supplement);a.output.write_text(json.dumps(r,indent=2,sort_keys=True,allow_nan=False)+'\n');print('PASS diagnostic reduction; no campaign acceptance claim')
