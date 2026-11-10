#!/usr/bin/env python3
"""Offline numerical evidence; unknown accumulated errors never become zero."""
import csv,json,math,pathlib,sys
from decimal import Decimal,localcontext

POINTS=[1,2,4,8,16,32,64,128,256,512,960,1440,1919]
ATOLS=[1e-24,1e-26,5e-27,1e-27,1e-28,1e-29]
def rows(path):
 with pathlib.Path(path).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def envelope(candidate,reference):
 """Measured discrepancy plus independent refinement spread, not a roundoff bound."""
 finest=reference[(0,16)]
 u=2*max(abs(reference[(m,8)]-reference[(m,4)]) for m in [0,1])
 u=max(u,2*max(abs(reference[(m,16)]-reference[(m,8)]) for m in [0,1]),2*abs(reference[(0,16)]-reference[(1,16)]))
 return abs(candidate-finest)+u
def assess(root,coefficients):
 root=pathlib.Path(root);candidates=rows(root/'candidates.tsv');witnesses=rows(root/'witnesses.tsv')
 complete=len(candidates)==len(POINTS)*len(ATOLS) and len(witnesses)==len(POINTS)*6
 result={'schema':'echem-refinement-assessment-v1','complete':complete,'candidates':[],
  'criterion':'unchanged diagnostic U <= 0.20 F; feasibility margin utilization <=0.80',
  'accumulated_roundoff_bound':'UNESTABLISHED','stage_A_authorized_by_assessment':False}
 baseline=next((r for r in candidates if int(r['index'])==8 and float(r['atol'])==1e-24),None)
 result['original_failure_reproduced']=bool(baseline and float(baseline['utilization'])==8.077630769832355)
 coeff=rows(coefficients)[0]
 with localcontext() as ctx:
  ctx.prec=100
  D=lambda x:Decimal.from_float(float(x))
  z00,z01,z10,z11=[D(coeff[k]) for k in ['Z00','Z01','Z10','Z11']]
  det=z00*z11-z01*z10
  def energy(e,m):return (z11*e*e-(z01+z10)*e*m+z00*m*m)/(2*det)
  for atol in ATOLS:
   selected=[r for r in candidates if float(r['atol'])==atol];points=[]
   for r in selected:
    index=int(r['index']);folder=root/f'point-{index}'/f'atol-{atol:.0e}'
    states=rows(folder/'state.tsv');diags=rows(folder/'diagnostics.tsv')
    refs={(int(w['method']),int(w['subdivisions'])):w for w in witnesses if int(w['index'])==index}
    echem=next(x for x in diags if x['channel']=='Echem')
    E=energy(D(r['e']),D(r['mu']));radius=[D(s['propagated_state_radius']) for s in states]
    # Match the already executed binary64 corner construction, then evaluate its
    # existing algebra at 100 digits. This is not a new RHS or trajectory.
    corners=[energy(D(float(r['e'])+se*float(radius[1])),D(float(r['mu'])+sm*float(radius[2]))) for se in [-1,1] for sm in [-1,1]]
    p={'index':index,'time_s':float(r['t_s']),'U_MeV':float(r['U']),'F_MeV':float(r['F']),
       'utilization':float(r['utilization']),'relative_estimated_uncertainty':float(r['U'])/float(r['Echem']),
       'state_qualified':all(s['pass']=='1' for s in states),'all_diagnostics_qualified':all(d['pass']=='1' for d in diags),
       'Echem_MeV':float(r['Echem']),'O1_O2_disagreement_MeV':float(echem['d_O']),
       'propagated_state_MeV':float(echem['propagated_state_term']),
       'witness_or_roundoff_MeV':float(echem['witness_or_roundoff_term']),
       'decimal_energy_MeV':str(E),'energy_algebra_rounding_MeV':str(abs(D(r['Echem'])-E)),
       'decimal_corner_envelope_MeV':str(2*max(abs(c-E) for c in corners)),
       'witness_matrix_complete':len(refs)==6}
    if len(refs)==6:
     p['measured_Echem_envelope_utilization']=envelope(float(r['Echem']),{k:float(w['Echem']) for k,w in refs.items()})/(.2*float(r['F']))
     p['measured_Echem_envelope_over_U']=envelope(float(r['Echem']),{k:float(w['Echem']) for k,w in refs.items()})/float(r['U'])
     p['measured_state_envelope_over_U']=[envelope(float(r[k]),{key:float(w[k]) for key,w in refs.items()})/float(s['propagated_state_radius']) for k,s in zip(['x','e','mu'],states)]
     p['maximum_witness_relative_Echem_difference']=max(abs(float(w['Echem'])-float(r['Echem'])) for w in refs.values())/abs(float(r['Echem']))
     p['witness_Echem_refinement_differences']={str(m):[abs(float(refs[m,8]['Echem'])-float(refs[m,4]['Echem'])),abs(float(refs[m,16]['Echem'])-float(refs[m,8]['Echem']))] for m in [0,1]}
    points.append(p)
   result['candidates'].append({'atol_MeV':atol,'points':points,'complete':len(points)==len(POINTS),
    'maximum_utilization':max((p['utilization'] for p in points),default=None),
    'configured_gate_and_margin_pass':len(points)==len(POINTS) and all(p['state_qualified'] and p['all_diagnostics_qualified'] and p['utilization']<=.8 for p in points),
    'measured_envelope_pass':len(points)==len(POINTS) and all(p.get('measured_Echem_envelope_over_U',math.inf)<=1 and max(p.get('measured_state_envelope_over_U',[math.inf]))<=1 for p in points)})
 result['state_dependence']=[{'index':i,'distinct_O2_states':len({(r['x'],r['e'],r['mu']) for r in candidates if int(r['index'])==i})} for i in POINTS]
 return result

if __name__=='__main__':
 result=assess(sys.argv[1],sys.argv[2]);pathlib.Path(sys.argv[3]).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
 print(json.dumps({k:result[k] for k in ['complete','original_failure_reproduced','accumulated_roundoff_bound','stage_A_authorized_by_assessment']},indent=2))
 for c in result['candidates']:print(c['atol_MeV'],c['maximum_utilization'],c['configured_gate_and_margin_pass'],c['measured_envelope_pass'])
