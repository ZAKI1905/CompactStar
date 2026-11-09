#!/usr/bin/env python3
"""Assess persisted budgets; Decimal checks only the existing Echem algebra, not an ODE."""
import csv,json,math,pathlib,sys
from decimal import Decimal,localcontext

run=pathlib.Path(sys.argv[1]);coefficients=pathlib.Path(sys.argv[2]);out=pathlib.Path(sys.argv[3])
def rows(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
folder=run/'replay/output/checkpoint-1'
meta=json.loads((folder/'checkpoint.json').read_text())
diags=rows(folder/'diagnostics.tsv');states=rows(folder/'state.tsv')
failures=[r for r in diags if r['pass']=='0']
result=dict(checkpoint=meta,failed_channels=failures,state_budgets=states,
 observed_reference_agreement={r['channel']:dict(
   range=max(float(r['reference_'+str(i)]) for i in range(5))-min(float(r['reference_'+str(i)]) for i in range(5)),
   maximum_difference_from_O2=max(abs(float(r['reference_'+str(i)])-float(r['O2'])) for i in range(5))) for r in diags},
 interpretation='Local data diagnose estimator contributions; no global accuracy or accumulated-roundoff certification is inferred.')
result['local_witness_settings']=dict(source='PassiveCheckpointOutput.cpp Reconstruct, unchanged library',
 refinement=dict(method='rk8pd',rtol=1e-14,atol=[1e-19,1e-25,1e-25],cap='local span / 2'),
 independent=dict(method='RKF45',rtol=1e-14,atol=[1e-19,1e-25,1e-25],cap='local span / 4'))
rhs=rows(run/'replay/output/local-rhs.tsv')
result['local_rhs_evaluations_by_context']={key:sum(r['context']==key for r in rhs) for key in sorted({r['context'] for r in rhs})}
result['state_event_decomposition']={r['channel']:dict(
 ulp_floor=64*(math.nextafter(abs(float(r['O2'])),math.inf)-abs(float(r['O2']))),
 twice_independent_discrepancy=float(r['U_event'])-64*(math.nextafter(abs(float(r['O2'])),math.inf)-abs(float(r['O2'])))) for r in states}
with localcontext() as ctx:
 ctx.prec=100
 D=lambda x:Decimal.from_float(float(x))
 coeff=rows(coefficients)[0];z00,z01,z10,z11=[D(coeff[k]) for k in ['Z00','Z01','Z10','Z11']]
 determinant=z00*z11-z01*z10
 def energy(e,m):return (z11*e*e-(z01+z10)*e*m+z00*m*m)/(2*determinant)
 e,m=[D(x) for x in meta['O2_state'][1:]];center=energy(e,m)
 trace=[r for r in rows(run/'replay/output/local-diagnostics.tsv') if float(r['t_s'])==meta['time_s']]
 assert len(trace)==10, 'Expected original O1/O2 plus eight existing diagnostic corners'
 assert [float(trace[1][k]) for k in ['x','e','mu']]==meta['O2_state']
 corners=[energy(D(r['e']),D(r['mu'])) for r in trace[2:]]
 energy_row=next(r for r in diags if r['channel']=='Echem')
 precise_propagated=2*max(abs(v-center) for v in corners)
 derivative_e=(z11*e-(z01+z10)*m/2)/determinant
 derivative_m=(z00*m-(z01+z10)*e/2)/determinant
 a_e,a_m=[D(x) for x in meta['O2_atol'][1:]];rtol=D(meta['O2_rtol'])
 # For this checkpoint F_O=D_O2, U_O=2 D_O2, and the corner envelope adds a factor 2.
 floor_linear=4*(abs(derivative_e)*a_e+abs(derivative_m)*a_m)
 relative_linear=4*rtol*(abs(derivative_e*e)+abs(derivative_m*m))
 limit=D(energy_row['threshold'])
 state_radius=[D(states[i]['propagated_state_radius']) for i in [1,2]]
 full_linear=2*(abs(derivative_e)*state_radius[0]+abs(derivative_m)*state_radius[1])
 max_common_atol=(limit-relative_linear)/(4*(abs(derivative_e)+abs(derivative_m)))
 algebra=dict(precision_digits=100,decimal_operands='exact binary64 coefficients/states reconstructed from serialized decimals',
   Echem_MeV=str(center),recorded_Echem_MeV=energy_row['O2'],
   center_binary64_absolute_error_MeV=str(abs(D(energy_row['O2'])-center)),
   recorded_propagated_MeV=energy_row['propagated_state_term'],decimal_propagated_MeV=str(precise_propagated),
   corner_envelope_binary64_difference_MeV=str(abs(D(energy_row['propagated_state_term'])-precise_propagated)),
   maximum_corner_algebra_roundoff_MeV=str(max(abs(D(r['Echem'])-energy(D(r['e']),D(r['mu']))) for r in trace)),
   quadratic_numerator_condition=str((abs(z11*e*e)+abs((z01+z10)*e*m)+abs(z00*m*m))/abs(2*determinant*center)),
   determinant_condition=str((abs(z00*z11)+abs(z01*z10))/abs(determinant)),
   linear_propagation_from_atol_MeV=str(floor_linear),linear_propagation_from_rtol_MeV=str(relative_linear),
   linear_propagation_full_box_MeV=str(full_linear),
   first_checkpoint_common_O2_chemical_atol_upper_bound_MeV=str(max_common_atol),
   bound_caveat='Linearized necessary budget estimate only; omits witness/roundoff/event contributions and is not an accepted tolerance change.')
 result['Echem_algebra']=algebra
 result['Echem_measured_vs_estimated']=dict(
   propagated_fraction=float(energy_row['propagated_state_term'])/float(energy_row['U']),
   atol_linear_fraction_of_propagated=float(floor_linear)/float(energy_row['propagated_state_term']),
   reference_disagreement_utilization=result['observed_reference_agreement']['Echem']['maximum_difference_from_O2']/float(energy_row['threshold']),
   primary_witness_disagreement_utilization=float(energy_row['d_O'])/float(energy_row['threshold']),
   atol_linear_utilization=float(floor_linear/limit),rtol_linear_utilization=float(relative_linear/limit))
 # Optional retained reference: no fresh solve, only a necessary floor-budget forecast.
 if len(sys.argv)>4:
  reference=rows(pathlib.Path(sys.argv[4]));first=reference[1]
  ee,mm=D(first['e']),D(first['mu']);energy_first=energy(ee,mm)
  ge=(z11*ee-(z01+z10)*mm/2)/determinant;gm=(z00*mm-(z01+z10)*ee/2)/determinant
  first_limit=D(.20)*D(1e-11)*energy_first
  first_relative=4*rtol*(abs(ge*ee)+abs(gm*mm))
  bound=(first_limit-first_relative)/(4*(abs(ge)+abs(gm)))
  result['finest_grid_first_checkpoint_necessary_budget']=dict(time_s=float(first['t']),
   common_O2_chemical_atol_upper_bound_MeV=str(bound),current_O2_chemical_atol_MeV=1e-24,
   caveat='Linearized necessary bound from retained reference, before witness/event/roundoff costs; no checkpoint replay or tolerance change performed.')
out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps({k:result[k] for k in ['failed_channels','Echem_algebra','Echem_measured_vs_estimated']},indent=2))
