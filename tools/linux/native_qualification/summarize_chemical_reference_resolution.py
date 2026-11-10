#!/usr/bin/env python3
"""Report retained D1 resolution results; never performs an IVP integration."""
import argparse
import csv
import json
from decimal import Decimal as D, localcontext
from pathlib import Path


def generate(result, output):
    points = result['points']
    records = []
    for p in points:
        for i, c in enumerate(p['chemical_state']):
            records.append(dict(point=p['point'], component=c['component'],
                status=p['status'], eta_ref_MeV=p['eta_ref_MeV'][i],
                actual_O2_error_MeV=c['O2_error_MeV'], actual_O2_error_ulp=c['O2_error_ulp'],
                u_discretization_MeV=p['u_ref_discretization_MeV'][i],
                u_offset_MeV=p['u_ref_thermal_MeV'][i],
                uniform_thermal_bound_MeV=p['thermal_sensitivity']['uniform_thermal_bound_MeV'],
                u_supported_MeV=p['u_supported_MeV'][i], recorded_state_radius_MeV=c['recorded_radius_MeV'],
                actual_error_over_radius=D(c['O2_error_MeV'])/D(c['recorded_radius_MeV']),
                supported_uncertainty_over_radius=D(p['u_supported_MeV'][i])/D(c['recorded_radius_MeV']),
                C2_original_ratio=c['ratio'], C2_supported_ratio=p['stronger_chemical_state'][i]['ratio'],
                rho=p['thermal']['rho'], C1=p['C1'], C2=p['C2'], C3=p['C3'], C4_report_only=p['C4']))
    energy = [dict(point=p['point'], **p['stronger_chemical_energy'],
                   reference_MeV=p['chemical_energy']['reference_MeV'],
                   original_u_E_MeV=p['chemical_energy']['u_E_MeV'],
                   original_C3_ratio=p['chemical_energy']['ratio'],
                   recorded_F_MeV=p['chemical_energy']['recorded_F_MeV'],
                   recorded_threshold_MeV=p['chemical_energy']['recorded_threshold_MeV'],
                   original_budget_utilization=p['chemical_energy']['original_budget_utilization']) for p in points]
    original = []
    for p in result['original_failure_analysis']:
        for c in p['components']:
            for index, (a,b) in enumerate(((64,128),(128,256),(256,512))):
                original.append(dict(point=p['point'], component=c['component'], offset=c['offset'],
                    coarse_N=a, fine_N=b, signed_difference_MeV=c['signed_differences_MeV'][index],
                    absolute_difference_MeV=c['absolute_differences_MeV'][index],
                    recorded_state_radius_MeV=c['recorded_radius_MeV'],
                    absolute_difference_over_radius=c['absolute_differences_over_radius'][index],
                    endpoint_Decimal_spacing_MeV=c['endpoint_Decimal_spacing_MeV'],
                    difference_over_Decimal_spacing=c['absolute_differences_over_Decimal_spacing'][index]))
    maxima = {
        'C1_contraction': max(D(c['ratio']) for p in points for cs in p['convergence'].values()
                              for v in cs for c in v['checks'] if c['ratio'] is not None),
        'C2_original': max(D(c['C2_original_ratio']) for c in records),
        'C2_supported': max(D(c['C2_supported_ratio']) for c in records),
        'C3_original': max(D(e['original_C3_ratio']) for e in energy),
        'C3_supported': max(D(e['ratio']) for e in energy),
        'C4_report_only': max(D(v) for p in points for a in p['audit_validation'] for v in a['ratios']),
        'O2_error_ulp': max(D(c['actual_O2_error_ulp']) for c in records),
        'O2_error_MeV': max(D(c['actual_O2_error_MeV']) for c in records),
        'u_discretization_MeV': max(D(c['u_discretization_MeV']) for c in records),
        'u_offset_MeV': max(D(c['u_offset_MeV']) for c in records),
        'uniform_thermal_bound_MeV': max(D(c['uniform_thermal_bound_MeV']) for c in records),
        'rho': max(D(c['rho']) for c in records),
        'LH': max(D(p['thermal_sensitivity']['LH']) for p in points),
    }
    summary = dict(status=result['status'], original_D1_disposition=result['original_D1_disposition'],
        predeclared_execution_commit=result['predeclared_execution_commit'], execution=result['execution'],
        integrations=len(points)*3*4, analytic_tests=17, maxima=maxima,
        rhs_fidelity=result['rhs_fidelity'], chemical_states=records, chemical_energy=energy,
        assumptions=result['assumptions'], exclusions=result['exclusions'])
    (output/'summary.json').write_text(json.dumps(summary,indent=2,default=str)+'\n')
    for name, data in [('chemical-states.tsv', records), ('chemical-energy.tsv', energy),
                       ('original-differences.tsv', original)]:
        with (output/name).open('w') as f:
            writer = csv.DictWriter(f, fieldnames=list(data[0]), delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(data)
    lines = ['# Generated tables — no numerical rerun', '',
             '## Original central differences', '',
             '| Point | Component | N pair | Signed difference (MeV) | Absolute difference / radius |',
             '| ---: | --- | --- | ---: | ---: |']
    for r in original:
        if D(r['offset'])==0:
            lines.append(f"| {r['point']} | {r['component']} | {r['coarse_N']}→{r['fine_N']} | {D(r['signed_difference_MeV']):.12E} | {D(r['absolute_difference_over_radius']):.12E} |")
    lines += ['', '## Aligned state comparison', '',
        '| Point | Component | Actual error (MeV) | Discretization allowance (MeV) | Offset probe (MeV) | Uniform thermal bound (MeV) | Original radius (MeV) | Supported C2 ratio |',
        '| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for r in records:
        vals=[r[k] for k in ('actual_O2_error_MeV','u_discretization_MeV','u_offset_MeV',
                            'uniform_thermal_bound_MeV','recorded_state_radius_MeV','C2_supported_ratio')]
        lines.append(f"| {r['point']} | {r['component']} | "+' | '.join(f'{D(v):.8E}' for v in vals)+' |')
    lines += ['', '## Aligned energy comparison', '',
        '| Point | Actual error (MeV) | Supported reference uncertainty (MeV) | Original U_Echem (MeV) | Supported C3 ratio |',
        '| ---: | ---: | ---: | ---: | ---: |']
    for r in energy:
        vals=[r[k] for k in ('actual_error_MeV','reference_uncertainty_MeV','original_budget_MeV','ratio')]
        lines.append(f"| {r['point']} | "+' | '.join(f'{D(v):.8E}' for v in vals)+' |')
    (output/'TABLES.md').write_text('\n'.join(lines)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    with localcontext() as ctx:
        ctx.prec = 110
        generate(json.loads(args.results.read_text()), args.output)
