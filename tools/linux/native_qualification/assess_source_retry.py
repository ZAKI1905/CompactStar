#!/usr/bin/env python3
"""Read-only assessment of the authorized source retry; never run a solver.

Native collection and raw agreement do not replace the independent uncertainty
requirements. Report every unassessed gate explicitly and retain partial results.
"""
import argparse
import csv
import gzip
import io
import json
import math
import re
from decimal import Decimal as D, localcontext
from pathlib import Path

from assess_stage_a_recovery import assess, lift, rows, sha


def equilibrium_count_roundoff(path, normalization):
    """Bound only the recorded linear-history equilibrium-energy arithmetic.

    B=fl(B0+fl(Bdot*t)); the close endpoint subtraction is exact (Sterbenz).
    Bound each multiplication/addition by half an ulp under the retained binary64
    nearest-rounding, -ffp-contract=off execution. No thermal/state bound follows.
    """
    data = rows(path);a,b = data[0],data[-1]
    assert float(a['t_s'])==0 and float(b['t_s'])==31557600
    assert all(r['Bdot_count_s']==a['Bdot_count_s'] for r in data)
    B0,Bdot,t = float(a['B_count']),float(a['Bdot_count_s']),float(b['t_s'])
    Bend = float(b['B_count'])
    assert B0+Bdot*t==Bend, 'recorded count differs from bound source expression'
    mev_to_erg = 1.602176634e33/1e39
    assert all(r['mu_B_inf_MeV']==a['mu_B_inf_MeV'] for r in data)
    mu = float(a['mu_B_inf_MeV'])
    scale = mev_to_erg*mu
    delta_B = Bend-B0
    recorded = scale*delta_B
    with localcontext() as ctx:
        ctx.prec = 200
        exact_scale = lift(mev_to_erg)*lift(mu)
        exact_delta = lift(Bdot)*lift(t)
        exact_equilibrium = exact_scale*exact_delta
        assert lift(delta_B)==lift(Bend)-lift(B0), 'endpoint subtraction not exact'
        uB = (lift(math.ulp(Bend))+lift(math.ulp(Bdot*t)))/2
        bound = abs(exact_scale)*uB+abs(lift(delta_B))*lift(math.ulp(scale))/2+lift(math.ulp(recorded))/2
        error = lift(recorded)-exact_equilibrium
        assert abs(error)<=bound
        return dict(diagnostics=str(path),time_s=b['t_s'],recorded_delta_B_count=lift(delta_B),
            exact_linear_delta_B_count=exact_delta,recorded_delta_Eeq_erg=lift(recorded),
            exact_linear_delta_Eeq_erg=exact_equilibrium,signed_arithmetic_error_erg=error,
            supported_component_bound_erg=bound,component_bound_over_N=bound/normalization,
            source_count_expression_reproduced=True,endpoint_count_subtraction_exact=True,
            scope='Only linear source count and equilibrium-energy arithmetic, with exact binary64 inputs. '
                  'Does not bound thermal RHS, accumulated main/reference errors, thermal quadrature or full R20. '
                  'Recorded residual is retained; this is not a residual correction.')


def analyze(run):
    result = assess(run)
    source = run/'run/trajectory/source'
    inputs = [run/'stage.log']
    root = Path(__file__).resolve().parents[3]
    prior_manifest = root/'docs/validation/results-first-bnv-pilot/stage-a-recovery/evidence/external-SHA256SUMS'
    expected = {name:digest for digest,name in (line.split('  ',1) for line in prior_manifest.read_text().splitlines())}
    prior = root.parents[1]/'stage-a-recovery-20261109-r1/source-stage/run/trajectory/source'
    regression = []
    names = [f'reference-{i}-{kind}.tsv' for i in range(5) for kind in ('history','observations')]
    names += ['tier-0-cap-0/main.accepted_states.tsv','tier-0-cap-0/endpoint-schedule-repeat/main.accepted_states.tsv']
    for name in names:
        current,old = source/name,prior/name
        if current.exists():
            assert sha(old)==expected[str(old)],('prior input authentication failed',str(old))
            regression.append(dict(file=name,prior_sha256=sha(old),actual_sha256=sha(current),
                                   identical=current.read_bytes()==old.read_bytes()))
    result['retained_trajectory_regression'] = dict(records=regression,
        scope='Byte comparison only; not a main-prefix/reference uncertainty bound.')
    inputs.append(prior_manifest)
    completion = run/'completion.txt'
    result['slurm'] = {}
    if completion.exists():
        inputs.append(completion)
        text = completion.read_text()
        result['slurm'] = dict(re.findall(r'\b(JobId|JobState|ExitCode|RunTime|TimeLimit|StartTime|EndTime|MinMemoryNode|NumCPUs|CPUs/Task)=([^\s]+)',text))
    if not result['failures'] and not result['native_stop_reasons']:
        state = result['slurm'].get('JobState')
        if state == 'TIMEOUT':
            result['status'] = 'INCOMPLETE_WALL_LIMIT'
        elif state in ('FAILED','CANCELLED','OUT_OF_MEMORY','NODE_FAIL'):
            result['status'] = 'INCOMPLETE_SLURM_'+state
    # Necessary lower-bound tests only: omitted positive uncertainties are never
    # set to zero for qualification. Exceeding these limits is already a failure.
    lower_bounds = []
    for main in result['mains']:
        for r in main['R20']:
            norm = lift(r['N'])
            assert norm > 0
            residual = abs(lift(r['residual']))
            reconstruction = lift(r['reconstruction_uncertainty'])
            assert reconstruction >= 0
            lower_bounds.append(dict(main=main['path'],intervals=int(r['intervals']),
                residual_erg=residual,N_erg=norm,reconstruction_uncertainty_erg=reconstruction,
                known_total_lower_bound_ratio=(residual+reconstruction)/norm,
                total_limit=D('2e-4'),reconstruction_lower_bound_ratio=reconstruction/norm,
                reconstruction_event_limit=D('5e-6'),
                fails_necessary_test=(residual+reconstruction)/norm>D('2e-4') or reconstruction/norm>D('5e-6'),
                qualification='UNQUALIFIED: main, event and quadrature uncertainty assessment required'))
    result['R20_necessary_lower_bound_tests'] = lower_bounds
    count_roundoff = []
    for row in lower_bounds:
        diagnostics = Path(row['main']).parent/f"diagnostics-{row['intervals']}.tsv"
        item = equilibrium_count_roundoff(diagnostics,row['N_erg'])
        item['recorded_residual_plus_count_bound_and_reconstruction_over_N'] = (
            row['residual_erg']+row['reconstruction_uncertainty_erg']+item['supported_component_bound_erg'])/row['N_erg']
        count_roundoff.append(item)
    result['equilibrium_count_arithmetic_component'] = count_roundoff
    count_sources = [root/p for p in (
        'tests/bnv/controlled_neutron_sink_fixture.hpp',
        'CompactStar/Physics/Rotochemical/ChemicalImbalanceState.hpp',
        'CompactStar/Units.hpp',
        'CompactStar/Physics/BNV/src/PassiveCheckpointOutput.cpp')]
    result['equilibrium_count_arithmetic_source_hashes'] = {str(p):sha(p) for p in count_sources}
    inputs.extend(count_sources)
    scientific_failures = []
    if result['maximum_R18_ratio'] and result['maximum_R18_ratio']['ratio'] > D('1e-10'):
        scientific_failures.append(dict(gate='R18',**result['maximum_R18_ratio']))
    scientific_failures.extend(dict(gate='R20',**r) for r in lower_bounds if r['fails_necessary_test'])
    result['additional_scientific_failures'] = scientific_failures
    if scientific_failures and not result['failures']:
        result['status'] = 'FAIL_SCIENTIFIC_NECESSARY_TEST'

    packet_rows = []
    for path in sorted(source.glob('tier-*/checkpoints-*/point-*/checkpoint.json')):
        try:
            q = json.loads(path.read_text())
        except json.JSONDecodeError:
            continue
        states = rows(path.parent/'state.tsv')
        diagnostics = rows(path.parent/'diagnostics.tsv')
        energy = next((r for r in diagnostics if r['channel']=='Echem'),None)
        packet_rows.append(dict(path=str(path.parent),time_s=q['time_s'],index=q['checkpoint_index'],
            main_step=q['main_step_ordinal'],state=q['state'],
            propagated_state_radii=[r['propagated_state_radius'] for r in states],
            state_pass=all(r['pass']=='1' for r in states),
            diagnostics_pass=all(r['pass']=='1' for r in diagnostics),
            Echem=energy['O2'] if energy else None,U_Echem=energy['U'] if energy else None,
            Echem_threshold=energy['threshold'] if energy else None,
            Echem_utilization=energy['utilization'] if energy else None,
            high_precision_qualification='NOT_ESTABLISHED'))
        inputs.extend([path,path.parent/'state.tsv',path.parent/'diagnostics.tsv',path.parent/'failure.txt'])
    result['checkpoint_inventory'] = packet_rows
    result['grid_counts'] = {}
    for row in packet_rows:
        key = str(Path(row['path']).parent.relative_to(source))
        result['grid_counts'][key] = result['grid_counts'].get(key,0)+1
    result['completed_scientific_grids'] = sum(len(m['R20']) for m in result['mains'])

    domain = []
    for path in source.glob('tier-*/diagnostics-*.tsv'):
        inputs.append(path)
        for row in rows(path):
            temperature = lift(row['Tinf_K'])
            utilization = lift(row['max_frozen_utilization'])
            depletion = abs(lift(row['DeltaB_over_B0']))
            domain.append(dict(path=str(path),time_s=row['t_s'],Tinf_K=temperature,
                frozen_utilization=utilization,depletion=depletion,
                passed=row['valid_through_sample']=='1' and D('1e6')<=temperature<=D('1e9')
                       and utilization<=1 and depletion<=D('1e-6')))
    result['published_domain_checks'] = dict(samples=len(domain),
        failures=[r for r in domain if not r['passed']],
        maximum_frozen_utilization=max((r['frozen_utilization'] for r in domain),default=None),
        maximum_depletion=max((r['depletion'] for r in domain),default=None))
    if result['published_domain_checks']['failures']:
        result['additional_scientific_failures'].extend(
            dict(gate='DOMAIN',**r) for r in result['published_domain_checks']['failures'])
        result['status'] = 'FAIL_SCIENTIFIC_NECESSARY_TEST'

    reference_path = source/'reference-4-observations.tsv'
    if reference_path.exists():
        reference = rows(reference_path)
        if reference and len(reference)==1921 and float(reference[-1]['t'])==31557600:
            inputs.append(reference_path)
            initial, final = reference[0],reference[-1]
            integrals = {}
            for key in ('P_dir_actual','LH','DeltaLnu','Lnu_full','Lgamma','Lother','Pnet'):
                integrals[key+'_erg'] = sum((lift(b[key])+lift(a[key]))*(lift(b['t'])-lift(a['t']))/2
                                            for a,b in zip(reference,reference[1:]))
            result['diagnostic_only_scientific_context'] = dict(reference='4 (finest independent RKF45)',
                qualification='UNQUALIFIED; descriptive retained outputs, not uncertainty bounds',
                endpoint=final,source_initial_temperature_change_K=lift(final['Tinf'])-lift(initial['Tinf']),
                LH_over_direct_power=lift(final['LH'])/lift(final['P_dir_actual']),
                original_1920_grid_trapezoids=integrals,
                scope='Direct BNV deposition test; does not establish measurable chemical heating. No control comparison.')
            by_time = {float(r['t']):r for r in reference}
            maximum = [None,None,None]
            for item in packet_rows:
                q = json.loads((Path(item['path'])/'checkpoint.json').read_text())
                ref = by_time.get(q['time_s'])
                if ref is None:
                    continue
                for i,key in enumerate(('x','e','mu')):
                    state, reference_state = lift(q['state'][i]),lift(ref[key])
                    error = abs(state-reference_state)
                    output_radius = lift(item['propagated_state_radii'][i])
                    budget = lift(q['main_atol'][i])+lift(q['main_rtol'])*max(abs(state),abs(reference_state))
                    ratio = (error+output_radius)/budget
                    r = dict(component=key,path=item['path'],time_s=q['time_s'],raw_error=error,
                        output_radius=output_radius,budget=budget,known_lower_bound_ratio=ratio)
                    if maximum[i] is None or ratio>maximum[i]['known_lower_bound_ratio']:
                        maximum[i] = r
            result['main_state_necessary_lower_bound_tests'] = dict(maximum_by_component=maximum,
                reference='4; independently justified U_ref is still required and is not set to zero',
                qualification='UNQUALIFIED: raw agreement does not bound main-prefix or accumulated roundoff')
            for r in maximum:
                if r and r['known_lower_bound_ratio']>1:
                    result['additional_scientific_failures'].append(dict(gate='STATE_REFERENCE_LOWER_BOUND',**r))
                    result['status'] = 'FAIL_SCIENTIFIC_NECESSARY_TEST'

    restarts = []
    for path in source.glob('tier-*/restart.txt'):
        lines = path.read_text().splitlines()
        fields = lines[1].split();t = float(fields[0]);state = tuple(map(float,fields[1:]))
        ordinal,rejected = map(int,lines[2].split())
        main = rows(path.parent/'main.accepted_states.tsv')
        selected = [r for r in main if 0<float(r['t_after_s'])<=31557600/2][-1]
        exact = (int(selected['sequence'])==ordinal and float(selected['t_after_s'])==t
                 and state==tuple(float(selected[k]) for k in ('x_state','eta_e_MeV','eta_mu_MeV'))
                 and rejected==int(selected['cumulative_rejected']))
        tail_path = path.parent/'restart-tail/main.accepted_states.tsv'
        tail = rows(tail_path) if tail_path.exists() else []
        restarts.append(dict(path=str(path),ordinal=ordinal,time_s=t,tail_duration_s=31557600-t,
            exact_selected_binary64_state=exact,tail_collected=bool(tail) and float(tail[-1]['t_after_s'])==31557600,
            numerical_accuracy='UNQUALIFIED pending reference, prefix and roundoff assessment'))
        inputs.append(path)
    result['actual_restarts'] = restarts
    result['required_source_matrix'] = dict(references=5,mains=12,schedule_repeats=12,restart_tails=3,
        scientific_grids_per_main=[240,480,960,1920])
    if result['status'] != 'SOURCE_COLLECTION_COMPLETE_UNQUALIFIED':
        result['remaining_accuracy_obligations']['chemical_high_precision_actual_checkpoints'] = 'NOT_REACHED: source collection stopped before full downstream qualification'
    result['original_D1_scope'] = 'Prior 13 intervals only; not substituted for actual-output qualification.'
    result['remaining_accuracy_obligations']['event_coverage'] = (
        'UNQUALIFIED: no actual event crossing in completed main histories'
        if not any(m['actual_events'] for m in result['mains']) else
        'PENDING complete event-neighbor and uncertainty assessment')
    result['resource_counter_caveat'] = ('Slurm accounting disabled; live sstat reports an invalid CPU counter and no RSS. '
        'Allocation wall time is distinct from measured CPU time. See resource-counter-query.txt.')
    result['source_qualified'] = False
    inputs.extend(p for p in source.glob('reference-*.tsv'))
    inputs.extend(p for p in source.glob('tier-*/*.tsv'))
    result['input_hashes'] = {str(p):sha(p) for p in sorted(set(inputs)) if p.is_file()}
    result['assessor_sha256'] = sha(Path(__file__))
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('run',type=Path);p.add_argument('output',type=Path)
    a = p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    with localcontext() as ctx:
        ctx.prec = 110
        result = analyze(a.run)
    (a.output/'assessment.json').write_text(json.dumps(result,indent=2,default=str)+'\n')
    summary = {k:v for k,v in result.items() if k not in ('input_hashes','checkpoint_inventory')}
    summary['full_assessment_sha256'] = sha(a.output/'assessment.json')
    (a.output/'summary.json').write_text(json.dumps(summary,indent=2,default=str)+'\n')
    inventory = result['checkpoint_inventory']
    if inventory:
        buffer = io.StringIO()
        writer = csv.DictWriter(buffer,fieldnames=list(inventory[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(inventory)
        (a.output/'checkpoint-inventory.tsv.gz').write_bytes(gzip.compress(buffer.getvalue().encode(),mtime=0))
    print(result['status'], 'packets',result['actual_checkpoint_packets'],
          'native failures',len(result['failures']),'completed grids',result['completed_scientific_grids'])
