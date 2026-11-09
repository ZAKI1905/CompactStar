#!/usr/bin/env python3
"""Read retained Stage A evidence and plot diagnostic-only completed source paths.

No solver, RHS, interpolation, or checkpoint reconstruction is called here.
Reference spread is descriptive, never an accumulated-roundoff/error bound.
"""
import argparse
import csv
import hashlib
import json
import math
import pathlib


def rows(path):
    with path.open() as stream:
        return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(stream, delimiter='\t')]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run', type=pathlib.Path)
    parser.add_argument('output', type=pathlib.Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    source = args.run / 'run/trajectory/source'
    log = (args.run / 'stage.log').read_text()
    assert 'STOP native checkpoint diagnostics unresolved' in log
    assert not (args.run / 'run/trajectory/control').exists()
    references = [rows(source / f'reference-{i}-observations.tsv') for i in range(5)]
    for reference in references:
        assert len(reference) == 1921 and reference[0]['t'] == 0 and reference[-1]['t'] == 31557600
        assert all(math.isfinite(v) for row in reference for v in row.values())
        assert all(1e6 <= row['Tinf'] <= 1e9 for row in reference)
    times = [r['t'] for r in references[0]]
    assert all([r['t'] for r in ref] == times for ref in references)
    main_dir = source / 'tier-0-cap-0'
    accepted = rows(main_dir / 'main.accepted_states.tsv')
    repeated = rows(main_dir / 'endpoint-schedule-repeat/main.accepted_states.tsv')
    assert accepted == repeated
    assert (main_dir / 'main.accepted_states.tsv').read_bytes() == (main_dir / 'endpoint-schedule-repeat/main.accepted_states.tsv').read_bytes()
    assert accepted[-1]['t_after_s'] == 31557600
    assert all(r['cumulative_rejected'] == 0 for r in accepted)
    errors = rows(main_dir / 'uncertainties-240.tsv')
    assert len(errors) == 1 and errors[0]['t'] == 0
    assert rows(main_dir / 'R20-grids.tsv') == []
    event_count = sum(r['knot'] >= 0 for r in rows(main_dir / 'events.tsv'))
    ref = references[4]
    spreads = {key: [max(r[j][key] for r in references) - min(r[j][key] for r in references)
                     for j in range(len(times))] for key in ref[0] if key != 't'}
    final = accepted[-1]
    endpoint_error = {key: abs(final[column] - ref[-1][key]) for key, column in
                      [('x', 'x_state'), ('e', 'eta_e_MeV'), ('mu', 'eta_mu_MeV')]}
    summary = dict(
        status='FAIL_CHECKPOINT_DIAGNOSTIC_QUALIFICATION', stage_A_ran=True, stage_A_passed=False,
        stage_B_ran=False, stage_B_passed=None, job=19174, candidate=True,
        label='DIAGNOSTIC ONLY — nonlinear qualification failed; control not run',
        source_references_completed=5, source_main_completed=1, source_schedule_repeats_completed=1,
        control_trajectories_completed=0, restart_tails_completed=0,
        reference_observations_each=1921, accepted_main_steps=len(accepted), main_rejected_steps=0,
        actual_main_events=event_count, schedule_accepted_states_byte_identical=True,
        failure=dict(role='source', tier='BASELINE', policy='uncapped', grid_intervals=240,
                     checkpoint_index=1, time_s=31557600/240,
                     gate='diagnostic U <= 0.20 F for every diagnostic',
                     state_reconstruction='passed before diagnostic refusal',
                     failing_channels='not serialized by unchanged collector'),
        selected_diagnostic_reference='reference-4: RKF45, rtol=1e-15, atol=(1e-20,1e-26,1e-26), cap=15408.984375 s',
        reference_initial=ref[0], reference_final=ref[-1],
        source_temperature_increase_K=ref[-1]['Tinf']-ref[0]['Tinf'],
        maximum_raw_reference_spread={k:max(v) for k,v in spreads.items()},
        main_endpoint_absolute_difference_from_reference4=endpoint_error,
        acceptance=dict(reference_accuracy='UNQUALIFIED: no independent accumulated-roundoff bound or reference interpolation bound',
                        R18='No published dynamic checkpoint packet; static preflight passed',
                        R20='NOT_EVALUATED: original grid file has header only',
                        restart='NOT_RUN', source_minus_control='UNAVAILABLE: control not run',
                        actual_event_coverage='NOT_COVERED: zero events on completed main'),
        resource_usage=json.loads(log.split('RESOURCE_USAGE_JSON ', 1)[1]),
        data_sha256={str(p.relative_to(args.run)): hashlib.sha256(p.read_bytes()).hexdigest()
                     for p in sorted((args.run/'run/trajectory').rglob('*')) if p.is_file()})
    (args.output/'assessment.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    with (args.output/'source-endpoints.tsv').open('w') as f:
        fields=list(ref[0]); writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows([ref[0],ref[-1]])

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    year=np.array(times)/31557600
    fig,axes=plt.subplots(4,1,figsize=(9,12),sharex=True,layout='constrained')
    axes[0].plot(year,[r['Tinf'] for r in ref],label='Source reference 4')
    axes[0].plot([r['t_after_s']/31557600 for r in accepted],
                 [1e8*math.exp(r['x_state']) for r in accepted],'o',ms=3,label='BASELINE main endpoints')
    axes[0].set_ylabel(r'$T_\infty$ [K]');axes[0].ticklabel_format(axis='y',style='plain',useOffset=False)
    axes[0].legend(fontsize=8)
    for key,label in [('e',r'$\eta_e$'),('mu',r'$\eta_\mu$')]:
        axes[1].plot(year,[r[key] for r in ref],label=label)
    axes[1].set_ylabel('Chemical imbalance [MeV]');axes[1].legend()
    for key,label in [('P_dir_actual','BNV direct deposition'),('Lgamma','Photon luminosity'),('Lnu_full','Neutrino luminosity')]:
        axes[2].plot(year,[r[key] for r in ref],label=label)
    axes[2].set_yscale('log');axes[2].set_ylabel('Power at infinity [erg/s]');axes[2].legend(fontsize=8)
    for key,label in [('LH','Reaction heating'),('DeltaLnu','Excess neutrino power')]:
        power=np.array([r[key] for r in ref]);axes[3].plot(year,np.where(power>0,power,np.nan),label=label)
    axes[3].set_yscale('log');axes[3].set_ylabel('Power at infinity [erg/s]');axes[3].legend(fontsize=8)
    axes[3].set_xlabel('Time [Julian years]')
    for ax in axes:ax.grid(alpha=.25)
    fig.suptitle('Mathematical P2 source: completed 1-year paths\nDIAGNOSTIC ONLY — checkpoint qualification failed; control not run',fontsize=11)
    for ext in ('png','pdf'):fig.savefig(args.output/('source-diagnostic-curves.'+ext),dpi=150)
    plt.close(fig)
    fig,axes=plt.subplots(3,1,figsize=(9,8),sharex=True,layout='constrained')
    for ax,key,label in zip(axes,('x','e','mu'),('x spread',r'$\eta_e$ spread [MeV]',r'$\eta_\mu$ spread [MeV]')):
        ax.plot(year,spreads[key],linewidth=.7);ax.set_ylabel(label);ax.grid(alpha=.25)
    axes[-1].set_xlabel('Time [Julian years]')
    fig.suptitle('Five source references: raw maximum-minus-minimum spread\nDIAGNOSTIC ONLY — not an uncertainty bound or convergence PASS',fontsize=11)
    for ext in ('png','pdf'):fig.savefig(args.output/('reference-spread.'+ext),dpi=150)
    plt.close(fig)
    print(json.dumps({k:summary[k] for k in ('status','failure','source_temperature_increase_K','maximum_raw_reference_spread','main_endpoint_absolute_difference_from_reference4')},indent=2))


if __name__ == '__main__':
    main()
