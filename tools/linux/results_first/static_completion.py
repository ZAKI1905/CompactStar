#!/usr/bin/env python3
"""ADR-0022 retained-data ledger. No solver, RHS, or certificate-writing path."""
import argparse
import copy
import csv
import gzip
import hashlib
import importlib.util
import io
import json
import math
import os
from pathlib import Path
import subprocess

import mpmath as mp
import numpy as np

REPO = Path(__file__).resolve().parents[3]
DOC = REPO / 'docs/validation/results-first-bnv-pilot/static-completion'
GROUPS = ('envelope', 'Tsurface_inf')
START = '30f2ab39572abdf33f51f0ef9609d4fef9b90642'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path, precise=False):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as f:
        return json.load(f, parse_float=str if precise else float)


def authenticate(registry, root=REPO):
    required = dict(status='ACCEPTED_STATIC_PREPARATION_ONLY', starting_commit=START,
                    scope='RF-BNV-P2-NATIVE-BOUNDED-v1', amended_groups=list(GROUPS),
                    pointwise_multiplier=10, worst_case_increment_multiplier=2,
                    paired_trajectories_authorized=False, historical_certificate_promoted=False,
                    new_stellar_calculations_authorized=False, missing_terms_must_refuse=True)
    if any(registry.get(k) != v for k, v in required.items()):
        raise ValueError('wrong or broadened static authority')
    for name, digest in registry['files'].items():
        if sha(root / name) != digest:
            raise ValueError('changed authority input: ' + name)
    spec = importlib.util.spec_from_file_location('original_cref', root/'tools/linux/results_first/cref_completion.py')
    original = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(original)
    original.authenticate(read(root/'docs/validation/results-first-bnv-pilot/preparation-registry.json'))
    return original


def point_allowance(name, detail):
    measured = max(abs(v) for row in detail['signed_error'] for v in row)
    return 10 * (9.1e-12 if name == 'Ltilde_Me' else measured)


def revised_partial(A, b, e, fit, *, origin=False, absolute=False):
    if any(not math.isfinite(v) or v < 0 for v in (A, b, e, fit)) or b + e >= 1:
        raise ValueError('invalid budget or unseparated denominator')
    added = (0 if origin else 2*e) + fit
    return A + added if absolute else (A*(1-b)+added)/(1-b-e)


def thermal_projection(x, thermal, comparison, digits):
    """Independent high-precision normal equations, with historical binary64 weights."""
    result = {}
    with mp.workdps(digits):
        xx = [mp.mpf(float(v)) for v in x]
        L = mp.mpf('0.000001')
        V = mp.matrix([[1, v/L] for v in xx])
        Vf = V[10:21, :]
        for name, value_key, error_key in (
                ('envelope', 'Tsurface_local_K', None),
                ('Tsurface_inf', 'Tsurface_inf_K', 'u_Tsurface_inf_K')):
            d = comparison[name]
            fnorm = xfit = efit = null = affine = mp.mpf(0)
            for component in range(13):
                rows = thermal[component::13]
                if [int(r['index']) for r in rows] != list(range(21)):
                    raise ValueError('thermal row alignment')
                y = [float(r[value_key]) for r in rows]
                u = [float(r[error_key]) if error_key else 0. for r in rows]
                sigma = [max(u[i], 64*np.finfo(float).eps*max(1., abs(y[i]))) for i in range(10, 21)]
                W = mp.diag([1/mp.mpf(s)**2 for s in sigma])
                H = V * (Vf.T*W*Vf)**-1 * Vf.T*W
                M = mp.eye(21)
                for i in range(21):
                    for j in range(11):
                        M[i, j+10] -= H[i, j]
                fnorm = max(fnorm, max(mp.fsum(abs(M[i, j]) for j in range(21)) for i in range(21)))
                z = mp.matrix([mp.mpf(float(row[component])) for row in d['interpolated']])
                error = mp.matrix([mp.mpf(float(row[component])) for row in d['signed_error']])
                xfit = max(xfit, max(abs(v) for v in M*z))
                efit = max(efit, max(abs(v) for v in M*error))
                null = max(null, max(abs(v) for v in M*mp.matrix([1]*21)),
                           max(abs(v) for v in M*mp.matrix(xx)))
                e = mp.mpf(point_allowance(name, d))
                varying = mp.matrix([e*(1+2*v/L) for v in xx])
                affine = max(affine, max(abs(v) for v in M*varying))
            curvature = fnorm*mp.mpf(d['K'])*L**2/8
            F = curvature+xfit
            if null > mp.mpf(10)**(-digits+8) or affine > mp.mpf(10)**(-digits+8):
                raise ValueError('fit fails affine identities')
            if abs(float(fnorm)-d['fit_operator_norm']) > 1e-13:
                raise ValueError('historical fit norm mismatch')
            if abs(float(efit)-d['fit_deterministic']) > 1e-18:
                raise ValueError('historical deterministic fit mismatch')
            result[name] = dict(f=float(fnorm), K=d['K'], L=1e-6, fit_reference=float(xfit),
                                measured_fit_error=float(efit), curvature_fit=float(curvature),
                                F=math.nextafter(float(F), math.inf),
                                affine_identity_residual=mp.nstr(null, 12),
                                varying_full_allowance_identity=mp.nstr(affine, 12),
                                digits=digits)
    return result


def tangent_witness(paths, stored, digits):
    result = {}
    with mp.workdps(digits):
        for target in (0, 10, 20):
            raw = {(g, sign): read(paths[f't{target}_n{g}_{sign}'], True)
                   for g in (512, 1024, 2048) for sign in (-1, 1)}
            for name, species in (('t_n', 0), ('t_e', 2), ('t_mu', 3)):
                a = {g: mp.mpf(raw[g, 1]['N'][species])-mp.mpf(raw[g, -1]['N'][species])
                     for g in (512, 1024, 2048)}
                d = {g: mp.mpf(raw[g, 1]['primary'][0])-mp.mpf(raw[g, -1]['primary'][0]) for g in a}
                t = {g: a[g]/d[g] for g in a}
                delta1, delta2 = t[1024]-t[512], t[2048]-t[1024]
                joint = ((a[2048]-a[1024])-t[1024]*(d[2048]-d[1024]))/d[2048]
                rec = stored['records'][f'{target}:{name}']
                order = mp.log(abs(delta1/delta2), 2)
                if delta1*delta2 <= 0 or not mp.mpf('3.5') <= order <= mp.mpf('4.5'):
                    raise ValueError('tangent asymptotic gate')
                if max(abs(t[g]-mp.mpf(rec['t'][str(g)])) for g in t) > mp.mpf('1e-55'):
                    raise ValueError('raw tangent does not reproduce retained ratio')
                if abs(joint-delta2) > mp.mpf(10)**(-digits+8):
                    raise ValueError('correlated tangent identity')
                if not 3 <= float(rec['stencil_ratio']) <= 5:
                    raise ValueError('tangent width gate')
                result[f'{target}:{name}'] = dict(t=float(t[2048]), tail=float(abs(delta2)/15),
                    order=float(order), denominator=float(d[2048]), denominator_relative_allowance=1e-14,
                    separated=bool(abs(d[2048])*mp.mpf('1e-14') < abs(d[2048])),
                    correlated_identity_residual=mp.nstr(abs(joint-delta2), 12), digits=digits)
    return result


def term(name, value, kind, conditions, derivation, sources, **extra):
    if value is not None and (not math.isfinite(value) or value < 0):
        raise ValueError('invalid error contribution: '+name)
    if value is None and kind != 'unresolved':
        raise ValueError('missing value must remain unresolved')
    return dict(name=name, value=value, evidence_type=kind, conditions=conditions,
                derivation=derivation, source_ids=sources, **extra)


def finish_ledger(historical_drift, terms, baseline, *, absolute=False):
    """Unknown is not zero. Return known partial and a distinct nullable full total."""
    if any(t['value'] is not None and (not math.isfinite(t['value']) or t['value'] < 0) for t in terms):
        raise ValueError('invalid ledger value')
    known = historical_drift + math.fsum(t['value'] for t in terms if t['value'] is not None)
    b = 0. if absolute else math.fsum(v['value'] for v in baseline if v['value'] is not None)
    if not 0 <= b < 1:
        raise ValueError('unseparated baseline')
    missing = [t['name'] for t in terms if t['value'] is None]
    missing += ['baseline:'+t['name'] for t in baseline if t['value'] is None and not absolute]
    partial = known if absolute else known/(1-b)
    normalization = 0. if absolute else None if missing else known*b/(1-b)
    return dict(known_numerator=known, known_baseline=b, partial_budget=partial,
                normalization=normalization, missing=missing,
                full_budget=None if missing else partial)


def retained_term(t):
    name = t['name']
    conditions = []
    if t['evidence_type'] == 'validated_numerical_estimate':
        conditions = ['S', 'C-TAN' if name == 'tangent_stencil' else 'C-REF']
    elif t['evidence_type'] == 'conditional_estimate':
        if name in ('reference_interpolation', 'physical_inter_sample_remainder'):
            conditions = ['C-CONT', 'C-CUT']
        elif name == 'X2_common_B_alignment':
            conditions = ['C-CONT', 'C-REF']
        elif name == 'tangent_stencil':
            conditions = ['C-TAN']
        elif name == 'fit_reference_inputs':
            conditions = ['C-REF', 'C-TAN', 'C-CONT', 'C-INPUT']
        elif name == 'equilibrium_profile_definition':
            conditions = ['C-MUN']
        else:
            raise ValueError('unmapped conditional: '+name)
    return term(name, t['value'], t['evidence_type'], conditions,
                t.get('derivation') or 'Retained finite entry; original provenance and estimate unchanged.',
                ['matrix'], historical_dependencies=t['dependencies'])


def assemble(matrix, comparison, m1, m3, m4, tangents, diagnostics, projection):
    rows = []
    partial_violations = []
    fullstage_expectations = []
    missing_cut = ('C-CUT fixes delta ln Ps <=2.1e-12 but no validated observable-specific '
                   'finite-perturbation response at fixed B is retained. Along-path K is not that input.')
    for old in matrix['rows']:
        n, i = old['quantity'], old['index']
        h = old['historical']
        new = copy.deepcopy(old)
        if n == 'metric_structure':
            new['static_completion'] = dict(status='CONSUMER_SCOPE_ONLY', full_utilization=None,
                nodewise_qualified=False, reason='Accepted P2 consumer limitation; original nodewise unknowns preserved.')
            rows.append(new)
            continue
        if n in ('P2_average', 'species_support'):
            if h['drift_bound'] != 0 or any(t['value'] not in (None, 0) for t in old['terms']):
                raise ValueError('identity/topology historical mismatch')
            new['static_completion'] = dict(status='EXACT_IDENTITY' if n == 'P2_average' else 'CONDITIONALLY_QUALIFIED',
                full_utilization=0., partial_utilization=0., baseline=[], full_budget=0., missing=[],
                terms=[term('exact_P2' if n == 'P2_average' else 'support_continuation', 0.,
                            'analytic_bound' if n == 'P2_average' else 'conditional_estimate',
                            [] if n == 'P2_average' else ['C-SUPPORT', 'C-CONT'],
                            'Exact P2 null.' if n == 'P2_average' else 'No crossing under accepted topology and positive-B slope.',
                            ['matrix', 'continuation'])])
            rows.append(new)
            continue
        d = comparison[n]
        absolute = n.startswith('t_')
        base = 'mu_B' if n == 'mu_n_profile' else n
        e = point_allowance(n, d)
        f = d['fit_operator_norm']
        A = h['drift_bound'] + old['finite_added_estimates']
        F = projection[n]['F'] if n in GROUPS else f*e
        partial = revised_partial(A, d['base_error_estimate'], e, F, origin=i == 0, absolute=absolute)
        delta = partial/h['threshold'] - old['partial_utilization']
        if delta > 1e-3:
            partial_violations.append(dict(index=i, quantity=n, increase=delta,
                                           amended_partial_utilization=partial/h['threshold']))
        entries = [retained_term(t) for t in old['terms']
                   if t['value'] is not None and t['name'] != 'full_numerator_normalization']
        baseline = [] if absolute else [term('retained_base', d['base_error_estimate'], 'conditional_estimate',
                    ['C-REF', 'C-INPUT'], 'Retained reference error at B0, in original component scales.', ['comparison'])]
        coefficient = f + (0 if i == 0 else 2)
        # mu_n_profile shares the scalar mu_B identity, rescaled to its own nodewise scales.
        scale = (max(abs(a) for row in d['reference_values'] for a in row) /
                 max(abs(a) for row in comparison['mu_B']['reference_values'] for a in row)) if n == 'mu_n_profile' else 1.
        tail = (max(float(tangents['records'][f'{t}:{n}']['fine_tail']) for t in (0, 10, 20)) if absolute else
                max(m1['reference'][str(t)][base]['integration'] for t in (0, 10, 20))*scale)
        stage_joint = 15*tail
        for name in ('reference_stage_remainder', 'reference_arithmetic_remainder'):
            u = stage_joint/2
            entries.append(term(name, coefficient*u, 'conditional_estimate', ['C-REF', 'S'],
                                'Half of the single finest-increment stage+arithmetic allowance; paired+fit response.',
                                ['m1', 'tangents'], pointwise_allowance=u))
            if not absolute:
                baseline.append(term(name, u, 'conditional_estimate', ['C-REF', 'S'],
                                     'Same joint allowance at B0.', ['m1']))
        entries.append(term('unseen_reference_bias_increment', 0. if i == 0 else 2*e, 'conditional_estimate',
                            ['C-REF'], 'Full ten-times pointwise allowance; worst-case paired increment.',
                            ['comparison'], pointwise_allowance=e))
        entries.append(term('total_reference_fit_response' if n in GROUPS else 'unseen_reference_fit_response',
                            F, 'conditional_estimate', ['C-REF', 'C-CONT', 'ADR-0022'] if n in GROUPS else ['C-REF'],
                            'f*K*L^2/8+norm(M z) bounds TOTAL reference fit; prior finite entries retained conservatively.'
                            if n in GROUPS else 'Unchanged all-row operator norm times full pointwise allowance.',
                            ['comparison', 'projection']))
        if not absolute:
            baseline.append(term('unseen_reference_bias', e, 'conditional_estimate', ['C-REF'],
                                 'Full original pointwise allowance also enters the denominator.', ['comparison']))
        entries.append(term('reference_parameter_remainders', 0., 'conditional_estimate', ['C-REF', 'S'],
                            'Accepted convergence models make charged EOS/centre estimates complete; avoid duplicate remainder.', ['m1']))
        point_extra = 0.
        if absolute:
            stencil = max(float(tangents['records'][f'{t}:{n}']['stencil_tail_observation']) for t in (0, 10, 20))
            cross = max(abs(float(v['residual'])) for v in diagnostics['tangents'][n]['cross_geometry'].values())
            eos = max(m1['reference'][str(t)][n]['eos_equilibrium'] for t in (0, 10, 20))
            center = max(m1['reference'][str(t)][n]['center'] for t in (0, 10, 20))
            U = tail+2*stencil+cross+eos+center
            correction = U*1e-14/(1-1e-14)
            point_extra = stencil+cross+correction
            for name, u, why in [('tangent_second_width_tail', stencil, 'Second narrow-width tail required by C-TAN.'),
                                 ('tangent_moving_cut_path', cross, 'Accepted cross-geometry fixed-cut/path residual.'),
                                 ('tangent_denominator_remainder', correction, 'Conservative U*delta/(1-delta), delta=1e-14; ratio-level J=U*abs(D).')]:
                entries.append(term(name, coefficient*u, 'conditional_estimate', ['C-TAN'], why,
                                    ['tangents', 'diagnostics', 'm1'], pointwise_allowance=u))
        gaps = m3['curvatures'][base]['intervals']
        measured = max(abs(v) for row in d['signed_error'] for v in row)
        nonlinear = max((g['coordinate_estimate']/g['achieved_spacing'])**2*measured for g in gaps)
        entries.append(term('fit_nonlinear_input_remainder', nonlinear, 'conditional_estimate', ['C-INPUT'],
                            'Accepted second-order design range; coordinate first-order and historical weights remain separate.', ['m3', 'comparison']))
        entries.append(term('moving_cut_input_remainder', None, 'unresolved', ['C-CUT'], missing_cut,
                            ['m1', 'm3', 'continuation', 'diagnostics'], delta_log_pressure=2.1e-12))
        if not absolute:
            baseline.append(term('moving_cut_at_B0', None, 'unresolved', ['C-CUT'], missing_cut, ['continuation']))
        # Endpoint slope bound includes full known endpoint error, including C-REF.
        sliver = max(0., 1e-6+float(m4['achieved_x'][-1]))
        width = abs(m4['achieved_x'][20]-m4['achieved_x'][10])
        secant = max(abs(a-b)/width for a, b in zip(d['reference_values'][2], d['reference_values'][1]))
        known = tail + stage_joint + point_extra + e + max(d['root_each'].values())
        known += max(m1['reference'][str(t)][base]['eos_equilibrium'] + m1['reference'][str(t)][base]['center']
                     + m1['reference'][str(t)][base]['stencil'] for t in (0, 10, 20))*scale
        endpoint = (secant+2*known/width+d['K']*width/2)*sliver+d['K']*sliver**2/2 if i == 20 else 0.
        entries.append(term('endpoint_sliver_known_response', endpoint, 'conditional_estimate', ['C-CONT', 'C-REF'],
                            'Last secant + two known endpoint bands/spacing + K*spacing/2; derivative times sliver + K*sliver^2/2.',
                            ['m4', 'comparison', 'm1'], endpoint_sliver=sliver))
        if i == 20:
            entries.append(term('endpoint_sliver_cut_response', None, 'unresolved', ['C-CUT', 'C-CONT'],
                                'Endpoint derivative still depends on missing finite-cut endpoint bands.', ['continuation']))
        totals = finish_ledger(h['drift_bound'], entries, baseline, absolute=absolute)
        norm = term('full_numerator_normalization', totals['normalization'],
                    'unresolved' if totals['normalization'] is None else 'analytic_bound', ['C-INPUT'],
                    'Full numerator and full baseline required; unknown cut contribution is never zero.', ['matrix', 'comparison'])
        status = 'UNQUALIFIED' if totals['missing'] else 'CONDITIONALLY_QUALIFIED' if totals['full_budget'] <= h['threshold'] else 'FAIL'
        if totals['partial_budget'] > h['threshold']:
            status = 'FAIL_PARTIAL_BUDGET'
        new['static_completion'] = dict(status=status, terms=entries+[norm], baseline=baseline,
            **totals, full_utilization=None if totals['full_budget'] is None else totals['full_budget']/h['threshold'],
            partial_utilization=totals['partial_budget']/h['threshold'],
            amended_cref_partial_utilization=partial/h['threshold'],
            amended_cref_partial_status='SATISFIES_NECESSARY_TEST' if partial <= h['threshold'] else 'FAIL',
            original_cref_failure_preserved=n in GROUPS,
            qualification_wording='conditional on named estimates and assumptions; not a rigorous enclosure')
        increase = totals['partial_budget']/h['threshold']-old['partial_utilization']
        if increase > 1e-3:
            fullstage_expectations.append(dict(index=i, quantity=n, increase=increase))
        rows.append(new)
    stripped = copy.deepcopy(rows)
    for r in stripped:
        del r['static_completion']
    if stripped != matrix['rows']:
        raise ValueError('historical matrix changed')
    return rows, partial_violations, fullstage_expectations


def run(output):
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('EKU Slurm execution required')
    original = authenticate(read(DOC/'acceptance-registry.json'))
    manifest = read(DOC/'inputs.json')
    paths = {}
    for key, rec in manifest.items():
        p = Path(rec['path'])
        p = p if p.is_absolute() else REPO/p
        if sha(p) != rec['sha256']:
            raise ValueError('changed retained input: '+key)
        paths[key] = p
    data = {key: read(paths[key]) for key in ('matrix', 'comparison', 'm1', 'm3', 'm4', 'tangents', 'diagnostics', 'continuation')}
    old = original.assess(data['matrix'], data['comparison'])['restricted_summary']
    saved = read(paths['original_summary'])
    if old['failed_rows'] != saved['failed_rows'] or old['review_expectation_failures'] != saved['review_expectation_failures']:
        raise ValueError('original failure evidence not reproduced')
    if len(old['failed_rows']) != 42:
        raise ValueError('original 42 failures missing')
    with paths['thermal'].open() as f:
        thermal = list(csv.DictReader(f, delimiter='\t'))
    projection = [thermal_projection(data['m4']['achieved_x'], thermal, data['comparison'], digits) for digits in (70, 90)]
    tangent = [tangent_witness(paths, data['tangents'], digits) for digits in (70, 90)]
    for name in GROUPS:
        for key in ('f', 'F', 'fit_reference', 'measured_fit_error', 'curvature_fit'):
            if projection[0][name][key] != projection[1][name][key]:
                raise ValueError('projection precision instability')
    for name in tangent[0]:
        for key in ('t', 'tail', 'order', 'denominator'):
            if tangent[0][name][key] != tangent[1][name][key]:
                raise ValueError('tangent precision instability')
    rows, violations, fullstage = assemble(
        data['matrix'], data['comparison'], data['m1'], data['m3'], data['m4'], data['tangents'], data['diagnostics'], projection[1])
    if len(violations) != 353:
        raise ValueError('353 amended-partial expectation violations not reproduced')
    amended_failures = [r for r in rows if r['quantity'] in GROUPS and r['static_completion']['amended_cref_partial_status'] != 'SATISFIES_NECESSARY_TEST']
    if amended_failures:
        raise ValueError('amended temperature necessary test failed')
    if any(r['static_completion'].get('missing') and r['static_completion']['full_utilization'] is not None for r in rows):
        raise ValueError('unknown error was promoted to full qualification')
    counts = {s: sum(r['static_completion']['status'] == s for r in rows)
              for s in sorted({r['static_completion']['status'] for r in rows})}
    incomplete = [dict(index=r['index'], quantity=r['quantity'], missing=r['static_completion'].get('missing', []))
                  for r in rows if r['static_completion']['status'] in ('UNQUALIFIED', 'FAIL', 'FAIL_PARTIAL_BUDGET')]
    group_summary = {}
    for name in sorted({r['quantity'] for r in rows}):
        rr = [r['static_completion'] for r in rows if r['quantity'] == name]
        group_summary[name] = dict(statuses=sorted({r['status'] for r in rr}),
            maximum_supported_partial_utilization=max((r.get('partial_utilization', 0) for r in rr), default=0) if name != 'metric_structure' else None,
            maximum_full_utilization=None if any(r['full_utilization'] is None for r in rr) else max(r['full_utilization'] for r in rr),
            amended_cref_partial_maximum=max((r.get('amended_cref_partial_utilization', 0) for r in rr), default=0) if name in data['comparison'] else None)
    summary = dict(status='BLOCKED_UNSUPPORTED_CUT_RESPONSE' if incomplete else 'STATIC_QUALIFICATION_COMPLETE',
        execution_commit=subprocess.check_output(['git', '-C', str(REPO), 'rev-parse', 'HEAD'], text=True).strip(),
        slurm_job_id=os.environ['SLURM_JOB_ID'], counts=counts, rows=len(rows), original_42_failures_preserved=True,
        amended_42_partial_satisfy=True, amended_partial_expectation_violation_count=len(violations),
        completion_partial_expectation_violation_count=len(fullstage), groups=group_summary, blockers=incomplete,
        operational_certificate_issued=False, historical_certificate_promoted=False, trajectories_authorized=False,
        stellar_calculations=0, trajectories_run=0,
        scientific_limit='Finite-cut input response unsupported; partial sums are not full qualification.')
    output.mkdir(exist_ok=False)
    def write(name, value):
        raw = (json.dumps(value, indent=2, allow_nan=False)+'\n').encode()
        (output/name).write_bytes(gzip.compress(raw, mtime=0) if name.endswith('.gz') else raw)
    write('summary.json', summary)
    write('ledger.json.gz', dict(summary=summary, rows=rows, original_cref_summary=saved,
                               accepted_registry=read(DOC/'acceptance-registry.json'), inputs=manifest,
                               derived_sources={'projection': 'precision-witness.json'}))
    write('precision-witness.json', dict(projection=projection, tangents=tangent, precision_repeat='identical binary64 reporting values'))
    write('amended-expectation-violations.json', violations)
    write('completion-expectation-violations.json', fullstage)
    write('checks.json', dict(authority=True, all_inputs_authenticated=True, original_42_reproduced=True,
        original_355_reproduced=True, amended_353_reproduced=True, historical_rows_unchanged=True,
        temperature_fit_affine_identities=True, temperature_fit_norms=True, precision_70_90_agree=True,
        raw_tangent_ratios=True, correlated_tangent_identities=True, no_full_pass_with_missing_cut=True))
    stream = io.StringIO()
    writer = csv.writer(stream, delimiter='\t')
    writer.writerow(['index', 'quantity', 'status', 'threshold', 'historical_utilization', 'amended_cref_partial',
                     'supported_partial_utilization', 'full_utilization', 'missing'])
    for r in rows:
        s = r['static_completion']
        writer.writerow([r['index'], r['quantity'], s['status'], r['historical']['threshold'], r['historical']['utilization'],
            s.get('amended_cref_partial_utilization'), s.get('partial_utilization'), s['full_utilization'], ';'.join(s.get('missing', []))])
    (output/'matrix.tsv.gz').write_bytes(gzip.compress(stream.getvalue().encode(), mtime=0))
    print(json.dumps({k: summary[k] for k in ('status', 'counts', 'rows', 'amended_partial_expectation_violation_count')}))
    return 2 if incomplete else 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    raise SystemExit(run(args.output))
