#!/usr/bin/env python3
"""Predeclared D1 knot-aligned offline reference; scientific candidate.

The original verifier and its FAIL evidence are immutable. Reuse only its
source-bound mathematical model, input checks and reporting definitions.
Never import the production solver or the roundoff-audit recursion.
"""
import argparse
import collections
import json
import os
import platform
import resource
import subprocess
import sys
import time
import traceback
from decimal import Decimal, localcontext, ROUND_CEILING, ROUND_FLOOR
from pathlib import Path

from verify_chemical_reference import (
    AUDIT, Chemistry, DELTA_X, GateFailure, NS, PI, POINTS, PRECISION,
    REFINEMENT, RHO_LIMIT, authenticate, convergence, corner_uncertainty,
    digest, fidelity, lift, make_path, require, rows, save, ulp,
)

START = '871dbd77230f9e21b54ea0129473fe658236fd56'
REVIEW_SHA = '46201727575270055d237178b2ad0433ddc22837c0cc1fa646a8d34eddcb0066'
RESOLUTION = REFINEMENT / 'd1-resolution'
KS = (8, 16, 32, 64)


def segment_value(a, b, s):
    """Hermite at normalized s; durations stay in exact local seconds."""
    h = b[0]-a[0]
    return ((2*s**3-3*s**2+1)*a[1] + (s**3-2*s**2+s)*h*a[2]
            + (-2*s**3+3*s**2)*b[1] + (s**3-s**2)*h*b[2])


def aligned_rk4(rhs, nodes, initial, k, offset=Decimal(0)):
    """K classical RK4 steps on EACH segment, including tiny terminal pieces.

    rhs takes x, eta, and returns MeV/s. The grid uses integer stage indices
    divided by 2*K so its final stage is exactly s=1, without rounded t+h.
    No IVP integration or evaluation of x outside its owning segment occurs.
    """
    require(k in KS, 'method', f'undeclared steps per segment: {k}')
    y = tuple(initial)
    for a, b in zip(nodes, nodes[1:]):
        duration = b[0]-a[0]
        require(duration > 0, 'method', 'nonpositive segment duration')
        h = duration/k
        xs = [segment_value(a, b, Decimal(j)/(2*k))+offset for j in range(2*k+1)]
        for j in range(k):
            k1 = rhs(xs[2*j], y)
            k2 = rhs(xs[2*j+1], tuple(v+h*d/2 for v, d in zip(y, k1)))
            k3 = rhs(xs[2*j+1], tuple(v+h*d/2 for v, d in zip(y, k2)))
            k4 = rhs(xs[2*j+2], tuple(v+h*d for v, d in zip(y, k3)))
            y = tuple(v+h*(a1+2*a2+2*a3+a4)/6
                      for v, a1, a2, a3, a4 in zip(y, k1, k2, k3, k4))
    return y


def sensitivity_bound(model, nodes, initial):
    """Positive upper bounds for any uniform-band thermal perturbation.

    Implements the polynomial/first-exit/Gronwall derivation in PREDECLARATION.
    Decimal exp rounds to nearest regardless of context; next_plus encloses it.
    Floor-rounded positive denominators ensure division is an upper bound.
    This does NOT bound RK4 discretization or prove the thermal band itself.
    """
    with localcontext() as ctx:
        ctx.prec = PRECISION
        ctx.rounding = ROUND_FLOOR
        denominators = [model.kberg*11513*model.kb**(2*k+1)*PI**(2*k+2) for k in range(4)]
    with localcontext() as ctx:
        ctx.prec = PRECISION
        ctx.rounding = ROUND_CEILING
        H = nodes[-1][0]-nodes[0][0]
        y0 = max(abs(v) for v in initial)
        # Bound the defining source expression, not a possibly rounded signed sum.
        q = max(sum(abs(model.z[i][j])*abs(model.sigma[j]).next_plus() for j in range(2))
                for i in range(2))
        Y = 2*(y0+H*q)
        x_upper = max(n[1] for n in nodes) + max(
            (b[0]-a[0])*max(abs(a[2]), abs(b[2]))/3 for a, b in zip(nodes, nodes[1:]))
        T = Decimal('1e8')*(x_upper+DELTA_X).exp().next_plus()
        aa = [[abs(model.luminosity[j])*c/denominators[k]
               for k, c in enumerate((14680, 7560, 840, 24))] for j in range(2)]
        rate = [sum(aa[j][k]*Y**(2*k+1)*T**(6-2*k) for k in range(4)) for j in range(2)]
        deta = [sum((2*k+1)*aa[j][k]*Y**(2*k)*T**(6-2*k) for k in range(4)) for j in range(2)]
        dx = [sum((6-2*k)*aa[j][k]*Y**(2*k+1)*T**(6-2*k) for k in range(4)) for j in range(2)]

        def norm(v):
            return max(sum(abs(model.z[i][j])*v[j] for j in range(2)) for i in range(2))

        beta, L, G = norm(rate), norm(deta), norm(dx)
        exit_upper = y0+H*(q+beta)
        bound = DELTA_X*G*H*(L*H).exp().next_plus()
        return dict(status='PASS' if exit_upper < Y else 'FAIL',
                    chemical_domain_abs_MeV=Y, first_exit_upper_MeV=exit_upper,
                    x_upper=x_upper, T_upper_K=T, source_bound_MeV_s=q,
                    reaction_bound_MeV_s=beta, lipschitz_per_s=L,
                    thermal_derivative_MeV_s=G, LH=L*H,
                    uniform_thermal_bound_MeV=bound,
                    scope='All time-dependent |delta_x(t)|<=1e-12; assumes uniform band; excludes RK4 error.')


def verify_initial(checkpoint):
    result = []
    for name, value in zip(('eta_e', 'eta_mu'), checkpoint['main_left_state'][1:]):
        value = float(value)
        numerator, denominator = value.as_integer_ratio()
        rational = Decimal(numerator)/Decimal(denominator)
        exact = lift(value)
        require(rational == exact, 'initial-state', f'exact double conversion mismatch {name}')
        result.append(dict(component=name, binary64_hex=value.hex(), numerator=numerator,
                           denominator=denominator, exact_Decimal_MeV=exact, status='PASS'))
    return result


def old_diagnosis(old, model):
    """Only analyze retained outputs; no old-grid IVP rerun."""
    nodes = [tuple(Decimal(v) for v in row) for row in old['thermal']['nodes_tau_s_x_dx']]
    span = nodes[-1][0]
    eta = tuple(map(Decimal, old['eta_ref_MeV']))
    radii = [Decimal(r['recorded_radius_MeV']) for r in old['chemical_state']]
    segments = [b[0]-a[0] for a, b in zip(nodes, nodes[1:])]
    grids = {}
    for n in NS:
        h = span/n
        knots = []
        for t, x, dx in nodes[1:-1]:
            grid_index = int((t/h).to_integral_value())
            offset = t-span*grid_index/n
            knots.append(dict(tau_s=t, nearest_step=grid_index, offset_s=offset,
                              offset_over_h=offset/h, aligned=(offset == 0)))
        grids[n] = dict(h_s=h, knot_offsets=knots, terminal_segment_over_h=segments[-1]/h,
                        final_midpoint_inside_terminal_segment=(h/2 < segments[-1]))
    a, b = nodes[-3:-1]
    extrapolated = segment_value(a, b, (span-a[0])/(b[0]-a[0]))
    components = []
    for offset, values in old['integrations_eta_MeV'].items():
        vals = {int(k): tuple(map(Decimal, v)) for k, v in values.items()}
        center = vals[512]
        actual = model.rhs(nodes[-1][1]+Decimal(offset), center)
        extrap = model.rhs(extrapolated+Decimal(offset), center)
        for i, name in enumerate(('eta_e', 'eta_mu')):
            d = [vals[b][i]-vals[a][i] for a, b in zip(NS, NS[1:])]
            A = (d[1]-d[0]/16)/(Decimal('0.5')-Decimal(1)/16)
            B = d[0]-A
            predicted = A/4+B/256
            jump = actual[i]-extrap[i]
            apred = -span*jump/(64*12)
            spacing = Decimal(10)**(center[i].adjusted()-PRECISION+1)
            components.append(dict(offset=offset, component=name,
                signed_differences_MeV=d, absolute_differences_MeV=list(map(abs, d)),
                recorded_radius_MeV=radii[i], absolute_differences_over_radius=[abs(v)/radii[i] for v in d],
                contraction_ratios=[abs(b/a) if a else None for a, b in zip(d, d[1:])],
                endpoint_Decimal_spacing_MeV=spacing,
                absolute_differences_over_Decimal_spacing=[abs(v)/spacing for v in d],
                fitted_first_order_A_MeV=A, fitted_fourth_order_B_MeV=B,
                predicted_third_difference_MeV=predicted,
                third_prediction_relative_residual=(d[2]-predicted)/d[2] if d[2] else None,
                endpoint_rhs_jump_MeV_s=jump, endpoint_sampling_predicted_A_MeV=apred,
                A_prediction_relative_residual=(A-apred)/A if A else None))
    return dict(point=old['point'], original_status=old['status'], span_s=span,
                segment_durations_s=segments, original_grids=grids,
                endpoint_x_minus_regular_extrapolation=nodes[-1][1]-extrapolated,
                components=components,
                interpretation='Grid-sampling diagnostic, not a discretization enclosure or Decimal roundoff estimate.')


def verify_preserved(repo, result):
    """Authenticate all original D1 files/scripts against START, not a new hash."""
    paths = subprocess.check_output(['git', '-C', str(repo), 'ls-tree', '-r', '--name-only', START,
                str(REFINEMENT/'d1'), 'tools/linux/native_qualification/verify_chemical_reference.py',
                'tools/linux/native_qualification/test_chemical_reference.py',
                'tools/linux/native_qualification/chemical_reference.sbatch']).decode().splitlines()
    import hashlib
    checks = []
    for relative in paths:
        blob = subprocess.check_output(['git', '-C', str(repo), 'show', f'{START}:{relative}'])
        expected = hashlib.sha256(blob).hexdigest()
        actual = digest(repo/relative)
        checks.append(dict(path=relative, expected_sha256=expected, actual_sha256=actual,
                           status='PASS' if actual == expected else 'FAIL'))
    result['preserved_original_D1'] = dict(starting_commit=START, files=checks,
        status='PASS' if all(r['status']=='PASS' for r in checks) else 'FAIL')
    require(paths and result['preserved_original_D1']['status']=='PASS', 'original-evidence', checks)


def stronger_checks(out, model):
    """Additional conservative thermal test; never reduce original u_ref."""
    B = out['thermal_sensitivity']['uniform_thermal_bound_MeV']
    u = tuple(out['u_ref_discretization_MeV'][i]+max(out['u_ref_thermal_MeV'][i], B) for i in range(2))
    out['u_supported_MeV'] = u
    out['stronger_chemical_state'] = [dict(component=r['component'],
        ratio=(r['O2_error_MeV']+u[i])/r['recorded_radius_MeV'],
        actual_error_MeV=r['O2_error_MeV'], reference_uncertainty_MeV=u[i],
        original_budget_MeV=r['recorded_radius_MeV']) for i, r in enumerate(out['chemical_state'])]
    for r in out['stronger_chemical_state']:
        r['status'] = 'PASS' if r['ratio']<=1 else 'FAIL'
    u_E, corners = corner_uncertainty(model.energy, out['eta_ref_MeV'], u)
    erow = out['chemical_energy']
    ratio = (erow['O2_error_MeV']+u_E)/erow['recorded_U_Echem_MeV']
    out['stronger_chemical_energy'] = dict(actual_error_MeV=erow['O2_error_MeV'],
        reference_uncertainty_MeV=u_E, original_budget_MeV=erow['recorded_U_Echem_MeV'],
        four_corner_energies_MeV=corners, ratio=ratio, status='PASS' if ratio<=1 else 'FAIL')
    out['supported_budget_check'] = 'PASS' if ratio<=1 and all(r['status']=='PASS' for r in out['stronger_chemical_state']) else 'FAIL'


# Original C2/C3/C4 input joins and gates retained verbatim from START.
def verify_point(point, position, study, model, candidates, witnesses, steps, contexts, audits):
    directory = study / f'point-{point}/atol-5e-27'
    checkpoint = json.loads((directory/'checkpoint.json').read_text())
    state_rows = {r['channel']: r for r in rows(directory/'state.tsv')}
    diagnostics = {r['channel']: r for r in rows(directory/'diagnostics.tsv')}
    candidate = [r for r in candidates if int(r['index']) == point and float(r['atol']) == 5e-27]
    require(len(candidate) == 1, 'input-structure', f'candidate selection point {point}')
    candidate = candidate[0]
    o2 = tuple(lift(v) for v in checkpoint['O2_state'][1:])
    require(o2 == tuple(lift(candidate[k]) for k in ('e', 'mu')), 'input-structure', f'O2 state disagreement point {point}')
    require(checkpoint['O2_atol'][1:] == [5e-27, 5e-27], 'input-structure', f'O2 atol mismatch point {point}')
    left, obs = lift(checkpoint['main_left_time_s']), lift(checkpoint['time_s'])
    span = obs-left
    require(span > 0 and obs == lift(candidate['t_s']), 'input-structure', f'interval mismatch point {point}')
    initial = tuple(lift(v) for v in checkpoint['main_left_state'][1:])
    out = dict(point=point, time_origin_s=left, observation_time_s=obs, span_s=span,
               initial_eta_MeV=initial, O2_eta_MeV=o2, input_hashes={str(p): digest(p) for p in sorted(directory.iterdir()) if p.is_file()},
               status='RUNNING', C1='NOT_RUN', C2='NOT_RUN', C3='NOT_RUN', C4='NOT_RUN')
    path, out['thermal'] = make_path(point, checkpoint, steps, contexts, position)
    if out['thermal']['status'] != 'PASS':
        out['status'] = 'FAIL'
        out['failure_gate'] = 'thermal-band'
        return out
    out['initial_state_conversion'] = verify_initial(checkpoint)
    out['thermal_sensitivity'] = sensitivity_bound(model, path.nodes, initial)
    require(out['thermal_sensitivity']['status']=='PASS', 'sensitivity-domain', dict(point=point, detail=out['thermal_sensitivity']))
    out['steps_per_segment'] = KS
    out['total_steps'] = {8*k: k*(len(path.nodes)-1) for k in KS}
    out['segment_count'] = len(path.nodes)-1
    integrations = {}
    c1 = {}
    for offset in (-DELTA_X, Decimal(0), DELTA_X):
        values = {8*k: aligned_rk4(model.rhs, path.nodes, initial, k, offset) for k in KS}
        integrations[str(offset)] = values
        c1[str(offset)] = convergence(values)
    out['integrations_eta_MeV'] = integrations
    out['convergence'] = c1
    out['C1'] = 'PASS' if all(r['status']=='PASS' for rs in c1.values() for r in rs) else 'FAIL'
    ref = integrations['0'][512]
    discretization = tuple(2*abs(ref[i]-integrations['0'][256][i]) for i in range(2))
    thermal = tuple(max(abs(v[512][i]-ref[i]) for v in integrations.values()) for i in range(2))
    uncertainty = tuple(discretization[i]+thermal[i] for i in range(2))
    out.update(eta_ref_MeV=ref, u_ref_discretization_MeV=discretization,
               u_ref_thermal_MeV=thermal, u_ref_MeV=uncertainty)
    c2 = []
    for i, channel in enumerate(('eta_e', 'eta_mu')):
        row = state_rows[channel]
        require(lift(row['O2']) == o2[i], 'input-structure', f'state.tsv O2 mismatch point {point} {channel}')
        radius = lift(row['propagated_state_radius'])
        error = abs(o2[i]-ref[i])
        ratio = (error+uncertainty[i])/radius
        c2.append(dict(component=channel, O2_error_MeV=error, O2_error_ulp=error/ulp(o2[i]),
                       recorded_radius_MeV=radius, ratio=ratio, status='PASS' if ratio<=1 else 'FAIL'))
    out['chemical_state'] = c2
    out['C2'] = 'PASS' if all(r['status']=='PASS' for r in c2) else 'FAIL'
    erow = diagnostics['Echem']
    require(lift(erow['O2']) == lift(candidate['Echem']) and lift(erow['U']) == lift(candidate['U']),
            'input-structure', f'energy record disagreement point {point}')
    e_ref = model.energy(ref)
    u_e, corners = corner_uncertainty(model.energy, ref, uncertainty)
    error = abs(lift(erow['O2'])-e_ref)
    allowance = lift(erow['U'])
    ratio = (error+u_e)/allowance
    out['chemical_energy'] = dict(reference_MeV=e_ref, recorded_O2_MeV=lift(erow['O2']),
                                  O2_error_MeV=error, u_E_MeV=u_e, four_corner_energies_MeV=corners,
                                  recorded_U_Echem_MeV=allowance, recorded_F_MeV=lift(erow['F']),
                                  recorded_threshold_MeV=lift(erow['threshold']),
                                  original_budget_utilization=lift(erow['utilization']), ratio=ratio,
                                  status='PASS' if ratio<=1 else 'FAIL')
    out['C3'] = out['chemical_energy']['status']
    local_witnesses = [r for r in witnesses if int(r['index']) == point]
    require(len(local_witnesses) == 6, 'input-structure', f'witness count point {point}')
    targets = [('O2', 42*position+14, o2)]
    for row in local_witnesses:
        method, subdivisions = int(row['method']), int(row['subdivisions'])
        context = 42*position+37+3*method+(4, 8, 16).index(subdivisions)
        targets.append((f'{("rk8pd", "rkf45")[method]}/{subdivisions}', context, tuple(lift(row[k]) for k in ('e', 'mu'))))
    c4 = []
    for name, context, eta in targets:
        last = contexts[context][-1]
        require(tuple(lift(last[k]) for k in ('e', 'mu')) == eta and lift(last['t_s']) == obs,
                'input-structure', dict(point=point, context=context, reason='C4 endpoint/context mismatch'))
        audit = audits[context]
        require(audit['supported_trace'], 'input-structure', f'unsupported C4 audit context {context}')
        # These audit outputs are already decimal mathematical estimates, not doubles.
        enclosure = tuple(Decimal(v) for v in audit['chemical_roundoff_enclosure_MeV'])
        errors = tuple(abs(eta[i]-ref[i]) for i in range(2))
        ratios = tuple(errors[i]/(enclosure[i]+uncertainty[i]) for i in range(2))
        c4.append(dict(name=name, context=context, error_MeV=errors,
                       error_ulp=tuple(errors[i]/ulp(eta[i]) for i in range(2)),
                       enclosure_ctx_MeV=enclosure, ratios=ratios,
                       status='PASS' if max(ratios)<=1 else 'FAIL', report_only=True))
    out['audit_validation'] = c4
    out['C4'] = 'PASS' if all(r['status']=='PASS' for r in c4) else 'FAIL'
    stronger_checks(out, model)
    out['status'] = 'PASS' if all(out[k]=='PASS' for k in ('C1', 'C2', 'C3', 'supported_budget_check')) else 'FAIL'
    return out


def execute(args, result):
    repo = Path(__file__).resolve().parents[3]
    verify_preserved(repo, result)
    head = subprocess.check_output(['git', '-C', str(repo), 'rev-parse', 'HEAD']).decode().strip()
    result['predeclared_execution_commit'] = head
    result['implementation_hashes'] = {}
    # Require these bytes to have been committed BEFORE the job starts.
    for rel in [RESOLUTION/'PREDECLARATION.md',
                Path('tools/linux/native_qualification/resolve_chemical_reference.py'),
                Path('tools/linux/native_qualification/test_chemical_reference_resolution.py'),
                Path('tools/linux/native_qualification/chemical_reference_resolution.sbatch')]:
        import hashlib
        committed = subprocess.check_output(['git', '-C', str(repo), 'show', f'{head}:{rel}'])
        require(hashlib.sha256(committed).hexdigest() == digest(repo/rel), 'predeclaration', str(rel))
        result['implementation_hashes'][str(rel)] = digest(repo/rel)
    study, inputs = authenticate(repo, repo/REFINEMENT/'d1/review/BNV-THERMO-01_FABLE_ADJUDICATION.md', REVIEW_SHA, result)
    model = Chemistry(rows(inputs['coefficients.tsv'])[0], rows(inputs['tangent.tsv'])[0])
    result['constants'] = dict(kb_MeV_K=model.kb, erg_per_MeV=model.erg, kb_erg_K=model.kberg,
        Bdot_count_s=model.bdot, pi=PI, Z=model.z, Ltilde=model.luminosity,
        sigma_count_s=model.sigma, source_MeV_s=model.source)
    old_records = {p: json.loads((repo/REFINEMENT/f'd1/evidence/point-{p}.json').read_text()) for p in POINTS}
    result['original_failure_analysis'] = [old_diagnosis(old_records[p], model) for p in (1, 64, 1440)]
    save(args.output/'original-failure-analysis.json', result['original_failure_analysis'])
    trace = rows(study/'local-rhs.tsv')
    result['rhs_fidelity'] = fidelity(model, trace)
    save(args.output/'results.json', result)
    require(result['rhs_fidelity']['status']=='PASS', 'rhs-fidelity', result['rhs_fidelity']['failures'])
    contexts = collections.defaultdict(list)
    for row in trace:
        contexts[int(row['context'])].append(row)
    audits = {r['context']: r for r in json.loads((repo/REFINEMENT/'evidence/conditional-roundoff.json').read_text())['records']}
    candidates, witnesses, steps = [rows(study/name) for name in ('candidates.tsv', 'witnesses.tsv', 'witness-steps.tsv')]
    require(tuple(sorted({int(r['index']) for r in candidates})) == POINTS, 'input-structure', 'point list changed')
    for position, point in enumerate(POINTS):
        record = verify_point(point, position, study, model, candidates, witnesses, steps, contexts, audits)
        record['input_hashes'].update({r['path']: r['actual_sha256'] for r in result['authentication']['files']
            if Path(r['path']).name in ('candidates.tsv', 'witnesses.tsv', 'witness-steps.tsv',
                'local-rhs.tsv', 'coefficients.tsv', 'tangent.tsv', 'conditional-roundoff.json', AUDIT.name)})
        old = old_records[point]
        require(record['initial_eta_MeV'] == tuple(map(Decimal, old['initial_eta_MeV'])),
                'initial-state', f'original/revised left state mismatch point {point}')
        if 'eta_ref_MeV' in record:
            record['original_to_aligned_reference_shift_MeV'] = [record['eta_ref_MeV'][i]-Decimal(old['eta_ref_MeV'][i]) for i in range(2)]
            require(record['thermal']['nodes_tau_s_x_dx'] == [tuple(map(Decimal, n)) for n in old['thermal']['nodes_tau_s_x_dx']],
                    'thermal-band', f'original/revised path mismatch point {point}')
        result['points'].append(record)
        if record['status'] == 'FAIL':
            for gate in ('C1', 'C2', 'C3', 'supported_budget_check'):
                if record.get(gate) == 'FAIL':
                    key = {'C1': 'convergence', 'C2': 'chemical_state', 'C3': 'chemical_energy',
                           'supported_budget_check': 'stronger_chemical_state'}[gate]
                    result['failure_details'].append(dict(gate=gate, point=point, detail=record[key],
                        stronger_energy=record.get('stronger_chemical_energy')))
            if record.get('failure_gate'):
                result['failure_details'].append(dict(gate=record['failure_gate'], point=point, detail=record['thermal']))
        save(args.output/f'point-{point}.json', record)
        save(args.output/'results.json', result)
        print(f'point={point} status={record["status"]} C1={record["C1"]} C2={record["C2"]} C3={record["C3"]} C4={record["C4"]} supported={record.get("supported_budget_check")}', flush=True)
    result['status'] = 'PASS_CONDITIONAL' if all(r['status']=='PASS' for r in result['points']) else 'FAIL'
    failures = [r['path'] for r in result['authentication']['files'] if digest(r['path']) != r['actual_sha256']]
    require(not failures, 'authentication', dict(reason='input changed during execution', paths=failures))
    verify_preserved(repo, result)
    result['post_execution_input_hashes'] = 'PASS'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    result = dict(schema='bnv-thermo-01-d1-resolution-v1', status='RUNNING', starting_commit=START,
        original_D1_disposition='FAIL_PRESERVED', precision_digits=PRECISION,
        K_per_segment=KS, nominal_N=NS, delta_x=DELTA_X, rho_limit=RHO_LIMIT,
        points=[], failure_details=[],
        units=dict(time='s', eta='MeV', Echem='MeV', x='ln(Tinf/1e8 K)'),
        assumptions=['Exact retained binary64 left state; no accumulated main-prefix error.',
          'Sampled thermal agreement does not prove a uniform continuous enclosure.',
          'Constant-offset runs are sensitivity probes, not extrema for arbitrary time-dependent perturbations.',
          'Uniform sensitivity bound assumes |x_true-xhat|<=1e-12 and the checked first-exit domain.',
          'Twice the last refinement difference is a convergence-supported allowance, conditional on continued contraction beyond tested levels; not a rigorous IVP enclosure.',
          'C4 is report-only: ratio>1 invalidates the audit enclosure as a bound, not independently D1.'],
        exclusions=['all Stage A checkpoints', 'main-prefix error from t=0', 'thermal RHS roundoff',
          'R18/R20 and complete paired-run restart', 'chemical heating in one-year temperature', 'Stage B readiness'],
        production_change_applied=False, stage_A_launched=False)
    try:
        require(os.environ.get('SLURM_JOB_ID') and len(os.sched_getaffinity(0)) == 1,
                'execution', 'requires Slurm allocation and one logical execution CPU')
        with localcontext() as ctx:
            ctx.prec = PRECISION
            execute(args, result)
    except GateFailure as error:
        result['status'] = 'FAIL'
        result['failure_details'].append(dict(gate=error.gate, detail=error.detail))
    except Exception as error:
        result['status'] = 'FAIL'
        result['failure_details'].append(dict(gate='implementation', detail=str(error), traceback=traceback.format_exc()))
    finally:
        completed = {r['point'] for r in result['points']}
        for point in POINTS:
            if point not in completed:
                r = dict(point=point, status='NOT_RUN', reason='Earlier mandatory gate failed.')
                result['points'].append(r)
                save(args.output/f'point-{point}.json', r)
        result['execution'] = dict(job_id=os.environ.get('SLURM_JOB_ID'), hostname=platform.node(),
            python=sys.version, affinity=sorted(os.sched_getaffinity(0)), elapsed_s=time.monotonic()-start,
            max_rss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            script_sha256=digest(__file__))
        save(args.output/'results.json', result)
        save(args.output/'failure-details.json', result['failure_details'])
        print(result['status'], 'failure gates:', [r['gate'] for r in result['failure_details']], flush=True)
    return 0 if result['status']=='PASS_CONDITIONAL' else 1


if __name__ == '__main__':
    sys.exit(main())

