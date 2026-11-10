#!/usr/bin/env python3
"""BNV-THERMO-01-D1: bounded, offline Decimal reference (scientific candidate).

No production solver is imported or executed. Mathematical definitions are
source-bound to audit_chemical_roundoff.py at REVIEWED, not its recursion.
Retained floating records are decoded as binary64 BEFORE exact Decimal lifting.
The adjudication and its independently supplied digest are mandatory inputs.
"""
import argparse
import bisect
import collections
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import resource
import subprocess
import sys
import time
import traceback
from decimal import Decimal, localcontext

REVIEWED = 'ad4e91aee3837b8b07e7ed723afdd40b765011fe'
POINTS = (1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 960, 1440, 1919)
NS = (64, 128, 256, 512)
PRECISION = 110
DELTA_X = Decimal('1e-12')
RHO_LIMIT = Decimal('1e-13')
PI = Decimal('3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068')
BASE = Path('docs/validation/results-first-bnv-pilot')
REFINEMENT = BASE / 'echem-refinement'
AUDIT = Path('tools/linux/native_qualification/audit_chemical_roundoff.py')


def lift(value):
    """Exact value of the retained double, NOT Decimal(printed decimal)."""
    value = float(value)
    if not math.isfinite(value):
        raise ValueError('non-finite retained binary64 input')
    return Decimal.from_float(value)


def ulp(value):
    return lift(math.ulp(float(value)))


def rows(path):
    with Path(path).open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def save(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, default=str, allow_nan=False) + '\n')
    temporary.replace(path)


class GateFailure(Exception):
    def __init__(self, gate, detail):
        super().__init__(f'{gate}: {detail}')
        self.gate, self.detail = gate, detail


def require(condition, gate, detail):
    if not condition:
        raise GateFailure(gate, detail)


def authenticate(repo, adjudication, expected_review_hash, result):
    """Recheck every external byte; no prior review's check is assumed."""
    auth = result['authentication'] = {'status': 'RUNNING', 'files': [], 'failures': []}

    def git(*args):
        return subprocess.check_output(['git', '-C', str(repo), *args])

    def committed(relative):
        content = git('show', f'{REVIEWED}:{relative}')
        path = repo / relative
        # Reports may acquire addenda; consumed definitions/inventories may not.
        check(path, hashlib.sha256(content).hexdigest(), f'git:{REVIEWED}:{relative}')
        return content

    def check(path, expected, authority):
        path = Path(path)
        actual = digest(path) if path.is_file() else None
        item = dict(path=str(path), expected_sha256=expected, actual_sha256=actual,
                    authority=str(authority), status='PASS' if actual == expected else 'FAIL')
        auth['files'].append(item)
        if item['status'] == 'FAIL':
            auth['failures'].append(item)
        return item

    require(git('rev-parse', REVIEWED + '^{commit}').decode().strip() == REVIEWED,
            'authentication', 'reviewed commit cannot be resolved exactly')
    auth['reviewed_commit'] = REVIEWED
    auth['execution_head'] = git('rev-parse', 'HEAD').decode().strip()
    require(subprocess.run(['git', '-C', str(repo), 'merge-base', '--is-ancestor',
                            REVIEWED, 'HEAD']).returncode == 0,
            'authentication', 'HEAD is not a descendant of reviewed commit')
    auth['reviewed_commit_object_sha256'] = hashlib.sha256(git('cat-file', 'commit', REVIEWED)).hexdigest()
    manifest_path = REFINEMENT / 'evidence/external-SHA256SUMS'
    manifest = committed(manifest_path).decode()
    inventory = {}
    for line in manifest.splitlines():
        sha, name = line.split(maxsplit=1)
        require(name not in inventory, 'authentication', f'duplicate manifest path {name}')
        inventory[name] = sha
        check(name, sha, manifest_path)
    require(len(inventory) == 354, 'authentication', f'external file count {len(inventory)} != 354')
    auth['external_files_checked'] = len(inventory)
    # Discover the study from authenticated committed paths, never a guessed run.
    candidates = [Path(p) for p in inventory if p.endswith('/study/output/candidates.tsv')
                  and '/pre-reconstruction-attempt-' not in p]
    require(len(candidates) == 1, 'authentication', 'ambiguous study candidates path')
    study = candidates[0].parent
    required = [study / f for f in ('candidates.tsv', 'witnesses.tsv', 'witness-steps.tsv', 'local-rhs.tsv')]
    required += [study / f'point-{i}/atol-5e-27/{name}' for i in POINTS
                 for name in ('checkpoint.json', 'state.tsv', 'diagnostics.tsv')]
    for path in required:
        require(str(path) in inventory, 'authentication', f'unmanifested required input: {path}')
    committed(AUDIT)
    committed(REFINEMENT / 'evidence/conditional-roundoff.json')
    # These two physical tables are outside the 354-file study inventory.
    # Their reviewed-commit authentication record supplies the prior hashes.
    physical_record = BASE / 'native-regeneration/evidence/authentication.json'
    physical = json.loads(committed(physical_record))
    inputs = {}
    for name in ('coefficients.tsv', 'tangent.tsv'):
        entries = [r for r in physical['checks'] if Path(r['path']).name == name]
        require(len(entries) == 1, 'authentication', f'ambiguous physical input {name}')
        entry = entries[0]
        check(entry['path'], entry['sha256'], physical_record)
        inputs[name] = Path(entry['path'])
    if not adjudication or not expected_review_hash:
        auth['failures'].append(dict(input='BNV-THERMO-01_FABLE_ADJUDICATION.md',
                                     reason='Missing adjudication path and/or independently supplied SHA-256; no substitute review invented.'))
    else:
        require(len(expected_review_hash) == 64 and all(c in '0123456789abcdef' for c in expected_review_hash),
                'authentication', 'invalid adjudication SHA-256')
        check(adjudication, expected_review_hash, 'owner-supplied review identity')
    auth['status'] = 'FAIL' if auth['failures'] else 'PASS'
    result['study_directory'] = str(study)
    require(not auth['failures'], 'authentication', auth['failures'])
    return study, inputs


class Chemistry:
    """Real-valued RHS; x=ln(Tinf/1e8 K), eta in MeV, time in seconds."""
    def __init__(self, coefficients, tangent):
        self.z = [[lift(coefficients[f'Z{i}{j}']) for j in range(2)] for i in range(2)]
        self.luminosity = [lift(coefficients[k]) for k in ('Ltilde_Me', 'Ltilde_Mmu')]
        self.bdot = lift(-2.413652026364133e37)  # count/s, source-bound double
        self.sigma = [-lift(tangent[k]) * self.bdot for k in ('t_e', 't_mu')]
        self.source = [-sum(self.z[i][j] * self.sigma[j] for j in range(2)) for i in range(2)]
        # Preserve the binary64 constant expressions in the authenticated audit.
        self.kb = lift(8.61733326214518e-5 * 1e-6)  # MeV/K
        self.erg = lift(1.602176634e33 / 1e39)       # erg/MeV
        self.kberg = lift(float(self.kb) * float(self.erg))  # erg/K

    def reaction(self, temperature, eta):
        xi = eta / (self.kb * temperature)
        u2 = (xi / PI) ** 2
        return xi / PI**2 * (14680 + u2 * (7560 + u2 * (840 + 24 * u2))) / 11513

    def rhs(self, x, eta):
        temperature = Decimal('1e8') * x.exp()
        rates = [self.luminosity[j] / self.kberg * temperature**7 * self.reaction(temperature, eta[j])
                 for j in range(2)]
        return tuple(self.source[i] - sum(self.z[i][j] * rates[j] for j in range(2)) for i in range(2))

    def energy(self, eta):
        return energy(self.z, eta)


def energy(z, eta):
    e, mu = eta
    det = z[0][0] * z[1][1] - z[0][1] * z[1][0]
    return (z[1][1] * e*e - (z[0][1] + z[1][0]) * e*mu + z[0][0] * mu*mu) / (2 * det)


def corner_uncertainty(evaluate, eta, radius):
    center = evaluate(eta)
    corners = [evaluate((eta[0] + se*radius[0], eta[1] + sm*radius[1]))
               for se in (-1, 1) for sm in (-1, 1)]
    return max(abs(e - center) for e in corners), corners


def rk4(rhs, left, right, initial, n):
    """Classical RK4, uniform rational grid; no adaptive or production stepping."""
    h = (right - left) / n
    y = tuple(initial)
    for k in range(n):
        t = left + h*k
        k1 = rhs(t, y)
        k2 = rhs(t + h/2, tuple(v + h*a/2 for v, a in zip(y, k1)))
        k3 = rhs(t + h/2, tuple(v + h*a/2 for v, a in zip(y, k2)))
        k4 = rhs(t + h, tuple(v + h*a for v, a in zip(y, k3)))
        y = tuple(v + h*(a + 2*b + 2*c + d)/6 for v, a, b, c, d in zip(y, k1, k2, k3, k4))
    return y


class Hermite:
    def __init__(self, nodes):
        self.nodes = nodes  # (tau_s, x, dx/dt), tau measured from exact t_left
        self.times = [n[0] for n in nodes]
        require(len(nodes) >= 2 and all(b > a for a, b in zip(self.times, self.times[1:])),
                'thermal-band', 'Hermite nodes are not strictly ordered')

    def __call__(self, t):
        require(self.times[0] <= t <= self.times[-1], 'thermal-band', f'out-of-domain thermal query {t}')
        i = min(bisect.bisect_right(self.times, t) - 1, len(self.nodes)-2)
        a, x, dx = self.nodes[i]
        b, y, dy = self.nodes[i+1]
        h = b-a
        s = (t-a)/h
        return (2*s**3-3*s**2+1)*x + (s**3-2*s**2+s)*h*dx + (-2*s**3+3*s**2)*y + (s**3-s**2)*h*dy


def convergence(values):
    records = []
    for i, component in enumerate(('eta_e', 'eta_mu')):
        differences = [abs(values[b][i]-values[a][i]) for a, b in zip(NS, NS[1:])]
        floor = Decimal('1e-60') * abs(values[512][i])
        checks = []
        for old, new in zip(differences, differences[1:]):
            ratio = new/old if old else None
            zero_case = 'both_zero' if old == new == 0 else ('zero_denominator_nonzero_numerator' if not old else None)
            passed = new == 0 or new < floor or (old != 0 and ratio <= Decimal('0.5'))
            checks.append(dict(ratio=ratio, zero_case=zero_case, status='PASS' if passed else 'FAIL'))
        records.append(dict(component=component, differences_MeV=differences, floor_MeV=floor,
                            checks=checks, status='PASS' if all(c['status']=='PASS' for c in checks) else 'FAIL'))
    return records


def fidelity(model, trace):
    failures, maxima = [], [Decimal(0), Decimal(0)]
    for line, row in enumerate(trace, 2):
        reference = model.rhs(lift(row['x']), tuple(lift(row[k]) for k in ('e', 'mu')))
        for i, key in enumerate(('de', 'dmu')):
            recorded = lift(row[key])
            difference = abs(recorded-reference[i])
            spacing = ulp(recorded)
            maxima[i] = max(maxima[i], difference/spacing)
            if difference > 2*spacing:
                failures.append(dict(line=line, context=int(row['context']), component=key,
                                     t_s=lift(row['t_s']), row=row, recorded_MeV_s=recorded,
                                     reference_MeV_s=reference[i], difference_MeV_s=difference,
                                     limit_MeV_s=2*spacing, difference_ulp=difference/spacing))
    return dict(status='FAIL' if failures else 'PASS', rows=len(trace), components_checked=2*len(trace),
                maximum_difference_ulp=maxima, failures=failures)


def make_path(point, checkpoint, steps, contexts, position):
    left = lift(checkpoint['main_left_time_s'])
    obs = lift(checkpoint['time_s'])
    context = 42*position + 38  # six tolerances * six contexts, then rk8pd /4,/8,/16
    trace = contexts[context]

    def derivative(t, state):
        matches = {lift(r['dx']) for r in trace if lift(r['t_s']) == t
                   and tuple(lift(r[k]) for k in ('x', 'e', 'mu')) == tuple(state)}
        require(len(matches) == 1, 'thermal-band', dict(point=point, context=context, time_s=t,
                                                       reason='missing or ambiguous matching node derivative'))
        return matches.pop()

    state = tuple(lift(v) for v in checkpoint['main_left_state'])
    nodes = [(Decimal(0), state[0], derivative(left, state))]
    selected = [r for r in steps if int(r['index']) == point and r['method']=='0' and r['subdivisions']=='8']
    previous = left
    for row in selected:
        require(lift(row['left']) == previous, 'thermal-band', f'noncontiguous witness steps point {point}')
        right = lift(row['right'])
        state = tuple(lift(row[k]) for k in ('x', 'e', 'mu'))
        nodes.append((right-left, state[0], derivative(right, state)))
        previous = right
    require(previous == obs, 'thermal-band', f'witness endpoint mismatch point {point}')
    path = Hermite(nodes)
    comparisons = []
    for row in steps:
        if int(row['index']) == point and row['subdivisions']=='16':
            t = lift(row['right'])
            residual = abs(path(t-left)-lift(row['x']))
            comparisons.append(dict(method='rk8pd' if row['method']=='0' else 'rkf45', time_s=t, residual=residual))
    comparisons.append(dict(method='O2', time_s=obs, residual=abs(path(obs-left)-lift(checkpoint['O2_state'][0]))))
    rho = max(r['residual'] for r in comparisons)
    return path, dict(rho=rho, limit=RHO_LIMIT, delta_x=DELTA_X,
                      status='PASS' if rho <= RHO_LIMIT else 'FAIL',
                      nodes_tau_s_x_dx=nodes, comparisons=comparisons)


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
    integrations = {}
    c1 = {}
    for offset in (-DELTA_X, Decimal(0), DELTA_X):
        # Cache only the prescribed thermal interpolant, never the chemical RHS.
        thermal_cache = {}

        def rhs(t, y):
            if t not in thermal_cache:
                thermal_cache[t] = path(t) + offset
            return model.rhs(thermal_cache[t], y)

        values = {n: rk4(rhs, Decimal(0), span, initial, n) for n in NS}
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
    out['status'] = 'PASS' if all(out[k]=='PASS' for k in ('C1', 'C2', 'C3')) else 'FAIL'
    return out


def execute(args, result):
    repo = Path(__file__).resolve().parents[3]
    study, inputs = authenticate(repo, args.adjudication, args.adjudication_sha256, result)
    model = Chemistry(rows(inputs['coefficients.tsv'])[0], rows(inputs['tangent.tsv'])[0])
    result['constants'] = dict(kb_MeV_K=model.kb, erg_per_MeV=model.erg, kb_erg_K=model.kberg,
                               Bdot_count_s=model.bdot, pi=PI, Z=model.z, Ltilde=model.luminosity,
                               sigma_count_s=model.sigma, source_MeV_s=model.source)
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
                                                                  'local-rhs.tsv', 'coefficients.tsv', 'tangent.tsv',
                                                                  'conditional-roundoff.json', AUDIT.name)})
        result['points'].append(record)
        if record['status'] == 'FAIL':
            for gate in ('C1', 'C2', 'C3'):
                if record[gate] == 'FAIL':
                    detail_key = {'C1': 'convergence', 'C2': 'chemical_state', 'C3': 'chemical_energy'}[gate]
                    result['failure_details'].append(dict(gate=gate, point=point, detail=record[detail_key]))
            if record.get('failure_gate'):
                result['failure_details'].append(dict(gate=record['failure_gate'], point=point, detail=record['thermal']))
        save(args.output/f'point-{point}.json', record)
        save(args.output/'results.json', result)
        print(f'point={point} status={record["status"]} C1={record["C1"]} C2={record["C2"]} C3={record["C3"]} C4={record["C4"]}', flush=True)
    result['status'] = 'PASS' if all(r['status']=='PASS' for r in result['points']) else 'FAIL'
    # Reauthenticate after reading/integration: original inputs must remain intact.
    failures = [r['path'] for r in result['authentication']['files'] if digest(r['path']) != r['actual_sha256']]
    require(not failures, 'authentication', dict(reason='input changed during execution', paths=failures))
    result['post_execution_input_hashes'] = 'PASS'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--adjudication', type=Path)
    parser.add_argument('--adjudication-sha256')
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    result = dict(schema='bnv-thermo-01-d1-v1', status='RUNNING', reviewed_commit=REVIEWED,
                  precision_digits=PRECISION, N=NS, delta_x=DELTA_X, rho_limit=RHO_LIMIT,
                  points=[], failure_details=[], units=dict(time='s', eta='MeV', Echem='MeV', x='ln(Tinf/1e8 K)'),
                  assumptions=['Exact retained binary64 left chemical state; no main-prefix error included.',
                               'Sampled thermal agreement is not a proof of uniform continuous enclosure.',
                               'Constant offsets +/-1e-12 are conditional sensitivity probes, not proved extrema over arbitrary time-dependent perturbations.',
                               'C4 is report-only; ratio >1 invalidates the audit enclosure as a bound but is not an independent D1 gate.'],
                  exclusions=['every Stage A checkpoint', 'accumulated error from t=0', 'thermal RHS roundoff qualification',
                              'R18/R20 and complete paired-run restart qualification', 'chemical heating in one-year temperature', 'Stage B readiness'],
                  production_change_applied=False, stage_A_launched=False)
    try:
        require(os.environ.get('SLURM_JOB_ID') and len(os.sched_getaffinity(0)) == 1,
                'execution', 'requires Slurm allocation and exactly one logical execution CPU')
        with localcontext() as context:
            context.prec = PRECISION
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
                record = dict(point=point, status='NOT_RUN', C1='NOT_RUN', C2='NOT_RUN', C3='NOT_RUN', C4='NOT_RUN',
                              reason='An earlier mandatory gate failed; no numerical result is claimed.',
                              eta_ref_MeV=None, contraction_ratios=None, u_ref_discretization_MeV=None,
                              u_ref_thermal_MeV=None, rho=None, O2_error_MeV=None, O2_error_ulp=None,
                              recorded_state_radii_MeV=None, recorded_U_Echem_MeV=None,
                              input_hashes={r['path']: r['actual_sha256'] for r in result.get('authentication', {}).get('files', [])
                                            if f'/point-{point}/atol-5e-27/' in r['path']})
                result['points'].append(record)
                save(args.output/f'point-{point}.json', record)
        result['execution'] = dict(job_id=os.environ.get('SLURM_JOB_ID'), hostname=platform.node(),
                                   python=sys.version, affinity=sorted(os.sched_getaffinity(0)),
                                   elapsed_s=time.monotonic()-start, max_rss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                                   script_sha256=digest(__file__))
        save(args.output/'results.json', result)
        if result['failure_details']:
            save(args.output/'failure-details.json', result['failure_details'])
        print(result['status'], 'failure gates:', [r['gate'] for r in result['failure_details']], flush=True)
    return 0 if result['status']=='PASS' else 1


if __name__ == '__main__':
    sys.exit(main())
