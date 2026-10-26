#!/usr/bin/env python3
"""Bounded, write-once execution of the predeclared evidence experiment."""
import decimal
import hashlib
import json
import math
import os
import pathlib
import resource
import subprocess
import sys
import time

S = pathlib.Path(__file__).resolve().parents[3]
R = S.parents[1]
O = R / 'ba13-evidence-20261026'
OLD = R / 'ba13-closure-20261019'
PRE = '93c28c3bd85afea8a893a34a8a6dd92ce06c52be'
D = decimal.Decimal
decimal.getcontext().prec = 60


def read(path):
    return json.loads(path.read_text(), parse_float=D)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cpu():
    return sum(x.ru_utime + x.ru_stime for x in (
        resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)))


def write(path, obj):
    with path.open('x') as f:
        json.dump(obj, f, indent=2, allow_nan=False,
                  default=lambda v: str(v) if isinstance(v, D) else None)
        f.write('\n')


def tangent_report():
    records = {}
    old_m1 = read(OLD / 'assessment-final/m1.json')
    for i in (0, 10, 20):
        raw = {n: {s: read((OLD if n == 512 else O) /
                (f't{i}-h0.000025-{s}.json' if n == 512 else f't{i}-n{n}-offset-{s}.json'))
                   for s in (-1, 1)} for n in (512, 1024, 2048)}
        for k, name in ((0, 't_n'), (2, 't_e'), (3, 't_mu')):
            ts = {n: (v[1]['N'][k] - v[-1]['N'][k]) /
                  (v[1]['primary'][0] - v[-1]['primary'][0]) for n, v in raw.items()}
            d1, d2 = ts[1024] - ts[512], ts[2048] - ts[1024]
            p = (abs(d1 / d2)).ln() / D(2).ln() if d1 and d2 and d1*d2 > 0 else None
            order_ok = p is not None and D('3.5') <= p <= D('4.5')
            prior = old_m1['reference'][str(i)][name]
            width_ok = prior['stencil_status'] == 'ORDER_COMPATIBLE'
            # Exact difference of correlated ratios, including the B denominator.
            bn = {n: raw[n][1]['primary'][0]-raw[n][-1]['primary'][0] for n in raw}
            nn = {n: raw[n][1]['N'][k]-raw[n][-1]['N'][k] for n in raw}
            correlated = ((nn[2048]-nn[1024])-ts[1024]*(bn[2048]-bn[1024]))/bn[2048]
            # Common-H sample movement is measured, not charged as a fictitious
            # independent endpoint error. Root uncertainty remains a separate term.
            bp = {n: (raw[n][1]['primary'][0]+raw[n][-1]['primary'][0])/2 for n in raw}
            records[f'{i}:{name}'] = dict(
                target=i, quantity=name, t={str(n): v for n, v in ts.items()},
                increments=[d1, d2], order=p, order_ok=order_ok,
                stencil_ratio=prior['stencil_ratio'], width_ok=width_ok,
                stencil_tail_observation=prior['stencil'],
                correlated_ratio_increment=correlated,
                correlated_identity_residual=correlated-d2,
                B_span={str(n): v for n, v in bn.items()},
                B_midpoint={str(n): v for n, v in bp.items()},
                fine_tail=abs(d2)/15 if order_ok else None,
                coarse_tail=abs(ts[2048]-ts[512])+abs(d2)/15 if order_ok else None,
                evidence_type='validated_numerical_estimate' if order_ok else 'unresolved',
                standard_id='offset-H-GL4-three-resolution-v1', authority=None,
                historical_uncorrelated_estimate=prior['integration'],
                status='PASS_EXPERIMENT' if order_ok and width_ok else 'INCONCLUSIVE')
    failed = [k for k, v in records.items() if v['status'] != 'PASS_EXPERIMENT']
    return dict(records=records, failed=failed,
                numerical_expansion_allowed=not failed,
                status='STOP' if failed else 'PASS_EXPERIMENT',
                qualification='UNQUALIFIED',
                reason='failed frozen tangent order/width gate' if failed else
                'experiment passes; reference remainders and authority remain unqualified')


def main():
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Slurm required')
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    O.mkdir(exist_ok=True)
    phase = sys.argv[1]
    previous = sum(float(json.loads(p.read_text())['numerical_cpu_seconds'])
                   for p in O.glob('receipt-*.json'))
    begin, wall = cpu(), time.monotonic()
    records, status = [], 'ERROR'
    decl = S / 'docs/validation/ba13-evidence/PREDECLARATION.md'
    committed = subprocess.check_output(['git', '-C', str(S), 'show',
                                        PRE+':docs/validation/ba13-evidence/PREDECLARATION.md'])
    if hashlib.sha256(committed).hexdigest() != sha(decl):
        raise RuntimeError('predeclaration changed')

    def call(label, command, env=None):
        remaining = 180 - previous - (cpu()-begin)
        if remaining < 15:
            raise RuntimeError('numerical reserve stop')
        cap = min(30, math.floor(remaining-5))
        def limit():
            resource.setrlimit(resource.RLIMIT_CPU, (cap, cap))
        a, w = cpu(), time.monotonic()
        stdout, stderr = O/(label+'.json'), O/(label+'.err')
        with stdout.open('x') as f, stderr.open('x') as e:
            result = subprocess.run(list(map(str, command)), stdout=f, stderr=e,
                                    env=env, preexec_fn=limit, timeout=60)
        records.append(dict(label=label, command=list(map(str, command)),
                            cpu_seconds=cpu()-a, wall_seconds=time.monotonic()-w,
                            returncode=result.returncode, stdout_sha256=sha(stdout),
                            stderr_sha256=sha(stderr)))
        if result.returncode:
            raise RuntimeError(label+' failed')
        return read(stdout)

    try:
        exe = OLD/'assessment'
        retained = json.loads((OLD/'receipt-m1.json').read_text())
        if sha(exe) != retained['executable_sha256']:
            raise RuntimeError('retained executable hash mismatch')
        if phase == 'tangent':
            targets = read(S/'docs/validation/x1-oracle-anchor/targets.json')
            for i in (0, 10, 20):
                hc = read(R/f'x2-fixed-step/t{i}-n4096.json')['values'][8]
                for n in (1024, 2048):
                    for sign in (-1, 1):
                        h = hc*(1+D('0.000025')*sign)
                        call(f't{i}-n{n}-offset-{sign}',
                             [exe, h, targets[str(i)]['surfaceP_cgs'], n, '1e-4', 96, 0])
            result = tangent_report()
            write(O/'tangents.json', result)
            status = result['status']
        elif phase == 'freeze':
            if not read(O/'tangents.json')['numerical_expansion_allowed']:
                raise RuntimeError('tangent stop')
            call('frozen-holdouts', [sys.executable, '-B', S/'tools/linux/ba13_evidence/freeze.py'])
            status = 'PASS'
        elif phase == 'holdouts':
            if not read(O/'tangents.json')['numerical_expansion_allowed']:
                raise RuntimeError('tangent stop')
            frozen = read(O/'frozen-holdouts.json')
            # Numeric criteria must be committed, not merely written, first.
            frozen_doc = S/'docs/validation/ba13-evidence/evidence/frozen-holdouts.json'
            if sha(frozen_doc) != sha(O/'frozen-holdouts.json'):
                raise RuntimeError('frozen criterion mismatch')
            blob = subprocess.check_output(['git','-C',str(S),'show',
                'HEAD:docs/validation/ba13-evidence/evidence/frozen-holdouts.json'])
            if hashlib.sha256(blob).hexdigest() != sha(frozen_doc):
                raise RuntimeError('criterion not committed')
            from reference import load_reference, normalized_fields
            comp = load_reference()[0]
            root_exe = R/'x2-fixed-step/build/extended'
            manifest = json.loads((S/'docs/validation/x2-fixed-step/evidence/external-manifest.json').read_text())
            if sha(root_exe) != manifest['files']['build/extended']['sha256']:
                raise RuntimeError('root executable hash mismatch')
            scores = {}
            status = 'PASS_EXPERIMENT'
            for i in (5,15):
                t = frozen['targets'][str(i)]
                root = call(f'holdout-{i}-root', [root_exe,t['B'],t['H_guess'],t['surfaceP_cgs'],
                            2048,'1e-4',96,0,0,0])
                result = call(f'holdout-{i}-functional', [exe,root['values'][8],t['surfaceP_cgs'],
                              2048,'1e-4',96,0])
                if abs(result['primary'][0]/D(t['B'])-1)>D('2e-16'):
                    raise RuntimeError('holdout alignment violates frozen root residual')
                values = normalized_fields(result,comp)
                score = {}
                for n,gate in t['groups'].items():
                    diff=max(abs(float(v)-float(p)) for v,p in zip(values[n],gate['prediction']))
                    tol=float(gate['criterion'])
                    score[n]=dict(discrepancy=diff,criterion=tol,ratio=diff/tol,
                                  status='PASS' if diff<=tol else 'FAIL')
                scores[str(i)]=dict(groups=score,pressure=t['surfaceP_cgs'],B_target=t['B'],
                    B_achieved=result['primary'][0],root_shots=root['shots'])
                if any(v['status']=='FAIL' for v in score.values()):
                    status='STOP_HOLDOUT_FAIL'
                    break
            write(O/'holdouts.json',dict(status=status,targets=scores,
                  criterion_sha256=sha(frozen_doc),continuous_theorem=False,
                  qualification='UNQUALIFIED',reference_transfer_allowed=status=='PASS_EXPERIMENT'))
        elif phase == 'historical-audit':
            mirror = O/'historical-audit'
            mirror.mkdir()
            for path in OLD.glob('t*.json'):
                (mirror/path.name).symlink_to(path)
            (mirror/'assessment-final').symlink_to(OLD/'assessment-final', target_is_directory=True)
            env = dict(os.environ, BA13_CLOSURE_OUTPUT_ROOT=str(mirror),
                       BA13_CLOSURE_AUDIT_NAME='audit-completion.json')
            call('historical-audit-stdout',
                 [sys.executable, '-B', R/'worktrees/ba13-closure/tools/linux/ba13_closure/audit.py'], env)
            status = 'PASS'
        else:
            raise ValueError(phase)
    finally:
        write(O/f'receipt-{phase}.json', dict(
            phase=phase, status=status, job=os.environ['SLURM_JOB_ID'],
            affinity=sorted(os.sched_getaffinity(0)), prior_cpu_seconds=previous,
            numerical_cpu_seconds=cpu()-begin, wall_seconds=time.monotonic()-wall,
            peak_rss_kib=max(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                             resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss),
            predeclaration_commit=PRE, predeclaration_sha256=sha(decl),
            source_commit=subprocess.check_output(['git','-C',str(S),'rev-parse','HEAD'],text=True).strip(),
            runner_sha256=sha(pathlib.Path(__file__)), executable_sha256=sha(OLD/'assessment'),
            solver_sha256=sha(S/'tests/bnv/x2_fixed_step.cpp'),
            passive_source_sha256=sha(S/'tests/bnv/x2_uncertainty.cpp'), runs=records))
    print(json.dumps(dict(phase=phase, status=status, numerical_cpu_seconds=cpu()-begin)))


if __name__ == '__main__':
    main()
