#!/usr/bin/env python3
"""Write-once fixed-Hc M1 evaluations; no shooting and no production solve."""
import decimal, hashlib, json, os, pathlib, resource, subprocess, sys, time

S = pathlib.Path(__file__).resolve().parents[3]
R = S.parents[1]
O = R / 'ba13-closure-20261019'
TC = R / 'toolchain'
PREDECLARATION = 'fdffaaf'

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def cpu():
    return sum(v.ru_utime + v.ru_stime for v in
               (resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)))

def main():
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Slurm required')
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    phase = sys.argv[1]
    records = []
    start, wall = cpu(), time.monotonic()
    previous = sum(json.loads(p.read_text())['numerical_cpu_seconds'] for p in O.glob('receipt-*.json'))
    build_cpu = 0.
    status = 'RUNNING'
    def call(label, command, build=False, expected_rc=0):
        nonlocal build_cpu
        used = previous + cpu() - start - build_cpu
        if not build and used + 90 > 900:
            raise RuntimeError('initial numerical CPU ceiling reserve reached')
        def cap():
            resource.setrlimit(resource.RLIMIT_CPU, (90, 91))
        a, w = cpu(), time.monotonic()
        out, err = O / (label + '.json'), O / (label + '.err')
        with out.open('x') as f, err.open('x') as e:
            q = subprocess.run(list(map(str, command)), stdout=f, stderr=e,
                               timeout=120, preexec_fn=cap)
        elapsed = cpu() - a
        if build:
            build_cpu += elapsed
        records.append(dict(label=label, command=list(map(str, command)),
                            cpu_seconds=elapsed, wall_seconds=time.monotonic()-w,
                            returncode=q.returncode, expected_returncode=expected_rc,
                            stdout_sha256=sha(out), stderr_sha256=sha(err), build=build))
        if q.returncode != expected_rc:
            raise RuntimeError(label + ' failed')
        return json.loads(out.read_text(), parse_float=decimal.Decimal) if not build and expected_rc == 0 else None
    try:
        if phase == 'benchmark':
            call('build', [TC/'bin/g++', '-O3', '-std=c++17', '-ffp-contract=off',
                          '-Wall', '-Wextra', '-isystem', TC/'install/Release/Zaki/include',
                          S/'tests/bnv/x2_uncertainty.cpp', '-o', O/'assessment',
                          TC/'install/Release/Zaki/lib/libZaki.a'], build=True)
        exe = O/'assessment'
        def solve(i, label, n=512, width=None, sign=0, retired=0):
            target=json.loads((S/'docs/validation/x1-oracle-anchor/targets.json').read_text())[str(i)]
            with decimal.localcontext() as ctx:
                ctx.prec=60
                hc=json.loads((R/f'x2-fixed-step/t{i}-n4096.json').read_text(),
                              parse_float=decimal.Decimal)['values'][8]
                if width is not None:
                    hc *= 1 + decimal.Decimal(width)*sign
            return call(label, [exe, str(hc), target['surfaceP_cgs'], n, '1e-4', 96, retired],
                        expected_rc=1 if retired else 0)
        if phase == 'benchmark':
            replay=solve(0,'t0-replay512')
            retained=json.loads((R/'ba13-headroom/t0-n512.json').read_text(),parse_float=decimal.Decimal)
            if replay != retained:
                raise RuntimeError('retained decimal replay differs')
            solve(0,'retired-defect-negative',retired=1)
            solve(0,'t0-n1024',1024)
        elif phase == 'm1':
            estimate=json.loads((O/'receipt-benchmark.json').read_text())
            if estimate['status']!='PASS':
                raise RuntimeError('benchmark gate')
            for i in (0,10,20):
                if i:
                    solve(i,f't{i}-n1024',1024)
                solve(i,f't{i}-n2048',2048)
                for sign in (-1,1):
                    solve(i,f't{i}-h0.000025-{sign}',512,'0.000025',sign)
        else:
            raise ValueError(phase)
        status='PASS'
    finally:
        receipt=dict(phase=phase,status=status,job=os.environ['SLURM_JOB_ID'],
                     affinity=sorted(os.sched_getaffinity(0)),
                     numerical_cpu_seconds=cpu()-start-build_cpu,build_cpu_seconds=build_cpu,
                     wall_seconds=time.monotonic()-wall,
                     peak_rss_kib=max(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                                      resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss),
                     runs=records,source={str(p.relative_to(S)):sha(p) for p in
                     (pathlib.Path(__file__),S/'tests/bnv/x2_uncertainty.cpp',S/'tests/bnv/x2_fixed_step.cpp',
                      S/'docs/validation/ba13-closure/PREDECLARATION.md')},
                     executable_sha256=sha(O/'assessment') if (O/'assessment').exists() else None)
        with (O/f'receipt-{phase}.json').open('x') as f:
            json.dump(receipt,f,indent=2,allow_nan=False)
        print(json.dumps(receipt,indent=2))

if __name__=='__main__':
    main()
