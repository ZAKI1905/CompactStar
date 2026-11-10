#!/usr/bin/env python3
"""Read-only necessary tests of completed source evidence, never an IVP solve.

Exact binary64 lifts, 110-digit arithmetic. Missing bounds remain UNQUALIFIED.
The original native-assembly/NONLINEAR_PLAN.md supplies every threshold.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
from decimal import Decimal as D, localcontext
from pathlib import Path


def rows(p):
    with p.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def lift(s):
    return D.from_float(float(s))


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()


def analyze(run):
    source = run/'run/trajectory/source'
    files = [source/f'reference-{i}-observations.tsv' for i in range(5)]
    refs = [rows(p) for p in files]
    assert all(len(r) == 1921 for r in refs)
    assert all([r['t'] for r in a] == [r['t'] for r in refs[0]] for a in refs)
    comparisons = ((0, 1), (1, 2), (3, 4), (2, 4))
    out = {'schema': 'completed-source-necessary-tests-v1', 'precision': 110,
           'source_qualified': False, 'control_executed': False,
           'criteria': 'native-assembly/NONLINEAR_PLAN.md:95-125',
           'reference_comparisons': comparisons, 'reference_maxima': {},
           'reference_noncontracting_counts': {}, 'restarts': [], 'quadrature': [],
           'failures': [], 'missing_bounds': [
               'Independent actual-path accumulated reference and main-prefix roundoff',
               'Thermal RHS evaluation roundoff and continuous thermal enclosure',
               'Actual-output knot-aligned chemical reference: witness nodes and RHS not serialized',
               'Reference interpolation error at non-reference accepted endpoints',
               'Full R20 U_main, U_event and independently justified quadrature rounding bound',
               'Restart reference/prefix uncertainty; small raw differences are not a bound'],
           'event_scope': 'No actual main crossing; event coverage remains UNQUALIFIED'}
    for k, atol in zip(('x', 'e', 'mu'), ('1e-16', '1e-22', '1e-22')):
        maximum = None
        noncontracting = 0
        for j in range(1921):
            v = [lift(r[j][k]) for r in refs]
            diffs = [abs(v[a]-v[b]) for a, b in comparisons]
            floor = 64*lift(math.ulp(float(max(map(abs, v)))))
            lower = 2*max(diffs)+floor
            budget = D(atol)+D('1e-11')*max(map(abs, v))
            ratio = lower/(D('.05')*budget)
            contraction = diffs[1]/diffs[0] if diffs[0] else None
            noncontracting += (diffs[1] > diffs[0]/2)
            row = dict(component=k, time_s=refs[0][j]['t'], differences=diffs,
                       reference_values=v, rounding_floor_64ulp=floor,
                       reference_uncertainty_necessary_minimum=lower,
                       ULTRA_state_budget=budget, reference_allocation=D('.05')*budget,
                       allocation_ratio=ratio, rk8pd_contraction=contraction,
                       contraction_status='PASS' if diffs[1] <= diffs[0]/2 else
                       'UNQUALIFIED: independent rounding bound required')
            if maximum is None or ratio > maximum['allocation_ratio']:
                maximum = row
            if ratio > 1:
                out['failures'].append(dict(gate='REFERENCE_ALLOCATION_NECESSARY', **row))
        out['reference_maxima'][k] = maximum
        out['reference_noncontracting_counts'][k] = dict(samples=1921, count=noncontracting,
            interpretation='No failure inferred from ratio alone; independent roundoff bound is missing.')
    for path in sorted(source.glob('tier-*/restart-comparison.tsv')):
        vals = rows(path)
        assert len(vals) == 2
        tier = int(path.parent.name.split('-')[1])
        rt = [D('1e-7'), D('1e-9'), D('1e-11')][tier]
        at = [(D('1e-12'), D('1e-18')), (D('1e-14'), D('1e-20')),
              (D('1e-16'), D('1e-22'))][tier]
        components = []
        main_path = path.parent/'main.accepted_states.tsv'
        tail_path = path.parent/'restart-tail/main.accepted_states.tsv'
        main_final, tail_final = rows(main_path)[-1], rows(tail_path)[-1]
        for k, field in zip(('x', 'e', 'mu'), ('x_state', 'eta_e_MeV', 'eta_mu_MeV')):
            assert lift(vals[0][k]) == lift(tail_final[field]), 'first row must be restarted tail'
            assert lift(vals[1][k]) == lift(main_final[field]), 'second row must be uninterrupted main'
        files.extend((main_path, tail_path))
        for i, k in enumerate(('x', 'e', 'mu')):
            # stage.cpp writes tail first, uninterrupted main second.
            b, a, ref = lift(vals[0][k]), lift(vals[1][k]), lift(refs[4][-1][k])
            budget = at[0 if i == 0 else 1]+rt*max(abs(b), abs(ref))
            error = abs(b-ref)
            row = dict(component=k, uninterrupted=a, restarted=b, reference=ref,
                       restart_minus_uninterrupted=b-a, absolute_reference_error=error,
                       raw_error_ulp=error/lift(math.ulp(float(ref))), budget=budget,
                       raw_error_over_budget=error/budget,
                       status='UNQUALIFIED: reference and prefix bounds missing')
            components.append(row)
            if error > budget:
                out['failures'].append(dict(gate='RESTART_RAW_NECESSARY', path=str(path), **row))
        out['restarts'].append(dict(path=str(path), components=components))
        files.append(path)
    for path in sorted(source.glob('tier-*/R20-grids.tsv')):
        r = rows(path)
        assert [int(x['intervals']) for x in r] == [240, 480, 960, 1920]
        channels = {}
        for k in ('residual', 'delta_U', 'delta_Echem', 'outgoing_integral'):
            v = [lift(x[k]) for x in r]
            d = [abs(b-a) for a, b in zip(v, v[1:])]
            channels[k] = dict(values=v, absolute_adjacent_differences=d,
                final_ratio=d[-1]/d[-2] if d[-2] else None,
                twice_maximum_adjacent_difference=2*max(d),
                contraction_status='PASS' if d[-1] <= d[-2]/2 else
                'UNQUALIFIED: independent rounding bound required')
        out['quadrature'].append(dict(path=str(path), channels=channels,
            note='Original trapezoids retained; differences are not an independent quadrature enclosure.'))
        files.append(path)
    out['actual_output_trace_inventory'] = {
        name: len(list(source.rglob(name))) for name in
        ('local-rhs.tsv', 'witness-steps.tsv', 'witnesses.tsv')}
    out['status'] = 'FAIL_NECESSARY_ACCURACY_TEST' if out['failures'] else 'UNQUALIFIED_MISSING_REQUIRED_BOUNDS'
    out['input_hashes'] = {str(p): sha(p) for p in files}
    out['assessor_sha256'] = sha(Path(__file__))
    return out


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--skip-manifest', action='store_true', help='Reuse a separately preserved immutable source manifest')
    args = parser.parse_args()
    with localcontext() as ctx:
        ctx.prec = 110
        result = analyze(args.run)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output/'completed-tests.json').write_text(json.dumps(result, indent=2, default=str)+'\n')
    print(result['status'], 'necessary failures', len(result['failures']), flush=True)
    # Preserve all source bytes with a compressed external manifest. No mutable
    # analysis output lies below the source-run directory being inventoried.
    if args.skip_manifest:
        raise SystemExit(0)
    count = 0
    with (args.output/'external-SHA256SUMS.gz').open('wb') as raw:
        with gzip.GzipFile(fileobj=raw, mode='wb', mtime=0) as zipped:
            for p in sorted(args.run.rglob('*')):
                if p.is_file():
                    zipped.write(f'{sha(p)}  {p}\n'.encode())
                    count += 1
    (args.output/'manifest-receipt.json').write_text(json.dumps(dict(
        files=count, manifest_sha256=sha(args.output/'external-SHA256SUMS.gz'),
        source=str(args.run), scope='All completed source files; original bytes unmodified'), indent=2)+'\n')
    print('authenticated completed source inventory', count, flush=True)
