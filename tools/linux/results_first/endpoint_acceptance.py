#!/usr/bin/env python3
"""Restricted conditional C-CUT accounting; no stellar or evolution entrypoint."""
import argparse
import copy
from collections import Counter
from decimal import Decimal, localcontext
import gzip
import hashlib
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT/'docs/validation/results-first-bnv-pilot'
DOC = BASE/'endpoint-acceptance'
WITNESS = BASE/'surface-transfer/evidence/assessment'
CONDITIONS = ['S', 'C-REF', 'C-TAN', 'C-CONT', 'C-INPUT', 'C-SUPPORT',
              'C-CUT', 'C-MUN', 'ADR-0022', 'C-ENDPOINT', 'METRIC-CONSUMER-SCOPE']

def read(path):
    with (gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)) as f:
        return json.load(f)

def D(value):
    return Decimal(str(value))

def require(ok, message):
    if not ok:
        raise ValueError(message)

def authenticate(registry, bounds):
    expected = dict(schema='compactstar.restricted-endpoint-acceptance.v1',
        status='ACCEPTED_CONDITIONAL_STATIC_ONLY', scope='RF-BNV-P2-20261108',
        starting_commit='afd4b06074a5f597ef4cbe38ce3c9cb04386f556',
        condition='C-ENDPOINT', maximum_extension_x='3.022e-26',
        physical_interval_x=['-1e-6', '0'], rigorous_endpoint_bound=False,
        paired_trajectories_authorized=False, historical_certificate_promoted=False,
        new_stellar_calculations_authorized=False)
    require(all(registry.get(k) == v for k, v in expected.items()), 'wrong or broadened authority')
    for p, digest in registry['files'].items():
        require(hashlib.sha256((ROOT/p).read_bytes()).hexdigest() == digest, 'authority hash: '+p)
    require(0 < D(bounds['eta']) <= D(registry['maximum_extension_x']), 'extension exceeds authority')
    require(bounds['endpoint_extension_supported'] is False, 'historical endpoint result changed')
    inputs = read(DOC/'inputs.json')
    for p, digest in inputs.items():
        require(hashlib.sha256((ROOT/p).read_bytes()).hexdigest() == digest, 'input hash: '+p)

def transition_check(bounds, retained, continuation):
    """Falsification from existing samples, explicitly not an endpoint proof."""
    states = []
    secants = []
    for t in ('0', '20'):
        states += list(retained['centers'][t].values())
        for kind in ('widths', 'grids'):
            for pair in retained[kind][t].values():
                states += list(pair.values())
                lo, hi = pair['-1']['primary'], pair['1']['primary']
                require(D(hi[8]) > D(lo[8]), 'reversed central coordinate')
                secants.append((D(hi[0])-D(lo[0]))/(D(hi[8])-D(lo[8])))
    margins = []
    for s in states:
        hs, hn, hb, hm, hc = [D(s[k]) for k in ('Hs', 'Hn', 'Hb', 'Hmu')] + [D(s['primary'][8])]
        require(0 < hs < hn < hb < hm < hc, 'known endpoint onset/topology contradiction')
        margins.append(dict(terminal_neutron=hn-hs, central_muon=hc-hm))
    require(min(secants) > 0, 'nonpositive sampled baryon derivative')
    for k in ('pressure_denominator', 'fixed_B_derivative_denominator', 'chemical_denominator',
              'horizon_denominator', 'pressure_lower', 'topology_margin'):
        require(D(bounds[k]) > 0, 'terminal domain contradiction: '+k)
    w = continuation['window']
    require(all(v['constant'] for v in w['selection'].values()), 'known sampled selection transition')
    require(w['neutron_edge_sample_monotone'] and w['envelope_sample_monotone'], 'sample monotonicity contradiction')
    return dict(status='NO_KNOWN_CONTRADICTION_IN_RETAINED_EVIDENCE',
        endpoint_states_checked=len(states), positive_baryon_secants=len(secants),
        minimum_sampled_B_H=str(min(secants)),
        minimum_sampled_terminal_neutron_margin=str(min(m['terminal_neutron'] for m in margins)),
        minimum_sampled_central_muon_margin=str(min(m['central_muon'] for m in margins)),
        terminal_tube_topology_margin=bounds['topology_margin'],
        continuous_selection_proven=w['continuous_selection_proven'],
        rigorous_endpoint_bound=False, endpoint_grid_diagnostics='SIX_INCONCLUSIVE_RETAINED',
        interpretation='Retained samples and terminal tube show no known transition; unsampled family extension is C-ENDPOINT, not a numerical theorem.')

def complete(old, bounds, comparison, m4):
    new = copy.deepcopy(old)
    maximum = Decimal(0)
    max_delta = Decimal(0)
    group_max = {}
    for row in new['rows']:
        name, i = row['quantity'], row['index']
        prior = row['static_completion']
        if name not in bounds['allowances']:
            row['endpoint_acceptance'] = dict(status=row['surface_transfer']['status'],
                full_utilization=row['surface_transfer']['full_utilization'],
                conditions=CONDITIONS, rigorous_bound=False)
            continue
        U = D(bounds['allowances'][name]['interior_point_bound'])
        require(U > 0 and U.is_finite(), 'unsupported or nonpositive cut allowance')
        fit = D(comparison[name]['fit_operator_norm'])
        extra = (fit+(0 if i == 0 else 2))*U
        endpoint = Decimal(0)
        if i == 20:
            xx = list(map(D, m4['achieved_x']))
            endpoint = 2*U*max(Decimal(0), D('1e-6')+xx[20])/abs(xx[20]-xx[10])
        baseline = Decimal(0) if name.startswith('t_') else U
        N, b = D(prior['known_numerator']), D(prior['known_baseline'])
        denominator = 1-b-baseline
        require(denominator > 0, 'unseparated baseline')
        budget = (N+extra+endpoint)/denominator
        utilization = budget/D(row['historical']['threshold'])
        delta = (budget-N/(1-b))/D(row['historical']['threshold'])
        row['endpoint_acceptance'] = dict(
            status='CONDITIONALLY_QUALIFIED' if utilization <= 1 else 'FAIL_CONDITIONAL_BUDGET',
            conditions=CONDITIONS, rigorous_bound=False,
            full_interval_point_bound=str(U), baseline_cut_bound='NOT_APPLICABLE' if name.startswith('t_') else str(U),
            additional_numerator=str(extra), additional_endpoint=str(endpoint), additional_baseline=str(baseline),
            full_denominator=str(denominator), full_budget=str(budget), full_utilization=str(utilization),
            utilization_increase=str(delta), evidence='validated local witness plus owner-accepted C-ENDPOINT')
        maximum = max(maximum, utilization)
        group_max[name] = str(max(D(group_max.get(name, '0')), utilization))
        if delta > max_delta:
            max_delta = delta
            impact = dict(quantity=name, index=i)
    counts = dict(Counter(r['endpoint_acceptance']['status'] for r in new['rows']))
    new['endpoint_acceptance_summary'] = dict(
        status='CONDITIONAL_STATIC_PASS' if counts == {'CONDITIONALLY_QUALIFIED':420, 'EXACT_IDENTITY':21, 'CONSUMER_SCOPE_ONLY':21} else 'FAIL_STATIC',
        row_statuses=counts, maximum_utilization=str(maximum), maximum_cut_utilization_increase=str(max_delta),
        maximum_cut_increase_row=impact, group_max_utilization=group_max,
        conditions=CONDITIONS, rigorous_bound=False, historical_CREF_failures=42,
        original_expectation_violations=355, amended_expectation_violations=353,
        historical_ledger_preserved=True, historical_certificate_promoted=False,
        nonlinear_trajectory_qualified=False, trajectories_authorized=False,
        all_restricted_pilot_gates_pass=False, actual_path_uncertainty=None)
    return new

def main():
    require(bool(os.environ.get('SLURM_JOB_ID')), 'EKU Slurm required')
    p = argparse.ArgumentParser(); p.add_argument('--output', type=Path, required=True); args=p.parse_args()
    args.output.mkdir()
    bounds = read(WITNESS/'bounds.json')
    authenticate(read(DOC/'acceptance-registry.json'), bounds)
    with localcontext() as ctx:
        # Delta utilization subtracts two O(1) budgets to resolve effects down
        # to O(1e-26). Keep enough guard digits for the independent 95-digit
        # relative arithmetic check; this does not change a scientific budget.
        ctx.prec = 160
        transition = transition_check(bounds, read(BASE/'cut-assessment/retained-inputs.json.gz'),
                                      read(ROOT/'docs/validation/ba13-evidence/evidence/continuation.json'))
        result = complete(read(WITNESS/'ledger.json.gz'), bounds,
                          read(ROOT/'docs/validation/ba13-evidence/evidence/comparison.json.gz'),
                          read(ROOT/'docs/validation/ba13-closure/evidence/m4.json'))
    for name, data in [('ledger.json.gz',result), ('summary.json',result['endpoint_acceptance_summary']),
                       ('endpoint-transition-audit.json',transition)]:
        data = (json.dumps(data, indent=2)+'\n').encode()
        (args.output/name).write_bytes(gzip.compress(data, mtime=0) if name.endswith('.gz') else data)
    with (args.output/'matrix.tsv').open('w') as f:
        f.write('index\tquantity\tstatus\tfull_utilization\tconditions\n')
        for row in result['rows']:
            s=row['endpoint_acceptance']
            f.write('\t'.join(map(str,(row['index'],row['quantity'],s['status'],s['full_utilization'],','.join(s['conditions']))))+'\n')
    print(json.dumps(result['endpoint_acceptance_summary']))

if __name__ == '__main__':
    main()
