#!/usr/bin/env python3
"""Independent serialized-ledger audit; no assessor imports and no solver calls."""
import argparse
import copy
from decimal import Decimal, localcontext
import gzip
import json
import os
from pathlib import Path


def load(path):
    with (gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)) as f:
        return json.load(f)


def validate(ledger, original, comparison, violations):
    checks = 0
    def require(value, why):
        nonlocal checks
        checks += 1
        if not value:
            raise ValueError(why)
    require(len(ledger['rows']) == 462, 'row count')
    require(len({(r['index'], r['quantity']) for r in ledger['rows']}) == 462, 'row uniqueness')
    require(len(violations) == 353, 'expectation count')
    require(len({(r['index'], r['quantity']) for r in violations}) == 353, 'expectation uniqueness')
    expected_violations = set()
    for r, old in zip(ledger['rows'], original['rows']):
        preserved = {k: v for k, v in r.items() if k != 'static_completion'}
        require(preserved == old, 'historical record changed')
        n, i = r['quantity'], r['index']
        s = r['static_completion']
        if n not in comparison:
            require(s['status'] == {'metric_structure': 'CONSUMER_SCOPE_ONLY', 'P2_average': 'EXACT_IDENTITY',
                                    'species_support': 'CONDITIONALLY_QUALIFIED'}[n], 'scope/identity disposition')
            continue
        terms = {t['name']: t for t in s['terms']}
        require(len(terms) == len(s['terms']), 'duplicate ledger term')
        sources = set(ledger['inputs']) | set(ledger['derived_sources'])
        require(all(set(t['source_ids']) <= sources for t in s['terms']+s['baseline']), 'unresolved provenance')
        require(terms['moving_cut_input_remainder']['value'] is None, 'unsupported cut assigned a number')
        require(s['full_utilization'] is None and s['full_budget'] is None, 'partial promoted to full')
        require(s['status'] == 'UNQUALIFIED', 'missing cut did not refuse')
        d = comparison[n]
        e = 10*(9.1e-12 if n == 'Ltilde_Me' else max(abs(v) for a in d['signed_error'] for v in a))
        require(terms['unseen_reference_bias_increment']['value'] == (0. if i == 0 else 2*e), '2e changed')
        base = {t['name']: t for t in s['baseline']}
        if n.startswith('t_'):
            require(not base and terms['full_numerator_normalization']['value'] == 0, 'absolute tangent normalization')
            require(all(k in terms for k in ('tangent_second_width_tail', 'tangent_moving_cut_path',
                                            'tangent_denominator_remainder')), 'C-TAN incomplete')
        else:
            require(base['unseen_reference_bias']['value'] == e, 'baseline e changed')
            require(base['moving_cut_at_B0']['value'] is None, 'missing baseline treated as zero')
            require(terms['full_numerator_normalization']['value'] is None, 'incomplete normalization populated')
        with localcontext() as context:
            context.prec = 60
            N = Decimal.from_float(r['historical']['drift_bound']) + sum(
                (Decimal.from_float(t['value']) for k, t in terms.items()
                 if t['value'] is not None and k != 'full_numerator_normalization'), Decimal(0))
            b = sum((Decimal.from_float(t['value']) for t in base.values() if t['value'] is not None), Decimal(0))
            partial = N/(1-b)
            expected = float(partial/Decimal.from_float(r['historical']['threshold']))
        require(abs(expected-s['partial_utilization']) <= 8e-15*max(abs(expected), 1e-300), 'full-numerator accounting mismatch')
        if s['amended_cref_partial_utilization']-r['partial_utilization'] > 1e-3:
            expected_violations.add((i, n))
        require(s['partial_utilization'] < 1, 'supported allowance exceeds unchanged threshold')
        if n in ('envelope', 'Tsurface_inf'):
            require(s['amended_cref_partial_status'] == 'SATISFIES_NECESSARY_TEST', 'amended partial failed')
            require('total_reference_fit_response' in terms and 'unseen_reference_fit_response' not in terms, 'wrong fit response')
        else:
            require(terms['unseen_reference_fit_response']['value'] == d['fit_operator_norm']*e, 'other fit changed')
        require('fit_nonlinear_input_remainder' in terms, 'C-INPUT missing')
    require(expected_violations == {(v['index'], v['quantity']) for v in violations}, '353 row identities mismatch')
    require(ledger['summary']['operational_certificate_issued'] is False, 'certificate issued')
    return checks


def main():
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Slurm required')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--assessment', type=Path, required=True)
    a = parser.parse_args()
    repo = Path(__file__).resolve().parents[3]
    ledger = load(a.assessment/'ledger.json.gz')
    original = load(repo/'docs/validation/ba13-final-decision/evidence/matrix.json.gz')
    comparison = load(repo/'docs/validation/ba13-evidence/evidence/comparison.json.gz')
    violations = load(a.assessment/'amended-expectation-violations.json')
    count = validate(ledger, original, comparison, violations)
    # Mutations are memory-only copies: no retained file is modified.
    mutations = {}
    for label in ('drop-cut', 'half-increment', 'drop-denominator', 'promote-partial', 'historical-threshold'):
        bad = copy.deepcopy(ledger)
        row = next(r for r in bad['rows'] if r['index'] == 20 and r['quantity'] == 'envelope')
        s = row['static_completion']
        if label == 'drop-cut':
            next(t for t in s['terms'] if t['name'] == 'moving_cut_input_remainder')['value'] = 0.
        elif label == 'half-increment':
            next(t for t in s['terms'] if t['name'] == 'unseen_reference_bias_increment')['value'] /= 2
        elif label == 'drop-denominator':
            next(t for t in s['baseline'] if t['name'] == 'unseen_reference_bias')['value'] = 0.
        elif label == 'promote-partial':
            s['full_utilization'] = s['partial_utilization']
        else:
            row['historical']['threshold'] *= 2
        try:
            validate(bad, original, comparison, violations)
        except ValueError as error:
            mutations[label] = str(error)
        else:
            raise RuntimeError('mutation escaped: '+label)
    result = dict(status='PASS_LEDGER_AUDIT_NOT_QUALIFICATION', checks=count,
                  mutations_refused=mutations, mutation_scope='in-memory copies only',
                  stellar_calculations=0, trajectories_run=0)
    (a.assessment/'independent-ledger-audit.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
