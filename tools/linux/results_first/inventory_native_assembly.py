#!/usr/bin/env python3
"""Read-only assembly-input inventory; never constructs a star or runs an RHS.

This is a reproducible evidence collector, NOT an experiment builder or a
qualification assessor. A matching hash does not reconstruct a semantic owner.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path


def sha(path):
    with path.open('rb') as stream:
        digest = hashlib.sha256()
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    evidence = args.source / 'docs/validation/results-first-bnv-pilot'
    pins = json.loads((evidence / 'evidence/preflight.json').read_text())['input_sha256']
    raw = Path(next(p for p in pins if p.endswith('/raw/coefficients.tsv'))).parent
    cluster = raw.parents[2]
    binding = raw.parent / 'fit/binding.json'
    checks = []

    def check(path, expected):
        actual = sha(path) if path.is_file() else None
        checks.append(dict(path=str(path), expected_sha256=expected,
                           actual_sha256=actual, matches=actual == expected))

    for path, expected in pins.items():
        check(Path(path), expected)
    # The original native fitter retained the input identities for its aggregate
    # tangent/number tables. Record that provenance without calling it a loader.
    bound = json.loads(binding.read_text())
    for name in ('tangent.tsv', 'target_stars.tsv', 'thermal_surface.tsv'):
        check(raw / name, bound['inputs'][name])

    directories = []
    roots = [cluster / 'ba13-recovery' / run / 'stars/tangent-sequences/target-0'
             for run in ('run1', 'repeat')]
    for run in ('native-current', 'native-verify'):
        owner = cluster / 'results-first-bnv-preparation-20261108' / run / 'native-fixture/owning-star'
        roots.extend((owner / 'central', owner / 'sequence'))
    for root in roots:
        files = sorted(p for p in root.rglob('*') if p.is_file()) if root.exists() else []
        directories.append(dict(path=str(root), exists=root.is_dir(),
                                child_directories=sum(p.is_dir() for p in root.rglob('*')),
                                files=[dict(path=str(p.relative_to(root)), sha256=sha(p)) for p in files]))

    tables = {}
    for name in ('coefficients.tsv', 'tangent.tsv', 'target_stars.tsv', 'thermal_surface.tsv'):
        with (raw / name).open() as stream:
            rows = csv.DictReader(stream, delimiter='\t')
            tables[name] = dict(columns=rows.fieldnames, first_row=next(rows))
    profile = raw / 'chemical-profiles/target-0/profile.tsv'
    with profile.open() as stream:
        tables['profile.tsv'] = dict(columns=stream.readline().strip().split('\t'),
                                     rows=sum(1 for _ in stream))
    summary_path = evidence / 'endpoint-acceptance/evidence/assessment/summary.json'
    summary = json.loads(summary_path.read_text())
    result = dict(
        schema='compactstar.native-assembly.input-inventory.v1',
        scope='Input authentication and retained serialization coverage only',
        input_hashes_match=all(c['matches'] for c in checks), checks=checks,
        fitter_binding=dict(path=str(binding), sha256=sha(binding)),
        owner_directories=directories, retained_tables=tables,
        existing_static_result=dict(path=str(summary_path), sha256=sha(summary_path),
                                    status=summary['status'], row_statuses=summary['row_statuses'],
                                    maximum_utilization=summary['maximum_utilization']),
        source_snapshot=None, control_snapshot=None,
        assembly_status=('BLOCKED_EMPTY_KNOWN_OWNER_DIRECTORIES'
                         if all(not row['files'] for row in directories)
                         else 'REVIEW_REQUIRED_OWNER_FILES_PRESENT'),
        limitation='Directory inventory is not proof that no other archive exists. '
                   'Aggregate tangent rows are not Phase-5B derivative snapshots.',
        stellar_solves=0, rhs_evaluations=0, thermal_trajectories=0)
    # Exclusive creation prevents overwriting retained evidence.
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(result['assembly_status'])
    print('input_hashes_match=' + str(result['input_hashes_match']))
    # Coverage is an observed blocked disposition, never a readiness PASS.
    return 2


if __name__ == '__main__':
    raise SystemExit(main())
