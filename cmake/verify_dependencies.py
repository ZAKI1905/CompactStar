#!/usr/bin/env python3
"""Authenticate the locally qualified Mac packages before CMake loads them.

The lock is evidence for this candidate, not a claim of portable archive hashes.
Rebuilt packages require explicit requalification and a reviewed lock update.
"""
import argparse
import hashlib
import json
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def manifest(prefix, paths):
    return hashlib.sha256(''.join(
        f'{digest(p)}  {p.relative_to(prefix)}\n' for p in sorted(paths) if p.is_file()
    ).encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=['Debug', 'Release'], required=True)
    parser.add_argument('--zaki', type=Path, required=True)
    parser.add_argument('--confind', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    lock = json.loads(Path(__file__).with_name('dependency-lock.json').read_text())
    records = []
    for name, supplied in [('Zaki', args.zaki), ('CONFIND', args.confind)]:
        if not supplied.is_absolute():
            raise SystemExit(f'{name}: an absolute explicit prefix is required')
        prefix = supplied.resolve(strict=True)
        expected = next(r for r in lock if r['name'] == name and r['mode'] == args.mode)
        actual = {
            'archive_sha256': digest(prefix / 'lib' / f'lib{name}.a'),
            'headers_sha256': manifest(prefix, (prefix / 'include').rglob('*')),
            'configs_sha256': manifest(prefix, (prefix / 'lib/cmake').rglob('*')),
        }
        for key, value in actual.items():
            if value != expected[key]:
                raise SystemExit(f'{name} {args.mode}: {key} mismatch: {value}')
        provenance = json.loads((prefix / 'provenance.json').read_text())
        for key in ('source_sha', 'mode', 'dirty'):
            if provenance[key] != expected[key]:
                raise SystemExit(f'{name}: provenance {key} mismatch')
        records.append({**expected, **actual, 'prefix': str(prefix)})
    args.output.write_text(json.dumps(records, indent=2) + '\n')


if __name__ == '__main__':
    main()
