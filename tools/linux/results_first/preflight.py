#!/usr/bin/env python3
"""Read-only identity/consumer gate audit. Never constructs a star or runs an ODE."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import time

SOURCE = 'bb0f9c5af8ebaf86b5a59aece5551bb9dafbd264'
ROOT = Path('/mnt/sdd/zaki/CompactStar/linux-transition-20261010')
REPO = Path(__file__).resolve().parents[3]


def git(*args):
    return subprocess.check_output(['git', '-C', str(REPO), *args], text=True).strip()


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    start = time.monotonic()
    review = ROOT / 'ba13-final-decision-independent-review-20261030'
    manifest_checks = []
    for line in (review / 'SHA256SUMS').read_text().splitlines():
        expected, name = line.split(maxsplit=1)
        path = review / name.lstrip('*')
        actual = digest(path)
        if actual != expected:
            raise RuntimeError(f'review manifest mismatch: {path}')
        manifest_checks.append(dict(path=str(path), sha256=actual))
    native = ROOT / 'ba13-recovery/run1'
    expected = {
        'fit/certificate.json': '4212cad05bb4244ab75e80b8dfea023e6d0fcee69086c9993bd48d26a40e41fb',
        'fit/certificate.tsv': '38bc85c62621b390ff2d7c8affee1f5ba4e9f0f44233047bb9959f2d685fd5ba',
    }
    identities = {}
    for name, pin in expected.items():
        actual = digest(native / name)
        if actual != pin:
            raise RuntimeError(f'native certificate mismatch: {name}')
        identities[str(native / name)] = actual
    with (native / 'fit/certificate.tsv').open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    if len(rows) != 462:
        raise RuntimeError('native certificate incomplete')
    b0 = float(rows[0]['B_solved_count'])
    context = git('show', SOURCE + ':CompactStar/Physics/Rotochemical/FrozenRotochemicalRunContext.hpp')
    mac_pins = {
        'profile.tsv': 'e9cd03b0b8449806f6c9883d75de1d3dff0cf1481675efc45a56519655d40890',
        'model.txt': '3ea70de79e15b70c5a6d68f48335d18047ff80e60b55a9acdb78084e9be4d6d4',
        'certificate.txt': '7fc892b50bfbcdd2c5d963e0ff866777f3628f40fc282e1ce5a0b0364200a453',
        'freegas.tsv': '7cd44c92e1e7206e0e68e3fed7e3f0ca68e79ab4517d02b96ff78b9be23d3f1a',
    }
    binding = []
    profile = native / 'raw/chemical-profiles/target-0'
    for name, pin in mac_pins.items():
        if pin not in context:
            raise RuntimeError('audited context pins changed')
        actual = digest(profile / name)
        identities[str(profile / name)] = actual
        binding.append(dict(file=name, native_sha256=actual, context_pin=pin,
                            matches=actual == pin))
    for path in sorted((profile / 'thermal').iterdir()):
        if path.is_file():
            identities[str(path)] = digest(path)
    identities[str(native / 'raw/coefficients.tsv')] = digest(native / 'raw/coefficients.tsv')
    branches = []
    for ref in git('for-each-ref', '--format=%(refname)', 'refs/heads').splitlines():
        code = subprocess.run(['git', '-C', str(REPO), 'merge-base', '--is-ancestor', ref, SOURCE]).returncode
        branches.append(dict(ref=ref, commit=git('rev-parse', ref), ancestor_of_source=code == 0))
    registry = json.loads(git('show', SOURCE + ':docs/validation/ba13-final-decision/acceptance-registry.json'))
    stat = os.statvfs(ROOT)
    result = dict(
        status='AUTHENTICATED_NOT_AUTHORIZED_FOR_TRAJECTORIES',
        scientific_disposition='UNQUALIFIED_PENDING_OWNER_AND_NATIVE_RUN_BINDING',
        host=platform.node(), system=platform.system(), repository=str(REPO),
        source_commit=SOURCE, audit_head=git('rev-parse', 'HEAD'),
        branch=git('branch', '--show-current'), branches=branches,
        worktrees=git('worktree', 'list', '--porcelain'),
        storage=dict(root=str(ROOT), available_bytes=stat.f_bavail * stat.f_frsize,
                     writable=os.access(ROOT, os.W_OK),
                     scriptorium_visible_on_this_host=Path('/Volumes/Scriptorium').exists(),
                     mac_mount_status='NOT_OBSERVABLE_FROM_THIS_LINUX_SESSION'),
        review_manifest_checks=manifest_checks, input_sha256=identities,
        native_context_hash_comparison=binding,
        native_B0_count=b0, native_Bdot_count_s=-1e-12 * b0 / (365.25 * 86400),
        historical_Bdot_literal_count_s=-2.4136520263641375e37,
        native_Bdot_matches_historical_literal=(-1e-12*b0/(365.25*86400) == -2.4136520263641375e37),
        bounded_duration_s=462269531250.0, bounded_duration_year=462269531250.0/(365.25*86400),
        predicted_fractional_depletion=1e-12*462269531250.0/(365.25*86400),
        historical_native_max_utilization=max(float(r['utilization']) for r in rows),
        existing_acceptance_registry=registry,
        trajectories_run=0, slurm_jobs_submitted=0,
        elapsed_wall_seconds=time.monotonic()-start,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as f:
        json.dump(result, f, indent=2)
        f.write('\n')
    print(json.dumps({k: result[k] for k in ('status','native_Bdot_count_s',
          'bounded_duration_year','predicted_fractional_depletion')}))
    print('Review manifest entries verified:', len(manifest_checks))
    print('Native files differing from fixed Mac context pins:', sum(not r['matches'] for r in binding))


if __name__ == '__main__':
    main()
