"""Read-only authentication of the onset pilot candidate and retained evidence."""
import collections
import glob
import hashlib
import json
import subprocess
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, PILOT, write

WT = ROOT / 'worktrees' / 'ba13-onset-pilot'
DOC = 'docs/validation/ba13-onset-pilot/'


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def git(*args):
    return subprocess.check_output(['git', '-C', str(WT), *args])


def main():
    out = {}
    out['head'] = git('rev-parse', 'HEAD').decode().strip()
    out['worktree_clean'] = git('status', '--porcelain') == b''
    out['commits'] = git('log', '--format=%H %cI %s', '4f3cf97..HEAD').decode().splitlines()
    out['frozen_files_unchanged'] = {}
    for commit, name in (('8bb72e2', 'PREDECLARATION.md'), ('e657efb', 'evidence/predictions.json'), ('650871e', 'evidence/eos-predictions.json')):
        blob = git('show', f'{commit}:{DOC}{name}')
        out['frozen_files_unchanged'][name] = dict(commit=commit, sha256=hashlib.sha256(blob).hexdigest(),
                                                   unchanged=blob == (WT / DOC / name).read_bytes())
    out['tools_changed_after_final_gate_650871e'] = [n for n in git('diff', '--name-only', '650871e', 'HEAD').decode().split()
                                                     if not n.startswith('docs/')]
    rec = [json.loads(Path(p).read_text()) for p in glob.glob(str(PILOT / 'runs' / '*.receipt.json'))]
    out['receipts'] = {k: dict(collections.Counter(r[k] for r in rec)) for k in
                       ('executable_sha256', 'predeclaration', 'predictions', 'eos_predictions', 'status', 'source_commit')}
    out['build_executable_sha256'] = sha(PILOT / 'build' / 'onset_pilot')
    m = json.loads((PILOT / 'PILOT_COMPLETE.json').read_text())
    out['pilot_complete_sha256'] = sha(PILOT / 'PILOT_COMPLETE.json')
    out['manifest_files'] = len(m['files'])
    out['manifest_bytes'] = sum(v['bytes'] for v in m['files'].values())
    out['manifest_mismatches'] = [n for n, v in m['files'].items() if sha(PILOT / n) != v['sha256']]
    out['source_files_mismatching_worktree'] = [n for n, h in m['source_files'].items() if sha(WT / n) != h]
    out['protected_mismatching_worktree'] = [n for n, h in m['protected'].items() if sha(WT / n) != h]
    out['compact_copies_identical'] = {n: sha(WT / DOC / 'evidence' / n) == sha(PILOT / n) for n in
                                       ('scores.json', 'high-precision.json', 'resource-summary.json', 'constants.json',
                                        'diagnosis.json', 'geometry-diagnostic.json')}
    hand = ROOT / 'ba13-onset-handoff-20261013'
    out['bundle_sha256'] = sha(hand / 'CompactStar-ba13-onset-incremental.bundle')
    out['bundle_sha_recorded'] = (hand / 'SHA256SUMS').read_text().split()[0]
    rev = ROOT / 'ba13-uncertainty-independent-review-20261012'
    out['previous_review_unchanged'] = all(sha(rev / ln.split('  ', 1)[1]) == ln.split('  ', 1)[0]
                                           for ln in (rev / 'SHA256SUMS').read_text().splitlines())
    times = {p: Path(p).stat().st_mtime for p in [str(PILOT / 'table-2048-11.log'), str(PILOT / 'discovery.json'),
                                                   str(PILOT / 'runs' / 'r640000-0.receipt.json'), str(PILOT / 'scores.json')]}
    out['first_controller_output_mtime'] = times[str(PILOT / 'table-2048-11.log')]
    out['final_gate_commit_time'] = git('log', '-1', '--format=%ct', '650871e').decode().strip()
    write('authentication.json', out)
    print(json.dumps({k: v for k, v in out.items() if k not in ('commits',)}, indent=1))


if __name__ == '__main__':
    main()
