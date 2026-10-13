"""Read-only authentication; output only in the descendant worktree."""
import hashlib,json,subprocess
from pathlib import Path
S=Path(__file__).resolve().parents[3]; R=S.parents[1]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
review=R/'ba13-uncertainty-independent-review-20261012'
review_hashes={}
for line in (review/'SHA256SUMS').read_text().splitlines():
 h,n=line.split('  ',1);assert sha(review/n)==h,n;review_hashes[n]=h
study=R/'ba13-uncertainty'; manifest=study/'STUDY_COMPLETE.json'
assert sha(manifest)=='b68c7a675c8bb2df0de24b37092767c75f0c1104e0ec01fc8aeafcbb642c4845'
d=json.loads(manifest.read_text())
for n,v in d['files'].items():assert sha(study/n)==v['sha256'],n
for n,h in d['source_files'].items():
 b=subprocess.check_output(['git','show',d['source_commit']+':'+n],cwd=S)
 assert hashlib.sha256(b).hexdigest()==h,n
protected=['CompactStar','cmake','tests/eos','tests/analysis','tests/bnv/compare_candidate.py','tests/bnv/frozen_sensitivity.cpp','tests/bnv/produce_frozen_certificates.py']
assert not subprocess.check_output(['git','diff','4f3cf97','--']+protected,cwd=S)
for branch in subprocess.check_output(['git','for-each-ref','--format=%(refname)','refs/heads'],cwd=S,text=True).splitlines():
 assert subprocess.run(['git','merge-base','--is-ancestor',branch,'HEAD'],cwd=S).returncode==0,branch
fitter=sha(S/'tests/bnv/compare_candidate.py');assert fitter=='9b3ef1ce17df7e1f8d2d7c40e111ee0844abfddc2cc54b8f93cb3f2ce1322370'
record=dict(status='PASS',source_entry='4f3cf97462ca030f715e26808938f781a2a40945',review=review_hashes,review_manifest=sha(review/'SHA256SUMS'),study_manifest=sha(manifest),verified_study_files=len(d['files']),study_source_commit=d['source_commit'],fitter=fitter,protected_paths=protected,certificate=sha(R/'ba13-recovery/run1/fit/certificate.json'),dependencies={n:subprocess.check_output(['git','-C',str(R/'sources'/n),'rev-parse','HEAD'],text=True).strip() for n in ['Zaki','CONFIND']},dependency_lock=sha(R/'toolchain/dependency-lock-linux.json'))
(S/'docs/validation/ba13-onset-pilot/evidence/authentication.json').write_text(json.dumps(record,indent=2)+'\n')
print('PASS',len(d['files']),'retained files; review unchanged; all other branches ancestors')
