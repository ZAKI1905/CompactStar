import pathlib, subprocess, hashlib, json
S=pathlib.Path('/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/ba13-final-decision'); R=S.parents[1]
def git(*a):return subprocess.check_output(['git','-C',str(S),*a],text=True).strip()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
base='84e1062db78b91e4fb0988074941b1aea6019917'
branches=[]
for line in git('for-each-ref','--format=%(refname:short) %(objectname)','refs/heads').splitlines():
 n,h=line.split(); rc=subprocess.run(['git','-C',str(S),'merge-base','--is-ancestor',h,base]).returncode
 branches.append(dict(branch=n,sha=h,ancestor_of_entry=rc==0));assert rc==0,line
records=[]
def manifest(path):
 d=json.loads(path.read_text()); root=pathlib.Path(d['root']);bad=[]
 for rel,v in d['files'].items():
  p=root/rel; h=v['sha256'] if isinstance(v,dict) else v
  if not p.exists() or sha(p)!=h:bad.append(rel)
 records.append(dict(path=str(path),sha256=sha(path),entries=len(d['files']),bad=bad));assert not bad,bad
for group in ['ba13-evidence','ba13-closure','ba13-headroom','x2-fixed-step','x1-oracle-anchor']:
 manifest(S/f'docs/validation/{group}/evidence/external-manifest.json')
inputs=S/'docs/validation/ba13-closure/evidence/inputs.json'
for name,d in json.loads(inputs.read_text()).items():
 if not isinstance(d,dict) or 'files' not in d:continue
 bad=[n for n,v in d['files'].items() if not (pathlib.Path(d['root'])/n).exists() or sha(pathlib.Path(d['root'])/n)!=(v['sha256'] if isinstance(v,dict) else v)]
 records.append(dict(path=str(inputs),group=name,entries=len(d['files']),sha256=sha(inputs),bad=bad));assert not bad,bad
for root,p in [(R/'ba13-closure-independent-review-20261020',R/'ba13-closure-independent-review-20261020/SHA256SUMS'),(S,S/'docs/validation/ba13-evidence/evidence/SHA256SUMS'),(S/'docs/validation/x1-oracle-anchor',S/'docs/validation/x1-oracle-anchor/SHA256SUMS')]:
 lines=p.read_text().splitlines();bad=[]
 for line in lines:
  h,n=line.split(maxsplit=1); n=n.lstrip('*');q=root/n
  if not q.exists() or sha(q)!=h:bad.append(n)
 records.append(dict(path=str(p),sha256=sha(p),entries=len(lines),bad=bad));assert not bad,bad
cert=R/'ba13-recovery/run1/fit/certificate.json'; fitter=S/'tests/bnv/compare_candidate.py'
assert sha(fitter)=='9b3ef1ce17df7e1f8d2d7c40e111ee0844abfddc2cc54b8f93cb3f2ce1322370'
checks={str(p):sha(p) for p in [cert,fitter,S/'docs/validation/PHASE6A1_CONTROLLED_BNV_IMPLEMENTATION_PREFLIGHT.md',S/'docs/adr/ADR-0021-ba13-static-uncertainty-completion.md',R/'ba13-closure-independent-review-20261020/REVIEW.md']}
out=dict(entry=base,entry_object_type=git('cat-file','-t',base),entry_tree=git('rev-parse',base+'^{tree}'),entry_matches_named_branch=git('rev-parse','candidate/linux-ba13-evidence-20261026')==base,worktree=str(S),branch=git('branch','--show-current'),initial_status=git('status','--porcelain'),branches=branches,manifests=records,identities=checks,authentication_scope='Git object/ancestry and SHA-256 content authentication; no cryptographic author signature asserted',latest_available_independent_review='2026-10-20 review of 093d121; no independent review of 84e1062 located')
p=S/'docs/validation/ba13-final-decision/evidence';p.mkdir(parents=True)
(p/'authentication.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(dict(manifests=len(records),entries=sum(r['entries'] for r in records),all_verified=True,identities=checks),indent=2))
