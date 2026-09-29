from pathlib import Path
import subprocess,json,hashlib,re
q=Path(__file__).resolve().parent
src=Path('/Users/keeper/Documents/CompactStar/worktrees/CompactStar-external-deps-taskmanager-migration')
canonical='812463ac9ed374f64ac9cadd500066ab723d3a6c'
def git(*args):return subprocess.check_output(['git',*args],cwd=src)
def digest(b):return hashlib.sha256(b).hexdigest()
rows=[]
paths=git('ls-tree','-r','--name-only',canonical,'tests/baselines','dependencies').decode().splitlines()
paths+=git('ls-tree','-r','--name-only',canonical,'docs/adr').decode().splitlines()
paths=[p for p in paths if p.startswith(('tests/baselines/','dependencies/')) or 'ADR-0018' in p]
for rel in paths:
 old=git('show',canonical+':'+rel);actual=(src/rel).read_bytes()
 rows.append(dict(path=rel,expected=digest(old),actual=digest(actual)))
 assert old==actual,rel
plan=(src/'docs/validation/PHASE5D1_FROZEN_CONTEXT_PROVENANCE_PLAN.md').read_text()
protected=re.findall(r'^\| `([0-9a-f]{64})` \| `([^`]+)` \|$',plan,re.M)
assert len(protected)==33
for expected,rel in protected:
 actual=digest((src/rel).read_bytes());assert actual==expected,rel
 rows.append(dict(path=rel,expected=expected,actual=actual))
baseline=json.loads((src/'tests/baselines/phase5d1_controlled_evolution.json').read_text())
for rel,expected in baseline['source_provenance']['scientific_production_source_hashes'].items():
 actual=digest((src/rel).read_bytes());assert actual==expected,rel
 rows.append(dict(path=rel,expected=expected,actual=actual))
p='docs/validation/PHASE6_EXTERNAL_DEPS_TASKMANAGER_MIGRATION_PREDECLARATION.md'
original=git('show','6ab1783b8003295224dcf2d720d500a174f9ef9e:'+p)
assert (src/p).read_bytes().startswith(original)
result=dict(candidate_head=git('rev-parse','HEAD').decode().strip(),canonical=canonical,predeclaration_original_body_sha256=digest(original),predeclaration_body_unchanged=True,rows=rows,all_equal=True)
(q/'evidence/protected-authority.json').write_text(json.dumps(result,indent=2)+'\n')
print('PASS',len(rows),'protected entries; original predeclaration unchanged')
