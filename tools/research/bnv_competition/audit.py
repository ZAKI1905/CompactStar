#!/usr/bin/env python3
"""Final source, documentation and evidence authentication (Python 3.11)."""
from pathlib import Path
import subprocess,hashlib,json,re,csv,ast,sys
assert sys.version_info>=(3,11), 'Use the declared Python 3.11 study environment'
root=Path(__file__).resolve().parents[3];d=root/'docs/research/bnv-rotochemical-competition';e=d/'evidence'
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(1<<20),b''):h.update(block)
 return h.hexdigest()
checks={}
stage='docs/validation/results-first-bnv-pilot/stage-b/CONDITIONAL_ACCEPTANCE.md'
old=subprocess.check_output(['git','show','81b808f:'+stage],text=True,cwd=root)
new=(root/stage).read_text()
conditions=lambda s:[x for x in s.splitlines() if re.match(r'\| C[1-7] ',x)]
checks['C1_C7_byte_identical_to_entry']=conditions(old)==conditions(new) and len(conditions(new))==7
checks['correct_reviewed_commit']='49253200a801ad9b336789dffb4e93d535a8e903' in new
checks['exact_accepted_endpoint']='+15,325,803.204209456 K' in new
changes=subprocess.check_output(['git','diff','--name-only','81b808f','--','docs/validation'],text=True,cwd=root).splitlines()
checks['historical_validation_evidence_untouched']=sorted(changes)==sorted([stage,'docs/validation/results-first-bnv-pilot/stage-b/REPORT.md'])
source=json.loads((e/'short-invocation.json').read_text())
checks['native_run_production_source_matches_delivered']=all(digest(root/k)==v for k,v in source['source_sha256'].items())
final=json.loads((e/'final-source-sha256.json').read_text())
checks['final_source_manifest']=all(digest(root/k)==v for k,v in final.items())
summary=json.loads((e/'summary.json').read_text())
checks['static_solver_bytes_match_map']=all(digest(root/k)==summary['source_sha256'][k] for k in ['tools/research/bnv_competition/model.py','tools/research/bnv_competition/scan.py'])
checks['copied_input_hashes']=all(digest(d/'inputs'/k)==v['sha256'] for k,v in json.loads((d/'inputs/manifest.json').read_text()).items())
checks['static_independent_pass']=json.loads((e/'verification.json').read_text())['status']=='PASS'
checks['short_independent_pass']=json.loads((e/'short/independent-check.json').read_text())['status']=='PASS'
checks['row_counts']=sum(1 for _ in csv.DictReader((d/'data/map.csv').open()))==5580 and sum(1 for _ in csv.DictReader((d/'data/bnv-only.csv').open()))==5580
bad=[]
for p in d.glob('*.md'):
 for link in re.findall(r'\]\(([^)]+)\)',p.read_text()):
  if not link.startswith(('https:','http:','#')) and not (p.parent/link.split('#')[0]).exists():bad.append([str(p),link])
checks['research_links_exist']=not bad
checks['six_pdf_png_pairs']=all((d/'figures'/(name+ext)).exists() for ext in ['.pdf','.png'] for name in ['A_competition','B_temperature','C_channels','D_applicability','E_reachability','F_scaling_deposition'])
for p in (root/'tools/research/bnv_competition').glob('*.py'):ast.parse(p.read_text(),filename=str(p))
checks['python_syntax']=True
worktrees=subprocess.check_output(['git','worktree','list','--porcelain'],text=True,cwd=root)
checks['nine_worktrees_preserved']=worktrees.count('worktree ')==9
receipt=json.loads((e/'external-receipt.json').read_text())
checks['external_manifest_matches']=digest(Path(receipt['manifest']))==receipt['sha256']==digest(e/'external-SHA256SUMS')
checks['external_files_match']=all(digest(Path(line.split('  ',1)[1]))==line.split('  ',1)[0] for line in (e/'external-SHA256SUMS').read_text().splitlines())
result=dict(status='PASS' if all(checks.values()) else 'FAIL',checks=checks,broken_links=bad,scope='Packaging/source/provenance audit, not a new scientific certificate.')
(e/'packaging-check.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2));assert all(checks.values())
