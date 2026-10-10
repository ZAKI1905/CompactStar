#!/usr/bin/env python3
"""Copy available scientific inputs only after independent historical-hash checks."""
import argparse,hashlib,json,pathlib,re,shutil
from recover_mac import tree,sha
R=pathlib.Path(__file__).resolve().parents[3]
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--workspace',type=pathlib.Path,required=True);ap.add_argument('--output',type=pathlib.Path,required=True);a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=False)
 auth=json.loads((R/'docs/validation/evidence/obs32-forensics/authentication.json').read_text());records=[]
 for key,source in auth['input_paths'].items():
  p=pathlib.Path(source);expected=auth['input_hashes'][key];assert tree(p)==expected
  relative=pathlib.Path('phase6')/(key if p.is_dir() else key+p.suffix);dest=a.output/relative;dest.parent.mkdir(parents=True,exist_ok=True)
  if p.is_dir():shutil.copytree(p,dest)
  else:shutil.copyfile(p,dest)
  assert tree(dest)==expected
  records.append({'key':key,'original':source,'relative_path':relative.as_posix(),'sha256':expected,'authority':'docs/validation/evidence/obs32-forensics/authentication.json','status':'PASS'})
 text=(R/'docs/validation/PASSIVE_COOLING_BASELINE.md').read_text();cold=text.split('**Raw CompOSE cold**')[1].split('**Raw CompOSE finite-`T`**')[0];hot=text.split('**Raw CompOSE finite-`T`**')[1].split('Supplied via')[0]
 for dirname,part in [('DS-CMF-1-with-crust',cold),('DNS-CMF-Hadronic-with-electrons',hot)]:
  expected=dict((n,h) for h,n in re.findall(r'\| `([0-9a-f]{64})` \| `([^`]+)` \|',part))
  if dirname=='DS-CMF-1-with-crust':
   expected.update({'DS(CMF)-1_with_crust.eos':'5747dd73256c0c28bc56be337cbb96d0918a54bc9ed9fc40984c5befd47ae5dd','eos.mr':'98e38994b4191b1024f5f5edbada0e6c1abe9e0cb6d81de94e47be08e3e8efa8'})
  for n,h in expected.items():
   src=a.workspace/'data/compose'/dirname/n;assert sha(src)==h,str(src)
   relative=pathlib.Path('compose')/dirname/n;dest=a.output/relative;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dest);assert sha(dest)==h
   records.append({'key':'compose/'+dirname+'/'+n,'original':str(src),'relative_path':relative.as_posix(),'sha256':h,'authority':'docs/validation/PASSIVE_COOLING_BASELINE.md:268-300 and TOV_REFERENCE.md:124-136','status':'PASS'})
 review=a.workspace/'data/phase6-breakpoint-candidate/obs32-review-evidence'
 for line in (review/'SHA256SUMS').read_text().splitlines():
  h,n=line.split(maxsplit=1);n=n.lstrip('*');src=review/n;assert sha(src)==h,n
  if n.startswith('inputs/'):
   relative=pathlib.Path('obs32')/n;dest=a.output/relative;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dest);assert sha(dest)==h
   records.append({'key':'obs32/'+n,'original':str(src),'relative_path':relative.as_posix(),'sha256':h,'authority':'Original obs32-review-evidence/SHA256SUMS (entire original manifest verified before input subset copy)','status':'PASS'})
 (a.output/'input-map.json').write_text(json.dumps(records,indent=2)+'\n')
 print('PASS copied',len(records),'input identities; bytes',sum(p.stat().st_size for p in a.output.rglob('*') if p.is_file()))
if __name__=='__main__':main()
