#!/usr/bin/env python3
"""Assemble a deterministic source/input handoff from clean committed sources."""
import argparse,hashlib,json,pathlib,shutil,subprocess,zipfile
HERE=pathlib.Path(__file__).resolve().parent
ROOT=HERE.parents[2]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--workspace',type=pathlib.Path,required=True);p.add_argument('--staging',type=pathlib.Path,required=True);p.add_argument('--output',type=pathlib.Path,required=True);a=p.parse_args()
 assert subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True)=='','source must be clean and committed'
 assert subprocess.check_output(['git','branch','--show-current'],cwd=ROOT,text=True).strip()=='docs/phase6-cluster-preparation'
 subprocess.run(['git','diff','--exit-code','5433fc2','--','CompactStar','cmake','CMakeLists.txt','tests','dependencies'],cwd=ROOT,check=True)
 a.output.mkdir(parents=True,exist_ok=False);payload=a.output/'CompactStar-EKU-handoff-20261010';payload.mkdir()
 subprocess.run(['python3',str(HERE/'package_sources.py'),'--workspace',str(a.workspace),'--compactstar',str(ROOT),'--output',str(payload/'bundles')],check=True)
 for directory in ['scientific','upstream']:shutil.copytree(a.staging/directory,payload/directory)
 shutil.copyfile(a.staging/'unavailable-required-inputs.json',payload/'unavailable-required-inputs.json')
 for name in ['README.md','import_sources.py','verify_transfer.py']:shutil.copyfile(HERE/name,payload/name)
 shutil.copyfile(ROOT/'docs/validation/PHASE6_CLUSTER_TRANSFER_PREPARATION.md',payload/'PREPARATION_REPORT.md')
 shutil.copytree(HERE/'evidence',payload/'evidence')
 identities={name:json.loads((payload/'bundles'/(name+'.json')).read_text())['refs'] for name in ['CompactStar','Zaki','CONFIND']}
 (payload/'source-identities.json').write_text(json.dumps({'source_refs':identities,'mac_dependency_lock_sha256':sha(ROOT/'cmake/dependency-lock.json'),'source_production_diff_from_readiness':'EMPTY','prepared_on':'2026-10-10','status':'SOURCE HANDOFF; LINUX UNQUALIFIED; NO CAMPAIGN AUTHORIZATION'},indent=2)+'\n')
 files=sorted(p for p in payload.rglob('*') if p.is_file());assert not any(p.is_symlink() for p in payload.rglob('*'))
 (payload/'SHA256SUMS').write_text(''.join(sha(p)+'  '+p.relative_to(payload).as_posix()+'\n' for p in files))
 subprocess.run(['python3',str(payload/'verify_transfer.py'),'--root',str(payload)],check=True)
 archive=a.output/(payload.name+'.zip')
 with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=9,allowZip64=True) as z:
  for p in sorted(p for p in payload.rglob('*') if p.is_file()):
   info=zipfile.ZipInfo(p.relative_to(a.output).as_posix(),date_time=(2026,10,10,0,0,0));info.create_system=3;info.external_attr=(0o100644<<16);info.compress_type=zipfile.ZIP_DEFLATED
   with z.open(info,'w',force_zip64=True) as dst,p.open('rb') as src:shutil.copyfileobj(src,dst,1048576)
 h=sha(archive);(a.output/(archive.name+'.sha256')).write_text(h+'  '+archive.name+'\n')
 receipt={'zip':str(archive),'sha256':h,'bytes':archive.stat().st_size,'payload_sha256sums_sha256':sha(payload/'SHA256SUMS'),'source_refs':identities,'source_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()}
 (a.output/'transfer-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
