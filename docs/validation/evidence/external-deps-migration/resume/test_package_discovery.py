from pathlib import Path
import subprocess,json,shutil,tempfile
q=Path(__file__).resolve().parent;src=Path('/Users/keeper/Documents/CompactStar/worktrees/CompactStar-external-deps-taskmanager-migration');work=Path(tempfile.mkdtemp(prefix='csm-g0-',dir='/private/tmp'));records=[]
base=['/opt/homebrew/bin/cmake','-S',str(src),'-DBUILD_TESTING=OFF','-DCMAKE_BUILD_TYPE=Debug','-DPython3_EXECUTABLE=/Users/keeper/miniforge3/bin/python3','-DGSL_ROOT_DIR=/opt/local','-DGSL_CONFIG_EXECUTABLE=/opt/local/bin/gsl-config','-DGSL_INCLUDE_DIR=/opt/local/include','-DGSL_LIBRARY=/opt/local/lib/libgsl.dylib','-DGSL_CBLAS_LIBRARY=/opt/local/lib/libgslcblas.dylib']
z=str(q/'Debug/Zaki');c=str(q/'Debug/CONFIND')
def run(name,args,success):
 command=base+['-B',str(work/name)]+args
 p=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
 (q/'evidence'/f'G0-{name}.log').write_text(p.stdout)
 records.append(dict(case=name,expected_success=success,returncode=p.returncode,passed=(p.returncode==0)==success,command=command))
 assert records[-1]['passed'],name
run('no-prefixes',['-DZaki_DIR='+z+'/lib/cmake/Zaki','-DCONFIND_DIR='+c+'/lib/cmake/CONFIND'],False)
run('wrong-mode',['-DCOMPACTSTAR_ZAKI_PREFIX='+str(q/'Release/Zaki'),'-DCOMPACTSTAR_CONFIND_PREFIX='+c],False)
run('explicit-pinned-over-stale-cache',['-DCOMPACTSTAR_ZAKI_PREFIX='+z,'-DCOMPACTSTAR_CONFIND_PREFIX='+c,'-DZaki_DIR=/usr/local/lib/cmake/Zaki','-DCONFIND_DIR=/opt/homebrew/lib/cmake/CONFIND'],True)
copy=work/'modified-Zaki';shutil.copytree(z,copy)
header=copy/'include/Zaki/Version.hpp';header.write_text(header.read_text()+'\n// modified\n')
run('changed-header',['-DCOMPACTSTAR_ZAKI_PREFIX='+str(copy),'-DCOMPACTSTAR_CONFIND_PREFIX='+c],False)
shutil.copyfile(Path(z)/'include/Zaki/Version.hpp',header)
archive=copy/'lib/libZaki.a'
with archive.open('ab') as f:f.write(b'changed')
run('changed-archive',['-DCOMPACTSTAR_ZAKI_PREFIX='+str(copy),'-DCOMPACTSTAR_CONFIND_PREFIX='+c],False)
shutil.copyfile(Path(z)/'lib/libZaki.a',archive)
prov=copy/'provenance.json';data=json.loads(prov.read_text());data['source_sha']='c8c68131b04e9d216673725075bef38df81e6041';prov.write_text(json.dumps(data))
run('wrong-source',['-DCOMPACTSTAR_ZAKI_PREFIX='+str(copy),'-DCOMPACTSTAR_CONFIND_PREFIX='+c],False)
shutil.copyfile(Path(z)/'provenance.json',prov)
injected=work/'injected-target.cmake'
injected.write_text('add_library(Zaki::Zaki STATIC IMPORTED GLOBAL)\nset_target_properties(Zaki::Zaki PROPERTIES IMPORTED_LOCATION_DEBUG /usr/local/lib/libZaki.a INTERFACE_INCLUDE_DIRECTORIES /usr/local/include)\n')
run('preexisting-wrong-target',['-DCOMPACTSTAR_ZAKI_PREFIX='+z,'-DCOMPACTSTAR_CONFIND_PREFIX='+c,'-DCMAKE_PROJECT_INCLUDE='+str(injected)],False)
(q/'evidence/package-discovery-controls.json').write_text(json.dumps(records,indent=2)+'\n');print(len(records),'controls passed')
