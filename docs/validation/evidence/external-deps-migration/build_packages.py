#!/usr/bin/env python3
import subprocess, os, json, hashlib, pathlib, time
root=pathlib.Path(__file__).resolve().parent
env=dict(os.environ);env.pop('CONDA_PREFIX',None);env['PYTHONDONTWRITEBYTECODE']='1'
cmake='/opt/homebrew/bin/cmake'
records=[]
for mode in ('Debug','Release'):
 for name,src,sha in [('Zaki','/Users/keeper/Documents/CompactStar/worktrees/ZakiLib-2.0.1-fp-preservation','e263a6e180c5c417198e7778bd21fc9c0a32dc33'),('CONFIND','/Users/keeper/Documents/CompactStar/external/CONFIND','b0cbd510fd3fd0c772fa50499cd749287cb39e7b')]:
  assert subprocess.check_output(['git','-C',src,'rev-parse','HEAD'],text=True).strip()==sha
  assert not subprocess.check_output(['git','-C',src,'status','--porcelain'],text=True).strip()
  build=root/f'build-{name}-{mode}';prefix=root/mode/name
  args=[cmake,'-S',src,'-B',str(build),f'-DCMAKE_BUILD_TYPE={mode}','-DCMAKE_C_COMPILER=/usr/bin/clang','-DCMAKE_CXX_COMPILER=/usr/bin/clang++','-DGSL_ROOT_DIR=/opt/local','-DGSL_CONFIG_EXECUTABLE=/opt/local/bin/gsl-config','-DGSL_INCLUDE_DIR=/opt/local/include','-DGSL_LIBRARY=/opt/local/lib/libgsl.dylib','-DGSL_CBLAS_LIBRARY=/opt/local/lib/libgslcblas.dylib',f'-DCMAKE_INSTALL_PREFIX={prefix}','-DCMAKE_EXPORT_COMPILE_COMMANDS=ON']
  if name=='Zaki':args+=['-DZAKI_BUILD_TESTS=ON','-DZAKI_BUILD_EXAMPLES=OFF']
  else:args+=['-DCONFIND_BUILD_TESTS=ON',f'-DZaki_DIR={root/mode/"Zaki/lib/cmake/Zaki"}']
  for stage,cmd in [('configure',args),('build',[cmake,'--build',str(build),'-j','6']),('test',['/opt/homebrew/bin/ctest','--test-dir',str(build),'--output-on-failure','-j','1']),('install',[cmake,'--install',str(build)])]:
   with (root/'evidence'/f'{name}-{mode}-{stage}.log').open('w') as log:
    r=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,env=env)
   print(name,mode,stage,r.returncode,flush=True)
   if r.returncode:raise SystemExit(r.returncode)
  def hashfile(p):return hashlib.sha256(p.read_bytes()).hexdigest()
  def manifest(folder):
   return ''.join(hashfile(f)+'  '+str(f.relative_to(prefix))+'\n' for f in sorted((prefix/folder).rglob('*')) if f.is_file())
  headers=manifest('include');configs=manifest('lib/cmake')
  (prefix/'headers.sha256').write_text(headers);(prefix/'configs.sha256').write_text(configs)
  record=dict(name=name,mode=mode,source=src,source_sha=sha,dirty=False,prefix=str(prefix),archive_sha256=hashfile(prefix/'lib'/f'lib{name}.a'),headers_sha256=hashlib.sha256(headers.encode()).hexdigest(),configs_sha256=hashlib.sha256(configs.encode()).hexdigest(),compiler=subprocess.check_output(['/usr/bin/clang++','--version'],text=True),cmake=subprocess.check_output([cmake,'--version'],text=True),cache_sha256=hashfile(build/'CMakeCache.txt'))
  (prefix/'provenance.json').write_text(json.dumps(record,indent=2)+'\n');records.append(record)
  (root/'evidence/packages.json').write_text(json.dumps(records,indent=2)+'\n')
