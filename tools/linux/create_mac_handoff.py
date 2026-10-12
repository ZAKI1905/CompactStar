#!/usr/bin/env python3
"""Source-only independent AppleClang arithmetic receipt bundle."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import zipfile
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[2]
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
a.output.mkdir(parents=True,exist_ok=False)
files=['CompactStar/EOS/src/LocalThermodynamics.cpp','CompactStar/EOS/src/TrackRFreeGasThermodynamics.cpp',
       'tests/relativity/relativistic_unit_boundary.cpp','tests/eos/rotochemical_trackr_pe.cpp',
       'tests/linux/arithmetic_probe.cpp','tests/linux/arithmetic_oracle.py',
       'docs/validation/linux-corrections/WIDTH_AUDIT.md','tests/eos/structure1/table.hpp']
for name in files:
 target=a.output/'sources'/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes((ROOT/name).read_bytes())
spec=importlib.util.spec_from_file_location('oracle',ROOT/'tests/linux/arithmetic_oracle.py');oracle=importlib.util.module_from_spec(spec);spec.loader.exec_module(oracle)
constants=[939.56542052,938.27208816,0.51099895000,5.067730716156395e-3]
expected={'constants_hex':[x.hex() for x in constants],'decimal80':str(oracle.calculate(constants,80)),
          'decimal120':str(oracle.calculate(constants,120)),'rounded_hex':float(oracle.calculate(constants,120)).hex(),
          'derivation':'Exact binary constants, Kallen momentum squared times sqrt, k^3/(3*pi^2); Machin pi with guarded alternating-series remainder, 80 and 120 digits.'}
(a.output/'expected.json').write_text(json.dumps(expected,indent=2)+'\n')
# Standalone probe faithfully transcribes the changed arithmetic, not a linked
# CompactStar test. Receipt distinguishes these scopes explicitly.
(a.output/'portable_probe.cpp').write_text(r'''#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <stdexcept>
int main(){
const double mn=939.56542052,mp=938.27208816,me=0.51099895000,k=5.067730716156395e-3;
auto density=[&](double mu){const double p=std::sqrt((mu-me)*(mu+me));const double pi=std::acos(-1.0);const double ratio=p/(1.0/k);return ratio*ratio*ratio/(3.0*pi*pi);};
const double mu=std::fma(mn-mp,mn+mp,me*me)/(2.0*mn);
std::cout<<std::hexfloat<<"constants\t"<<mn<<'\t'<<mp<<'\t'<<me<<'\t'<<k<<'\n';
std::cout<<"onset\t"<<density(mu)<<'\n';
for(int offset:{-32,-8,-1,0,1,8,32}){double x=mu;for(int j=0;j<std::abs(offset);++j)x=std::nextafter(x,offset<0?0.0:2.0);std::cout<<"inverse\t"<<x<<'\t'<<density(x)<<'\n';}
for(double x:{0.8125,1.0625,1.75,2.125,8.5})std::cout<<"inverse\t"<<x<<'\t'<<density(x)<<'\n';
constexpr double c=29979245800.,G=6.673e-8;
const double boundary=1e10*(G/(c*c)*1e10),below=std::nextafter(boundary,0.0),above=std::nextafter(boundary,INFINITY);
for(double value:{below,boundary,above}){double col[]={above,value,below};int index=-1;for(int i=2;i>=0;--i)if(col[i]>=boundary){index=i;break;}if(index!=(value<boundary?0:1))throw std::runtime_error("U11 rule");}
std::cout<<"U11_probe\tPASS\t"<<below<<'\t'<<boundary<<'\t'<<above<<'\n';
std::cout<<"types\t"<<sizeof(float)<<'\t'<<std::numeric_limits<float>::digits<<'\t'<<sizeof(double)<<'\t'<<std::numeric_limits<double>::digits<<'\t'<<sizeof(long double)<<'\t'<<std::numeric_limits<long double>::digits<<'\n';
}
''')
(a.output/'run_mac.py').write_text(r'''#!/usr/bin/env python3
import hashlib,json,platform,subprocess
from pathlib import Path
root=Path(__file__).resolve().parent
if platform.system()!='Darwin':raise SystemExit('Refusing Mac receipt: this is not actual Darwin execution')
compiler=subprocess.check_output(['xcrun','--find','clang++'],text=True).strip()
version=subprocess.check_output([compiler,'--version'],text=True)
if 'Apple clang' not in version:raise SystemExit('Actual AppleClang required')
r={'actual_mac_execution':True,'scope':'standalone arithmetic/U11 comparison probes; linked production tests still required','platform':platform.platform(),'compiler':version,'runs':[],'linked_production_tests':'PENDING'}
for mode,flags in [('default',[]),('contraction-off',['-ffp-contract=off'])]:
 exe=root/('probe-'+mode);cmd=[compiler,'-std=c++17','-O0',*flags,str(root/'portable_probe.cpp'),'-o',str(exe)]
 subprocess.run(cmd,check=True)
 raw=subprocess.check_output([str(exe)],text=True);out=root/(mode+'.tsv');out.write_text(raw)
 subprocess.run(['python3',str(root/'sources/tests/linux/arithmetic_oracle.py'),'--standalone','--input',str(out),'--output',str(root/(mode+'-oracle.json'))],check=True)
 r['runs'].append({'mode':mode,'command':cmd,'raw_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'libraries':subprocess.check_output(['otool','-L',str(exe)],text=True)})
(root/'mac-receipt.json').write_text(json.dumps(r,indent=2)+'\n')
''')
(a.output/'README.txt').write_text('''Actual Mac confirmation; no cluster result is Mac execution.
1. Verify: shasum -a 256 -c SHA256SUMS
2. Run: python3 run_mac.py
   Requires native Darwin, xcrun AppleClang and Python 3 (stdlib only).
   Runs default contraction and explicit off standalone arithmetic probes;
   changes no qualified Mac flags. Expected values and derivation in expected.json
   and sources/tests/linux/arithmetic_oracle.py.
3. Production confirmation in a disposable Mac checkout: apply arithmetic.patch,
   copy the supplied tests/linux sources and append the supplied CMake fragment.
   Configure with the owner's authenticated Mac dependency stack and unchanged
   qualified flags. MAC_BUILD must name that independently configured build.
   cmake --build "$MAC_BUILD" --target rotochemical_trackr_pe relativistic_unit_boundary linux_arithmetic_probe
   ctest --test-dir "$MAC_BUILD" -R '^(rotochemical_trackr_pe|relativistic_unit_boundary|candidate_arithmetic_oracle)$' --output-on-failure
   Also rerun chemical_production_contract, phase5d_response and
   phase5d_independent_oracles for retained width-dependent quantities.
   Record actual repository commit, dependency lock hash, compile_commands hash,
   compiler/OS/CPU, flags and raw logs in production-results-template.json.
   Repeat production confirmation for both Debug and Release, keeping a
   separately named completed results receipt and raw logs for each mode.
4. Return mac-receipt.json, filled production-results-template.json and raw logs
   through the office transfer workflow. Do not infer production qualification
   from a standalone probe or a simulated double-width run.
These files are candidate changes, not baseline promotion or owner acceptance.
''')
paths=files[:4]+['tests/eos/structure1/table.hpp']
patch=subprocess.check_output(['git','diff','f4e7c18','--',*paths],cwd=ROOT)
(a.output/'arithmetic.patch').write_bytes(patch)
s=(ROOT/'tests/CMakeLists.txt').read_text();(a.output/'CMake-fragment.txt').write_text(s[s.index('# ADR-0020 candidate independent'):])
(a.output/'production-results-template.json').write_text(json.dumps(dict(status='NOT_RUN',actual_mac_execution=False,source_commit=None,source_hashes={},dependency_lock_sha256=None,compile_commands_sha256=None,compiler=None,cpu=None,fp_flags=None,tests={},raw_logs=[]),indent=2)+'\n')
identities=dict(compactstar_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
    starting_commit='f4e7c18f568045febc5bf65c6e95200752196dad',Zaki='e263a6e180c5c417198e7778bd21fc9c0a32dc33',CONFIND='b0cbd510fd3fd0c772fa50499cd749287cb39e7b',
    source_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files})
(a.output/'identities.json').write_text(json.dumps(identities,indent=2)+'\n')
(a.output/'SHA256SUMS').write_text(''.join(hashlib.sha256(f.read_bytes()).hexdigest()+'  '+f.relative_to(a.output).as_posix()+'\n' for f in sorted(a.output.rglob('*')) if f.is_file()))
zip_path=a.output.with_suffix('.zip')
with zipfile.ZipFile(zip_path,'x',zipfile.ZIP_DEFLATED) as z:
 for f in sorted(a.output.rglob('*')):
  if f.is_file():z.write(f,a.output.name+'/'+f.relative_to(a.output).as_posix())
zip_path.with_suffix('.zip.sha256').write_text(hashlib.sha256(zip_path.read_bytes()).hexdigest()+'  '+zip_path.name+'\n')
print(zip_path)
