#!/usr/bin/env python3
"""Authenticated declarative native pair snapshots. The stage subcommand requires separate owner authorization."""
import argparse,copy,hashlib,json,pathlib,subprocess
from runtime import runtime_binding,verify_runtime
ROOT=pathlib.Path(__file__).resolve().parents[3]
def digest(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def encoded(x):return (json.dumps(x,indent=2,sort_keys=True)+'\n').encode()
def load_pair(directory):
 directory=pathlib.Path(directory)
 sums=json.loads((directory/'identities.json').read_text())
 pair=[]
 for role in ('source','control'):
  p=directory/(role+'.json')
  if digest(p)!=sums[role]:raise ValueError('snapshot identity mismatch: '+role)
  s=json.loads(p.read_text());pair.append(s)
  for item in s['files']:
   if digest(item['path'])!=item['sha256']:raise ValueError('missing/mismatched authenticated dependency: '+item['path'])
 if pair[0]['source']!={'mode':'P2-controlled-uniform-proper-neutron-sink','fractional_Bdot_per_Julian_year':-1e-12}:raise ValueError('source role mismatch')
 if pair[1]['source']!={'mode':'exact-zero','fractional_Bdot_per_Julian_year':0}:raise ValueError('control role mismatch')
 a,b=copy.deepcopy(pair[0]),copy.deepcopy(pair[1]);del a['source'];del b['source']
 if a!=b:raise ValueError('source/control common snapshot mismatch')
 if a['spin_rad_s']!=0 or a['B0_count']!=7.616906518718875e56 or a['Julian_year_s']!=31557600 or a['initial_time_s']!=0:raise ValueError('native physical configuration mismatch')
 if a['thermal']!={'Tmin_for_derivative_MeV':0,'clamp_to_domain':False}:raise ValueError('thermal configuration mismatch')
 if a['integration']!={'method':'uninterrupted breakpoint-aware RKF45','initial_step_s':1,'relative_tiers':[1e-7,1e-9,1e-11],'absolute_tiers':[[1e-12,1e-18,1e-18],[1e-14,1e-20,1e-20],[1e-16,1e-22,1e-22]],'checkpoint':'passive; existing isolated two-level self-qualified rk8pd'}:raise ValueError('integration configuration mismatch')
 bound={x['path']:x['sha256'] for x in a['files']}
 for key in ('archive','verify_binary','stage_binary','entry','frozen','coefficients'):
  if a[key] not in bound:raise ValueError('unbound runtime input: '+key)
 if a['archive_sha256']!=bound[a['archive']]:raise ValueError('archive identity mismatch')
 if a['schema']!='RF-BNV-P2-native-owned-snapshot-v1' or a['initial_state']!=[0,0,0] or a['Tref_K']!=1e8:raise ValueError('initial-state/schema mismatch')
 if set(a.get('runtime_bindings',{}))!={a['verify_binary'],a['stage_binary']}:raise ValueError('missing native runtime bindings')
 verify_runtime(a['runtime_bindings'],bound)
 return pair

def publish(run,out):
 run=pathlib.Path(run).resolve();out=pathlib.Path(out);out.mkdir(parents=True,exist_ok=False)
 if 'PASS native pair static checks' not in (run/'verify.log').read_text():raise ValueError('native static receipt missing')
 if '100% tests passed' not in (run/'contracts.log').read_text():raise ValueError('retained contracts receipt missing')
 # EOS path spelling is part of the retained Phase-5B identity; do not canonicalize it.
 native=ROOT.parents[1]/'toolchain/../ba13-recovery/run1'
 archive=run/'capture/sequence.archive'
 inventory=json.loads((ROOT/'docs/validation/results-first-bnv-pilot/native-assembly/evidence/input-inventory.json').read_text())
 paths=[pathlib.Path(x['path']) for x in inventory['checks']]
 paths += [native/p.relative_to(ROOT.parents[1]/'ba13-recovery/run1') for p in list(paths)]
 paths += [archive,run/'build/tests/native_p2_verify',run/'build/tests/native_p2_stage']
 paths += [ROOT.parents[1]/'toolchain/dependency-lock-linux.json']
 bindings={str(run/'build/tests'/binary):runtime_binding(run/'build/tests'/binary) for binary in ('native_p2_verify','native_p2_stage')}
 for binding in bindings.values():
  paths += [pathlib.Path(item['resolved']) for item in binding['dependencies']]
 paths += [ROOT/'docs/validation/results-first-bnv-pilot/preparation-registry.json']
 side=ROOT/'docs/validation/results-first-bnv-pilot/endpoint-acceptance'
 for line in (side/'SHA256SUMS').read_text().splitlines():
  expected,name=line.split(maxsplit=1);p=ROOT/name
  if digest(p)!=expected:raise ValueError('conditional sidecar dependency mismatch: '+name)
  paths.append(p)
 summary=json.loads((side/'evidence/assessment/summary.json').read_text())
 if summary['status']!='CONDITIONAL_STATIC_PASS' or float(summary['maximum_utilization'])>1 or summary['historical_certificate_promoted']:raise ValueError('conditional sidecar applicability mismatch')
 paths += [side/'evidence/assessment/summary.json',side/'evidence/assessment/ledger.json.gz',side/'acceptance-registry.json',side/'OWNER_APPROVAL.md',side/'inputs.json']
 paths += [ROOT/p for p in subprocess.check_output(['git','ls-files','CompactStar','tests/native','tests/bnv','tests/rotochemical','tests/analysis','tools/linux/native_pair'],cwd=ROOT,text=True).splitlines() if (ROOT/p).is_file()]
 files=[{'path':str(p),'sha256':digest(p)} for p in sorted(set(paths))]
 s={'schema':'RF-BNV-P2-native-owned-snapshot-v1','experiment':'RF-BNV-P2-20261108','candidate':True,
 'assembly_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
 'files':files,'runtime_bindings':bindings,'archive':str(archive),'archive_sha256':digest(archive),'central_contributor':2,
 'initial_time_s':0,'initial_state':[0,0,0],'state_units':['ln(Tinf/Tref)','MeV','MeV'],'Tref_K':1e8,
 'spin_rad_s':0,'B0_count':7.616906518718875e56,'Julian_year_s':31557600,
 'photon_boundary':{'surface_model':'EnvelopeTbTs','composition':'Iron','rho_b_g_cm3':1e10,'xi':0,'radiating_fraction':1,'global_scale':1,
 'law':'Potekhin2003 Iron dispatch to 1997 fit; Ts^4=1e24*g14*(1.81*Tb/1e8)^2.42; Lgamma_inf=4*pi*(R_km*1e5)^2*sigma*Ts^4*exp(2*nu_surface)'},
 'thermal':{'Tmin_for_derivative_MeV':0,'clamp_to_domain':False},
 'integration':{'method':'uninterrupted breakpoint-aware RKF45','initial_step_s':1,'relative_tiers':[1e-7,1e-9,1e-11],
 'absolute_tiers':[[1e-12,1e-18,1e-18],[1e-14,1e-20,1e-20],[1e-16,1e-22,1e-22]],
 'checkpoint':'passive; existing isolated two-level self-qualified rk8pd'},
 'certificate_status':'conditional static applicability only; historical runtime certificate unchanged; nonlinear actual path unqualified',
 'profile':str(native/'raw/chemical-profiles/target-0'),'thermal_path':str(native/'raw/chemical-profiles/target-0/thermal'),
 'entry':str(ROOT/'docs/validation/results-first-bnv-pilot/preparation-registry.json'),
 'frozen':str(native/'fit/certificate.tsv'),'coefficients':str(native/'raw/coefficients.tsv'),
 'verify_binary':str(run/'build/tests/native_p2_verify'),'stage_binary':str(run/'build/tests/native_p2_stage')}
 ids={}
 for role,source in [('source',{'mode':'P2-controlled-uniform-proper-neutron-sink','fractional_Bdot_per_Julian_year':-1e-12}),('control',{'mode':'exact-zero','fractional_Bdot_per_Julian_year':0})]:
  s['source']=source;p=out/(role+'.json');p.write_bytes(encoded(s));ids[role]=digest(p)
 (out/'identities.json').write_bytes(encoded(ids));load_pair(out);print(json.dumps(ids,indent=2))

def main():
 p=argparse.ArgumentParser();sub=p.add_subparsers(dest='mode',required=True)
 q=sub.add_parser('publish');q.add_argument('run');q.add_argument('directory')
 q=sub.add_parser('check');q.add_argument('directory')
 q=sub.add_parser('stage');q.add_argument('directory');q.add_argument('years',type=int,choices=[1,2048]);q.add_argument('fresh_work')
 q=sub.add_parser('verify');q.add_argument('directory');q.add_argument('fresh_work')
 a=p.parse_args()
 if a.mode=='publish':publish(a.run,a.directory);return
 pair=load_pair(a.directory)
 if a.mode=='check':print('PASS authenticated native matched snapshots');return
 s=pair[0];w=pathlib.Path(a.fresh_work);w.mkdir(parents=True,exist_ok=False)
 binary=s['stage_binary'] if a.mode=='stage' else s['verify_binary']
 args=[binary,s['profile'],s['thermal_path'],s['entry'],s['frozen'],s['coefficients'],s['archive'],s['archive_sha256'],str(w/'assembly'),*([str(a.years),str(w/'trajectory')] if a.mode=='stage' else [str(w/'static-rhs.tsv')])]
 subprocess.run(args,cwd=ROOT,check=True)
if __name__=='__main__':main()
