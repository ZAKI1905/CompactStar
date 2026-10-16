#!/usr/bin/env python3
"""Post-trace, pre-result support/endpoint/shell-coordinate mechanism test."""
import hashlib,json,os,pathlib,resource,subprocess,sys,time
SRC=pathlib.Path(__file__).resolve().parents[3];ROOT=SRC.parents[1];OUT=ROOT/'x2-fixed-step';DOC=SRC/'docs/validation/x2-fixed-step';TC=ROOT/'toolchain'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
if len(sys.argv)>1 and sys.argv[1]=='build':
 base=(OUT/'build/x1-regularized.cpp').read_text();files=[pathlib.Path(__file__),DOC/'SHELL_DIAGNOSTIC_PREDECLARATION.md',OUT/'build/x1-regularized.cpp']
 for variant in ['support','endpoint','shell-s']:
  s=base.replace('bool transformed=false;','bool transformed=false;bool neutron_segment=true;').replace('transformed=reg;onset=on;','transformed=reg;onset=on;neutron_segment=high>hn;').replace('f[3]=H>o.hn ?', 'f[3]=o.neutron_segment ?')
  if variant!='support':s=s.replace('auto q=o.at_H(H); double u=', 'auto q=o.at_H(H);if(!o.neutron_segment && q.ns[0]>0){q.ns[0]=0;q=o.eos.from_ns(q.ns);} double u=')
  if variant=='shell-s':s=s.replace('bool reg=low>=hn;double on=(low==hm)?hm:hn;', 'bool reg=true;double on=(low==hm)?hm:((low==hs)?0:hn);')
  path=OUT/f'build/shell-{variant}.cpp';path.write_text(s);exe=OUT/f'build/shell-{variant}';files += [path,exe]
  subprocess.run([str(TC/'bin/g++'),'-O3','-std=c++17','-ffp-contract=off','-isystem',str(TC/'install/Release/Zaki/include'),'-isystem',str(TC/'install/gsl-2.7.1/include'),str(path),'-o',str(exe),str(TC/'install/Release/Zaki/lib/libZaki.a'),str(TC/'install/gsl-2.7.1/lib/libgsl.a'),str(TC/'install/gsl-2.7.1/lib/libgslcblas.a')],check=True)
 (DOC/'shell-gate.json').write_text(json.dumps({'files':{str(p):sha(p) for p in files}},indent=2)+'\n');sys.exit()
os.sched_setaffinity(0,{min(os.sched_getaffinity(0))});assert os.getenv('SLURM_JOB_ID')
for p,h in json.loads((DOC/'shell-gate.json').read_text())['files'].items():assert sha(pathlib.Path(p))==h
q=json.loads((SRC/'docs/validation/x1-oracle-anchor/targets.json').read_text())['0']
x1=json.loads((SRC/'docs/validation/x1-oracle-anchor/evidence/runs.json').read_text());hc=str(next(r['values']['Hc'] for r in x1 if r['target']=='0' and r['variant']=='regularized' and r['tol']=='2e-13'))
def cpu():return sum(r.ru_utime+r.ru_stime for r in [resource.getrusage(resource.RUSAGE_SELF),resource.getrusage(resource.RUSAGE_CHILDREN)])
prior=sum(json.loads(p.read_text())['cpu_seconds'] for p in OUT.glob('receipt-*.json'));assert prior<750
start=cpu();runs=[]
try:
 for variant in ['support','endpoint','shell-s']:
  for tol in ['2e-10','2e-11','2e-12','2e-13']:
   for fixed in [False,True]:
    label=f'shell-{variant}-'+('fixed' if fixed else 'shoot')+'-'+tol;out=OUT/(label+'.out');err=OUT/(label+'.err');assert not out.exists()
    cmd=[str(OUT/f'build/shell-{variant}'),q['B'],hc if fixed else q['rho_guess'],q['surfaceP_cgs'],tol];t=cpu()
    p=subprocess.run(cmd,capture_output=True,text=True,env=dict(os.environ,**({'X2_FIXED_H':'1'} if fixed else {})),timeout=30);out.write_text(p.stdout);err.write_text(p.stderr)
    runs.append(dict(label=label,command=cmd,fixed_H=fixed,cpu_seconds=cpu()-t,returncode=p.returncode,stdout_sha256=sha(out),stderr_sha256=sha(err)))
    if p.returncode:raise RuntimeError(label+': '+p.stderr)
finally:
 (OUT/'receipt-shell.json').write_text(json.dumps(dict(phase='shell',job=os.getenv('SLURM_JOB_ID'),affinity=sorted(os.sched_getaffinity(0)),cpu_seconds=cpu()-start,prior_cpu_seconds=prior,runs=runs,gate_sha256=sha(DOC/'shell-gate.json')),indent=2)+'\n')
