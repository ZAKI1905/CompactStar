#!/usr/bin/env python3
"""Receipt-bearing, CPU-limited static work. Does not launch stellar executables."""
import hashlib,json,os,pathlib,resource,subprocess,sys,time
S=pathlib.Path(__file__).resolve().parents[3];R=S.parents[1]
O=pathlib.Path(os.environ['BA13_FINAL_OUTPUT_ROOT']);PRE='d5b40b3'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def cpu():
    a=resource.getrusage(resource.RUSAGE_SELF);b=resource.getrusage(resource.RUSAGE_CHILDREN)
    return a.ru_utime+a.ru_stime+b.ru_utime+b.ru_stime
if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Slurm required')
os.sched_setaffinity(0,{min(os.sched_getaffinity(0))});phase=sys.argv[1]
O.mkdir(parents=True,exist_ok=True)
prior=sum(json.loads(p.read_text())['cpu_seconds'] for p in O.glob('receipt-*.json'))
if prior>=570:raise RuntimeError('600 CPU-second ceiling reserve stop')
resource.setrlimit(resource.RLIMIT_CPU,(int(600-prior),int(600-prior)))
decl=S/'docs/validation/ba13-final-decision/PREDECLARATION.md'
assert hashlib.sha256(subprocess.check_output(['git','-C',str(S),'show',PRE+':'+str(decl.relative_to(S))])).hexdigest()==sha(decl)
start=cpu();wall=time.monotonic();records=[];status='ERROR'
def call(label,args,env=None,expected=0):
    remaining=600-prior-(cpu()-start)
    if remaining<15:raise RuntimeError('resource stop')
    cap=min(60,int(remaining-5))
    def limit():resource.setrlimit(resource.RLIMIT_CPU,(cap,cap))
    label=(phase+'-'+label) if phase=='controls-v2' else label
    a=cpu();w=time.monotonic()
    with (O/(label+'.out')).open('x') as f,(O/(label+'.err')).open('x') as e:
        p=subprocess.run(args,stdout=f,stderr=e,env=env,timeout=90,preexec_fn=limit)
    records.append(dict(label=label,command=args,returncode=p.returncode,cpu_seconds=cpu()-a,wall_seconds=time.monotonic()-w,
                        stdout_sha256=sha(O/(label+'.out')),stderr_sha256=sha(O/(label+'.err'))))
    if expected=='nonzero':
        assert p.returncode!=0 and 'FAILED' in (O/(label+'.err')).read_text(),label
    else:assert p.returncode==expected,label
try:
    if phase=='benchmark':
        call('assessment',[sys.executable,'-B',str(S/'tools/linux/ba13_final_decision/assess.py')])
    elif phase in ('validate','controls-v2'):
        benchmark=json.loads((O/'receipt-benchmark.json').read_text())
        assert benchmark['status']=='PASS' and benchmark['cpu_seconds']<30,'benchmark expansion gate'
        for p in (sorted((S/'tests/bnv').glob('test_ba13_*.py')) if phase=='validate' else [S/'tests/bnv/test_ba13_final_decision.py']):
            call(p.stem,[sys.executable,'-B',str(p)])
        review=R/'ba13-closure-independent-review-20261020/scripts/verify_reference.py'
        if phase=='validate':call('independent-reference-replay',[sys.executable,'-I','-B',str(review),str(O/'independent-reference-replay.json')])
        # Mutations operate only on an external copy, with the same unchanged test suite.
        mutations={
          'missing-error':("if v is None or kind=='unresolved' or t.get('missing'):","if v is None: v=0\n        if kind=='unresolved' or t.get('missing'):"),
          'authority-bypass':("if not authority or registry.get(table,{}).get(t.get(field))!=authority:",'if False:'),
          'correlation-cancelled':('return (a - t*d) / denominator','return a*0'),
          'wrong-reference':("if reference_kind != 'X2':",'if False:'),
          'sampled-curvature-promoted':("if evidence_scope not in ('analytic_supremum','accepted_conditional_supremum'):",'if False:'),
          'normalization-omitted':('coefficient=2+fit_norm+(0 if absolute else threshold)','coefficient=2+fit_norm'),
          'dependency-omitted':("if any(set(t.get('dependencies',()))-names for t in terms):return None",'if False:return None')}
        code=(S/'tools/linux/ba13_final_decision/model.py').read_text()
        for label,(old,new) in mutations.items():
            assert code.count(old)==1
            p=O/(phase+'-'+label+'.py');p.write_text(code.replace(old,new))
            env=dict(os.environ,BA13_FINAL_MODEL=str(p))
            call('mutation-'+label,[sys.executable,'-B',str(S/'tests/bnv/test_ba13_final_decision.py')],env,expected='nonzero')
    else:raise ValueError('unknown phase')
    status='PASS'
finally:
    receipt=dict(status=status,phase=phase,slurm_job=os.environ['SLURM_JOB_ID'],source_commit=subprocess.check_output(['git','-C',str(S),'rev-parse','HEAD'],text=True).strip(),
        predeclaration_sha256=sha(decl),cpu_seconds=cpu()-start,wall_seconds=time.monotonic()-wall,prior_cpu_seconds=prior,
        ceiling_cpu_seconds=600,stellar_runs=0,runs=records,
        source_sha256={str(p.relative_to(S)):sha(p) for p in list((S/'tools/linux/ba13_final_decision').glob('*.py'))+[S/'tests/bnv/test_ba13_final_decision.py']})
    with (O/('receipt-'+phase+'.json')).open('x') as f:json.dump(receipt,f,indent=1);f.write('\n')
    print(json.dumps(dict(status=status,cpu_seconds=receipt['cpu_seconds'],phase=phase)))
