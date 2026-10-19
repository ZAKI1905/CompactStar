#!/usr/bin/env python3
"""Independent arithmetic/provenance checks; never imports the closure assessor."""
import csv,decimal,gzip,hashlib,json,os,pathlib,resource,subprocess,sys,time
import numpy as np
S=pathlib.Path(__file__).resolve().parents[3];R=S.parents[1]
O=pathlib.Path(os.environ.get('BA13_CLOSURE_OUTPUT_ROOT',R/'ba13-closure-20261019'));F=O/'assessment-final';REV=R/'ba13-headroom-independent-review-20261018'
D=decimal.Decimal

def read(p,precise=False):return json.loads(p.read_text(),parse_float=D if precise else float)
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def main():
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Slurm required')
    os.sched_setaffinity(0,{min(os.sched_getaffinity(0))})
    start=time.process_time();wall=time.monotonic();checks={};tests=[]
    children_start=resource.getrusage(resource.RUSAGE_CHILDREN)
    def check(name,ok):
        checks[name]=bool(ok)
        if not ok:raise AssertionError(name)
    for script in sorted((S/'tests/bnv').glob('test_ba13_*.py')):
        p=subprocess.run([sys.executable,'-B',str(script)],capture_output=True,text=True)
        (O/(script.stem+'.out')).write_text(p.stdout+p.stderr)
        tests.append(dict(script=str(script.relative_to(S)),returncode=p.returncode,
                          output_sha256=sha(O/(script.stem+'.out'))))
        check(script.stem,p.returncode==0)
    for line in (REV/'SHA256SUMS').read_text().splitlines():
        digest,rel=line.split(maxsplit=1)
        check('review '+rel,sha(REV/rel)==digest)
    for rel,root in [('docs/validation/ba13-headroom/evidence/external-manifest.json',R/'ba13-headroom'),
                     ('docs/validation/x2-fixed-step/evidence/external-manifest.json',R/'x2-fixed-step'),
                     ('docs/validation/ba13-recovery/evidence/manifests/run1-raw.json',R/'ba13-recovery/run1/raw')]:
        entries=read(S/rel)['files']
        check('manifest '+str(root),all(sha(root/p)==(v['sha256'] if isinstance(v,dict) else v) for p,v in entries.items()))
    mx=read(F/'matrix.json');rows=mx['rows'];cert=read(R/'ba13-recovery/run1/fit/certificate.json')
    check('462 unique rows',len(rows)==462 and len({(r['index'],r['quantity']) for r in rows})==462)
    check('frozen historical identity',all(r['historical']=={k:cert['samples'][r['index']][k][r['quantity']] for k in ('drift_bound','threshold','utilization')} for r in rows))
    check('headrooms unchanged',all(r['headroom']==r['historical']['threshold']-r['historical']['drift_bound'] for r in rows))
    check('only exact P2 passes',all((r['status']=='PASS')==(r['quantity']=='P2_average') for r in rows))
    check('no incomplete totals',all(r['proposed_total_utilization'] is None for r in rows if r['status']=='UNQUALIFIED'))
    check('all missing values named',all(k in r['missing'] for r in rows for k,v in r['contributions'].items() if v['value'] is None))
    check('no fabricated owner authority',all(not t.get('authority') for r in rows for t in r['contributions'].values()))
    check('units of tangents absolute',all('absolute' in r['units'] and r['contributions']['normalization']['value']==0 for r in rows if r['quantity'].startswith('t_')))
    check('no historical estimator consumption','defect_transport_uv' not in (S/'tools/linux/ba13_closure/assess.py').read_text())
    # Independent Decimal calculation of the three central N_mu errors/orders.
    with decimal.localcontext() as ctx:
        ctx.prec=60
        central={i:{n:read((R/'ba13-headroom' if n==512 else O)/f't{i}-n{n}.json',True) for n in (512,1024,2048)} for i in (0,10,20)}
        raw=list(csv.DictReader((R/'ba13-recovery/run1/raw/target_stars.tsv').open(),delimiter='\t'))
        m1=read(F/'m1.json')
        p0=D(float(raw[0]['N_mu']))
        for i in (0,10,20):
            delta=[central[i][1024]['N'][3]-central[i][512]['N'][3],central[i][2048]['N'][3]-central[i][1024]['N'][3]]
            order=float((abs(delta[0]/delta[1])).ln()/D(2).ln())
            check(f'independent N_mu order {i}',abs(order-m1['reference'][str(i)]['N_mu']['order'][0])<1e-5)
        for i in (10,20):
            dp=(D(float(raw[i]['N_mu']))-p0)/p0
            dx=(central[i][2048]['N'][3]-central[0][2048]['N'][3])/p0
            check(f'independent N_mu paired {i}',abs(float(abs(dp-dx))-m1['paired']['N_mu'][str(i)]['paired_production_minus_X2'])<2e-16)
        for i in (0,10,20):
            hs=central[i][512]['primary'][8];tt=[]
            for w in ('0.0001','0.00005','0.000025'):
                folder=O if w=='0.000025' else R/'ba13-headroom'
                lo=read(folder/f't{i}-h{w}--1.json',True);hi=read(folder/f't{i}-h{w}-1.json',True)
                tt.append((hi['N'][3]-lo['N'][3])/(hi['primary'][0]-lo['primary'][0]))
            ratio=float((tt[1]-tt[0])/(tt[2]-tt[1]))
            check(f'independent third-width t_mu {i}',3<ratio<5 and abs(ratio-m1['reference'][str(i)]['t_mu']['stencil_ratio'])<.001)
    # Independent Neumann SUM (not the enclosure.py linear solve), applied to
    # every changing 3x3 Jacobian. The omitted series tail is negligible for
    # this verification, whose tolerance is not a scientific uncertainty.
    inv=read(F/'m5.json');L=np.array([[-1.,-1.],[1.,0.],[0.,1.]])
    worst=[0.,0.]
    for i in range(1,21):
        total=np.zeros((2,2));nonlinear=np.zeros((2,2))
        for j in (0,i):
            p=read(R/f'ba13-headroom/profile-v2-{j}.json');g=np.array(p['G_fixed_profile'])
            E=sum(np.array(p[k]) for k in ('background_inherited','center_inherited','tail_error'))
            a=np.linalg.inv(g);K=abs(a)@E;term=K@abs(a);series=term.copy();linear=term.copy()
            for _ in range(8):term=K@term;series+=term
            total+=abs(L).T@series@abs(L);nonlinear+=abs(L).T@(series-linear)@abs(L)
        recorded=np.array(inv[str(i)]['total_Z_matrix'])
        check(f'Neumann independent sum {i}',np.max(abs(total-recorded)/recorded)<5e-13)
        check(f'both nonlinear tails {i}',np.all(np.array(inv[str(i)]['nonlinear_Z_matrix'])>0))
        for row,n in enumerate(('Z_row_npe','Z_row_npmu')):worst[row]=max(worst[row],inv[str(i)]['row_majorants'][n])
    review=read(REV/'z-independent-endpoint.json')
    check('review Z majorant reproduced',all(abs(worst[k]/review[n]['independent_total']-1)<1e-8 for k,n in enumerate(('Z_row_npe','Z_row_npmu'))))
    m2=read(F/'m2.json');m3=read(F/'m3.json');m4=read(F/'m4.json')
    check('grid counterexample localized to surface',not m2['grid_identity'] and m2['interior_grid_identity'] and all(g['first_different_node']==g['nodes']-1 for g in m2['grids'][1:]))
    check('no transfer promoted',all(not v['transfer_validated'] for v in m2['transfer'].values()))
    check('all 20 intervals',all(len(v['intervals'])==20 for v in m3['curvatures'].values()))
    check('curvature remains conditional',all(not v['validated_supremum'] for v in m3['curvatures'].values()))
    check('all three tangent common B terms',all(max(m4['common_B'][n]['paired'])>0 for n in ('t_n','t_e','t_mu')))
    check('no fit covariance duplicate',all(not v['fit_covariance_added'] for v in m4['fit'].values()))
    git=R/'toolchain/bin/git'
    changed=subprocess.check_output([str(git),'-C',str(S),'diff','--name-only','286aafb63002093b3052e23a3972903011f36035','HEAD'],text=True).splitlines()
    allowed=('tools/linux/ba13_headroom/','tools/linux/ba13_closure/','docs/validation/ba13-closure/')
    exact={'tests/bnv/x2_uncertainty.cpp','tests/bnv/test_ba13_closure.py','docs/adr/ADR-0021-ba13-static-uncertainty-completion.md'}
    check('protected source/failure trees unchanged',all(p.startswith(allowed) or p in exact for p in changed))
    adr='docs/adr/ADR-0021-ba13-static-uncertainty-completion.md'
    old_adr=subprocess.check_output([str(git),'-C',str(S),'show','286aafb63002093b3052e23a3972903011f36035:'+adr])
    check('ADR historical bytes preserved as prefix',(S/adr).read_bytes().startswith(old_adr))
    result=dict(status='PASS',checks=checks,tests=tests,check_count=len(checks),scientific_qualification=False,
                method='independent Decimal functionals/stencils; Neumann series sum; complete source/manifest/decision audit',
                cpu_seconds=time.process_time()-start,
                numerical_cpu_seconds=time.process_time()-start+resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime+resource.getrusage(resource.RUSAGE_CHILDREN).ru_stime-children_start.ru_utime-children_start.ru_stime,
                wall_seconds=time.monotonic()-wall,
                peak_rss_kib=max(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss),
                job=os.environ['SLURM_JOB_ID'],source_sha256=sha(pathlib.Path(__file__)))
    with (O/'audit.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
