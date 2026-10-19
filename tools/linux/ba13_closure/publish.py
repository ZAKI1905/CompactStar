#!/usr/bin/env python3
"""Package existing assessment evidence; never scores a new scientific run."""
import gzip,hashlib,json,pathlib,shutil,subprocess
S=pathlib.Path(__file__).resolve().parents[3];R=S.parents[1]
O=R/'ba13-closure-20261019';F=O/'assessment-final';DOC=S/'docs/validation/ba13-closure';E=DOC/'evidence'
def read(p):return json.loads(p.read_text())
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def save(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')

def main():
    assert read(O/'audit.json')['status']=='PASS'
    E.mkdir(exist_ok=True)
    for name in ('summary.json','m1.json','m2.json','m3.json','m4.json','m5.json'):
        shutil.copyfile(F/name,E/name)
    shutil.copyfile(O/'audit.json',E/'audit.json')
    for p in O.glob('test_ba13_*.out'):shutil.copyfile(p,E/p.name)
    for name in ('matrix.json','matrix.tsv'):
        (E/(name+'.gz')).write_bytes(gzip.compress((F/name).read_bytes(),mtime=0))
    records=[]
    for p in sorted(O.rglob('receipt-*.json')):
        data=read(p)
        records.append(dict(path=str(p.relative_to(O)),numerical_cpu_seconds=data['numerical_cpu_seconds'],
                            build_cpu_seconds=data.get('build_cpu_seconds',0),wall_seconds=data['wall_seconds'],
                            peak_rss_kib=data['peak_rss_kib'],job=data.get('job'),status=data.get('status','completed'),sha256=sha(p)))
    audit=read(O/'audit.json')
    records.append(dict(path='audit.json',numerical_cpu_seconds=audit['numerical_cpu_seconds'],build_cpu_seconds=0,
                        wall_seconds=audit['wall_seconds'],peak_rss_kib=audit['peak_rss_kib'],job=audit['job'],status=audit['status'],sha256=sha(O/'audit.json')))
    resources=dict(records=records,numerical_cpu_seconds=sum(r['numerical_cpu_seconds'] for r in records),
                   build_cpu_seconds=sum(r['build_cpu_seconds'] for r in records),
                   measured_wall_seconds_sum=sum(r['wall_seconds'] for r in records),
                   peak_rss_kib=max(r['peak_rss_kib'] for r in records),
                   initial_numerical_cpu_ceiling=900,total_numerical_cpu_ceiling=1800,
                   extension_numerical_cpu_seconds=0,new_blanket_sequence_solves=0,
                   scope='all newly scored stellar/static/audit runs, including stopped schema replay and both static assemblies',
                   excluded='initial historical fitter/authentication replay, login-node source edits/hash/packaging/startup and Slurm connectivity probes were not metered; no exact whole-session CPU claim',
                   stopped_run='receipt-benchmark.json says RUNNING because its then-current finally block recorded that before the schema exception; final-slurm-18965.txt and traceback retain actual FAILED state',
                   allocation='one pinned logical core; Slurm may allocate both hardware threads of that core; see allocation receipts')
    save(E/'resources.json',resources);save(O/'resources.json',resources)
    m1=read(F/'m1.json');m2=read(F/'m2.json');m3=read(F/'m3.json');m5=read(F/'m5.json')
    summary=read(F/'summary.json')
    disp=dict(**summary,implementation=dict(D1='CORRECTED',D2='RETIRED',D3='CORRECTED',D4='CORRECTED',D5='CORRECTED',D6='CORRECTED'),
              milestones={
              'M1':dict(status='INCONCLUSIVE',per_quantity=m1['paired'],reason='tangent endpoint reference dominates measured paired error; all reference bands lack named remainders'),
              'M2':dict(status='INCONCLUSIVE',full_grid_identity=m2['grid_identity'],interior_grid_identity=m2['interior_grid_identity'],reason='moving final node plus unvalidated all-target structural magnitude'),
              'M3':dict(status='CONDITIONAL',all_intervals_evaluated=20,narrower_width_test_satisfied=all(v['narrower_width_prediction_satisfied'] for v in m3['curvatures'].values()),reason='factor 10 hypothesis is tested, not a uniform derivative enclosure'),
              'M4':dict(status='PARTIAL_EVIDENCE',reason='functional/common-B/denominator/design terms evaluated; numerical and residual-input remainders not established'),
              'M5':dict(status='ALGEBRA_VERIFIED_CONDITIONAL_INPUT',majorants={n:max(v['row_majorants'][n] for v in m5.values()) for n in ('Z_row_npe','Z_row_npmu')},reason='per-target inherited background magnitude hypothesis still needs validation/acceptance'),
              'M6':dict(status='PASS_IMPLEMENTATION',existing_tests=17,new_tests=13,independent_checks=audit['check_count'],scientific_acceptance=False)},
              owner_decisions=['inter-sample remainder','accepted evidence classes/standard','physical response versus production output','sidecar versus successor runtime record','magnitude transfer across targets'],
              observer=dict(feasible_in_principle=True,implemented=False,reason='M1-M3 prerequisite/transfer stop',
                            method='passive existing step observer; retain collocation polynomials and evaluate original support coordinates; verify dense-output order',
                            estimated_three_target_512_1024_2048_solver_cpu_seconds=12,observer_overhead_benchmarked=False),
              no_promotion=True)
    save(E/'dispositions.json',disp);save(O/'dispositions.json',disp)
    code=list((S/'tools/linux/ba13_closure').glob('*.py'))+list((S/'tools/linux/ba13_closure').glob('*.sbatch'))
    code += [S/p for p in ('tools/linux/ba13_headroom/enclosure.py','tools/linux/ba13_headroom/assess.py','tools/linux/ba13_refit.py',
                          'tests/bnv/x2_fixed_step.cpp','tests/bnv/x2_uncertainty.cpp','tests/bnv/compare_candidate.py',
                          'tests/bnv/test_ba13_closure.py','tests/bnv/test_ba13_headroom.py','tests/bnv/test_ba13_recovery.py',
                          'docs/adr/ADR-0021-ba13-static-uncertainty-completion.md')]
    code += list(DOC.glob('*.md'))
    source=dict(entry='286aafb63002093b3052e23a3972903011f36035',
                files={str(p.relative_to(S)):sha(p) for p in sorted(code)},
                scientific_assessor_sha256=read(F/'receipt-static.json')['source_sha256'],
                final_audit_sha256=audit['source_sha256'])
    source['receipt_source_note']='Execution receipts bind the exact earlier Git revisions. Later output-directory configuration changes only reproduction I/O, not equations/scoring; final source files are hashed separately here.'
    git=R/'toolchain/bin/git'
    source['changed_paths']=sorted(set(subprocess.check_output([str(git),'-C',str(S),'diff','--name-only',source['entry']],text=True).splitlines()+
        subprocess.check_output([str(git),'-C',str(S),'ls-files','--others','--exclude-standard'],text=True).splitlines()))
    save(E/'source-hashes.json',source);save(O/'source-hashes.json',source)
    manifests=[('headroom','docs/validation/ba13-headroom/evidence/external-manifest.json',R/'ba13-headroom'),
               ('x2','docs/validation/x2-fixed-step/evidence/external-manifest.json',R/'x2-fixed-step'),
               ('nominal','docs/validation/ba13-recovery/evidence/manifests/run1-raw.json',R/'ba13-recovery/run1/raw')]
    inputs={}
    for name,rel,root in manifests:
        inputs[name]=dict(manifest_source=rel,manifest_sha256=sha(S/rel),root=str(root),files=read(S/rel)['files'])
    inputs['review']=dict(root=str(R/'ba13-headroom-independent-review-20261018'),manifest_sha256=sha(R/'ba13-headroom-independent-review-20261018/SHA256SUMS'))
    constants=R/'toolchain/install/Release/Zaki/include/Zaki/Physics/Constants.hpp'
    inputs['toolchain']=dict(constants_sha256=sha(constants),Zaki_library_sha256=sha(R/'toolchain/install/Release/Zaki/lib/libZaki.a'),
                             compiler_wrapper_sha256=sha(R/'toolchain/bin/g++'),python=str(R/'toolchain/python/bin/python'))
    save(E/'inputs.json',inputs);save(O/'inputs.json',inputs)
    mx=read(F/'matrix.json');names=list(dict.fromkeys(r['quantity'] for r in mx['rows']))
    lines=['# Complete frozen-threshold BA13 matrix','',
           'All 462 rows are in [matrix.json.gz](evidence/matrix.json.gz) and [matrix.tsv.gz](evidence/matrix.tsv.gz).',
           'The JSON carries terms, evidence classes, hypotheses, units, scales, source pointers, correlations and refusal reasons.',
           'Component values remain in the external component-data.json.gz, authenticated by the manifest.',
           'Finite partial utilization is a diagnostic; it is never a qualified total.','',
           '| Quantity | PASS | FAIL | UNQUALIFIED | Minimum headroom | Largest partial utilization |',
           '|---|---:|---:|---:|---:|---:|']
    for n in names:
        rows=[r for r in mx['rows'] if r['quantity']==n]
        counts=[sum(r['status']==v for r in rows) for v in ('PASS','FAIL','UNQUALIFIED')]
        lines.append(f'| {n} | {counts[0]} | {counts[1]} | {counts[2]} | {min(r["headroom"] for r in rows):.9g} | {max(r["partial_utilization"] for r in rows):.9g} |')
    lines+=['','Total: **21 PASS / 0 FAIL / 441 UNQUALIFIED**. No approximate group has complete evidence.',
            'Muon-sector numerical rows: **84 UNQUALIFIED**; exact support/domain rows: **21 UNQUALIFIED**.',
            'Disposition **E**; [scientific report](REPORT.md) and [owner decisions](OWNER_DECISIONS.md).','']
    (DOC/'MATRIX.md').write_text('\n'.join(lines))
    # Include generated MATRIX in the exact source map after creating it.
    source['files']['docs/validation/ba13-closure/MATRIX.md']=sha(DOC/'MATRIX.md')
    save(E/'source-hashes.json',source);save(O/'source-hashes.json',source)
    files={str(p.relative_to(O)):dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(O.rglob('*'))
           if p.is_file() and p.name not in ('external-manifest.json','COMPLETE.json')}
    manifest=dict(schema='compactstar.ba13.closure.artifacts.v1',root=str(O),files=files,scientific_qualification=False)
    save(E/'external-manifest.json',manifest);shutil.copyfile(E/'external-manifest.json',O/'external-manifest.json')
    complete=dict(external_manifest_sha256=sha(E/'external-manifest.json'),matrix_sha256=sha(F/'matrix.json'),
                  status='ASSESSMENT_COMPLETE_SCIENTIFICALLY_UNQUALIFIED',counts=summary['counts'])
    save(O/'COMPLETE.json',complete);save(E/'COMPLETE.json',complete)
    # Verify every external byte now, and compact matrix identity separately.
    assert all(sha(O/p)==v['sha256'] for p,v in files.items())
    assert all(gzip.decompress((E/(n+'.gz')).read_bytes())==(F/n).read_bytes() for n in ('matrix.json','matrix.tsv'))
    print(json.dumps(dict(external_files=len(files),numerical_cpu_seconds=resources['numerical_cpu_seconds'],
                         build_cpu_seconds=resources['build_cpu_seconds'],peak_rss_kib=resources['peak_rss_kib'],counts=summary['counts']),indent=2))

if __name__=='__main__':main()
