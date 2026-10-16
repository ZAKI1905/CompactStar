#!/usr/bin/env python3
"""Publish compact derived evidence and authenticated external manifests."""
import gzip,hashlib,json,pathlib,shutil,subprocess
S=pathlib.Path(__file__).resolve().parents[3];R=S.parents[1];O=R/'ba13-headroom';D=S/'docs/validation/ba13-headroom/evidence';D.mkdir(exist_ok=True)
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
receipts=[read(p) for p in sorted(O.glob('receipt-*.json'))]
profiles=[read(p) for p in sorted(O.glob('profile-v2-*.json'))]
supp=read(O/'supplement.json');summary=read(O/'summary.json');audit=read(O/'audit.json')
# Failed serializer jobs had no completed process receipt. Durations were read
# from scontrol in the session (13,12,12 s), before records expired. Preserve
# this as a transcription, not invented raw scheduler output or exact CPU use.
failures=dict(source='scontrol show job output observed during session; raw final records expired before archival',jobs=[dict(job=17038,runtime_seconds=13),dict(job=17045,runtime_seconds=12),dict(job=17047,runtime_seconds=12)],charge='entire two-CPU allocations; conservative 74 CPU seconds',cpu_charge=74)
resources=dict(solver_controller_cpu_seconds=sum(r['cpu_seconds'] for r in receipts if not r['phase'].startswith('build')),
 build_cpu_seconds=sum(r['cpu_seconds'] for r in receipts if r['phase'].startswith('build')),
 nominal_profile_cpu_seconds=sum(r['cpu_seconds'] for r in profiles),supplement_cpu_seconds=supp['cpu_seconds'],
 final_assessment_cpu_seconds=summary['cpu_seconds'],final_audit_cpu_seconds=audit['cpu_seconds'],failed_profile_jobs=failures,
 receipts=receipts,profile_jobs=[{k:r[k] for k in ['index','job','cpu_seconds','wall_seconds','peak_rss_kib']} for r in profiles],
 max_recorded_rss_kib=max([r['peak_rss_kib'] for r in receipts+profiles]+[supp['peak_rss_kib'],summary['peak_rss_kib']]),
 numerical_ceiling_cpu_seconds=1800,substantial_work='single pinned logical CPU per Slurm task; arrays for independent work',
 controls='static executable 180 CPU seconds/200 wall seconds; profile/supplement jobs Slurm 10 minute wall limit; benchmarked 3-target profile task before 18-profile expansion',
 accounting_limit='CPU is retained measurements plus explicit failed-allocation charge, not exact whole-session CPU. Earlier local assembly/audits, source hashing and interpreter startup outside timers are not included.',
 failed_checks=['17036: fixed-H diagnostics differ from shooting-history maxima; primary physical outputs agree; original-history replay passes',
 '17038/17045/17047: NumPy int64 serialization after calculation; partial files retained, corrected write validates JSON before opening final output',
 'signed N=128 defect correction has wrong sign for R; not used as an error bound'])
resources['retained_numerical_cpu_seconds_plus_failed_charge']=sum(resources[k] for k in ['solver_controller_cpu_seconds','nominal_profile_cpu_seconds','supplement_cpu_seconds','final_assessment_cpu_seconds','final_audit_cpu_seconds'])+74
save(O/'resources.json',resources)
for name in ['summary.json','audit.json','resources.json','pilots.json','inner-half.json','tests.txt']:
 shutil.copyfile(O/name,D/name)
for name in ['matrix.json','matrix.tsv']:
 (D/(name+'.gz')).write_bytes(gzip.compress((O/name).read_bytes(),mtime=0))
# Compact human-readable 21-by-22 status table; complete numeric rows in matrix.tsv.gz.
m=read(O/'matrix.json');names=list(dict.fromkeys(r['quantity'] for r in m['rows']));lines=['# Proposed BA13 qualification matrix','','All thresholds and historical values are unchanged. P = exact PASS; U = UNQUALIFIED.','Full numeric rows, contributions, units, evidence types and individual gaps: `evidence/matrix.json.gz` and `evidence/matrix.tsv.gz`.','','| Target | '+' | '.join(names)+' |','|---|'+'---|'*len(names)]
for i in range(21):lines.append('| '+str(i)+' | '+' | '.join('P' if next(r for r in m['rows'] if r['index']==i and r['quantity']==n)['status']=='PASS' else 'U' for n in names)+' |')
(D.parent/'MATRIX.md').write_text('\n'.join(lines)+'\n')
manifest=dict(root=str(O),source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=S,text=True).strip(),files={str(p.relative_to(O)):dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(O.rglob('*')) if p.is_file() and p.name not in ['COMPLETE.json']})
save(D/'external-manifest.json',manifest)
save(O/'COMPLETE.json',dict(schema_id='compactstar.ba13.headroom-assessment-completion.v1',qualified=False,status='COMPLETE_ASSESSMENT_WITH_EXPLICIT_GAPS',recommendation='B',rows=462,external_manifest_sha256=sha(D/'external-manifest.json'),source_commit=manifest['source_commit']))
print(json.dumps({k:v for k,v in resources.items() if k not in ['receipts','profile_jobs']},indent=2))
