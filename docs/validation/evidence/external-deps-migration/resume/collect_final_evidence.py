from pathlib import Path
import subprocess,sys,shutil,gzip,hashlib,json
q=Path(__file__).resolve().parent
repo=Path('/Users/keeper/Documents/CompactStar/worktrees/CompactStar-external-deps-taskmanager-migration')
dest=repo/'docs/validation/evidence/external-deps-migration/resume'
subprocess.run([sys.executable,q/'summarize_final_science.py'],check=True)
subprocess.run([sys.executable,q/'summarize_test_suites.py'],check=True)
subprocess.run([sys.executable,q/'audit_protected_authority.py'],check=True)
subprocess.run([sys.executable,q/'audit_final_install.py'],check=True)
subprocess.run([sys.executable,q/'collect_static_evidence.py'],check=True)
selected=['canonical-final-identities.json','install-archive-diagnostic.json','significant-commits.json','final-install-audit.json','ADR0017-process-timings.json','test-suite-matrix.json','runtime-library-versions.json','install-probe-preservation.json','ADR0017-same-mode-comparison.json','Phase5D1-complete-matrix.json','protected-authority.json','paired-performance.json','final-build-checks.json','final-image-changes.json','final-all-images.json','final-concurrent-heat-capacity.json','final-Release-applicable-suite.json','Release-debug-authority-diagnostics.json','Release-extra-reference-comparison.json','old-Release-repaired-structural.json','phase5d1-install-provenance-diagnostic.json','test-registration-summary.json','heat-concurrency-repair.json','T1-without-external-EOS.json']
for name in selected:shutil.copyfile(q/'evidence'/name,dest/name)
for pattern in ['*-full-suite.log','*phase5d1-repaired.log','final-*-tests.log','final-*-suite.log','final-*-configure.log','final-*-build.log','final-*-install.log','*same-mode.log','*structural-repaired.log','*heat-capacity-serial.log']:
 for f in sorted((q/'evidence').glob(pattern)):(dest/(f.name+'.gz')).write_bytes(gzip.compress(f.read_bytes(),mtime=0))
for stack in ['old','new']:
 for mode in ['Debug','Release']:
  for i,f in enumerate(sorted((q/f'{stack}-{mode}').rglob('regression-evidence.json'))):
   label=f'phase5d1-{stack}-{mode}-{i}'
   shutil.copyfile(f,dest/(label+'-receipt.json'))
   receipt=json.loads(f.read_text());sidecar=Path(receipt['execution_sidecar'])
   shutil.copyfile(sidecar,dest/(label+'-execution-sidecar.json'))
   artifact=Path(receipt['scratch_root'])/'governed-artifact.json'
   assert hashlib.sha256(artifact.read_bytes()).hexdigest()==receipt['baseline_sha256']
   (dest/(label+'-governed-artifact.json.gz')).write_bytes(gzip.compress(artifact.read_bytes(),mtime=0))
for stack,name in [('old','adr0017-old-Debug'),('new','adr0017-new-Debug-attempt2')]:
 root=q/name
 for rel in ['result.json','commands.json','auth/authentication.json']:
  f=root/rel
  if f.exists():shutil.copyfile(f,dest/f'ADR0017-{stack}-{f.name}')
 for f in sorted((root/'output').iterdir()):
  if f.is_file():(dest/f'ADR0017-{stack}-{f.name}.gz').write_bytes(gzip.compress(f.read_bytes(),mtime=0))
for name in ['collect_static_evidence.py','collect_final_evidence.py','summarize_final_science.py','summarize_test_suites.py','audit_protected_authority.py','audit_final_install.py','qualify_release_extra_references.py','audit_images.py','finish_build_matrix.py','finish_test_matrix.py']:
 shutil.copyfile(q/name,dest/name)
readme=dest/'README.md'
readme.write_text(readme.read_text()+'''\nThe initial Release suite logs intentionally retain failed cross-mode Debug\nreference checks and the repaired infrastructure/provenance failures. The\nfinal applicable suite and separate complete Phase-5D1 rerun are authoritative;\nsee the report for their exact inventory. OLD/NEW Debug complete-suite logs\ninclude every long test. Phase-5D1 receipts authenticate byte-identical\ngoverned artifacts and ten fail-closed controls per run. ADR-0017 files retain\nwhole output tables, with timing fields explicitly separated from science.\n''')
(dest/'SHA256SUMS').write_text(''.join(f'{hashlib.sha256(f.read_bytes()).hexdigest()}  {f.relative_to(dest)}\n' for f in sorted(dest.rglob('*')) if f.is_file() and f.name!='SHA256SUMS'))
print(len(list(dest.iterdir())),'files',sum(f.stat().st_size for f in dest.iterdir()),'bytes')
