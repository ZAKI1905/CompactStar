from pathlib import Path
import shutil,gzip,hashlib,json
q=Path(__file__).resolve().parent
repo=Path('/Users/keeper/Documents/CompactStar/worktrees/CompactStar-external-deps-taskmanager-migration')
dest=repo/'docs/validation/evidence/external-deps-migration/resume';dest.mkdir(exist_ok=True)
selected=['package-discovery-controls.json','integrated-zaki-characterization.json','unpatched-replay-control.json','final-provider-audit.json','final-link-audit.json','taskmanager-final-image-inventory.json','active-source-audit.json','Release-baseline-comparison.json','new-baseline-reproduction.json','running-suite-binary-hashes.json']
for name in selected:shutil.copyfile(q/'evidence'/name,dest/name)
for f in sorted((q/'evidence').glob('G0-*.log')):shutil.copyfile(f,dest/f.name)
for f in sorted((q/'evidence').glob('final-*.gz')):shutil.copyfile(f,dest/f.name)
for mode in ['Debug','Release']:
 for fixture in ['T1','T2']:
  f=q/f'old-plot-free-{mode}/{fixture}-result.json';shutil.copyfile(f,dest/f'old-plot-free-{mode}-{fixture}.json')
 for f in (q/f'new-{mode}/tests/dependency_migration').glob('*-result.json'):shutil.copyfile(f,dest/f'new-{mode}-{f.name}')
 for name in ['compile-link-commands.json']:
  for stack in ['old','new']:
   f=q/'evidence'/f'{stack}-{mode}-{name}';(dest/(f.name+'.gz')).write_bytes(gzip.compress(f.read_bytes(),mtime=0))
 for f in [q/f'provider-{mode}/old.map',q/f'provider-{mode}/old-disassembly.txt',q/f'provider-{mode}/old-nm.txt']:
  (dest/f'{mode}-{f.name}.gz').write_bytes(gzip.compress(f.read_bytes(),mtime=0))
for short in ['D','R']:
 for label in ['old','t1-old','ob-old']:
  f=Path(f'/private/tmp/csm-{label}-{short}.json');shutil.copyfile(f,dest/f.name)
for name in ['audit_final_links.py','audit_final_providers.py','qualify_plot_free_old.py','qualify_release_baselines.py','qualify_zaki_integrated.py','test_package_discovery.py','check_new_baselines.py','run_adr0017.py','run_adr0017_old.py','compare_performance.py']:
 shutil.copyfile(q/name,dest/name)
(dest/'README.md').write_text('''# Resume implementation evidence

These records supplement the immutable earlier stop/probe evidence in the
parent directory. The owner resume authorizes scientific/functional gates,
not universal object-code or symbol-provider identity. Final suite/ADR0017
results and candidate disposition are recorded in the migration report.

Scripts are local-Mac qualification recipes. They intentionally identify
this task's source repositories and exact pinned packages. Copy the scripts
into the external execution root before running them, so `Path(__file__).parent`
resolves to `/Users/keeper/Documents/CompactStar/external/qualification/compactstar-migration/e263a6e-b0cbd510`.
They never download or update source repositories. Build OLD canonical and
NEW candidate in separate `old-{Debug,Release}` and `new-{Debug,Release}`
directories; configure NEW as documented in `docs/build/EXTERNAL_PACKAGES_MAC.md`.
The package-build and historical probe recipes are in the parent directory.

`csm-old-*` and `csm-t1-old-*` contain original OLD artifact manifests;
`old-plot-free-*` compare the two mechanical plotting-deletion TUs against
those references. `new-*-result.json` contain the actual candidate comparisons.
The committed test fixtures retain compressed OLD numerical files, so replay
does not depend on temporary successful NEW output directories.

`integrated-zaki-characterization.json` records every fixed/stress comparison
and full output hashes; the large TSVs stay in the external execution root.
Maps, symbols and disassembly are compressed text, not generated binary trees.
The raw build/install logs and package trees remain in that same external root.
A SHA256SUMS file authenticates the compact committed evidence.
''')
(dest/'SHA256SUMS').write_text(''.join(f'{hashlib.sha256(f.read_bytes()).hexdigest()}  {f.relative_to(dest)}\n' for f in sorted(dest.rglob('*')) if f.is_file() and f.name!='SHA256SUMS'))
print(len(list(dest.iterdir())),'files',sum(f.stat().st_size for f in dest.iterdir()),'bytes')
