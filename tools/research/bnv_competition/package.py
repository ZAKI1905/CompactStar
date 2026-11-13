#!/usr/bin/env python3
"""Package compact results and authenticate external evidence; no scientific solve."""
import argparse,csv,hashlib,json,shutil
from pathlib import Path

ROOT=Path(__file__).resolve().parents[3]
p=argparse.ArgumentParser();p.add_argument('out',type=Path);a=p.parse_args()
dest=ROOT/'docs/research/bnv-rotochemical-competition'
for name in ['data','figures','evidence/short','evidence/logs']:
 (dest/name).mkdir(parents=True,exist_ok=True)
for name in ['map.csv','controls.csv','fiducial.csv','cancellation-zoom.csv','representative-periods.csv']:
 (dest/'data'/name).write_bytes((a.out/name).read_bytes().replace(b'\r\n',b'\n'))
# Preserve every BNV-only matched state in compact form, in the same row order
# as map.csv; full stationary outputs remain outside Git.
keys=['Gamma_yr','D_count_s','Tinf_K','Tsurface_inf_K','eta_e_MeV','eta_mu_MeV',
      'Lgamma_erg_s','LH_erg_s','DeltaLnu_erg_s','net_weak_erg_s','Lnu_eq_erg_s',
      'P_direct_erg_s','A_B_erg_s','A_spin_erg_s','Echem_erg','Lescape_star_erg_s',
      'tau_slow_yr','build_time_yr','depletion_build','stable','chemical_residual',
      'thermal_residual','subsystem_residual']
with (a.out/'bnv-only.csv').open() as f,(dest/'data/bnv-only.csv').open('w') as g:
 w=csv.DictWriter(g,keys,lineterminator='\n');w.writeheader()
 for r in csv.DictReader(f):w.writerow({k:r[k] for k in keys})
for f in (a.out/'figures').iterdir():
 if f.suffix in ['.pdf','.png']:shutil.copyfile(f,dest/'figures'/f.name)
for name in ['summary.json','verification.json','supplement.json','short-invocation.json',
             'short-completion.json','plot-receipt.json']:
 shutil.copyfile(a.out/name,dest/'evidence'/name)
for f in (a.out/'short').iterdir():
 if f.is_file():shutil.copyfile(f,dest/'evidence/short'/f.name)
for name in ['contracts.log','short.log','check-19423.log','short-19427.log',
             'supplement-19428.log']:
 shutil.copyfile(a.out/name,dest/'evidence/logs'/name)
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(1<<20),b''):h.update(block)
 return h.hexdigest()
paths=[]
for run in ['r1','r2','r3']:
 directory=a.out.parent/('bnv-rotochemical-competition-20261112-'+run)
 def visit(path):
  for f in sorted(path.iterdir()):
   if f.name in ['python-deps','matplotlib-cache','build','tmp','__pycache__']:continue
   if f.is_dir():visit(f)
   elif f.is_file():paths.append(f)
 visit(directory)
for name in ['bnv_competition_short','bnv_competition_contract']:
 paths.append(a.out/'build/tests'/name)
manifest=a.out.parent/'bnv-rotochemical-competition-20261112-external-SHA256SUMS'
manifest.write_text(''.join(digest(f)+'  '+str(f)+'\n' for f in sorted(paths)))
receipt=dict(manifest=str(manifest),sha256=digest(manifest),files=len(paths),
             scope='All r1/r2/r3 scientific outputs and logs plus new native executables; excludes reproducible build intermediates and dependency/cache directories.',
             accepted_static_directory=str(a.out),
             retained_archive=dict(path=str(a.out.parent/'native-p2-regeneration-20261108-r1/capture/sequence.archive'),
               sha256='26dd2b060ad8ed10d26ec46b67a9b957e3517fe6718708ba91d09c5f204158cf'),
             scientific_status='candidate; no broader certificate or human scientific review')
(dest/'evidence/external-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
shutil.copyfile(manifest,dest/'evidence/external-SHA256SUMS')
source=[]
for directory in ['CompactStar/Physics/BNV','CompactStar/Physics/Rotochemical','tools/research/bnv_competition']:
 source += [f for f in (ROOT/directory).rglob('*') if f.is_file() and '__pycache__' not in f.parts]
source += [ROOT/'tests/native/competition.cpp',ROOT/'tests/bnv/competition_contract.cpp',ROOT/'tests/CMakeLists.txt']
(dest/'evidence/final-source-sha256.json').write_text(json.dumps({str(f.relative_to(ROOT)):digest(f) for f in sorted(source)},indent=2)+'\n')
print(json.dumps(receipt,indent=2))
