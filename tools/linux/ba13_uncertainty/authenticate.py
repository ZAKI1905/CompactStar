from pathlib import Path
import sys,json,hashlib,subprocess,gzip
s=Path('/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/ba13-uncertainty');r=s.parents[1];a=r/'ba13-recovery'
sys.path.insert(0,str(s/'tools/linux'));from ba13_stage import verify,digest
subprocess.run(['git','-C',str(s),'merge-base','--is-ancestor','8e6ac589d8d34e78447ee893e43e5673490eb728','HEAD'],check=True)
# Scientific owners must remain byte-identical to the authenticated entry.
assert not subprocess.check_output(['git','-C',str(s),'diff','8e6ac58','--','CompactStar','cmake','tests/eos','tests/analysis','tests/bnv/compare_candidate.py','tests/bnv/frozen_sensitivity.cpp','tests/bnv/produce_frozen_certificates.py'],text=True)
for run in ['run1','repeat']:
 verify(a/run/'raw');verify(a/run/'fit')
assert all((a/'run1/raw'/n).read_bytes()==(a/'repeat/raw'/n).read_bytes() for n in ['target_stars.tsv','tangent.tsv','coefficients.tsv','thermal_surface.tsv','profile_metrics.tsv'])
for n in ['certificate.json','certificate.tsv']:assert (a/'run1/fit'/n).read_bytes()==(a/'repeat/fit'/n).read_bytes()
for name in ['native','mac']:
 d=json.loads((a/(name+'-independent-refit.json')).read_text());assert d['fit_count']==10303 and not d['failures']
 with gzip.open(a/(name+'-independent-refit.json.fits.json.gz'),'rt') as f:assert len(json.load(f))==10303
archive=r.parents[0]/'incoming/BA13-original-evidence-20261012.zip'
assert digest(archive)=='990d3708a24c71e518d4bb59e6ad3c23ae47f337ba924007b083ca598e6b8e2b'
h=a/'historical-package/BA13-original-evidence-20261012'
for line in (h/'SHA256SUMS').read_text().splitlines():
 sha,name=line.split('  ',1);assert digest(h/name)==sha
lock=r/'toolchain/dependency-lock-linux.json';assert digest(lock)==json.loads((a/'build-receipts/dependency-provenance.json').read_text())['lock_sha256']
for x in json.loads((s/'docs/validation/ba13-recovery/recovered-identities.json').read_text()):assert digest(s/x['path'])==x['sha256']
heads={n:subprocess.check_output(['git','-C',str(r/'sources'/n),'rev-parse','HEAD'],text=True).strip() for n in ['CompactStar','Zaki','CONFIND']}
record=dict(status='PASS_REAUTHENTICATED_STUDY_INPUTS',heads=heads,base='8e6ac589d8d34e78447ee893e43e5673490eb728',dependency_lock_sha256=digest(lock),historical_archive_sha256=digest(archive),refit_components_each=10303,files={str(p.relative_to(a)):digest(p) for p in [a/'run1/raw/COMPLETE.json',a/'repeat/raw/COMPLETE.json',a/'run1/fit/COMPLETE.json',a/'repeat/fit/COMPLETE.json',a/'native-independent-refit.json',a/'mac-independent-refit.json',a/'native-independent-refit.json.fits.json.gz',a/'mac-independent-refit.json.fits.json.gz']})
(s.parents[1]/'ba13-uncertainty/reauthentication.json').write_text(json.dumps(record,indent=2)+'\n')
print(record['status'])
