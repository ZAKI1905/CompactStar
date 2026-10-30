#!/usr/bin/env python3
"""Deterministic packaging of bounded static evidence; no scientific computation."""
import gzip,hashlib,json,os,pathlib,shutil,subprocess
S=pathlib.Path(__file__).resolve().parents[3];R=S.parents[1]
O=pathlib.Path(os.environ.get('BA13_FINAL_OUTPUT_ROOT',R/'ba13-final-decision-20261030'))
P=S/'docs/validation/ba13-final-decision';E=P/'evidence'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,d):p.write_text(json.dumps(d,indent=1,allow_nan=False)+'\n')
E.mkdir(exist_ok=True)
receipts=[json.loads(p.read_text()) for p in sorted(O.glob('receipt-*.json'))]
assert all(v['status']=='PASS' for v in receipts)
dump(O/'resources.json',dict(measured_numerical_cpu_seconds=sum(v['cpu_seconds'] for v in receipts),
    ceiling_cpu_seconds=600,stellar_evaluations=0,slurm_jobs=[v['slurm_job'] for v in receipts],
    budget_includes='static numeric assessment, tests, independent reference replay, source mutations; hash-only packaging excluded',
    receipts={p.name:sha(p) for p in O.glob('receipt-*.json')},
    scheduler_accounting='Slurm accounting storage is disabled; process CPU and allocation receipts available',
    initial_submission='sandbox DNS resolution of ekucn1 failed; retry loop interrupted, no job submitted; permitted scheduler submission then succeeded'))
# Raw large files stay external. Compact gzip payloads are reproducible (mtime=0).
for n in ('matrix.json','matrix.tsv','conditional-intervals.json'):
    (E/(n+'.gz')).write_bytes(gzip.compress((O/n).read_bytes(),mtime=0))
for p in sorted(O.iterdir()):
    if p.is_file() and p.name not in ('matrix.json','matrix.tsv','conditional-intervals.json','external-manifest.json') and p.suffix!='.py':
        shutil.copyfile(p,E/p.name)
manifest=dict(root=str(O),files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(O.iterdir()) if p.is_file() and p.name!='external-manifest.json'})
dump(E/'external-manifest.json',manifest);dump(O/'external-manifest.json',manifest)
summary=json.loads((O/'summary.json').read_text())
lines=['# BA13 final decision matrix','','Updated inventory: **462 rows (21 stars × 22 groups)**. Full [JSON](evidence/matrix.json.gz) and [TSV](evidence/matrix.tsv.gz).',
       'Historical values, thresholds, norms, terms, classes and provenance are retained; all missing full uncertainties stay null.',
       'This successor adds common-blocker identifiers and a joint required-accuracy region. It does not rescore old outcomes.','',
       '| Group | Historical max utilization | Max finite addition / headroom | Minimum partial remaining | Required pointwise-error cap* |',
       '|---|---:|---:|---:|---:|']
for n,g in summary['groups'].items():
    cap=g['min_required_pointwise_cap'];lines.append(f"| {n} | {g['historical_utilization']:.9g} | {100*g['worst_finite_fraction']:.6g}% | {g['partial_remaining']:.8g} | {cap:.8g} |" if cap is not None else f"| {n} | {g['historical_utilization']:.9g} | {100*g['worst_finite_fraction']:.6g}% | {g['partial_remaining']:.8g} | — |")
lines += ['', '*Caps are single-axis feasibility limits, not estimates: additional gap/fit errors are set to zero for this display, and inherited finite terms are assumed valid. The full joint inequality is in every applicable JSON row and METHODS.md. Tangents are absolute dN_i/dB; other numeric caps use the historical scaled component norm. Metric structure has only P−R, not an independent reference discrepancy. Species support is categorical.',
          '', 'Disposition: **21 exact-P2 PASS, 0 FAIL, 441 UNQUALIFIED**. Approximate full utilization remains null. No finite sum is a complete qualification. Source: `evidence/summary.json:1`; formulas: `METHODS.md`, final section.']
(P/'MATRIX.md').write_text('\n'.join(lines)+'\n')
print(json.dumps(dict(external_entries=len(manifest['files']),cpu=sum(v['cpu_seconds'] for v in receipts),compressed_matrix_bytes=(E/'matrix.json.gz').stat().st_size)))
