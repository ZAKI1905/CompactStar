#!/usr/bin/env python3
"""Package existing evidence only. No solver, re-scoring or scientific calculation."""
import gzip
import hashlib
import json
import pathlib
import shutil
import subprocess

S=pathlib.Path(__file__).resolve().parents[3];R=S.parents[1]
O=R/'ba13-evidence-20261026';DOC=S/'docs/validation/ba13-evidence';EV=DOC/'evidence'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def write(p,data):p.write_text(json.dumps(data,indent=1,allow_nan=False)+'\n')


def main():
    receipt={p.name:read(p) for p in sorted(O.glob('receipt-*.json'))}
    resources=dict(numerical_cpu_seconds=sum(v['numerical_cpu_seconds'] for v in receipt.values()),
        initial_plan_cpu_seconds=[40,60],ceiling_cpu_seconds=180,
        within_initial_plan=40<=sum(v['numerical_cpu_seconds'] for v in receipt.values())<=60,
        builds_cpu_seconds=0,execution='authenticated retained executables',
        phases={k:{f:v[f] for f in ('job','status','numerical_cpu_seconds','wall_seconds','affinity')} for k,v in receipt.items()},
        scheduler_accounting='sacct failed: Slurm accounting storage is disabled',
        available_receipts='scontrol job allocations and per-process self+children getrusage, tied to Slurm job IDs',
        stopped_static_attempt='19135: audit SyntaxError before any test or solver; preserved; corrected static run 19136 passed',
        packaging='hashing, compression, copying and report formatting only; outside numerical scope')
    if resources['numerical_cpu_seconds']>180:raise RuntimeError('resource ceiling exceeded')
    write(O/'resources.json',resources)
    write(O/'scheduler-accounting.json',dict(command='sacct -j 19130,19131,19132,19133,19134,19135,19136,19137 --parsable2 --units=K --format=JobIDRaw,JobName,State,ExitCode,ElapsedRaw,TotalCPU,UserCPU,SystemCPU,MaxRSS,AllocCPUS,NodeList',
          returncode=1,stderr='Slurm accounting storage is disabled',fabricated_accounting=False))
    write(O/'github-authentication.json',dict(noninteractive_credential_available=False,credential_fill_returncode=128,
          credential_values_logged=False,push_attempted=False,delivery='verified incremental Git bundle'))
    mx=read(O/'matrix.json')
    lines=['# BA13 complete matrix','',
      'Full 21-by-22 data: [matrix.tsv.gz](evidence/matrix.tsv.gz) and [matrix.json.gz](evidence/matrix.json.gz).',
      'Every row retains historical drift, threshold, utilization and assessment status. Unknown total justified error and physical headroom are null.',
      'The table below shows the row with largest finite-added/headroom fraction for each quantity. Finite sums include unaccepted numerical and conditional estimates; the remaining column is partial headroom after those terms, not a complete physical enclosure.','',
      '| Quantity | Sample | Historical utilization | Finite added terms | Justified total added error | Partial remaining headroom | Status |',
      '|---|---:|---:|---:|---:|---:|---|']
    names=list(dict.fromkeys(r['quantity'] for r in mx['rows']))
    for n in names:
        r=max((r for r in mx['rows'] if r['quantity']==n),key=lambda r:r['finite_added_estimates']/r['headroom'])
        justified='UNRESOLVED' if r['added_justified_error'] is None else f"{r['added_justified_error']:.6g}"
        lines.append(f"| {n} | {r['index']} | {r['historical']['utilization']:.9g} | {r['finite_added_estimates']:.6g} | {justified} | {r['remaining_headroom_after_finite_terms']:.6g} | {r['status']} |")
    lines+=['','For every individual sample, the full files list each finite term, evidence class, source hashes, accepted-authority field, missing remainder and conditional dependency.',
       'All approximate X2 rows retain unresolved reference stage/arithmetic, unseen bias, finite-sensitivity remainders, moving-cut input and nonlinear fit/input terms. Interpolation, positive-B alignment, stencil/normalization and continuous curvature are explicit conditionals. The registry accepts none.',
       'mu_n_profile additionally records the equilibrium-definition condition; metric_structure lacks R-Q nodewise evidence; species_support lacks continuous no-crossing acceptance. Only P2 is an exact closed identity.','']
    (DOC/'MATRIX.md').write_text('\n'.join(lines))
    # External manifest includes all new regular files, never dereferences the
    # historical audit mirror's symlinks or mutates its immutable inputs.
    files={}
    for p in sorted(O.rglob('*')):
        if not p.is_file() or p.is_symlink() or p.name in ('external-manifest.json','COMPLETE.json'):continue
        files[str(p.relative_to(O))]=dict(sha256=sha(p),bytes=p.stat().st_size)
    manifest=dict(root=str(O),files=files,symlinks='historical inputs excluded; bound by inherited manifests')
    write(O/'external-manifest.json',manifest)
    write(O/'COMPLETE.json',dict(matrix_sha256=sha(O/'matrix.json'),external_manifest_sha256=sha(O/'external-manifest.json'),
          row_count=462,owner_ratified=False,certificate_promoted=False,scientific_qualification='UNQUALIFIED'))
    EV.mkdir(exist_ok=True)
    for rel in [*files,'external-manifest.json','COMPLETE.json']:
        p=O/rel;target=EV/rel;target.parent.mkdir(parents=True,exist_ok=True)
        if rel in ('matrix.json','matrix.tsv','comparison.json'):
            target=target.with_suffix(target.suffix+'.gz')
            target.write_bytes(gzip.compress(p.read_bytes(),mtime=0))
        else:shutil.copyfile(p,target)
    paths=sorted((S/'tools/linux/ba13_evidence').glob('*'))+[
        S/'tools/linux/ba13_closure/assess.py',S/'tests/bnv/test_ba13_evidence.py',
        S/'tests/bnv/x2_fixed_step.cpp',S/'tests/bnv/x2_uncertainty.cpp',S/'tests/bnv/compare_candidate.py',
        S/'tools/linux/ba13_refit.py',S/'docs/adr/ADR-0021-ba13-static-uncertainty-completion.md']+list(DOC.glob('*.md'))+[DOC/'acceptance-registry.json']
    write(EV/'source-hashes.json',dict(files={str(p.relative_to(S)):sha(p) for p in paths if p.is_file()},
          assembly_source_commit='2cac8b0',final_reference_change='output-root override only',
          packaging_source_commit=subprocess.check_output(['git','-C',str(S),'rev-parse','HEAD'],text=True).strip(),
          review_manifest_sha256=sha(R/'ba13-closure-independent-review-20261020/SHA256SUMS')))
    entries=[]
    for p in sorted(EV.rglob('*')):
        if p.is_file() and p.name!='SHA256SUMS':entries.append(f'{sha(p)}  {p.relative_to(S)}')
    (EV/'SHA256SUMS').write_text('\n'.join(entries)+'\n')
    print(json.dumps(dict(external_files=len(files),compact_files=len(entries),numerical_cpu_seconds=resources['numerical_cpu_seconds'])))


if __name__=='__main__':main()
