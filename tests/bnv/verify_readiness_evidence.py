#!/usr/bin/env python3
"""Verify completed diagnostic output and seal its executable/source identities."""
import csv,hashlib,json,pathlib,subprocess
ROOT=pathlib.Path(__file__).resolve().parents[2]
OUT=ROOT/'docs/validation/evidence/phase6-campaign-readiness'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    record={'commands':[],'files':{},'accepted':'0c3f9269c884161c9d89163660a47b141a93ace6'}
    for name,source in [('breakpoint_contract','tests/bnv/breakpoint_contract.cpp'),('linear_probe','tests/bnv/readiness_linear_probe.cpp'),('edge_probe','tests/bnv/readiness_edge_probe.cpp')]:
        exe='build/readiness-standalone/'+name
        cmd=['/usr/bin/clang++','-std=c++17','-g','-I.','-I/opt/local/include',source,'-L/opt/local/lib','-lgsl','-lgslcblas','-o',exe]
        # Record the build invocation, but do not relink a sealed run executable.
        record.setdefault('build_commands',[]).append(cmd)
        for p in [source,exe]:record['files'][p]=sha(ROOT/p)
        if name!='linear_probe':
            r=subprocess.run([str(ROOT/exe)],cwd=ROOT,text=True,capture_output=True)
            record['commands'].append({'command':[exe],'returncode':r.returncode,'stdout':r.stdout,'stderr':r.stderr});assert r.returncode==0
            if name=='edge_probe':assert r.stdout==(OUT/'edge-probe.tsv').read_text()
    record['files']['CompactStar/Physics/BNV/BreakpointStepper.hpp']=sha(ROOT/'CompactStar/Physics/BNV/BreakpointStepper.hpp')
    assessment=json.loads((OUT/'linear-stability/assessment.json').read_text())
    assert assessment['probe_sha256']==record['files']['build/readiness-standalone/linear_probe']
    assert assessment['script_sha256']==sha(ROOT/'tests/bnv/readiness_stability.py')
    for p in (ROOT/'tests/bnv').glob('*readiness*'):record['files'][str(p.relative_to(ROOT))]=sha(p)
    rows=list(csv.DictReader((OUT/'linear-stability/linear.tsv').open(),delimiter='\t'))
    assert len({(r['tier'],r['axis'],r['sign'],r['policy']) for r in rows})==72
    earlier=list(csv.DictReader((OUT/'linear-stability-build1/linear.tsv').open(),delimiter='\t'))
    assert len(rows)==len(earlier)
    assert all({k:v for k,v in a.items() if k!='wall_s'}=={k:v for k,v in b.items() if k!='wall_s'} for a,b in zip(rows,earlier))
    record['fresh_relinked_linear_numerical_rows_identical']=True
    record['linear_rows']=len(rows);record['edge_rows']=len(list(csv.DictReader((OUT/'edge-probe.tsv').open(),delimiter='\t')));assert record['edge_rows']==14
    record['source_diff']=subprocess.check_output(['git','diff',record['accepted'],'--','CompactStar','cmake','CMakeLists.txt','dependencies','tests/baselines'],cwd=ROOT,text=True);assert record['source_diff']==''
    record['master']=subprocess.check_output(['git','rev-parse','master'],cwd=ROOT,text=True).strip();assert record['master']=='9a5c1eca40758e4c39f9f80db7324e12ad23e505'
    (OUT/'verification.json').write_text(json.dumps(record,indent=2)+'\n')
    print('PASS',record['linear_rows'],'linear rows;',record['edge_rows'],'edge rows; production bytes unchanged')
if __name__=='__main__':main()
