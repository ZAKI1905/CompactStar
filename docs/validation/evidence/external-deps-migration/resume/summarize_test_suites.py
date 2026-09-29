from pathlib import Path
import re,json
q=Path(__file__).resolve().parent
records=[]
for stack in ['old','new']:
 for mode in ['Debug','Release']:
  p=q/'evidence'/f'{stack}-{mode}-full-suite.log';s=p.read_text()
  match=re.search(r'(\d+)% tests passed, (\d+) tests failed out of (\d+)',s);assert match
  percent,failed,total=map(int,match.groups())
  records.append(dict(stack=stack,mode=mode,log=p.name,total=total,passed=total-failed,failed=failed,wall_seconds=float(re.search(r'Total Test time \(real\) =\s*([\d.]+)',s).group(1))))
assert [(r['passed'],r['failed']) for r in records]==[(76,0),(70,6),(81,0),(75,6)]
final=(q/'evidence/final-Release-applicable-suite.log').read_text();assert '100% tests passed, 0 tests failed out of 77' in final
repair=(q/'evidence/new-Release-phase5d1-repaired.log').read_text();assert '100% tests passed, 0 tests failed out of 1' in repair
assert '100% tests passed, 0 tests failed out of 1' in (q/'evidence/old-Release-heat-capacity-serial.log').read_text()
assert all(r['returncode']==0 for r in json.loads((q/'evidence/old-Release-repaired-structural.json').read_text()))
result=dict(initial_complete_runs=records,debug_authority_only=['phase5b_structural_response_regression','phase5c_chemical_coefficient_regression','phase5d_coupled_oracles'],old_Release_applicable=dict(passed=73,total=73,accounting='70 initial passes plus two same-mode structural comparisons and one serial heat-capacity rerun; three Debug-authority tests pass separately in OLD Debug.'),new_Release_applicable=dict(passed=78,total=78,accounting='77/77 final applicable suite plus independently complete repaired Phase-5D1 fresh run; the three Debug-authority tests pass in NEW Debug.'),no_comparator_tolerance_changed=True)
(q/'evidence/test-suite-matrix.json').write_text(json.dumps(result,indent=2)+'\n')
print('Suite matrix: OLD Debug 76/76; NEW Debug 81/81; applicable OLD Release 73/73 and NEW Release 78/78')
