"""Authority, algebra, refusal and historical preservation, without stellar work."""
import copy
import importlib.util
import json
from pathlib import Path
import unittest

REPO=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('completion',REPO/'tools/linux/results_first/cref_completion.py')
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)

class CompletionTests(unittest.TestCase):
    def test_joint_inequality_boundary(self):
        A,b,f,T=.2,.01,5.,1.
        cap=(T-A)*(1-b)/(2+f+T)
        self.assertAlmostEqual(c.partial_budget(A,b,f,cap),T)
        self.assertGreater(c.partial_budget(A,b,f,cap*1.01),T)
    def test_missing_or_bad_error_refused(self):
        for e in (-1,float('nan'),float('inf'),1):
            with self.assertRaises(ValueError):c.partial_budget(.1,.01,5.,e)
    def test_origin_and_absolute(self):
        self.assertEqual(c.partial_budget(.1,0,5,.02,absolute=True),.1+7*.02)
        self.assertEqual(c.partial_budget(.1,0,5,.02,absolute=True,origin=True),.1+5*.02)
    def test_owner_authority(self):
        r=json.loads((c.DOC/'preparation-registry.json').read_text());c.authenticate(r)
        for key in r['accepted_conditions']:
            bad=copy.deepcopy(r);del bad['accepted_conditions'][key]
            with self.assertRaises(ValueError):c.authenticate(bad)
        for key,value in [('scope','foreign'),('paired_trajectories_authorized',True),('approval_sha256','0'*64)]:
            bad=copy.deepcopy(r);bad[key]=value
            with self.assertRaises(ValueError):c.authenticate(bad)
    def test_failure_preserves_history_and_does_not_issue_certificate(self):
        old={'rows':[dict(index=1,quantity='envelope',historical={'threshold':1e-4,'drift_bound':1e-6},
            finite_added_estimates=1e-7,partial_utilization=.011,status='UNQUALIFIED')]}
        detail={'envelope':dict(signed_error=[[2.55e-6]],base_error_estimate=3e-6,fit_operator_norm=5.18)}
        out=c.assess(old,detail)
        self.assertEqual(out['rows'][0]['status'],'UNQUALIFIED')
        self.assertGreater(out['rows'][0]['restricted_preparation']['utilization_lower_bound'],1)
        self.assertFalse(out['restricted_summary']['operational_certificate_issued'])
        self.assertNotIn('restricted_preparation',old['rows'][0])

if __name__=='__main__':unittest.main()
