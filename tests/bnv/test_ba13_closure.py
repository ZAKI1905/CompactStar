"""Independent fault controls for closure qualification and paired quantities."""
import copy
import pathlib
import sys
import unittest
import numpy as np

sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[2]/'tools/linux/ba13_headroom'))
from enclosure import FRACTIONS,L,inverse_majorant,qualify,root_error

class ClosureTests(unittest.TestCase):
    def setUp(self):
        self.history=dict(drift_bound=2e-5,threshold=1e-4,utilization=.2)
        self.terms={k:dict(value=0.,evidence_type='analytic_bound') for k in FRACTIONS}
    def test_conditional_without_missing_flag_refuses(self):
        self.terms['integration']=dict(value=0.,evidence_type='conditional_estimate',
                                      condition_id='unproved-remainder',missing_estimate='')
        q=qualify(self.history,self.terms)
        self.assertEqual(q['status'],'UNQUALIFIED')
        self.assertIsNone(q['proposed_total_utilization'])
    def test_forged_authority_cannot_accept_itself(self):
        self.terms['integration']=dict(value=1e-10,evidence_type='conditional_estimate',
                                      condition_id='foo',authority='invented owner acceptance')
        self.assertEqual(qualify(self.history,self.terms)['status'],'UNQUALIFIED')
    def test_explicit_acceptance_and_remaining_gap(self):
        self.terms['integration']=dict(value=1e-10,evidence_type='conditional_estimate',
                                      condition_id='test-only',authority='test-fixture')
        accepted={'test-only':'test-fixture'}
        self.assertEqual(qualify(self.history,self.terms,accepted_authorities=accepted)['status'],'PASS')
        self.terms['center']['missing_estimate']='unbounded center remainder'
        self.assertEqual(qualify(self.history,self.terms,accepted_authorities=accepted)['status'],'UNQUALIFIED')
    def test_numerical_standard_requires_acceptance(self):
        self.terms['integration']=dict(value=1e-10,evidence_type='validated_numerical_estimate',
                                      standard_id='observed-order',authority='review-proposal')
        self.assertEqual(qualify(self.history,self.terms)['status'],'UNQUALIFIED')
    def test_each_missing_contribution_refuses(self):
        for k in FRACTIONS:
            t=copy.deepcopy(self.terms);t[k]['value']=None
            self.assertEqual(qualify(self.history,t)['status'],'UNQUALIFIED',k)
        self.terms.pop('common_B')
        with self.assertRaises(ValueError):qualify(self.history,self.terms)
    def test_nonfinite_and_negative_refuse(self):
        for v in (float('nan'),float('inf'),-1.):
            t=copy.deepcopy(self.terms);t['center']['value']=v
            with self.assertRaises(ValueError):qualify(self.history,t)
    def test_exact_does_not_use_tolerance(self):
        for drift in (1e-30,2e-5):
            h=dict(self.history,drift_bound=drift)
            self.assertEqual(qualify(h,self.terms,exact=True)['status'],'FAIL')
        self.assertEqual(qualify(dict(self.history,drift_bound=0),self.terms,exact=True)['status'],'PASS')
    def test_threshold_mutation_refuses(self):
        h=dict(self.history,threshold=.1)
        with self.assertRaises(ValueError):qualify(h,self.terms,frozen_historical=self.history)
    def test_wrong_reference_counterexample(self):
        # X2 self-error can vanish while the fixed-production-profile structure
        # is wrong. P-R alone is zero and cannot describe P-Q.
        p=np.array([1.,1.00004]);r=p.copy();x=np.array([1.,1.00001])
        rep=(p[1]-r[1])-(p[0]-r[0])
        structural=(r[1]-x[1])-(r[0]-x[0])
        paired=(p[1]-x[1])-(p[0]-x[0])
        self.assertEqual(rep,0.)
        self.assertAlmostEqual(rep+structural,paired)
        self.assertGreater(abs(paired),2.9e-5)
    def test_common_B_tangent_absolute_units(self):
        sensitivity=.127;coordinate=2.32e-12;t0=.0033776
        absolute=sensitivity*coordinate
        self.assertAlmostEqual(absolute,2.9464e-13)
        self.assertGreater(absolute/t0,absolute)
        self.assertGreater(root_error(0.,2.5e-15,2.),0.)
    def test_independent_inverse_both_nonlinear_tails(self):
        g0=np.diag([1.,2.,3.]);gi=np.diag([1.2,2.1,2.8]);e=.01*np.eye(3)
        a0,t0,r0,_=inverse_majorant(g0,e);ai,ti,ri,_=inverse_majorant(gi,e)
        bound=abs(L).T@(t0+ti)@abs(L)
        actual=abs(L.T@(np.linalg.inv(gi-e)-ai-np.linalg.inv(g0+e)+a0)@L)
        self.assertTrue(np.all(actual<=bound))
        self.assertTrue(np.any((abs(L).T@(r0+ri)@abs(L))>0))

if __name__=='__main__':unittest.main(verbosity=2)
