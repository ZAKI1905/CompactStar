#!/usr/bin/env python3
import copy,importlib.util,os,pathlib,unittest
from fractions import Fraction as F
S=pathlib.Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('final_model',os.environ.get('BA13_FINAL_MODEL',S/'tools/linux/ba13_final_decision/model.py'))
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class DecisionTests(unittest.TestCase):
    def term(self,**kw):
        t=dict(name='a',value=0.1,evidence_type='conditional_estimate',condition_id='C-REF',authority=None,
               provenance=[dict(path='fixture',sha256='a'*64)],dependencies=[],missing=None)
        t.update(kw);return t
    def test_missing_error_is_not_zero(self):
        self.assertIsNone(m.complete_total([self.term(value=None)],{}))
    def test_authority_must_match_exactly(self):
        r={'accepted_conditions':{'C-REF':'owner:reviewed-commit'}}
        self.assertIsNone(m.complete_total([self.term(authority='some-other-review')],r))
    def test_general_standard_does_not_accept_condition(self):
        self.assertIsNone(m.complete_total([self.term(authority='owner')],{'accepted_standards':{'C-REF':'owner'}}))
    def test_missing_dependency(self):
        self.assertIsNone(m.complete_total([self.term(evidence_type='measured_difference',dependencies=['absent'])],{}))
    def test_independent_ratio_algebra(self):
        A,D,a,d=F(2,3),F(7,5),F(1,100),F(-1,200)
        t=A/D
        self.assertEqual(m.correlated_increment(a,d,t,D+d),(A+a)/(D+d)-t)
    def test_correlated_error_is_not_blanket_cancellation(self):
        self.assertEqual(m.correlated_increment(F(1,10),F(1,20),F(1,2),F(2)),F(3,80))
        self.assertEqual(m.correlated_increment(F(1,10),F(1,20),F(2),F(2)),0)
    def test_joint_ratio_bound(self):
        A,D,a,d=F(3),F(4),F(1,10),F(-1,5);t=A/D
        J=abs(a-t*d)
        self.assertEqual(m.joint_ratio_bound(J,abs(d),D),abs((A+a)/(D+d)-t))
        with self.assertRaises(ValueError):m.joint_ratio_bound(None,F(1),F(4))
        with self.assertRaises(ValueError):m.joint_ratio_bound(F(1),F(4),F(4))
    def test_feasibility_is_not_estimate(self):
        cap=m.feasible_cap(.001,.0002,.0001,5.18,.00001)
        e=cap['single_axis_pointwise_cap'];b=.00001
        # Independent finite-normalizer calculation, including existing finite numerator.
        result=((.0002+.0001)*(1-b)+(2+5.18)*e)/(1-b-e)
        self.assertAlmostEqual(result,.001,places=14)
        self.assertIsNone(cap['qualifying_value'])
    def test_tangents_remain_absolute(self):
        cap=m.feasible_cap(.0001,.00002,.00001,5.18,.01,absolute=True)
        self.assertAlmostEqual(cap['single_axis_pointwise_cap']*(2+5.18),.00007,places=18)
    def test_wrong_reference_refused(self):
        with self.assertRaises(ValueError):m.continuous_remainder(1,1,reference_kind='reconstruction',evidence_scope='analytic_supremum')
    def test_sampled_curvature_refused(self):
        for scope in ('three_stars','two_holdouts','smooth','numerically_observed'):
            with self.assertRaises(ValueError):m.continuous_remainder(1,1,reference_kind='X2',evidence_scope=scope)
        self.assertEqual(m.continuous_remainder(2,4,reference_kind='X2',evidence_scope='accepted_conditional_supremum'),4)
    def test_between_sample_counterexample(self):
        # C2 compact bump on one unobserved gap. Zero with first two derivatives
        # at both endpoints; arbitrary interior amplitude survives finite sampling.
        def bump(u,A):return A*64*u**3*(1-u)**3 if 0<=u<=1 else 0
        self.assertTrue(all(bump(u,100)==0 for u in (-1,0,1,2)))
        self.assertEqual(bump(.5,100),100)
        with self.assertRaises(ValueError):m.continuous_remainder(0,1,reference_kind='X2',evidence_scope='two_holdouts')
if __name__=='__main__':unittest.main()
