"""ADR-0022 semantics and fail-closed completion. Runs in the static Slurm job."""
import copy
import importlib.util
import math
from pathlib import Path
import unittest

REPO = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('static_completion', REPO/'tools/linux/results_first/static_completion.py')
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)


class StaticTests(unittest.TestCase):
    def test_original_pointwise_model_is_special_case(self):
        A, b, e, f = .2, .01, .02, 5.
        self.assertEqual(c.revised_partial(A, b, e, f*e), (A*(1-b)+(2+f)*e)/(1-b-e))
        self.assertEqual(c.revised_partial(A, b, e, f*e, absolute=True), A+(2+f)*e)

    def test_origin_still_charges_fit_and_baseline(self):
        self.assertEqual(c.revised_partial(.1, .01, .02, .003, origin=True), (.1*.99+.003)/.97)

    def test_affine_error_drift_not_cancelled(self):
        # The full pointwise amplitude survives as 2e even when its fit response is zero.
        self.assertEqual(c.revised_partial(0., 0., .1, 0.), .2/.9)
        self.assertEqual(c.revised_partial(0., 0., .1, 0., origin=True), 0.)

    def test_separation_and_nonfinite_refusals(self):
        for value in (None, -1., math.nan, math.inf):
            with self.assertRaises((ValueError, TypeError)):
                c.revised_partial(.1, .01, value, .02)
        with self.assertRaises(ValueError):
            c.revised_partial(.1, .6, .4, .01)

    def test_full_allowance_and_LMe_exception(self):
        d = {'signed_error': [[-3e-6], [2e-6]]}
        self.assertEqual(c.point_allowance('envelope', d), 10*3e-6)
        self.assertEqual(c.point_allowance('Ltilde_Me', d), 10*9.1e-12)

    def test_missing_cut_cannot_become_zero(self):
        terms = [dict(name='finite', value=.01), dict(name='cut', value=None)]
        baseline = [dict(name='known', value=.02), dict(name='cut', value=None)]
        out = c.finish_ledger(.1, terms, baseline)
        self.assertIsNone(out['full_budget'])
        self.assertIsNone(out['normalization'])
        self.assertEqual(out['missing'], ['cut', 'baseline:cut'])
        self.assertAlmostEqual(out['partial_budget'], .11/.98)

    def test_full_numerator_normalization(self):
        out = c.finish_ledger(.1, [dict(name='fit', value=.03), dict(name='drift', value=.02)],
                              [dict(name='base', value=.04)])
        self.assertAlmostEqual(out['full_budget'], .15/.96)
        self.assertAlmostEqual(out['normalization'], .15*.04/.96)

    def test_absolute_tangent_still_refuses_missing_cut(self):
        out = c.finish_ledger(.1, [dict(name='cut', value=None)], [], absolute=True)
        self.assertEqual(out['normalization'], 0.)
        self.assertIsNone(out['full_budget'])

    def test_authority_is_exact_and_narrow(self):
        r = c.read(c.DOC/'acceptance-registry.json')
        c.authenticate(r)
        for k, v in [('pointwise_multiplier', 1), ('worst_case_increment_multiplier', 0),
                     ('amended_groups', ['envelope']), ('new_stellar_calculations_authorized', True),
                     ('paired_trajectories_authorized', True), ('historical_certificate_promoted', True)]:
            bad = copy.deepcopy(r)
            bad[k] = v
            with self.assertRaises(ValueError):
                c.authenticate(bad)
        bad = copy.deepcopy(r)
        bad['files'][next(iter(bad['files']))] = '0'*64
        with self.assertRaises(ValueError):
            c.authenticate(bad)


if __name__ == '__main__':
    unittest.main()
