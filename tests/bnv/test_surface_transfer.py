#!/usr/bin/env python3
"""Independent terminal interval, finite-map and qualification refusal tests."""
import importlib.util
from pathlib import Path
import unittest
import mpmath as mp

ROOT=Path(__file__).resolve().parents[2]
s=importlib.util.spec_from_file_location('transfer',ROOT/'tools/linux/results_first/surface_transfer.py')
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)


class SurfaceTransfer(unittest.TestCase):
    def setUp(self):
        mp.mp.dps=90;m.iv.dps=90;m.EPS=mp.mpf('2.1e-12')

    def inside(self,value,interval):
        self.assertLessEqual(m.low(interval),value)
        self.assertLessEqual(value,m.high(interval))

    def test_interval_fermi_series_and_derivative(self):
        for hs in ('0.000012','0.0000130736','0.000014'):
            H=mp.mpf(hs);q=m.eos_jet(m.Jet(m.I(H),[m.I(1),m.I(0),m.I(0)]))
            ref=m.cut.eos(H)
            for name in ('e','p','n','cv','chi'):
                self.inside(ref[name],q[name].v)
                derivative=mp.diff(lambda h:m.cut.eos(h)[name],H)
                self.inside(derivative,q[name].d[0])

    def test_structural_jacobian_with_independent_differentiation(self):
        y=list(map(mp.mpf,('.000013','12.76','.92')));scales=[mp.mpf(1)]*3
        B0=mp.mpf('7.6e56');v=m.structural([m.I(x-mp.mpf('1e-25'),x+mp.mpf('1e-25')) for x in y],scales,B0)
        for name,key in (('H','alpha'),('R','R'),('M','M'),('b','B')):
            for j in range(3):
                val=mp.diff(lambda h,r,mm:m.cut.terminal(h,r,mm,0)[key],tuple(y),tuple(1 if k==j else 0 for k in range(3)))
                if name=='b':val/=B0
                self.inside(val,v[name].d[j])

    def test_finite_chemical_inverse_identity(self):
        G=mp.matrix([[2,mp.mpf('-.2')],[mp.mpf('-.2'),1]])
        U=G**-1;Z=U+mp.ones(2)/3
        for sign in (-1,1):
            ell=sign*mp.mpf('2e-12');dk=sign*mp.mpf('3e-15')
            actual=mp.exp(dk)*((G+mp.matrix([[ell,0],[0,0]]))**-1+mp.ones(2)/3)
            formula=mp.matrix(2)
            for i in range(2):
                for j in range(2):formula[i,j]=mp.exp(dk)*(Z[i,j]-ell*U[i,0]*U[0,j]/(1+ell*U[0,0]))
            self.assertLess(max(abs(actual[i,j]-formula[i,j]) for i in range(2) for j in range(2)),mp.mpf('1e-80'))

    def test_local_count_identity_and_nonlinear_outputs(self):
        inputs=m.cut.read(m.cut.DOC/'retained-inputs.json.gz');raw=inputs['centers']['0']['2048'];y=m.state(raw)
        z=m.flow(y,mp.mpf('2.1e-12'),2);out=m.output(raw,z)
        self.assertLess(abs((z[4]-y[4])-(z[10]-y[10])),mp.mpf('1e-30'))
        self.assertEqual(out['N_n'],mp.mpf(raw['N'][0]))
        self.assertEqual(out['N_mu'],mp.mpf(raw['N'][3]))
        self.assertNotEqual(out['Tsurface_inf'],m.cut.fields(raw)[0]['Tsurface_inf'])

    def test_endcap_unknown_survives_tiny_allowance_and_large_headroom(self):
        old={'rows':[dict(quantity='radius',index=0,historical={'threshold':'1'},
            static_completion={'known_numerator':'.01','known_baseline':'.001'})]}
        b={'allowances':{'radius':{'interior_point_bound':mp.mpf('1e-99')}}}
        result=m.successor(old,b,{'radius':{'fit_operator_norm':'5'}},{})
        row=result['rows'][0]['surface_transfer']
        self.assertIsNone(row['full_budget']);self.assertIsNone(row['baseline_cut_bound'])
        self.assertFalse(row['endpoint_extension_scenario']['authorized'])
        self.assertEqual(result['rows'][0]['static_completion'],old['rows'][0]['static_completion'])

    def test_series_refuses_outside_convergence_domain(self):
        with self.assertRaises(ValueError):m.series(m.Jet(m.I('1','1.1')),mp.mpf('.5'),3)


if __name__=='__main__':unittest.main()
