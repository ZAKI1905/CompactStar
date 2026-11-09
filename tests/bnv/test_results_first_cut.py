#!/usr/bin/env python3
"""Independent EOS, terminal and implicit-response algebra checks."""
import importlib.util
from pathlib import Path
import unittest
import mpmath as mp

ROOT=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('cut',ROOT/'tools/linux/results_first/cut_response.py')
cut=importlib.util.module_from_spec(spec)
spec.loader.exec_module(cut)


class CutResponse(unittest.TestCase):
    def setUp(self):
        mp.mp.dps=70

    def test_eos_against_independent_momentum_integrals(self):
        q=cut.eos(mp.mpf('0.000013073642058'))
        k,hc=q['k'],q['hc']
        energy=sum(mp.quad(lambda z:z*z*mp.sqrt(m*m+z*z),[0,k]) for m in q['masses'])/(mp.pi**2*hc**3)
        pressure=sum(mp.quad(lambda z:z**4/mp.sqrt(m*m+z*z),[0,k]) for m in q['masses'])/(3*mp.pi**2*hc**3)
        self.assertLess(abs(q['e']/energy-1),mp.mpf('1e-45'))
        self.assertLess(abs(q['p']/pressure-1),mp.mpf('1e-45'))

    def test_pressure_enthalpy_and_density_derivatives(self):
        H=mp.mpf('.000013')
        q=cut.eos(H)
        self.assertLess(abs(mp.diff(lambda h:cut.eos(h)['p'],H)/q['p_H']-1),mp.mpf('1e-45'))
        self.assertLess(abs(mp.diff(lambda h:cut.eos(h)['n'],H)/q['n_H']-1),mp.mpf('1e-45'))
        self.assertLess(abs((q['e']+q['p'])/(q['n']*q['chemical_enthalpy'])-1),mp.mpf('1e-45'))

    def test_terminal_redshift_identity(self):
        t=cut.terminal(mp.mpf('.000013'),mp.mpf('12.7'),mp.mpf('.92'),mp.mpf('-.078'))
        self.assertLess(abs(t['K_identity_residual']),mp.mpf('1e-60'))
        # The independently widened binary64 h0 and masses leave a real source
        # convention defect; neither set it to zero nor call it new cut error.
        self.assertNotEqual(t['thermodynamic_residual'],0)
        self.assertLess(t['R'],0)
        self.assertLess(t['B'],0)
        self.assertGreater(t['K'],0)

    def test_inverse_rank_one_response(self):
        g=mp.matrix([[2,mp.mpf('-.2')],[mp.mpf('-.2'),1]])
        k,l=mp.mpf('.03'),mp.mpf('-.04')
        def value(s,i,j):
            gs=g*mp.exp(-k*s)+mp.matrix([[l*s,0],[0,0]])
            return (gs**-1)[i,j]+mp.exp(k*s)/3
        z=g**-1+mp.ones(2)/3
        inv=g**-1
        for i in range(2):
            for j in range(2):
                self.assertLess(abs(mp.diff(lambda s:value(s,i,j),0)-(k*z[i,j]-l*inv[i,0]*inv[0,j])),mp.mpf('1e-60'))

    def test_fixed_B_mixed_tangent_by_implicit_root(self):
        H=mp.mpf('.7'); B0=mp.mpf('2'); a=mp.mpf('.4')
        b=lambda h:2+h+h*h
        u=lambda h:mp.exp(h)/10
        n=lambda h:mp.sin(h)
        for shell in (0,1):
            B=lambda h,s:b(h)+u(h)*s
            N=lambda h,s:n(h)+shell*u(h)*s
            bh=mp.diff(b,H); bs=u(H)
            out=cut.mixed_response(bh,bs,mp.diff(u,H),0,mp.diff(n,H),mp.diff(n,H,2),mp.diff(b,H,2),shell,B0,a)
            def path(s):
                h=mp.findroot(lambda h:B(h,s)-b(H),H)
                t=mp.diff(lambda v:N(v,s),h)/mp.diff(lambda v:B(v,s),h)
                return t+a/B0*(shell-t)*u(h)
            self.assertLess(abs(mp.diff(path,0)-out['path_s']),mp.mpf('1e-55'))
            def onlaw(x):
                s=a*x
                h=mp.findroot(lambda h:B(h,s)-(b(H)+B0*x),H)
                t=mp.diff(lambda v:N(v,s),h)/mp.diff(lambda v:B(v,s),h)
                return t+a/B0*(shell-t)*u(h),u(h)
            Tx=mp.diff(lambda x:onlaw(x)[0],0)
            Bsx=mp.diff(lambda x:onlaw(x)[1],0)
            alternate=cut.path_response_from_onlaw(out['path'],Tx,bs,Bsx,B0,a,shell)
            self.assertLess(abs(alternate-mp.diff(path,0)),mp.mpf('1e-55'))

    def test_wrong_order_and_zero_are_not_estimates(self):
        for values in ((1,1,1),(1,2,1),(1,2,3)):
            result=cut.diagnostic(list(map(mp.mpf,values)),4)
            self.assertEqual(result['status'],'INCONCLUSIVE')
            self.assertIsNone(result['formal_tail'])

    def test_unequal_stencil(self):
        c={'primary':[0]*8+['2']};lo={'primary':[0]*8+['1.8']};hi={'primary':[0]*8+['2.3']}
        d1,d2=cut.stencil(c,lo,hi)
        f=lambda x:mp.mpf(x)**2+3*mp.mpf(x)
        self.assertLess(abs(d1(f('2'),f('1.8'),f('2.3'))-7),mp.mpf('1e-60'))
        self.assertLess(abs(d2(f('2'),f('1.8'),f('2.3'))-2),mp.mpf('1e-60'))


if __name__=='__main__':
    unittest.main()
