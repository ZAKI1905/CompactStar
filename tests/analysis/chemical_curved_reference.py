"""Independent expected mathematics; NEVER import an integration implementation here.

Interior Schwarzschild, R=12 km, M=9/5 km; C=1 and C=1+(r/R)^2.
The constant primitive follows r=sin(theta)/sqrt(k), and polynomial division
sin²(theta)/(a-cos(theta)) = a+cos(theta)+(1-a²)/(a-cos(theta)).
The variable reference uses arbitrary precision theta quadrature, not the
double radial rule under test and not any production Geometry helper.
"""
import mpmath as mp


def references(dps):
    with mp.workdps(dps):
        r, m = mp.mpf(12), mp.mpf(9)/5
        k = 2*m/r**3
        a = 3*mp.sqrt(1-2*m/r)
        end = mp.asin(mp.sqrt(k)*r)
        primitive = 8*mp.pi/k**mp.mpf('1.5')*(a*end+mp.sin(end)+(1-a*a)*
            2/mp.sqrt(a*a-1)*mp.atan(mp.sqrt((a+1)/(a-1))*mp.tan(end/2)))
        def theta_integrand(t, variable):
            shape = 1+(mp.sin(t)**2/(k*r*r) if variable else 0)
            return 8*mp.pi/k**mp.mpf('1.5')*mp.sin(t)**2/(a-mp.cos(t))*shape
        constant_check = mp.quad(lambda t:theta_integrand(t,False), [0,end])
        variable = mp.quad(lambda t:theta_integrand(t,True), [0,end])
        assert abs(primitive-constant_check) < mp.mpf(10)**(-dps+8)
        return [mp.nstr(primitive,dps),mp.nstr(variable,dps)]
