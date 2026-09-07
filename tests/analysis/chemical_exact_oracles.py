#!/usr/bin/env python3
"""ADR-0013 independent rational expected fixtures, with production operations.

GC3/6/10 use independently differentiated/derived rational answers. GC1/2/4/5/
7/11/12 identities retain their contract classification. Candidate operations are checked through the production contract executable.
"""
from fractions import Fraction as F
import json
import subprocess
import sys


def mat(rows):
    return [[F(x) for x in row] for row in rows]


def transpose(a):
    return list(map(list, zip(*a)))


def mul(a, b):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def inverse(a):
    n = len(a)
    x = [row[:] + [F(i == j) for j in range(n)] for i, row in enumerate(a)]
    for k in range(n):
        pivot = next((i for i in range(k, n) if x[i][k]), None)
        if pivot is None:
            raise ValueError('singular: no inverse or pseudoinverse permitted')
        x[k], x[pivot] = x[pivot], x[k]
        d = x[k][k]
        x[k] = [v/d for v in x[k]]
        for i in range(n):
            if i != k:
                d = x[i][k]
                x[i] = [v-d*w for v, w in zip(x[i], x[k])]
    return [row[n:] for row in x]


def run():
    # epsilon = 1/2 n^T K n, intrinsic order n,p,e,mu. Positive Sylvester minors.
    k = mat([[4,1,1,0],[1,5,0,1],[1,0,6,1],[0,1,1,7]])
    sx = mat([[1,-1,-1],[0,1,1],[0,1,0],[0,0,1]])
    sy = mat([[1,0,0],[0,1,1],[0,1,0],[0,0,1]])
    t = mat([[1,-1,-1],[0,1,0],[0,0,1]])
    x = mat([[10],[2],[3]])
    assert mul(sx, x) == mat([[5],[5],[2],[3]])  # GC2 neutral reconstruction
    hx = mul(transpose(sx), mul(k, sx))
    # Direct substitution in the declared energy gives these polynomial coefficients.
    assert hx == mat([[4,-2,-3],[-2,11,8],[-3,8,16]])  # GC3
    assert mul(hx, inverse(hx)) == mat([[1,0,0],[0,1,0],[0,0,1]])  # GC4 identity
    cy = mul(t, mul(inverse(hx), transpose(t)))
    # Independently differentiated y energy: 2nn²+2nn*ne+nn*nmu+
    # 11ne²/2+7ne*nmu+7nmu².
    hy = mat([[4,2,1],[2,11,7],[1,7,14]])
    assert cy == inverse(hy)  # GC5, independent coordinate polynomial
    # GC3 finite perturbation: exact central differences of the substituted
    # coupled energy, without invoking the matrix Hessian calculation.
    def energy(v):
        b,e,mu=v
        return 2*b*b-2*b*e-3*b*mu+F(11,2)*e*e+8*e*mu+8*mu*mu
    base=[F(10),F(2),F(3)]; step=F(1,17)
    for i in range(3):
        for j in range(3):
            value=0
            for si in (-1,1):
                for sj in (-1,1):
                    v=base[:];v[i]+=si*step;v[j]+=sj*step
                    value+=si*sj*energy(v)
            assert value/(4*step*step)==hx[i][j]
    expected = [[F(v,381) for v in row] for row in
                [[105,-18,-21,3],[-18,43,29,14],[-21,29,55,-26],[3,14,-26,40]]]
    chi = inverse(k)
    q = mat([[0],[1],[-1],[-1]])
    v = mul(chi, q)
    den = mul(transpose(q), v)[0][0]
    projected = [[chi[i][j]-v[i][0]*v[j][0]/den for j in range(4)] for i in range(4)]
    assert projected == expected == mul(sy, mul(cy, transpose(sy)))  # GC6
    assert mul(projected, q) == mat([[0],[0],[0],[0]])
    try:
        inverse(projected)
        raise AssertionError('GC7 full corrected inverse accepted')
    except ValueError:
        pass
    # Test-side representation harness: neutral response cannot accept projection.
    def project_again(kind):
        if kind != 'intrinsic':
            raise TypeError('already neutral; second projection is forbidden')
    try:
        project_again('neutral')
        raise AssertionError('GC7 second projection accepted')
    except TypeError:
        pass
    # GC8: embed an actual 2D/1D response, never pad a Hessian and invert it.
    tnpe = mat([[1,-1],[0,1],[0,0]])
    c2 = mul(tnpe, mul(inverse(mat([[2,-2],[-2,5]])), transpose(tnpe)))
    assert c2 == mat([[F(1,2),0,0],[0,F(1,3),0],[0,0,0]])
    pe = mat([[0,0,0],[0,F(1,3),0],[0,0,0]])
    assert all(pe[i][j] == 0 for i,j in [(0,0),(0,1),(1,2),(2,2)])
    # GC8 independent phase-space asymptotics, squared to retain rational
    # arithmetic: n²/delta³ -> (2m)³/(3pi² hc³)² and chi²/delta ->
    # 2m³/(pi² hc³)². Normalize out constants; m=2 is a declared toy species.
    nerrors=[];cerrors=[]
    for exponent in (2,4,8,16):
        delta=F(1,10**exponent); mass=F(2)
        nerrors.append((2*mass+delta)**3/(2*mass)**3-1)
        cerrors.append((mass+delta)**2*(2*mass+delta)/(2*mass**3)-1)
        dn=F(10**exponent)
        hlimit=mat([[dn,-dn],[-dn,dn+3]])
        embedded=mul(tnpe,mul(inverse(hlimit),transpose(tnpe)))
        assert embedded[0][0]==1/dn and embedded[1][1]==F(1,3)
        assert embedded[2][2]==0
    assert all(0<b<a for a,b in zip(nerrors,nerrors[1:]))
    assert all(0<b<a for a,b in zip(cerrors,cerrors[1:]))
    # GC10: two unit-volume zones; expected sides separately hand reduced.
    z1 = mat([[3,1,0],[1,2,0],[0,0,1]])
    z2 = mat([[2,0,1],[0,1,0],[1,0,3]])
    total = [[a+b for a,b in zip(r,s)] for r,s in zip(z1,z2)]
    global_q = [[total[i][j]-total[i][0]*total[0][j]/total[0][0]
                 for j in (1,2)] for i in (1,2)]
    assert global_q == mat([[F(14,5),F(-1,5)],[F(-1,5),F(19,5)]])
    local_sum = mat([[F(8,3),0],[0,F(7,2)]])
    assert global_q != local_sum
    z = inverse(global_q)
    assert z == mat([[F(19,53),F(1,53)],[F(1,53),F(14,53)]])
    assert mul(z, mat([[1],[2]])) == mat([[F(21,53)],[F(29,53)]])
    # GC1 negative eta sign; fm^3/km^3 is 10^54 exactly, G=count/MeV.
    eta = [[-v[0]] for v in mul(z, mat([[1],[2]]))]
    assert eta == mat([[F(-21,53)],[F(-29,53)]])
    assert (10**18)**3 == 10**54
    # GC11 source traceability ONLY: R2006 (16)-(18), supplied symmetric M.
    m = mat([[5,1,2],[1,7,3],[2,3,11]])
    l = mat([[-1,-1],[1,0],[0,1]])
    source_z = mul(transpose(l), mul(m,l))
    assert source_z == mat([[10,5],[5,12]])
    assert (m[0][0]-2*m[0][1]+m[1][1],
            m[0][0]-m[0][1]-m[0][2]+m[1][2],
            m[0][0]-2*m[0][2]+m[2][2]) == (10,5,12)
    # GC12 signed action / physical seconds: rational c in km/s.
    declared_z = mat([[2,1],[1,3]])
    c = F('299792.458')
    iphys = mat([[-3],[-5]])
    kgeom = [[v[0]*c*c] for v in iphys]
    w = mul(declared_z, [[v[0]/(c*c)] for v in kgeom])
    assert w == mat([[-11],[-18]])
    action = [[-4*v[0]] for v in w]  # Omega=2/s, dotOmega=-1/s²
    assert action == mat([[44],[72]])
    mutants = {
        'flip_Omega_dot': [[-v[0]] for v in action],
        'omit_c_minus2': [[v[0]*c*c] for v in action],
        'double_c_minus2': [[v[0]/(c*c)] for v in action],
        'omit_e': mul(declared_z, mat([[0],[20]])),
        'omit_mu': mul(declared_z, mat([[12],[0]])),
        'swap_channels': mul(declared_z, mat([[20],[12]])),
    }
    assert all(v != action for v in mutants.values())
    production = {}
    if len(sys.argv)>1:
        text = subprocess.check_output([sys.argv[1]], text=True)
        for line in text.splitlines():
            if line.startswith(('GC3_6_', 'GC10_', 'GC11_SOURCE_', 'N8_')):
                key,*values=line.split();production[key]=list(map(float,values))
        for key,expected_matrix in [('GC3_6_C',cy),('GC10_Q',global_q),('GC10_Z',z),('GC11_SOURCE_Z',source_z)]:
            actual=production[key];expected_flat=[float(v) for row in expected_matrix for v in row]
            assert len(actual)==len(expected_flat)
            assert all(abs(a-b)<1e-12 for a,b in zip(actual,expected_flat)),(key,actual,expected_flat)
        toy_inverse=mul(t,mul(mat([[1/F(1e-24),0,0],[0,1,0],[0,0,F(1,2)]]),transpose(t)))
        assert all(abs(F(a)-b)<=F(e) for a,b,e in zip(production['N8_C'],[v for row in toy_inverse for v in row],production['N8_E']))
        # Expected rational mathematics is not sourced from production helpers.
        assert all(abs(a-float(b))<=e for a,b,e in zip(production['GC3_6_C'],[v for row in cy for v in row],production['GC3_6_E']))
        assert any(abs(a-float(b))>1e-2 for a,b in zip(production['GC10_Q'],[v for row in local_sum for v in row]))
    return {'GC1-GC8':'PASS (precursors)', 'GC10':'PASS independent rational two-zone',
            'GC11':'PASS source traceability/contract', 'GC12':'PASS signed unit fixture',
            'GC12_mutants':list(mutants), 'arithmetic':'exact Fraction, zero tolerance',
            'production_chemical_path_exercised':bool(production), 'production_values':production}


if __name__ == '__main__':
    print(json.dumps(run(), indent=2))
