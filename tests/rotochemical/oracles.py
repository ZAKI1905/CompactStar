"""Independent R1995 Fermi integrals and RE10b analytic/high precision oracles."""
import subprocess,sys
import mpmath as m
m.mp.dps=45
exe=sys.argv[1]
for line in subprocess.check_output([exe,'functions'],text=True).splitlines():
    x,fd,hd,fm,hm=map(m.mpf,line.split())
    if abs(x)>10 or abs(x)<m.mpf('.1'): continue
    for modified,F,H in [(False,fd,hd),(True,fm,hm)]:
        P=(lambda z:9*m.pi**4+10*m.pi**2*z*z+z**4) if modified else (lambda z:m.pi**2+z*z)
        def I(power,eta):
            return m.quad(lambda y:y**power*P(y-eta)/(1+m.exp(y-eta)),[0,10,30,100,m.inf])
        norm=2*I(3,0)
        ef=(I(3,x)+I(3,-x))/norm
        eh=(I(2,x)-I(2,-x))/norm
        assert abs(F/ef-1)<m.mpf('3e-14'),('F',modified,x,F,ef)
        assert abs(H/eh-1)<m.mpf('3e-14'),('H',modified,x,H,eh)
print('PASS RE3 independent Fermi convolutions')
mutants=0
for line in subprocess.check_output([exe,'integrals'],text=True).splitlines():
    if not line.startswith('INTEGRAL '): continue
    _,q,mutant,value=line.split();q=int(q);mutant=int(mutant);value=m.mpf(value)
    a=m.mpf('.4')-m.mpf('.1')*q
    bounds=[(m.mpf('.2'),m.mpf('.45')),(m.mpf('.7'),m.mpf('.9'))]
    # Independent entire exponential power-series antiderivative, not production GL.
    def primitive(x):
        return m.fsum(a**n/m.factorial(n)*(2*x**(2*n+3)/(2*n+3)+x**(2*n+5)/(2*n+5)) for n in range(100))
    integ=m.fsum(primitive(hi)-primitive(lo) for lo,hi in bounds)
    quad=m.fsum(m.quad(lambda x:x*x*(2+x*x)*m.exp(a*x*x),[lo,hi]) for lo,hi in bounds)
    assert abs(integ/quad-1)<m.mpf('1e-40')
    expected=4*m.pi*m.mpf('1e-25')*m.exp(m.mpf('.4')*(q-2))*integ
    residual=abs(value/expected-1)
    if mutant==0: assert residual<m.mpf('1e-10'),(q,value,expected,residual)
    else:
        assert residual>m.mpf('1e-10'),('surviving mutation',q,mutant,residual)
        mutants+=1
    print('RE10b',q,mutant,'relative_residual',float(residual))
print('PASS RE10b independent integrator oracle; transformed-input mutant cases',mutants,'(9 families, 2 q values)')
