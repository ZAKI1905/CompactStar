#!/usr/bin/env python3
"""Generate separate instrumented X1 copies; never edit or rescore X1."""
import hashlib,json,pathlib,subprocess
SRC=pathlib.Path(__file__).resolve().parents[3];ROOT=SRC.parents[1];OUT=ROOT/'x2-fixed-step/build';TC=ROOT/'toolchain'

def replace(s,a,b):
 assert a in s,a
 return s.replace(a,b)
for variant,name in [('regularized','x1_regularized_oracle'),('unregularized','onset_structure_oracle')]:
 p=SRC/f'tests/bnv/{name}.cpp';s=p.read_text()
 s=replace(s,'#include "../eos/structure1/tov_oracle.hpp"',f'#include "{SRC}/tests/eos/structure1/tov_oracle.hpp"' + '\n#include <cstdlib>\n#include <vector>\nstatic double center_factor(){const char*p=getenv("X2_CENTER");return p?std::stod(p):1e-10;}\nstatic bool corrected(){return getenv("X2_CORRECT_CENTER");}\nstatic bool fixed(){return getenv("X2_FIXED_H");}')
 s=replace(s,'delta=Hc*1e-10','delta=Hc*center_factor()')
 if variant=='regularized':
  s=replace(s,'void add(double x) {','void add(double x) { if(getenv("X2_PLAIN")){sum+=x;return;}')
  s=replace(s,'unsigned accepted=0, rejected=0;','unsigned accepted=0, rejected=0;\n double passive[8]={};int calls=0;Kahan masssum;')
  s=replace(s,'if(jac==0) {','if(jac==0) { o.passive[o.calls++]=0;')
  s=replace(s,'f[5]=measure*q.ns[3]*jac;','f[5]=measure*q.ns[3]*jac;\n  if(o.calls>=8)throw std::runtime_error("passive stage overflow");\n  o.passive[o.calls++]=4*M_PI*q.e*o.to_geo*r*A*jac;')
  s=replace(s,'int rc=gsl_odeiv2_step_apply','calls=0;\n   int rc=gsl_odeiv2_step_apply')
  s=replace(s,'for(int j=0;j<6;++j)y[j]=proposed[j];','if(calls!=7)throw std::runtime_error("passive stage count");\n   masssum.add(h*(16./135*passive[0]+6656./12825*passive[2]+28561./56430*passive[3]-9./50*passive[4]+2./55*passive[5]));\n   for(int j=0;j<6;++j)y[j]=proposed[j];')
  s=replace(s,'Kahan sum[4];','if(corrected()){double v0=4*M_PI*ec/3,Q=v0+4*M_PI*pc,a=2/Q,cc=-4*M_PI*(ec+pc)*q.slope/5,bb=a/2*(-2*a*v0-(cc-4*M_PI*(ec+pc))/Q);y[0]+=bb*delta*delta;double cv=0;for(int j=0;j<4;j++){double k=eos.b*std::cbrt(3*M_PI*M_PI*q.ns[j]);cv+=k*std::hypot(k,eos.masses[j])/(3*std::pow(eos.b,3));}y[4]=vol*std::exp(Hc)*cv;}\n  masssum={};masssum.sum=y[1]*std::pow(y[0],1.5);\n  Kahan sum[4];')
  s=replace(s,'if(low==hn)rn=', 'std::cerr<<std::setprecision(17)<<"SEG "<<Hc<<" "<<low<<" "<<std::sqrt(y[0])<<" "<<y[1]*std::pow(y[0],1.5)<<" "<<y[2]<<" "<<y[3]<<" "<<y[4]<<" "<<masssum.sum<<" "<<accepted<<" "<<rejected<<"\\n";\n   if(low==hn)rn=')
  start=s.index('void checks()');s=s[:start]
  cls='Regularized'
 else:
  s=replace(s,'struct Oracle:structure1::EnthalpyOracle {',r"""struct Oracle:structure1::EnthalpyOracle {
 std::vector<double> passive;long double masssum=0;
 int advance(gsl_odeiv2_driver*d,double*x,double end,double*y){
  d->n=0;
  while(*x<end){
   auto q=at_H(Hc-*x);double r=std::sqrt(y[0]),A=(1-2*y[0]*y[1])/(y[1]+4*M_PI*q.p*to_geo);
   double first=4*M_PI*q.e*to_geo*r*A;passive.clear();
   int rc=gsl_odeiv2_evolve_apply(d->e,d->c,d->s,d->sys,x,end,&d->h,y);if(rc)return rc;
   if(passive.size()<6)throw std::runtime_error("passive stage count");size_t k=passive.size()-6;
   masssum+=d->e->last_step*(16./135*first+6656./12825*passive[k+1]+28561./56430*passive[k+2]-9./50*passive[k+3]+2./55*passive[k+4]);
   if(d->nmax>0&&d->n>d->nmax)return GSL_EMAXITER;
   if(d->h>d->hmax)d->h=d->hmax;if(d->h<d->hmin)return GSL_ENOPROG;d->n++;
  }return GSL_SUCCESS;
 }
""")
  s=replace(s,'f[5]=v*q.ns[3];','f[5]=v*q.ns[3];o.passive.push_back(4*M_PI*q.e*o.to_geo*r*A);')
  s=replace(s,'gsl_odeiv2_driver_apply(d,&x,Hc-h,y)','advance(d,&x,Hc-h,y)')
  s=replace(s,'gsl_odeiv2_system sys={rhs7','masssum=y[1]*std::pow(y[0],1.5);gsl_odeiv2_system sys={rhs7')
  s=replace(s,'gsl_odeiv2_system sys={rhs7' , 'if(corrected()){double v0=4*M_PI*ec/3,Q=v0+4*M_PI*pc,a=2/Q,cc=-4*M_PI*(ec+pc)*q.slope/5,bb=a/2*(-2*a*v0-(cc-4*M_PI*(ec+pc))/Q);y[0]+=bb*delta*delta;double cv=0;for(int j=0;j<4;j++){double k=eos.b*std::cbrt(3*M_PI*M_PI*q.ns[j]);cv+=k*std::hypot(k,eos.masses[j])/(3*std::pow(eos.b,3));}y[4]=vol*std::exp(Hc)*cv;}\n gsl_odeiv2_system sys={rhs7')
  s=replace(s,'if(h==hn)rn=', 'std::cerr<<std::setprecision(17)<<"SEG "<<Hc<<" "<<h<<" "<<std::sqrt(y[0])<<" "<<y[1]*std::pow(y[0],1.5)<<" "<<y[2]<<" "<<y[3]<<" "<<y[4]<<" "<<masssum<<"\\n";\n  if(h==hn)rn=')
  s=s[:s.index('int main(')];cls='Oracle'
 s+='''
int main(int argc,char**argv){try{
 gsl_set_error_handler_off();if(argc!=5)throw std::runtime_error("B rho_or_H P tol");
 CLASS o;double target=std::stod(argv[1]),hc=fixed()?std::stod(argv[2]):o.H_for_rho(std::stod(argv[2])),ps=std::stod(argv[3]),tol=std::stod(argv[4]);
 auto a=o.run(hc,ps,tol);double hp=hc,bp=a[0];
 if(!fixed()) { hc*=1.00001;for(int k=0;k<8;k++){a=o.run(hc,ps,tol);if(std::abs(a[0]/target-1)<ALIGN)break;double next=hc+(target-a[0])*(hc-hp)/(a[0]-bp);hp=hc;bp=a[0];hc=next;} }
 std::cout<<std::setprecision(17);for(double x:a)std::cout<<x<<' '; EXTRA std::cout<<'\\n';
 }catch(const std::exception&e){std::cerr<<e.what()<<'\\n';return 1;}}
'''.replace('CLASS',cls).replace('ALIGN','5e-13' if variant=='regularized' else 'std::max(5e-13,tol*.1)').replace('EXTRA','' if variant=='regularized' else "std::cout<<hc<<' ';")
 output=OUT/f'x1-{variant}.cpp';output.write_text(s)
 cmd=[str(TC/'bin/g++'),'-O3','-std=c++17','-ffp-contract=off','-isystem',str(TC/'install/Release/Zaki/include'),'-isystem',str(TC/'install/gsl-2.7.1/include'),str(output),'-o',str(OUT/f'x1-{variant}'),str(TC/'install/Release/Zaki/lib/libZaki.a'),str(TC/'install/gsl-2.7.1/lib/libgsl.a'),str(TC/'install/gsl-2.7.1/lib/libgslcblas.a')]
 subprocess.run(cmd,check=True)
