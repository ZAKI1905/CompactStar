// Independent validation: analytic momentum-integral EOS and enthalpy TOV.
// Reuses authenticated independent local oracle, no production table/interpolator.
#include "../eos/structure1/tov_oracle.hpp"
#include <iostream>
#include <iomanip>
struct Oracle:structure1::EnthalpyOracle {
 static int rhs7(double x,const double*y,double*f,void*ptr){
  auto&o=*static_cast<Oracle*>(ptr);auto q=o.at_H(o.Hc-x);double r=std::sqrt(y[0]),metric=1-2*y[0]*y[1];
  double A=metric/(y[1]+4*M_PI*q.p*o.to_geo),v=4*M_PI*r*A/std::sqrt(metric);
  f[0]=2*A;f[1]=(4*M_PI*q.e*o.to_geo-3*y[1])*A/y[0];f[2]=v*q.n*.001;
  f[3]=q.ns[0]>0?v*std::exp(6*(o.Hc-x)):0;
  double slope=0;for(int i=0;i<4;++i){double k=o.eos.b*std::cbrt(3*M_PI*M_PI*q.ns[i]);slope+=k*std::hypot(k,o.eos.masses[i])/(3*std::pow(o.eos.b,3));}
  f[4]=v*std::exp(o.Hc-x)*slope;f[5]=v*q.ns[3];
  for(int j=0;j<6;++j)if(!std::isfinite(f[j]))return GSL_EBADFUNC;return GSL_SUCCESS;
 }
 std::array<double,8> run(double hc,double ps,double tol){
 Hc=hc;double hs=H_for_pressure(ps),delta=Hc*1e-10,x=delta;auto q=at_H(Hc);
 double ec=q.e*to_geo,pc=q.p*to_geo,r2=3*delta/(2*M_PI*(ec+3*pc)),r=std::sqrt(r2),vol=4*M_PI*r*r*r/3;
 double y[6]={r2,4*M_PI*ec/3-4*M_PI*(ec+pc)*q.slope*delta/5,vol*q.n*.001,vol*std::exp(6*Hc),0,vol*q.ns[3]};
 gsl_odeiv2_system sys={rhs7,nullptr,6,this};auto*d=gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_rkf45,delta,tol*1e-5,tol);
 double hn=std::log(eos.masses[0]/h0);double hp=std::sqrt(eos.masses[3]*eos.masses[3]-eos.masses[2]*eos.masses[2]);
 double hm=std::log((std::hypot(eos.masses[1],hp)+eos.masses[3])/h0);
 double ba=0,bb=.3;for(int k=0;k<70;++k){double mid=(ba+bb)/2;if(at_H(mid).e*1.602176634e33/(c*c)<1e10)ba=mid;else bb=mid;}
 double hb=(ba+bb)/2,rn=0,rb=0;
 for(double h:{hm,hb,hn,hs}){ // actual ordering hb > hn is checked below
  if(!(h<Hc&&h>=hs))continue;
  if(Hc-h<x)throw std::runtime_error("oracle boundary ordering");
  if(gsl_odeiv2_driver_apply(d,&x,Hc-h,y)!=GSL_SUCCESS)throw std::runtime_error("oracle integration failed");
  if(h==hn)rn=std::sqrt(y[0]);if(h==hb)rb=std::sqrt(y[0]);
 }
 gsl_odeiv2_driver_free(d);r=std::sqrt(y[0]);double m=y[1]*r*r*r,nu=.5*std::log(1-2*m/r),muinf=h0*std::exp(nu+hs);
 return {y[2]*1e57,r,m,rn,y[3]*std::exp(-6*(nu+hs))*1e-36,y[4]*std::exp(-(nu+hs))*1e15*std::pow(8.61733326214518e-11,2)*1.602176634e33,rb,muinf};
 }
};
int main(int argc,char**argv){try{
 gsl_set_error_handler_off();if(argc!=5)throw std::runtime_error("targetB rho_guess surfaceP_cgs tolerance");
 Oracle o;double target=std::stod(argv[1]),hc=o.H_for_rho(std::stod(argv[2])),ps=std::stod(argv[3]),tol=std::stod(argv[4]);
 auto a=o.run(hc,ps,tol);double hprev=hc,bprev=a[0];hc*=1.00001;std::array<double,8> b{};
 for(int k=0;k<8;++k){b=o.run(hc,ps,tol);if(std::abs(b[0]/target-1)<std::max(5e-13,tol*.1))break;
 double hn=hc+(target-b[0])*(hc-hprev)/(b[0]-bprev);hprev=hc;bprev=b[0];hc=hn;}
 if(std::abs(b[0]/target-1)>std::max(5e-13,tol*.1))throw std::runtime_error("oracle target alignment failed");
 std::cout<<std::setprecision(17);for(double v:b)std::cout<<v<<' ';std::cout<<hc<<'\n';
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
