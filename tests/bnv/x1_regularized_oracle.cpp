// VALIDATION ONLY / UNVERIFIED SCIENTIFIC CANDIDATE. No production linkage.
// Equations and numerical limitations: docs/validation/x1-oracle-anchor/MATHEMATICS.md.
#include "../eos/structure1/tov_oracle.hpp"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>

struct Kahan {
 double sum=0, correction=0;
 void add(double x) { double a=x-correction, b=sum+a; correction=(b-sum)-a; sum=b; }
};

struct Regularized : structure1::EnthalpyOracle {
 bool transformed=false;
 double onset=0, hn=0, hm=0;
 double max_thermo=0, max_charge=0, max_beta=0, max_measure=0, max_mass=0;
 unsigned accepted=0, rejected=0;
 Regularized() {
  hn=std::log(eos.masses[0]/h0);
  double p=std::sqrt(eos.masses[3]*eos.masses[3]-eos.masses[2]*eos.masses[2]);
  hm=std::log((std::hypot(eos.masses[1],p)+eos.masses[3])/h0);
 }
 // t=-s increases outward; x=Hc-H, dx/dt=-2t >=0.
 static int rhs(double t,const double*y,double*f,void*ptr) {
  auto&o=*static_cast<Regularized*>(ptr);
  double H=o.transformed ? o.onset+t*t : o.Hc-t;
  double jac=o.transformed ? -2*t : 1;
  if(jac==0) { std::fill(f,f+6,0.); return GSL_SUCCESS; }
  auto q=o.at_H(H); double u=y[0],v=y[1],r=std::sqrt(u),metric=1-2*u*v;
  if(!(u>0 && metric>0 && jac>0))return GSL_EBADFUNC;
  double A=metric/(v+4*M_PI*q.p*o.to_geo);
  double measure=4*M_PI*r*A/std::sqrt(metric);
  f[0]=2*A*jac;
  f[1]=(4*M_PI*q.e*o.to_geo-3*v)*A/u*jac;
  f[2]=measure*q.n*.001*jac;
  f[3]=H>o.hn ? measure*std::exp(6*H)*jac : 0;
  Kahan slope;
  for(int i=0;i<4;++i) {
   double k=o.eos.b*std::cbrt(3*M_PI*M_PI*q.ns[i]);
   slope.add(k*std::hypot(k,o.eos.masses[i])/(3*std::pow(o.eos.b,3)));
  }
  f[4]=measure*std::exp(H)*slope.sum*jac;
  f[5]=measure*q.ns[3]*jac;
  for(int j=0;j<6;++j)if(!std::isfinite(f[j]))return GSL_EBADFUNC;
  return GSL_SUCCESS;
 }
 void check_state(double H,const double*y) {
  auto q=at_H(H);
  max_thermo=std::max(max_thermo,std::abs((q.e+q.p)/(q.n*h0*std::exp(H))-1));
  max_charge=std::max(max_charge,std::abs((q.ns[1]-q.ns[2]-q.ns[3])/q.n));
  if(q.ns[0]>0) max_beta=std::max(max_beta,std::abs((eos.potential(q.ns[1],eos.masses[1])+eos.potential(q.ns[2],eos.masses[2]))/(h0*std::exp(H))-1));
  double r=std::sqrt(y[0]),m=y[1]*r*r*r,metric=1-2*m/r;
  double A=metric/(y[1]+4*M_PI*q.p*to_geo),drdx=A/r;
  double via_radial=4*M_PI*r*r*drdx/std::sqrt(metric);
  double via_H=4*M_PI*r*A/std::sqrt(metric);
  max_measure=std::max(max_measure,std::abs(via_radial/via_H-1));
  // Chain rule dm/dx = r^3 dv/dx + 3vr^2 dr/dx =4pi r^2 epsilon dr/dx.
  double dvdx=(4*M_PI*q.e*to_geo-3*y[1])*A/y[0];
  double dm=r*r*r*dvdx+3*y[1]*r*r*drdx,ref=4*M_PI*r*r*q.e*to_geo*drdx;
  // Absolute cancellation-normalized residual, not divided by tiny surface density.
  max_mass=std::max(max_mass,std::abs(dm-ref)/(std::abs(r*r*r*dvdx)+std::abs(3*y[1]*r*r*drdx)+std::abs(ref)));
 }
 void segment(double high,double low,bool reg,double on,double tol,double*y,Kahan*sum) {
  transformed=reg;onset=on;
  double t=reg ? -std::sqrt(std::max(0.,high-on)) : Hc-high;
  double end=reg ? -std::sqrt(std::max(0.,low-on)) : Hc-low;
  double step=std::min((end-t)/16,1e-6);
  gsl_odeiv2_system sys={rhs,nullptr,6,this};
  auto*s=gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkf45,6);
  auto*c0=gsl_odeiv2_control_y_new(tol*1e-5,tol);
  while(t<end) {
   if(accepted+rejected>2000000)throw std::runtime_error("step ceiling");
   double h=std::min(step,end-t),trial[6]={y[0],y[1],0,0,0,0},err[6],yp[6],proposed[6];
   if(t+h==t)throw std::runtime_error("step underflow");
   // Integrals do not enter any RHS. Resetting their local base to zero makes
   // trial[2:6] accurate increments, avoiding subtraction of two large totals.
   int rc=gsl_odeiv2_step_apply(s,t,h,trial,err,nullptr,yp,&sys);
   if(rc!=GSL_SUCCESS)throw std::runtime_error("RKF45 failure");
   proposed[0]=trial[0];proposed[1]=trial[1];
   Kahan next[4];for(int j=0;j<4;++j){next[j]=sum[j];next[j].add(trial[j+2]);proposed[j+2]=next[j].sum;}
   double next_h=h;
   int adjust=gsl_odeiv2_control_hadjust(c0,s,proposed,err,yp,&next_h);
   if(adjust==GSL_ODEIV_HADJ_DEC) {++rejected;step=next_h;gsl_odeiv2_step_reset(s);continue;}
   for(int j=0;j<6;++j)y[j]=proposed[j];
   for(int j=0;j<4;++j)sum[j]=next[j];
   t=(h==end-t)?end:t+h;step=next_h;++accepted;
   check_state(reg ? on+t*t : Hc-t,y);
  }
  gsl_odeiv2_control_free(c0);gsl_odeiv2_step_free(s);
 }
 std::array<double,9> run(double hc,double ps,double tol) {
  Hc=hc;double hs=H_for_pressure(ps),delta=Hc*1e-10;auto q=at_H(Hc);
  double ec=q.e*to_geo,pc=q.p*to_geo,u=3*delta/(2*M_PI*(ec+3*pc));
  double r=std::sqrt(u),vol=4*M_PI*r*r*r/3;
  double y[6]={u,4*M_PI*ec/3-4*M_PI*(ec+pc)*q.slope*delta/5,vol*q.n*.001,vol*std::exp(6*Hc),0,vol*q.ns[3]};
  // Preserve the established center prescription, including zero initial C.
  Kahan sum[4];for(int j=0;j<4;++j)sum[j].sum=y[j+2];
  double a=0,b=.3;for(int k=0;k<70;++k){double h=(a+b)/2;if(at_H(h).e*1.602176634e33/(c*c)<1e10)a=h;else b=h;}
  double hb=(a+b)/2,rn=0,rb=0,high=Hc-delta;
  if(!(Hc>hm && hm>hb && hb>hn && hn>hs))throw std::runtime_error("onset ordering");
  for(double low:{hm,hb,hn,hs}) {
   bool reg=low>=hn;double on=(low==hm)?hm:hn;
   segment(high,low,reg,on,tol,y,sum);
   if(low==hn)rn=std::sqrt(y[0]);if(low==hb)rb=std::sqrt(y[0]);high=low;
  }
  r=std::sqrt(y[0]);double m=y[1]*r*r*r,nu=.5*std::log(1-2*m/r),muinf=h0*std::exp(nu+hs);
  return {y[2]*1e57,r,m,rn,y[3]*std::exp(-6*(nu+hs))*1e-36,y[4]*std::exp(-(nu+hs))*1e15*std::pow(8.61733326214518e-11,2)*1.602176634e33,rb,muinf,hc};
 }
};

void checks() {
 Regularized o;
 // Independent Gaussian moments: integrate (H-Hon)^(p/2) dH via s,
 // exact result 2*S^(p+2)/(p+2), including nonsmooth p=1 heat-capacity onset.
 auto*g=gsl_integration_glfixed_table_alloc(16);double maxmoment=0;
 for(int p=0;p<=5;++p){double val=0;for(int j=0;j<16;++j){double s,w;gsl_integration_glfixed_point(0,.2,j,&s,&w,g);val+=w*2*s*std::pow(s,p);}double exact=2*std::pow(.2,p+2)/(p+2);maxmoment=std::max(maxmoment,std::abs(val/exact-1));}
 gsl_integration_glfixed_table_free(g);
 double maxdp=0,maxmu=0,endpoint=0,transform=0;
 for(double H:{.0002,.001,.02,.08,.13}) {
  double d=1e-7*H;auto q=o.at_H(H);double diff=(o.at_H(H+d).p-o.at_H(H-d).p)/(2*d);
  maxdp=std::max(maxdp,std::abs(diff/(q.e+q.p)-1));
  double nu=.1-H;maxmu=std::max(maxmu,std::abs(q.h*std::exp(nu)/(o.h0*std::exp(.1))-1));
 }
 // Compare implemented transformed RHS to an independent x-coordinate RHS
 // evaluation at an arbitrary physical state, for both onset coordinates.
 double y[6]={100,.001,0,0,0,0};o.Hc=.14;
 for(double on:{o.hn,o.hm}) {
  double s=.01,H=on+s*s,fx[6],fs[6];o.transformed=false;Regularized::rhs(o.Hc-H,y,fx,&o);
  o.transformed=true;o.onset=on;Regularized::rhs(-s,y,fs,&o);
  for(int j=0;j<6;++j)transform=std::max(transform,std::abs(fs[j]-2*s*fx[j])/(std::abs(fs[j])+1e-100));
  Regularized::rhs(0,y,fs,&o);for(double v:fs)endpoint=std::max(endpoint,std::abs(v));
 }
 auto belowN=o.at_H(o.hn-1e-10),aboveN=o.at_H(o.hn+1e-10),belowM=o.at_H(o.hm-1e-10),aboveM=o.at_H(o.hm+1e-10);
 bool support=belowN.ns[0]==0 && aboveN.ns[0]>0 && belowM.ns[3]==0 && aboveM.ns[3]>0;
 Kahan k;k.add(1);for(int j=0;j<10000;++j)k.add(1e-16);double sumerr=std::abs(k.sum-(1+1e-12));
 bool pass=maxmoment<1e-13 && maxdp<1e-5 && maxmu<1e-13 && transform<1e-9 && endpoint==0 && support && sumerr<3e-16;
 std::cout<<std::setprecision(17)<<"{\"pass\":"<<(pass?"true":"false")<<",\"moment_rel\":"<<maxmoment<<",\"dP_dH_rel\":"<<maxdp<<",\"mu_redshift_rel\":"<<maxmu<<",\"transform_rel\":"<<transform<<",\"endpoint\":"<<endpoint<<",\"species_support\":"<<(support?"true":"false")<<",\"kahan_abs\":"<<sumerr<<"}\n";
 if(!pass)throw std::runtime_error("independent checks failed");
}
int main(int argc,char**argv) {try {
 gsl_set_error_handler_off();if(argc==2 && std::string(argv[1])=="--check"){checks();return 0;}
 if(argc!=5)throw std::runtime_error("targetB rho_guess surfaceP_cgs tolerance");
 Regularized o;double target=std::stod(argv[1]),hc=o.H_for_rho(std::stod(argv[2])),ps=std::stod(argv[3]),tol=std::stod(argv[4]);
 auto a=o.run(hc,ps,tol);double hp=hc,bp=a[0];hc*=1.00001;std::array<double,9>b{};int shots=1;
 for(int k=0;k<8;++k){b=o.run(hc,ps,tol);++shots;if(std::abs(b[0]/target-1)<5e-13)break;double next=hc+(target-b[0])*(hc-hp)/(b[0]-bp);hp=hc;bp=b[0];hc=next;}
 if(std::abs(b[0]/target-1)>5e-13)throw std::runtime_error("oracle target alignment failed");
 std::cout<<std::setprecision(17);for(double v:b)std::cout<<v<<' ';std::cout<<'\n';
 std::cerr<<std::setprecision(17)<<"{\"hn\":"<<o.hn<<",\"hm\":"<<o.hm<<",\"shots\":"<<shots<<",\"accepted\":"<<o.accepted<<",\"rejected\":"<<o.rejected<<",\"thermo_rel\":"<<o.max_thermo<<",\"charge_rel\":"<<o.max_charge<<",\"beta_rel\":"<<o.max_beta<<",\"measure_rel\":"<<o.max_measure<<",\"mass_identity_residual\":"<<o.max_mass<<"}\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
