// Validation only. UNVERIFIED SCIENTIFIC CANDIDATE. No production/old ODE includes.
#include <Zaki/Physics/Constants.hpp> // authenticated physical constants only
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef X2_DOUBLE
using Real=double;
#else
using Real=long double;
#endif
using V2=std::array<Real,2>;
const Real pi=acos(Real(-1));
struct Sum { Real s=0,c=0; bool plain=false; void add(Real x){if(plain){s+=x;return;}Real a=x-c,b=s+a;c=(b-s)-a;s=b;} };
struct Rule {
 std::vector<Real> x,w;
 explicit Rule(int n):x(n),w(n){
  for(int i=0;i<(n+1)/2;i++){
   Real z=cos(pi*(i+Real(.75))/(n+Real(.5))),p=0,dp=0;
   for(int it=0;it<50;it++){
    Real p0=1,p1=z;
    for(int k=2;k<=n;k++){p=(Real(2*k-1)*z*p1-Real(k-1)*p0)/k;p0=p1;p1=p;}
    p=p1;dp=n*(z*p-p0)/(z*z-1);Real d=p/dp;z-=d;
    if(abs(d)<4*std::numeric_limits<Real>::epsilon())break;
    if(it==49)throw std::runtime_error("quadrature nodes");
   }
   // Re-evaluate derivative at the final root.
   Real p0=1,p1=z;for(int k=2;k<=n;k++){p=(Real(2*k-1)*z*p1-Real(k-1)*p0)/k;p0=p1;p1=p;}
   dp=n*(z*p1-p0)/(z*z-1);
   x[i]=(1-z)/2;x[n-1-i]=(1+z)/2;w[i]=w[n-1-i]=1/((1-z*z)*dp*dp);
  }
 }
};
struct Point {Real e=0,p=0,n=0,cv=0,eH=0,h=0;std::array<Real,4> ns{},k{},mu{};};
struct EOS {
 std::array<Real,4> m={Zaki::Physics::NEUTRON_M_MEV,Zaki::Physics::PROTON_M_MEV,Zaki::Physics::ELECTRON_M_MEV,Zaki::Physics::MUON_M_MEV};
 // Deliberately retain X1's binary64 physical unit constants, then widen.
 Real hc=double(1/Zaki::Physics::MEV_2_INV_FM),h0=double(Zaki::Physics::PROTON_M_MEV+Zaki::Physics::ELECTRON_M_MEV);
 Real light=double(2.99792458e10),G=double(6.673e-8);
 Real geo=double(1.602176634e33*6.673e-8/std::pow(2.99792458e10,4)*1e10);
 Rule rule;Real hn,hm;
 explicit EOS(int n):rule(n){hn=log(m[0]/h0);hm=log((sqrt(m[1]*m[1]+m[3]*m[3]-m[2]*m[2])+m[3])/h0);}
 Real dens(Real mu,Real mass)const{return mu<=mass?0:pow((mu-mass)*(mu+mass),Real(1.5))/(3*pi*pi*hc*hc*hc);}
 Point at(Real H)const{
  Point q;Real dh=h0*expm1(H),h=h0+dh;q.h=h;
  if(H<=0)return q;
  Real t=((h-m[1])*(h+m[1])+m[2]*m[2])/(2*h);
  if(h<m[0]){
   Real k2=dh*(h+h0)*(h-m[1]+m[2])*(h+m[1]-m[2])/(4*h*h);
   q.ns[1]=q.ns[2]=pow(k2,Real(1.5))/(3*pi*pi*hc*hc*hc);
  }else{
   if(t>m[3]){
    Real lo=m[3],hi=h-m[1];
    for(int j=0;j<80;j++){
     t=(lo+hi)/2;if(t==lo||t==hi)break;
     Real np=dens(t,m[2])+dens(t,m[3]),kp=hc*cbrt(3*pi*pi*np);
     if(sqrt(m[1]*m[1]+kp*kp)+t<h)lo=t;else hi=t;
    }t=(lo+hi)/2;
   }
   q.ns={dens(h,m[0]),dens(t,m[2])+dens(t,m[3]),dens(t,m[2]),dens(t,m[3])};
  }
  std::array<Real,4> susceptibility{};
  for(int i=0;i<4;i++){
   Real k=hc*cbrt(3*pi*pi*q.ns[i]);q.k[i]=k;q.mu[i]=sqrt(m[i]*m[i]+k*k);
   susceptibility[i]=q.mu[i]*k/(pi*pi*hc*hc*hc);
   q.cv+=k*q.mu[i]/(3*hc*hc*hc);
   if(k==0)continue;
   Sum energy,pressure;
   for(size_t j=0;j<rule.x.size();j++){
    Real z=rule.x[j],a=sqrt(m[i]*m[i]+z*z*k*k);
    energy.add(rule.w[j]*z*z*a);pressure.add(rule.w[j]*z*z*z*z/a);
   }
   q.e+=3*q.ns[i]*energy.s;q.p+=q.ns[i]*k*k*pressure.s;
  }
  q.n=q.ns[0]+q.ns[1];Real l=susceptibility[2]+susceptibility[3];
  Real chi=susceptibility[0]+susceptibility[1]*l/(susceptibility[1]+l);
  q.eH=h*h*chi;return q;
 }
 Real cut(Real value,bool pressure)const{
  Real a=0,b=Real(.3);for(int k=0;k<100;k++){
   Real h=(a+b)/2;if(h==a||h==b)break;auto q=at(h);
   Real v=pressure?q.p*Real(1.602176634e33):q.e*Real(1.602176634e33)/(light*light);
   if(v<value)a=h;else b=h;
  }return (a+b)/2;
 }
};
struct Map {
 int mode;Real hc,on; // 0: center z, 1: onset -s, 2: x ablation
 std::pair<Real,Real> at(Real t)const{
  if(mode==0)return {hc-t*t,2*t};
  if(mode==1)return {on+t*t,-2*t};
  return {hc-t,1};
 }
};
struct Evaluation {V2 f;std::array<Real,4> J;std::array<Real,6> g;};
// Optional passive validation observer. Null in the historical executable.
// It cannot modify stages, step selection, shooting, or primary accumulators.
using X2Observer=void (*)(const EOS&,const Point&,Real,const V2&,Real);
X2Observer x2_observer=nullptr;
Evaluation evaluate(const EOS&e,const Point&q,Real H,Real jac,const V2&y){
 Real u=y[0],v=y[1],D=1-2*u*v,Q=v+4*pi*q.p*e.geo;
 if(!(u>0&&D>0&&Q>0))throw std::runtime_error("structure domain");
 Real r=sqrt(u),A=D/Q,E=4*pi*q.e*e.geo-3*v,Au=-2*v/Q,Av=(-2*u*Q-D)/(Q*Q);
 Real W=4*pi*r*A/sqrt(D)*jac;
 Evaluation f;f.f={2*A*jac,E*A/u*jac};
 f.J={2*Au*jac,2*Av*jac,E*(Au/u-A/(u*u))*jac,(-3*A+E*Av)/u*jac};
 f.g={W*q.n*Real(.001),H>e.hn?W*exp(6*H):0,W*exp(H)*q.cv,W*q.ns[3],4*pi*q.e*e.geo*r*A*jac,jac};
 return f;
}
std::array<Real,4> linear(std::array<std::array<Real,5>,4>a){
 for(int j=0;j<4;j++){
  int k=j;for(int i=j+1;i<4;i++)if(abs(a[i][j])>abs(a[k][j]))k=i;
  std::swap(a[k],a[j]);if(abs(a[j][j])<Real(1e-35))throw std::runtime_error("singular Newton");
  Real d=a[j][j];for(int l=j;l<5;l++)a[j][l]/=d;
  for(int i=0;i<4;i++)if(i!=j){d=a[i][j];for(int l=j;l<5;l++)a[i][l]-=d*a[j][l];}
 }
 return {a[0][4],a[1][4],a[2][4],a[3][4]};
}
using X2StepObserver=void (*)(const EOS&,const Map&,Real,Real,const V2&,const V2*,const Evaluation*);
X2StepObserver x2_step_observer=nullptr;
struct Solver {
 EOS eos;Rule audit{4};int N;Real cf;bool raw,plain;Real Hc,Hs,hb;
 Real thermo=0,beta=0,charge=0,stage=0,mudrift=0,defect=0;
 std::array<Sum,9> sums{}; // B,JL,JC,Nmu,Mstage,nu,Maudit,Baudit,nuaudit
 std::vector<std::array<Real,7>> segments;
 Solver(int n,int eq,Real center,bool r,bool p):eos(eq),N(n),cf(center),raw(r),plain(p){for(auto&s:sums)s.plain=p;}
 void step(const Map&map,Real t,Real h,V2&y){
  const Real a= sqrt(Real(3))/6,c[2]={Real(.5)-a,Real(.5)+a};
  const Real A[2][2]={{Real(.25),Real(.25)-a},{Real(.25)+a,Real(.25)}};
  Point q[2];Real H[2],jac[2];V2 Y[2];Evaluation f[2];
  auto [ht,jt]=map.at(t);auto ft=evaluate(eos,eos.at(ht),ht,jt,y);
  for(int i=0;i<2;i++){auto z=map.at(t+c[i]*h);H[i]=z.first;jac[i]=z.second;q[i]=eos.at(H[i]);for(int k=0;k<2;k++)Y[i][k]=y[k]+c[i]*h*ft.f[k];}
  Real residual=0;
  for(int it=0;it<16;it++){
   for(int i=0;i<2;i++)f[i]=evaluate(eos,q[i],H[i],jac[i],Y[i]);
   std::array<std::array<Real,5>,4> mat{};residual=0;
   for(int i=0;i<2;i++)for(int k=0;k<2;k++){
    int row=2*i+k;Real R=Y[i][k]-y[k];for(int j=0;j<2;j++)R-=h*A[i][j]*f[j].f[k];
    Real scale=std::max(abs(Y[i][k]),Real(1e-20));residual=std::max(residual,abs(R)/scale);mat[row][4]=-R;
    for(int j=0;j<2;j++)for(int l=0;l<2;l++)mat[row][2*j+l]=(i==j&&k==l?1:0)-h*A[i][j]*f[j].J[k*2+l];
   }
   if(residual<32*std::numeric_limits<Real>::epsilon())break;
   auto d=linear(mat);for(int i=0;i<2;i++)for(int k=0;k<2;k++)Y[i][k]+=d[2*i+k];
   if(it==15)throw std::runtime_error("Newton convergence");
  }
  stage=std::max(stage,residual);if(residual>Real(1e-15))throw std::runtime_error("stage residual");
  for(int j=0;j<2;j++){
   thermo=std::max(thermo,abs((q[j].e+q[j].p)/(q[j].n*q[j].h)-1));
   charge=std::max(charge,abs((q[j].ns[1]-q[j].ns[2]-q[j].ns[3])/q[j].n));
   beta=std::max(beta,abs((q[j].mu[1]+q[j].mu[2])/q[j].h-1));
  }
  // Independent off-collocation quadrature on the degree-two structure polynomial.
  // Uses radial dr/dt from its derivative, not A/r from the TOV RHS.
  for(size_t j=0;j<audit.x.size();j++){
   Real z=audit.x[j],l0=(z-c[1])/(c[0]-c[1]),l1=(z-c[0])/(c[1]-c[0]);
   Real b0=(z*z/2-c[1]*z)/(c[0]-c[1]),b1=(z*z/2-c[0]*z)/(c[1]-c[0]);
   V2 yp,yd;for(int k=0;k<2;k++){yp[k]=y[k]+h*(b0*f[0].f[k]+b1*f[1].f[k]);yd[k]=l0*f[0].f[k]+l1*f[1].f[k];}
   auto [H0,jac0]=map.at(t+h*z);auto qp=eos.at(H0);Real r=sqrt(yp[0]),m=yp[1]*r*r*r,dr=yd[0]/(2*r),D=1-2*m/r;
   Real dm=4*pi*r*r*qp.e*eos.geo*dr,dB=4*pi*r*r*qp.n*Real(.001)*dr/sqrt(D);
   Real dnu=(m+4*pi*r*r*r*qp.p*eos.geo)/(r*r*D)*dr;
   sums[6].add(h*audit.w[j]*dm);sums[7].add(h*audit.w[j]*dB);sums[8].add(h*audit.w[j]*dnu);
   auto rhs=evaluate(eos,qp,H0,jac0,yp);
   defect=std::max(defect,abs(yd[0]-rhs.f[0])/(abs(rhs.f[0])+Real(1e-30)));
  }
  if(x2_step_observer)x2_step_observer(eos,map,t,h,y,Y,f);
  for(int k=0;k<6;k++)sums[k].add(h*(f[0].g[k]+f[1].g[k])/2);
  if(x2_observer)for(int j=0;j<2;j++){
   Real D=1-2*Y[j][0]*Y[j][1],Q=Y[j][1]+4*pi*q[j].p*eos.geo;
   Real volume=h/2*4*pi*sqrt(Y[j][0])*sqrt(D)/Q*jac[j];
   x2_observer(eos,q[j],H[j],Y[j],volume);
  }
  for(int k=0;k<2;k++)y[k]+=h*(f[0].f[k]+f[1].f[k])/2;
  auto [hend,jend]=map.at(t+h);(void)jend;
  mudrift=std::max(mudrift,abs(sums[8].s-(Hc-hend)));
 }
 void segment(Map map,Real start,Real end,V2&y){
  Real c0=sums[2].s,l0=sums[1].s;
  for(int i=0;i<N;i++)step(map,start+(end-start)*Real(i)/N,(end-start)/N,y);
  segments.push_back({map.at(end).first,sqrt(y[0]),y[1]*pow(y[0],Real(1.5)),sums[0].s,sums[2].s-c0,sums[1].s-l0,sums[8].s});
 }
 std::array<Real,15> run(Real hc,Real ps){
  Hc=hc;Hs=eos.cut(ps,true);hb=eos.cut(Real(1e10),false);auto q=eos.at(Hc);
  Real e=q.e*eos.geo,p=q.p*eos.geo,eH=q.eH*eos.geo;
  Real v0=4*pi*e/3,Q=v0+4*pi*p,a=2/Q,c=-4*pi*eH/5,b=a/2*(-2*a*v0-(c-4*pi*(e+p))/Q);
  Real x0=Hc*cf/(N*Real(N)),u=a*x0+b*x0*x0,r=sqrt(u),vol=4*pi*r*r*r/3;
  V2 y={u,v0+c*x0};for(auto&s:sums){s.s=s.c=0;}
  Real mc=y[1]*r*r*r;
  std::array<Real,9> initial={vol*q.n*Real(.001),vol*exp(6*Hc),vol*exp(Hc)*q.cv,vol*q.ns[3],mc,x0,mc,vol*q.n*Real(.001),x0};
  for(int i=0;i<9;i++)sums[i].s=initial[i];
  if(x2_observer)x2_observer(eos,q,Hc,y,vol);
  segments.clear();
  Real mid=(Hc+eos.hm)/2;
  if(!(Hc>mid&&mid>eos.hm&&eos.hm>hb&&hb>eos.hn&&eos.hn>Hs))throw std::runtime_error("event ordering");
  segment({0,Hc,0},sqrt(x0),sqrt(Hc-mid),y);
  auto onsetseg=[&](Real hi,Real lo,Real on){if(raw)segment({2,Hc,on},Hc-hi,Hc-lo,y);else segment({1,Hc,on},-sqrt(hi-on),-sqrt(std::max(Real(0),lo-on)),y);};
  onsetseg(mid,eos.hm,eos.hm);onsetseg(eos.hm,hb,eos.hn);Real rb=sqrt(y[0]);onsetseg(hb,eos.hn,eos.hn);Real rn=sqrt(y[0]);
  // Surface segment also uses sqrt(H), regardless of onset ablation.
  segment({1,Hc,0},-sqrt(eos.hn),-sqrt(Hs),y);
  r=sqrt(y[0]);Real m=y[1]*r*r*r,nu=log1p(-2*m/r)/2,K=nu+Hs;
  Real cvfactor=Real(1e15)*Real(double(8.61733326214518e-11))*Real(double(8.61733326214518e-11))*Real(double(1.602176634e33));
  return {sums[0].s*Real(1e57),r,m,rn,sums[1].s*exp(-6*K)*Real(1e-36),sums[2].s*exp(-K)*cvfactor,rb,eos.h0*exp(K),Hc,(sums[4].s/m-1),(sums[6].s/m-1),(sums[7].s/sums[0].s-1),mudrift,stage,defect};
 }
};
void checks(){
 EOS e(96);Real thermo=0,dp=0,moment=0;Rule r(8);
 for(int p=0;p<=5;p++){Sum v;for(size_t j=0;j<r.x.size();j++)v.add(r.w[j]*2*pow(r.x[j],p+1));moment=std::max(moment,abs(v.s*Real(p+2)/2-1));}
 for(Real H:{Real(.0002),Real(.001),Real(.02),Real(.08),Real(.13)}){
  auto q=e.at(H);Real d=H*Real(1e-5);thermo=std::max(thermo,abs((q.e+q.p)/(q.n*q.h)-1));dp=std::max(dp,abs((e.at(H+d).p-e.at(H-d).p)/(2*d*(q.e+q.p))-1));
 }
 // Independent series residual: u equation O(x^2), v equation O(x).
 auto q=e.at(Real(.13));Real ec=q.e*e.geo,pc=q.p*e.geo,eh=q.eH*e.geo;
 Real v0=4*pi*ec/3,Q=v0+4*pi*pc,a=2/Q,c=-4*pi*eh/5,b=a/2*(-2*a*v0-(c-4*pi*(ec+pc))/Q);
 Real prevu=0,prevv=0,pu=0,pv=0;
 for(Real x:{Real(1e-5),Real(5e-6),Real(2.5e-6)}){
  auto f=evaluate(e,e.at(Real(.13)-x),Real(.13)-x,1,{a*x+b*x*x,v0+c*x});
  Real ru=abs((a+2*b*x-f.f[0])/a),rv=abs((c-f.f[1])/c);
  if(prevu){pu=log(prevu/ru)/log(Real(2));pv=log(prevv/rv)/log(Real(2));}prevu=ru;prevv=rv;
 }
 // Scalar y'=y: exact GL4 stability function, independently from tableau stages.
 Real err[3];for(int j=0;j<3;j++){int n=8<<j;Real h=Real(1)/n,y=1,aa=sqrt(Real(3))/6;
  for(int k=0;k<n;k++){Real d=1-h/4,det=d*d-h*h*(Real(.25)-aa)*(Real(.25)+aa);Real Y0=y*(d+h*(Real(.25)-aa))/det,Y1=y*(d+h*(Real(.25)+aa))/det;y+=h*(Y0+Y1)/2;}err[j]=abs(y-exp(Real(1)));}
 Real order=log(err[1]/err[2])/log(Real(2));
 bool support=e.at(e.hn-Real(1e-10)).ns[0]==0&&e.at(e.hn+Real(1e-10)).ns[0]>0&&e.at(e.hm-Real(1e-10)).ns[3]==0&&e.at(e.hm+Real(1e-10)).ns[3]>0;
 bool pass=thermo<Real(1e-13)&&dp<Real(1e-6)&&moment<Real(1e-14)&&pu>Real(1.9)&&pu<Real(2.1)&&pv>Real(.9)&&pv<Real(1.1)&&order>Real(3.8)&&order<Real(4.2)&&support;
 std::cout<<std::setprecision(21)<<"{\"pass\":"<<(pass?"true":"false")<<",\"thermo\":"<<thermo<<",\"dP\":"<<dp<<",\"moments\":"<<moment<<",\"center_u_order\":"<<pu<<",\"center_v_order\":"<<pv<<",\"GL4_order\":"<<order<<",\"support\":"<<(support?"true":"false")<<",\"digits\":"<<std::numeric_limits<Real>::digits<<"}\n";
 if(!pass)throw std::runtime_error("analytic checks failed");
}
int main(int argc,char**argv){try{
 if(argc==2&&std::string(argv[1])=="--check"){checks();return 0;}
 if(argc!=10)throw std::runtime_error("B Hc_guess surfaceP N center_factor eos_order raw plain fixed_H");
 Real target=std::stold(argv[1]),hc=std::stold(argv[2]),ps=std::stold(argv[3]);int n=std::stoi(argv[4]);
 Solver s(n,std::stoi(argv[6]),std::stold(argv[5]),std::stoi(argv[7]),std::stoi(argv[8]));
 auto v=s.run(hc,ps);int shots=1;Real hprev=hc,bprev=v[0],sensitivity=0;
 if(!std::stoi(argv[9])){
  hc*=Real(1.00001);for(int k=0;k<9;k++){
   v=s.run(hc,ps);shots++;sensitivity=(v[0]-bprev)/(hc-hprev);
   if(abs(v[0]/target-1)<=Real(2e-16))break;
   Real next=hc+(target-v[0])/sensitivity;hprev=hc;bprev=v[0];hc=next;
  }
  if(abs(v[0]/target-1)>Real(2e-16))throw std::runtime_error("common B root");
 }
 std::cout<<std::setprecision(21)<<"{\"values\":[";for(size_t i=0;i<v.size();i++)std::cout<<(i?",":"")<<v[i];
 std::cout<<"],\"shots\":"<<shots<<",\"dB_dH\":"<<sensitivity<<",\"thermo\":"<<s.thermo<<",\"charge\":"<<s.charge<<",\"beta\":"<<s.beta<<",\"segments\":[";
 for(size_t i=0;i<s.segments.size();i++){std::cout<<(i?",[":"[");for(int j=0;j<7;j++)std::cout<<(j?",":"")<<s.segments[i][j];std::cout<<"]";}
 std::cout<<"]}\n";
 return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
