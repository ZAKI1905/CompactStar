// Passive uncertainty assessment on X2's existing GL4 path, not a third solver.
// UNVERIFIED SCIENTIFIC CANDIDATE. See docs/validation/ba13-headroom/METHOD.md.
#define main x2_historical_main
#include "x2_fixed_step.cpp"
#undef main
std::array<Sum,12> extra;
Rule momentum(32);
// HISTORICAL INCONSISTENT estimator retained for failure provenance only.
// Its first-order intra-step propagator omits a leading O(h^5) contribution.
// Retired from execution and qualification; see ba13-closure/PREDECLARATION.md.
V2 transported{},majorant{},stage_effect{};
Real defect_quadrature_spread=0;
void transport(const EOS&e,const Map&map,Real t,Real h,const V2&y,const V2*Y,const Evaluation*f){
 const Real a=sqrt(Real(3))/6,c[2]={Real(.5)-a,Real(.5)+a};
 const Real A[2][2]={{Real(.25),Real(.25)-a},{Real(.25)+a,Real(.25)}};
 std::array<std::array<Real,5>,4> mat{};
 for(int i=0;i<2;i++)for(int k=0;k<2;k++)for(int j=0;j<2;j++)for(int l=0;l<2;l++)
  mat[2*i+k][2*j+l]=(i==j&&k==l?1:0)-h*A[i][j]*f[j].J[2*k+l];
 Real T[2][2]{};
 for(int l=0;l<2;l++){
  auto mm=mat;for(int i=0;i<2;i++)mm[2*i+l][4]=1;
  auto d=linear(mm);
  for(int k=0;k<2;k++){T[k][l]=(k==l);for(int j=0;j<2;j++)for(int z=0;z<2;z++)T[k][l]+=h/2*f[j].J[2*k+z]*d[2*j+z];}
 }
 V2 force{},absolute{};
 for(int order:{4,8}){
  Rule rule(order);V2 current{},magnitude{};
  for(int j=0;j<order;j++){
   Real z=rule.x[j],l0=(z-c[1])/(c[0]-c[1]),l1=(z-c[0])/(c[1]-c[0]);
   Real b0=(z*z/2-c[1]*z)/(c[0]-c[1]),b1=(z*z/2-c[0]*z)/(c[1]-c[0]);V2 yp,yd;
   for(int k=0;k<2;k++){yp[k]=y[k]+h*(b0*f[0].f[k]+b1*f[1].f[k]);yd[k]=l0*f[0].f[k]+l1*f[1].f[k];}
   auto [H,jac]=map.at(t+h*z);auto rhs=evaluate(e,e.at(H),H,jac,yp);V2 residual;
   for(int k=0;k<2;k++)residual[k]=rhs.f[k]-yd[k];
   for(int k=0;k<2;k++){
    Real v=residual[k];for(int l=0;l<2;l++)v+=(1-z)*h*rhs.J[2*k+l]*residual[l];
    current[k]+=h*rule.w[j]*v;magnitude[k]+=h*rule.w[j]*abs(v);
   }
  }
  if(order==8)for(int k=0;k<2;k++)defect_quadrature_spread=std::max(defect_quadrature_spread,abs(current[k]-force[k])/std::max(abs(y[k]),Real(1e-30)));
  force=current;absolute=magnitude;
 }
 auto mm=mat;
 for(int i=0;i<2;i++)for(int k=0;k<2;k++){
  Real residual=Y[i][k]-y[k];for(int j=0;j<2;j++)residual-=h*A[i][j]*f[j].f[k];mm[2*i+k][4]=-residual;
 }
 auto ds=linear(mm);V2 se{};
 for(int k=0;k<2;k++)for(int j=0;j<2;j++)for(int l=0;l<2;l++)se[k]+=abs(h/2*f[j].J[2*k+l]*ds[2*j+l]);
 V2 nt=force,nm=absolute,ns=se;
 for(int k=0;k<2;k++)for(int l=0;l<2;l++){nt[k]+=T[k][l]*transported[l];nm[k]+=abs(T[k][l])*majorant[l];ns[k]+=abs(T[k][l])*stage_effect[l];}
 transported=nt;majorant=nm;stage_effect=ns;
}
void observe(const EOS&e,const Point&q,Real H,const V2&,Real volume){
 for(int i=0;i<4;i++)extra[i].add(volume*q.ns[i]);
 std::array<Real,4>s{};
 for(int i=0;i<4;i++)s[i]=q.mu[i]*q.k[i]/(pi*pi*e.hc*e.hc*e.hc);
 Real den=s[1]+s[2]+s[3];
 std::array<Real,4>g={s[0],s[2]*(s[1]+s[3])/den,-s[2]*s[3]/den,s[3]*(s[1]+s[2])/den};
 for(int i=0;i<4;i++)extra[4+i].add(volume*exp(H)*g[i]);
 if(H>e.hm)extra[8].add(volume*exp(6*H));
 Real measure=volume*exp(-H)*q.ns[0];
 extra[9].add(measure);extra[10].add(measure*exp(-H)*q.mu[0]);
 // Historical P1 constants are deliberately retained, including its hbar*c.
 Real pf=Real(197.3269804)*cbrt(3*pi*pi*q.ns[0]);Sum avg;
 for(size_t j=0;j<momentum.x.size();j++){
  Real z=momentum.x[j];avg.add(3*momentum.w[j]*z*z*sqrt(Real(939.56542052)*Real(939.56542052)+pf*pf*z*z));
 }
 extra[11].add(measure*exp(-H)*avg.s);
}
int main(int argc,char**argv){try{
 if(argc!=7)throw std::runtime_error("Hc surfaceP N center_factor eos_order defect_transport");
 Real H=std::stold(argv[1]),ps=std::stold(argv[2]);
 Solver s(std::stoi(argv[3]),std::stoi(argv[5]),std::stold(argv[4]),false,false);
 x2_observer=observe;
 if(std::stoi(argv[6]))throw std::runtime_error("retired inconsistent defect estimator; qualification requires independent reference evidence");
 auto v=s.run(H,ps);x2_observer=nullptr;x2_step_observer=nullptr;
 Real K=log1p(-2*v[2]/v[1])/2+s.Hs;
 std::cout<<std::setprecision(21)<<"{\"primary\":[";
 for(size_t i=0;i<v.size();i++)std::cout<<(i?",":"")<<v[i];
 std::cout<<"],\"N\":[";
 for(int i=0;i<4;i++)std::cout<<(i?",":"")<<extra[i].s*Real(1e54);
 std::cout<<"],\"G_packed\":[";
 for(int i=0;i<4;i++)std::cout<<(i?",":"")<<extra[4+i].s*exp(-K)*Real(1e54);
 std::cout<<"],\"Lmu\":"<<extra[8].s*exp(-6*K)*Real(2e-36)
 <<",\"P0\":"<<exp(K)*extra[10].s/extra[9].s
 <<",\"P1\":"<<exp(K)*extra[11].s/extra[9].s
 <<",\"Hs\":"<<s.Hs<<",\"Hb\":"<<s.hb<<",\"Hn\":"<<s.eos.hn<<",\"Hmu\":"<<s.eos.hm
 <<",\"K\":"<<K<<",\"digits\":"<<std::numeric_limits<Real>::digits
 <<",\"defect_transport_uv\":["<<transported[0]<<","<<transported[1]<<"]"
 <<",\"defect_absolute_transport_uv\":["<<majorant[0]<<","<<majorant[1]<<"]"
 <<",\"stage_transport_uv\":["<<stage_effect[0]<<","<<stage_effect[1]<<"]"
 <<",\"defect_quadrature_spread\":"<<defect_quadrature_spread<<"}\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
