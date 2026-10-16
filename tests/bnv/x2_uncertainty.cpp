// Passive uncertainty assessment on X2's existing GL4 path, not a third solver.
// UNVERIFIED SCIENTIFIC CANDIDATE. See docs/validation/ba13-headroom/METHOD.md.
#define main x2_historical_main
#include "x2_fixed_step.cpp"
#undef main
std::array<Sum,12> extra;
Rule momentum(32);
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
 if(argc!=7)throw std::runtime_error("Hc surfaceP N center_factor eos_order output_profile_reserved");
 Real H=std::stold(argv[1]),ps=std::stold(argv[2]);
 Solver s(std::stoi(argv[3]),std::stoi(argv[5]),std::stold(argv[4]),false,false);
 x2_observer=observe;auto v=s.run(H,ps);x2_observer=nullptr;
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
 <<",\"K\":"<<K<<",\"digits\":"<<std::numeric_limits<Real>::digits<<"}\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
