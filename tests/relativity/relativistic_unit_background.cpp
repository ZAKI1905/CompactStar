// ADR-0012 U3/U4: transformed CGS RHS at every canonical node, with
// independently coded constants; U12 export only (no particle-response production).
#include "tests/eos/structure1/table.hpp"
#include "CompactStar/Core/TOVSolver.hpp"
#include "CompactStar/Core/NStar.hpp"
#include <gsl/gsl_errno.h>
#include "CompactStar/Physics/Evolution/StarContext.hpp"
#include "CompactStar/Physics/Driver/Thermal/Boundary/TbDefinition.hpp"
#include <array>
#include <chrono>
#include <iostream>
using namespace CompactStar::Core;
using namespace structure1;
struct Probe:TOVSolver {const auto& labels()const{return eos_tab.extra_labels;} using TOVSolver::PressureCutoff;};
void run(Probe& t,const std::filesystem::path& dir,const std::string& name,double x,int radial,bool check,double rho_c=1.10e15)
{
 t.SetRadialRes(radial);std::vector<TOVPoint> v;
 require(t.SingleStarSolveToTOVPoints(rho_c*std::exp(x),v)>0&&t.LastSolveStatus()==TOVSolveStatus::SURFACE_REACHED,"canonical completion");
 require(std::abs(v.front().e/(rho_c*std::exp(x))-1)<1e-12,"achieved center/clamp");
 NStar star(v,t.labels()); require(star.ComputeHartleMonopoleResponse(),"monopole completion");
 const auto &p=star.Profile();const auto& h=*star.MonopoleResponse();const auto& f=star.RotationResponse();
 if(check) {
  NStar append;append.InitFromTOVSolver(&t);for(const auto& q:v)append.Append(q);append.FinalizeSurface();
  const auto& a=append.Profile().Radial();const auto& b=p.Radial();
  require(a.Dim().size()==b.Dim().size(),"U4 column schema");
  for(size_t k=0;k<a.Dim().size();++k)for(size_t j=0;j<v.size();++j)require(a[int(k)][j]==b[int(k)][j],"U4 profile path bits");
  require(append.MassSurface()==v.back().m && star.MassSurface()==v.back().m,"U4 literal public mass bits");
  require(append.GetSequence().ec==v.front().e && append.GetSequence().pc==v.front().p,"U4 physical center bits");
  require(append.IsSurfaceFinalized() && star.IsSurfaceFinalized() && (*p.GetPressure()).Size()==v.size(),"U4 endpoint");
  std::cout<<"U4 PASS all radial/species/metric columns identical; physical public bits retained.\n";
  CompactStar::Physics::Evolution::StarContext context(p);double rho_error=0;
  for(size_t i=0;i<v.size();i++){
   require((*p.GetBaryonDensity())[i]==v[i].rho,"U10 whole-baryon input bits");
   for(size_t k=0;k<t.labels().size();k++){
    const double Y=(*p.GetSpeciesPtr(t.labels()[k]))[i];
    require(Y==v[i].rho_i[k],"U10 species fractions preserve physical input");
    require(Y*(*p.GetBaryonDensity())[i]==v[i].rho_i[k]*v[i].rho,"U10 n_i=Y_i n_B");
   }
   rho_error=std::max(rho_error,std::abs((*context.MassDensity_gcm3())[i]/v[i].e-1));
  }
  require(rho_error<16*std::numeric_limits<double>::epsilon(),"U11 physical density round trip");
  for(double threshold:{1e9,1e10,1e12,1e14}){
   size_t expected=v.size();for(size_t i=v.size();i-->0;)if(v[i].e>=threshold){expected=i;break;}
   require(expected<v.size(),"physical threshold reached");
   require(CompactStar::Physics::Driver::Thermal::Boundary::FindTbIndex(context,threshold)==expected,"U11 physical threshold index");
  }
  std::cout<<std::setprecision(17)<<"U10 raw species fraction/number-density/whole-baryon semantics PASS; U11 inverse max="<<rho_error<<" threshold indices PASS\n";

  // U3 derivative here is the independently transformed CGS ODE, not a naive
  // finite difference of samples. Segment integral closure is checked separately.
  constexpr double G=6.673e-8,c=29979245800.,Ms=1.98892e33;
  std::array<double,5> worst{};
  for(size_t i=0;i<v.size();++i) {
   double r=(*p.GetRadius())[i],m=(*p.GetMass())[i],e=(*p.GetEnergyDensity())[i],pr=(*p.GetPressure())[i],np=(*p.GetMetricNuPrime())[i];
   double rc=v[i].r*1e5,grams=v[i].m*Ms,rho=v[i].e,pc=v[i].p;
   double mprime=4*M_PI*rc*rc*rho * G/(c*c); // (g/cm)*(km/g)*(cm/km)
   double physical_np=G/(rc*rc*c*c)*(grams+4*M_PI*rc*rc*rc*pc/(c*c))/(1-2*G*grams/(rc*c*c));
   double physical_dp=-(rho*c*c+pc)*physical_np;
   double pprime=physical_dp*G/(c*c*c*c)*1e15;
   double F=(m+4*M_PI*r*r*r*pr)/(r*(r-2*m));
   std::array<double,5> lhs{std::exp(-2*(*p.GetMetricLambda())[i]),mprime,pprime,np,pprime};
   std::array<double,5> rhs{1-2*m/r,4*M_PI*r*r*e,-(e+pr)*F,F,-(e+pr)*np};
   for(size_t k=0;k<5;++k)worst[k]=std::max(worst[k],std::abs(lhs[k]/rhs[k]-1));
  }
  std::cout<<std::setprecision(17)<<"U3 metric mass hydrostatic nuprime pressure_nuprime:";
  for(double w:worst){std::cout<<' '<<w;require(w<8e-15,"U3 transformed CGS/GSL identity");}std::cout<<'\n';
 }
 std::ofstream o(dir/(name+".tsv"));o<<std::setprecision(17);
 o<<"r m eps p nb nu nup dedp s sp mh ph xi nn np ne nm lambda r_cgs m_cgs_recovered rho_cgs p_cgs tp_m tp_nup tp_nu\n";
 for(size_t i=0;i<v.size();++i){
 o<<(*p.GetRadius())[i]<<' '<<(*p.GetMass())[i]<<' '<<(*p.GetEnergyDensity())[i]<<' '<<(*p.GetPressure())[i]<<' '<<(*p.GetBaryonDensity())[i]<<' '<<(*p.GetMetricNu())[i]<<' '<<(*p.GetMetricNuPrime())[i]<<' '<<(*p.GetEosDEdP())[i]<<' '<<f.omega_bar_over_Omega[i]<<' '<<f.domega_bar_over_Omega_dr[i]<<' '<<h.m0_over_Omega2[i]<<' '<<h.p0star_over_Omega2[i]<<' '<<h.xi0_over_Omega2[i];
 for(const char* label:{"10","11","0","1"})o<<' '<<(*p.GetSpeciesPtr(label))[i]*(*p.GetBaryonDensity())[i];
 o<<' '<<(*p.GetMetricLambda())[i]<<' '<<v[i].r*1e5<<' '<<v[i].m*1.98892e33<<' '<<v[i].e<<' '<<v[i].p<<' '<<v[i].m<<' '<<v[i].nu_der<<' '<<v[i].nu<<'\n';}
 std::ofstream meta(dir/"meta.tsv",std::ios::app);meta<<std::setprecision(17)<<name<<' '<<x<<' '<<radial<<' '<<v.size()<<' '<<v.front().e<<' '<<t.GetInitPress()<<' '<<t.PressureCutoff()<<' '<<f.I<<' '<<h.deltaM_over_Omega2<<' '<<h.surface_xi0_over_Omega2<<' '<<0<<' '<<h.surface_shell_mass_over_Omega2<<' '<<star.GetSequence().b<<' '<<v.back().m<<'\n';
 std::cout<<name<<" done\n"<<std::flush;
}
int main(int argc,char**argv)
{
 try {
  gsl_set_error_handler_off();TrackRFreeGasThermodynamicProvider eos;
  auto dir=std::filesystem::absolute(argc>1?argv[1]:"unit-background-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
  require(!std::filesystem::exists(dir),"fresh artifact directory required");std::filesystem::create_directories(dir);
  if(argc>3 && std::string(argv[2])=="--cmf") {
   Probe t;t.ImportEOS(std::string(argv[3])+"/DS-CMF-1-with-crust/DS(CMF)-1_with_crust.eos",true);
   std::vector<TOVPoint> selected;std::vector<std::string> labels;
   require(t.SolveToProfile(1.6,selected,&labels)>0,"CMF target selection");
   const double rho=selected.front().e;
   for(int n:{10000,20000,40000,80000})run(t,dir,"cmf-r"+std::to_string(n),0,n,true,rho);
   return 0;
  }
  bool campaign=argc>2 && std::string(argv[2])=="--campaign";
  for(int g:{4096,8192}) {
   if(!campaign && g!=8192)continue;
   auto file=dir/("eos-"+std::to_string(g)+".tsv");generate(eos,file,g);Probe t;t.ImportEOS(file.string(),true);
   for(int n:{20000,40000}) {if(!campaign && n!=40000)continue;run(t,dir,"g"+std::to_string(g)+"r"+std::to_string(n),0,n,true);}
   if(campaign)for(int j:{-4,-2,-1,1,2,4})run(t,dir,"seq"+std::to_string(g)+"_"+std::to_string(j),j*.0005,40000,false);
  }
  return 0;
 }catch(const std::exception& e){std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}
}
