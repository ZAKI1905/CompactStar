// Validation-only control of existing owners. Historical source is unchanged.
#include "ba13_original_helpers.hpp"
struct StudySolver:Solver{void Center(double cm){r_min=cm;}};
std::shared_ptr<Core::NStar> StudyStar(const std::string& table,double rho,size_t radial,double center,const std::filesystem::path& work)
{
 if(center==1.)return solve(table,rho,radial,work);
 StudySolver s;s.SetWrkDir(work.string());s.ImportEOS(table,true);s.SetRadialRes(radial);s.Center(center);
 std::vector<Core::TOVPoint> pts;
 require(s.SingleStarSolveToTOVPoints(rho,pts)>0&&s.LastSolveStatus()==Core::TOVSolveStatus::SURFACE_REACHED,"study star solve refused");
 auto star=std::make_shared<Core::NStar>(pts,s.labels());require(star->ComputeHartleMonopoleResponse(),"study monopole refused");return star;
}
int main(int argc,char** argv)
{
 try {
 gsl_set_error_handler_off();
 require(argc==10,"mode nominal-raw out radial eos center tau index lower-density");
 const std::string mode=argv[1];const std::filesystem::path raw=argv[2],out=argv[3];
 size_t radial=std::stoull(argv[4]),eos=std::stoull(argv[5]),index=std::stoull(argv[8]);
 double center=std::stod(argv[6]),taufrac=std::stod(argv[7]),lower=std::stod(argv[9]);
 const auto nominal=raw/"chemical-profiles"/("target-"+std::to_string(index));
 TrackRFreeGasThermodynamicProvider provider;
 if(mode=="star") {
 require(!std::filesystem::exists(out),"fresh study output required");std::filesystem::create_directories(out/"inputs");
 for(const auto& n:{"model.txt","windows.tsv"})std::filesystem::copy_file(nominal/n,out/"inputs"/n);
 if(eos==8192)std::filesystem::copy_file(nominal/"freegas.tsv",out/"inputs/freegas.tsv");
 else generate(provider,out/"inputs/freegas.tsv",eos);
 if(lower<1e-14){
   // Preserve all original EOS nodes; extend only below the nominal floor.
   const auto extra=generate(provider,out/"low-density.tsv",64,lower,1e-14);
   std::ifstream old(out/"inputs/freegas.tsv"),low(extra.file);std::string h,line,rest;std::getline(old,h);
   while(std::getline(old,line))rest+=line+'\n';old.close();std::getline(low,line);
   std::ofstream extended(out/"inputs/freegas.tsv");extended<<h<<'\n';
   while(std::getline(low,line)){std::istringstream rr(line);double en,pr,nb;rr>>en>>pr>>nb;if(nb<1e-14)extended<<line<<'\n';}
   extended<<rest;extended.close();require(bool(extended),"extended EOS write");
 }
 auto src=source(out/"inputs/freegas.tsv");
 std::ifstream targets(raw/"target_stars.tsv");std::string line;std::getline(targets,line);
 double rho0=0,target=0,B0=0;size_t idx=0;
 while(std::getline(targets,line)){std::istringstream row(line);double fraction,rho,B,Btarget;row>>idx>>fraction>>rho>>B>>Btarget;if(idx==0)B0=B;if(idx==index){rho0=rho;target=Btarget;}}
 require(rho0>0&&target>0,"target missing");size_t evals=0;
 auto at=[&](double rho){TargetStar x;x.rho_c=rho;x.target=target;x.star=StudyStar(src->table_path,rho,radial,center,out/("solve-"+std::to_string(evals++)));x.numbers=ParticleNumbers::Compute(input(*x.star,src));x.numbers.RequireCurrent();x.B=x.numbers.Values()[0]+x.numbers.Values()[1];x.residual=x.B-target;return x;};
 auto mid=at(rho0),lo=mid,hi=mid;double step=rho0*1e-4;
 if(mid.residual!=0){
  for(int n=0;n<10&&(lo.residual>0||hi.residual<0);++n){if(lo.residual>0)lo=at(rho0-step);if(hi.residual<0)hi=at(rho0+step);step*=2;}
  require(lo.residual<=0&&hi.residual>=0,"common-B bracket unavailable");
  for(int n=0;n<80&&hi.B-lo.B>taufrac*B0;++n){double rho=lo.rho_c+(hi.rho_c-lo.rho_c)*(-lo.residual)/(hi.residual-lo.residual);
   const double margin=64*std::numeric_limits<double>::epsilon()*rho0;
   if(!(rho>lo.rho_c+margin&&rho<hi.rho_c-margin)||n%6==5)rho=(lo.rho_c+hi.rho_c)/2;
   auto v=at(rho);if(v.residual<=0)lo=std::move(v);else hi=std::move(v);
  }
 }
 require(hi.B-lo.B<=taufrac*B0,"common-B bracket precision failure");
 auto best=std::abs(lo.residual)<=std::abs(hi.residual)?lo:hi;
 WriteChemicalProfile(*best.star,provider,out/"inputs",out/"profile");
 std::ofstream f(out/"star.tsv");f<<std::setprecision(17)<<"index\tradial\teos\tcenter_cm\trho_c\tB_target\tB_solved\tresidual\tbracket\tR\tM\tnu_surface\tN_n\tN_e\tN_mu\tu_N_n\tu_N_e\tu_N_mu\trho_lo\trho_hi\tB_lo\tB_hi\tR_lo\tR_hi\tM_lo\tM_hi\tevaluations\n";
 const auto& p=best.star->Profile();const auto& n=best.numbers.Values();const auto& u=best.numbers.Errors();
 f<<index<<'\t'<<radial<<'\t'<<eos<<'\t'<<center<<'\t'<<best.rho_c<<'\t'<<target<<'\t'<<best.B<<'\t'<<best.residual<<'\t'<<hi.B-lo.B<<'\t'<<(*p.GetRadius())[-1]<<'\t'<<(*p.GetMass())[-1]<<'\t'<<(*p.GetMetricNu())[-1]<<'\t'<<n[0]<<'\t'<<n[2]<<'\t'<<n[3]<<'\t'<<u[0]<<'\t'<<u[2]<<'\t'<<u[3]<<'\t'<<lo.rho_c<<'\t'<<hi.rho_c<<'\t'<<lo.B<<'\t'<<hi.B<<'\t'<<(*lo.star->Profile().GetRadius())[-1]<<'\t'<<(*hi.star->Profile().GetRadius())[-1]<<'\t'<<(*lo.star->Profile().GetMass())[-1]<<'\t'<<(*hi.star->Profile().GetMass())[-1]<<'\t'<<evals<<'\n';f.close();require(bool(f),"study star write failure");
 std::cout<<"STUDY_COMMON_B_STAR PASS evaluations "<<evals<<'\n';return 0;
 }
 std::ifstream meta(out/"star.tsv");std::string line;std::getline(meta,line);std::getline(meta,line);std::istringstream row(line);double i,rad,eo,ce,rho,target,B;row>>i>>rad>>eo>>ce>>rho>>target>>B;
 auto star=StudyStar((out/"profile/freegas.tsv").string(),rho,radial,center,out/(mode+"-owning-star"));
 auto src=source(out/"profile/freegas.tsv");auto nums=ParticleNumbers::Compute(input(*star,src));require(nums.Values()[0]+nums.Values()[1]==B,"study observable star identity differs");
 if(mode=="thermal"){
  CompactStar::EOS::CompOSE_Thermo::Options op;op.Tmin_for_derivative_MeV=0;op.clamp_to_domain=false;
  RC::FrozenThermalSource thermal(out/"profile/thermal",op,"unchanged Phase-5D adapter: study candidate");
  EV::StarContext sc(star->Profile());EV::GeometryCache geo(sc);EV::Config config;
  TH::PhotonCooling::Options po;po.surface_model=TH::PhotonCooling::Options::SurfaceModel::EnvelopeTbTs;TH::PhotonCooling photon(po);
  std::ofstream f(out/"thermal.tsv");f<<std::setprecision(17)<<"index\tTinf_K\tCstar\tTsurface\tTsurface_inf\tg14\n";
  for(int k=0;k<=12;++k){double T=std::pow(10.,6+.25*k);CompactStar::Physics::State::ThermalState st;st.Resize(1);st.SetTinf(T);EV::StateVector state;state.Register(CompactStar::Physics::State::StateTag::Thermal,st);EV::DriverContext ctx;ctx.star=&sc;ctx.geo=&geo;ctx.cfg=&config;ctx.thermo=&thermal.Table();auto d=TH::Detail::ComputeDerived(photon,state,ctx);require(d.ok,"thermal study refused");f<<index<<'\t'<<T<<'\t'<<d.C_star_erg_K<<'\t'<<d.Tsurf_K<<'\t'<<d.Tsurf_K*std::sqrt(d.exp2nu_surf)<<'\t'<<TH::Boundary::SurfaceGravity_g14(sc,&geo)<<'\n';}
  f.close();require(bool(f),"thermal write");return 0;
 }
 if(mode=="events"){
  std::ofstream f(out/"events.tsv");f<<std::setprecision(17)<<"order\tmu_B\tP0\tP1\tP2\tneutron_measure\n";
  for(unsigned order:{16u,32u,64u,96u}){auto d=DirectReferenceAverages(*star,provider,order);f<<order<<'\t'<<d.mu_B_inf_MeV<<'\t'<<d.p0_inf_MeV<<'\t'<<d.p1_inf_MeV<<'\t'<<d.p2_inf_MeV<<'\t'<<d.neutron_measure<<'\n';}
  return 0;
 }
 throw std::runtime_error("unknown study mode");
 }catch(const std::exception& e){std::cerr<<"STOP "<<e.what()<<'\n';return 1;}
}
