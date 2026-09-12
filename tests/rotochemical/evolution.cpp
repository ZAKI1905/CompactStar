#include "coupled_oracles.hpp"
#include "trajectory.hpp"
int main(int argc,char**argv){try{
 gsl_set_error_handler_off(); // Repository executable policy: transport numerical failures, never abort.
 require(argc==7,"profile certificate thermal output mode entry-manifest");std::cout<<std::setprecision(17)<<std::unitbuf;
 auto f=Fixture(argv[1],argv[2],std::filesystem::path(argv[4])/"owning-star",80000);auto c=Channels(f);auto fine=Channels(f,32,true);
 for(auto p:{RC::UrcaProcess::Me,RC::UrcaProcess::Mmu}){double a=c->LuminosityCoefficient(p),b=fine->LuminosityCoefficient(p);Near(a,b,1e-6,"Ltilde quadrature refinement failed");std::cout<<"LTILDE "<<RC::ProcessIndex(p)<<' '<<a<<' '<<b<<" relative "<<std::abs(a-b)/b<<'\n';}
 EOS::CompOSE_Thermo::Options options;options.Tmin_for_derivative_MeV=0;options.clamp_to_domain=false;auto thermal=std::make_shared<const RC::FrozenThermalSource>(argv[3],options,"controlled mathematical fixed-background free-gas entropy; qualified radial80000");
 RC::RunQualification q;q.radial_resolution=80000;q.eos_resolution=8192;q.rho_c_g_cm3=1.10e15;q.profile_path=(std::filesystem::path(argv[1])/"profile.tsv").string();q.model_path=(std::filesystem::path(argv[1])/"model.txt").string();q.eos_path=(std::filesystem::path(argv[1])/"freegas.tsv").string();q.certificate_path=argv[2];q.entry_manifest_path=argv[6];q.entry_manifest_sha256=RC::FrozenSource(argv[6]).sha256;
 if(std::string(argv[5])=="oracles")CoupledOracles(f,c,thermal,q);
 auto token=std::make_shared<RC::RunDependencyToken>();auto spin=std::make_shared<const RC::PrescribedDipoleHistory>(token);auto context=Context(f,c,thermal,spin,q,token);FixtureReport(*context);
 if(std::string(argv[5])=="trajectory"){
  RadialLtildeComparison(f,*context,std::filesystem::path(argv[1])/"freegas.tsv",std::filesystem::path(argv[4])/"radial20000-comparison");
  std::cout<<"PRE_TRAJECTORY_READY"<<std::endl;std::string gate;require(bool(std::cin>>gate)&&gate=="PROTECTED_GOVERNED_GATES_PASS","pre-trajectory gate did not authorize execution");context->RequireFullCurrent();
  Trajectory(context,std::filesystem::path(argv[4])/"trajectory.tsv");Trajectory(context,std::filesystem::path(argv[4])/"refined.tsv",true);
  Trajectory(context,std::filesystem::path(argv[4])/"initial-T1e7.tsv",false,1e7);Trajectory(context,std::filesystem::path(argv[4])/"initial-T1e9.tsv",false,1e9);
  Trajectory(context,std::filesystem::path(argv[4])/"initial-xi1.tsv",false,1e8,1);Trajectory(context,std::filesystem::path(argv[4])/"initial-xi20.tsv",false,1e8,20);
 }
 // Destroy every caller semantic handle; context retains the complete lifetime chain.
 std::weak_ptr<Core::NStar> weak_central=f.central;f.central.reset();f.provider.reset();f.g.reset();f.z.reset();f.fixed.reset();f.metric={};f.interpolate={};c.reset();fine.reset();thermal.reset();spin.reset();token.reset();context->RequireFullCurrent();RunState retained(context);RC::SecularEvolutionDriver retained_driver(context);auto value=retained_driver.Evaluate(0,retained.state,retained.ctx);require(std::isfinite(value.x_dot_s),"context failed retained lifetime");std::cout<<"CALLER_LIFETIME PASS all caller handles destroyed\n";
 if(std::string(argv[5])=="oracles"){const_cast<Core::StarProfile&>(weak_central.lock()->Profile()).Touch(); // Test mutation of an actually non-const owning star, matching governed refusal controls.
retained.rhs.AddTo(Tag::Thermal,0,135);retained.rhs.AddTo(Tag::Chem,0,246);retained.rhs.AddTo(Tag::Chem,1,357);MustRefuse([&]{retained_driver.AccumulateRHS(0,retained.state,retained.rhs,retained.ctx);},"actual upstream profile version changed");require(retained.rhs.Peek(Tag::Thermal,0)==135&&retained.rhs.Peek(Tag::Chem,0)==246&&retained.rhs.Peek(Tag::Chem,1)==357,"profile change altered derivative sentinel");}
 return 0;
 }catch(const std::exception&e){std::cerr<<"STOP "<<e.what()<<'\n';return 1;}}
