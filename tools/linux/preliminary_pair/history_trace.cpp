// One-step residual evidence on retained histories. Never reruns a source trajectory.
#define NATIVE_STAGE_NO_MAIN
#include "../../../tests/native/stage.cpp"
#undef NATIVE_STAGE_NO_MAIN
#include "retained_helpers.inc"
struct MapRow{double a,b,h;State before,after;};
std::vector<MapRow> maps(const std::filesystem::path& path,bool native){
 std::ifstream f(path);require(bool(f),"missing retained history");std::string line;std::getline(f,line);
 std::vector<MapRow> result;State y{{0,0,0}};
 while(std::getline(f,line)){
  std::istringstream in(line);MapRow r;r.before=y;
  if(native){std::size_t ordinal,rejected;double next;in>>ordinal>>r.a>>r.b>>r.h>>next>>rejected;}
  else{in>>r.a>>r.b;r.h=r.b-r.a;}
  for(double& v:r.after)in>>v;require(bool(in),"bad retained history row");y=r.after;result.push_back(r);
 }return result;
}
int main(int argc,char** argv){
 try{
  gsl_set_error_handler_off();std::cout<<std::setprecision(17)<<std::unitbuf;
  require(argc==12,"profile thermal entry frozen coefficients archive sha work histories output role");
  Inputs in;in.profile=argv[1];in.certificate=in.profile/"certificate.txt";in.thermal=argv[2];in.entry=argv[3];in.frozen=argv[4];in.coefficients=argv[5];in.work=argv[8];
  auto sequence=Core::StellarArchive::Load(argv[6],argv[7]);auto pair=NativePair(in,sequence);
  std::string role=argv[11];require(role=="source"||role=="control","invalid role");auto p=role=="source"?pair.source:pair.control;
  const std::filesystem::path retained=argv[9],out=argv[10];require(!std::filesystem::exists(out),"fresh traces required");std::filesystem::create_directories(out);
  NativeQualification::ProcessPool(6,6,out/"workers",[&](std::size_t task,std::size_t){
   const bool native=task==5;const auto input=native?retained/"tier-0-cap-0/main.accepted_states.tsv":retained/("reference-"+std::to_string(task)+"-history.tsv");
   const auto records=maps(input,native);const auto dir=out/(native?"main":"reference-"+std::to_string(task));std::filesystem::create_directories(dir);
   std::ofstream rhs(dir/"local-rhs.tsv"),diagnostics(dir/"diagnostics.tsv"),mapping(dir/"maps.tsv");
   rhs<<std::setprecision(17)<<"context\tt_s\tcell\tx\te\tmu\tdx\tde\tdmu\n";mapping<<std::setprecision(17)<<"context\tleft\tright\th\tx\te\tmu\treplay_x\treplay_e\treplay_mu\tsteps\trejected\n";
   TraceContext context(p,0,rhs,diagnostics);
   for(const auto& r:records){
    ++context.id;auto y=r.before;double t=r.a,h=r.b-r.a;std::size_t count=0;
    Stepper solver([&](double tt,const State& yy,State& d,int cell){context.SmoothDerivative(tt,yy,d,cell);},context.Breakpoints(),task<3?gsl_odeiv2_step_rk8pd:gsl_odeiv2_step_rkf45,native?tiers[0].r:(task==0?tight.r:tighter.r),native?tiers[0].a:(task==0?tight.a:tighter.a),r.b-r.a);
    while(t<r.b){require(++count<=16,"unexpected one-step reconstruction refinement");solver.Advance(t,y,h,r.b);}
    mapping<<context.id<<'\t'<<r.a<<'\t'<<r.b<<'\t'<<r.h<<'\t'<<values(r.after)<<'\t'<<values(y)<<'\t'<<count<<'\t'<<solver.Stats().rejected<<'\n';
   }
   rhs.flush();mapping.flush();require(bool(rhs)&&bool(mapping),"trace write failed");
  });
  std::cout<<"PASS six retained histories mapped to local arithmetic traces\n";return 0;
 }catch(const std::exception& e){std::cerr<<"STOP "<<e.what()<<'\n';return 1;}
}
