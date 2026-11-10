// Read-only numerical evidence from retained Stage B histories.
#define NATIVE_STAGE_NO_MAIN
#include "../../../tests/native/stage.cpp"
#undef NATIVE_STAGE_NO_MAIN
#include "trace_context.inc"
struct MapRow{double a,b;State before,after;};
std::vector<MapRow> ReadMaps(const std::filesystem::path& path,bool native){
 std::ifstream f(path);require(bool(f),"missing retained history");std::string line;std::getline(f,line);
 std::vector<MapRow> result;State y{{0,0,0}};
 while(std::getline(f,line)){
  std::istringstream in(line);MapRow r;r.before=y;
  if(native){std::size_t ordinal,rejected;double h,next;in>>ordinal>>r.a>>r.b>>h>>next>>rejected;}
  else in>>r.a>>r.b;
  for(double& v:r.after)in>>v;require(bool(in),"bad retained history row");y=r.after;result.push_back(r);
 }return result;
}
int main(int argc,char** argv){try{
 gsl_set_error_handler_off();std::cout<<std::setprecision(17)<<std::unitbuf;
 require(argc==12,"profile thermal entry frozen coefficients archive sha work histories output role");
 Inputs in;in.profile=argv[1];in.certificate=in.profile/"certificate.txt";in.thermal=argv[2];in.entry=argv[3];in.frozen=argv[4];in.coefficients=argv[5];in.work=argv[8];
 auto sequence=Core::StellarArchive::Load(argv[6],argv[7]);auto pair=NativePair(in,sequence);
 std::string role=argv[11];require(role=="source"||role=="control","invalid role");auto p=role=="source"?pair.source:pair.control;
 const std::filesystem::path retained=argv[9],out=argv[10];require(!std::filesystem::exists(out),"fresh supplemental evidence required");std::filesystem::create_directories(out);
 const auto workers=StageWorkers();
 std::vector<Solution> refs(5);
 for(int k=0;k<5;++k)for(const auto& r:ReadMaps(retained/("reference-"+std::to_string(k)+"-history.tsv"),false)){
  Stepper::Step s{};s.left=r.a;s.right=r.b;s.before=r.before;s.after=r.after;refs[k].steps.push_back(s);
 }
 // Reference maps are independent; local steps are reconstructed from each exact
 // retained left state. No original-t0 reference or main trajectory is rerun.
 NativeQualification::ProcessPool(5,workers,out/"trace-workers",[&](std::size_t task,std::size_t){
  const auto dir=out/("reference-"+std::to_string(task));std::filesystem::create_directories(dir);
  std::ofstream rhs(dir/"local-rhs.tsv"),diagnostics(dir/"diagnostics.tsv"),mapping(dir/"maps.tsv");
  rhs<<std::setprecision(17)<<"context\tt_s\tcell\tx\te\tmu\tdx\tde\tdmu\n";mapping<<std::setprecision(17)<<"context\tleft\tright\th\tx\te\tmu\treplay_x\treplay_e\treplay_mu\tsteps\trejected\n";
  TraceContext context(p,0,rhs,diagnostics);
  for(const auto& r:refs[task].steps){
   ++context.id;auto y=r.before;double t=r.left,h=r.right-r.left;std::size_t count=0;
   Stepper solver([&](double tt,const State& yy,State& d,int cell){context.SmoothDerivative(tt,yy,d,cell);},context.Breakpoints(),task<3?gsl_odeiv2_step_rk8pd:gsl_odeiv2_step_rkf45,task==0?tight.r:tighter.r,task==0?tight.a:tighter.a,r.right-r.left);
   while(t<r.right){require(++count<=16,"unexpected one-step reconstruction refinement");solver.Advance(t,y,h,r.right);}
   mapping<<context.id<<'\t'<<r.left<<'\t'<<r.right<<'\t'<<r.right-r.left<<'\t'<<values(r.after)<<'\t'<<values(y)<<'\t'<<count<<'\t'<<solver.Stats().rejected<<'\n';
  }
  rhs.flush();mapping.flush();require(bool(rhs)&&bool(mapping),"trace write failed");
 });
 std::set<double> unique;
 for(int tier=0;tier<3;++tier)for(int cap=0;cap<4;++cap){
  auto dir=retained/("tier-"+std::to_string(tier)+"-cap-"+std::to_string(cap));
  for(const auto& r:ReadMaps(dir/"main.accepted_states.tsv",true))unique.insert(r.b);
  if(cap==3)for(const auto& r:ReadMaps(dir/"restart-tail/main.accepted_states.tsv",true))unique.insert(r.b);
 }
 std::vector<double> times(unique.begin(),unique.end());
 NativeQualification::ProcessPool(times.size(),workers,out/"accepted-workers",[&](std::size_t i,std::size_t){
  const double t=times[i];ActualCheckpointContext context(p,1);
  std::ofstream f(out/("accepted-"+std::to_string(i)+".tsv"));f<<std::setprecision(17)<<"reference\tt\tx\te\tmu";
  for(const auto& kv:DiagnosticValues(context.EvaluateDiagnostics(0,{{0,0,0}})))f<<'\t'<<kv.first;f<<'\n';
  for(int k=0;k<5;++k){auto y=observe(p,refs[k],t);f<<k<<'\t'<<t<<'\t'<<values(y);for(const auto& kv:DiagnosticValues(context.EvaluateDiagnostics(t,y)))f<<'\t'<<kv.second;f<<'\n';}
  f.flush();require(bool(f),"accepted observation write failure");
 });
 std::cout<<"PASS retained references traced and "<<times.size()<<" unique main/restart endpoints observed\n";return 0;
}catch(const std::exception& e){std::cerr<<"STOP "<<e.what()<<'\n';return 1;}}
