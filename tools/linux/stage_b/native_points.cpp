// Native diagnostics at retained main/restart accepted states, no evolution.
#define NATIVE_STAGE_NO_MAIN
#include "../../../tests/native/stage.cpp"
#undef NATIVE_STAGE_NO_MAIN
struct NativePoint {double t;State y;std::size_t ordinal;};
std::vector<NativePoint> ReadPoints(const std::filesystem::path& path){
 std::ifstream f(path);require(bool(f),"missing complete accepted history");std::string line;std::getline(f,line);std::vector<NativePoint> result;
 while(std::getline(f,line)){std::istringstream r(line);NativePoint p;double left,h,next;std::size_t rejected;r>>p.ordinal>>left>>p.t>>h>>next>>rejected;for(auto& v:p.y)r>>v;require(bool(r),"bad accepted history");result.push_back(p);}return result;
}
int main(int argc,char** argv){try{
 gsl_set_error_handler_off();require(argc==12,"profile thermal entry frozen coefficients archive sha work histories output role");
 Inputs in;in.profile=argv[1];in.certificate=in.profile/"certificate.txt";in.thermal=argv[2];in.entry=argv[3];in.frozen=argv[4];in.coefficients=argv[5];in.work=argv[8];
 auto sequence=Core::StellarArchive::Load(argv[6],argv[7]);auto pair=NativePair(in,sequence);std::string role=argv[11];require(role=="source"||role=="control","invalid role");auto p=role=="source"?pair.source:pair.control;
 std::filesystem::path retained=argv[9],out=argv[10];require(!std::filesystem::exists(out),"fresh native diagnostic evidence required");std::filesystem::create_directories(out);
 NativeQualification::ProcessPool(15,StageWorkers(),out/"workers",[&](std::size_t task,std::size_t){
  const bool restart=task>=12;int tier=restart?task-12:task/4,cap=restart?3:task%4;
  auto name="tier-"+std::to_string(tier)+"-cap-"+std::to_string(cap);auto path=retained/name;
  if(restart)path/="restart-tail";
  auto points=ReadPoints(path/"main.accepted_states.tsv");ActualCheckpointContext context(p,task+1);
  std::ofstream f(out/(name+(restart?"-restart":"")+".tsv"));f<<std::setprecision(17)<<"ordinal\tt\tx\te\tmu";
  for(const auto& kv:DiagnosticValues(context.EvaluateDiagnostics(0,{{0,0,0}})))f<<'\t'<<kv.first;f<<'\n';
  for(const auto& q:points){auto d=context.EvaluateDiagnostics(q.t,q.y);require(d.valid_through_sample,"native accepted diagnostic guard");f<<q.ordinal<<'\t'<<q.t<<'\t'<<values(q.y);for(const auto& kv:DiagnosticValues(d))f<<'\t'<<kv.second;f<<'\n';}
  f.flush();require(bool(f),"native diagnostic write failure");
 });return 0;
}catch(const std::exception& e){std::cerr<<"STOP "<<e.what()<<'\n';return 1;}}
