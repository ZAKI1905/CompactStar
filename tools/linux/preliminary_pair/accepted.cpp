// Observe retained references at every main accepted endpoint, without t0 solves.
#define NATIVE_STAGE_NO_MAIN
#include "../../../tests/native/stage.cpp"
#undef NATIVE_STAGE_NO_MAIN
#include "retained_helpers.inc"
int main(int argc,char** argv){try{
 gsl_set_error_handler_off();require(argc==12,"profile thermal entry frozen coefficients archive sha work source control output");
 Inputs in;in.profile=argv[1];in.certificate=in.profile/"certificate.txt";in.thermal=argv[2];in.entry=argv[3];in.frozen=argv[4];in.coefficients=argv[5];in.work=argv[8];
 auto sequence=Core::StellarArchive::Load(argv[6],argv[7]);auto pair=NativePair(in,sequence);
 std::filesystem::path out=argv[11];require(!std::filesystem::exists(out),"fresh accepted evidence required");std::filesystem::create_directories(out);
 for(int role=0;role<2;++role){
  auto p=role==0?pair.source:pair.control;std::filesystem::path retained=argv[9+role],dir=out/(role==0?"source":"control");std::filesystem::create_directories(dir);
  auto main=ReadRetainedMain(retained/"tier-0-cap-0");std::vector<Solution> refs(5);
  for(int k=0;k<5;++k){std::ifstream f(retained/("reference-"+std::to_string(k)+"-history.tsv"));std::string line;std::getline(f,line);State y{{0,0,0}};
   while(std::getline(f,line)){Stepper::Step s{};s.before=y;std::istringstream r(line);r>>s.left>>s.right;for(auto& v:s.after)r>>v;require(bool(r),"bad history");y=s.after;refs[k].steps.push_back(s);}
  }
  NativeQualification::ProcessPool(main.accepted_steps.size(),24,dir/"workers",[&](std::size_t i,std::size_t){
   auto step=main.accepted_steps[i];const double t=step.t_right_s;ActualCheckpointContext context(p,1);
   std::ofstream f(dir/("point-"+std::to_string(i)+".tsv"));f<<std::setprecision(17)<<"reference\tt\tx\te\tmu";
   for(const auto& kv:DiagnosticValues(context.EvaluateDiagnostics(t,step.y_right)))f<<'\t'<<kv.first;f<<'\n';
   for(int k=-1;k<5;++k){auto y=k<0?step.y_right:observe(p,refs[k],t);f<<k<<'\t'<<t<<'\t'<<values(y);for(const auto& kv:DiagnosticValues(context.EvaluateDiagnostics(t,y)))f<<'\t'<<kv.second;f<<'\n';}
   f.flush();require(bool(f),"accepted observation write failure");
  });
 }
 std::cout<<"PASS both roles observed at all accepted endpoints\n";return 0;
}catch(const std::exception& e){std::cerr<<"STOP "<<e.what()<<'\n';return 1;}}
