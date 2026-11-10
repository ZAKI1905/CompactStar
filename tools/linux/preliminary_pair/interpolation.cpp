// Recover actual reference-observation arithmetic without a t0 trajectory solve.
#define NATIVE_STAGE_NO_MAIN
#include "../../../tests/native/stage.cpp"
#undef NATIVE_STAGE_NO_MAIN
#include "retained_helpers.inc"
struct ReferenceStep{double left,right;State before,after;};
int main(int argc,char** argv){
 try{
  gsl_set_error_handler_off();require(argc==12,"profile thermal entry frozen coefficients archive sha work histories output role");
  Inputs in;in.profile=argv[1];in.certificate=in.profile/"certificate.txt";in.thermal=argv[2];in.entry=argv[3];in.frozen=argv[4];in.coefficients=argv[5];in.work=argv[8];
  auto sequence=Core::StellarArchive::Load(argv[6],argv[7]);auto pair=NativePair(in,sequence);
  std::string role=argv[11];require(role=="source"||role=="control","invalid role");auto p=role=="source"?pair.source:pair.control;
  std::filesystem::path retained=argv[9],out=argv[10];require(!std::filesystem::exists(out),"fresh interpolation output required");std::filesystem::create_directories(out);
  std::vector<std::vector<ReferenceStep>> histories(5);std::vector<std::vector<State>> observations(5);std::vector<double> times;
  for(int k=0;k<5;++k){
   std::ifstream f(retained/("reference-"+std::to_string(k)+"-history.tsv"));std::string line;std::getline(f,line);State y{{0,0,0}};
   while(std::getline(f,line)){ReferenceStep s;s.before=y;std::istringstream r(line);r>>s.left>>s.right;for(auto& v:s.after)r>>v;require(bool(r),"bad reference history");y=s.after;histories[k].push_back(s);}
   std::ifstream obs(retained/("reference-"+std::to_string(k)+"-observations.tsv"));std::getline(obs,line);std::size_t i=0;
   while(std::getline(obs,line)){std::istringstream r(line);double t;State s;r>>t;for(auto& v:s)r>>v;require(bool(r),"bad observation");if(k==0)times.push_back(t);else require(t==times.at(i),"reference time mismatch");observations[k].push_back(s);++i;}
   require(observations[k].size()==1921,"reference observation count");std::filesystem::create_directories(out/("reference-"+std::to_string(k)));
  }
  NativeQualification::ProcessPool(5*1921,24,out/"workers",[&](std::size_t task,std::size_t){
   const auto k=task/1921,index=task%1921;const double target=times[index];
   const auto dir=out/("reference-"+std::to_string(k))/("point-"+std::to_string(index));std::filesystem::create_directories(dir);
   auto it=std::find_if(histories[k].begin(),histories[k].end(),[&](const auto& s){return target<=s.right;});require(it!=histories[k].end(),"observation outside history");
   if(target==it->left||target==it->right){const State& y=target==it->left?it->before:it->after;require(y==observations[k][index],"exact reference endpoint mismatch");std::ofstream f(dir/"endpoint.txt");f<<std::setprecision(17)<<target<<'\t'<<values(y)<<'\n';return;}
   std::ofstream rhs(dir/"local-rhs.tsv"),diag(dir/"diagnostics.tsv"),mapping(dir/"maps.tsv");
   rhs<<std::setprecision(17)<<"context\tt_s\tcell\tx\te\tmu\tdx\tde\tdmu\n";
   mapping<<std::setprecision(17)<<"context\tleft\tright\th\tx\te\tmu\treplay_x\treplay_e\treplay_mu\tsteps\trejected\n";
   TraceContext context(p,1,rhs,diag);auto y=it->before;double t=it->left,h=target-t;std::size_t count=0;
   Stepper solver([&](double tt,const State& yy,State& d,int cell){context.SmoothDerivative(tt,yy,d,cell);},context.Breakpoints(),gsl_odeiv2_step_rk8pd,tighter.r,tighter.a,(target-it->left)/4);
   while(t<target){require(++count<=1000,"interpolation work bound");solver.Advance(t,y,h,target);}
   mapping<<1<<'\t'<<it->left<<'\t'<<target<<'\t'<<target-it->left<<'\t'<<values(observations[k][index])<<'\t'<<values(y)<<'\t'<<count<<'\t'<<solver.Stats().rejected<<'\n';
   rhs.flush();mapping.flush();require(bool(rhs)&&bool(mapping),"interpolation trace write failure");
   require(y==observations[k][index],"deterministic reference interpolation differs");
  });
  std::cout<<"PASS all 9605 reference observations reproduced exactly\n";return 0;
 }catch(const std::exception& e){std::cerr<<"STOP "<<e.what()<<'\n';return 1;}
}
