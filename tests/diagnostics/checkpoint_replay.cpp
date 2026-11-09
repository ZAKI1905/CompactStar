// One source-only diagnostic replay; preserves every original snapshot-bound file.
#include "../native/pair.hpp"
#define main historical_breakpoint_qualification_main
#include "../bnv/breakpoint_qualification.cpp"
#undef main
#include "checkpoint_evidence.hpp"
#include "retained_native_main.inc" // verbatim function extracted from authenticated stage.cpp

struct TraceContext final:BNV::ICheckpointEvaluationContext {
 ActualCheckpointContext context;std::size_t id;std::ostream& rhs;std::ostream& diagnostics;
 TraceContext(Ctx p,std::size_t n,std::ostream& r,std::ostream& d):context(p,n),id(n),rhs(r),diagnostics(d){}
 void RequireFullCurrent()const override{context.RequireFullCurrent();}
 std::string Identity()const override{return context.Identity();}
 bool HasExactZeroChemicalSolution()const override{return context.HasExactZeroChemicalSolution();}
 std::vector<double> Breakpoints()const override{return context.Breakpoints();}
 void Derivative(double t,const State& y,State& out)override{context.Derivative(t,y,out);Record(t,y,out,-1);}
 void SmoothDerivative(double t,const State& y,State& out,int cell)override{context.SmoothDerivative(t,y,out,cell);Record(t,y,out,cell);}
 void Record(double t,const State& y,const State& d,int cell){rhs<<id<<'\t'<<t<<'\t'<<cell<<'\t'<<values(y)<<'\t'<<values(d)<<'\n';}
 BNV::BnvDiagnostics EvaluateDiagnostics(double t,const State& y)override{
  auto d=context.EvaluateDiagnostics(t,y);diagnostics<<id<<'\t'<<t<<'\t'<<values(y);
  for(const auto& kv:CheckpointEvidence::Values(d))diagnostics<<'\t'<<kv.second;
  diagnostics<<'\n';diagnostics.flush();return d;
 }
};
std::string Bytes(const std::filesystem::path& p){std::ifstream f(p);require(bool(f),"missing retained evidence");return {std::istreambuf_iterator<char>(f),{}};}
std::vector<std::string> Fields(const std::string& line){std::stringstream s(line);std::string v;std::vector<std::string> r;while(std::getline(s,v,'\t'))r.push_back(v);return r;}
std::map<std::string,double> Reference(const std::filesystem::path& p,double time){
 std::ifstream f(p);require(bool(f),"missing retained reference");std::string line;std::getline(f,line);auto names=Fields(line);
 while(std::getline(f,line)){auto fields=Fields(line);require(fields.size()==names.size(),"reference row width");if(std::stod(fields[0])!=time)continue;
 std::map<std::string,double> row;for(std::size_t i=0;i<names.size();++i)row[names[i]]=std::stod(fields[i]);return row;}
 throw std::runtime_error("reference does not contain exact checkpoint time");
}
int main(int argc,char** argv){
 std::filesystem::path output;
 try{
  gsl_set_error_handler_off();std::cout<<std::setprecision(17)<<std::unitbuf;
  require(argc==11,"profile thermal entry frozen coefficients archive sha fresh-work retained-source output");
  output=argv[10];require(!std::filesystem::exists(output),"fresh diagnostic output required");std::filesystem::create_directories(output);
  Inputs in;in.profile=argv[1];in.certificate=in.profile/"certificate.txt";in.thermal=argv[2];in.entry=argv[3];in.frozen=argv[4];in.coefficients=argv[5];in.work=argv[8];
  auto sequence=Core::StellarArchive::Load(argv[6],argv[7]);auto pair=NativePair(in,sequence);
  const std::filesystem::path retained=argv[9];
  std::vector<double> grid;for(int i=0;i<=1920;++i)grid.push_back(Year*i/1920);
  auto main=native_main(pair.source,0,{{0,0,0}},Year,tiers[0],INFINITY,grid,output/"main");
  require(Bytes(output/"main/main.accepted_states.tsv")==Bytes(retained/"tier-0-cap-0/main.accepted_states.tsv"),"replayed main differs from original accepted states");
  std::ofstream rhs(output/"local-rhs.tsv"),diag(output/"local-diagnostics.tsv");rhs<<std::setprecision(17);diag<<std::setprecision(17);
  rhs<<"context\tt_s\tcell\tx\te\tmu\tdx\tde\tdmu\n";
  diag<<"context\tt_s\tx\te\tmu";for(const auto& kv:CheckpointEvidence::Values(BNV::BnvDiagnostics{}))diag<<'\t'<<kv.first;diag<<'\n';
  std::size_t count=0;BNV::CheckpointContextFactory factory=[&]{return std::make_unique<TraceContext>(pair.source,++count,rhs,diag);};
  BNV::Rk8pdCheckpointReconstructor reconstruction;
  grid.clear();for(int i=0;i<=240;++i)grid.push_back(Year*i/240);
  BNV::PassiveObservationSchedule schedule(grid,"native-R20-grid");schedule.Replay(main.initial_state,main.accepted_steps);
  for(const auto& bracket:schedule.Brackets()){
   require(bracket.observation_index<=1,"diagnostic replay cannot pass first nonzero checkpoint");
   auto q=reconstruction.Reconstruct(bracket,main.accepted_steps,factory);
   if(q.status==BNV::CheckpointStatus::Qualified)reconstruction.EvaluateDiagnostics(q,factory);
   std::vector<std::map<std::string,double>> refs;for(int i=0;i<5;++i)refs.push_back(Reference(retained/("reference-"+std::to_string(i)+"-observations.tsv"),bracket.t_observation_s));
   const auto& step=main.accepted_steps.at(bracket.initial?0:bracket.right_accepted_ordinal-1);
   auto dir=output/("checkpoint-"+std::to_string(bracket.observation_index));
   CheckpointEvidence::Write(dir,q,step,main.configuration,refs);rhs.flush();diag.flush();
   require(bool(rhs)&&bool(diag),"trace write failed");
   if(q.status!=BNV::CheckpointStatus::Qualified||!q.diagnostics_evaluated){
    std::cout<<Bytes(dir/"failure.txt");throw std::runtime_error("diagnostic replay stopped at first unresolved checkpoint; evidence persisted");}
   if(bracket.observation_index==1)throw std::runtime_error("original failure did not reproduce; diagnostic stop, no further checkpoints authorized");
  }
  throw std::runtime_error("diagnostic replay ended without expected checkpoint");
 }catch(const std::exception& e){
  if(!output.empty()&&std::filesystem::exists(output)){std::ofstream f(output/"stop.txt");f<<e.what()<<'\n';f.flush();}
  std::cerr<<"STOP "<<e.what()<<'\n';return 1;
 }
}
