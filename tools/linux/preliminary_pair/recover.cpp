// Isolated deterministic reconstructions from immutable accepted states only.
#define NATIVE_STAGE_NO_MAIN
#include "../../../tests/native/stage.cpp"
#undef NATIVE_STAGE_NO_MAIN
#include "retained_helpers.inc"

std::string bytes(const std::filesystem::path& p){std::ifstream f(p);require(bool(f),"missing comparison file");return {std::istreambuf_iterator<char>(f),{}};}
int main(int argc,char** argv) {
 try {
  gsl_set_error_handler_off();std::cout<<std::setprecision(17)<<std::unitbuf;
  require(argc==12,"profile thermal entry frozen coefficients archive sha fresh-work retained-main output indices");
  Inputs in;in.profile=argv[1];in.certificate=in.profile/"certificate.txt";in.thermal=argv[2];in.entry=argv[3];in.frozen=argv[4];in.coefficients=argv[5];in.work=argv[8];
  const std::filesystem::path retained=argv[9],out=argv[10];
  require(!std::filesystem::exists(out),"fresh recovery output required");std::filesystem::create_directories(out);
  auto sequence=Core::StellarArchive::Load(argv[6],argv[7]);auto pair=NativePair(in,sequence);
  auto main=ReadRetainedMain(retained);
  std::vector<double> times;for(int i=0;i<=1920;++i)times.push_back(Year*i/1920);
  BNV::PassiveObservationSchedule schedule(times,"recovered-local-evidence");schedule.Replay(main.initial_state,main.accepted_steps);
  std::ifstream index_file(argv[11]);std::vector<std::size_t> indices;std::size_t index;
  while(index_file>>index){require(index<=1920,"index out of range");indices.push_back(index);}require(!indices.empty(),"empty reconstruction request");
  NativeQualification::ProcessPool(indices.size(),24,out/"workers",[&](std::size_t task,std::size_t){
   const auto i=indices[task];const auto bracket=schedule.Brackets().at(i);
   const auto dir=out/("point-"+std::to_string(i));std::filesystem::create_directories(dir);
   std::ofstream rhs(dir/"local-rhs.tsv"),diag(dir/"local-diagnostics.tsv"),steps(dir/"witness-steps.tsv"),witness(dir/"witnesses.tsv");
   for(auto* f:{&rhs,&diag,&steps,&witness})*f<<std::setprecision(17);
   rhs<<"context\tt_s\tcell\tx\te\tmu\tdx\tde\tdmu\n";
   diag<<"context\tt_s\tx\te\tmu";for(const auto& kv:CheckpointEvidence::Values(BNV::BnvDiagnostics{}))diag<<'\t'<<kv.first;diag<<'\n';
   steps<<"index\tmethod\tsubdivisions\tleft\tright\tx\te\tmu\terr_x\terr_e\terr_mu\n";
   witness<<"index\tmethod\tsubdivisions\tcontext\tx\te\tmu\taccepted\trejected\trhs\n";
   std::size_t count=0;BNV::CheckpointContextFactory factory=[&]{return std::make_unique<TraceContext>(pair.source,++count,rhs,diag);};
   BNV::Rk8pdCheckpointReconstructor reconstruction;auto q=reconstruction.Reconstruct(bracket,main.accepted_steps,factory);
   if(q.status==BNV::CheckpointStatus::Qualified)reconstruction.EvaluateDiagnostics(q,factory);
   const auto& step=main.accepted_steps.at(bracket.initial?0:bracket.right_accepted_ordinal-1);
   CheckpointEvidence::Write(dir/"packet",q,step,main.configuration);
   require(q.status==BNV::CheckpointStatus::Qualified&&q.diagnostics_evaluated,"reconstruction refused");
   for(const auto* name:{"checkpoint.json","state.tsv","diagnostics.tsv","failure.txt"})
    require(bytes(dir/"packet"/name)==bytes(retained/"checkpoints-1920"/("point-"+std::to_string(i))/name),(std::string("retained packet differs: ")+name).c_str());
   if(!bracket.initial&&bracket.t_observation_s!=step.t_right_s){
    const double t=bracket.t_observation_s,span=t-step.t_left_s;
    for(int method=0;method<2;++method)for(int n:{4,8,16}){
     TraceContext ctx(pair.source,++count,rhs,diag);
     Stepper solver([&](double tt,const State& y,State& d,int cell){ctx.SmoothDerivative(tt,y,d,cell);},ctx.Breakpoints(),method?gsl_odeiv2_step_rkf45:gsl_odeiv2_step_rk8pd,1e-15,{{1e-20,1e-30,1e-30}},span/n);
     auto y=step.y_left;double tt=step.t_left_s,h=span;std::size_t accepted=0;
     while(tt<t){require(++accepted<=1000,"witness work bound");const auto s=solver.Advance(tt,y,h,t);
      steps<<i<<'\t'<<method<<'\t'<<n<<'\t'<<s.left<<'\t'<<s.right<<'\t'<<values(y)<<'\t'<<values(s.rk_error)<<'\n';}
     const auto stats=solver.Stats();witness<<i<<'\t'<<method<<'\t'<<n<<'\t'<<count<<'\t'<<values(y)<<'\t'<<accepted<<'\t'<<stats.rejected<<'\t'<<stats.rhs<<'\n';
     require(stats.rejected<=.1*(accepted+stats.rejected),"witness rejection ceiling");
    }
   }
   for(auto* f:{&rhs,&diag,&steps,&witness}){f->flush();require(bool(*f),"evidence write failure");}
   std::ofstream done(dir/"consistency.txt");done<<"PASS four packet files byte-identical; no main integration\n";
  });
  std::cout<<"PASS reconstructed "<<indices.size()<<" exact retained packets and collected local witnesses\n";return 0;
 }catch(const std::exception& e){std::cerr<<"STOP "<<e.what()<<'\n';return 1;}
}
