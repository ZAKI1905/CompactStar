// Retained-history reconstruction only: no main integration, reference trajectory,
// restart tail or stellar solve. Exercise the exact shared grid collector.
#define main stage_entry_not_invoked
#include "stage.cpp"
#undef main

BNV::MainTrajectoryResult ReadRetainedMain(const std::filesystem::path& dir) {
 BNV::MainTrajectoryResult main;main.final_time_s=Year;main.configuration.relative=tiers[0].r;
 main.configuration.absolute=tiers[0].a;main.configuration.initial_step_s=1;main.configuration.maximum_step_s=INFINITY;
 std::ifstream states(dir/"main.accepted_states.tsv"),events(dir/"events.tsv"),internal(dir/"main.internal_steps.tsv");
 require(bool(states)&&bool(events)&&bool(internal),"retained history missing");
 std::string line,event,detail;std::getline(states,line);std::getline(events,event);std::getline(internal,detail);
 State left{{0,0,0}};
 while(std::getline(states,line)) {
  BNV::AcceptedStepRecord s;std::istringstream r(line);
  r>>s.ordinal>>s.t_left_s>>s.t_right_s>>s.accepted_step_s>>s.suggested_next_h_s>>s.cumulative_rejected;
  for(auto& v:s.y_right)r>>v;s.y_left=left;left=s.y_right;
  require(bool(r)&&s.ordinal==main.accepted_steps.size()+1,"bad retained history");
  require(bool(std::getline(events,event))&&bool(std::getline(internal,detail)),"missing retained metadata");
  std::size_t ordinal;std::istringstream e(event);e>>ordinal>>s.event.time>>s.event.knot>>s.event.direction>>s.event.residual>>s.event.time_low>>s.event.time_high;
  for(auto& v:s.local_rk_error)e>>v;require(bool(e)&&ordinal==s.ordinal&&s.event.knot==-1,"unexpected retained event");
  std::istringstream d(detail);double ceiling,t0,t1,h,next,x0,x1;std::size_t cp,rejected;
  d>>ordinal>>cp>>ceiling>>t0>>t1>>h>>next>>rejected>>x0>>x1>>s.cstar_cell_left>>s.cstar_cell_right>>s.cstar_knots_crossed;
  require(bool(d)&&ordinal==s.ordinal&&t0==s.t_left_s&&t1==s.t_right_s,"retained internal mismatch");
  main.accepted_steps.push_back(s);
 }
 require(main.accepted_steps.size()==12&&main.accepted_steps.back().t_right_s==Year,"unexpected retained main");
 main.final_state=main.accepted_steps.back().y_right;return main;
}
void WriteLedger(const std::filesystem::path& dir,const BNV::CheckpointR20Result& r) {
 std::ofstream f(dir/"R20-grids.tsv");f<<std::setprecision(17)<<"intervals\tresidual\tN\treconstruction_uncertainty\tdelta_U\tdelta_Echem\toutgoing_integral\n"
 <<240<<'\t'<<r.R20_residual_erg<<'\t'<<r.N_R20_erg<<'\t'<<r.reconstruction_uncertainty_erg<<'\t'<<r.delta_Uth_erg<<'\t'<<r.delta_Echem_erg<<'\t'<<r.outgoing_integral_erg<<'\n';
}
int main(int argc,char** argv) {
 try {
  gsl_set_error_handler_off();std::cout<<std::unitbuf;
  require(argc==11,"profile thermal entry frozen coefficients archive sha fresh-work retained-main output");
  require(NativeQualification::PhysicalCpus().size()>=24,"24-core regression allocation required");
  Inputs inputs;inputs.profile=argv[1];inputs.certificate=inputs.profile/"certificate.txt";inputs.thermal=argv[2];inputs.entry=argv[3];inputs.frozen=argv[4];inputs.coefficients=argv[5];inputs.work=argv[8];
  const std::filesystem::path out=argv[10];require(!std::filesystem::exists(out),"fresh regression output required");
  auto sequence=Core::StellarArchive::Load(argv[6],argv[7]);auto pair=NativePair(inputs,sequence);
  const auto main=ReadRetainedMain(argv[9]);
  for(std::size_t workers:{1,24}) {
   const auto dir=out/("workers-"+std::to_string(workers));std::filesystem::create_directories(dir);
   auto count=std::make_shared<std::size_t>(0);
   BNV::CheckpointContextFactory factory=[&]{return std::make_unique<ActualCheckpointContext>(pair.source,++*count);};
   const auto start=std::chrono::steady_clock::now();auto r=CollectGrid(pair.source,main,240,dir,factory,workers);
   WriteLedger(dir,r);
   std::ofstream timing(dir/"timing.txt");timing<<std::setprecision(17)<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
   std::cout<<"COLLECTED retained 240-grid with workers="<<workers<<'\n';
  }
  std::cout<<"PASS retained-history reconstruction collected; external exact-byte comparison required; zero main/stellar solves\n";return 0;
 }catch(const std::exception& e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}
}
