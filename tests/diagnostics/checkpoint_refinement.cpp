// Bounded local-IVP study only: no call to native_main or a t0 trajectory.
#include "../native/pair.hpp"
#define main historical_breakpoint_qualification_main
#include "../bnv/breakpoint_qualification.cpp"
#undef main
#include "checkpoint_evidence.hpp"
#include "retained_trace_helpers.inc"

std::vector<BNV::AcceptedStepRecord> History(const std::filesystem::path& path) {
 std::ifstream f(path);require(bool(f),"missing retained main history");std::string line;std::getline(f,line);
 std::vector<BNV::AcceptedStepRecord> result;State left{{0,0,0}};
 while(std::getline(f,line)) {auto v=Fields(line);require(v.size()==9,"main history width");
  BNV::AcceptedStepRecord s;s.ordinal=std::stoul(v[0]);s.t_left_s=std::stod(v[1]);s.t_right_s=std::stod(v[2]);
  s.accepted_step_s=std::stod(v[3]);s.suggested_next_h_s=std::stod(v[4]);s.cumulative_rejected=std::stoul(v[5]);
  s.y_left=left;for(int k=0;k<3;++k)s.y_right[k]=std::stod(v[k+6]);left=s.y_right;
  require(s.ordinal==result.size()+1&&s.t_right_s>s.t_left_s,"invalid retained history");result.push_back(s);
 }require(result.size()==12&&result.back().t_right_s==Year,"unexpected source history");return result;
}

int main(int argc,char** argv) {std::filesystem::path output;
 try {
  gsl_set_error_handler_off();std::cout<<std::setprecision(17)<<std::unitbuf;
  require(argc==11,"profile thermal entry frozen coefficients archive sha fresh-work retained-source output");
  output=argv[10];require(!std::filesystem::exists(output),"fresh refinement output required");std::filesystem::create_directories(output);
  Inputs in;in.profile=argv[1];in.certificate=in.profile/"certificate.txt";in.thermal=argv[2];in.entry=argv[3];in.frozen=argv[4];in.coefficients=argv[5];in.work=argv[8];
  auto sequence=Core::StellarArchive::Load(argv[6],argv[7]);auto pair=NativePair(in,sequence);
  const std::filesystem::path retained=argv[9];const auto history=History(retained/"tier-0-cap-0/main.accepted_states.tsv");
  std::ofstream summary(output/"candidates.tsv"),witness(output/"witnesses.tsv"),steps(output/"witness-steps.tsv"),rhs(output/"local-rhs.tsv"),diagnostics(output/"local-diagnostics.tsv"),physics(output/"chemistry.tsv");
  for(auto* f:{&summary,&witness,&steps,&rhs,&diagnostics,&physics})*f<<std::setprecision(17);
  summary<<"index\tt_s\tatol\tx\te\tmu\tEchem\tU\tF\tutilization\tstatus\tcontexts\n";
  witness<<"index\tmethod\tsubdivisions\tx\te\tmu\tEchem\taccepted\trejected\trhs\n";
  steps<<"index\tmethod\tsubdivisions\tleft\tright\tx\te\tmu\terr_x\terr_e\terr_mu\n";
  rhs<<"context\tt_s\tcell\tx\te\tmu\tdx\tde\tdmu\n";
  diagnostics<<"context\tt_s\tx\te\tmu";for(const auto& kv:CheckpointEvidence::Values(BNV::BnvDiagnostics{}))diagnostics<<'\t'<<kv.first;diagnostics<<'\n';
  physics<<"index\tt_s\tx\te\tmu\tTinf\tsigma_e\tsigma_mu\tsource_e\tsource_mu\tbeta_e\tbeta_mu\tEchem\n";
  std::size_t count=0;BNV::CheckpointContextFactory factory=[&]{return std::make_unique<TraceContext>(pair.source,++count,rhs,diagnostics);};
  std::vector<double> times;for(int i:{1,2,4,8,16,32,64,128,256,512,960,1440,1919})times.push_back(Year*i/1920);
  BNV::PassiveObservationSchedule schedule(times,"predeclared-reconstruction-study");schedule.Replay({{0,0,0}},history);
  ActualCheckpointContext evaluation(pair.source,0);
  for(auto bracket:schedule.Brackets()) {
   bracket.observation_index=std::llround(bracket.t_observation_s/Year*1920);
   const auto& step=history.at(bracket.right_accepted_ordinal-1);const double t=bracket.t_observation_s,span=t-step.t_left_s;
   std::cout<<"STUDY index="<<bracket.observation_index<<" time="<<t<<" left="<<step.t_left_s<<'\n';
   std::vector<std::map<std::string,double>> refs;for(int i=0;i<5;++i)refs.push_back(Reference(retained/("reference-"+std::to_string(i)+"-observations.tsv"),t));
   for(double atol:{1e-24,1e-26,5e-27,1e-27,1e-28,1e-29}) {
    BNV::Rk8pdCheckpointReconstructor reconstructor;
    // Legal only because the object and underlying configuration are nonconst.
    // Test parameter injection preserves the authenticated library and every formula.
    auto& configuration=const_cast<BNV::Rk8pdConfiguration&>(reconstructor.O2());
    configuration.absolute[1]=configuration.absolute[2]=atol;configuration.Validate();
    const auto before=count;auto q=reconstructor.Reconstruct(bracket,history,factory);
    if(q.status==BNV::CheckpointStatus::Qualified)reconstructor.EvaluateDiagnostics(q,factory);
    std::ostringstream name;name<<std::scientific<<std::setprecision(0)<<atol;
    BNV::Rkf45Configuration config;config.relative=tiers[0].r;config.absolute=tiers[0].a;
    CheckpointEvidence::Write(output/("point-"+std::to_string(bracket.observation_index))/("atol-"+name.str()),q,step,config,refs);
    require(q.diagnostic_U_O.count("Echem"),"state refusal before Echem evaluation");
    summary<<bracket.observation_index<<'\t'<<t<<'\t'<<atol<<'\t'<<values(q.O2_state)<<'\t'<<q.O2_diagnostics.Echem_MeV<<'\t'<<q.diagnostic_U_O.at("Echem")<<'\t'<<q.diagnostic_F_P.at("Echem")<<'\t'<<q.diagnostic_U_O.at("Echem")/(.20*q.diagnostic_F_P.at("Echem"))<<'\t'<<(q.status==BNV::CheckpointStatus::Qualified)<<'\t'<<count-before<<'\n';summary.flush();
    if(atol==1e-24) {auto d=q.O2_diagnostics;
     physics<<bracket.observation_index<<'\t'<<t<<'\t'<<values(q.O2_state)<<'\t'<<d.Tinf_K<<'\t'<<d.sigma_count_s[0]<<'\t'<<d.sigma_count_s[1]<<'\t'<<d.eta_dot_from_sigma_MeV_s[0]<<'\t'<<d.eta_dot_from_sigma_MeV_s[1]<<'\t'<<d.eta_dot_from_beta_MeV_s[0]<<'\t'<<d.eta_dot_from_beta_MeV_s[1]<<'\t'<<d.Echem_MeV<<'\n';physics.flush();
    }
   }
   for(int method=0;method<2;++method)for(int n:{4,8,16}) {
    TraceContext ctx(pair.source,++count,rhs,diagnostics);
    Stepper solver([&](double tt,const State& y,State& d,int cell){ctx.SmoothDerivative(tt,y,d,cell);},ctx.Breakpoints(),method?gsl_odeiv2_step_rkf45:gsl_odeiv2_step_rk8pd,1e-15,{{1e-20,1e-30,1e-30}},span/n);
    auto y=step.y_left;double tt=step.t_left_s,h=span;std::size_t accepted=0;
    while(tt<t) {require(++accepted<=1000,"witness work bound");const auto s=solver.Advance(tt,y,h,t);
     steps<<bracket.observation_index<<'\t'<<method<<'\t'<<n<<'\t'<<s.left<<'\t'<<s.right<<'\t'<<values(y)<<'\t'<<values(s.rk_error)<<'\n';}
    const auto d=evaluation.EvaluateDiagnostics(t,y);auto st=solver.Stats();
    witness<<bracket.observation_index<<'\t'<<method<<'\t'<<n<<'\t'<<values(y)<<'\t'<<d.Echem_MeV<<'\t'<<accepted<<'\t'<<st.rejected<<'\t'<<st.rhs<<'\n';witness.flush();steps.flush();
    require(st.rejected<=.1*(accepted+st.rejected),"witness rejection ceiling");
   }
   rhs.flush();diagnostics.flush();
  }
  for(auto* f:{&summary,&witness,&steps,&rhs,&diagnostics,&physics}){f->flush();require(bool(*f),"study evidence write failure");}
  std::cout<<"COLLECTION COMPLETE; feasibility requires independent offline assessment; no Stage A executed\n";return 0;
 }catch(const std::exception& e){if(!output.empty()&&std::filesystem::exists(output)){std::ofstream f(output/"stop.txt");f<<e.what()<<'\n';}std::cerr<<"STOP "<<e.what()<<'\n';return 1;}
}
