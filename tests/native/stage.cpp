// Bounded actual-path evidence collection. Compile now; execution requires separate owner authorization.
#include "pair.hpp"
#define main historical_breakpoint_qualification_main
#include "../bnv/breakpoint_qualification.cpp"
#undef main

BNV::MainTrajectoryResult native_main(Ctx p,double a,State y,double b,Tol tol,double cap,
 const std::vector<double>& times,const std::filesystem::path& out) {
    std::filesystem::create_directories(out);
    RunState state(p->OrdinaryContext());
    auto driver=std::make_shared<BNV::ControlledBnvSecularDriver>(p,BNV::ControlledBnvSecularDriver::Mode::Coupled);
    EV::EvolutionSystem system(state.ctx,state.state,state.rhs,state.layout,{driver});
    BNV::PassiveObservationSchedule schedule(times,"RF-native-passive");
    BNV::Rkf45Configuration config;config.relative=tol.r;config.absolute=tol.a;config.initial_step_s=1;config.maximum_step_s=cap;
    BNV::UninterruptedBnvTrajectory solver(system,state.layout,p->OrdinaryContext(),config,
      {p->RunCardIdentity(),"native snapshot mathematical source","full/cheap unchanged guards",schedule.Identity()});
    auto result=solver.Integrate(a,b,y,schedule,true);
    WriteMainEvidence(out,result,schedule,0);
    std::ofstream events(out/"events.tsv");events<<std::setprecision(17)<<"ordinal\ttime\tknot\tdirection\tresidual\tlow\thigh\tx_error\te_error\tmu_error\n";
    for(const auto& s:result.accepted_steps)events<<s.ordinal<<'\t'<<s.event.time<<'\t'<<s.event.knot<<'\t'<<s.event.direction<<'\t'<<s.event.residual<<'\t'<<s.event.time_low<<'\t'<<s.event.time_high<<'\t'<<values(s.local_rk_error)<<'\n';
    require(result.final_time_s==b,"incomplete bounded endpoint");
    const auto& st=result.statistics;
    require(st.rejected_steps<=.1*(st.accepted_steps+st.rejected_steps),"bounded rejection ceiling exceeded");
    return result;
}
void stage_role(Ctx p,double end,const std::filesystem::path& out) {
    std::filesystem::create_directories(out);
    std::vector<double> times;for(int i=0;i<=1920;++i)times.push_back(end*i/1920);
    // Five independent, original-t0 references; retain their complete histories.
    const double H=end/256;std::vector<Solution> refs;
    for(int i=0;i<5;++i) {
      const double caps[]={H,H/2,H/4,H/4,H/8};
      refs.push_back(integrate(p,0,{{0,0,0}},end,1,i==0?tight:tighter,i<3?gsl_odeiv2_step_rk8pd:gsl_odeiv2_step_rkf45,caps[i]));
      const auto& ref=refs.back();std::ofstream f(out/("reference-"+std::to_string(i)+"-history.tsv"));
      f<<std::setprecision(17)<<"left\tright\tx\te\tmu\n";
      for(const auto& step:ref.steps)f<<step.left<<'\t'<<step.right<<'\t'<<values(step.after)<<'\n';
      require(ref.stats.rejected<=.1*(ref.steps.size()+ref.stats.rejected),"reference rejection ceiling exceeded");
      ActualCheckpointContext evaluation(p,100+i);
      std::ofstream obs(out/("reference-"+std::to_string(i)+"-observations.tsv"));
      obs<<std::setprecision(17)<<"t\tx\te\tmu";
      for(const auto& kv:DiagnosticValues(evaluation.EvaluateDiagnostics(0,{{0,0,0}})))obs<<'\t'<<kv.first;obs<<'\n';
      for(double t:times) {
        auto y=observe(p,ref,t);auto d=evaluation.EvaluateDiagnostics(t,y);
        require(d.valid_through_sample&&d.Tinf_K>=1e6&&d.Tinf_K<=1e9,"reference scientific domain refusal");
        obs<<t<<'\t'<<values(y);for(const auto& kv:DiagnosticValues(d))obs<<'\t'<<kv.second;obs<<'\n';
      }
    }
    const double caps[]={INFINITY,38723807757.0,19361903878.5,9680951939.25};
    for(int tier=0;tier<3;++tier)for(int policy=0;policy<4;++policy) {
      auto dir=out/("tier-"+std::to_string(tier)+"-cap-"+std::to_string(policy));
      auto main=native_main(p,0,{{0,0,0}},end,tiers[tier],caps[policy],times,dir);
      auto repeat=native_main(p,0,{{0,0,0}},end,tiers[tier],caps[policy],{0,end},dir/"endpoint-schedule-repeat");
      require(main.final_state==repeat.final_state&&main.accepted_steps.size()==repeat.accepted_steps.size(),"schedule changed native main");
      for(std::size_t i=0;i<main.accepted_steps.size();++i) {
        const auto& a=main.accepted_steps[i];const auto& b=repeat.accepted_steps[i];
        require(a.y_right==b.y_right&&a.t_right_s==b.t_right_s&&a.suggested_next_h_s==b.suggested_next_h_s&&a.cumulative_rejected==b.cumulative_rejected,"schedule changed accepted history");
      }
      auto count=std::make_shared<std::size_t>(0);
      BNV::CheckpointContextFactory factory=[p,count]{return std::make_unique<ActualCheckpointContext>(p,++*count);};
      BNV::Rk8pdCheckpointReconstructor reconstruction;
      for(const auto& step:main.accepted_steps)if(step.event.knot>=0) {
        BNV::ObservationBracket exact{0,step.t_right_s,step.ordinal-1,step.ordinal,step.t_left_s,step.t_right_s,false};
        require(reconstruction.Reconstruct(exact,main.accepted_steps,{}).state==step.y_right,"event endpoint bypass changed state");
        exact.t_observation_s=std::nextafter(step.t_right_s,step.t_left_s);
        require(reconstruction.Reconstruct(exact,main.accepted_steps,factory).status==BNV::CheckpointStatus::Qualified,"left event neighbor unresolved");
        if(step.ordinal<main.accepted_steps.size()) {
          const auto& next=main.accepted_steps.at(step.ordinal);
          BNV::ObservationBracket right{0,std::nextafter(step.t_right_s,next.t_right_s),step.ordinal,next.ordinal,next.t_left_s,next.t_right_s,false};
          require(reconstruction.Reconstruct(right,main.accepted_steps,factory).status==BNV::CheckpointStatus::Qualified,"right event neighbor unresolved");
        }
      }
      std::ofstream r20(dir/"R20-grids.tsv");r20<<std::setprecision(17)<<"intervals\tresidual\tN\treconstruction_uncertainty\tdelta_U\tdelta_Echem\toutgoing_integral\n";
      for(int n:{240,480,960,1920}) {
        std::vector<double> grid;for(int j=0;j<=n;++j)grid.push_back(end*j/n);
        BNV::PassiveObservationSchedule schedule(grid,"native-R20-grid");schedule.Replay(main.initial_state,main.accepted_steps);
        std::vector<BNV::CheckpointOutput> outputs;
        BNV::ControlledBnvSecularDriver driver(p,BNV::ControlledBnvSecularDriver::Mode::Coupled);
        Capture diagnostics(dir/("diagnostics-"+std::to_string(n)+".tsv"),driver);
        std::ofstream errors(dir/("uncertainties-"+std::to_string(n)+".tsv"));
        errors<<std::setprecision(17)<<"t\tx\te\tmu\tUx\tUe\tUmu\tEx\tEe\tEmu\n";
        for(const auto& bracket:schedule.Brackets()) {
          auto q=reconstruction.Reconstruct(bracket,main.accepted_steps,factory);
          require(q.status==BNV::CheckpointStatus::Qualified,"native checkpoint state unresolved");
          reconstruction.EvaluateDiagnostics(q,factory);
          require(q.status==BNV::CheckpointStatus::Qualified&&q.diagnostics_evaluated,"native checkpoint diagnostics unresolved");
          diagnostics.SaveSnapshot(q.state[0],q.diagnostics);
          errors<<q.t_observation_s<<'\t'<<values(q.state)<<'\t'<<values(q.reconstruction_error_contribution)<<'\t'<<values(q.qualification.U_event)<<'\n';
          outputs.push_back(std::move(q));
        }
        diagnostics.Publish();auto r=BNV::ComputeCheckpointR20(outputs);
        r20<<n<<'\t'<<r.R20_residual_erg<<'\t'<<r.N_R20_erg<<'\t'<<r.reconstruction_uncertainty_erg<<'\t'<<r.delta_Uth_erg<<'\t'<<r.delta_Echem_erg<<'\t'<<r.outgoing_integral_erg<<'\n';
      }
      if(policy==3) {
        auto it=std::find_if(main.accepted_steps.begin(),main.accepted_steps.end(),[&](const auto& s){return s.t_right_s>=end/2;});
        require(it!=main.accepted_steps.end()&&it->t_right_s<end,"no interior restart endpoint");
        const auto checkpoint=dir/"restart.txt";
        {std::ofstream f(checkpoint);f<<std::setprecision(17)<<std::quoted(p->RunCardIdentity())<<'\n'<<it->t_right_s<<'\t'<<values(it->y_right)<<'\n'<<it->ordinal<<'\t'<<it->cumulative_rejected<<'\n';}
        std::ifstream f(checkpoint);std::string role;double t;State y;std::size_t prefix,rejected;
        f>>std::quoted(role)>>t;for(auto& v:y)f>>v;f>>prefix>>rejected;
        require(bool(f)&&role==p->RunCardIdentity()&&prefix==it->ordinal&&rejected==it->cumulative_rejected,"restart prefix/identity mismatch");
        auto tail=native_main(p,t,y,end,tiers[tier],caps[policy],{t,end},dir/"restart-tail");
        std::ofstream result(dir/"restart-comparison.tsv");result<<std::setprecision(17)<<"x\te\tmu\n"<<values(tail.final_state)<<'\n'<<values(main.final_state)<<'\n';
        // Full original prefix remains in parent main.accepted_states.tsv; no claim of bitwise step continuation.
      }
    }
}
int main(int argc,char** argv) {
 try {
    gsl_set_error_handler_off();std::cout<<std::setprecision(17)<<std::unitbuf;
    require(argc==11,"profile thermal entry frozen coefficients archive sha fresh-work years output");
    const int years=std::stoi(argv[9]);require(years==1||years==2048,"only separately authorized 1/2048-year bounds supported");
    Inputs inputs;inputs.profile=argv[1];inputs.certificate=inputs.profile/"certificate.txt";inputs.thermal=argv[2];inputs.entry=argv[3];inputs.frozen=argv[4];inputs.coefficients=argv[5];inputs.work=argv[8];
    require(!std::filesystem::exists(argv[10]),"fresh stage output required");
    auto sequence=Core::StellarArchive::Load(argv[6],argv[7]);auto pair=NativePair(inputs,sequence);
    stage_role(pair.source,years*Year,std::filesystem::path(argv[10])/"source");
    stage_role(pair.control,years*Year,std::filesystem::path(argv[10])/"control");
    std::cout<<"COLLECTION COMPLETE; nonlinear acceptance requires all predeclared error terms and paired assessment.\n";
    return 0;
 } catch(const std::exception& e){std::cerr<<"STOP "<<e.what()<<'\n';return 1;}
}
