// Bounded actual-path evidence collection. Compile now; execution requires separate owner authorization.
#include "pair.hpp"
#define main historical_breakpoint_qualification_main
#include "../bnv/breakpoint_qualification.cpp"
#undef main
#include "../diagnostics/checkpoint_evidence.hpp"
#include <cstdlib>
#include "restart_schedule.hpp"
#include "process_pool.hpp"
#include "checkpoint_transport.hpp"

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
std::size_t StageWorkers() {
    const char* value=std::getenv("RF_STAGE_WORKERS");
    if(!value||std::string(value)=="1")return 1;
    require(std::string(value)=="24","only 1 or 24 stage workers supported");
    require(NativeQualification::PhysicalCpus().size()>=24,"24 distinct allocated physical cores required");
    return 24;
}
BNV::CheckpointR20Result CollectGrid(Ctx p,const BNV::MainTrajectoryResult& main,int n,
 const std::filesystem::path& dir,const BNV::CheckpointContextFactory& factory,std::size_t workers) {
    BNV::Rk8pdCheckpointReconstructor reconstruction;
        std::vector<double> grid;for(int j=0;j<=n;++j)grid.push_back(main.final_time_s*j/n);
        BNV::PassiveObservationSchedule schedule(grid,"native-R20-grid");schedule.Replay(main.initial_state,main.accepted_steps);
        std::vector<BNV::CheckpointOutput> outputs;
        BNV::ControlledBnvSecularDriver driver(p,BNV::ControlledBnvSecularDriver::Mode::Coupled);
        Capture diagnostics(dir/("diagnostics-"+std::to_string(n)+".tsv"),driver);
        std::ofstream errors(dir/("uncertainties-"+std::to_string(n)+".tsv"));
        errors<<std::setprecision(17)<<"t\tx\te\tmu\tUx\tUe\tUmu\tEx\tEe\tEmu\n";
        auto calculate=[&](std::size_t index,const BNV::CheckpointContextFactory& local_factory) {
          const auto& bracket=schedule.Brackets().at(index);
          auto q=reconstruction.Reconstruct(bracket,main.accepted_steps,local_factory);
          if(q.status==BNV::CheckpointStatus::Qualified)reconstruction.EvaluateDiagnostics(q,local_factory);
          // Serialize the actual packet, including rejected channels, before any
          // refusal. This is evidence only; it never changes qualification.
          const auto& step=main.accepted_steps.at(bracket.initial?0:bracket.right_accepted_ordinal-1);
          auto evidence=dir/("checkpoints-"+std::to_string(n))/("point-"+std::to_string(bracket.observation_index));
          CheckpointEvidence::Write(evidence,q,step,main.configuration);
          if(q.status!=BNV::CheckpointStatus::Qualified||!q.diagnostics_evaluated) {
            if(workers==1){diagnostics.Publish();errors.flush();}
            std::cout<<"FAILED_CHECKPOINT "<<evidence<<'\n';
            std::ifstream failure(evidence/"failure.txt");std::cout<<failure.rdbuf();
          }
          if(workers>1)NativeQualification::Transport::Write(evidence/"packet.bin",q);
          require(q.status==BNV::CheckpointStatus::Qualified&&q.diagnostics_evaluated,"native checkpoint diagnostics unresolved");
          return q;
        };
        auto consume=[&](BNV::CheckpointOutput q) {
          diagnostics.SaveSnapshot(q.state[0],q.diagnostics);
          errors<<q.t_observation_s<<'\t'<<values(q.state)<<'\t'<<values(q.reconstruction_error_contribution)<<'\t'<<values(q.qualification.U_event)<<'\n';
          outputs.push_back(std::move(q));
        };
        if(workers==1) {
          for(std::size_t i=0;i<schedule.Brackets().size();++i)consume(calculate(i,factory));
        } else {
          auto packet=[&](std::size_t i){return dir/("checkpoints-"+std::to_string(n))/("point-"+std::to_string(i))/"packet.bin";};
          try {
            NativeQualification::ProcessPool(schedule.Brackets().size(),workers,
              dir/("workers-"+std::to_string(n)),[&](std::size_t i,std::size_t){
                // Stable provenance independent of completion order; fresh private contexts.
                std::size_t identity=(std::size_t(n)<<32)+(i<<16);
                BNV::CheckpointContextFactory local=[&]{return std::make_unique<ActualCheckpointContext>(p,++identity);};
                calculate(i,local);
              });
          } catch(...) {
            // Preserve every completed packet. No ledger is formed from a partial grid.
            for(std::size_t i=0;i<schedule.Brackets().size();++i)if(std::filesystem::exists(packet(i))) {
              auto q=NativeQualification::Transport::Read(packet(i));
              if(q.status==BNV::CheckpointStatus::Qualified&&q.diagnostics_evaluated)consume(std::move(q));
            }
            diagnostics.Publish();errors.flush();throw;
          }
          for(std::size_t i=0;i<schedule.Brackets().size();++i) {
            auto q=NativeQualification::Transport::Read(packet(i));
            require(q.observation_index==i&&q.t_observation_s==schedule.Brackets()[i].t_observation_s,
                    "worker returned wrong checkpoint identity");
            consume(std::move(q));
          }
        }
        diagnostics.Publish();errors.flush();
        require(bool(errors),"checkpoint uncertainty evidence write failed");
        return BNV::ComputeCheckpointR20(outputs);
}

void stage_role(Ctx p,double end,const std::filesystem::path& out) {
    std::filesystem::create_directories(out);
    std::vector<double> times;for(int i=0;i<=1920;++i)times.push_back(end*i/1920);
    // Five independent, original-t0 references; retain their complete histories.
    const double H=end/256;std::vector<Solution> refs;
    for(int i=0;i<5;++i) {
      std::cout<<"REFERENCE "<<p->RunCardIdentity()<<" index="<<i<<'\n';
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
      std::cout<<"MAIN "<<p->RunCardIdentity()<<" tier="<<tier<<" policy="<<policy<<'\n';
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
        auto r=CollectGrid(p,main,n,dir,factory,StageWorkers());
        r20<<n<<'\t'<<r.R20_residual_erg<<'\t'<<r.N_R20_erg<<'\t'<<r.reconstruction_uncertainty_erg<<'\t'<<r.delta_Uth_erg<<'\t'<<r.delta_Echem_erg<<'\t'<<r.outgoing_integral_erg<<'\n';
        r20.flush();std::cout<<"GRID_COMPLETE tier="<<tier<<" policy="<<policy<<" n="<<n<<'\n';
      }
      if(policy==3) {
        // Stage A coverage policy; Stage B retains its original rule and authorization boundary.
        auto it=end==Year
          ? NativeQualification::RestartAtOrBeforeMidpoint(main.accepted_steps,0.,end)
          : std::find_if(main.accepted_steps.cbegin(),main.accepted_steps.cend(),[&](const auto& s){return s.t_right_s>=end/2;});
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
#ifndef NATIVE_STAGE_NO_MAIN
int main(int argc,char** argv) {
 try {
    gsl_set_error_handler_off();std::cout<<std::setprecision(17)<<std::unitbuf;
    require(argc==11,"profile thermal entry frozen coefficients archive sha fresh-work years output");
    const char* requested_role=std::getenv("RF_STAGE_ROLE");
    require(requested_role!=nullptr,"explicit source/control role required; no automatic paired execution");
    const std::string role=requested_role;
    require(role=="source"||role=="control","invalid thermal role");
    const int years=std::stoi(argv[9]);require(years==1||years==2048,"only separately authorized 1/2048-year bounds supported");
    const auto workers=StageWorkers();
    require(workers==1||(years==1&&role=="source"),"parallel execution is source Stage A only");
    Inputs inputs;inputs.profile=argv[1];inputs.certificate=inputs.profile/"certificate.txt";inputs.thermal=argv[2];inputs.entry=argv[3];inputs.frozen=argv[4];inputs.coefficients=argv[5];inputs.work=argv[8];
    require(!std::filesystem::exists(argv[10]),"fresh stage output required");
    auto sequence=Core::StellarArchive::Load(argv[6],argv[7]);auto pair=NativePair(inputs,sequence);
    stage_role(role=="source"?pair.source:pair.control,years*Year,std::filesystem::path(argv[10])/role);
    std::cout<<"ROLE COLLECTION COMPLETE "<<role<<"; accuracy assessment is required before any matched control.\n";
    return 0;
 } catch(const std::exception& e){std::cerr<<"STOP "<<e.what()<<'\n';return 1;}
}
#endif
