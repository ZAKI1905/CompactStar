#include "breakpoint_fixture.hpp"
#include "checked_evidence.hpp"
using namespace BreakpointTest;
using Stepper=BNV::BreakpointStepper;
struct Tol{double r;State a;};
const Tol tiers[]={{1e-7,{{1e-12,1e-18,1e-18}}},{1e-9,{{1e-14,1e-20,1e-20}}},{1e-11,{{1e-16,1e-22,1e-22}}}};
const Tol tight{1e-14,{{1e-19,1e-25,1e-25}}};
const Tol tighter{1e-15,{{1e-20,1e-26,1e-26}}};
std::string values(const State& y){std::ostringstream s;s<<std::setprecision(17)<<y[0]<<'\t'<<y[1]<<'\t'<<y[2];return s.str();}
std::map<std::string,double> DiagnosticValues(const BNV::BnvDiagnostics& d)
{
    return {
      {"P_dir_eq",d.P_dir_eq_erg_s},{"P_dir_actual",d.P_dir_actual_erg_s},
      {"LH",d.LH_erg_s},{"DeltaLnu",d.DeltaLnu_erg_s},{"DeltaPbeta",d.DeltaPbeta_erg_s},
      {"Lnu_eq",d.Lnu_eq_erg_s},{"Lnu_full",d.Lnu_full_erg_s},{"Lgamma",d.Lgamma_erg_s},
      {"Lother",d.Lother_erg_s},{"Pnet",d.Pnet_erg_s},{"mu_B",d.mu_B_inf_MeV},
      {"mu_n_actual",d.mu_n_actual_inf_MeV},{"sigma_e",d.sigma_count_s[0]},
      {"sigma_mu",d.sigma_count_s[1]},{"Echem",d.Echem_MeV},
      {"Cstar",d.Cstar_erg_K},{"Tinf",d.Tinf_K},{"Lout",d.L_out_fluid_inf_erg_s}
    };
}

struct Solution{State y;std::vector<Stepper::Step> steps;Stepper::Statistics stats;};
Solution integrate(Ctx p,double a,State y,double b,double h,Tol tol,const gsl_odeiv2_step_type* method,double cap=INFINITY)
{
    ActualCheckpointContext ctx(p,1);auto knots=ctx.Breakpoints();
    Stepper solver([&](double t,const State& in,State& out,int cell){ctx.SmoothDerivative(t,in,out,cell);},knots,method,tol.r,tol.a,cap);
    double t=a;Solution result;
    while(t<b){require(result.steps.size()<100000,"qualification work limit");result.steps.push_back(solver.Advance(t,y,h,b));}
    result.y=y;result.stats=solver.Stats();return result;
}
State observe(Ctx p,const Solution& s,double t)
{
    if(t==s.steps.front().left)return s.steps.front().before;
    for(const auto& r:s.steps)if(t<=r.right){if(t==r.right)return r.after;
      return integrate(p,r.left,r.before,t,t-r.left,tighter,gsl_odeiv2_step_rk8pd,(t-r.left)/4).y;}
    throw std::runtime_error("missing reference observation");
}
BNV::MainTrajectoryResult main_run(Ctx p,double a,State y,double b,double h,Tol tol,const std::vector<double>& times,bool historical=false)
{
    RunState state(p->OrdinaryContext());
    auto driver=std::make_shared<BNV::ControlledBnvSecularDriver>(p,BNV::ControlledBnvSecularDriver::Mode::Coupled);
    EV::EvolutionSystem system(state.ctx,state.state,state.rhs,state.layout,{driver});
    BNV::PassiveObservationSchedule schedule(times,"bounded-passive");
    BNV::Rkf45Configuration config;config.relative=tol.r;config.absolute=tol.a;config.initial_step_s=h;config.breakpoint_aware=!historical;
    BNV::UninterruptedBnvTrajectory solver(system,state.layout,p->OrdinaryContext(),config,
      {p->RunCardIdentity(),"bounded fixture","full/cheap frozen currentness",schedule.Identity()});
    return solver.Integrate(a,b,y,schedule);
}
void window(const std::string& name,Ctx p,double a,State y,double b,double h,const std::filesystem::path& out,int n)
{
    std::cout<<"WINDOW "<<name<<' '<<a<<' '<<b<<'\n';
    std::vector<double> times;for(int i=0;i<=n;++i)times.push_back(a+(b-a)*double(i)/n);
    std::vector<BNV::MainTrajectoryResult> mains;
    for(const auto& tol:tiers)mains.push_back(main_run(p,a,y,b,h,tol,times));
    auto invariant=main_run(p,a,y,b,h,tiers[0],{a,b});
    require(invariant.final_state==mains[0].final_state&&invariant.accepted_steps.size()==mains[0].accepted_steps.size(),"schedule changed main");
    for(std::size_t i=0;i<invariant.accepted_steps.size();++i){auto& x=invariant.accepted_steps[i];auto& z=mains[0].accepted_steps[i];
      require(x.y_right==z.y_right&&x.t_right_s==z.t_right_s&&x.suggested_next_h_s==z.suggested_next_h_s&&x.cumulative_rejected==z.cumulative_rejected,"schedule changed history");}
    const double cap=(b-a)/256;
    auto ref=integrate(p,a,y,b,h,tight,gsl_odeiv2_step_rk8pd,cap);
    auto fine=integrate(p,a,y,b,h,tighter,gsl_odeiv2_step_rk8pd,cap/2);
    auto independent=integrate(p,a,y,b,h,tighter,gsl_odeiv2_step_rkf45,cap/2);
    std::ostringstream rows,stats,events,history,uncertainty,envelopes;
    envelopes<<std::setprecision(17)<<"tier\tobservation\tquantity\tbound_ratio\treference_uncertainty\n";std::size_t envelope_rows=0;
    history<<std::setprecision(17)<<"tier\tordinal\tt_left\tt_right\th_next\tcell\tboundary_knot\tx\te\tm\trejected\n";
    uncertainty<<std::setprecision(17)<<"tier\tobservation\tcomponent\tU_O\tU_event\tF_i\td_O\tD_O1\n";
    std::size_t history_rows=0;rows<<std::setprecision(17)<<"time\ttier\tx\te\tm\taccuracy_x\taccuracy_e\taccuracy_m\treconstruction_qualified\tR18\n";
    stats<<std::setprecision(17)<<"tier\taccepted\trejected\trhs\tevents\tevent_trials\tdiscarded\tmin_step\tmax_step\twall_s\tevent_rejected\tboundary_limited_fraction\n";
    events<<std::setprecision(17)<<"tier\tstep\ttime\tknot\tdirection\tresidual\ttime_low\ttime_high\ttime_ulp\n";
    std::vector<std::vector<State>> observed(3);std::size_t event_count=0;double max_accuracy=0,max_reference_uncertainty=0;std::size_t failed=0;
    auto count=std::make_shared<std::size_t>(0);
    BNV::CheckpointContextFactory factory=[p,count]{return std::make_unique<ActualCheckpointContext>(p,++*count);};
    BNV::Rk8pdCheckpointReconstructor reconstruction;
    std::vector<std::vector<BNV::CheckpointOutput>> all_outputs(3);
    for(int tier=0;tier<3;++tier)
    {
      const auto& m=mains[tier];auto& st=m.statistics;
      for(const auto& r:m.accepted_steps){history<<tier<<'\t'<<r.ordinal<<'\t'<<r.t_left_s<<'\t'<<r.t_right_s<<'\t'<<r.suggested_next_h_s<<'\t'<<r.cstar_cell_left<<'\t'<<r.event.knot<<'\t'<<values(r.y_right)<<'\t'<<r.cumulative_rejected<<'\n';++history_rows;}std::size_t event_n=0;
      for(const auto& step:m.accepted_steps)if(step.event.knot>=0){++event_n;++event_count;const auto& e=step.event;
        events<<tier<<'\t'<<step.ordinal<<'\t'<<e.time<<'\t'<<e.knot<<'\t'<<e.direction<<'\t'<<e.residual<<'\t'<<e.time_low<<'\t'<<e.time_high<<'\t'<<e.time_ulp<<'\n';}
      stats<<tier<<'\t'<<st.accepted_steps<<'\t'<<st.rejected_steps<<'\t'<<st.rhs_evaluations<<'\t'<<event_n<<'\t'<<st.breakpoint.event_trials<<'\t'<<st.breakpoint.discarded<<'\t'<<st.minimum_step_s<<'\t'<<st.maximum_step_s<<'\t'<<st.wall_seconds<<'\t'<<st.breakpoint.event_rejected<<'\t'<<double(event_n)/st.accepted_steps<<'\n';
      BNV::PassiveObservationSchedule schedule(times,"bounded-passive");schedule.Replay(y,m.accepted_steps);
      for(const auto& bracket:schedule.Brackets())
      {
        if(bracket.observation_index%16==0)std::cout<<"CHECKPOINT "<<name<<" tier "<<tier<<" index "<<bracket.observation_index<<'\n';
        auto q=reconstruction.Reconstruct(bracket,m.accepted_steps,factory);
        if(q.status==BNV::CheckpointStatus::Qualified)reconstruction.EvaluateDiagnostics(q,factory);
        if(q.source==BNV::CheckpointSource::Rk8pdReconstructed&&q.status==BNV::CheckpointStatus::Qualified){
          const auto& left=m.accepted_steps.at(bracket.right_accepted_ordinal-1);
          const double span=bracket.t_observation_s-left.t_left_s;
          const auto local=integrate(p,left.t_left_s,left.y_left,bracket.t_observation_s,span,tighter,gsl_odeiv2_step_rk8pd,span/8).y;
          const auto other_local=integrate(p,left.t_left_s,left.y_left,bracket.t_observation_s,span,tighter,gsl_odeiv2_step_rkf45,span/16).y;
          for(int i=0;i<3;++i){const double u=2*std::abs(local[i]-other_local[i]);
            const double ratio=(std::abs(q.O2_state[i]-local[i])+u)/q.reconstruction_error_contribution[i];
            envelopes<<tier<<'\t'<<bracket.observation_index<<"\tstate"<<i<<'\t'<<ratio<<'\t'<<u<<'\n';++envelope_rows;}
          ActualCheckpointContext local_context(p,1);
          const auto lv=DiagnosticValues(local_context.EvaluateDiagnostics(bracket.t_observation_s,local));
          const auto iv=DiagnosticValues(local_context.EvaluateDiagnostics(bracket.t_observation_s,other_local));
          const auto qv=DiagnosticValues(q.diagnostics);
          for(const auto& [key,value]:lv){const double u=2*std::abs(value-iv.at(key));
            const double error=std::abs(qv.at(key)-value)+u;
            const double bound=q.diagnostic_U_O.at(key);
            const double ratio=bound>0?error/bound:(error==0?0:INFINITY);
            envelopes<<tier<<'\t'<<bracket.observation_index<<'\t'<<key<<'\t'<<ratio<<'\t'<<u<<'\n';++envelope_rows;}
        }
        const auto state=q.source==BNV::CheckpointSource::MainEndpoint?q.state:q.O2_state;
        const auto reference=observe(p,fine,bracket.t_observation_s),coarse=observe(p,ref,bracket.t_observation_s),other=observe(p,independent,bracket.t_observation_s);
        rows<<bracket.t_observation_s<<'\t'<<tier<<'\t'<<values(state);
        for(int i=0;i<3;++i){const double floor=64*(std::nextafter(std::abs(reference[i]),INFINITY)-std::abs(reference[i]));
          const double uncertainty=2*std::max(std::abs(reference[i]-coarse[i]),std::abs(reference[i]-other[i]))+floor;
          max_reference_uncertainty=std::max(max_reference_uncertainty,uncertainty/(tiers[2].a[i]+tiers[2].r*std::abs(reference[i])));
          double ratio=(std::abs(state[i]-reference[i])+uncertainty)/(tiers[tier].a[i]+tiers[tier].r*std::max(std::abs(state[i]),std::abs(reference[i])));
          rows<<'\t'<<ratio;max_accuracy=std::max(max_accuracy,ratio);}
        for(int i=0;i<3;++i)uncertainty<<tier<<'\t'<<bracket.observation_index<<'\t'<<i<<'\t'<<q.qualification.U_O[i]<<'\t'<<q.qualification.U_event[i]<<'\t'<<q.qualification.F_i[i]<<'\t'<<q.qualification.d_O[i]<<'\t'<<q.qualification.D_O1[i]<<'\n';
        const bool ok=q.status==BNV::CheckpointStatus::Qualified&&q.diagnostics_evaluated;
        rows<<'\t'<<ok<<'\t'<<(ok?q.diagnostics.R18_residual_erg_s:0)<<'\n';failed+=!ok;
        observed[tier].push_back(state);all_outputs[tier].push_back(q);
      }
      // Exact endpoints bypass factory; nextafter-adjacent points must be valid brackets.
      for(const auto& step:m.accepted_steps)if(step.event.knot>=0){
        BNV::ObservationBracket exact{0,step.t_right_s,step.ordinal-1,step.ordinal,step.t_left_s,step.t_right_s,false};
        auto q=reconstruction.Reconstruct(exact,m.accepted_steps,{});require(q.state==step.y_right,"endpoint changed");
        exact.t_observation_s=std::nextafter(step.t_right_s,step.t_left_s);
        auto adjacent=reconstruction.Reconstruct(exact,m.accepted_steps,factory);
        require(adjacent.status==BNV::CheckpointStatus::Qualified,"adjacent boundary reconstruction unresolved");
        if(step.ordinal<m.accepted_steps.size()){
          const auto& next=m.accepted_steps.at(step.ordinal);
          BNV::ObservationBracket right{0,std::nextafter(step.t_right_s,next.t_right_s),step.ordinal,next.ordinal,next.t_left_s,next.t_right_s,false};
          auto right_adjacent=reconstruction.Reconstruct(right,m.accepted_steps,factory);
          require(right_adjacent.status==BNV::CheckpointStatus::Qualified,"right-adjacent boundary reconstruction unresolved");}}

    }
    std::size_t c1=0,c2=0,c3=0;double contraction=0;
    for(int j=0;j<=n;++j)for(int i=0;i<3;++i){double B=observed[0][j][i],R=observed[1][j][i],U=observed[2][j][i];
      double dBR=std::abs(B-R),dRU=std::abs(R-U),M=std::max({std::abs(B),std::abs(R),std::abs(U)});
      double F=std::max(tiers[2].a[i]+tiers[2].r*M,64*(std::nextafter(M,INFINITY)-M));
      c1+=dRU>tiers[1].a[i]+tiers[1].r*M;
      if(dBR>10*F){c2+=dRU/dBR>.1;contraction=std::max(contraction,dRU/dBR);}else c3+=dRU>10*F;}
    std::size_t ledger_c1=0,ledger_c2=0,ledger_c3=0;double max_ledger_accuracy=0,max_r18=0;
    const auto ledger=[](const BNV::BnvDiagnostics& d){return std::array<double,11>{{d.P_dir_eq_erg_s,d.P_dir_actual_erg_s,d.LH_erg_s,d.DeltaLnu_erg_s,d.DeltaPbeta_erg_s,d.Lnu_eq_erg_s,d.Lnu_full_erg_s,d.L_out_fluid_inf_erg_s,d.Lgamma_erg_s,d.Lother_erg_s,d.Pnet_erg_s}};};
    const auto gross=[](const BNV::BnvDiagnostics& d){return std::max(1.,std::abs(d.P_dir_actual_erg_s)+std::abs(d.LH_erg_s)+std::abs(d.DeltaLnu_erg_s)+std::abs(d.Lnu_eq_erg_s)+std::abs(d.Lgamma_erg_s)+std::abs(d.Lother_erg_s)+std::abs(d.L_out_fluid_inf_erg_s));};
    for(int j=0;j<=n;++j){
      std::array<std::array<double,11>,3> val;double G=1;
      for(int tier=0;tier<3;++tier){const auto& q=all_outputs[tier][j];val[tier]=ledger(q.diagnostics);G=std::max(G,gross(q.diagnostics));
        if(q.diagnostics_evaluated)max_r18=std::max(max_r18,std::abs(q.diagnostics.R18_residual_erg_s)/gross(q.diagnostics));}
      for(int k=0;k<11;++k){double dBR=std::abs(val[0][k]-val[1][k]),dRU=std::abs(val[1][k]-val[2][k]);
        double M=std::max({std::abs(val[0][k]),std::abs(val[1][k]),std::abs(val[2][k])});
        double F=std::max(1e-11*G,64*(std::nextafter(M,INFINITY)-M));ledger_c1+=dRU/G>1e-9;
        if(dBR>10*F)ledger_c2+=dRU/dBR>.1;else ledger_c3+=dRU>F;}
      ActualCheckpointContext dc(p,1);
      auto fd=dc.EvaluateDiagnostics(times[j],observe(p,fine,times[j]));
      auto cd=dc.EvaluateDiagnostics(times[j],observe(p,ref,times[j]));
      auto id=dc.EvaluateDiagnostics(times[j],observe(p,independent,times[j]));
      const auto lv=ledger(fd),lc=ledger(cd),li=ledger(id);
      for(int tier=1;tier<3;++tier)for(int k=0;k<11;++k){
        const double u=2*std::max(std::abs(lv[k]-lc[k]),std::abs(lv[k]-li[k]));
        max_ledger_accuracy=std::max(max_ledger_accuracy,(std::abs(val[tier][k]-lv[k])+u)/(1e-9*std::max(G,gross(fd))));}
    }
    std::ostringstream ledger_out;ledger_out<<std::setprecision(17)<<"legacy_c1\tlegacy_c2\tlegacy_c3\tmax_accuracy_ratio\tmax_R18_normalized\n"<<ledger_c1<<'\t'<<ledger_c2<<'\t'<<ledger_c3<<'\t'<<max_ledger_accuracy<<'\t'<<max_r18<<'\n';
    PublishTsv(out/(name+"-ledger.tsv"),ledger_out.str(),5,1);
    std::ostringstream r20;r20<<std::setprecision(17)<<"tier\tR20_normalized\treconstruction_over_N\tquadrature_grid_points\n";
    if(failed==0){for(int tier=0;tier<3;++tier){auto result=BNV::ComputeCheckpointR20(all_outputs[tier]);
      r20<<tier<<'\t'<<result.R20_normalized<<'\t'<<result.reconstruction_uncertainty_erg/result.N_R20_erg<<'\t'<<all_outputs[tier].size()<<'\n';}
      PublishTsv(out/(name+"-r20.tsv"),r20.str(),4,3);}
    std::ostringstream summary;summary<<std::setprecision(17)<<"max_accuracy_ratio\tunresolved_checkpoints\tlegacy_c1\tlegacy_c2\tlegacy_c3\tmax_contraction\tschedule_invariant\tmax_reference_uncertainty_over_DU\n"<<max_accuracy<<'\t'<<failed<<'\t'<<c1<<'\t'<<c2<<'\t'<<c3<<'\t'<<contraction<<"\t1\t"<<max_reference_uncertainty<<'\n';
    PublishTsv(out/(name+"-envelopes.tsv"),envelopes.str(),5,envelope_rows);
    PublishTsv(out/(name+"-history.tsv"),history.str(),11,history_rows);
    PublishTsv(out/(name+"-uncertainty.tsv"),uncertainty.str(),8,9*(n+1));
    PublishTsv(out/(name+"-observations.tsv"),rows.str(),10,3*(n+1));PublishTsv(out/(name+"-statistics.tsv"),stats.str(),12,3);
    PublishTsv(out/(name+"-events.tsv"),events.str(),9,event_count);PublishTsv(out/(name+"-summary.tsv"),summary.str(),8,1);
    std::cout<<"RESULT "<<name<<" accuracy "<<max_accuracy<<" unresolved "<<failed<<" old "<<c1<<' '<<c2<<' '<<c3<<'\n';
}
int main(int argc,char** argv)
{
 try{
    gsl_set_error_handler_off();std::cout<<std::setprecision(17)<<std::unitbuf;
    require(argc==11,"profile certificate thermal work entry frozen coefficients archived-inputs output mode");
    Inputs in;in.profile=argv[1];in.certificate=argv[2];in.thermal=argv[3];in.work=argv[4];in.entry=argv[5];in.frozen=argv[6];in.coefficients=argv[7];in.output=argv[9];
    require(!std::filesystem::exists(in.work)&&!std::filesystem::exists(in.output),"fresh qualification output required");std::filesystem::create_directories(in.output);
    auto p=BuildBoth(in);
    {const auto dc=p.source->OrdinaryContext()->DriverContext();const auto cells=dc.star->ReadHeatCapacityCells(*dc.thermo,dc.geo);
      for(std::size_t i=1;i+1<cells.temperature_MeV.size();++i)for(int side:{-1,1}){
        const std::size_t hint=side<0?i-1:i;
        for(double T:{std::nextafter(cells.temperature_MeV[i],0.),cells.temperature_MeV[i],std::nextafter(cells.temperature_MeV[i],INFINITY)}){
          dc.star->HeatCapacityStar_Tinf(std::sqrt(cells.temperature_MeV[hint]*cells.temperature_MeV[hint+1]),*dc.thermo,dc.geo);
          require(dc.star->HeatCapacityStar_Tinf(T,*dc.thermo,dc.geo)==cells.SmoothCapacity(T,cells.Cell(T,hint)),"typed cache interface differs from production");}}
      std::cout<<"CACHE exact interior/knots/adjacent values and both hint directions PASS\n";}
auto source=ReadAccepted(std::filesystem::path(argv[8])/"source-accepted.tsv");auto control=ReadAccepted(std::filesystem::path(argv[8])/"control-accepted.tsv");
    const std::string mode=argv[10];
    if(mode=="all"||mode=="events"){
      auto historical=main_run(p.source,0,{{0,0,0}},462269531250.,1,tiers[2],{0,462269531250.},true);
      require(historical.final_state==State{{0.49240008824076903,-2.5123474256442210e-7,-4.7906773046561003e-7}}
        &&historical.statistics.accepted_steps==232&&historical.statistics.rejected_steps==60,"old Arm E failed exact reproduction");
      std::ifstream cases("docs/validation/evidence/obs32-forensics/cases.tsv");std::string line;std::getline(cases,line);
      std::ostringstream reproduction;reproduction<<std::setprecision(17)<<"observation\told_d_over_D\told_U_over_budget\tnew_qualified\tsmooth_bitwise\n";
      std::size_t count=0,event_rows=0;BNV::Rk8pdCheckpointReconstructor old(false),candidate;std::size_t id=0;
      std::ostringstream event_evidence;event_evidence<<std::setprecision(17)<<"observation\tlevel\tknot\ttime\tresidual\ttime_ulp\tx\te\tm\n";
      BNV::CheckpointContextFactory factory=[&]{return std::make_unique<ActualCheckpointContext>(p.source,++id);};
      while(std::getline(cases,line)){auto v=Split(line);require(v.size()==13,"forensic case schema");
        BNV::AcceptedStepRecord step;step.ordinal=1;step.t_left_s=std::stod(v[1]);step.t_right_s=std::stod(v[3]);
        step.y_left={{std::stod(v[4]),std::stod(v[5]),std::stod(v[6])}};step.y_right={{std::stod(v[7]),std::stod(v[8]),std::stod(v[9])}};
        BNV::ObservationBracket bracket{0,std::stod(v[2]),0,1,step.t_left_s,step.t_right_s,false};
        auto a=old.Reconstruct(bracket,{step},factory),b=candidate.Reconstruct(bracket,{step},factory);
        const std::vector<std::pair<std::string,std::vector<Stepper::Event>>> levels={
          {"O1",b.provenance.O1_events},{"O2",b.provenance.O2_events},
          {"half-rk8pd",b.provenance.refinement_events},{"quarter-rkf45",b.provenance.independent_events}};
        for(const auto& level:levels)for(const auto& ev:level.second){event_evidence<<v[0]<<'\t'<<level.first<<'\t'<<ev.knot<<'\t'<<ev.time<<'\t'<<ev.residual<<'\t'<<ev.time_ulp<<'\t'<<values(ev.state)<<'\n';++event_rows;}
        if(std::stoi(v[0])==32){require(a.O1_state==State{{0.12472490938574241,-2.0287605965623057e-7,-5.0871537810150439e-7}}
          &&a.O2_state==State{{0.12472490938614615,-2.0287605965622805e-7,-5.0871537810149634e-7}}
          &&a.status==BNV::CheckpointStatus::NumericallyUnresolved,"obs32 failure reproduction changed");}
        else require(a.O1_state==b.O1_state&&a.O2_state==b.O2_state,"smooth interval changed");
        reproduction<<v[0]<<'\t'<<a.qualification.d_O[0]/a.qualification.D_O1[0]<<'\t'<<a.qualification.U_O[0]/(.2*a.qualification.F_i[0])<<'\t'<<(b.status==BNV::CheckpointStatus::Qualified)<<'\t'<<(a.O1_state==b.O1_state&&a.O2_state==b.O2_state)<<'\n';++count;
      }
      PublishTsv(in.output/"historical-reproduction.tsv",reproduction.str(),5,count);
      PublishTsv(in.output/"reconstruction-events.tsv",event_evidence.str(),9,event_rows);
      std::cout<<"HISTORICAL Arm E exact 232/60; obs32 old states exact; smooth intervals exact\n";
    }

    if(mode=="source"||mode=="all")window("source",p.source,source[19].t1,source[19].y,1e12,source[19].hnext,in.output,64);
    if(mode=="control"||mode=="all")window("control",p.control,control[48].t1,control[48].y,9e12,control[48].hnext,in.output,64);
    if(mode=="bounded"||mode=="all")window("bounded",p.source,0,{{0,0,0}},462269531250.,1,in.output,240);
    if(mode=="tail"||mode=="all"){
      const auto i=std::find_if(source.begin(),source.end(),[](const Archived& a){return a.t1>1.3e13;});require(i!=source.end(),"missing tail state");
      window("tail",p.source,i->t1,i->y,std::min(i->t1+5e11,source.back().t1),i->hnext,in.output,32);}
    if(mode=="crossings"||mode=="all"||mode=="events"){
      std::ifstream input(std::filesystem::path(argv[8])/"steps.tsv");std::string line;std::getline(input,line);std::ostringstream out;
      out<<std::setprecision(17)<<"run\tstep\ttier\tcomponent\taccuracy_ratio\tevents\tevent_accuracy_ratio\n";std::size_t count=0,crossings=0;double max=0;
      while(std::getline(input,line)){auto v=Split(line);require(v.size()==12,"archived crossing schema");if(v[11]!="crossing")continue;++crossings;
        auto ctx=v[0]=="source"?p.source:p.control;double a=std::stod(v[2]),b=std::stod(v[7]),h=std::stod(v[6]);State y{{std::stod(v[3]),std::stod(v[4]),std::stod(v[5])}};
        auto ref=integrate(ctx,a,y,b,h,tight,gsl_odeiv2_step_rk8pd,(b-a)/16);
        auto fine=integrate(ctx,a,y,b,h,tighter,gsl_odeiv2_step_rk8pd,(b-a)/32);
        auto other=integrate(ctx,a,y,b,h,tighter,gsl_odeiv2_step_rkf45,(b-a)/32);
        for(int tier=0;tier<3;++tier){auto candidate=integrate(ctx,a,y,b,h,tiers[tier],gsl_odeiv2_step_rkf45);
          for(int i=0;i<3;++i){double uncertainty=2*std::max(std::abs(ref.y[i]-fine.y[i]),std::abs(other.y[i]-fine.y[i]));
            double ratio=(std::abs(candidate.y[i]-fine.y[i])+uncertainty)/(tiers[tier].a[i]+tiers[tier].r*std::max(std::abs(candidate.y[i]),std::abs(fine.y[i])));
            double event_ratio=0;
            for(const auto& segment:candidate.steps)if(segment.event.knot>=0){
              const auto er=observe(ctx,fine,segment.right),ec=observe(ctx,ref,segment.right),ei=observe(ctx,other,segment.right);
              const double eu=2*std::max(std::abs(er[i]-ec[i]),std::abs(er[i]-ei[i]));
              event_ratio=std::max(event_ratio,(std::abs(segment.after[i]-er[i])+eu)/(tiers[tier].a[i]+tiers[tier].r*std::max(std::abs(segment.after[i]),std::abs(er[i]))));}
            max=std::max(max,ratio);out<<v[0]<<'\t'<<v[1]<<'\t'<<tier<<'\t'<<i<<'\t'<<ratio<<'\t'<<candidate.stats.discarded<<'\t'<<event_ratio<<'\n';++count;}}
      }
      require(crossings==37,"archived crossing count");PublishTsv(in.output/"crossings.tsv",out.str(),7,count);std::cout<<"CROSSINGS "<<crossings<<" max_accuracy "<<max<<'\n';
    }
    return 0;
 }catch(const std::exception&e){std::cerr<<"UNRESOLVED "<<e.what()<<'\n';return 1;}
}
