// Uses the accepted actual-fixture wrappers and evidence serializers.
// The historical bounded qualification main is not executed.
#define PHASE6A1_CAMPAIGN_SUPPORT_ONLY
#include "adr0017_production_qualification.cpp"

namespace {
std::shared_ptr<const BNV::FrozenControlledBnvRunContext> BuildCampaignContext(const Inputs& inputs,bool control)
{
    const auto& card=Campaign::RunCards().back();
    require(card.identity=="CPL-P2-LINEAR-QSS-v1"&&card.partition=="P2"
      &&card.fractional_drive_per_year==-1.0e-12&&card.duration_year==5.0e5
      &&card.checkpoints==8193&&!card.reaction_free,"bounded P2 card changed");
    auto fixture=Fixture(inputs.profile,inputs.certificate.string(),inputs.work/"owning-star",80000);
    auto channels=Channels(fixture);
    EOS::CompOSE_Thermo::Options options;options.Tmin_for_derivative_MeV=0;options.clamp_to_domain=false;
    auto thermal=std::make_shared<const RC::FrozenThermalSource>(inputs.thermal,options,
      "controlled mathematical fixed-background free-gas entropy; qualified radial80000");
    auto tangent=Campaign::Tangent(fixture,inputs.profile,inputs.work);
    auto bnv_token=std::make_shared<BNV::BnvDependencyToken>();
    auto monitor=Phase6A1Test::LoadFrozenMonitor(inputs.frozen,tangent,bnv_token);
    auto run_token=std::make_shared<RC::RunDependencyToken>();
    auto spin=std::make_shared<const BNV::StaticZeroSpinHistory>(run_token);
    auto ordinary=Context(fixture,channels,thermal,spin,
      Campaign::Qualification(inputs.profile,inputs.certificate,inputs.entry),run_token);
    const double mu_B=Campaign::Column(inputs.coefficients,"mu_B_inf");
    const double Bdot=card.fractional_drive_per_year*tangent->B0Count()/Year;
    require(Bdot==-2.4136520263641375e37,"bounded P2 Bdot changed");
    const std::string fate_id="phase6a1-P2-terminal-fate-v1";
    BNV::ProductFateLedger fate(fate_id,{{"controlled-terminal-product",BNV::TerminalProductFate::BoundInert,1,
      "abstract-neutron-disappearance"}});
    auto partition=Partition(card,fixture,bnv_token,fate_id);
    auto history=std::make_shared<const Phase6A1Test::UniformProperNeutronSinkHistory>(fixture.central,
      tangent->B0Count(),Bdot,tangent->DomainIdentity(),tangent->StarIdentity(),tangent->SequenceStateIdentity(),
      partition->Identity(),fate_id,bnv_token);
    if(control)
    {
        auto zero_history=std::make_shared<const Campaign::ExactZeroHistory>(tangent->B0Count(),
          tangent->DomainIdentity(),tangent->StarIdentity(),tangent->SequenceStateIdentity(),bnv_token);
        auto zero_partition=std::make_shared<const Campaign::ExactZeroPartition>(bnv_token);
        BNV::ProductFateLedger zero_fate("zero-control-fate",{{"none",BNV::TerminalProductFate::BoundInert,1,"zero-control"}});
        return std::make_shared<const BNV::FrozenControlledBnvRunContext>(ordinary,tangent,zero_history,
          zero_partition,zero_fate,monitor,channels,mu_B,
          "frozen Structure-1 B0 equilibrium mu_B plus governed moving-reference actual-potential correction",
          card.identity+"-MATCHED-CONTROL");
    }
    return std::make_shared<const BNV::FrozenControlledBnvRunContext>(ordinary,tangent,history,partition,
      fate,monitor,channels,mu_B,
      "frozen Structure-1 B0 equilibrium mu_B plus governed moving-reference actual-potential correction",
      card.identity);
}

void WriteCampaignSummary(const Inputs& inputs,const BNV::MainTrajectoryResult& main,
  const std::vector<BNV::CheckpointOutput>& outputs,double assembly_wall,double assembly_cpu,
  double output_wall,std::size_t contexts,const std::string& tier,const std::string& mode,bool pass)
{
    std::ofstream out(inputs.output/"performance.tsv");out<<std::setprecision(17)<<"key\tvalue\n"
      <<"mode\t"<<mode<<"\ntier\t"<<tier<<"\npass\t"<<pass
      <<"\nphysics_context_wall_s\t"<<assembly_wall<<"\nphysics_context_cpu_s\t"<<assembly_cpu
      <<"\nmain_wall_s\t"<<main.statistics.wall_seconds<<"\nmain_cpu_s\t"<<main.statistics.cpu_seconds
      <<"\ncheckpoint_output_wall_s\t"<<output_wall<<"\ncontext_count\t"<<contexts
      <<"\nmain_accepted\t"<<main.statistics.accepted_steps<<"\nmain_rejected\t"<<main.statistics.rejected_steps
      <<"\nmain_rhs\t"<<main.statistics.rhs_evaluations<<"\nmain_run_count\t1"
      <<"\nunique_positive_t1_targets\t"<<main.statistics.distinct_positive_t1_targets
      <<"\npositive_t1_s\t"<<main.statistics.positive_t1_s
      <<"\nrequested_checkpoints\t8193\nprocessed_checkpoints\t"<<outputs.size();
    std::size_t interior=0,qualified=0,invocations=0;
    double o1=0,o2=0,max_d=0,max_u=0;
    for(const auto& q:outputs)
    {
        interior+=q.source==BNV::CheckpointSource::Rk8pdReconstructed;
        qualified+=q.status==BNV::CheckpointStatus::Qualified&&q.diagnostics_evaluated;
        invocations+=q.provenance.rk8pd_invocations;
        o1+=q.performance.O1_solve_wall_seconds;o2+=q.performance.O2_solve_wall_seconds;
        if(q.source==BNV::CheckpointSource::Rk8pdReconstructed)
            for(std::size_t i=0;i<3;++i){max_d=std::max(max_d,q.qualification.d_O[i]/q.qualification.D_O1[i]);
                max_u=std::max(max_u,q.qualification.U_O[i]/(.20*q.qualification.F_i[i]));}
    }
    out<<"\nstrict_interior\t"<<interior<<"\nqualified_checkpoints\t"<<qualified
      <<"\nrk8pd_invocations\t"<<invocations<<"\nO1_solve_wall_s\t"<<o1<<"\nO2_solve_wall_s\t"<<o2
      <<"\nmax_d_over_D_O1\t"<<max_d<<"\nmax_U_over_point2F\t"<<max_u<<'\n';
    if(pass)
    {
        const auto r=BNV::ComputeCheckpointR20(outputs);
        out<<"R20_residual_erg\t"<<r.R20_residual_erg<<"\nN_R20_erg\t"<<r.N_R20_erg
          <<"\nR20_normalized\t"<<r.R20_normalized<<"\nR20_reconstruction_uncertainty_erg\t"<<r.reconstruction_uncertainty_erg<<'\n';
    }
}
}

int main(int argc,char** argv)
{
    try
    {
        gsl_set_error_handler_off();std::cout<<std::setprecision(17)<<std::unitbuf;
        require(argc==12,"profile certificate thermal work entry frozen coefficients output mode tier prerequisite-flag");
        Inputs inputs;inputs.profile=argv[1];inputs.certificate=argv[2];inputs.thermal=argv[3];
        inputs.work=argv[4];inputs.entry=argv[5];inputs.frozen=argv[6];inputs.coefficients=argv[7];inputs.output=argv[8];
        const std::string mode=argv[9],tier=argv[10];
        require(mode=="source"||mode=="control","invalid campaign mode");
        BNV::Rkf45Configuration config;
        if(tier=="baseline"){config.relative=1e-7;config.absolute={{1e-12,1e-18,1e-18}};}
        else if(tier=="refined"){config.relative=1e-9;config.absolute={{1e-14,1e-20,1e-20}};}
        else require(tier=="ultra","invalid campaign tier");
        {std::ifstream flag(argv[11]);std::string value((std::istreambuf_iterator<char>(flag)),{});
         require(value=="BNV_RESUME_PREREQUISITES_PASS\n","prerequisites not authenticated");}
        require(!std::filesystem::exists(inputs.work)&&!std::filesystem::exists(inputs.output),"fresh work/output required");
        std::filesystem::create_directories(inputs.output);
        const auto begin=Clock::now();const auto cpu=std::clock();
        auto physics=BuildCampaignContext(inputs,mode=="control");
        const double assembly_wall=std::chrono::duration<double>(Clock::now()-begin).count();
        const double assembly_cpu=double(std::clock()-cpu)/CLOCKS_PER_SEC;
        RunState state(physics->OrdinaryContext());state.thermal.SetTinf(1e8);RC::ChemicalImbalanceState(0,0).Store(state.chem);
        auto driver=std::make_shared<BNV::ControlledBnvSecularDriver>(physics,BNV::ControlledBnvSecularDriver::Mode::Coupled);
        EV::EvolutionSystem system(state.ctx,state.state,state.rhs,state.layout,{driver});
        const auto& card=Campaign::RunCards().back();const double end=card.duration_year*Year;
        std::vector<double> times;for(std::size_t i=0;i<card.checkpoints;++i)
            times.push_back(end*double(i)/double(card.checkpoints-1));
        BNV::PassiveObservationSchedule schedule(times,"CPL-P2-8193-uniform-passive-v1");
        BNV::MainIntegrationProvenance provenance{physics->RunCardIdentity(),mode,
          "full-boundary-cheap-per-RHS-frozen-validity",schedule.Identity()};
        BNV::UninterruptedBnvTrajectory solver(system,state.layout,physics->OrdinaryContext(),config,provenance);
        // No observers evaluate diagnostics on the authoritative evolving context.
        const auto main=solver.Integrate(0,end,{{0,0,0}},schedule,true);
        require(main.statistics.distinct_positive_t1_targets==1,"multiple main ceilings");
        WriteMainEvidence(inputs.output,main,schedule,0);
        std::cout<<"MAIN_COMPLETE "<<mode<<' '<<tier<<" accepted "<<main.statistics.accepted_steps
          <<" rejected "<<main.statistics.rejected_steps<<'\n';
        auto count=std::make_shared<std::size_t>(0);
        BNV::CheckpointContextFactory factory=[physics,count]{return std::make_unique<ActualCheckpointContext>(physics,++*count);};
        BNV::Rk8pdCheckpointReconstructor reconstructor;
        Capture capture(inputs.output/"trajectory.tsv",*driver);
        std::vector<BNV::CheckpointOutput> outputs;std::vector<MatrixRow> rows;
        const auto output_begin=Clock::now();bool pass=true;
        for(const auto& bracket:schedule.Brackets())
        {
            if(!bracket.initial){MatrixRow row;row.knots=main.accepted_steps.at(bracket.right_accepted_ordinal-1).cstar_knots_crossed;
              row.category=row.knots==0?'A':row.knots==1?'B':'C';rows.push_back(row);}
            auto q=reconstructor.Reconstruct(bracket,main.accepted_steps,factory);
            if(q.status==BNV::CheckpointStatus::Qualified)reconstructor.EvaluateDiagnostics(q,factory);
            pass=q.status==BNV::CheckpointStatus::Qualified&&q.diagnostics_evaluated;
            outputs.push_back(std::move(q));
            if(!pass){std::cout<<"NUMERICALLY_UNRESOLVED index "<<bracket.observation_index<<" time "<<bracket.t_observation_s<<'\n';break;}
            capture.SaveSnapshot(outputs.back().state[0],outputs.back().diagnostics);
            if(bracket.observation_index%1024==0)std::cout<<"CHECKPOINT "<<bracket.observation_index<<'\n';
        }
        WriteCheckpointEvidence(inputs.output/"checkpoints.tsv",outputs,rows);
        WriteCampaignSummary(inputs,main,outputs,assembly_wall,assembly_cpu,
          std::chrono::duration<double>(Clock::now()-output_begin).count(),*count,tier,mode,pass);
        require(pass&&outputs.size()==card.checkpoints,"ADR-0017 reconstruction gate failed; no retry or fallback");
        std::cout<<"CAMPAIGN_TRAJECTORY_COMPLETE "<<mode<<' '<<tier<<'\n';return 0;
    }
    catch(const std::exception& e){std::cerr<<"STOP "<<e.what()<<'\n';return 1;}
}
