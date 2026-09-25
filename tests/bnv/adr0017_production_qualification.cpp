#define PHASE6A1_BA12R_ULTRA_SUPPORT_ONLY
#include "coupled_trajectory.cpp"

#include <CompactStar/Physics/BNV/PassiveCheckpointOutput.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace
{
using Clock=std::chrono::steady_clock;

struct MatrixRow
{
    std::size_t observation=0,left_index=0,right_index=0;
    double t_obs=0,t_left=0,t_right=0;
    BNV::CheckpointState left{},right{};
    char category='A';int knots=0;bool deep=false,exact=false,strict=false;
};

struct Inputs
{
    std::filesystem::path matrix,profile,certificate,thermal,work,entry,frozen,coefficients,output,oracle_flag;
};

std::vector<std::string> Split(const std::string& line)
{
    std::vector<std::string> values;std::string item;std::istringstream input(line);
    while(std::getline(input,item,'\t'))values.push_back(item);return values;
}

std::vector<MatrixRow> ReadMatrix(const std::filesystem::path& path)
{
    std::ifstream input(path);require(bool(input),"cannot open solve matrix");
    std::string line;require(bool(std::getline(input,line)),"empty solve matrix");
    const auto header=Split(line);require(header.size()==26&&header.front()=="observation_index"
      &&header.back()=="local_integrations","solve matrix schema changed");
    std::vector<MatrixRow> rows;std::size_t A=0,B=0,C=0,deep=0,exact=0,strict=0;
    while(std::getline(input,line))
    {
        if(line.empty())continue;const auto v=Split(line);require(v.size()==26,"solve matrix row shape changed");
        MatrixRow row;row.observation=std::stoull(v[0]);row.t_obs=std::stod(v[1]);row.left_index=std::stoull(v[2]);
        row.t_left=std::stod(v[3]);row.left={{std::stod(v[4]),std::stod(v[5]),std::stod(v[6])}};
        row.right_index=std::stoull(v[7]);row.t_right=std::stod(v[8]);
        row.right={{std::stod(v[9]),std::stod(v[10]),std::stod(v[11])}};
        require(v[12].size()==1,"invalid knot category");row.category=v[12][0];row.knots=std::stoi(v[13]);
        row.deep=std::stoi(v[14])!=0;row.exact=std::stoi(v[15])!=0;row.strict=std::stoi(v[16])!=0;
        require(row.observation==rows.size()+1&&row.right_index==row.left_index+1,"matrix observation identity changed");
        require(row.strict==(row.t_left<row.t_obs&&row.t_obs<row.t_right),"matrix interior semantics changed");
        require(row.exact==(row.t_obs==row.t_right),"matrix endpoint semantics changed");
        require(std::stoull(v[25])==(row.strict?6u:0u),"matrix solve authority changed");
        A+=row.category=='A';B+=row.category=='B';C+=row.category=='C';deep+=row.deep;exact+=row.exact;strict+=row.strict;
        rows.push_back(row);
    }
    require(rows.size()==240&&A==237&&B==3&&C==0&&deep==81&&exact==1&&strict==239,
            "matrix population changed");
    return rows;
}

std::shared_ptr<const BNV::FrozenControlledBnvRunContext> BuildContext(const Inputs& inputs)
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
    return std::make_shared<const BNV::FrozenControlledBnvRunContext>(ordinary,tangent,history,partition,
      fate,monitor,channels,mu_B,
      "frozen Structure-1 B0 equilibrium mu_B plus governed moving-reference actual-potential correction",
      card.identity);
}

class ActualCheckpointContext final:public BNV::ICheckpointEvaluationContext
{
  public:
    ActualCheckpointContext(std::shared_ptr<const BNV::FrozenControlledBnvRunContext> context,std::size_t identity)
      :context_(std::move(context)),state_(context_->OrdinaryContext()),
       driver_(std::make_shared<BNV::ControlledBnvSecularDriver>(context_,BNV::ControlledBnvSecularDriver::Mode::Coupled)),
       system_(state_.ctx,state_.state,state_.rhs,state_.layout,{driver_}),identity_(identity){}
    void RequireFullCurrent() const override{context_->RequireFullCurrent();}
    void Derivative(double t,const BNV::CheckpointState& state,BNV::CheckpointState& output) override
    {
        if(system_(t,state.data(),output.data())!=0)throw std::runtime_error("production checkpoint RHS failure");
        for(double value:output)require(std::isfinite(value),"nonfinite production checkpoint RHS");
    }
    BNV::BnvDiagnostics EvaluateDiagnostics(double t,const BNV::CheckpointState& packed) override
    {
        EV::UnpackStateVector(state_.state,state_.layout,packed.data());
        return driver_->Evaluate(t,state_.state,state_.ctx).diagnostics;
    }
    std::string Identity() const override
    {return "CPL-P2-LINEAR-QSS-v1:isolated-production-context:"+std::to_string(identity_);}
  private:
    std::shared_ptr<const BNV::FrozenControlledBnvRunContext> context_;
    RunState state_;
    std::shared_ptr<BNV::ControlledBnvSecularDriver> driver_;
    EV::EvolutionSystem system_;
    std::size_t identity_;
};

void ValidateAgainstMatrix(const BNV::MainTrajectoryResult& main,
  const BNV::PassiveObservationSchedule& schedule,const std::vector<MatrixRow>& rows)
{
    require(main.final_state==BNV::CheckpointState{{0.49240008824076903,-2.5123474256442210e-7,
      -4.7906773046561003e-7}},"Arm-E final state mismatch");
    require(main.statistics.accepted_steps==232&&main.statistics.rejected_steps==60,
            "Arm-E accepted/rejected mismatch");
    require(main.statistics.distinct_positive_t1_targets==1
      &&main.statistics.positive_t1_s==462269531250.0,"main target authority changed");
    require(main.accepted_steps.size()==232&&schedule.Brackets().size()==241,"main evidence count changed");
    for(const auto& row:rows)
    {
        const auto& step=main.accepted_steps.at(row.right_index-1);
        const auto& bracket=schedule.Brackets().at(row.observation);
        require(step.ordinal==row.right_index&&step.t_left_s==row.t_left&&step.t_right_s==row.t_right
          &&step.y_left==row.left&&step.y_right==row.right,"production main history differs from matrix");
        require(bracket.observation_index==row.observation&&bracket.t_observation_s==row.t_obs
          &&bracket.left_accepted_ordinal==row.left_index&&bracket.right_accepted_ordinal==row.right_index
          &&bracket.t_left_s==row.t_left&&bracket.t_right_s==row.t_right,"production passive bracket differs");
        require(step.cstar_knots_crossed==row.knots,"production knot classification differs");
    }
}

void WriteMainEvidence(const std::filesystem::path& root,const BNV::MainTrajectoryResult& main,
  const BNV::PassiveObservationSchedule& schedule,std::size_t trajectory_rows)
{
    std::ofstream schedule_out(root/"schedule.tsv");schedule_out<<std::setprecision(17)<<"index\tt_s\n";
    for(std::size_t i=0;i<schedule.RequestedTimes().size();++i)
        schedule_out<<i<<'\t'<<schedule.RequestedTimes()[i]<<'\n';
    std::ofstream accepted(root/"main.accepted_states.tsv");accepted<<std::setprecision(17)
      <<"sequence\tt_before_s\tt_after_s\tstep_s\tsuggested_next_h_s\tcumulative_rejected\tx_state\teta_e_MeV\teta_mu_MeV\n";
    for(const auto& step:main.accepted_steps)accepted<<step.ordinal<<'\t'<<step.t_left_s<<'\t'<<step.t_right_s<<'\t'
      <<step.accepted_step_s<<'\t'<<step.suggested_next_h_s<<'\t'<<step.cumulative_rejected<<'\t'
      <<step.y_right[0]<<'\t'<<step.y_right[1]<<'\t'<<step.y_right[2]<<'\n';
    std::ofstream internal(root/"main.internal_steps.tsv");internal<<std::setprecision(17)
      <<"sequence\tarchived_checkpoint_index\tceiling_s\tt_before_s\tt_after_s\tstep_s\tsuggested_next_h_s"
      <<"\tcumulative_rejected\tx_before\tx_after\tcstar_cell_before\tcstar_cell_after\tcstar_knots_crossed\n";
    for(const auto& step:main.accepted_steps)internal<<step.ordinal<<"\t240\t"<<main.final_time_s<<'\t'
      <<step.t_left_s<<'\t'<<step.t_right_s<<'\t'<<step.accepted_step_s<<'\t'<<step.suggested_next_h_s<<'\t'
      <<step.cumulative_rejected<<'\t'<<step.y_left[0]<<'\t'<<step.y_right[0]<<'\t'
      <<step.cstar_cell_left<<'\t'<<step.cstar_cell_right<<'\t'<<step.cstar_knots_crossed<<'\n';
    std::ofstream observations(root/"main.observations.tsv");observations<<std::setprecision(17)
      <<"observation_index\trequested_t_s\tprevious_accepted_step\tnew_accepted_step\tt_previous_s\tt_new_s\tkind\n";
    for(const auto& bracket:schedule.Brackets())observations<<bracket.observation_index<<'\t'
      <<bracket.t_observation_s<<'\t'<<bracket.left_accepted_ordinal<<'\t'<<bracket.right_accepted_ordinal<<'\t'
      <<bracket.t_left_s<<'\t'<<bracket.t_right_s<<'\t'<<(bracket.initial?"initial":"bracket")<<'\n';
    const auto& stats=main.statistics;const auto& last=main.accepted_steps.back();
    std::ofstream steps(root/"main.steps");steps<<std::setprecision(17)
      <<"accepted\trejected\trhs\tminimum_step_s\tmaximum_step_s\trows\n"
      <<stats.accepted_steps<<'\t'<<stats.rejected_steps<<'\t'<<stats.rhs_evaluations<<'\t'
      <<stats.minimum_step_s<<'\t'<<stats.maximum_step_s<<'\t'<<trajectory_rows
      <<"\nt_s\tcumulative_accepted\tcumulative_rejected\tlast_step_s\tminimum_step_since_previous_output_s\n"
      <<main.final_time_s<<'\t'<<stats.accepted_steps<<'\t'<<stats.rejected_steps<<'\t'
      <<last.accepted_step_s<<'\t'<<stats.minimum_step_s<<'\n';
    std::ofstream audit(root/"main.audit.tsv");audit<<std::setprecision(17)
      <<"gsl_apply_calls\tunique_positive_t1_targets\tt1_s\tintermediate_observation_t1_matches"
      <<"\tobserver_callbacks\trequested_observations\tobservation_records\tpositive_brackets\n"
      <<stats.gsl_apply_calls<<'\t'<<stats.distinct_positive_t1_targets<<'\t'<<stats.positive_t1_s
      <<"\t0\t"<<schedule.AcceptedCallbacks()<<'\t'<<schedule.RequestedTimes().size()<<'\t'
      <<schedule.Brackets().size()<<'\t'<<schedule.Brackets().size()-1<<'\n';
}

void WriteCheckpointEvidence(const std::filesystem::path& path,
  const std::vector<BNV::CheckpointOutput>& outputs,const std::vector<MatrixRow>& rows)
{
    std::ofstream out(path);require(bool(out),"cannot create checkpoint evidence");out<<std::setprecision(17);
    out<<"observation_index\tt_obs_s\tcategory\tknots\tsource\tstatus\trk8pd_invocations\tself_qualified"
      <<"\tx_O1\teta_e_O1\teta_mu_O1\tx_O2\teta_e_O2\teta_mu_O2"
      <<"\td_x\td_eta_e\td_eta_mu\tD_O1_x\tD_O1_eta_e\tD_O1_eta_mu"
      <<"\tU_x\tU_eta_e\tU_eta_mu\tF_x\tF_eta_e\tF_eta_mu"
      <<"\tdiagnostic_self_qualified\tP_dir_eq_O1\tP_dir_eq_O2\tP_dir_actual_O1\tP_dir_actual_O2"
      <<"\tLH_O1\tLH_O2\tDeltaLnu_O1\tDeltaLnu_O2\tDeltaPbeta_O1\tDeltaPbeta_O2"
      <<"\tLnu_eq_O1\tLnu_eq_O2\tLnu_full_O1\tLnu_full_O2\tLgamma_O1\tLgamma_O2"
      <<"\tLother_O1\tLother_O2\tPnet_O1\tPnet_O2\tmu_B_O1\tmu_B_O2"
      <<"\tmu_n_actual_O1\tmu_n_actual_O2\tsigma_e_O1\tsigma_e_O2\tsigma_mu_O1\tsigma_mu_O2"
      <<"\tEchem_O1\tEchem_O2\tCstar_O1\tCstar_O2\tTinf_O1\tTinf_O2\tB_O1\tB_O2"
      <<"\trun_card\tsource_identity\tdomain_identity\trevision_identity\tpartition_identity\tproduct_fate_identity"
      <<"\tO1_context_wall_s\tO1_context_cpu_s\tO1_solve_wall_s\tO1_solve_cpu_s"
      <<"\tO2_context_wall_s\tO2_context_cpu_s\tO2_solve_wall_s\tO2_solve_cpu_s"
      <<"\tdiagnostic_context_wall_s\tdiagnostic_context_cpu_s\tdiagnostic_wall_s\tdiagnostic_cpu_s\n";
    for(std::size_t i=0;i<outputs.size();++i)
    {
        const auto& q=outputs[i];const char category=i?rows[i-1].category:'A';const int knots=i?rows[i-1].knots:0;
        const auto source=q.source==BNV::CheckpointSource::MainEndpoint?"MAIN_ENDPOINT":"RK8PD_RECONSTRUCTED";
        const auto status=q.status==BNV::CheckpointStatus::Qualified?"QUALIFIED":"NUMERICALLY_UNRESOLVED";
        const auto& a=q.O1_diagnostics;const auto& b=q.O2_diagnostics;const auto& p=q.performance;
        out<<q.observation_index<<'\t'<<q.t_observation_s<<'\t'<<category<<'\t'<<knots<<'\t'<<source<<'\t'<<status
          <<'\t'<<q.provenance.rk8pd_invocations<<'\t'<<q.qualification.self_qualified;
        for(double v:q.O1_state)out<<'\t'<<v;for(double v:q.O2_state)out<<'\t'<<v;
        for(double v:q.qualification.d_O)out<<'\t'<<v;for(double v:q.qualification.D_O1)out<<'\t'<<v;
        for(double v:q.qualification.U_O)out<<'\t'<<v;
        for(double v:q.qualification.F_i)out<<'\t'<<v;
        out<<'\t'<<q.diagnostic_self_qualified
          <<'\t'<<a.P_dir_eq_erg_s<<'\t'<<b.P_dir_eq_erg_s<<'\t'<<a.P_dir_actual_erg_s<<'\t'<<b.P_dir_actual_erg_s
          <<'\t'<<a.LH_erg_s<<'\t'<<b.LH_erg_s<<'\t'<<a.DeltaLnu_erg_s<<'\t'<<b.DeltaLnu_erg_s
          <<'\t'<<a.DeltaPbeta_erg_s<<'\t'<<b.DeltaPbeta_erg_s<<'\t'<<a.Lnu_eq_erg_s<<'\t'<<b.Lnu_eq_erg_s
          <<'\t'<<a.Lnu_full_erg_s<<'\t'<<b.Lnu_full_erg_s<<'\t'<<a.Lgamma_erg_s<<'\t'<<b.Lgamma_erg_s
          <<'\t'<<a.Lother_erg_s<<'\t'<<b.Lother_erg_s<<'\t'<<a.Pnet_erg_s<<'\t'<<b.Pnet_erg_s
          <<'\t'<<a.mu_B_inf_MeV<<'\t'<<b.mu_B_inf_MeV<<'\t'<<a.mu_n_actual_inf_MeV<<'\t'<<b.mu_n_actual_inf_MeV
          <<'\t'<<a.sigma_count_s[0]<<'\t'<<b.sigma_count_s[0]<<'\t'<<a.sigma_count_s[1]<<'\t'<<b.sigma_count_s[1]
          <<'\t'<<a.Echem_MeV<<'\t'<<b.Echem_MeV<<'\t'<<a.Cstar_erg_K<<'\t'<<b.Cstar_erg_K
          <<'\t'<<a.Tinf_K<<'\t'<<b.Tinf_K<<'\t'<<a.B_count<<'\t'<<b.B_count
          <<'\t'<<b.run_card_identity<<'\t'<<b.source_identity<<'\t'<<b.domain_identity<<'\t'<<b.revision_identity
          <<'\t'<<b.partition_identity<<'\t'<<b.product_fate_identity
          <<'\t'<<p.O1_context_wall_seconds<<'\t'<<p.O1_context_cpu_seconds
          <<'\t'<<p.O1_solve_wall_seconds<<'\t'<<p.O1_solve_cpu_seconds
          <<'\t'<<p.O2_context_wall_seconds<<'\t'<<p.O2_context_cpu_seconds
          <<'\t'<<p.O2_solve_wall_seconds<<'\t'<<p.O2_solve_cpu_seconds
          <<'\t'<<p.diagnostic_context_wall_seconds<<'\t'<<p.diagnostic_context_cpu_seconds
          <<'\t'<<p.diagnostic_evaluation_wall_seconds<<'\t'<<p.diagnostic_evaluation_cpu_seconds<<'\n';
    }
}

void WriteSummary(const std::filesystem::path& path,const BNV::MainTrajectoryResult& main,
  const std::vector<BNV::CheckpointOutput>& outputs,const BNV::CheckpointR20Result& R20,
  double physics_context_wall,double physics_context_cpu,std::size_t context_count)
{
    double O1_wall=0,O1_cpu=0,O2_wall=0,O2_cpu=0,context_wall=0,context_cpu=0,diagnostic_wall=0,diagnostic_cpu=0;
    for(const auto& q:outputs)
    {
        O1_wall+=q.performance.O1_solve_wall_seconds;O1_cpu+=q.performance.O1_solve_cpu_seconds;
        O2_wall+=q.performance.O2_solve_wall_seconds;O2_cpu+=q.performance.O2_solve_cpu_seconds;
        context_wall+=q.performance.O1_context_wall_seconds+q.performance.O2_context_wall_seconds
          +q.performance.diagnostic_context_wall_seconds;
        context_cpu+=q.performance.O1_context_cpu_seconds+q.performance.O2_context_cpu_seconds
          +q.performance.diagnostic_context_cpu_seconds;
        diagnostic_wall+=q.performance.diagnostic_evaluation_wall_seconds;
        diagnostic_cpu+=q.performance.diagnostic_evaluation_cpu_seconds;
    }
    const double checkpoint_wall=O1_wall+O2_wall+context_wall+diagnostic_wall;
    std::ofstream out(path);out<<std::setprecision(17)<<"key\tvalue\n"
      <<"main_wall_s\t"<<main.statistics.wall_seconds<<"\nmain_cpu_s\t"<<main.statistics.cpu_seconds
      <<"\nphysics_context_wall_s\t"<<physics_context_wall<<"\nphysics_context_cpu_s\t"<<physics_context_cpu
      <<"\nO1_wall_s\t"<<O1_wall<<"\nO1_cpu_s\t"<<O1_cpu
      <<"\nO2_wall_s\t"<<O2_wall<<"\nO2_cpu_s\t"<<O2_cpu
      <<"\ncontext_wall_s\t"<<context_wall<<"\ncontext_cpu_s\t"<<context_cpu
      <<"\ndiagnostic_wall_s\t"<<diagnostic_wall<<"\ndiagnostic_cpu_s\t"<<diagnostic_cpu
      <<"\ncheckpoint_output_wall_s\t"<<checkpoint_wall<<"\ncontext_count\t"<<context_count
      <<"\nprocess_concurrency\t1\nstrict_interior\t239\nmain_run_count\t1"
      <<"\nR20_residual_erg\t"<<R20.R20_residual_erg<<"\nN_R20_erg\t"<<R20.N_R20_erg
      <<"\nR20_normalized\t"<<R20.R20_normalized
      <<"\nR20_reconstruction_uncertainty_erg\t"<<R20.reconstruction_uncertainty_erg
      <<"\n8192_serial_wall_estimate_s\t"<<checkpoint_wall*8191.0/239.0
      <<"\n8192_process_parallel_2_wall_estimate_s\t"<<checkpoint_wall*8191.0/(239.0*2.0)<<'\n';
}
}

int main(int argc,char** argv)
{
    try
    {
        gsl_set_error_handler_off();
        require(argc==11,"matrix profile certificate thermal work entry frozen coefficients output oracle-qualified-flag");
        Inputs inputs{argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],argv[9],argv[10]};
        require(!std::filesystem::exists(inputs.work),"fresh qualification context work root required");
        require(!std::filesystem::exists(inputs.output),"fresh qualification output root required");
        {std::ifstream flag(inputs.oracle_flag);std::string text((std::istreambuf_iterator<char>(flag)),{});
         require(text=="ADR0017_ORACLE_AUTHENTICATED\n","historical oracle is not authenticated");}
        const auto rows=ReadMatrix(inputs.matrix);std::vector<double> times{0};
        for(const auto& row:rows)times.push_back(row.t_obs);
        std::filesystem::create_directories(inputs.output);
        const auto physics_wall_begin=Clock::now();const auto physics_cpu_begin=std::clock();
        auto physics=BuildContext(inputs);
        const double physics_context_wall=std::chrono::duration<double>(Clock::now()-physics_wall_begin).count();
        const double physics_context_cpu=static_cast<double>(std::clock()-physics_cpu_begin)/CLOCKS_PER_SEC;

        RunState main_state(physics->OrdinaryContext());main_state.thermal.SetTinf(1.0e8);
        RC::ChemicalImbalanceState(0,0).Store(main_state.chem);
        auto driver=std::make_shared<BNV::ControlledBnvSecularDriver>(physics,BNV::ControlledBnvSecularDriver::Mode::Coupled);
        EV::EvolutionSystem system(main_state.ctx,main_state.state,main_state.rhs,main_state.layout,{driver});
        auto capture=std::make_shared<Capture>(inputs.output/"main.tsv",*driver);system.AddObserver(capture);
        BNV::PassiveObservationSchedule schedule(times,
          "43ec23ada72bfa59c4e89672ae2987e9927164db765f05afb7b45c924cc0c67e");
        BNV::MainIntegrationProvenance provenance{"CPL-P2-LINEAR-QSS-v1",
          "Phase-6A-1 controlled abstract uniform proper neutron sink v1",
          "phase6a1-frozen-currentness:full-before-and-after;cheap-per-step",
          schedule.Identity()};
        BNV::UninterruptedBnvTrajectory integrator(system,main_state.layout,physics->OrdinaryContext(),{},provenance);
        const auto main=integrator.Integrate(0,462269531250.0,{{0,0,0}},schedule,true);
        require(capture->rows==2,"bounded main trajectory did not emit start/final only");
        ValidateAgainstMatrix(main,schedule,rows);
        BNV::PassiveObservationSchedule alternate({0,times[1],times[120],times.back()},"alternate-postprocess-only");
        alternate.Replay(main.initial_state,main.accepted_steps);
        require(alternate.AcceptedCallbacks()==main.accepted_steps.size(),"alternate schedule replay incomplete");
        WriteMainEvidence(inputs.output,main,schedule,capture->rows);

        auto context_count=std::make_shared<std::size_t>(0);
        BNV::CheckpointContextFactory factory=[physics,context_count]
        {return std::make_unique<ActualCheckpointContext>(physics,++*context_count);};
        BNV::Rk8pdCheckpointReconstructor reconstructor;
        std::vector<BNV::CheckpointOutput> outputs;outputs.reserve(schedule.Brackets().size());
        for(const auto& bracket:schedule.Brackets())
        {
            auto checkpoint=reconstructor.Reconstruct(bracket,main.accepted_steps,factory);
            require(checkpoint.status==BNV::CheckpointStatus::Qualified,"production state reconstruction unresolved");
            reconstructor.EvaluateDiagnostics(checkpoint,factory);
            require(checkpoint.status==BNV::CheckpointStatus::Qualified&&checkpoint.diagnostics_evaluated,
                    "production diagnostic reconstruction unresolved");
            outputs.push_back(std::move(checkpoint));
        }
        require(std::count_if(outputs.begin()+1,outputs.end(),[](const auto& q)
          {return q.source==BNV::CheckpointSource::MainEndpoint;})==1,"positive exact-endpoint count changed");
        require(std::count_if(outputs.begin(),outputs.end(),[](const auto& q)
          {return q.source==BNV::CheckpointSource::Rk8pdReconstructed;})==239,"strict-interior count changed");
        const auto R20=BNV::ComputeCheckpointR20(outputs);
        WriteCheckpointEvidence(inputs.output/"checkpoints.tsv",outputs,rows);
        WriteSummary(inputs.output/"performance.tsv",main,outputs,R20,
          physics_context_wall,physics_context_cpu,*context_count);
        std::cout<<std::setprecision(17)
          <<"ADR0017_PRODUCTION_EXECUTION PASS main_accepted "<<main.statistics.accepted_steps
          <<" main_rejected "<<main.statistics.rejected_steps<<" contexts "<<*context_count
          <<" R20 "<<R20.R20_residual_erg<<" N_R20 "<<R20.N_R20_erg<<'\n';
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<"STOP "<<error.what()<<'\n';return 1;
    }
}
