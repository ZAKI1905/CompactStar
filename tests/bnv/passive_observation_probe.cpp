#define PHASE6A1_BA12R_ULTRA_SUPPORT_ONLY
#include "coupled_trajectory.cpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace
{
struct Endpoint
{
    std::size_t archived_index;
    double time_s;
};

struct InternalStep
{
    std::size_t sequence;
    std::size_t archived_checkpoint_index;
    double ceiling_s;
    double t_before_s;
    double t_after_s;
    double step_s;
    double suggested_next_h_s;
    std::size_t cumulative_rejected;
    double x_before;
    double x_after;
    int cstar_cell_before;
    int cstar_cell_after;
    int cstar_knots_crossed;
};

struct AcceptedState
{
    std::size_t sequence;
    double t_before_s;
    double t_after_s;
    double step_s;
    double suggested_next_h_s;
    std::size_t cumulative_rejected;
    double x_state;
    double eta_e_MeV;
    double eta_mu_MeV;
};

struct ObservationRecord
{
    std::size_t observation_index;
    double requested_t_s;
    std::size_t previous_accepted_step;
    std::size_t new_accepted_step;
    double t_previous_s;
    double t_new_s;
    const char* kind;
};

struct ProbeSummary
{
    std::size_t accepted=0,rejected=0,rhs=0,rows=0;
    double minimum_step_s=std::numeric_limits<double>::infinity();
    double maximum_step_s=0;
    std::vector<RC::StepOutput> outputs;
    std::vector<InternalStep> internal_steps;
    std::vector<AcceptedState> accepted_states;
};

int CstarCell(double x)
{
    const double T_MeV=1.0e8*std::exp(x)*RC::BoltzmannMeVPerK;
    const double u=(std::log(T_MeV)-std::log(1.0e-5))/(std::log(1.0)-std::log(1.0e-5))*159.0;
    return std::max(0,std::min(158,static_cast<int>(std::floor(u))));
}

std::vector<Endpoint> ReadEndpoints(const std::filesystem::path& path)
{
    std::ifstream in(path);
    if(!in)throw std::runtime_error("cannot open authenticated observation schedule");
    std::string line;
    if(!std::getline(in,line)||line!="index\tt_s")throw std::runtime_error("invalid observation-schedule header");
    std::vector<Endpoint> result;
    while(std::getline(in,line))
    {
        if(line.empty())continue;
        std::istringstream row(line);
        Endpoint endpoint{};
        if(!(row>>endpoint.archived_index>>endpoint.time_s))throw std::runtime_error("invalid observation-schedule row");
        result.push_back(endpoint);
    }
    if(result.size()!=241||result.front().archived_index!=0||result.front().time_s!=0||
       result.back().archived_index!=240||result.back().time_s!=462269531250.0)
        throw std::runtime_error("authenticated observation schedule is not indices 0..240");
    double previous=-1;
    for(std::size_t i=0;i<result.size();++i)
    {
        if(result[i].archived_index!=i||!(result[i].time_s>previous))
            throw std::runtime_error("observation schedule ordering changed");
        previous=result[i].time_s;
    }
    return result;
}

class PassiveObservationSchedule final
{
  public:
    explicit PassiveObservationSchedule(std::vector<Endpoint> schedule):schedule_(std::move(schedule))
    {
        if(schedule_.size()!=241||schedule_.front().time_s!=0)
            throw std::runtime_error("invalid passive observation schedule");
        records_.reserve(schedule_.size());
        records_.push_back({0,0,0,0,0,0,"initial"});
    }

    void NotifyAccepted(double previous,double current,std::size_t accepted_step)
    {
        ++callbacks_;
        if(callbacks_!=accepted_step||!(current>previous))
            throw std::runtime_error("invalid accepted-step notification");
        while(next_<schedule_.size()&&schedule_[next_].time_s<=current)
        {
            const auto& observation=schedule_[next_];
            if(!(observation.time_s>previous))
                throw std::runtime_error("passive observer missed or duplicated an observation");
            records_.push_back({observation.archived_index,observation.time_s,
              accepted_step-1,accepted_step,previous,current,"bracket"});
            ++next_;
        }
    }

    void RequireComplete()const
    {
        if(next_!=schedule_.size()||records_.size()!=schedule_.size())
            throw std::runtime_error("passive observation schedule incomplete");
    }

    const std::vector<ObservationRecord>& Records()const{return records_;}
    std::size_t CallbackCount()const{return callbacks_;}
    std::size_t RequestedCount()const{return schedule_.size();}
    std::size_t PositiveBracketCount()const{return records_.size()-1;}

  private:
    const std::vector<Endpoint> schedule_;
    std::vector<ObservationRecord> records_;
    std::size_t next_=1;
    std::size_t callbacks_=0;
};

class GslTargetAudit final
{
  public:
    GslTargetAudit(double expected,const std::vector<Endpoint>& schedule)
      :expected_(expected),schedule_(schedule){}

    void Record(double target)
    {
        if(target!=expected_)throw std::runtime_error("unexpected GSL t1 target");
        ++apply_calls_;
    }

    void RequireComplete()const
    {
        if(apply_calls_==0||UniquePositiveTargetCount()!=1||IntermediateScheduleMatches()!=0)
            throw std::runtime_error("GSL target audit failed");
    }

    std::size_t ApplyCalls()const{return apply_calls_;}
    std::size_t UniquePositiveTargetCount()const{return apply_calls_?1:0;}
    double Target()const{return expected_;}
    std::size_t IntermediateScheduleMatches()const
    {
        std::size_t matches=0;
        for(std::size_t i=1;i+1<schedule_.size();++i)if(schedule_[i].time_s==expected_)++matches;
        return matches;
    }

  private:
    const double expected_;
    const std::vector<Endpoint>& schedule_;
    std::size_t apply_calls_=0;
};

class PassiveAuditedRKF45 final
{
  public:
    PassiveAuditedRKF45(EV::EvolutionSystem& system,const EV::StateLayout& layout,
      std::shared_ptr<const RC::FrozenRotochemicalRunContext> context,RC::ComponentTolerances tolerances)
      :system_(system),context_(std::move(context)),tolerances_(tolerances)
    {
        tolerances_.Validate();
        if(!context_||layout.TotalSize()!=3||layout.BlockSize(Tag::Thermal)!=1||
           layout.BlockSize(Tag::Chem)!=2||layout.Offset(Tag::Thermal)!=0||layout.Offset(Tag::Chem)!=1)
            throw std::runtime_error("passive audit RKF45 requires thermal,Npe,NpMu layout");
    }

    void Derivative(double t,const double* y,double* out)
    {
        context_->RequireCheapCurrent();
        if(stats_)++stats_->rhs;
        std::array<double,3> local;
        const int rc=system_(t,y,local.data());
        if(rc)throw std::runtime_error("EvolutionSystem RHS failure");
        context_->RequireCheapCurrent();
        for(double value:local)if(!std::isfinite(value))throw std::runtime_error("nonfinite derivative");
        std::copy(local.begin(),local.end(),out);
    }

    void Integrate(double start,double final_time_s,double* y,ProbeSummary& stats,
      PassiveObservationSchedule& observer,GslTargetAudit& target_audit,bool enforce_thermal_domain=true)
    {
        context_->RequireFullCurrent();
        if(!std::isfinite(start)||!std::isfinite(final_time_s)||!(final_time_s>start))
            throw std::runtime_error("invalid passive integration interval");
        stats={};stats.internal_steps.reserve(512);stats.accepted_states.reserve(512);
        stats_=&stats;failure_=nullptr;
        std::unique_ptr<gsl_odeiv2_step,decltype(&gsl_odeiv2_step_free)> step(
          gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkf45,3),gsl_odeiv2_step_free);
        std::unique_ptr<gsl_odeiv2_control,decltype(&gsl_odeiv2_control_free)> control(
          tolerances_.Allocate(),gsl_odeiv2_control_free);
        std::unique_ptr<gsl_odeiv2_evolve,decltype(&gsl_odeiv2_evolve_free)> evolve(
          gsl_odeiv2_evolve_alloc(3),gsl_odeiv2_evolve_free);
        if(!step||!evolve)throw std::bad_alloc();
        gsl_odeiv2_system sys{Callback,nullptr,3,this};
        double t=start,h=std::min(1.0,(final_time_s-start)*1.0e-3),last=0;
        if(h!=1.0)throw std::runtime_error("initial h is not exactly 1 s");
        try
        {
            ValidateAccepted(y,enforce_thermal_domain);
            system_.NotifyStart(start,final_time_s,y);
            std::size_t count=0;
            double interval_minimum_step=std::numeric_limits<double>::infinity();
            while(t<final_time_s)
            {
                if(++count>100000)throw std::runtime_error("RKF45 internal step budget exceeded");
                context_->RequireCheapCurrent();
                const double before=t;
                const double x_before=y[0];
                target_audit.Record(final_time_s);
                const int rc=gsl_odeiv2_evolve_apply(evolve.get(),control.get(),step.get(),&sys,
                  &t,final_time_s,&h,y);
                stats.rejected=evolve->failed_steps;
                if(failure_)std::rethrow_exception(failure_);
                if(rc!=GSL_SUCCESS)throw std::runtime_error(std::string("RKF45 failure: ")+gsl_strerror(rc));
                last=t-before;
                if(!(last>0)||!std::isfinite(last))throw std::runtime_error("RKF45 failed progress");
                ++stats.accepted;
                stats.minimum_step_s=std::min(stats.minimum_step_s,last);
                stats.maximum_step_s=std::max(stats.maximum_step_s,last);
                interval_minimum_step=std::min(interval_minimum_step,last);
                ValidateAccepted(y,enforce_thermal_domain);
                std::array<double,3> endpoint_rhs;
                Derivative(t,y,endpoint_rhs.data());
                const int before_cell=CstarCell(x_before);
                const int after_cell=CstarCell(y[0]);
                stats.internal_steps.push_back({stats.accepted,240,final_time_s,before,t,last,h,
                  stats.rejected,x_before,y[0],before_cell,after_cell,std::abs(after_cell-before_cell)});
                stats.accepted_states.push_back({stats.accepted,before,t,last,h,stats.rejected,
                  y[0],y[1],y[2]});
                observer.NotifyAccepted(before,t,stats.accepted);
            }
            stats.outputs.push_back({t,stats.accepted,stats.rejected,last,interval_minimum_step});
            system_.NotifySample(t,y,0);
            system_.NotifyFinish(t,y,true);
        }
        catch(...)
        {
            auto original=std::current_exception();stats_=nullptr;
            try{system_.NotifyFinish(t,y,false);}catch(...){}
            std::rethrow_exception(original);
        }
        stats_=nullptr;
        observer.RequireComplete();
        target_audit.RequireComplete();
    }

  private:
    static int Callback(double t,const double* y,double* out,void* p)noexcept
    {
        auto& self=*static_cast<PassiveAuditedRKF45*>(p);
        try{self.Derivative(t,y,out);return GSL_SUCCESS;}
        catch(...){self.failure_=std::current_exception();return GSL_EBADFUNC;}
    }

    void ValidateAccepted(const double* y,bool thermal)const
    {
        context_->RequireCheapCurrent();
        for(std::size_t i=0;i<3;++i)if(!std::isfinite(y[i]))throw std::runtime_error("nonfinite accepted state");
        const double T=1.0e8*std::exp(y[0]);
        if(!(T>0)||!std::isfinite(T))throw std::runtime_error("invalid accepted temperature");
        if(thermal&&!(T*RC::BoltzmannMeVPerK>1.0e-5&&T*RC::BoltzmannMeVPerK<1.0))
            throw std::runtime_error("accepted state outside predeclared thermal cache domain");
    }

    EV::EvolutionSystem& system_;
    const std::shared_ptr<const RC::FrozenRotochemicalRunContext> context_;
    const RC::ComponentTolerances tolerances_;
    ProbeSummary* stats_=nullptr;
    std::exception_ptr failure_;
};

void WriteEvidence(const std::filesystem::path& output,const ProbeSummary& summary,
  const PassiveObservationSchedule& observer,const GslTargetAudit& target_audit)
{
    std::ofstream steps(output.string()+".steps");
    steps<<std::setprecision(17)
      <<"accepted\trejected\trhs\tminimum_step_s\tmaximum_step_s\trows\n"
      <<summary.accepted<<'\t'<<summary.rejected<<'\t'<<summary.rhs<<'\t'<<summary.minimum_step_s<<'\t'
      <<summary.maximum_step_s<<'\t'<<summary.rows
      <<"\nt_s\tcumulative_accepted\tcumulative_rejected\tlast_step_s\tminimum_step_since_previous_output_s\n";
    for(const auto& sample:summary.outputs)
        steps<<sample.time_s<<'\t'<<sample.accepted<<'\t'<<sample.rejected<<'\t'<<sample.last_step_s<<'\t'
             <<sample.minimum_step_since_previous_output_s<<'\n';
    if(!steps)throw std::runtime_error("step-summary serialization failure");

    std::ofstream internal(output.string()+".internal_steps.tsv");
    internal<<std::setprecision(17)
      <<"sequence\tarchived_checkpoint_index\tceiling_s\tt_before_s\tt_after_s\tstep_s\tsuggested_next_h_s"
      <<"\tcumulative_rejected\tx_before\tx_after\tcstar_cell_before\tcstar_cell_after\tcstar_knots_crossed\n";
    for(const auto& sample:summary.internal_steps)
        internal<<sample.sequence<<'\t'<<sample.archived_checkpoint_index<<'\t'<<sample.ceiling_s<<'\t'
          <<sample.t_before_s<<'\t'<<sample.t_after_s<<'\t'<<sample.step_s<<'\t'<<sample.suggested_next_h_s<<'\t'
          <<sample.cumulative_rejected<<'\t'<<sample.x_before<<'\t'<<sample.x_after<<'\t'
          <<sample.cstar_cell_before<<'\t'<<sample.cstar_cell_after<<'\t'<<sample.cstar_knots_crossed<<'\n';
    if(!internal)throw std::runtime_error("internal-step serialization failure");

    std::ofstream accepted(output.string()+".accepted_states.tsv");
    accepted<<std::setprecision(17)
      <<"sequence\tt_before_s\tt_after_s\tstep_s\tsuggested_next_h_s\tcumulative_rejected"
      <<"\tx_state\teta_e_MeV\teta_mu_MeV\n";
    for(const auto& sample:summary.accepted_states)
        accepted<<sample.sequence<<'\t'<<sample.t_before_s<<'\t'<<sample.t_after_s<<'\t'<<sample.step_s<<'\t'
          <<sample.suggested_next_h_s<<'\t'<<sample.cumulative_rejected<<'\t'<<sample.x_state<<'\t'
          <<sample.eta_e_MeV<<'\t'<<sample.eta_mu_MeV<<'\n';
    if(!accepted)throw std::runtime_error("accepted-state serialization failure");

    std::ofstream observations(output.string()+".observations.tsv");
    observations<<std::setprecision(17)
      <<"observation_index\trequested_t_s\tprevious_accepted_step\tnew_accepted_step\tt_previous_s\tt_new_s\tkind\n";
    for(const auto& sample:observer.Records())
        observations<<sample.observation_index<<'\t'<<sample.requested_t_s<<'\t'
          <<sample.previous_accepted_step<<'\t'<<sample.new_accepted_step<<'\t'
          <<sample.t_previous_s<<'\t'<<sample.t_new_s<<'\t'<<sample.kind<<'\n';
    if(!observations)throw std::runtime_error("observation serialization failure");

    std::ofstream audit(output.string()+".audit.tsv");
    audit<<std::setprecision(17)
      <<"gsl_apply_calls\tunique_positive_t1_targets\tt1_s\tintermediate_observation_t1_matches"
      <<"\tobserver_callbacks\trequested_observations\tobservation_records\tpositive_brackets\n"
      <<target_audit.ApplyCalls()<<'\t'<<target_audit.UniquePositiveTargetCount()<<'\t'<<target_audit.Target()<<'\t'
      <<target_audit.IntermediateScheduleMatches()<<'\t'<<observer.CallbackCount()<<'\t'
      <<observer.RequestedCount()<<'\t'<<observer.Records().size()<<'\t'<<observer.PositiveBracketCount()<<'\n';
    if(!audit)throw std::runtime_error("audit serialization failure");
}

ProbeSummary RunProbe(const std::shared_ptr<const BNV::FrozenControlledBnvRunContext>& context,
    const std::vector<Endpoint>& schedule,const std::filesystem::path& output,
    const RC::ComponentTolerances& tolerances)
{
    RunState state(context->OrdinaryContext());
    state.thermal.SetTinf(1.0e8);
    RC::ChemicalImbalanceState(0,0).Store(state.chem);
    auto y=state.Pack();
    auto driver=std::make_shared<BNV::ControlledBnvSecularDriver>(context,BNV::ControlledBnvSecularDriver::Mode::Coupled);
    EV::EvolutionSystem system(state.ctx,state.state,state.rhs,state.layout,{driver});
    auto capture=std::make_shared<Capture>(output,*driver);
    system.AddObserver(capture);
    PassiveObservationSchedule observer(schedule);
    GslTargetAudit target_audit(schedule.back().time_s,schedule);
    PassiveAuditedRKF45 solver(system,state.layout,context->OrdinaryContext(),tolerances);
    ProbeSummary summary;
    solver.Integrate(0,schedule.back().time_s,y.data(),summary,observer,target_audit,true);
    summary.rows=capture->rows;
    if(summary.rows!=2)throw std::runtime_error("passive trajectory row count changed");
    WriteEvidence(output,summary,observer,target_audit);
    return summary;
}
}

int main(int argc,char** argv)
{
 try
 {
    gsl_set_error_handler_off();
    require(argc==11,"schedule profile certificate thermal work entry frozen-certificate coefficients output pretrajectory-record");
    std::cout<<std::setprecision(17)<<std::unitbuf;
    const std::filesystem::path schedule_path=argv[1],profile=argv[2],certificate=argv[3],thermal_path=argv[4],
      work=argv[5],entry=argv[6],frozen=argv[7],coefficients=argv[8],output=argv[9],pretrajectory=argv[10];
    const std::array<std::string,5> suffixes{{"",".steps",".internal_steps.tsv",".accepted_states.tsv",".observations.tsv"}};
    require(!std::filesystem::exists(work),"fresh passive-observation assembly directory required");
    for(const auto& suffix:suffixes)require(!std::filesystem::exists(output.string()+suffix),"fresh passive-observation output path required");
    require(!std::filesystem::exists(output.string()+".audit.tsv"),"fresh passive-observation audit path required");
    {std::ifstream in(pretrajectory);std::string text((std::istreambuf_iterator<char>(in)),{});
     require(text.find("PRETRAJECTORY PASS")!=std::string::npos&&text.find("no BNV trajectory had been generated")!=std::string::npos,
       "committed pretrajectory authorization missing");}

    const auto& card=Campaign::RunCards().back();
    require(card.identity=="CPL-P2-LINEAR-QSS-v1"&&card.partition=="P2"&&
      card.fractional_drive_per_year==-1.0e-12&&card.duration_year==5.0e5&&
      card.checkpoints==8193&&!card.reaction_free,"immutable passive-observation P2 card changed");
    const RC::ComponentTolerances ultra{1.0e-11,{1.0e-16,1.0e-22,1.0e-22}};
    ultra.Validate();
    const auto schedule=ReadEndpoints(schedule_path);

    std::filesystem::create_directories(output.parent_path());
    auto fixture=Fixture(profile,certificate.string(),work/"owning-star",80000);
    auto channels=Channels(fixture);
    EOS::CompOSE_Thermo::Options options;options.Tmin_for_derivative_MeV=0;options.clamp_to_domain=false;
    auto thermal=std::make_shared<const RC::FrozenThermalSource>(thermal_path,options,
      "controlled mathematical fixed-background free-gas entropy; qualified radial80000");
    auto tangent=Campaign::Tangent(fixture,profile,work);
    auto bnv_token=std::make_shared<BNV::BnvDependencyToken>();
    auto monitor=Phase6A1Test::LoadFrozenMonitor(frozen,tangent,bnv_token);
    auto run_token=std::make_shared<RC::RunDependencyToken>();
    auto spin=std::make_shared<const BNV::StaticZeroSpinHistory>(run_token);
    auto ordinary=Context(fixture,channels,thermal,spin,Campaign::Qualification(profile,certificate,entry),run_token);
    const double mu_B=Campaign::Column(coefficients,"mu_B_inf");
    const std::string potential="frozen Structure-1 B0 equilibrium mu_B plus governed moving-reference actual-potential correction";
    const double Bdot=card.fractional_drive_per_year*tangent->B0Count()/Year;
    require(Bdot==-2.4136520263641375e37,"immutable passive-observation Bdot changed");
    const std::string fate_id="phase6a1-P2-terminal-fate-v1";
    BNV::ProductFateLedger fate(fate_id,{{"controlled-terminal-product",BNV::TerminalProductFate::BoundInert,1,
      "abstract-neutron-disappearance"}});
    auto partition=Partition(card,fixture,bnv_token,fate_id);
    auto history=std::make_shared<const Phase6A1Test::UniformProperNeutronSinkHistory>(fixture.central,
      tangent->B0Count(),Bdot,tangent->DomainIdentity(),tangent->StarIdentity(),tangent->SequenceStateIdentity(),
      partition->Identity(),fate_id,bnv_token);
    auto run_context=std::make_shared<const BNV::FrozenControlledBnvRunContext>(ordinary,tangent,history,partition,
      fate,monitor,channels,mu_B,potential,card.identity);

    const auto summary=RunProbe(run_context,schedule,output,ultra);
    std::cout<<"PASSIVE_OBSERVATION_PROBE COMPLETE rows "<<summary.rows
             <<" accepted "<<summary.accepted<<" rejected "<<summary.rejected
             <<" rhs "<<summary.rhs<<" minimum_step_s "<<summary.minimum_step_s
             <<" maximum_step_s "<<summary.maximum_step_s<<'\n';
    return 0;
 }
 catch(const std::exception& e){std::cerr<<"STOP "<<e.what()<<'\n';return 1;}
}
