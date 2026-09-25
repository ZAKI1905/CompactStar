#include <CompactStar/Physics/BNV/PassiveCheckpointOutput.hpp>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>
#include <utility>

namespace CompactStar::Physics::BNV
{
namespace
{
using Clock=std::chrono::steady_clock;

void Need(bool condition,const char* message)
{
    if(!condition)throw std::runtime_error(message);
}

double WallSeconds(const Clock::time_point& begin)
{
    return std::chrono::duration<double>(Clock::now()-begin).count();
}

int CstarCell(double x)
{
    const double T_MeV=1.0e8*std::exp(x)*Rotochemical::BoltzmannMeVPerK;
    const double u=(std::log(T_MeV)-std::log(1.0e-5))
      /(std::log(1.0)-std::log(1.0e-5))*159.0;
    return std::max(0,std::min(158,static_cast<int>(std::floor(u))));
}

double Ulp(double magnitude)
{
    Need(std::isfinite(magnitude)&&magnitude>=0,"invalid ulp magnitude");
    const double next=std::nextafter(magnitude,std::numeric_limits<double>::infinity());
    const double ulp=next-magnitude;
    return ulp>0?ulp:std::numeric_limits<double>::denorm_min();
}

struct LocalRk8pdCallback
{
    ICheckpointEvaluationContext* context=nullptr;
    std::exception_ptr failure;
    static int Call(double t,const double* y,double* dydt,void* opaque) noexcept
    {
        auto& self=*static_cast<LocalRk8pdCallback*>(opaque);
        try
        {
            CheckpointState in{{y[0],y[1],y[2]}},out{};
            self.context->Derivative(t,in,out);
            for(std::size_t i=0;i<out.size();++i)
            {
                if(!std::isfinite(out[i]))throw std::runtime_error("nonfinite reconstruction derivative");
                dydt[i]=out[i];
            }
            return GSL_SUCCESS;
        }
        catch(...)
        {
            self.failure=std::current_exception();
            return GSL_EBADFUNC;
        }
    }
};

double DiagnosticScale(const BnvDiagnostics& a,const BnvDiagnostics& b)
{
    const auto power_sum=[](const BnvDiagnostics& d)
    {
        return std::abs(d.P_dir_actual_erg_s)+std::abs(d.LH_erg_s)
          +std::abs(d.DeltaLnu_erg_s)+std::abs(d.Lnu_eq_erg_s)
          +std::abs(d.Lgamma_erg_s)+std::abs(d.Lother_erg_s)
          +std::abs(d.L_out_fluid_inf_erg_s);
    };
    return std::max({1.0,power_sum(a),power_sum(b)});
}

std::map<std::string,double> DiagnosticValues(const BnvDiagnostics& d)
{
    return {
      {"P_dir_eq",d.P_dir_eq_erg_s},{"P_dir_actual",d.P_dir_actual_erg_s},
      {"LH",d.LH_erg_s},{"DeltaLnu",d.DeltaLnu_erg_s},{"DeltaPbeta",d.DeltaPbeta_erg_s},
      {"Lnu_eq",d.Lnu_eq_erg_s},{"Lnu_full",d.Lnu_full_erg_s},{"Lgamma",d.Lgamma_erg_s},
      {"Lother",d.Lother_erg_s},{"Pnet",d.Pnet_erg_s},{"mu_B",d.mu_B_inf_MeV},
      {"mu_n_actual",d.mu_n_actual_inf_MeV},{"sigma_e",d.sigma_count_s[0]},
      {"sigma_mu",d.sigma_count_s[1]},{"Echem",d.Echem_MeV}
    };
}

struct R20Terms
{
    double delta_eq=0,delta_chem=0,delta_uth=0,outgoing=0,residual=0,normalizer=0;
};

R20Terms R20From(const std::vector<CheckpointOutput>& checkpoints,bool O1)
{
    const auto diagnostic=[O1](const CheckpointOutput& c)->const BnvDiagnostics&
    {return O1?c.O1_diagnostics:c.O2_diagnostics;};
    const auto& first=diagnostic(checkpoints.front());
    const auto& last=diagnostic(checkpoints.back());
    R20Terms out;
    out.delta_eq=Rotochemical::MeVToErg*first.mu_B_inf_MeV*(last.B_count-first.B_count);
    out.delta_chem=Rotochemical::MeVToErg*(last.Echem_MeV-first.Echem_MeV);
    double outgoing_abs=0;
    for(std::size_t i=1;i<checkpoints.size();++i)
    {
        const auto& a=diagnostic(checkpoints[i-1]);
        const auto& b=diagnostic(checkpoints[i]);
        const double dt=checkpoints[i].t_observation_s-checkpoints[i-1].t_observation_s;
        Need(dt>0&&std::isfinite(dt),"unordered R20 checkpoints");
        out.delta_uth+=0.5*(a.Cstar_erg_K+b.Cstar_erg_K)*(b.Tinf_K-a.Tinf_K);
        const double pa=a.L_out_fluid_inf_erg_s+a.Lnu_full_erg_s+a.Lgamma_erg_s+a.Lother_erg_s;
        const double pb=b.L_out_fluid_inf_erg_s+b.Lnu_full_erg_s+b.Lgamma_erg_s+b.Lother_erg_s;
        out.outgoing+=0.5*dt*(pa+pb);
        const double aa=std::abs(a.Lnu_full_erg_s)+std::abs(a.Lgamma_erg_s)+std::abs(a.Lother_erg_s);
        const double ab=std::abs(b.Lnu_full_erg_s)+std::abs(b.Lgamma_erg_s)+std::abs(b.Lother_erg_s);
        outgoing_abs+=0.5*dt*(aa+ab);
    }
    out.residual=out.delta_eq+out.delta_chem+out.delta_uth+out.outgoing;
    out.normalizer=std::max({1.0,std::abs(out.delta_uth),std::abs(out.delta_chem),outgoing_abs});
    return out;
}
} // namespace

void Rkf45Configuration::Validate() const
{
    Need(relative>0&&std::isfinite(relative),"invalid RKF45 relative tolerance");
    Need(initial_step_s>0&&std::isfinite(initial_step_s),"invalid RKF45 initial step");
    for(double value:absolute)Need(value>0&&std::isfinite(value),"invalid RKF45 absolute tolerance");
}

void Rk8pdConfiguration::Validate() const
{
    Need(!level.empty(),"missing rk8pd level identity");
    Need(relative>0&&std::isfinite(relative),"invalid rk8pd relative tolerance");
    for(double value:absolute)Need(value>0&&std::isfinite(value),"invalid rk8pd absolute tolerance");
}

PassiveObservationSchedule::PassiveObservationSchedule(
  std::vector<double> requested_times_s,std::string identity)
  :requested_times_s_(std::move(requested_times_s)),identity_(std::move(identity))
{
    Need(!requested_times_s_.empty()&&!identity_.empty(),"missing observation schedule");
    Need(std::isfinite(requested_times_s_.front()),"nonfinite first observation");
    for(std::size_t i=1;i<requested_times_s_.size();++i)
        Need(std::isfinite(requested_times_s_[i])&&requested_times_s_[i]>requested_times_s_[i-1],
             "observation schedule must be strictly ordered without duplicates");
}

void PassiveObservationSchedule::Reset()
{
    brackets_.clear();
    next_=0;
    accepted_callbacks_=0;
}

void PassiveObservationSchedule::NotifyInitial(const CheckpointState& initial_state)
{
    Reset();
    for(double value:initial_state)Need(std::isfinite(value),"nonfinite initial observation state");
    Need(requested_times_s_.front()==0.0,"observation schedule must begin at the main initial time");
    brackets_.push_back({0,requested_times_s_.front(),0,0,requested_times_s_.front(),requested_times_s_.front(),true});
    next_=1;
}

void PassiveObservationSchedule::NotifyAccepted(const AcceptedStepRecord& step)
{
    ++accepted_callbacks_;
    Need(step.ordinal==accepted_callbacks_,"accepted-step notifications are not contiguous");
    Need(step.t_right_s>step.t_left_s,"nonpositive accepted step");
    while(next_<requested_times_s_.size()&&requested_times_s_[next_]<=step.t_right_s)
    {
        const double requested=requested_times_s_[next_];
        Need(requested>step.t_left_s,"passive observer missed or duplicated an observation");
        brackets_.push_back({next_,requested,step.ordinal-1,step.ordinal,
                             step.t_left_s,step.t_right_s,false});
        ++next_;
    }
}

void PassiveObservationSchedule::RequireComplete() const
{
    Need(next_==requested_times_s_.size(),"passive observation schedule incomplete");
    Need(brackets_.size()==requested_times_s_.size(),"observation bracket count mismatch");
}

void PassiveObservationSchedule::Replay(
  const CheckpointState& initial_state,const std::vector<AcceptedStepRecord>& history)
{
    NotifyInitial(initial_state);
    for(const auto& step:history)NotifyAccepted(step);
    RequireComplete();
}

UninterruptedBnvTrajectory::UninterruptedBnvTrajectory(
  Evolution::EvolutionSystem& system,const Evolution::StateLayout& layout,
  std::shared_ptr<const Rotochemical::FrozenRotochemicalRunContext> currentness_owner,
  Rkf45Configuration configuration,MainIntegrationProvenance provenance)
  :system_(system),layout_(layout),currentness_owner_(std::move(currentness_owner)),
   configuration_(std::move(configuration)),provenance_(std::move(provenance))
{
    configuration_.Validate();
    Need(currentness_owner_!=nullptr,"missing main currentness owner");
    Need(layout_.TotalSize()==3&&layout_.BlockSize(State::StateTag::Thermal)==1
      &&layout_.BlockSize(State::StateTag::Chem)==2&&layout_.Offset(State::StateTag::Thermal)==0
      &&layout_.Offset(State::StateTag::Chem)==1,"main integration requires thermal,Npe,NpMu layout");
    Need(!provenance_.run_card_identity.empty()&&!provenance_.source_identity.empty()
      &&!provenance_.currentness_identity.empty()&&!provenance_.observation_schedule_identity.empty(),
      "incomplete main integration provenance");
}

int UninterruptedBnvTrajectory::Callback(double t,const double* y,double* out,void* opaque) noexcept
{
    auto& self=*static_cast<UninterruptedBnvTrajectory*>(opaque);
    try{self.Derivative(t,y,out);return GSL_SUCCESS;}
    catch(...){self.failure_=std::current_exception();return GSL_EBADFUNC;}
}

void UninterruptedBnvTrajectory::Derivative(double t,const double* y,double* out)
{
    currentness_owner_->RequireCheapCurrent();
    if(statistics_)++statistics_->rhs_evaluations;
    CheckpointState local{};
    const int rc=system_(t,y,local.data());
    Need(rc==0,"EvolutionSystem RHS failure");
    currentness_owner_->RequireCheapCurrent();
    for(std::size_t i=0;i<local.size();++i)
    {
        Need(std::isfinite(local[i]),"nonfinite main derivative");
        out[i]=local[i];
    }
}

void UninterruptedBnvTrajectory::ValidateAccepted(const double* y,bool thermal) const
{
    currentness_owner_->RequireCheapCurrent();
    for(std::size_t i=0;i<3;++i)Need(std::isfinite(y[i]),"nonfinite accepted main state");
    const double T=1.0e8*std::exp(y[0]);
    Need(T>0&&std::isfinite(T),"invalid accepted main temperature");
    if(thermal)Need(T*Rotochemical::BoltzmannMeVPerK>1.0e-5
      &&T*Rotochemical::BoltzmannMeVPerK<1.0,"accepted state outside thermal cache domain");
}

MainTrajectoryResult UninterruptedBnvTrajectory::Integrate(
  double start,double final_time,CheckpointState initial_state,
  PassiveObservationSchedule& observation_schedule,bool enforce_thermal_domain)
{
    currentness_owner_->RequireFullCurrent();
    configuration_.Validate();
    Need(std::isfinite(start)&&std::isfinite(final_time)&&final_time>start,"invalid main interval");
    Need(observation_schedule.Identity()==provenance_.observation_schedule_identity,
         "main/schedule identity mismatch");
    Need(observation_schedule.RequestedTimes().front()==start
      &&observation_schedule.RequestedTimes().back()==final_time,"schedule/main endpoints mismatch");
    MainTrajectoryResult result;
    result.start_time_s=start;result.final_time_s=final_time;
    result.initial_state=initial_state;result.configuration=configuration_;result.provenance=provenance_;
    result.accepted_steps.reserve(512);
    observation_schedule.NotifyInitial(initial_state);
    auto& statistics=result.statistics;
    statistics.positive_t1_s=final_time;
    statistics.distinct_positive_t1_targets=1;
    statistics.minimum_step_s=std::numeric_limits<double>::infinity();
    statistics_=&statistics;failure_=nullptr;
    std::unique_ptr<gsl_odeiv2_step,decltype(&gsl_odeiv2_step_free)> step(
      gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkf45,3),gsl_odeiv2_step_free);
    std::unique_ptr<gsl_odeiv2_control,decltype(&gsl_odeiv2_control_free)> control(
      gsl_odeiv2_control_scaled_new(1.0,configuration_.relative,1.0,0.0,
                                    configuration_.absolute.data(),3),gsl_odeiv2_control_free);
    std::unique_ptr<gsl_odeiv2_evolve,decltype(&gsl_odeiv2_evolve_free)> evolve(
      gsl_odeiv2_evolve_alloc(3),gsl_odeiv2_evolve_free);
    Need(step&&control&&evolve,"failed to allocate main GSL objects");
    gsl_odeiv2_system gsl_system{Callback,nullptr,3,this};
    double t=start,h=configuration_.initial_step_s;
    auto y=initial_state;
    const auto wall_begin=Clock::now();
    try
    {
        ValidateAccepted(y.data(),enforce_thermal_domain);
        system_.NotifyStart(start,final_time,y.data());
        while(t<final_time)
        {
            Need(statistics.accepted_steps<100000,"main RKF45 accepted-step budget exceeded");
            currentness_owner_->RequireCheapCurrent();
            const double before=t;
            const CheckpointState y_before=y;
            ++statistics.gsl_apply_calls;
            const int rc=gsl_odeiv2_evolve_apply(evolve.get(),control.get(),step.get(),
              &gsl_system,&t,final_time,&h,y.data());
            statistics.rejected_steps=evolve->failed_steps;
            if(failure_)std::rethrow_exception(failure_);
            Need(rc==GSL_SUCCESS,"main RKF45 GSL failure");
            const double accepted_step=t-before;
            Need(accepted_step>0&&std::isfinite(accepted_step),"main RKF45 failed progress");
            ++statistics.accepted_steps;
            statistics.minimum_step_s=std::min(statistics.minimum_step_s,accepted_step);
            statistics.maximum_step_s=std::max(statistics.maximum_step_s,accepted_step);
            ValidateAccepted(y.data(),enforce_thermal_domain);
            CheckpointState endpoint_rhs{};
            Derivative(t,y.data(),endpoint_rhs.data());
            const int left_cell=CstarCell(y_before[0]),right_cell=CstarCell(y[0]);
            AcceptedStepRecord record;
            record.ordinal=statistics.accepted_steps;record.t_left_s=before;record.y_left=y_before;
            record.t_right_s=t;record.y_right=y;record.accepted_step_s=accepted_step;
            record.suggested_next_h_s=h;record.cumulative_rejected=statistics.rejected_steps;
            record.source_identity=provenance_.source_identity;
            record.currentness_identity=provenance_.currentness_identity;
            record.cstar_cell_left=left_cell;record.cstar_cell_right=right_cell;
            record.cstar_knots_crossed=std::abs(right_cell-left_cell);
            result.accepted_steps.push_back(record);
            observation_schedule.NotifyAccepted(result.accepted_steps.back());
        }
        system_.NotifySample(t,y.data(),0);
        system_.NotifyFinish(t,y.data(),true);
        observation_schedule.RequireComplete();
        currentness_owner_->RequireFullCurrent();
    }
    catch(...)
    {
        statistics_=nullptr;
        try{system_.NotifyFinish(t,y.data(),false);}catch(...){}
        throw;
    }
    statistics.wall_seconds=WallSeconds(wall_begin);
    statistics_=nullptr;
    result.final_state=y;
    return result;
}

Rk8pdCheckpointReconstructor::Rk8pdCheckpointReconstructor()
  :O1_{"O1",1.0e-12,{{1.0e-17,1.0e-23,1.0e-23}}},
   O2_{"O2",1.0e-13,{{1.0e-18,1.0e-24,1.0e-24}}}
{
    O1_.Validate();O2_.Validate();
}

Rk8pdCheckpointReconstructor::SolveResult Rk8pdCheckpointReconstructor::Integrate(
  double start,const CheckpointState& initial,double target,const Rk8pdConfiguration& configuration,
  const CheckpointContextFactory& context_factory)
{
    configuration.Validate();
    Need(target>start&&std::isfinite(start)&&std::isfinite(target),"invalid strict-interior interval");
    Need(static_cast<bool>(context_factory),"missing reconstruction context factory");
    SolveResult result;
    auto begin=Clock::now();
    auto context=context_factory();
    result.context_wall_seconds=WallSeconds(begin);
    Need(context!=nullptr,"reconstruction factory returned null");
    context->RequireFullCurrent();
    result.context_identity=context->Identity();
    Need(!result.context_identity.empty(),"missing reconstruction context identity");
    std::unique_ptr<gsl_odeiv2_step,decltype(&gsl_odeiv2_step_free)> step(
      gsl_odeiv2_step_alloc(gsl_odeiv2_step_rk8pd,3),gsl_odeiv2_step_free);
    std::unique_ptr<gsl_odeiv2_control,decltype(&gsl_odeiv2_control_free)> control(
      gsl_odeiv2_control_scaled_new(1.0,configuration.relative,1.0,0.0,
                                    configuration.absolute.data(),3),gsl_odeiv2_control_free);
    std::unique_ptr<gsl_odeiv2_evolve,decltype(&gsl_odeiv2_evolve_free)> evolve(
      gsl_odeiv2_evolve_alloc(3),gsl_odeiv2_evolve_free);
    Need(step&&control&&evolve,"failed to allocate isolated rk8pd objects");
    LocalRk8pdCallback callback{context.get(),nullptr};
    gsl_odeiv2_system system{LocalRk8pdCallback::Call,nullptr,3,&callback};
    double t=start,h=target-start;
    auto y=initial;
    begin=Clock::now();
    std::size_t calls=0;
    while(t<target)
    {
        Need(++calls<=100000,"rk8pd local step budget exceeded");
        const int rc=gsl_odeiv2_evolve_apply(evolve.get(),control.get(),step.get(),
          &system,&t,target,&h,y.data());
        if(callback.failure)std::rethrow_exception(callback.failure);
        Need(rc==GSL_SUCCESS,"rk8pd reconstruction GSL failure");
        for(double value:y)Need(std::isfinite(value),"nonfinite rk8pd reconstructed state");
    }
    result.solve_wall_seconds=WallSeconds(begin);
    context->RequireFullCurrent();
    result.state=y;
    return result;
}

ReconstructionQualification Rk8pdCheckpointReconstructor::Qualify(
  const AcceptedStepRecord& step,const CheckpointState& O1,const CheckpointState& O2,
  const Rk8pdConfiguration& C1,const Rk8pdConfiguration& C2)
{
    ReconstructionQualification result;
    result.self_qualified=true;
    const Rkf45Configuration ultra;
    for(std::size_t i=0;i<3;++i)
    {
        const double M_O=std::max(std::abs(O1[i]),std::abs(O2[i]));
        result.d_O[i]=std::abs(O2[i]-O1[i]);
        result.D_O1[i]=C1.absolute[i]+C1.relative*M_O;
        const double D_O2=C2.absolute[i]+C2.relative*M_O;
        result.F_O[i]=std::max(D_O2,64.0*Ulp(M_O));
        const double M_i=std::max({std::abs(step.y_left[i]),std::abs(step.y_right[i]),
                                   std::abs(O1[i]),std::abs(O2[i])});
        const double D_U=ultra.absolute[i]+ultra.relative*M_i;
        result.F_i[i]=std::max(D_U,64.0*Ulp(M_i));
        result.U_O[i]=2.0*std::max(result.d_O[i],result.F_O[i]);
        if(!(result.d_O[i]<=result.D_O1[i]&&result.U_O[i]<=0.20*result.F_i[i]))
            result.self_qualified=false;
    }
    return result;
}

CheckpointOutput Rk8pdCheckpointReconstructor::Reconstruct(
  const ObservationBracket& bracket,const std::vector<AcceptedStepRecord>& history,
  const CheckpointContextFactory& context_factory) const
{
    Need(!history.empty(),"empty accepted-step history");
    CheckpointOutput output;
    output.observation_index=bracket.observation_index;
    output.t_observation_s=bracket.t_observation_s;
    output.provenance.O1=O1_;output.provenance.O2=O2_;
    output.provenance.left_accepted_ordinal=bracket.left_accepted_ordinal;
    output.provenance.right_accepted_ordinal=bracket.right_accepted_ordinal;
    output.provenance.t_left_s=bracket.t_left_s;output.provenance.t_right_s=bracket.t_right_s;
    output.provenance.t_observation_s=bracket.t_observation_s;
    if(bracket.initial)
    {
        Need(bracket.t_observation_s==history.front().t_left_s,"initial checkpoint/history mismatch");
        output.source=CheckpointSource::MainEndpoint;output.status=CheckpointStatus::Qualified;
        output.state=history.front().y_left;output.O1_state=output.state;output.O2_state=output.state;
        output.qualification.self_qualified=true;output.provenance.method="main-endpoint";
        return output;
    }
    Need(bracket.right_accepted_ordinal>0&&bracket.right_accepted_ordinal<=history.size(),
         "checkpoint references a missing accepted step");
    const auto& step=history[bracket.right_accepted_ordinal-1];
    Need(step.ordinal==bracket.right_accepted_ordinal&&step.t_left_s==bracket.t_left_s
      &&step.t_right_s==bracket.t_right_s,"checkpoint bracket/history mismatch");
    if(bracket.t_observation_s==step.t_left_s||bracket.t_observation_s==step.t_right_s)
    {
        output.source=CheckpointSource::MainEndpoint;output.status=CheckpointStatus::Qualified;
        output.state=bracket.t_observation_s==step.t_left_s?step.y_left:step.y_right;
        output.O1_state=output.state;output.O2_state=output.state;
        output.qualification.self_qualified=true;output.provenance.method="main-endpoint";
        return output;
    }
    Need(bracket.t_observation_s>step.t_left_s&&bracket.t_observation_s<step.t_right_s,
         "checkpoint is neither exact endpoint nor strict interior");
    output.source=CheckpointSource::Rk8pdReconstructed;
    output.provenance.method="rk8pd-two-level";
    const auto O1=Integrate(step.t_left_s,step.y_left,bracket.t_observation_s,O1_,context_factory);
    ++output.provenance.rk8pd_invocations;
    const auto O2=Integrate(step.t_left_s,step.y_left,bracket.t_observation_s,O2_,context_factory);
    ++output.provenance.rk8pd_invocations;
    output.O1_state=O1.state;output.O2_state=O2.state;
    output.performance.O1_context_wall_seconds=O1.context_wall_seconds;
    output.performance.O1_solve_wall_seconds=O1.solve_wall_seconds;
    output.performance.O2_context_wall_seconds=O2.context_wall_seconds;
    output.performance.O2_solve_wall_seconds=O2.solve_wall_seconds;
    output.provenance.context_identity_O1=O1.context_identity;
    output.provenance.context_identity_O2=O2.context_identity;
    output.qualification=Qualify(step,O1.state,O2.state,O1_,O2_);
    output.reconstruction_error_contribution=output.qualification.U_O;
    if(output.qualification.self_qualified)
    {
        output.status=CheckpointStatus::Qualified;
        output.state=O2.state;
    }
    else output.status=CheckpointStatus::NumericallyUnresolved;
    return output;
}

void Rk8pdCheckpointReconstructor::EvaluateDiagnostics(
  CheckpointOutput& output,const CheckpointContextFactory& context_factory) const
{
    Need(output.status==CheckpointStatus::Qualified,"diagnostics refuse unresolved checkpoint");
    Need(static_cast<bool>(context_factory),"missing diagnostic context factory");
    if(output.source==CheckpointSource::MainEndpoint)
    {
        auto begin=Clock::now();auto context=context_factory();
        output.performance.diagnostic_context_wall_seconds=WallSeconds(begin);
        Need(context!=nullptr,"diagnostic factory returned null");context->RequireFullCurrent();
        begin=Clock::now();output.diagnostics=context->EvaluateDiagnostics(output.t_observation_s,output.state);
        output.performance.diagnostic_evaluation_wall_seconds=WallSeconds(begin);
        context->RequireFullCurrent();
        output.O1_diagnostics=output.diagnostics;output.O2_diagnostics=output.diagnostics;
        output.diagnostics_evaluated=true;output.diagnostic_self_qualified=true;
        return;
    }
    auto begin=Clock::now();auto context1=context_factory();
    output.performance.diagnostic_context_wall_seconds=WallSeconds(begin);
    Need(context1!=nullptr,"O1 diagnostic factory returned null");context1->RequireFullCurrent();
    begin=Clock::now();output.O1_diagnostics=context1->EvaluateDiagnostics(output.t_observation_s,output.O1_state);
    output.performance.diagnostic_evaluation_wall_seconds=WallSeconds(begin);context1->RequireFullCurrent();
    begin=Clock::now();auto context2=context_factory();
    output.performance.diagnostic_context_wall_seconds+=WallSeconds(begin);
    Need(context2!=nullptr,"O2 diagnostic factory returned null");context2->RequireFullCurrent();
    begin=Clock::now();output.O2_diagnostics=context2->EvaluateDiagnostics(output.t_observation_s,output.O2_state);
    output.performance.diagnostic_evaluation_wall_seconds+=WallSeconds(begin);context2->RequireFullCurrent();
    output.diagnostics=output.O2_diagnostics;
    const auto values1=DiagnosticValues(output.O1_diagnostics),values2=DiagnosticValues(output.O2_diagnostics);
    const double G=DiagnosticScale(output.O1_diagnostics,output.O2_diagnostics);
    output.diagnostic_self_qualified=true;
    for(const auto& [name,value1]:values1)
    {
        const double value2=values2.at(name),magnitude=std::max(std::abs(value1),std::abs(value2));
        const double d=std::abs(value2-value1);
        const double floor=std::max(1.0e-11*std::max(1.0,magnitude),64.0*Ulp(magnitude));
        const double F=name=="P_dir_eq"||name=="P_dir_actual"||name=="LH"||name=="DeltaLnu"
          ||name=="DeltaPbeta"||name=="Lnu_eq"||name=="Lnu_full"||name=="Lgamma"
          ||name=="Lother"||name=="Pnet"?std::max(1.0e-11*G,floor):floor;
        const double U=2.0*std::max(d,64.0*Ulp(magnitude));
        output.diagnostic_d_O[name]=d;output.diagnostic_F_P[name]=F;output.diagnostic_U_O[name]=U;
        if(!(U<=0.20*F))output.diagnostic_self_qualified=false;
    }
    if(!output.diagnostic_self_qualified)
    {
        output.status=CheckpointStatus::NumericallyUnresolved;
        output.diagnostics_evaluated=false;
        return;
    }
    output.diagnostics_evaluated=true;
}

CheckpointR20Result ComputeCheckpointR20(const std::vector<CheckpointOutput>& checkpoints)
{
    Need(checkpoints.size()>=2,"R20 requires at least two checkpoints");
    for(const auto& checkpoint:checkpoints)
        Need(checkpoint.status==CheckpointStatus::Qualified&&checkpoint.diagnostics_evaluated
          &&checkpoint.diagnostic_self_qualified,"R20 refuses unresolved checkpoint diagnostics");
    const auto O2=R20From(checkpoints,false),O1=R20From(checkpoints,true);
    CheckpointR20Result result;
    result.delta_Eeq_erg=O2.delta_eq;result.delta_Echem_erg=O2.delta_chem;
    result.delta_Uth_erg=O2.delta_uth;result.outgoing_integral_erg=O2.outgoing;
    result.R20_residual_erg=O2.residual;result.N_R20_erg=O2.normalizer;
    result.R20_normalized=O2.residual/O2.normalizer;
    const double magnitude=std::max(std::abs(O1.residual),std::abs(O2.residual));
    result.reconstruction_uncertainty_erg=2.0*std::max(std::abs(O2.residual-O1.residual),64.0*Ulp(magnitude));
    return result;
}

} // namespace CompactStar::Physics::BNV
