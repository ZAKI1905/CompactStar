#pragma once

#include <CompactStar/Physics/BNV/BnvDiagnostics.hpp>
#include <CompactStar/Physics/Evolution/EvolutionSystem.hpp>
#include <CompactStar/Physics/Evolution/StateLayout.hpp>
#include <CompactStar/Physics/Rotochemical/FrozenRotochemicalRunContext.hpp>

#include <array>
#include <cstddef>
#include <exception>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace CompactStar::Physics::BNV
{

using CheckpointState=std::array<double,3>;

struct Rkf45Configuration
{
    double relative=1.0e-11;
    CheckpointState absolute{{1.0e-16,1.0e-22,1.0e-22}};
    double initial_step_s=1.0;
    void Validate() const;
};

struct Rk8pdConfiguration
{
    std::string level;
    double relative=0;
    CheckpointState absolute{};
    void Validate() const;
};

struct MainIntegrationProvenance
{
    std::string run_card_identity;
    std::string source_identity;
    std::string currentness_identity;
    std::string observation_schedule_identity;
};

struct AcceptedStepRecord
{
    std::size_t ordinal=0;
    double t_left_s=0;
    CheckpointState y_left{};
    double t_right_s=0;
    CheckpointState y_right{};
    double accepted_step_s=0;
    double suggested_next_h_s=0;
    std::size_t cumulative_rejected=0;
    std::string source_identity;
    std::string currentness_identity;
    int cstar_cell_left=0;
    int cstar_cell_right=0;
    int cstar_knots_crossed=0;
};

struct ObservationBracket
{
    std::size_t observation_index=0;
    double t_observation_s=0;
    std::size_t left_accepted_ordinal=0;
    std::size_t right_accepted_ordinal=0;
    double t_left_s=0;
    double t_right_s=0;
    bool initial=false;
};

class PassiveObservationSchedule final
{
  public:
    PassiveObservationSchedule(std::vector<double> requested_times_s,std::string identity);
    void Replay(const CheckpointState& initial_state,const std::vector<AcceptedStepRecord>& history);
    void NotifyInitial(const CheckpointState& initial_state);
    void NotifyAccepted(const AcceptedStepRecord& step);
    void RequireComplete() const;
    const std::vector<double>& RequestedTimes() const{return requested_times_s_;}
    const std::vector<ObservationBracket>& Brackets() const{return brackets_;}
    const std::string& Identity() const{return identity_;}
    std::size_t AcceptedCallbacks() const{return accepted_callbacks_;}
  private:
    void Reset();
    std::vector<double> requested_times_s_;
    std::string identity_;
    std::vector<ObservationBracket> brackets_;
    std::size_t next_=0;
    std::size_t accepted_callbacks_=0;
};

struct MainIntegrationStatistics
{
    std::size_t accepted_steps=0;
    std::size_t rejected_steps=0;
    std::size_t rhs_evaluations=0;
    std::size_t gsl_apply_calls=0;
    std::size_t distinct_positive_t1_targets=0;
    double positive_t1_s=0;
    double minimum_step_s=0;
    double maximum_step_s=0;
    double wall_seconds=0;
    double cpu_seconds=0;
};

struct MainTrajectoryResult
{
    double start_time_s=0;
    double final_time_s=0;
    CheckpointState initial_state{};
    CheckpointState final_state{};
    Rkf45Configuration configuration;
    MainIntegrationProvenance provenance;
    MainIntegrationStatistics statistics;
    std::vector<AcceptedStepRecord> accepted_steps;
};

class UninterruptedBnvTrajectory final
{
  public:
    UninterruptedBnvTrajectory(
      Evolution::EvolutionSystem& system,
      const Evolution::StateLayout& layout,
      std::shared_ptr<const Rotochemical::FrozenRotochemicalRunContext> currentness_owner,
      Rkf45Configuration configuration,
      MainIntegrationProvenance provenance);
    MainTrajectoryResult Integrate(
      double start_time_s,double final_time_s,CheckpointState initial_state,
      PassiveObservationSchedule& observation_schedule,bool enforce_thermal_domain=true);
  private:
    static int Callback(double,const double*,double*,void*) noexcept;
    void Derivative(double,const double*,double*);
    void ValidateAccepted(const double*,bool) const;
    Evolution::EvolutionSystem& system_;
    const Evolution::StateLayout& layout_;
    std::shared_ptr<const Rotochemical::FrozenRotochemicalRunContext> currentness_owner_;
    Rkf45Configuration configuration_;
    MainIntegrationProvenance provenance_;
    MainIntegrationStatistics* statistics_=nullptr;
    std::exception_ptr failure_;
};

class ICheckpointEvaluationContext
{
  public:
    virtual ~ICheckpointEvaluationContext()=default;
    virtual void RequireFullCurrent() const=0;
    virtual void Derivative(double,const CheckpointState&,CheckpointState&)=0;
    virtual BnvDiagnostics EvaluateDiagnostics(double,const CheckpointState&)=0;
    virtual std::string Identity() const=0;
};

using CheckpointContextFactory=std::function<std::unique_ptr<ICheckpointEvaluationContext>()>;

enum class CheckpointSource { MainEndpoint,Rk8pdReconstructed };
enum class CheckpointStatus { Qualified,NumericallyUnresolved };

struct ReconstructionQualification
{
    CheckpointState d_O{};
    CheckpointState D_O1{};
    CheckpointState F_O{};
    CheckpointState F_i{};
    CheckpointState U_O{};
    bool self_qualified=false;
};

struct CheckpointProvenance
{
    std::string method;
    std::size_t left_accepted_ordinal=0;
    std::size_t right_accepted_ordinal=0;
    double t_left_s=0;
    double t_right_s=0;
    double t_observation_s=0;
    Rk8pdConfiguration O1;
    Rk8pdConfiguration O2;
    std::size_t rk8pd_invocations=0;
    std::size_t fallback_invocations=0;
    std::string context_identity_O1;
    std::string context_identity_O2;
};

struct CheckpointPerformance
{
    double O1_context_wall_seconds=0;
    double O1_context_cpu_seconds=0;
    double O1_solve_wall_seconds=0;
    double O1_solve_cpu_seconds=0;
    double O2_context_wall_seconds=0;
    double O2_context_cpu_seconds=0;
    double O2_solve_wall_seconds=0;
    double O2_solve_cpu_seconds=0;
    double diagnostic_context_wall_seconds=0;
    double diagnostic_context_cpu_seconds=0;
    double diagnostic_evaluation_wall_seconds=0;
    double diagnostic_evaluation_cpu_seconds=0;
};

struct CheckpointOutput
{
    std::size_t observation_index=0;
    double t_observation_s=0;
    CheckpointSource source=CheckpointSource::MainEndpoint;
    CheckpointStatus status=CheckpointStatus::NumericallyUnresolved;
    CheckpointState state{};
    CheckpointState O1_state{};
    CheckpointState O2_state{};
    ReconstructionQualification qualification;
    CheckpointProvenance provenance;
    CheckpointPerformance performance;
    BnvDiagnostics diagnostics;
    BnvDiagnostics O1_diagnostics;
    BnvDiagnostics O2_diagnostics;
    bool diagnostics_evaluated=false;
    bool diagnostic_self_qualified=false;
    std::map<std::string,double> diagnostic_d_O;
    std::map<std::string,double> diagnostic_U_O;
    std::map<std::string,double> diagnostic_F_P;
    CheckpointState reconstruction_error_contribution{};
};

class Rk8pdCheckpointReconstructor final
{
  public:
    Rk8pdCheckpointReconstructor();
    CheckpointOutput Reconstruct(
      const ObservationBracket& bracket,const std::vector<AcceptedStepRecord>& history,
      const CheckpointContextFactory& context_factory) const;
    void EvaluateDiagnostics(CheckpointOutput&,const CheckpointContextFactory&) const;
    const Rk8pdConfiguration& O1() const{return O1_;}
    const Rk8pdConfiguration& O2() const{return O2_;}
  private:
    struct SolveResult
    {
        CheckpointState state{};
        std::string context_identity;
        double context_wall_seconds=0;
        double context_cpu_seconds=0;
        double solve_wall_seconds=0;
        double solve_cpu_seconds=0;
    };
    static SolveResult Integrate(
      double,const CheckpointState&,double,const Rk8pdConfiguration&,
      const CheckpointContextFactory&);
    static ReconstructionQualification Qualify(
      const AcceptedStepRecord&,const CheckpointState&,const CheckpointState&,
      const Rk8pdConfiguration&,const Rk8pdConfiguration&);
    Rk8pdConfiguration O1_;
    Rk8pdConfiguration O2_;
};

struct CheckpointR20Result
{
    double delta_Eeq_erg=0;
    double delta_Echem_erg=0;
    double delta_Uth_erg=0;
    double outgoing_integral_erg=0;
    double R20_residual_erg=0;
    double N_R20_erg=0;
    double R20_normalized=0;
    double reconstruction_uncertainty_erg=0;
};

CheckpointR20Result ComputeCheckpointR20(const std::vector<CheckpointOutput>& checkpoints);

} // namespace CompactStar::Physics::BNV
