#include <CompactStar/Physics/BNV/PassiveCheckpointOutput.hpp>

#include <cmath>
#include <cstring>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace BNV=CompactStar::Physics::BNV;
namespace RC=CompactStar::Physics::Rotochemical;

namespace
{
void Require(bool condition,const char* message)
{
    if(!condition)throw std::runtime_error(message);
}

template<class Function>
bool Throws(Function&& function)
{
    try{function();return false;}catch(const std::exception&){return true;}
}

struct FactoryAudit
{
    std::size_t constructions=0;
    std::vector<std::size_t> live_ids;
    BNV::CheckpointState derivative{{0,0,0}};
    bool varying_derivative=false;
    BNV::CheckpointState last_diagnostic_state{};
    std::size_t diagnostic_calls=0;
};

class SyntheticContext final:public BNV::ICheckpointEvaluationContext
{
  public:
    SyntheticContext(std::shared_ptr<FactoryAudit> audit,std::size_t id)
      :audit_(std::move(audit)),id_(id){audit_->live_ids.push_back(id_);}
    void RequireFullCurrent() const override{}
    void Derivative(double,const BNV::CheckpointState&,BNV::CheckpointState& out) override
    {
        out=audit_->derivative;
        if(audit_->varying_derivative)for(double& value:out)value*=static_cast<double>(id_-1);
    }
    BNV::BnvDiagnostics EvaluateDiagnostics(double t,const BNV::CheckpointState& state) override
    {
        ++audit_->diagnostic_calls;audit_->last_diagnostic_state=state;
        BNV::BnvDiagnostics d;
        d.t_s=t;d.B_count=10.0-t;d.mu_B_inf_MeV=2;d.mu_n_actual_inf_MeV=3+state[1];
        d.sigma_count_s={{state[1],state[2]}};d.Echem_MeV=1+state[1]*state[1]+state[2]*state[2];
        d.P_dir_eq_erg_s=10+state[0];d.P_dir_actual_erg_s=11+state[0];
        d.LH_erg_s=2+state[0];d.DeltaLnu_erg_s=3+state[0];d.DeltaPbeta_erg_s=4+state[0];
        d.Lnu_eq_erg_s=5+state[0];d.Lnu_full_erg_s=8+2*state[0];
        d.Lgamma_erg_s=6+state[0];d.Lother_erg_s=7+state[0];d.Pnet_erg_s=state[0];
        d.L_out_fluid_inf_erg_s=1;d.Cstar_erg_K=4;d.Tinf_K=5+state[0];
        d.valid_through_sample=true;return d;
    }
    std::string Identity() const override{return "synthetic-isolated-"+std::to_string(id_);}
  private:
    std::shared_ptr<FactoryAudit> audit_;
    std::size_t id_;
};

BNV::CheckpointContextFactory Factory(const std::shared_ptr<FactoryAudit>& audit)
{
    return [audit]
    {
        const auto id=++audit->constructions;
        return std::make_unique<SyntheticContext>(audit,id);
    };
}

BNV::AcceptedStepRecord Step()
{
    BNV::AcceptedStepRecord step;
    step.ordinal=1;step.t_left_s=0;step.t_right_s=1;step.accepted_step_s=1;
    step.y_left={{0,0,0}};step.y_right={{1,0,0}};step.suggested_next_h_s=2;
    step.source_identity="synthetic-source";step.currentness_identity="synthetic-current";
    return step;
}

BNV::ObservationBracket Interior()
{
    return {1,0.5,0,1,0,1,false};
}

BNV::ObservationBracket Endpoint()
{
    return {2,1.0,0,1,0,1,false};
}

void P1PassiveAuthority()
{
    auto step=Step();const auto before=step;
    BNV::PassiveObservationSchedule schedule({0,0.25,0.5,1},"schedule-a");
    schedule.Replay(step.y_left,{step});
    Require(schedule.Brackets().size()==4&&schedule.AcceptedCallbacks()==1,"P1 bracket count");
    Require(std::memcmp(&before.y_left,&step.y_left,sizeof(step.y_left))==0
      &&std::memcmp(&before.y_right,&step.y_right,sizeof(step.y_right))==0,"P1 mutated history");
}

void P2ExactEndpointBypass()
{
    auto audit=std::make_shared<FactoryAudit>();
    BNV::Rk8pdCheckpointReconstructor reconstruction;
    const auto result=reconstruction.Reconstruct(Endpoint(),{Step()},Factory(audit));
    Require(result.source==BNV::CheckpointSource::MainEndpoint,"P2 source");
    Require(result.status==BNV::CheckpointStatus::Qualified,"P2 status");
    Require(result.state==Step().y_right,"P2 exact main state");
    Require(result.provenance.rk8pd_invocations==0&&audit->constructions==0,"P2 invoked reconstruction");
}

void P3StrictInteriorExactlyTwoLevels()
{
    auto audit=std::make_shared<FactoryAudit>();audit->derivative={{1,0,0}};
    BNV::Rk8pdCheckpointReconstructor reconstruction;
    const auto result=reconstruction.Reconstruct(Interior(),{Step()},Factory(audit));
    Require(result.provenance.rk8pd_invocations==2&&audit->constructions==2,"P3 not exactly O1/O2");
    Require(result.provenance.context_identity_O1!=result.provenance.context_identity_O2,"P3 contexts shared");
}

void P4SelfQualificationPass()
{
    auto audit=std::make_shared<FactoryAudit>();audit->derivative={{1,0,0}};
    BNV::Rk8pdCheckpointReconstructor reconstruction;
    const auto result=reconstruction.Reconstruct(Interior(),{Step()},Factory(audit));
    Require(result.qualification.self_qualified&&result.status==BNV::CheckpointStatus::Qualified,"P4 pass refused");
    Require(result.state==result.O2_state,"P4 did not select O2");
}

void P5FailClosedSynthetic()
{
    auto audit=std::make_shared<FactoryAudit>();audit->derivative={{1.0e-4,1.0e-4,1.0e-4}};
    audit->varying_derivative=true;
    BNV::Rk8pdCheckpointReconstructor reconstruction;
    auto result=reconstruction.Reconstruct(Interior(),{Step()},Factory(audit));
    Require(!result.qualification.self_qualified,"P5 synthetic failure passed");
    Require(result.status==BNV::CheckpointStatus::NumericallyUnresolved,"P5 did not fail closed");
    Require(Throws([&]{reconstruction.EvaluateDiagnostics(result,Factory(audit));}),"P5 downstream did not refuse");
}

void P6NoFallbackOrRetuning()
{
    auto audit=std::make_shared<FactoryAudit>();audit->derivative={{1.0e-4,1.0e-4,1.0e-4}};
    audit->varying_derivative=true;
    BNV::Rk8pdCheckpointReconstructor reconstruction;
    const auto O1=reconstruction.O1(),O2=reconstruction.O2();
    const auto result=reconstruction.Reconstruct(Interior(),{Step()},Factory(audit));
    Require(result.provenance.rk8pd_invocations==2&&result.provenance.fallback_invocations==0,"P6 fallback/retry");
    Require(reconstruction.O1().relative==O1.relative&&reconstruction.O1().absolute==O1.absolute
      &&reconstruction.O2().relative==O2.relative&&reconstruction.O2().absolute==O2.absolute,"P6 tolerance change");
}

void P7ContextIsolation()
{
    auto audit=std::make_shared<FactoryAudit>();audit->derivative={{1,0,0}};
    BNV::Rk8pdCheckpointReconstructor reconstruction;
    const auto result=reconstruction.Reconstruct(Interior(),{Step()},Factory(audit));
    Require(audit->live_ids.size()==2&&audit->live_ids[0]!=audit->live_ids[1],"P7 context reuse");
    Require(result.O1_state==result.O2_state,"P7 isolated deterministic contexts disagree");
}

void P8DiagnosticsUseReconstructedState()
{
    auto audit=std::make_shared<FactoryAudit>();audit->derivative={{1,0,0}};
    BNV::Rk8pdCheckpointReconstructor reconstruction;
    auto result=reconstruction.Reconstruct(Interior(),{Step()},Factory(audit));
    reconstruction.EvaluateDiagnostics(result,Factory(audit));
    Require(result.diagnostics_evaluated&&result.diagnostic_self_qualified,"P8 diagnostic qualification");
    Require(audit->last_diagnostic_state==result.O2_state,"P8 diagnostics used another state");
    Require(result.diagnostics.Pnet_erg_s==result.O2_state[0],"P8 diagnostic packet mismatch");
}

BNV::CheckpointOutput R20Point(std::size_t index,double t,double B,double T,double Echem)
{
    BNV::CheckpointOutput point;
    point.observation_index=index;point.t_observation_s=t;
    point.status=BNV::CheckpointStatus::Qualified;point.qualification.self_qualified=true;
    point.diagnostics_evaluated=true;point.diagnostic_self_qualified=true;
    BNV::BnvDiagnostics d;
    d.t_s=t;d.B_count=B;d.mu_B_inf_MeV=2;d.Echem_MeV=Echem;
    d.Cstar_erg_K=4;d.Tinf_K=T;d.L_out_fluid_inf_erg_s=1;d.valid_through_sample=true;
    point.diagnostics=d;point.O1_diagnostics=d;point.O2_diagnostics=d;
    return point;
}

void P9R20CheckpointDiagnostics()
{
    const auto result=BNV::ComputeCheckpointR20({R20Point(0,0,10,5,1),R20Point(1,2,9,7,3)});
    Require(result.delta_Uth_erg==8&&result.outgoing_integral_erg==2,"P9 trapezoid changed");
    const double expected_eq=-2*RC::MeVToErg,expected_chem=2*RC::MeVToErg;
    Require(result.delta_Eeq_erg==expected_eq&&result.delta_Echem_erg==expected_chem,"P9 ledger changed");
    Require(result.R20_residual_erg==10&&result.N_R20_erg==8&&result.R20_normalized==1.25,"P9 R20 semantics");
    Require(result.reconstruction_uncertainty_erg>0,"P9 uncertainty not separate");
}

void P10ScheduleRefusalsAndReplay()
{
    Require(Throws([]{BNV::PassiveObservationSchedule duplicate({0,0.5,0.5,1},"duplicate");}),
            "P10 duplicate accepted");
    BNV::PassiveObservationSchedule missing({0,0.5,1},"missing");
    missing.NotifyInitial({{0,0,0}});
    Require(Throws([&]{missing.RequireComplete();}),"P10 missing observations accepted");
    const auto history=std::vector<BNV::AcceptedStepRecord>{Step()};
    BNV::PassiveObservationSchedule first({0,0.25,0.75,1},"first");
    BNV::PassiveObservationSchedule second({0,0.1,0.9,1},"second");
    first.Replay(history.front().y_left,history);second.Replay(history.front().y_left,history);
    Require(first.Brackets().size()==second.Brackets().size()&&history.front().ordinal==1,
            "P10 postprocessing schedule changed history");
}
}

int main()
{
    try
    {
        P1PassiveAuthority();std::cout<<"P1 PASS\n";
        P2ExactEndpointBypass();std::cout<<"P2 PASS\n";
        P3StrictInteriorExactlyTwoLevels();std::cout<<"P3 PASS\n";
        P4SelfQualificationPass();std::cout<<"P4 PASS\n";
        P5FailClosedSynthetic();std::cout<<"P5 PASS\n";
        P6NoFallbackOrRetuning();std::cout<<"P6 PASS\n";
        P7ContextIsolation();std::cout<<"P7 PASS\n";
        P8DiagnosticsUseReconstructedState();std::cout<<"P8 PASS\n";
        P9R20CheckpointDiagnostics();std::cout<<"P9 PASS\n";
        P10ScheduleRefusalsAndReplay();std::cout<<"P10 PASS\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<"FAIL: "<<error.what()<<'\n';
        return 1;
    }
}
