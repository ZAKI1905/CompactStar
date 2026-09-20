#pragma once

#include <CompactStar/Physics/BNV/BnvDiagnostics.hpp>
#include <CompactStar/Physics/BNV/DirectEnergyLedger.hpp>
#include <CompactStar/Physics/BNV/MovingReferenceSource.hpp>
#include <CompactStar/Physics/BNV/StaticZeroSpinHistory.hpp>
#include <CompactStar/Physics/Rotochemical/FrozenRotochemicalRunContext.hpp>
#include <optional>

namespace CompactStar::Physics::BNV
{

struct ControlledBnvEvaluation
{
    Rotochemical::SecularEvaluation ordinary;
    MovingReferenceSample moving;
    DirectEnergyResult direct;
    BnvDiagnostics diagnostics;
    std::array<double,2> eta_dot_MeV_s{};
    double x_dot_s=0;
};

class FrozenControlledBnvRunContext final
{
  public:
    FrozenControlledBnvRunContext(
      std::shared_ptr<const Rotochemical::FrozenRotochemicalRunContext> ordinary,
      std::shared_ptr<const EquilibriumBaryonTangent> tangent,
      std::shared_ptr<const OrdinaryMatterBnvHistory> history,
      std::shared_ptr<const IDirectBnvEnergyPartition> partition,
      ProductFateLedger fate,
      std::shared_ptr<const FrozenBnvValidityMonitor> validity,
      std::shared_ptr<const Rotochemical::GlobalUrcaChannelCoefficient> channels,
      double mu_B_inf_MeV,
      std::string actual_potential_provenance,
      std::string run_card_identity);

    void RequireFullCurrent() const;
    void RequireCheapCurrent() const;
    Evolution::DriverContext DriverContext() const {RequireCheapCurrent();return ordinary_->DriverContext();}
    ControlledBnvEvaluation Evaluate(double,const Evolution::StateVector&,const Evolution::DriverContext&) const;
    ControlledBnvEvaluation EvaluateReactionFree(double,const Evolution::StateVector&,const Evolution::DriverContext&) const;
    const std::shared_ptr<const Rotochemical::FrozenRotochemicalRunContext>& OrdinaryContext() const{return ordinary_;}
    const std::shared_ptr<const Rotochemical::PrescribedSpinHistory>& SpinOwner() const{return spin_;}
    const std::shared_ptr<const EquilibriumBaryonTangent>& Tangent() const{return tangent_;}
    const std::string& RunCardIdentity() const{return run_card_identity_;}
  private:
    static void Need(bool,const char*);
    ControlledBnvEvaluation EvaluateImpl(double,const Evolution::StateVector&,const Evolution::DriverContext&,bool) const;
    std::shared_ptr<const Rotochemical::FrozenRotochemicalRunContext> ordinary_;
    std::shared_ptr<const EquilibriumBaryonTangent> tangent_;
    std::shared_ptr<const OrdinaryMatterBnvHistory> history_;
    std::shared_ptr<const IDirectBnvEnergyPartition> partition_;
    ProductFateLedger fate_;
    std::shared_ptr<const FrozenBnvValidityMonitor> validity_;
    std::shared_ptr<const Rotochemical::GlobalUrcaChannelCoefficient> channels_;
    std::shared_ptr<const Rotochemical::PrescribedSpinHistory> spin_;
    double mu_B_inf_MeV_=0;
    bool spin_on_zero_source_regression_=false;
    std::string potential_provenance_,run_card_identity_,history_identity_,partition_identity_,spin_identity_;
    mutable std::optional<Rotochemical::SecularEvaluation> reaction_free_reference_;
    mutable double reaction_free_reference_Tinf_K_=0;
};

} // namespace CompactStar::Physics::BNV
