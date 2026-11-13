#pragma once
#include <CompactStar/Physics/BNV/DirectEnergyLedger.hpp>
#include <iomanip>
#include <sstream>

namespace CompactStar::Physics::BNV
{
// Controlled P-A mathematical partition. No microscopic rate/kinematics claim.
class ActualPotentialPartition final : public IDirectBnvEnergyPartition
{
  public:
    ActualPotentialPartition(double fraction,std::string fate,
                             std::shared_ptr<const BnvDependencyToken> token)
      :fraction_(fraction),fate_(std::move(fate)),token_(std::move(token)),generation_(token_?token_->generation:0)
    {
        if(!std::isfinite(fraction_)||fraction_<0||fraction_>1||fate_.empty())
            throw std::runtime_error("invalid P-A retained fraction");
        std::ostringstream s;s<<std::setprecision(17)<<"ADR-0025 P-A actual chemical budget fraction="<<fraction_
          <<"; all fluid exit promptly escapes; cold mathematical partition";
        identity_=s.str();RequireCurrent();
    }
    void RequireCurrent() const override
    {if(!token_||!token_->alive||token_->generation!=generation_)throw std::runtime_error("stale actual-potential partition");}
    const std::string& Identity() const override {RequireCurrent();return identity_;}
    const std::string& FiniteTemperatureWeightingClass() const override
    {
        static const std::string name="mathematical cold P-A; occupied-state weighting unspecified; diagnostic O(kBT) scale, not a physical bound";
        return name;
    }
    double OmittedFiniteTemperaturePowerErgPerSecond(const OrdinaryMatterBnvHistorySample& s,double T) const override
    {
        RequireCurrent();if(!(T>0)||!std::isfinite(T))throw std::runtime_error("invalid partition temperature");
        double count=0;for(const auto& e:s.events)count+=e.rate_count_s;
        return count*Rotochemical::BoltzmannErgPerK*T;
    }
    std::vector<DirectEventEnergy> Evaluate(const OrdinaryMatterBnvHistorySample& s,const ActualOrdinaryPotential& p) const override
    {
        RequireCurrent();s.Validate();std::vector<DirectEventEnergy> out;
        for(const auto& event:s.events)
        {
            if(event.delta_y!=std::array<double,3>{-1,0,0})throw std::runtime_error("controlled P-A requires one neutron sink per event");
            const double escape=(1-fraction_)*p.mu_n_actual_inf_MeV;
            DirectEventEnergy e{event.event_id,identity_,fate_,escape,escape,0,0,false,1,0};
            e.actual_budget_retained_fraction=fraction_;out.push_back(e);
        }
        return out;
    }
  private:
    const double fraction_;
    const std::string fate_;
    const std::shared_ptr<const BnvDependencyToken> token_;
    const std::uint64_t generation_;
    std::string identity_;
};
}
