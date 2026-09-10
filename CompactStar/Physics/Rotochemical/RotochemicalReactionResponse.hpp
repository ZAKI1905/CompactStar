#pragma once
#include <CompactStar/Physics/Rotochemical/GlobalUrcaChannelCoefficient.hpp>
#include <CompactStar/Physics/Rotochemical/UrcaImbalanceFunctions.hpp>

namespace CompactStar::Physics::Rotochemical
{
struct RotochemicalReactionResult
{
    std::array<double,2> rate_count_s{};
    std::array<double,4> equilibrium_erg_s{},increment_erg_s{},full_erg_s{};
    double chemical_power_MeV_s=0;
    double Rate(BetaChannel c) const { return rate_count_s[ChannelIndex(c)]; }
    double EquilibriumErgPerSecond() const { double x=0;for(double v:equilibrium_erg_s)x+=v;return x; }
    double IncrementErgPerSecond() const { double x=0;for(double v:increment_erg_s)x+=v;return x; }
};
// Layer 4: immutable response input, no Z/W/spin or thermal denominator.
class RotochemicalReactionResponse
{
  public:
    explicit RotochemicalReactionResponse(std::shared_ptr<const GlobalUrcaChannelCoefficient> c):coefficients_(std::move(c))
    { if(!coefficients_) throw std::runtime_error("Urca coefficients required");coefficients_->RequireCurrent(); }
    const std::shared_ptr<const GlobalUrcaChannelCoefficient>& Coefficients() const { coefficients_->RequireCurrent();return coefficients_; }
    RotochemicalReactionResult Evaluate(double Tinf_K,const ChemicalImbalanceState& eta) const
    {
        coefficients_->RequireCurrent();
        RotochemicalReactionResult result;
        for(auto p:UrcaProcesses) {
            auto l=LeptonChannel(p);const double xi=eta.Xi(l,Tinf_K);
            const auto i=ProcessIndex(p);const auto q=TemperatureExponent(p);
            const double L=coefficients_->LuminosityCoefficient(p);
            if(L==0) continue;
            const double increment=q==6?UrcaImbalanceFunctions::DirectIncrement(xi):UrcaImbalanceFunctions::ModifiedIncrement(xi);
            const double H=q==6?UrcaImbalanceFunctions::HD(xi):UrcaImbalanceFunctions::HM(xi);
            const double equilibrium=L*std::pow(Tinf_K,q);
            const double rate=L/BoltzmannErgPerK*std::pow(Tinf_K,q-1)*H;
            result.equilibrium_erg_s[i]=equilibrium;
            result.increment_erg_s[i]=equilibrium*increment;
            result.full_erg_s[i]=equilibrium+result.increment_erg_s[i];
            result.rate_count_s[ChannelIndex(l)]+=rate;
            result.chemical_power_MeV_s+=eta.InfinityMeV(l)*rate;
            if(!std::isfinite(result.full_erg_s[i])||!std::isfinite(rate)||!std::isfinite(result.chemical_power_MeV_s)) throw std::runtime_error("nonfinite rotochemical response");
        }
        return result;
    }
  private:
    std::shared_ptr<const GlobalUrcaChannelCoefficient> coefficients_;
};
// Thermal luminosity boundary: the sole public MeV/s -> erg/s conversion.
struct RotochemicalThermalPower
{
    double heating_erg_s, neutrino_increment_erg_s, incremental_beta_erg_s, full_beta_erg_s;
    static RotochemicalThermalPower From(const RotochemicalReactionResult& r)
    { double h=MeVToErg*r.chemical_power_MeV_s, d=r.IncrementErgPerSecond(); return {h,d,h-d,h-d-r.EquilibriumErgPerSecond()}; }
};
}
