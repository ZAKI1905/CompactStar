#pragma once
#include <CompactStar/Analysis/ChemicalResponse.hpp>
#include <CompactStar/Physics/State/ChemState.hpp>
#include <CompactStar/Units.hpp>
#include <Zaki/Physics/Constants.hpp>
#include <array>
#include <cmath>
#include <stdexcept>

namespace CompactStar::Physics::Rotochemical
{
using BetaChannel = Analysis::ImbalanceChannel;
inline const double BoltzmannMeVPerK = Zaki::Physics::K_BOLTZ_EV * 1e-6;
inline constexpr double MeVToErg = Units::MEV_FM3_TO_ERG_CM3 / 1e39;
inline const double BoltzmannErgPerK = BoltzmannMeVPerK * MeVToErg;
inline std::size_t ChannelIndex(BetaChannel c)
{
    switch(c) { case BetaChannel::Npe: return 0; case BetaChannel::NpMu: return 1; }
    throw std::runtime_error("invalid beta channel");
}
// Immutable typed value. Generic ChemState supplies storage only.
class ChemicalImbalanceState
{
  public:
    ChemicalImbalanceState(double eta_npe_infinity_MeV, double eta_npmu_infinity_MeV)
      : eta_{eta_npe_infinity_MeV,eta_npmu_infinity_MeV}
    { for(double v:eta_) if(!std::isfinite(v)) throw std::runtime_error("nonfinite eta"); }
    static ChemicalImbalanceState Read(const State::ChemState& s)
    { if(s.Size()!=2) throw std::runtime_error("rotochemical state requires Npe,NpMu"); return {s.Eta(0),s.Eta(1)}; }
    void Store(State::ChemState& s) const
    { if(s.Size()!=2) throw std::runtime_error("rotochemical state requires Npe,NpMu"); s.Eta(0)=eta_[0];s.Eta(1)=eta_[1]; }
    double InfinityMeV(BetaChannel c) const { return eta_[ChannelIndex(c)]; }
    double LocalMeV(BetaChannel c,double nu) const { return InfinityMeV(c)*std::exp(-nu); }
    double Xi(BetaChannel c,double Tinf_K) const
    { if(!(Tinf_K>0) || !std::isfinite(Tinf_K)) throw std::runtime_error("invalid Tinf K"); return InfinityMeV(c)/(BoltzmannMeVPerK*Tinf_K); }
  private:
    std::array<double,2> eta_;
};
}
