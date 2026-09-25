#include <CompactStar/Physics/BNV/ControlledBnvSecularDriver.hpp>
#include <cmath>
#include <stdexcept>

namespace CompactStar::Physics::BNV
{
ControlledBnvSecularDriver::ControlledBnvSecularDriver(std::shared_ptr<const FrozenControlledBnvRunContext> c,Mode mode)
  :context_(std::move(c)),mode_(mode){if(!context_)throw std::runtime_error("missing controlled BNV context");context_->RequireCheapCurrent();}
const std::vector<State::StateTag>& ControlledBnvSecularDriver::DependsOn() const
{static const std::vector<State::StateTag> tags{State::StateTag::Thermal,State::StateTag::Chem};return tags;}
ControlledBnvEvaluation ControlledBnvSecularDriver::Evaluate(double t,const Evolution::StateVector& s,const Evolution::DriverContext& c) const
{
    auto out=mode_==Mode::ReactionFreeControl?context_->EvaluateReactionFree(t,s,c):context_->Evaluate(t,s,c);
    if(mode_==Mode::ReactionFreeControl)
    {
        out.ordinary.reaction={};out.ordinary.beta={};out.ordinary.eta_dot_MeV_s={};
        out.eta_dot_MeV_s=out.diagnostics.eta_dot_from_sigma_MeV_s;out.x_dot_s=0;
        auto& d=out.diagnostics;
        d.R_count_s={};d.eta_dot_from_beta_MeV_s={};
        d.Echem_dot_reaction_MeV_s=0;d.Echem_dot_total_MeV_s=d.Echem_dot_source_MeV_s;
        d.LH_erg_s=0;d.DeltaLnu_erg_s=0;d.DeltaPbeta_erg_s=0;
        d.Lnu_eq_erg_s=0;d.Lnu_full_erg_s=0;
        d.Pnet_erg_s=d.P_dir_actual_erg_s;
    }
    return out;
}
void ControlledBnvSecularDriver::AccumulateRHS(double t,const Evolution::StateVector& state,Evolution::RHSAccumulator& rhs,const Evolution::DriverContext& ctx) const
{
    context_->RequireCheapCurrent();
    if(rhs.Block(State::StateTag::Thermal).size()!=1||rhs.Block(State::StateTag::Chem).size()!=2)
        throw std::runtime_error("malformed controlled BNV RHS blocks");
    const auto out=Evaluate(t,state,ctx);context_->RequireCheapCurrent();
    for(std::size_t i=0;i<2;++i)if(!std::isfinite(rhs.Block(State::StateTag::Chem)[i]+out.eta_dot_MeV_s[i]))
        throw std::runtime_error("invalid controlled BNV chemical accumulator");
    if(mode_==Mode::Coupled&&!std::isfinite(rhs.Block(State::StateTag::Thermal)[0]+out.x_dot_s))
        throw std::runtime_error("invalid controlled BNV thermal accumulator");
    if(mode_==Mode::Coupled)rhs.AddTo(State::StateTag::Thermal,0,out.x_dot_s);
    for(std::size_t i=0;i<2;++i)rhs.AddTo(State::StateTag::Chem,i,out.eta_dot_MeV_s[i]);
}
} // namespace CompactStar::Physics::BNV

