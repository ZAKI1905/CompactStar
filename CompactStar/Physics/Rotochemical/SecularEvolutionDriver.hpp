#pragma once
#include <CompactStar/Physics/Rotochemical/FrozenRotochemicalRunContext.hpp>
#include <CompactStar/Physics/Evolution/RHSAccumulator.hpp>
namespace CompactStar::Physics::Rotochemical
{
class SecularEvolutionDriver final:public IDriver
{
  public:
    enum class Mode {Coupled,ChemicalOnly};
    explicit SecularEvolutionDriver(std::shared_ptr<const FrozenRotochemicalRunContext> context,Mode mode=Mode::Coupled):context_(std::move(context)),spin_(context_?context_->SpinOwner():nullptr),mode_(mode){if(!context_)throw std::runtime_error("missing frozen context");context_->RequireCheapCurrent();}
    std::string Name()const override{return "ControlledSecularRotochemicalEvolution";}
    const std::vector<State::StateTag>& DependsOn()const override{static const std::vector<State::StateTag> tags{State::StateTag::Thermal,State::StateTag::Chem};return tags;}
    const std::vector<State::StateTag>& Updates()const override{return DependsOn();}
    SecularEvaluation Evaluate(double t,const Evolution::StateVector& state,const Evolution::DriverContext& ctx)const{return context_->Evaluate(t,state,ctx,spin_.get());}
    void AccumulateRHS(double t,const Evolution::StateVector& state,Evolution::RHSAccumulator& rhs,const Evolution::DriverContext& ctx)const override
    {
        context_->RequireOwners(ctx,spin_.get());
        if(rhs.Block(State::StateTag::Thermal).size()!=1||rhs.Block(State::StateTag::Chem).size()!=2)throw std::runtime_error("malformed RHS blocks");
        auto out=Evaluate(t,state,ctx);context_->RequireOwners(ctx,spin_.get());
        for(size_t i=0;i<2;++i)if(!std::isfinite(rhs.Block(State::StateTag::Chem)[i]+out.eta_dot_MeV_s[i]))throw std::runtime_error("invalid chemical accumulator");
        if(mode_==Mode::Coupled&&!std::isfinite(rhs.Block(State::StateTag::Thermal)[0]+out.x_dot_s))throw std::runtime_error("invalid thermal accumulator");
        // All fallible dependency and scientific work precedes the first update.
        if(mode_==Mode::Coupled)rhs.AddTo(State::StateTag::Thermal,0,out.x_dot_s);
        for(auto c:{BetaChannel::Npe,BetaChannel::NpMu})rhs.AddTo(State::StateTag::Chem,ChannelIndex(c),out.eta_dot_MeV_s[ChannelIndex(c)]);
    }
  private:
    const std::shared_ptr<const FrozenRotochemicalRunContext> context_;
    const std::shared_ptr<const PrescribedSpinHistory> spin_;
    const Mode mode_;
};
}
