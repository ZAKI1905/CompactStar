#pragma once

#include <CompactStar/Physics/BNV/FrozenControlledBnvRunContext.hpp>
#include <CompactStar/Physics/Evolution/RHSAccumulator.hpp>

namespace CompactStar::Physics::BNV
{

class ControlledBnvSecularDriver final:public IDriver
{
  public:
    enum class Mode { Coupled, ChemicalOnly, ReactionFreeControl };
    explicit ControlledBnvSecularDriver(std::shared_ptr<const FrozenControlledBnvRunContext>,Mode=Mode::Coupled);
    std::string Name() const override {return "Phase6A1ControlledBnvSecularEvolution";}
    const std::vector<State::StateTag>& DependsOn() const override;
    const std::vector<State::StateTag>& Updates() const override{return DependsOn();}
    ControlledBnvEvaluation Evaluate(double,const Evolution::StateVector&,const Evolution::DriverContext&) const;
    void AccumulateRHS(double,const Evolution::StateVector&,Evolution::RHSAccumulator&,const Evolution::DriverContext&) const override;
  private:
    std::shared_ptr<const FrozenControlledBnvRunContext> context_;
    Mode mode_;
};

} // namespace CompactStar::Physics::BNV
