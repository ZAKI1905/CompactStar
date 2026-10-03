#pragma once
#include <CompactStar/Physics/Evolution/DriverContext.hpp>
#include <CompactStar/Physics/Evolution/StarContext.hpp>
#include <CompactStar/Physics/Rotochemical/ChemicalImbalanceState.hpp>
#include <cmath>
#include <stdexcept>

namespace CompactStar::Physics::BNV
{
// Read-only adapter. The snapshot is bound to a currentness-checked RHS owner.
class CstarSegments
{
public:
    explicit CstarSegments(Evolution::DriverContext context):context_(context),
      cells_(context.star->ReadHeatCapacityCells(*context.thermo,context.geo)){}
    std::vector<double> Knots() const
    {
        std::vector<double> x;
        for(double T:cells_.temperature_MeV)
            x.push_back(std::log((T/Rotochemical::BoltzmannMeVPerK)/1e8));
        return x;
    }
    const Evolution::StarContext::HeatCapacityCells& Cells() const{return cells_;}
    void Extend(double x,int cell,double& xdot) const
    {
        if(cell<0)return;
        const double T=1e8*std::exp(x)*Rotochemical::BoltzmannMeVPerK;
        if(T>=cells_.temperature_MeV.at(cell)&&T<=cells_.temperature_MeV.at(cell+1))return;
        const double smooth=cells_.SmoothCapacity(T,cell);
        if(!(smooth>0&&std::isfinite(smooth)))throw std::runtime_error("invalid Cstar smooth extension");
        const double actual=context_.star->HeatCapacityStar_Tinf(T,*context_.thermo,context_.geo);
        xdot*=actual/smooth;
    }
private:
    Evolution::DriverContext context_;
    const Evolution::StarContext::HeatCapacityCells cells_;
};
}
