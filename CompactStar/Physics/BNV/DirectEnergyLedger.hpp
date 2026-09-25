#pragma once

#include <CompactStar/Physics/BNV/OrdinaryMatterSource.hpp>
#include <CompactStar/Physics/BNV/ProductFate.hpp>
#include <CompactStar/Physics/Rotochemical/ChemicalImbalanceState.hpp>
#include <CompactStar/Physics/BNV/EquilibriumBaryonTangent.hpp>
#include <array>
#include <memory>
#include <string>
#include <vector>

namespace CompactStar::Physics::BNV
{

struct ActualOrdinaryPotential
{
    std::array<double,3> equilibrium_inf_MeV{};
    std::array<double,3> actual_inf_MeV{};
    double mu_B_inf_MeV=0, mu_n_actual_inf_MeV=0;
    double reconstruction_residual_MeV=0;
    std::string provenance;
};

struct DirectEventEnergy
{
    std::string event_id, partition_identity, product_fate_identity;
    double Eesc_fluid_inf_MeV=0, Eesc_star_inf_MeV=0, EX_inf_MeV=0;
    double Ein_explicit_inf_MeV=0;
    bool local_frame_check=false;
    double exp_phi=1, Qlocal_MeV=0;
};

class IDirectBnvEnergyPartition
{
  public:
    virtual ~IDirectBnvEnergyPartition()=default;
    virtual std::vector<DirectEventEnergy> Evaluate(
        const OrdinaryMatterBnvHistorySample&,const ActualOrdinaryPotential&) const=0;
    virtual void RequireCurrent() const=0;
    virtual const std::string& Identity() const=0;
    virtual const std::string& FiniteTemperatureWeightingClass() const=0;
    virtual double OmittedFiniteTemperaturePowerErgPerSecond(
        const OrdinaryMatterBnvHistorySample&,double Tinf_K) const=0;
};

struct DirectEnergyResult
{
    ActualOrdinaryPotential potential;
    std::vector<double> Qeq_event_inf_MeV,Qactual_event_inf_MeV;
    double power_eq_MeV_s=0,power_actual_MeV_s=0;
    double power_eq_erg_s=0,power_actual_erg_s=0;
    double local_infinity_residual_MeV=0;
    double escape_fate_residual_MeV=0;
};

class BnvDirectEnergyLedger final
{
  public:
    static ActualOrdinaryPotential ActualPotential(
        double mu_B_inf_MeV,
        const Rotochemical::ChemicalImbalanceState&,
        const EquilibriumBaryonTangent&,
        std::string provenance);
    static ActualOrdinaryPotential ActualPotential(
        double mu_B_inf_MeV,
        const Rotochemical::ChemicalImbalanceState&,
        const ValidatedTangentSnapshot&,
        std::string provenance);
    static DirectEnergyResult Evaluate(
        const OrdinaryMatterBnvHistorySample&,
        const ActualOrdinaryPotential&,
        const std::vector<DirectEventEnergy>&,
        const ProductFateLedger&);
};

} // namespace CompactStar::Physics::BNV
