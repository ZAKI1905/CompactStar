#pragma once
#include <CompactStar/Physics/Rotochemical/UrcaChannelNormalization.hpp>

namespace CompactStar::Physics::Rotochemical
{
struct UrcaMetric { double nu, lambda; };
struct UrcaIntegrationRequest
{
    std::shared_ptr<const Analysis::GlobalChemicalNumberResponse> chemical_domain;
    std::function<UrcaMetric(double)> metric;
    // Explicit partition; includes regular centre and support endpoints as needed.
    std::vector<double> radial_partition_km;
    std::string domain_identity, metric_identity;
    UrcaProcessSelection selection = UrcaProcessSelection::ControlledModifiedOnly();
    std::vector<UrcaChannelNormalization> normalizations;
    unsigned order=16;
};
struct UrcaCoefficientEntry
{
    UrcaProcess process;
    double luminosity_erg_s_Kq=0;
    std::vector<UrcaSupportInterval> support;
    std::string normalization_identity;
};
// Layer 2: one immutable channel-resolved global normalization authority.
// No F/H, rates, thermal state or torque dependencies.
class GlobalUrcaChannelCoefficient
{
  public:
    GlobalUrcaChannelCoefficient(const GlobalUrcaChannelCoefficient&) = default;
    GlobalUrcaChannelCoefficient(GlobalUrcaChannelCoefficient&&) = default;
    GlobalUrcaChannelCoefficient& operator=(const GlobalUrcaChannelCoefficient&) = delete;
    GlobalUrcaChannelCoefficient& operator=(GlobalUrcaChannelCoefficient&&) = delete;
    static GlobalUrcaChannelCoefficient Compute(const UrcaIntegrationRequest&);
    void RequireCurrent() const;
    double LuminosityCoefficient(UrcaProcess p) const { RequireCurrent();return entries_[ProcessIndex(p)].luminosity_erg_s_Kq; }
    const UrcaCoefficientEntry& Entry(UrcaProcess p) const { RequireCurrent();return entries_[ProcessIndex(p)]; }
    const UrcaProcessSelection& Selection() const { RequireCurrent();return selection_; }
    const std::shared_ptr<const Analysis::GlobalChemicalNumberResponse>& ChemicalDomain() const { RequireCurrent();return domain_; }
    const std::string& DomainIdentity() const { RequireCurrent();return domain_identity_; }
    const std::string& MetricIdentity() const { RequireCurrent();return metric_identity_; }
    const std::vector<double>& PartitionKm() const { RequireCurrent();return partition_; }
    unsigned Order() const { return order_; }
  private:
    GlobalUrcaChannelCoefficient() = default;
    std::shared_ptr<const Analysis::GlobalChemicalNumberResponse> domain_;
    std::string domain_identity_,metric_identity_;
    std::vector<double> partition_;
    unsigned order_=16;
    UrcaProcessSelection selection_{ {} };
    std::array<UrcaCoefficientEntry,4> entries_{{{UrcaProcess::De,0,{},{}},{UrcaProcess::Dmu,0,{},{}},{UrcaProcess::Me,0,{},{}},{UrcaProcess::Mmu,0,{},{}}}};
};
}
