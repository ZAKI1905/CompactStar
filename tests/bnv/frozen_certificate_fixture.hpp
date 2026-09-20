#pragma once

#include <CompactStar/Physics/BNV/EquilibriumBaryonTangent.hpp>
#include <CompactStar/Physics/BNV/FrozenBnvValidity.hpp>
#include <array>
#include <fstream>
#include <sstream>

namespace Phase6A1Test
{
namespace BNV=CompactStar::Physics::BNV;
namespace AN=CompactStar::Analysis;

inline std::vector<BNV::FrozenValiditySample> ReadFrozenCertificateRows(
    const std::filesystem::path& certificate_tsv)
{
    std::ifstream input(certificate_tsv);if(!input)throw std::runtime_error("missing frozen certificate TSV");
    std::string line;std::getline(input,line);
    if(line!="index\tfractional_depletion\tB_solved_count\tB_target_residual_count\tfinal_bracket_width_count\tDeltaN_n_over_N_n\tDeltaN_e_over_N_e\tDeltaN_mu_over_N_mu\tquantity\tdrift_bound\tthreshold\tutilization")
        throw std::runtime_error("unexpected frozen certificate schema");
    std::vector<BNV::FrozenValiditySample> samples(21);std::array<std::size_t,21> counts{};
    while(std::getline(input,line))
    {
        std::istringstream row(line);std::size_t index;std::string quantity;
        double fraction,B,residual,width,dNn,dNe,dNmu,drift,threshold,utilization;
        if(!(row>>index>>fraction>>B>>residual>>width>>dNn>>dNe>>dNmu>>quantity>>drift>>threshold>>utilization)||index>=21)
            throw std::runtime_error("malformed frozen certificate row");
        auto& sample=samples[index];
        if(counts[index]++==0)
        {sample.fractional_depletion=fraction;sample.B_solved_count=B;sample.B_target_residual_count=residual;sample.final_bracket_width_count=width;sample.DeltaN_over_N={dNn,dNe,dNmu};}
        else if(sample.fractional_depletion!=fraction||sample.B_solved_count!=B||
                sample.B_target_residual_count!=residual||sample.final_bracket_width_count!=width||
                sample.DeltaN_over_N!=std::array<double,3>{dNn,dNe,dNmu})
            throw std::runtime_error("inconsistent frozen certificate row identity");
        if(sample.utilization.find(quantity)!=sample.utilization.end())throw std::runtime_error("duplicate frozen certificate quantity");
        sample.drift_bound[quantity]=drift;sample.threshold[quantity]=threshold;sample.utilization[quantity]=utilization;
    }
    for(auto count:counts)if(count!=22)throw std::runtime_error("incomplete frozen certificate quantity set");
    return samples;
}

inline std::shared_ptr<const BNV::FrozenBnvValidityMonitor> LoadFrozenMonitor(
    const std::filesystem::path& certificate_tsv,
    const std::shared_ptr<const BNV::EquilibriumBaryonTangent>& tangent,
    const std::shared_ptr<const BNV::BnvDependencyToken>& token)
{
    if(!tangent||!token)throw std::runtime_error("missing frozen monitor dependency");
    auto certificate=std::make_shared<const BNV::FrozenSensitivityCertificate>(
      tangent->B0Count(),tangent->StarIdentity(),tangent->DomainIdentity(),
      "Phase-6A-1 achieved-B 21-star error-aware frozen sensitivity certificate",
      token,ReadFrozenCertificateRows(certificate_tsv));
    return std::make_shared<const BNV::FrozenBnvValidityMonitor>(certificate);
}
}
