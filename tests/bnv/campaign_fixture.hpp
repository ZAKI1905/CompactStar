#pragma once

#include "../rotochemical/coupled.hpp"
#include "controlled_neutron_sink_fixture.hpp"
#include "energy_partition_fixtures.hpp"
#include "frozen_certificate_fixture.hpp"
#include <CompactStar/Analysis/EquilibriumBaryonTangent.hpp>
#include <CompactStar/Physics/BNV/ControlledBnvSecularDriver.hpp>
#include <CompactStar/Physics/BNV/StaticZeroSpinHistory.hpp>
#include <filesystem>
#include <fstream>
#include <map>

namespace Phase6A1Campaign
{
namespace BNV=CompactStar::Physics::BNV;
namespace AN=CompactStar::Analysis;

inline std::shared_ptr<const AN::EquilibriumBaryonTangent> Tangent(
    const ControlledFixture& f,const std::filesystem::path& profile,const std::filesystem::path& work)
{
    auto src=source(profile/"freegas.tsv");auto numbers=ParticleNumbers::Compute(input(*f.central,src));numbers.RequireCurrent();
    NumberSequenceRecipe recipe{src,species,whole,1.10e15,1.095e15,1.105e15,
      "Structure-1 monotone midpoint smooth central branch",{.001,.0005,.00025},80000,
      (work/"tangent-sequence").string(),tail_bound};
    recipe.tail_policy_identity="positive-source pe comparison inequalities";recipe.tail_policy_revision="PB13-comparison-v1";
    auto sequence=std::make_shared<const EquilibriumSequenceNumberDerivative>(EquilibriumSequenceNumberDerivative::Compute(recipe));
    const auto find=[&](const std::string& label){for(std::size_t i=0;i<species.size();++i)if(species[i].label==label)return i;throw std::runtime_error("missing number species");};
    const double B0=numbers.Values()[find("10")]+numbers.Values()[find("11")];
    return std::make_shared<const AN::EquilibriumBaryonTangent>(AN::EquilibriumBaryonTangent::Compute(
      sequence,B0,"Structure-1 rho_c=1.10e15 radial80000 EOS8192",AN::EquilibriumBaryonTangent::SerializeDomainIdentity(whole)));
}

inline RC::RunQualification Qualification(const std::filesystem::path& profile,const std::filesystem::path& certificate,
    const std::filesystem::path& entry)
{
    RC::RunQualification q;q.purpose=RC::RunPurpose::AnalyticControl;q.radial_resolution=80000;q.eos_resolution=8192;q.rho_c_g_cm3=1.10e15;
    q.profile_path=(profile/"profile.tsv").string();q.model_path=(profile/"model.txt").string();q.eos_path=(profile/"freegas.tsv").string();q.certificate_path=certificate.string();q.entry_manifest_path=entry.string();q.entry_manifest_sha256=RC::FrozenSource(entry).sha256;return q;
}

inline double Column(const std::filesystem::path& path,const std::string& name)
{
    std::ifstream in(path);if(!in)throw std::runtime_error("missing coefficient table");
    std::string header,row;if(!std::getline(in,header)||!std::getline(in,row))throw std::runtime_error("empty coefficient table");
    std::istringstream hs(header),rs(row);std::vector<std::string> names,values;std::string item;
    while(hs>>item)names.push_back(item);while(rs>>item)values.push_back(item);
    if(names.size()!=values.size())throw std::runtime_error("coefficient table shape mismatch");
    for(std::size_t i=0;i<names.size();++i)if(names[i]==name){const double v=std::stod(values[i]);if(!std::isfinite(v))throw std::runtime_error("nonfinite coefficient");return v;}
    throw std::runtime_error("coefficient column missing: "+name);
}

class ExactZeroHistory final:public BNV::OrdinaryMatterBnvHistory
{
  public:
    ExactZeroHistory(double B,std::string domain,std::string star,std::string sequence,
      std::shared_ptr<const BNV::BnvDependencyToken> token)
      :B_(B),domain_(std::move(domain)),star_(std::move(star)),sequence_(std::move(sequence)),token_(std::move(token)),generation_(token_?token_->generation:0){RequireCurrent();}
    BNV::OrdinaryMatterBnvHistorySample Sample(double t)const override
    {RequireCurrent();BNV::OrdinaryMatterBnvHistorySample s;s.epoch_s=t;s.B_count=B_;s.source_identity=identity_;s.channel_provenance="exact disabled BNV bundle";s.revision_identity="phase6a1-zero-source-v1";s.domain_identity=domain_;s.star_identity=star_;s.sequence_state_identity=sequence_;s.Validate();return s;}
    void RequireCurrent()const override{if(!token_||!token_->alive||token_->generation!=generation_)throw std::runtime_error("stale exact-zero source");}
    const std::string& Identity()const override{RequireCurrent();return identity_;}
  private:
    double B_;std::string domain_,star_,sequence_;std::shared_ptr<const BNV::BnvDependencyToken> token_;std::uint64_t generation_;
    const std::string identity_="Phase-6A-1 exact zero BNV source for matched controls";
};

class ExactZeroPartition final:public BNV::IDirectBnvEnergyPartition
{
  public:
    explicit ExactZeroPartition(std::shared_ptr<const BNV::BnvDependencyToken> token):token_(std::move(token)),generation_(token_?token_->generation:0){RequireCurrent();}
    std::vector<BNV::DirectEventEnergy> Evaluate(const BNV::OrdinaryMatterBnvHistorySample& s,const BNV::ActualOrdinaryPotential&)const override
    {RequireCurrent();if(!s.events.empty())throw std::runtime_error("zero partition received event");return {};}
    void RequireCurrent()const override{if(!token_||!token_->alive||token_->generation!=generation_)throw std::runtime_error("stale exact-zero partition");}
    const std::string& Identity()const override{RequireCurrent();return identity_;}
    const std::string& FiniteTemperatureWeightingClass()const override{static const std::string x="no events; zero BNV control";return x;}
    double OmittedFiniteTemperaturePowerErgPerSecond(const BNV::OrdinaryMatterBnvHistorySample&,double)const override{return 0;}
  private:
    std::shared_ptr<const BNV::BnvDependencyToken> token_;std::uint64_t generation_;const std::string identity_="P-OFF exact zero partition";
};

struct RunCard
{
    std::string identity,partition;
    double fractional_drive_per_year=0,duration_year=0;
    std::size_t checkpoints=0;
    bool reaction_free=false;
};

inline const std::array<RunCard,4>& RunCards()
{
    static const std::array<RunCard,4> cards{{
      {"RF-P0-TRANSIENT-v1","P0",-1e-13,1e6,1025,true},
      {"CPL-P0-TRANSIENT-v1","P0",-1e-13,1e5,2049,false},
      {"CPL-P1-TRANSIENT-v1","P1",-1e-13,1e5,2049,false},
      {"CPL-P2-LINEAR-QSS-v1","P2",-1e-12,5e5,8193,false}}};
    return cards;
}

} // namespace Phase6A1Campaign
