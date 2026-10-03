#pragma once
// Fixture assembly promoted from reviewed v3, no private cache access or scratch solvers.
#define PHASE6A1_CAMPAIGN_SUPPORT_ONLY
#include "adr0017_production_qualification.cpp"
namespace BreakpointTest {
using State=BNV::CheckpointState;
using Ctx=std::shared_ptr<const BNV::FrozenControlledBnvRunContext>;
struct Pair{Ctx source,control;};
Pair BuildBoth(const Inputs& inputs)
{
    const auto& card=Campaign::RunCards().back();
    require(card.identity=="CPL-P2-LINEAR-QSS-v1"&&card.partition=="P2"&&card.fractional_drive_per_year==-1.0e-12
      &&card.duration_year==5.0e5&&card.checkpoints==8193&&!card.reaction_free,"bounded P2 card changed");
    auto fixture=Fixture(inputs.profile,inputs.certificate.string(),inputs.work/"owning-star",80000);
    auto channels=Channels(fixture);
    EOS::CompOSE_Thermo::Options options;options.Tmin_for_derivative_MeV=0;options.clamp_to_domain=false;
    auto thermal=std::make_shared<const RC::FrozenThermalSource>(inputs.thermal,options,
      "controlled mathematical fixed-background free-gas entropy; qualified radial80000");
    auto tangent=Campaign::Tangent(fixture,inputs.profile,inputs.work);
    auto bnv_token=std::make_shared<BNV::BnvDependencyToken>();
    auto monitor=Phase6A1Test::LoadFrozenMonitor(inputs.frozen,tangent,bnv_token);
    auto run_token=std::make_shared<RC::RunDependencyToken>();
    auto spin=std::make_shared<const BNV::StaticZeroSpinHistory>(run_token);
    auto ordinary=Context(fixture,channels,thermal,spin,Campaign::Qualification(inputs.profile,inputs.certificate,inputs.entry),run_token);
    const double mu_B=Campaign::Column(inputs.coefficients,"mu_B_inf");
    const double Bdot=card.fractional_drive_per_year*tangent->B0Count()/Year;
    require(Bdot==-2.4136520263641375e37,"bounded P2 Bdot changed");
    const std::string fate_id="phase6a1-P2-terminal-fate-v1";
    BNV::ProductFateLedger fate(fate_id,{{"controlled-terminal-product",BNV::TerminalProductFate::BoundInert,1,"abstract-neutron-disappearance"}});
    auto partition=Partition(card,fixture,bnv_token,fate_id);
    auto history=std::make_shared<const Phase6A1Test::UniformProperNeutronSinkHistory>(fixture.central,tangent->B0Count(),Bdot,
      tangent->DomainIdentity(),tangent->StarIdentity(),tangent->SequenceStateIdentity(),partition->Identity(),fate_id,bnv_token);
    const std::string potential="frozen Structure-1 B0 equilibrium mu_B plus governed moving-reference actual-potential correction";
    auto source=std::make_shared<const BNV::FrozenControlledBnvRunContext>(ordinary,tangent,history,partition,fate,monitor,channels,mu_B,potential,card.identity);
    auto zero_history=std::make_shared<const Campaign::ExactZeroHistory>(tangent->B0Count(),tangent->DomainIdentity(),
      tangent->StarIdentity(),tangent->SequenceStateIdentity(),bnv_token);
    auto zero_partition=std::make_shared<const Campaign::ExactZeroPartition>(bnv_token);
    BNV::ProductFateLedger zero_fate("zero-control-fate",{{"none",BNV::TerminalProductFate::BoundInert,1,"zero-control"}});
    auto control=std::make_shared<const BNV::FrozenControlledBnvRunContext>(ordinary,tangent,zero_history,zero_partition,zero_fate,
      monitor,channels,mu_B,potential,card.identity+"-MATCHED-CONTROL");
    return {source,control};
}

struct Archived{int seq;double t0,t1,hnext;State y;std::size_t rej;};
std::vector<Archived> ReadAccepted(const std::filesystem::path& p)
{
    std::ifstream f(p);std::string line;std::getline(f,line);std::vector<Archived> v;
    while(std::getline(f,line)){auto s=Split(line);v.push_back({std::stoi(s[0]),std::stod(s[1]),std::stod(s[2]),std::stod(s[3]),
      {{std::stod(s[4]),std::stod(s[5]),std::stod(s[6])}},std::stoull(s[7])});}
    return v;
}

}
