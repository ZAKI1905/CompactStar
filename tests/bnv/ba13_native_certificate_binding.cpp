// Candidate-only static binding check. No BNV context, RHS or trajectory.
#include "../rotochemical/fixture.hpp"
#include "frozen_certificate_fixture.hpp"
#include <iostream>

int main(int argc,char** argv)
{
 try {
    require(argc==5,"native-profile chemical-certificate frozen-certificate fresh-work");
    const std::filesystem::path profile=argv[1],chemical=argv[2],frozen=argv[3],work=argv[4];
    auto fixture=Fixture(profile,chemical.string(),work,80000);
    fixture.g->RequireCurrent();fixture.z->RequireCurrent();
    const auto src=source(profile/"freegas.tsv");
    auto numbers=ParticleNumbers::Compute(input(*fixture.central,src));
    numbers.RequireCurrent();
    auto number=[&](const std::string& label){
      for(std::size_t i=0;i<species.size();++i)
        if(species[i].label==label)return numbers.Values()[i];
      throw std::runtime_error("missing governed number axis");
    };
    const double B0=number("10")+number("11");
    auto samples=Phase6A1Test::ReadFrozenCertificateRows(frozen);
    require(samples.front().B_solved_count==B0,"native frozen B0 differs from actual owning star");
    namespace BNV=CompactStar::Physics::BNV;
    auto token=std::make_shared<BNV::BnvDependencyToken>();
    const auto certificate=std::make_shared<const BNV::FrozenSensitivityCertificate>(
      B0,"native profile SHA256="+RC::FrozenSource(profile/"profile.tsv").sha256,
      BNV::EquilibriumBaryonTangent::SerializeDomainIdentity(whole),
      "unpromoted Linux BA13 historical-semantics static candidate",token,samples);
    BNV::FrozenBnvValidityMonitor monitor(certificate);
    require(monitor.RequireValid(B0).valid,"native frozen certificate refused owning star");
    std::cout<<"NATIVE_FROZEN_CERTIFICATE_EXACT_OWNING_PROFILE_AND_B0 PASS\n";
    ++token->generation;
    bool refused=false;try{certificate->RequireCurrent();}catch(const std::exception&){refused=true;}
    require(refused,"stale native frozen certificate survived");
    std::cout<<"NATIVE_FROZEN_CERTIFICATE_CURRENTNESS_REFUSAL PASS; no trajectory\n";
    return 0;
 } catch(const std::exception& e){std::cerr<<"STOP "<<e.what()<<'\n';return 1;}
}
