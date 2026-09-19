#include <CompactStar/Physics/BNV/FrozenBnvValidity.hpp>
#include "frozen_certificate_fixture.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace BNV=CompactStar::Physics::BNV;

void Need(bool value,const char* message){if(!value)throw std::runtime_error(message);}
template<class F>void MustRefuse(F&& function,const char* message)
{bool refused=false;try{function();}catch(...){refused=true;}Need(refused,message);}

int main(int argc,char** argv)
{
 try
 {
    Need(argc==2,"certificate.tsv");
    auto samples=Phase6A1Test::ReadFrozenCertificateRows(argv[1]);
    const double B0=samples.front().B_solved_count;
    auto token=std::make_shared<BNV::BnvDependencyToken>();
    auto certificate=std::make_shared<const BNV::FrozenSensitivityCertificate>(B0,
      "Structure-1 rho_c=1.10e15 radial80000 EOS8192",
      "WholeStar P=0; source support and active Track-R charts",
      "Phase-6A-1 achieved-B 21-star error-aware frozen sensitivity certificate",
      token,samples);
    BNV::FrozenBnvValidityMonitor monitor(certificate);
    for(const auto& sample:samples)
    {
        const auto result=monitor.RequireValid(sample.B_solved_count);
        Need(result.valid&&result.max_utilization<=1,"certified sample refused");
    }
    const double inside=B0*(1-1e-6+128*std::numeric_limits<double>::epsilon());
    const double outside=B0*(1-1e-6-128*std::numeric_limits<double>::epsilon());
    Need(monitor.RequireValid(inside).valid,"just-inside depletion ceiling refused");
    MustRefuse([&]{monitor.RequireValid(outside);},"M20 outside depletion ceiling survived");
    auto mutant_samples=samples;mutant_samples.back().drift_bound["N_mu"]=2*mutant_samples.back().threshold["N_mu"];
    mutant_samples.back().utilization["N_mu"]=2;
    MustRefuse([&]{BNV::FrozenSensitivityCertificate mutant(B0,certificate->StarIdentity(),certificate->DomainIdentity(),
      "M20 over-budget mutant",token,mutant_samples);},"M20 over-budget certificate survived");
    ++token->generation;
    MustRefuse([&]{certificate->RequireCurrent();},"stale frozen certificate survived");
    std::cout<<std::setprecision(17)<<"BA13_RUNTIME_MONITOR PASS just_inside_fraction "<<(inside-B0)/B0
      <<" just_outside_fraction "<<(outside-B0)/B0<<"\nMUTATION M20 DETECTED pre-RHS frozen validity refusal\n";
    return 0;
 }
 catch(const std::exception& error){std::cerr<<"STOP "<<error.what()<<'\n';return 1;}
}
