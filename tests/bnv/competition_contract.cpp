#include <CompactStar/Physics/BNV/ActualPotentialPartition.hpp>
#include <CompactStar/Physics/Rotochemical/ParameterizedSpinHistory.hpp>
#include <iostream>
#include <limits>

namespace BNV=CompactStar::Physics::BNV;
namespace RC=CompactStar::Physics::Rotochemical;
void need(bool x,const char* m){if(!x)throw std::runtime_error(m);}
template<class F>void refuses(F f,const char* m){bool ok=false;try{f();}catch(const std::exception&){ok=true;}need(ok,m);}
int main(){try{
    auto st=std::make_shared<RC::RunDependencyToken>();
    RC::ParameterizedDipoleHistory spin(.02,1e-19,st);
    const auto a=spin.Sample(0),b=spin.Sample(1e15);const double pi=std::acos(-1.);
    need(std::abs(a.omega_rad_s-2*pi/.02)<1e-12&&a.omega_dot_rad_s2<0,"spin sign and normalization");
    need(b.omega_rad_s<a.omega_rad_s&&b.omega_dot_rad_s2<0,"monotone spin-down");
    const double period=2*pi/b.omega_rad_s,pdot=-2*pi*b.omega_dot_rad_s2/std::pow(b.omega_rad_s,2);
    need(std::abs(period*pdot/(.02*1e-19)-1)<1e-14,"constant PPdot identity");
    refuses([&]{spin.Sample(-1);},"negative time guard");
    refuses([&]{RC::ParameterizedDipoleHistory bad(.02,-1e-19,st);},"spin-up not admitted");
    st->generation++;refuses([&]{spin.Sample(0);},"stale spin token");
    auto token=std::make_shared<BNV::BnvDependencyToken>();
    BNV::ActualOrdinaryPotential pot;pot.mu_B_inf_MeV=868.;pot.mu_n_actual_inf_MeV=868.+3e-12;
    pot.equilibrium_inf_MeV={868.,868.,868.};pot.actual_inf_MeV={pot.mu_n_actual_inf_MeV,867.,866.};
    pot.departure_inf_MeV=std::array<double,3>{3e-12,-1.,-2.};pot.provenance="independent algebra fixture";
    BNV::ProductFateLedger fate("escape",{{"escaped",BNV::TerminalProductFate::PromptEscape,1.,"neutron-sink"}});
    for(double f:{0.,1e-12,1e-6,1.}){
        BNV::ActualPotentialPartition part(f,"escape",token);
        BNV::OrdinaryMatterBnvHistorySample s;s.B_count=1e57;s.Bdot_count_s=-1e30;s.source_count_s={-1e30,0,0};
        s.source_identity="s";s.channel_provenance="test";s.revision_identity="r";s.domain_identity="d";s.star_identity="star";s.sequence_state_identity="seq";
        s.events.push_back({"e","neutron-sink","escape",part.Identity(),{-1,0,0},1e30});
        auto energies=part.Evaluate(s,pot);auto result=BNV::BnvDirectEnergyLedger::Evaluate(s,pot,energies,fate);
        const long double reference=static_cast<long double>(f)*pot.mu_n_actual_inf_MeV;
        need(std::abs(result.Qactual_event_inf_MeV[0]-reference)<=4*std::numeric_limits<double>::epsilon()*std::max(1e-300L,reference),"tiny-f direct budget loses precision");
        if(f==0){need(result.power_actual_erg_s==0,"P-A zero direct power");need(result.Qeq_event_inf_MeV[0]==-3e-12,"chemical departure lost to GeV cancellation");}
        need(result.escape_fate_residual_MeV==0,"escape closure");
        auto bad=energies;bad[0].Eesc_fluid_inf_MeV+=1;bad[0].Eesc_star_inf_MeV+=1;
        refuses([&]{BNV::BnvDirectEnergyLedger::Evaluate(s,pot,bad,fate);},"missing/altered escape energy survived");
    }
    std::cout<<"PASS prescribed history, stale/negative guards, P-A f=0/1/tiny, precise chemical departure, escape mutation\n";
    return 0;
}catch(const std::exception& e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}}
