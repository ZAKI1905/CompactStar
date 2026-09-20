#include <CompactStar/Physics/BNV/DirectEnergyLedger.hpp>
#include <iostream>
#include <limits>

namespace BNV=CompactStar::Physics::BNV;
namespace RC=CompactStar::Physics::Rotochemical;
void need(bool x,const char*m){if(!x)throw std::runtime_error(m);}
template<class F>void refuse(F f,const char*m){bool x=false;try{f();}catch(...){x=true;}need(x,m);}
int main(){try{
 BNV::OrdinaryMatterBnvHistorySample s;s.epoch_s=0;s.B_count=1e57;s.Bdot_count_s=-1e30;s.source_count_s={-1e30,0,0};s.proton_source_count_s=0;s.source_identity="direct fixture";s.channel_provenance="generic";s.revision_identity="v1";s.domain_identity="D";s.star_identity="S";s.sequence_state_identity="Q";s.events.push_back({"e","c","f","p",{-1,0,0},1e30});s.Validate();
 BNV::ActualOrdinaryPotential p;p.equilibrium_inf_MeV={900,900,900};p.actual_inf_MeV={901,899,898};p.mu_B_inf_MeV=900;p.mu_n_actual_inf_MeV=901;p.provenance="independent typed potential";
 BNV::ProductFateLedger fate("f",{{"x",BNV::TerminalProductFate::PromptEscape,.25,"c"},{"y",BNV::TerminalProductFate::BoundInert,.75,"c"}});
 std::vector<BNV::DirectEventEnergy> e{{"e","p","f",100,80,20,0,true,.8,1001.25}};
 auto r=BNV::BnvDirectEnergyLedger::Evaluate(s,p,e,fate);need(r.Qactual_event_inf_MeV[0]==801,"inclusive direct event");need(r.power_actual_erg_s==801e30*RC::MeVToErg,"exact direct conversion");
 need(r.power_actual_erg_s!=r.power_actual_MeV_s,"M13 omitted direct conversion survived");
 need(r.power_actual_erg_s!=r.power_actual_MeV_s*RC::MeVToErg*RC::MeVToErg,"M14 doubled direct conversion survived");
 auto bad=e;bad[0].EX_inf_MeV=21;refuse([&]{BNV::BnvDirectEnergyLedger::Evaluate(s,p,bad,fate);},"M15 fluid/star confusion");
 refuse([&]{BNV::ProductFateLedger("f",{{"x",BNV::TerminalProductFate::PromptEscape,.5,"c"},{"x",BNV::TerminalProductFate::BoundInert,.5,"c"}});},"M16 duplicate fate");
 bad=e;bad[0].Qlocal_MeV=1;refuse([&]{BNV::BnvDirectEnergyLedger::Evaluate(s,p,bad,fate);},"local/infinity mismatch");
 std::cout<<"DIRECT_ENERGY PASS inclusive event local/infinity and conversion\nMUTATION M13 DETECTED unit oracle\nMUTATION M14 DETECTED unit oracle\nMUTATION M15 DETECTED fate-energy closure\nMUTATION M16 DETECTED terminal duplicate\n";return 0;
 }catch(const std::exception&e){std::cerr<<"STOP "<<e.what()<<'\n';return 1;}}
