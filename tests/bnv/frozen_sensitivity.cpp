#include <CompactStar/Physics/BNV/FrozenBnvValidity.hpp>
#include <iostream>

namespace BNV=CompactStar::Physics::BNV;
void need(bool x,const char*m){if(!x)throw std::runtime_error(m);}
int main(){try{
 const double B0=7.6169065187188905e56;auto token=std::make_shared<BNV::BnvDependencyToken>();std::vector<BNV::FrozenValiditySample> rows;
 const std::vector<std::string> names{"t_n","t_e","t_mu","Z_row_npe","Z_row_npmu","Cstar","Ltilde_Me","Ltilde_Mmu","mu_B","mu_n_profile","P0_average","P1_average","P2_average","N_n","N_e","N_mu","species_support","metric_structure","radius","surface_gravity","envelope","Tsurface_inf"};
 for(int i=0;i<21;++i){BNV::FrozenValiditySample s;s.fractional_depletion=-5e-8*i;s.B_solved_count=B0*(1+s.fractional_depletion);s.B_target_residual_count=s.B_solved_count-B0*(1+s.fractional_depletion);s.final_bracket_width_count=0;for(const auto&n:names){s.threshold[n]=1;s.drift_bound[n]=.025*i;s.utilization[n]=s.drift_bound[n]/s.threshold[n];}rows.push_back(s);}
 auto cert=std::make_shared<const BNV::FrozenSensitivityCertificate>(B0,"S","D","synthetic monitor semantics only",token,rows);BNV::FrozenBnvValidityMonitor monitor(cert);
 auto below=monitor.RequireValid(B0*(1-9.999999e-7));need(below.valid&&below.max_utilization<=1,"below-limit validity");bool refused=false;try{monitor.RequireValid(B0*(1-1.000001e-6));}catch(...){refused=true;}need(refused,"M20 over-limit did not refuse");
 std::cout<<"FROZEN_MONITOR_SEMANTICS PASS 21 grid and trial-state stop\nMUTATION M20 DETECTED just-over ceiling refused\n";return 0;
 }catch(const std::exception&e){std::cerr<<"STOP "<<e.what()<<'\n';return 1;}}
