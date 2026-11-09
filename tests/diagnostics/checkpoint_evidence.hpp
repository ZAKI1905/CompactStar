#pragma once
// Observability only: serialize existing qualification fields; never change status.
#include <CompactStar/Physics/BNV/PassiveCheckpointOutput.hpp>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace CheckpointEvidence {
namespace BNV=CompactStar::Physics::BNV;
inline std::map<std::string,double> Values(const BNV::BnvDiagnostics& d) {
 return {{"P_dir_eq",d.P_dir_eq_erg_s},{"P_dir_actual",d.P_dir_actual_erg_s},{"LH",d.LH_erg_s},
 {"DeltaLnu",d.DeltaLnu_erg_s},{"DeltaPbeta",d.DeltaPbeta_erg_s},{"Lnu_eq",d.Lnu_eq_erg_s},
 {"Lnu_full",d.Lnu_full_erg_s},{"Lgamma",d.Lgamma_erg_s},{"Lother",d.Lother_erg_s},{"Pnet",d.Pnet_erg_s},
 {"mu_B",d.mu_B_inf_MeV},{"mu_n_actual",d.mu_n_actual_inf_MeV},{"sigma_e",d.sigma_count_s[0]},
 {"sigma_mu",d.sigma_count_s[1]},{"Echem",d.Echem_MeV},{"Cstar",d.Cstar_erg_K},{"Tinf",d.Tinf_K},
 {"Lout",d.L_out_fluid_inf_erg_s}};
}
inline double Gross(const BNV::BnvDiagnostics& d) {
 return std::abs(d.P_dir_actual_erg_s)+std::abs(d.LH_erg_s)+std::abs(d.DeltaLnu_erg_s)
 +std::abs(d.Lnu_eq_erg_s)+std::abs(d.Lgamma_erg_s)+std::abs(d.Lother_erg_s)+std::abs(d.L_out_fluid_inf_erg_s);
}
inline double Ulp(double x){return std::nextafter(x,INFINITY)-x;}
inline std::string Units(const std::string& n) {
 if(n=="Echem"||n=="mu_B"||n=="mu_n_actual")return "MeV";
 if(n=="sigma_e"||n=="sigma_mu")return "count/s";
 if(n=="Tinf")return "K";if(n=="Cstar")return "erg/K";return "erg/s";
}
inline std::string Quantity(const std::string& n) {
 if(n=="Echem")return "whole-star chemical free energy (not a power)";
 if(n=="Tinf")return "redshifted temperature";
 if(n=="Cstar")return "whole-star heat capacity";
 return n;
}
inline void Array(std::ostream& out,const BNV::CheckpointState& a) {
 out<<'['<<a[0]<<','<<a[1]<<','<<a[2]<<']';
}
// All rows, including passed channels, are durable before caller may refuse.
inline void Write(const std::filesystem::path& dir,const BNV::CheckpointOutput& q,
 const BNV::AcceptedStepRecord& step,const BNV::Rkf45Configuration& config,
 const std::vector<std::map<std::string,double>>& references={}) {
 std::filesystem::create_directories(dir);
 std::ofstream state(dir/"state.tsv"),diag(dir/"diagnostics.tsv"),meta(dir/"checkpoint.json"),reason(dir/"failure.txt");
 state<<std::setprecision(17);diag<<std::setprecision(17);meta<<std::setprecision(17);reason<<std::setprecision(17);
 reason<<"checkpoint="<<q.observation_index<<" time_s="<<q.t_observation_s<<" status="
 <<(q.status==BNV::CheckpointStatus::Qualified?"QUALIFIED":"NUMERICALLY_UNRESOLVED")<<'\n';
 state<<"channel\tunits\tO1\tO2\td_O\tD_O1\tF_O\tF_i\tU_O\tU_event\tpropagated_state_radius\tthreshold\tutilization\tpass\n";
 const char* names[]={"x","eta_e","eta_mu"};
 for(int i=0;i<3;++i){const auto& a=q.qualification;const double u=a.U_O[i]+a.U_event[i],limit=.20*a.F_i[i];
 const bool ok=a.d_O[i]<=a.D_O1[i]&&u<=limit;
 state<<names[i]<<'\t'<<(i?"MeV":"dimensionless")<<'\t'<<q.O1_state[i]<<'\t'<<q.O2_state[i]<<'\t'<<a.d_O[i]<<'\t'<<a.D_O1[i]<<'\t'<<a.F_O[i]<<'\t'<<a.F_i[i]<<'\t'<<a.U_O[i]<<'\t'<<a.U_event[i]<<'\t'<<q.reconstruction_error_contribution[i]<<'\t'<<limit<<'\t'<<(limit?u/limit:0)<<'\t'<<ok<<'\n';
 if(!ok)reason<<"FAIL state "<<names[i]<<": d_O <= D_O1 AND U_O+U_event <= 0.20 F_i; U="<<u<<" F="<<a.F_i[i]<<" utilization="<<(limit?u/limit:INFINITY)<<'\n';}
 const auto v1=Values(q.O1_diagnostics),v2=Values(q.O2_diagnostics);
 const double G=std::max({1.,Gross(q.O1_diagnostics),Gross(q.O2_diagnostics)});
 diag<<"channel\tquantity\tunits\tO1\tO2\tU\tF\tthreshold\tU_over_F\tutilization\td_O\t64ulp\twitness_or_roundoff_term\tpropagated_state_term\tG_erg_s\tpass";
 for(std::size_t i=0;i<references.size();++i)diag<<"\treference_"<<i;diag<<'\n';
 for(const auto& item:q.diagnostic_U_O){const auto& name=item.first;const double u=item.second,F=q.diagnostic_F_P.at(name),limit=.20*F;
 const double d=q.diagnostic_d_O.at(name),floor=64*Ulp(std::max(std::abs(v1.at(name)),std::abs(v2.at(name))));
 const auto it=q.propagated_state_uncertainty.find(name);const double propagated=it==q.propagated_state_uncertainty.end()?0:it->second;
 const bool ok=u<=limit;
 diag<<name<<'\t'<<Quantity(name)<<'\t'<<Units(name)<<'\t'<<v1.at(name)<<'\t'<<v2.at(name)<<'\t'<<u<<'\t'<<F<<'\t'<<limit<<'\t'<<u/F<<'\t'<<u/limit<<'\t'<<d<<'\t'<<floor<<'\t'<<2*std::max(d,floor)<<'\t'<<propagated<<'\t'<<G<<'\t'<<ok;
 for(const auto& ref:references)diag<<'\t'<<ref.at(name);diag<<'\n';
 if(!ok)reason<<"FAIL diagnostic "<<name<<" ("<<Quantity(name)<<", "<<Units(name)<<"): U <= 0.20 F; U="<<u<<" F="<<F<<" threshold="<<limit<<" U/F="<<u/F<<" utilization="<<u/limit<<" O1="<<v1.at(name)<<" O2="<<v2.at(name)<<" d="<<d<<" 64ulp="<<floor<<" propagated="<<propagated<<'\n';}
 meta<<"{\n\"checkpoint_index\":"<<q.observation_index<<",\n\"time_s\":"<<q.t_observation_s<<",\n\"method\":"<<std::quoted(q.provenance.method)<<",\n\"state\":";Array(meta,q.state);
 meta<<",\n\"O1_state\":";Array(meta,q.O1_state);meta<<",\n\"O2_state\":";Array(meta,q.O2_state);
 meta<<",\n\"main_step_ordinal\":"<<step.ordinal<<",\n\"main_left_time_s\":"<<step.t_left_s<<",\n\"main_right_time_s\":"<<step.t_right_s<<",\n\"main_left_state\":";Array(meta,step.y_left);meta<<",\n\"main_right_state\":";Array(meta,step.y_right);
 meta<<",\n\"main_local_rk_error\":";Array(meta,step.local_rk_error);
 meta<<",\n\"main_accepted_step_s\":"<<step.accepted_step_s<<",\n\"main_suggested_next_h_s\":"<<step.suggested_next_h_s<<",\n\"main_cumulative_rejected\":"<<step.cumulative_rejected;
 meta<<",\n\"main_rtol\":"<<config.relative<<",\n\"main_atol\":";Array(meta,config.absolute);
 meta<<",\n\"main_initial_step_s\":"<<config.initial_step_s<<",\n\"main_maximum_step_s\":";if(std::isfinite(config.maximum_step_s))meta<<config.maximum_step_s;else meta<<"null";
 meta<<",\n\"O1_rtol\":"<<q.provenance.O1.relative<<",\n\"O1_atol\":";Array(meta,q.provenance.O1.absolute);
 meta<<",\n\"O2_rtol\":"<<q.provenance.O2.relative<<",\n\"O2_atol\":";Array(meta,q.provenance.O2.absolute);
 meta<<",\n\"refinement_invocations\":"<<q.provenance.refinement_invocations<<",\n\"independent_method_invocations\":"<<q.provenance.independent_method_invocations;
 meta<<",\n\"normalization\":\"G=max(1 erg/s, gross power sum for O1 and O2); power F=max(1e-11 G, 1e-11 max(1,|O1|,|O2|),64ulp); nonpower F=max(1e-11 max(1,|O1|,|O2|),64ulp). Existing Lout branch is nonpower. Echem is MeV, not erg.\",\n\"criterion\":\"diagnostic U=2 max(|O2-O1|,64ulp)+propagated_state_uncertainty <= 0.20 F\",\n\"main_at_checkpoint\":\"No exact main state exists in strict interior; O2 is reconstructed candidate, main bracketing states are above\"\n}\n";
 for(auto* f:{&state,&diag,&meta,&reason}){f->flush();if(!*f)throw std::runtime_error("failure evidence write failed");}
}
} // namespace CheckpointEvidence
