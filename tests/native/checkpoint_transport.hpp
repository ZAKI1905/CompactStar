#pragma once
// Harness-local, same-build IPC. Preserve every field and binary64 bit pattern;
// never recompute diagnostics or sum ledger terms in workers.
#include <CompactStar/Physics/BNV/PassiveCheckpointOutput.hpp>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <type_traits>

namespace NativeQualification::Transport {
namespace B = CompactStar::Physics::BNV;
template<class A> void Fields(A& a,B::FrozenValidityResult& v) {
 a(v.fractional_depletion,v.max_utilization,v.DeltaN_over_N,v.limiting_quantity,
   v.drift_bound,v.threshold,v.utilization,v.valid);
}
template<class A> void Fields(A& a,B::BnvDiagnostics& v) {
 a(v.schema_id,v.run_card_identity,v.source_identity,v.domain_identity,v.revision_identity,
   v.partition_identity,v.product_fate_identity,v.actual_potential_provenance,v.finite_T_weighting_class,
   v.terminal_fate_branch_ids,v.terminal_fate_channel_ids,v.terminal_fate_branch_weights,
   v.t_s,v.B_count,v.Bdot_count_s,v.DeltaB_over_B0,v.S_count_s,v.t,v.t_error,
   v.sigma_count_s,v.eta_MeV,v.xi,v.R_count_s,v.eta_dot_from_sigma_MeV_s,v.eta_dot_from_beta_MeV_s,
   v.bSigma_residual_count_s,v.lift_residual_count_s,v.mu_B_inf_MeV,v.mu_n_actual_inf_MeV,v.g_actual_residual_MeV,
   v.Echem_MeV,v.Echem_dot_reaction_MeV_s,v.Echem_dot_source_MeV_s,v.Echem_dot_total_MeV_s,
   v.P_dir_eq_erg_s,v.P_dir_actual_erg_s,v.R18_residual_erg_s,v.Ra_Rb_residual_erg_s,v.Rb_Rc_residual_erg_s,
   v.Eesc_fluid_inf_MeV_s,v.Eesc_star_inf_MeV_s,v.EX_inf_MeV_s,v.L_out_fluid_inf_erg_s,v.L_esc_star_inf_erg_s,
   v.J_X_inf_erg_s,v.finite_T_omitted_floor_erg_s,v.LH_erg_s,v.DeltaLnu_erg_s,v.DeltaPbeta_erg_s,
   v.Lnu_eq_erg_s,v.Lnu_full_erg_s,v.Lgamma_erg_s,v.Lother_erg_s,v.Pnet_erg_s,
   v.Cstar_erg_K,v.Tinf_K,v.Tsurface_inf_K,v.DeltaTinf_K,v.DeltaTsurface_inf_K,v.DeltaLgamma_erg_s,v.DeltaU_th_erg,
   v.R20_residual_erg,v.N_R20_erg,v.R20_normalized,v.tau_relax_s,v.qss_evolution_time_s,v.qss_balance_ratio,
   v.regime_classification,v.DeltaN_over_N,v.frozen_drift_bound,v.frozen_threshold,
   v.sign_observable,v.sign_classification,v.sign_t0_s,v.sign_t1_s,v.sign_central,v.sign_lower_error,v.sign_upper_error,
   v.frozen,v.valid_through_sample);
}
template<class A> void Fields(A& a,B::BreakpointStepper::Event& v) {
 a(v.knot,v.direction,v.time,v.time_low,v.time_high,v.residual,v.time_ulp,v.state,v.rk_error);
}
template<class A> void Fields(A& a,B::Rk8pdConfiguration& v) {a(v.level,v.relative,v.absolute);}
template<class A> void Fields(A& a,B::ReconstructionQualification& v) {
 a(v.d_O,v.D_O1,v.F_O,v.F_i,v.U_O,v.U_event,v.discriminating_witness,v.self_qualified);
}
template<class A> void Fields(A& a,B::CheckpointProvenance& v) {
 a(v.method,v.left_accepted_ordinal,v.right_accepted_ordinal,v.t_left_s,v.t_right_s,v.t_observation_s,
   v.O1,v.O2,v.rk8pd_invocations,v.fallback_invocations,v.refinement_invocations,v.independent_method_invocations,
   v.context_identity_O1,v.context_identity_O2,v.O1_events,v.O2_events,v.refinement_events,v.independent_events,
   v.exact_zero_chemical_subspace);
}
template<class A> void Fields(A& a,B::CheckpointPerformance& v) {
 a(v.O1_context_wall_seconds,v.O1_context_cpu_seconds,v.O1_solve_wall_seconds,v.O1_solve_cpu_seconds,
   v.O2_context_wall_seconds,v.O2_context_cpu_seconds,v.O2_solve_wall_seconds,v.O2_solve_cpu_seconds,
   v.diagnostic_context_wall_seconds,v.diagnostic_context_cpu_seconds,
   v.diagnostic_evaluation_wall_seconds,v.diagnostic_evaluation_cpu_seconds);
}
template<class A> void Fields(A& a,B::CheckpointOutput& v) {
 a(v.observation_index,v.t_observation_s,v.source,v.status,v.state,v.O1_state,v.O2_state,
   v.qualification,v.provenance,v.performance,v.diagnostics,v.O1_diagnostics,v.O2_diagnostics,
   v.diagnostics_evaluated,v.diagnostic_self_qualified,v.diagnostic_d_O,v.diagnostic_U_O,
   v.diagnostic_F_P,v.propagated_state_uncertainty,v.reconstruction_error_contribution);
}
struct Writer {
 std::ostream& out;
 template<class T> void Put(const T& v) {
  if constexpr(std::is_arithmetic_v<T>||std::is_enum_v<T>)out.write(reinterpret_cast<const char*>(&v),sizeof v);
  else Fields(*this,const_cast<T&>(v)); // Fields only visits; Writer never mutates.
 }
 void Put(const std::string& v){Put(std::uint64_t(v.size()));out.write(v.data(),v.size());}
 template<class T,std::size_t N> void Put(const std::array<T,N>& v){for(const auto& x:v)Put(x);}
 template<class T> void Put(const std::vector<T>& v){Put(std::uint64_t(v.size()));for(const auto& x:v)Put(x);}
 template<class K,class V> void Put(const std::map<K,V>& v){Put(std::uint64_t(v.size()));for(const auto& x:v){Put(x.first);Put(x.second);}}
 template<class... T> void operator()(const T&... v){(Put(v),...);if(!out)throw std::runtime_error("checkpoint IPC write failed");}
};
struct Reader {
 std::istream& in;
 template<class T> void Get(T& v) {
  if constexpr(std::is_arithmetic_v<T>||std::is_enum_v<T>)in.read(reinterpret_cast<char*>(&v),sizeof v);
  else Fields(*this,v);
  if(!in)throw std::runtime_error("truncated checkpoint IPC");
 }
 std::size_t Count(){std::uint64_t n;Get(n);if(n>65536)throw std::runtime_error("checkpoint IPC length limit");return n;}
 void Get(std::string& v){v.resize(Count());in.read(v.data(),v.size());if(!in)throw std::runtime_error("truncated IPC string");}
 template<class T,std::size_t N> void Get(std::array<T,N>& v){for(auto& x:v)Get(x);}
 template<class T> void Get(std::vector<T>& v){v.resize(Count());for(auto& x:v)Get(x);}
 template<class K,class V> void Get(std::map<K,V>& v){auto n=Count();v.clear();for(std::size_t i=0;i<n;++i){K k;V x;Get(k);Get(x);if(!v.emplace(k,x).second)throw std::runtime_error("duplicate IPC key");}}
 template<class... T> void operator()(T&... v){(Get(v),...);}
};
inline void Write(const std::filesystem::path& path,const B::CheckpointOutput& q) {
 const auto temporary=path.string()+".partial";
 std::ofstream out(temporary,std::ios::binary);Writer w{out};
 w(std::string("CompactStar same-build checkpoint IPC v1"),q);out.close();
 if(!out)throw std::runtime_error("checkpoint IPC close failed");
 std::filesystem::rename(temporary,path);
}
inline B::CheckpointOutput Read(const std::filesystem::path& path) {
 if(std::filesystem::file_size(path)>4*1024*1024)throw std::runtime_error("oversize checkpoint IPC");
 std::ifstream in(path,std::ios::binary);Reader r{in};std::string magic;B::CheckpointOutput q;r(magic);
 if(magic!="CompactStar same-build checkpoint IPC v1")throw std::runtime_error("checkpoint IPC schema mismatch");
 r(q);if(in.peek()!=std::char_traits<char>::eof())throw std::runtime_error("trailing checkpoint IPC bytes");return q;
}
} // namespace NativeQualification::Transport
