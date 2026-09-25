#pragma once

#include <CompactStar/Physics/BNV/FrozenBnvValidity.hpp>
#include <array>
#include <map>
#include <string>
#include <vector>

namespace CompactStar::Physics::BNV
{

struct BnvDiagnostics
{
    std::string schema_id="compactstar.phase6a1.controlled-bnv-candidate.v1";
    std::string run_card_identity,source_identity,domain_identity,revision_identity;
    std::string partition_identity,product_fate_identity,actual_potential_provenance;
    std::string finite_T_weighting_class;
    std::vector<std::string> terminal_fate_branch_ids;
    std::vector<std::string> terminal_fate_channel_ids;
    std::vector<double> terminal_fate_branch_weights;
    double t_s=0,B_count=0,Bdot_count_s=0,DeltaB_over_B0=0;
    std::array<double,3> S_count_s{},t{},t_error{};
    std::array<double,2> sigma_count_s{},eta_MeV{},xi{},R_count_s{};
    std::array<double,2> eta_dot_from_sigma_MeV_s{},eta_dot_from_beta_MeV_s{};
    double bSigma_residual_count_s=0,lift_residual_count_s=0;
    double mu_B_inf_MeV=0,mu_n_actual_inf_MeV=0,g_actual_residual_MeV=0;
    double Echem_MeV=0,Echem_dot_reaction_MeV_s=0,Echem_dot_source_MeV_s=0,Echem_dot_total_MeV_s=0;
    double P_dir_eq_erg_s=0,P_dir_actual_erg_s=0,R18_residual_erg_s=0;
    double Ra_Rb_residual_erg_s=0,Rb_Rc_residual_erg_s=0;
    double Eesc_fluid_inf_MeV_s=0,Eesc_star_inf_MeV_s=0,EX_inf_MeV_s=0;
    double L_out_fluid_inf_erg_s=0,L_esc_star_inf_erg_s=0,J_X_inf_erg_s=0;
    double finite_T_omitted_floor_erg_s=0;
    double LH_erg_s=0,DeltaLnu_erg_s=0,DeltaPbeta_erg_s=0;
    double Lnu_eq_erg_s=0,Lnu_full_erg_s=0,Lgamma_erg_s=0,Lother_erg_s=0,Pnet_erg_s=0;
    double Cstar_erg_K=0,Tinf_K=0,Tsurface_inf_K=0;
    double DeltaTinf_K=0,DeltaTsurface_inf_K=0,DeltaLgamma_erg_s=0,DeltaU_th_erg=0;
    double R20_residual_erg=0,N_R20_erg=0,R20_normalized=0;
    std::array<double,2> tau_relax_s{},qss_evolution_time_s{},qss_balance_ratio{};
    std::array<std::string,2> regime_classification{{"NOT_CLASSIFIED","NOT_CLASSIFIED"}};
    std::array<double,3> DeltaN_over_N{};
    std::map<std::string,double> frozen_drift_bound,frozen_threshold;
    std::string sign_observable="NOT_CLASSIFIED",sign_classification="SIGN_UNRESOLVED";
    double sign_t0_s=0,sign_t1_s=0,sign_central=0,sign_lower_error=0,sign_upper_error=0;
    FrozenValidityResult frozen;
    bool valid_through_sample=false;
};

} // namespace CompactStar::Physics::BNV

