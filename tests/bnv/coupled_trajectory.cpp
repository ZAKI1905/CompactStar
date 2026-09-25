#include "campaign_fixture.hpp"
#include <CompactStar/Physics/Evolution/Observers/IObserver.hpp>
#include <CompactStar/Physics/Evolution/StatePacking.hpp>
#include <CompactStar/Physics/Rotochemical/ScaledRKF45.hpp>
#include <gsl/gsl_errno.h>
#include <iomanip>

namespace Campaign=Phase6A1Campaign;
namespace BNV=CompactStar::Physics::BNV;
namespace EV=CompactStar::Physics::Evolution;
namespace P=CompactStar::Physics;
using Tag=P::State::StateTag;

struct Capture final:EV::Observers::IObserver
{
    static const std::vector<std::string>& FrozenNames()
    {
        static const std::vector<std::string> names{
          "t_n","t_e","t_mu","Z_row_npe","Z_row_npmu","Cstar","Ltilde_Me","Ltilde_Mmu",
          "mu_B","mu_n_profile","P0_average","P1_average","P2_average","N_n","N_e","N_mu",
          "species_support","metric_structure","radius","surface_gravity","envelope","Tsurface_inf"};
        return names;
    }
    static std::string Join(const std::vector<std::string>& values)
    {
        std::ostringstream out;for(std::size_t i=0;i<values.size();++i){if(i)out<<';';out<<values[i];}return out.str();
    }
    static std::string Join(const std::vector<double>& values)
    {
        std::ostringstream out;out<<std::setprecision(17);for(std::size_t i=0;i<values.size();++i){if(i)out<<';';out<<values[i];}return out.str();
    }
    std::ofstream out;
    const BNV::ControlledBnvSecularDriver& driver;
    bool started=false;
    std::size_t rows=0;
    explicit Capture(const std::filesystem::path& path,const BNV::ControlledBnvSecularDriver& d):out(path),driver(d)
    {
        if(!out)throw std::runtime_error("cannot create trajectory output");
        out<<"t_s\tx_state\tB_count\tBdot_count_s\tDeltaB_over_B0\tS_n_count_s\tS_e_count_s\tS_mu_count_s\tt_n\tt_e\tt_mu\tt_error_n\tt_error_e\tt_error_mu\tsigma_e_count_s\tsigma_mu_count_s\teta_e_MeV\teta_mu_MeV\txi_e\txi_mu\tR_e_count_s\tR_mu_count_s\teta_dot_from_sigma_e_MeV_s\teta_dot_from_sigma_mu_MeV_s\teta_dot_from_beta_e_MeV_s\teta_dot_from_beta_mu_MeV_s\tmu_B_inf_MeV\tmu_n_actual_inf_MeV\tEchem_MeV\tEchem_dot_reaction_MeV_s\tEchem_dot_source_MeV_s\tEchem_dot_total_MeV_s\tP_dir_eq_erg_s\tP_dir_actual_erg_s\tLH_erg_s\tDeltaLnu_erg_s\tDeltaPbeta_erg_s\tLnu_eq_erg_s\tLnu_full_erg_s\tL_out_fluid_inf_erg_s\tL_esc_star_inf_erg_s\tJ_X_inf_erg_s\tLgamma_erg_s\tLother_erg_s\tPnet_erg_s\tCstar_erg_K\tTinf_K\tTsurface_inf_K\tR18_residual_erg_s\tRa_Rb_residual_erg_s\tRb_Rc_residual_erg_s\tfinite_T_omitted_floor_erg_s\tbSigma_residual_count_s\tlift_residual_count_s\tDeltaN_n_over_N_n\tDeltaN_e_over_N_e\tDeltaN_mu_over_N_mu\tmax_frozen_utilization\tfrozen_limiting_quantity\tvalid_through_sample\trun_card_identity\tpartition_identity\tsource_identity\tdomain_identity\trevision_identity\tproduct_fate_identity\tactual_potential_provenance\tfinite_T_weighting_class"
           <<"\tg_actual_residual_MeV\tEesc_fluid_inf_MeV_s\tEesc_star_inf_MeV_s\tEX_inf_MeV_s"
           <<"\tDeltaTinf_K\tDeltaTsurface_inf_K\tDeltaLgamma_erg_s\tDeltaU_th_erg"
           <<"\tR20_residual_erg\tN_R20_erg\tR20_normalized"
           <<"\ttau_relax_e_s\ttau_relax_mu_s\tqss_evolution_time_e_s\tqss_evolution_time_mu_s\tqss_balance_e\tqss_balance_mu"
           <<"\tregime_e\tregime_mu\tsign_observable\tsign_classification\tsign_t0_s\tsign_t1_s\tsign_central\tsign_lower_error\tsign_upper_error"
           <<"\tterminal_fate_branch_ids\tterminal_fate_channel_ids\tterminal_fate_branch_weights";
        for(const auto& name:FrozenNames())out<<'\t'<<name<<"_drift_bound\t"<<name<<"_threshold\t"<<name<<"_utilization";
        out<<'\n';
        out<<std::setprecision(17);
    }
    void Save(double t,const EV::StateVector& state,const EV::DriverContext& context)
    {
        const auto value=driver.Evaluate(t,state,context);const auto& d=value.diagnostics;
        const std::array<double,58> numbers{{
          d.t_s,state.GetThermal().LnTinfOverTref(),d.B_count,d.Bdot_count_s,d.DeltaB_over_B0,
          d.S_count_s[0],d.S_count_s[1],d.S_count_s[2],d.t[0],d.t[1],d.t[2],d.t_error[0],d.t_error[1],d.t_error[2],
          d.sigma_count_s[0],d.sigma_count_s[1],d.eta_MeV[0],d.eta_MeV[1],d.xi[0],d.xi[1],d.R_count_s[0],d.R_count_s[1],
          d.eta_dot_from_sigma_MeV_s[0],d.eta_dot_from_sigma_MeV_s[1],d.eta_dot_from_beta_MeV_s[0],d.eta_dot_from_beta_MeV_s[1],
          d.mu_B_inf_MeV,d.mu_n_actual_inf_MeV,d.Echem_MeV,d.Echem_dot_reaction_MeV_s,d.Echem_dot_source_MeV_s,d.Echem_dot_total_MeV_s,
          d.P_dir_eq_erg_s,d.P_dir_actual_erg_s,d.LH_erg_s,d.DeltaLnu_erg_s,d.DeltaPbeta_erg_s,d.Lnu_eq_erg_s,d.Lnu_full_erg_s,
          d.L_out_fluid_inf_erg_s,d.L_esc_star_inf_erg_s,d.J_X_inf_erg_s,d.Lgamma_erg_s,d.Lother_erg_s,d.Pnet_erg_s,d.Cstar_erg_K,d.Tinf_K,d.Tsurface_inf_K,
          d.R18_residual_erg_s,d.Ra_Rb_residual_erg_s,d.Rb_Rc_residual_erg_s,d.finite_T_omitted_floor_erg_s,d.bSigma_residual_count_s,d.lift_residual_count_s,
          d.DeltaN_over_N[0],d.DeltaN_over_N[1],d.DeltaN_over_N[2],d.frozen.max_utilization}};
        for(double x:numbers)if(!std::isfinite(x))throw std::runtime_error("nonfinite trajectory diagnostic");
        for(std::size_t i=0;i<numbers.size();++i){if(i)out<<'\t';out<<numbers[i];}
        out<<'\t'<<d.frozen.limiting_quantity<<'\t'<<(d.valid_through_sample?1:0)
           <<'\t'<<d.run_card_identity<<'\t'<<d.partition_identity<<'\t'<<d.source_identity
           <<'\t'<<d.domain_identity<<'\t'<<d.revision_identity<<'\t'<<d.product_fate_identity
           <<'\t'<<d.actual_potential_provenance<<'\t'<<d.finite_T_weighting_class
           <<'\t'<<d.g_actual_residual_MeV<<'\t'<<d.Eesc_fluid_inf_MeV_s<<'\t'<<d.Eesc_star_inf_MeV_s<<'\t'<<d.EX_inf_MeV_s
           <<'\t'<<d.DeltaTinf_K<<'\t'<<d.DeltaTsurface_inf_K<<'\t'<<d.DeltaLgamma_erg_s<<'\t'<<d.DeltaU_th_erg
           <<'\t'<<d.R20_residual_erg<<'\t'<<d.N_R20_erg<<'\t'<<d.R20_normalized
           <<'\t'<<d.tau_relax_s[0]<<'\t'<<d.tau_relax_s[1]<<'\t'<<d.qss_evolution_time_s[0]<<'\t'<<d.qss_evolution_time_s[1]
           <<'\t'<<d.qss_balance_ratio[0]<<'\t'<<d.qss_balance_ratio[1]<<'\t'<<d.regime_classification[0]<<'\t'<<d.regime_classification[1]
           <<'\t'<<d.sign_observable<<'\t'<<d.sign_classification<<'\t'<<d.sign_t0_s<<'\t'<<d.sign_t1_s<<'\t'<<d.sign_central<<'\t'<<d.sign_lower_error<<'\t'<<d.sign_upper_error
           <<'\t'<<Join(d.terminal_fate_branch_ids)<<'\t'<<Join(d.terminal_fate_channel_ids)<<'\t'<<Join(d.terminal_fate_branch_weights);
        for(const auto& name:FrozenNames())out<<'\t'<<d.frozen_drift_bound.at(name)<<'\t'<<d.frozen_threshold.at(name)<<'\t'<<d.frozen.utilization.at(name);
        out<<'\n';
        if(!out)throw std::runtime_error("trajectory serialization failure");++rows;
    }
    void OnStart(const EV::Observers::RunInfo&,const EV::StateVector& state,const EV::DriverContext& context)override
    {if(!started){Save(0,state,context);started=true;}}
    void OnSample(const EV::Observers::SampleInfo& sample,const EV::StateVector& state,const EV::DriverContext& context)override
    {Save(sample.t,state,context);}
};

struct RunSummary
{
    std::size_t accepted=0,rejected=0,rhs=0,rows=0;
    double min_step=std::numeric_limits<double>::infinity(),max_step=0;
    std::vector<RC::StepOutput> outputs;
};

RunSummary Run(const std::shared_ptr<const BNV::FrozenControlledBnvRunContext>& context,
    BNV::ControlledBnvSecularDriver::Mode mode,const Campaign::RunCard& card,
    const std::filesystem::path& path,const RC::ComponentTolerances& tolerances)
{
    RunState state(context->OrdinaryContext());state.thermal.SetTinf(1e8);RC::ChemicalImbalanceState(0,0).Store(state.chem);
    auto y=state.Pack();auto driver=std::make_shared<BNV::ControlledBnvSecularDriver>(context,mode);
    EV::EvolutionSystem system(state.ctx,state.state,state.rhs,state.layout,{driver});auto capture=std::make_shared<Capture>(path,*driver);system.AddObserver(capture);
    RC::ScaledRKF45 solver(system,state.layout,context->OrdinaryContext(),tolerances);
    std::vector<double> checkpoints;checkpoints.reserve(card.checkpoints-1);const double end=card.duration_year*Year;
    for(std::size_t i=1;i<card.checkpoints;++i)checkpoints.push_back(end*static_cast<double>(i)/static_cast<double>(card.checkpoints-1));
    RunSummary total;double start=0;
    for(std::size_t begin=0;begin<checkpoints.size();begin+=999)
    {
        const std::size_t stop=std::min(checkpoints.size(),begin+999);std::vector<double> batch(checkpoints.begin()+begin,checkpoints.begin()+stop);
        RC::IntegrationStatistics stats;solver.Integrate(start,batch,y.data(),stats,true);start=batch.back();
        total.accepted+=stats.accepted_steps;total.rejected+=stats.rejected_steps;total.rhs+=stats.rhs_evaluations;
        total.min_step=std::min(total.min_step,stats.minimum_step_s);total.max_step=std::max(total.max_step,stats.maximum_step_s);
        const std::size_t accepted_offset=total.accepted-stats.accepted_steps;
        const std::size_t rejected_offset=total.rejected-stats.rejected_steps;
        for(const auto& sample:stats.outputs)total.outputs.push_back({sample.time_s,
          accepted_offset+sample.accepted,rejected_offset+sample.rejected,sample.last_step_s});
    }
    total.rows=capture->rows;if(total.rows!=card.checkpoints)throw std::runtime_error("trajectory checkpoint count changed");
    std::ofstream steps(path.string()+".steps");steps<<std::setprecision(17)<<"accepted\trejected\trhs\tminimum_step_s\tmaximum_step_s\trows\n"<<total.accepted<<'\t'<<total.rejected<<'\t'<<total.rhs<<'\t'<<total.min_step<<'\t'<<total.max_step<<'\t'<<total.rows<<"\nt_s\tcumulative_accepted\tcumulative_rejected\tlast_step_s\n";
    for(const auto& sample:total.outputs)steps<<sample.time_s<<'\t'<<sample.accepted<<'\t'<<sample.rejected<<'\t'<<sample.last_step_s<<'\n';
    return total;
}

std::shared_ptr<const BNV::IDirectBnvEnergyPartition> Partition(const Campaign::RunCard& card,
    const ControlledFixture& fixture,const std::shared_ptr<const BNV::BnvDependencyToken>& token,const std::string& fate)
{
    if(card.partition=="P0")return std::make_shared<const Phase6A1Test::P0Partition>("P0-controlled-mathematical-v1",fate,1,token);
    if(card.partition=="P1")return std::make_shared<const Phase6A1Test::P1IntegratedUniformSeaPartition>("P1-controlled-relativistic-uniform-sea-v1",fate,fixture.central,fixture.provider,token);
    if(card.partition=="P2")return std::make_shared<const Phase6A1Test::P2Partition>("P2-controlled-full-retention-v1",fate,1,token);
    throw std::runtime_error("unknown run-card partition");
}

#ifndef PHASE6A1_BA12R_ULTRA_SUPPORT_ONLY
int main(int argc,char** argv)
{
 try
 {
    gsl_set_error_handler_off();require(argc==10,"profile certificate thermal work entry frozen-certificate coefficients output pretrajectory-record");
    std::cout<<std::setprecision(17)<<std::unitbuf;
    const std::filesystem::path profile=argv[1],certificate=argv[2],thermal_path=argv[3],work=argv[4],entry=argv[5],frozen=argv[6],coefficients=argv[7],output=argv[8],pretrajectory=argv[9];
    require(!std::filesystem::exists(work),"fresh campaign assembly directory required");require(!std::filesystem::exists(output),"fresh trajectory output directory required");
    {std::ifstream in(pretrajectory);std::string text((std::istreambuf_iterator<char>(in)),{});require(text.find("PRETRAJECTORY PASS")!=std::string::npos&&text.find("no BNV trajectory had been generated")!=std::string::npos,"committed pretrajectory authorization missing");}
    std::filesystem::create_directories(output);auto fixture=Fixture(profile,certificate.string(),work/"owning-star",80000);auto channels=Channels(fixture);
    EOS::CompOSE_Thermo::Options options;options.Tmin_for_derivative_MeV=0;options.clamp_to_domain=false;
    auto thermal=std::make_shared<const RC::FrozenThermalSource>(thermal_path,options,"controlled mathematical fixed-background free-gas entropy; qualified radial80000");
    auto tangent=Campaign::Tangent(fixture,profile,work);auto bnv_token=std::make_shared<BNV::BnvDependencyToken>();auto monitor=Phase6A1Test::LoadFrozenMonitor(frozen,tangent,bnv_token);
    auto run_token=std::make_shared<RC::RunDependencyToken>();auto spin=std::make_shared<const BNV::StaticZeroSpinHistory>(run_token);
    auto ordinary=Context(fixture,channels,thermal,spin,Campaign::Qualification(profile,certificate,entry),run_token);
    const double mu_B=Campaign::Column(coefficients,"mu_B_inf");
    const std::string potential="frozen Structure-1 B0 equilibrium mu_B plus governed moving-reference actual-potential correction";
    for(const auto& card:Campaign::RunCards())
    {
        const double Bdot=card.fractional_drive_per_year*tangent->B0Count()/Year;
        const std::string fate_id="phase6a1-"+card.partition+"-terminal-fate-v1";
        const auto terminal=card.partition=="P2"?BNV::TerminalProductFate::BoundInert:BNV::TerminalProductFate::PromptEscape;
        BNV::ProductFateLedger fate(fate_id,{{"controlled-terminal-product",terminal,1,"abstract-neutron-disappearance"}});
        auto partition=Partition(card,fixture,bnv_token,fate_id);
        auto history=std::make_shared<const Phase6A1Test::UniformProperNeutronSinkHistory>(fixture.central,tangent->B0Count(),Bdot,tangent->DomainIdentity(),tangent->StarIdentity(),tangent->SequenceStateIdentity(),partition->Identity(),fate_id,bnv_token);
        auto context=std::make_shared<const BNV::FrozenControlledBnvRunContext>(ordinary,tangent,history,partition,fate,monitor,channels,mu_B,potential,card.identity);
        const auto mode=card.reaction_free?BNV::ControlledBnvSecularDriver::Mode::ReactionFreeControl:BNV::ControlledBnvSecularDriver::Mode::Coupled;
        auto zero_history=std::make_shared<const Campaign::ExactZeroHistory>(tangent->B0Count(),tangent->DomainIdentity(),tangent->StarIdentity(),tangent->SequenceStateIdentity(),bnv_token);
        auto zero_partition=std::make_shared<const Campaign::ExactZeroPartition>(bnv_token);BNV::ProductFateLedger zero_fate("zero-control-fate",{{"none",BNV::TerminalProductFate::BoundInert,1,"zero-control"}});
        auto control=std::make_shared<const BNV::FrozenControlledBnvRunContext>(ordinary,tangent,zero_history,zero_partition,zero_fate,monitor,channels,mu_B,potential,card.identity+"-MATCHED-CONTROL");
        for(bool refined:{false,true})
        {
            const std::string level=refined?"refined":"baseline";
            const auto tolerance=refined?RC::ComponentTolerances::Refined():RC::ComponentTolerances{};
            const auto a=Run(context,mode,card,output/(card.identity+"."+level+".tsv"),tolerance);
            const auto b=Run(control,mode,card,output/(card.identity+".control."+level+".tsv"),tolerance);
            std::cout<<"TRAJECTORY PASS "<<card.identity<<' '<<level<<" source_rows "<<a.rows<<" control_rows "<<b.rows<<" accepted "<<a.accepted<<" control_accepted "<<b.accepted<<'\n';
        }
    }
    std::cout<<"BA11_RAW_TRAJECTORIES PASS four predeclared cards and matched controls\n";
    return 0;
 }
 catch(const std::exception& e){std::cerr<<"STOP "<<e.what()<<'\n';return 1;}
}
#endif
