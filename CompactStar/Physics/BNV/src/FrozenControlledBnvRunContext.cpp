#include <CompactStar/Physics/BNV/FrozenControlledBnvRunContext.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace CompactStar::Physics::BNV
{
namespace RC=Rotochemical;

void FrozenControlledBnvRunContext::Need(bool ok,const char* message)
{if(!ok)throw std::runtime_error(message);}

FrozenControlledBnvRunContext::FrozenControlledBnvRunContext(
  std::shared_ptr<const RC::FrozenRotochemicalRunContext> ordinary,
  std::shared_ptr<const EquilibriumBaryonTangent> tangent,
  std::shared_ptr<const OrdinaryMatterBnvHistory> history,
  std::shared_ptr<const IDirectBnvEnergyPartition> partition,
  ProductFateLedger fate,
  std::shared_ptr<const FrozenBnvValidityMonitor> validity,
  std::shared_ptr<const RC::GlobalUrcaChannelCoefficient> channels,
  double mu_B,std::string potential_provenance,std::string run_card)
  :ordinary_(std::move(ordinary)),tangent_(std::move(tangent)),history_(std::move(history)),
   partition_(std::move(partition)),fate_(std::move(fate)),validity_(std::move(validity)),channels_(std::move(channels)),
   spin_(ordinary_?ordinary_->SpinOwner():nullptr),mu_B_inf_MeV_(mu_B),potential_provenance_(std::move(potential_provenance)),
   run_card_identity_(std::move(run_card))
{
    Need(ordinary_&&tangent_&&history_&&partition_&&validity_&&channels_&&spin_,"missing controlled BNV dependency");
    const auto purpose=ordinary_->Qualification().purpose;
    if(purpose==RC::RunPurpose::AnalyticControl)
    {
        Need(dynamic_cast<const StaticZeroSpinHistory*>(spin_.get())!=nullptr,"controlled BNV requires StaticZeroSpinHistory");
        Need(spin_->Identity()=="static zero spin Omega=0 rad s^-1 OmegaDot=0 rad s^-2; Phase-6A-1 controlled abstract BNV","wrong zero-spin identity");
    }
    else
    {
        Need(purpose==RC::RunPurpose::ControlledTrajectory,"unsupported controlled BNV run purpose");
        Need(dynamic_cast<const RC::PrescribedDipoleHistory*>(spin_.get())!=nullptr,
             "spin-on regression requires governed PrescribedDipoleHistory");
        spin_on_zero_source_regression_=true;
    }
    const auto& selection=channels_->Selection();
    Need(selection.Enabled(RC::UrcaProcess::Me)&&selection.Enabled(RC::UrcaProcess::Mmu)&&
         !selection.Enabled(RC::UrcaProcess::De)&&!selection.Enabled(RC::UrcaProcess::Dmu),"controlled BNV process selection mismatch");
    Need(channels_->MetricIdentity()=="qualified Structure-1 radial80000 canonical nu/lambda","controlled BNV metric identity mismatch");
    Need(channels_->Entry(RC::UrcaProcess::Me).normalization_identity==
      "predeclared mathematical benchmark SMe=1e-51 erg cm^-3 s^-1 K^-8","controlled BNV Me normalization mismatch");
    Need(channels_->Entry(RC::UrcaProcess::Mmu).normalization_identity==
      "predeclared mathematical benchmark SMmu=2e-51 erg cm^-3 s^-1 K^-8","controlled BNV Mmu normalization mismatch");
    for(auto p:RC::UrcaProcesses)
    {
        Need(ordinary_->Ltilde(p)==channels_->LuminosityCoefficient(p),"foreign controlled BNV channel owner");
        Need(ordinary_->ChannelEntry(p).normalization_identity==channels_->Entry(p).normalization_identity,
             "controlled BNV normalization owner mismatch");
    }
    Need(std::isfinite(mu_B_inf_MeV_)&&mu_B_inf_MeV_>0&&!potential_provenance_.empty()&&!run_card_identity_.empty(),
         "invalid controlled BNV reference potential");
    Need(validity_->Certificate()->B0Count()==tangent_->B0Count(),"tangent/certificate B0 mismatch");
    Need(validity_->Certificate()->StarIdentity()==tangent_->StarIdentity(),"tangent/certificate star mismatch");
    Need(validity_->Certificate()->DomainIdentity()==tangent_->DomainIdentity(),"tangent/certificate domain mismatch");
    history_identity_=history_->Identity();partition_identity_=partition_->Identity();spin_identity_=spin_->Identity();
    const auto entry=history_->Sample(0);entry.Validate();
    if(spin_on_zero_source_regression_)
    {
        Need(entry.Bdot_count_s==0&&entry.source_count_s==std::array<double,3>{0,0,0}&&
             entry.proton_source_count_s==0&&entry.events.empty()&&entry.B_count==tangent_->B0Count(),
             "spin-on regression permits only exact zero BNV source");
    }
    Need(entry.B_count==tangent_->B0Count(),"history/tangent initial B mismatch");
    Need(entry.star_identity==tangent_->StarIdentity()&&entry.domain_identity==tangent_->DomainIdentity()&&
         entry.sequence_state_identity==tangent_->SequenceStateIdentity(),"history/tangent identity mismatch");
    Need(std::all_of(entry.events.begin(),entry.events.end(),[&](const auto& e){return e.partition_identity==partition_identity_&&e.product_fate_identity==fate_.Identity();}),
         "history partition/fate mismatch");
    RequireFullCurrent();
}

void FrozenControlledBnvRunContext::RequireCheapCurrent() const
{
    ordinary_->RequireCheapCurrent();history_->RequireCurrent();partition_->RequireCurrent();
    // The immutable channel coefficients are already bound into ordinary_, whose
    // cheap currentness check covers their chemical lifetime and source owners.
    // Full channel byte-currentness remains a construction/integration-boundary
    // gate through RequireFullCurrent(); do not reread the response artifact for
    // every RKF45 trial evaluation.
    tangent_->RequireCheapCurrent();validity_->Certificate()->RequireCurrent();spin_->RequireCurrent();
    Need(history_->Identity()==history_identity_&&partition_->Identity()==partition_identity_&&spin_->Identity()==spin_identity_,
         "changed controlled BNV identity");
}
void FrozenControlledBnvRunContext::RequireFullCurrent() const
{ordinary_->RequireFullCurrent();tangent_->RequireCurrent();RequireCheapCurrent();}

ControlledBnvEvaluation FrozenControlledBnvRunContext::Evaluate(
    double epoch,const Evolution::StateVector& state,const Evolution::DriverContext& ctx) const
{return EvaluateImpl(epoch,state,ctx,false);}

ControlledBnvEvaluation FrozenControlledBnvRunContext::EvaluateReactionFree(
    double epoch,const Evolution::StateVector& state,const Evolution::DriverContext& ctx) const
{return EvaluateImpl(epoch,state,ctx,true);}

ControlledBnvEvaluation FrozenControlledBnvRunContext::EvaluateImpl(
    double epoch,const Evolution::StateVector& state,const Evolution::DriverContext& ctx,bool reaction_free) const
{
    RequireCheapCurrent();
    const auto tangent=tangent_->SnapshotCheap();
    ControlledBnvEvaluation out;
    const auto source=history_->Sample(epoch);source.Validate();
    Need(source.source_identity==history_identity_,"source/history identity mismatch");
    if(spin_on_zero_source_regression_)
    {
        Need(source.Bdot_count_s==0&&source.source_count_s==std::array<double,3>{0,0,0}&&
             source.proton_source_count_s==0&&source.events.empty()&&source.B_count==tangent_->B0Count(),
             "spin-on regression source became nonzero");
    }
    out.diagnostics.frozen=validity_->RequireValid(source.B_count); // trial-state/time gate before any RHS publication
    out.moving=MovingReferenceSource::Project(source,tangent);
    if(reaction_free)
    {
        const double T=state.GetThermal().Tinf();
        if(!reaction_free_reference_)
        {reaction_free_reference_=ordinary_->Evaluate(epoch,state,ctx,spin_.get());reaction_free_reference_Tinf_K_=T;}
        Need(T==reaction_free_reference_Tinf_K_,"reaction-free control thermal state changed");
        out.ordinary=*reaction_free_reference_;out.ordinary.reaction={};out.ordinary.beta={};
        out.ordinary.eta_dot_MeV_s={};out.ordinary.x_dot_s=0;out.ordinary.equilibrium_x_dot_s=0;
    }
    else out.ordinary=ordinary_->Evaluate(epoch,state,ctx,spin_.get());
    const auto eta=RC::ChemicalImbalanceState::Read(state.GetChem());
    out.direct.potential=BnvDirectEnergyLedger::ActualPotential(mu_B_inf_MeV_,eta,tangent,potential_provenance_);
    const auto energies=partition_->Evaluate(source,out.direct.potential);
    out.direct=BnvDirectEnergyLedger::Evaluate(source,out.direct.potential,energies,fate_);
    const double T=state.GetThermal().Tinf();
    const double ee=eta.InfinityMeV(RC::BetaChannel::Npe),em=eta.InfinityMeV(RC::BetaChannel::NpMu);
    for(std::size_t i=0;i<2;++i)
    {
        out.eta_dot_MeV_s[i]=out.ordinary.eta_dot_MeV_s[i];
        for(std::size_t j=0;j<2;++j)out.eta_dot_MeV_s[i]-=ordinary_->Z(i==0?RC::BetaChannel::Npe:RC::BetaChannel::NpMu,
          j==0?RC::BetaChannel::Npe:RC::BetaChannel::NpMu)*out.moving.sigma_count_s[j];
    }
    out.x_dot_s=out.ordinary.x_dot_s+out.direct.power_actual_erg_s/(T*out.ordinary.Cstar_erg_K);
    auto& d=out.diagnostics;
    d.run_card_identity=run_card_identity_;d.source_identity=source.source_identity;d.domain_identity=source.domain_identity;
    d.revision_identity=source.revision_identity;d.partition_identity=partition_identity_;d.product_fate_identity=fate_.Identity();
    for(const auto& branch:fate_.Branches()){d.terminal_fate_branch_ids.push_back(branch.terminal_id);d.terminal_fate_channel_ids.push_back(branch.channel_id);d.terminal_fate_branch_weights.push_back(branch.weight);}
    d.actual_potential_provenance=potential_provenance_;d.finite_T_weighting_class=partition_->FiniteTemperatureWeightingClass();
    d.t_s=epoch;d.B_count=source.B_count;d.Bdot_count_s=source.Bdot_count_s;d.DeltaB_over_B0=(source.B_count-tangent.B0_count)/tangent.B0_count;
    d.S_count_s=source.source_count_s;d.t=tangent.closed;d.t_error=tangent.numerical_error;d.sigma_count_s=out.moving.sigma_count_s;
    d.bSigma_residual_count_s=out.moving.baryon_residual_count_s;d.lift_residual_count_s=out.moving.lift_residual_count_s;
    d.eta_MeV={ee,em};d.xi={eta.Xi(RC::BetaChannel::Npe,T),eta.Xi(RC::BetaChannel::NpMu,T)};
    d.R_count_s=out.ordinary.reaction.rate_count_s;
    for(std::size_t i=0;i<2;++i)
    {
        const auto row=i==0?RC::BetaChannel::Npe:RC::BetaChannel::NpMu;
        d.eta_dot_from_sigma_MeV_s[i]=-ordinary_->Z(row,RC::BetaChannel::Npe)*out.moving.sigma_count_s[0]-ordinary_->Z(row,RC::BetaChannel::NpMu)*out.moving.sigma_count_s[1];
        d.eta_dot_from_beta_MeV_s[i]=-ordinary_->Z(row,RC::BetaChannel::Npe)*d.R_count_s[0]-ordinary_->Z(row,RC::BetaChannel::NpMu)*d.R_count_s[1];
    }
    const double z00=ordinary_->Z(RC::BetaChannel::Npe,RC::BetaChannel::Npe),z01=ordinary_->Z(RC::BetaChannel::Npe,RC::BetaChannel::NpMu);
    const double z10=ordinary_->Z(RC::BetaChannel::NpMu,RC::BetaChannel::Npe),z11=ordinary_->Z(RC::BetaChannel::NpMu,RC::BetaChannel::NpMu);
    const double det=z00*z11-z01*z10;Need(det>0&&std::isfinite(det),"invalid frozen Z inverse");
    d.Echem_MeV=.5*(z11*ee*ee-(z01+z10)*ee*em+z00*em*em)/det;
    d.Echem_dot_reaction_MeV_s=-ee*d.R_count_s[0]-em*d.R_count_s[1];
    d.Echem_dot_source_MeV_s=-ee*d.sigma_count_s[0]-em*d.sigma_count_s[1];
    d.Echem_dot_total_MeV_s=d.Echem_dot_reaction_MeV_s+d.Echem_dot_source_MeV_s;
    d.mu_B_inf_MeV=mu_B_inf_MeV_;d.mu_n_actual_inf_MeV=out.direct.potential.mu_n_actual_inf_MeV;
    d.g_actual_residual_MeV=out.direct.potential.reconstruction_residual_MeV;
    d.P_dir_eq_erg_s=out.direct.power_eq_erg_s;d.P_dir_actual_erg_s=out.direct.power_actual_erg_s;
    d.R18_residual_erg_s=d.P_dir_actual_erg_s-d.P_dir_eq_erg_s-RC::MeVToErg*(ee*d.sigma_count_s[0]+em*d.sigma_count_s[1]);
    for(std::size_t i=0;i<source.events.size();++i)
    {
        d.Eesc_fluid_inf_MeV_s+=energies[i].Eesc_fluid_inf_MeV*source.events[i].rate_count_s;
        d.Eesc_star_inf_MeV_s+=energies[i].Eesc_star_inf_MeV*source.events[i].rate_count_s;
        d.EX_inf_MeV_s+=energies[i].EX_inf_MeV*source.events[i].rate_count_s;
    }
    d.L_out_fluid_inf_erg_s=RC::MeVToErg*d.Eesc_fluid_inf_MeV_s;d.L_esc_star_inf_erg_s=RC::MeVToErg*d.Eesc_star_inf_MeV_s;d.J_X_inf_erg_s=RC::MeVToErg*d.EX_inf_MeV_s;
    d.finite_T_omitted_floor_erg_s=partition_->OmittedFiniteTemperaturePowerErgPerSecond(source,T);
    d.frozen_drift_bound=d.frozen.drift_bound;d.frozen_threshold=d.frozen.threshold;
    d.DeltaN_over_N=d.frozen.DeltaN_over_N;
    d.LH_erg_s=out.ordinary.beta.heating_erg_s;d.DeltaLnu_erg_s=out.ordinary.beta.neutrino_increment_erg_s;
    d.DeltaPbeta_erg_s=out.ordinary.beta.incremental_beta_erg_s;d.Lnu_eq_erg_s=out.ordinary.reaction.EquilibriumErgPerSecond();
    d.Lnu_full_erg_s=d.Lnu_eq_erg_s+d.DeltaLnu_erg_s;d.Lgamma_erg_s=out.ordinary.Lgamma_erg_s;d.Lother_erg_s=out.ordinary.Lother_neutrino_erg_s;
    d.Pnet_erg_s=d.P_dir_actual_erg_s+d.DeltaPbeta_erg_s-d.Lnu_eq_erg_s-d.Lgamma_erg_s-d.Lother_erg_s;
    d.Cstar_erg_K=out.ordinary.Cstar_erg_K;d.Tinf_K=T;d.Tsurface_inf_K=out.ordinary.Tsurface_infinity_K;d.valid_through_sample=true;
    for(double v:{out.eta_dot_MeV_s[0],out.eta_dot_MeV_s[1],out.x_dot_s,d.Echem_MeV,d.Pnet_erg_s,d.R18_residual_erg_s})
        Need(std::isfinite(v),"nonfinite controlled BNV evaluation");
    RequireCheapCurrent();
    return out;
}

} // namespace CompactStar::Physics::BNV
