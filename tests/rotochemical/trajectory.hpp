#include "radial_diagnostic.hpp"
struct TrajectoryCapture final:EV::Observers::IObserver {
 std::ofstream out,jacobian;const RC::SecularEvolutionDriver& driver;
 TrajectoryCapture(const std::filesystem::path& path,const RC::SecularEvolutionDriver& d):out(path),jacobian(path.string()+".jacobian"),driver(d){out<<std::setprecision(17)<<"t_s Tinf_K Tsurface_inf_K eta_e_MeV eta_mu_MeV xi_e xi_mu Omega Omega_dot R_e R_mu Lnu_eq DeltaLnu LH DeltaPbeta Lgamma Lother_neutrino Cstar x_dot eta_dot_e eta_dot_mu Pnet Lnu_full\n";jacobian<<std::setprecision(17)<<"t_s J00 J01 J02 J10 J11 J12 J20 J21 J22\n";}
 void Save(double t,const EV::StateVector&s,const EV::DriverContext&ctx){
  const double T=s.GetThermal().Tinf();const auto eta=RC::ChemicalImbalanceState::Read(s.GetChem());auto v=driver.Evaluate(t,s,ctx);
  const double xe=eta.InfinityMeV(RC::BetaChannel::Npe),xm=eta.InfinityMeV(RC::BetaChannel::NpMu);
  std::vector<double> values{t,T,v.Tsurface_infinity_K,xe,xm,eta.Xi(RC::BetaChannel::Npe,T),eta.Xi(RC::BetaChannel::NpMu,T),v.spin.omega_rad_s,v.spin.omega_dot_rad_s2,v.reaction.rate_count_s[0],v.reaction.rate_count_s[1],v.reaction.EquilibriumErgPerSecond(),v.beta.neutrino_increment_erg_s,v.beta.heating_erg_s,v.beta.incremental_beta_erg_s,v.Lgamma_erg_s,v.Lother_neutrino_erg_s,v.Cstar_erg_K,v.x_dot_s,v.eta_dot_MeV_s[0],v.eta_dot_MeV_s[1],v.Pnet_erg_s,v.reaction.EquilibriumErgPerSecond()+v.beta.neutrino_increment_erg_s};
  for(double x:values)require(std::isfinite(x),"invalid trajectory checkpoint");for(size_t i=0;i<values.size();++i)out<<(i?" ":"")<<values[i];out<<std::endl;
  for(double age:{1e6,1e8,1e10})if(std::abs(t/Year/age-1)<1e-10){double j[3][3];for(size_t col=0;col<3;++col){P::State::ThermalState th;th.Resize(1);P::State::ChemState ch;ch.Resize(2);EV::StateVector local;local.Register(Tag::Thermal,th);local.Register(Tag::Chem,ch);double eps=col==0?1e-5:std::max(1e-10,std::abs(col==1?xe:xm)*1e-5);std::array<RC::SecularEvaluation,2> eval;
    for(size_t sign=0;sign<2;++sign){double step=sign==0?-eps:eps;th.SetTinf(T*std::exp(col==0?step:0));RC::ChemicalImbalanceState(xe+(col==1?step:0),xm+(col==2?step:0)).Store(ch);eval[sign]=driver.Evaluate(t,local,ctx);}j[0][col]=(eval[1].x_dot_s-eval[0].x_dot_s)/(2*eps);for(size_t row=1;row<3;++row)j[row][col]=(eval[1].eta_dot_MeV_s[row-1]-eval[0].eta_dot_MeV_s[row-1])/(2*eps);}
   jacobian<<t;for(const auto&row:j)for(double vj:row)jacobian<<' '<<vj;jacobian<<std::endl;
  }
 }
 void OnStart(const EV::Observers::RunInfo&,const EV::StateVector&s,const EV::DriverContext&c)override{Save(0,s,c);}
 void OnSample(const EV::Observers::SampleInfo&i,const EV::StateVector&s,const EV::DriverContext&c)override{Save(i.t,s,c);}
};
void Trajectory(std::shared_ptr<const RC::FrozenRotochemicalRunContext> context,const std::filesystem::path& output,bool refined=false,double T0=1e8,double xi0=0){
 require(context->Qualification().purpose==RC::RunPurpose::ControlledTrajectory,"analytic context cannot authorize physical trajectory");RunState r(context);r.thermal.SetTinf(T0);RC::ChemicalImbalanceState(xi0*RC::BoltzmannMeVPerK*T0,xi0*RC::BoltzmannMeVPerK*T0).Store(r.chem);
 auto driver=std::make_shared<RC::SecularEvolutionDriver>(context);EV::EvolutionSystem system(r.ctx,r.state,r.rhs,r.layout,{driver});system.AddObserver(std::make_shared<TrajectoryCapture>(output,*driver));RC::ScaledRKF45 solver(system,r.layout,context,refined?RC::ComponentTolerances::Refined():RC::ComponentTolerances{});RC::IntegrationStatistics stats;auto y=r.Pack();std::vector<double> checkpoints;for(unsigned i=0;i<=400;++i)checkpoints.push_back(Year*std::pow(10.,i/40.));checkpoints.back()=1e10*Year;
 auto save=[&]{std::ofstream out(output.string()+".steps");out<<std::setprecision(17)<<"accepted rejected RHS min_step_s max_step_s\n"<<stats.accepted_steps<<' '<<stats.rejected_steps<<' '<<stats.rhs_evaluations<<' '<<stats.minimum_step_s<<' '<<stats.maximum_step_s<<"\nt_s cumulative_accepted cumulative_rejected last_step_s\n";for(auto v:stats.outputs)out<<v.time_s<<' '<<v.accepted<<' '<<v.rejected<<' '<<v.last_step_s<<'\n';};
 try{solver.Integrate(0,checkpoints,y.data(),stats,true);}catch(...){save();throw;}save();std::cout<<"TRAJECTORY COMPLETE "<<output<<" accepted "<<stats.accepted_steps<<" rejected "<<stats.rejected_steps<<" RHS "<<stats.rhs_evaluations<<std::endl;
}
void FixtureReport(const RC::FrozenRotochemicalRunContext& context){
 for(auto c:{RC::BetaChannel::Npe,RC::BetaChannel::NpMu})std::cout<<"FROZEN_ROW "<<RC::ChannelIndex(c)<<" Z "<<context.Z(c,RC::BetaChannel::Npe)<<' '<<context.Z(c,RC::BetaChannel::NpMu)<<" W "<<context.W(c)<<" I "<<context.I(c)<<'\n';
 std::cout<<"RESULT W 2 "<<context.W(RC::BetaChannel::Npe)<<' '<<context.W(RC::BetaChannel::NpMu)<<"\nRESULT I 2 "<<context.I(RC::BetaChannel::Npe)<<' '<<context.I(RC::BetaChannel::NpMu)<<'\n';
 for(auto p:{RC::UrcaProcess::Me,RC::UrcaProcess::Mmu}){const auto&e=context.ChannelEntry(p);std::cout<<"FROZEN_LTILDE "<<RC::ProcessIndex(p)<<' '<<e.luminosity_erg_s_Kq<<" erg s^-1 K^-8 "<<e.normalization_identity;for(auto s:e.support)std::cout<<" support_km "<<s.left_km<<' '<<s.right_km;std::cout<<'\n';}
 for(const auto&s:context.ThermalSource().Sources())std::cout<<"THERMAL_SOURCE "<<s.sha256<<' '<<s.path<<'\n';
 std::cout<<"SPIN_AUTHORITY "<<context.SpinOwner()->Identity()<<"\nTHERMAL_AUTHORITY StarContext::HeatCapacityStar_Tinf; 160 log-T cache; PhotonCooling FR2005 equation (49) / PCY97 fully accreted envelope rho_b=1e10; non-controlled NeutrinoCooling DU/MU/PBF=false\n";
}
