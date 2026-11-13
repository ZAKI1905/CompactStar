// ADR-0025 S2: short simultaneous-driving candidate. No Stage B rerun.
#define NATIVE_STAGE_NO_MAIN
#include "stage.cpp"
#undef NATIVE_STAGE_NO_MAIN
#include <CompactStar/Physics/BNV/ActualPotentialPartition.hpp>

int main(int argc,char** argv){try{
    gsl_set_error_handler_off();std::cout<<std::setprecision(17)<<std::unitbuf;
    require(argc==10,"profile thermal entry frozen coefficients archive sha fresh-work output");
    Inputs in;in.profile=argv[1];in.certificate=in.profile/"certificate.txt";in.thermal=argv[2];
    in.entry=argv[3];in.frozen=argv[4];in.coefficients=argv[5];in.work=argv[8];
    const std::filesystem::path out=argv[9];std::filesystem::create_directories(out);
    auto sequence=Core::StellarArchive::Load(argv[6],argv[7]);
    auto fixture=NativeFixture(in.profile,in.certificate.string(),in.work/"owning-star",80000,sequence);
    auto channels=Channels(fixture);
    EOS::CompOSE_Thermo::Options options;options.Tmin_for_derivative_MeV=0;options.clamp_to_domain=false;
    auto thermal=std::make_shared<const RC::FrozenThermalSource>(in.thermal,options,"controlled mathematical fixed-background free-gas entropy; qualified radial80000");
    auto src=source(in.profile/"freegas.tsv");auto numbers=ParticleNumbers::Compute(input(*fixture.central,src));numbers.RequireCurrent();
    const double B0=numbers.Values()[0]+numbers.Values()[1];
    require(B0==Campaign::Column(in.coefficients,"B_solved"),"native B0 changed");
    auto tangent=std::make_shared<const BNV::EquilibriumBaryonTangent>(BNV::EquilibriumBaryonTangent::Compute(
        sequence,B0,"Structure-1 rho_c=1.10e15 radial80000 EOS8192",BNV::EquilibriumBaryonTangent::SerializeDomainIdentity(whole)));
    for(int j=0;j<2;++j)require(tangent->ClosedValues()[j+1]==Campaign::Column(in.coefficients.parent_path()/"tangent.tsv",j?"t_mu":"t_e"),"native tangent changed");
    for(int i=0;i<2;++i)for(int j=0;j<2;++j)require(fixture.z->Values()[i][j]==Campaign::Column(in.coefficients,"Z"+std::to_string(i)+std::to_string(j)),"native Z changed");
    auto token=std::make_shared<BNV::BnvDependencyToken>();auto run_token=std::make_shared<RC::RunDependencyToken>();
    auto monitor=Phase6A1Test::LoadFrozenMonitor(in.frozen,tangent,token);
    const double mu=Campaign::Column(in.coefficients,"mu_B_inf"),gamma=1e-12,period=.02,beta_target=1.05;
    const double D=gamma*B0/Year;
    auto q=Campaign::Qualification(in.profile,in.certificate,in.entry);q.background_inputs=RC::QualifiedBackgroundInputs::NativeStructure1;
    auto offspin=std::make_shared<const BNV::StaticZeroSpinHistory>(run_token);
    auto off=Context(fixture,channels,thermal,offspin,q,run_token);
    const double drive=D*tangent->ClosedValues()[1]/(2*std::abs(off->I(RC::BetaChannel::Npe))*beta_target);
    const double pdot=drive*std::pow(period,3)/(4*M_PI*M_PI);
    auto spin=std::make_shared<const RC::ParameterizedDipoleHistory>(period,pdot,run_token);
    auto combined_q=q;combined_q.purpose=RC::RunPurpose::PrescribedCombinedDriving;
    auto on=Context(fixture,channels,thermal,spin,combined_q,run_token);
    MustRefuse([&]{Context(fixture,channels,thermal,std::make_shared<const RC::ParameterizedDipoleHistory>(.0001,pdot,run_token),combined_q,run_token);},"fast-spin admission");
    MustRefuse([&]{Context(fixture,channels,thermal,std::make_shared<const TestSpin>(10,-1e-10),combined_q,run_token);},"unqualified spin adapter");
    const std::string fate_id="ADR-0025-P-A-prompt-escape";
    BNV::ProductFateLedger fate(fate_id,{{"all-products-escape",BNV::TerminalProductFate::PromptEscape,1,"abstract-neutron-disappearance"}});
    auto partition=std::make_shared<const BNV::ActualPotentialPartition>(0,fate_id,token);
    auto history=std::make_shared<const Phase6A1Test::UniformProperNeutronSinkHistory>(fixture.central,B0,-D,
        tangent->DomainIdentity(),tangent->StarIdentity(),tangent->SequenceStateIdentity(),partition->Identity(),fate_id,token);
    auto zero=std::make_shared<const Campaign::ExactZeroHistory>(B0,tangent->DomainIdentity(),tangent->StarIdentity(),tangent->SequenceStateIdentity(),token);
    auto make=[&](auto ordinary,std::shared_ptr<const BNV::OrdinaryMatterBnvHistory> h,
                  std::shared_ptr<const BNV::IDirectBnvEnergyPartition> p,const std::string& id){
        return std::make_shared<const BNV::FrozenControlledBnvRunContext>(ordinary,tangent,h,p,fate,monitor,channels,mu,
            "ADR-0015 cold quadratic actual potential; frozen Structure-1; no rotating global closure",id);};
    std::array<Ctx,4> cases{{make(off,zero,partition,"S2-00"),make(on,zero,partition,"S2-RC"),
                             make(off,history,partition,"S2-BNV"),make(on,history,partition,"S2-RC-BNV")}};
    require(cases[1]->SpinOwner()==cases[3]->SpinOwner(),"spin-on history not shared");
    RunState probe(on);RC::ChemicalImbalanceState(.001,-.002).Store(probe.chem);
    auto z=cases[1]->Evaluate(0,probe.state,probe.ctx),c=cases[3]->Evaluate(0,probe.state,probe.ctx);
    auto ordinary=on->Evaluate(0,probe.state,probe.ctx,spin.get());
    require(z.eta_dot_MeV_s==ordinary.eta_dot_MeV_s&&z.x_dot_s==ordinary.x_dot_s,"zero-source ordinary reduction");
    for(int j=0;j<2;++j){const double expected=ordinary.eta_dot_MeV_s[j]+c.diagnostics.eta_dot_from_sigma_MeV_s[j];
        require(c.eta_dot_MeV_s[j]==expected,"combined driving not included exactly once");
        require(std::abs(c.eta_dot_MeV_s[j]-(expected+c.diagnostics.eta_dot_from_spin_MeV_s[j]))>0,"doubled-spin mutation undetected");}
    require(c.direct.power_actual_erg_s==0&&c.diagnostics.L_esc_star_inf_erg_s>0,"zero deposition or escape accounting");
    const double chemscale=std::max({std::abs(c.diagnostics.A_B_erg_s),std::abs(c.diagnostics.A_spin_erg_s),c.diagnostics.LH_erg_s});
    require(std::abs(c.diagnostics.Echem_chain_residual_erg_s)<1e-11*chemscale,"chemical chain identity");
    require(std::abs(c.diagnostics.A_spin_erg_s)>1e-8*chemscale,"missing-spin-work mutation insensitive");
    const double thermal_expected=(c.diagnostics.P_dir_actual_erg_s+c.diagnostics.LH_erg_s-c.diagnostics.Lnu_full_erg_s-c.diagnostics.Lgamma_erg_s)/(c.diagnostics.Cstar_erg_K*c.diagnostics.Tinf_K);
    Near(c.x_dot_s,thermal_expected,1e-13,"thermal ledger");
    require(std::abs((c.x_dot_s-thermal_expected)*c.diagnostics.Cstar_erg_K*c.diagnostics.Tinf_K)<.001*std::abs(c.diagnostics.A_spin_erg_s),"extra-spin-heat mutation not discriminated");
    MustRefuse([&]{cases[3]->EvaluateReactionFree(0,probe.state,probe.ctx);},"combined reaction-free refusal");
    probe.thermal.SetTinf(9e5);MustRefuse([&]{cases[3]->Evaluate(0,probe.state,probe.ctx);},"unchanged lower temperature gate");
    probe.thermal.SetTinf(1.01e9);MustRefuse([&]{cases[3]->Evaluate(0,probe.state,probe.ctx);},"unchanged upper temperature gate");
    probe.thermal.SetTinf(1e8);MustRefuse([&]{cases[3]->Evaluate(2e6*Year,probe.state,probe.ctx);},"unchanged depletion gate");
    // Regression of the new partition to the existing P0 and P2 energy owners.
    for(int retained:{0,1}){
        auto pa=std::make_shared<const BNV::ActualPotentialPartition>(retained,fate_id,token);
        std::shared_ptr<const BNV::IDirectBnvEnergyPartition> old;
        if(retained)old=std::make_shared<const Phase6A1Test::P2Partition>("old-P2",fate_id,100,token);
        else old=std::make_shared<const Phase6A1Test::P0Partition>("old-P0",fate_id,100,token);
        auto hs=[&](const std::string& id){return std::make_shared<const Phase6A1Test::ConstantNeutronSinkHistory>(B0,-D,tangent->DomainIdentity(),tangent->StarIdentity(),tangent->SequenceStateIdentity(),id,fate_id,token);};
        auto cp=make(off,hs(pa->Identity()),pa,"partition-regression-new"),co=make(off,hs(old->Identity()),old,"partition-regression-old");
        RunState r(off);RC::ChemicalImbalanceState(-1e-8,-2e-8).Store(r.chem);
        auto ep=cp->Evaluate(0,r.state,r.ctx),eo=co->Evaluate(0,r.state,r.ctx);
        require(ep.x_dot_s==eo.x_dot_s&&ep.eta_dot_MeV_s==eo.eta_dot_MeV_s,"P-A endpoint partition regression");
    }
    std::ofstream inputs(out/"native-inputs.tsv");inputs<<std::setprecision(17)<<"quantity\tvalue\nP0_s\t"<<period<<"\nPdot0\t"<<pdot<<"\nGamma_yr\t"<<gamma<<"\nbeta_e\t"<<beta_target<<"\nB0\t"<<B0<<"\n";
    for(int j=0;j<2;++j){const auto ch=j?RC::BetaChannel::NpMu:RC::BetaChannel::Npe;
        inputs<<"I"<<j<<'\t'<<on->I(ch)<<"\nW"<<j<<'\t'<<on->W(ch)<<"\nN"<<j<<'\t'<<numbers.Values()[j+2]<<"\nrotation"<<j<<'\t'<<std::abs(on->I(ch))*std::pow(spin->Sample(0).omega_rad_s,2)/numbers.Values()[j+2]<<'\n';}
    inputs.close();
    std::ofstream records(out/"accepted.csv");records<<std::setprecision(17)
      <<"case,tier,t_s,x,eta_e_MeV,eta_mu_MeV,Tinf_K,Tsurface_inf_K,LH,DeltaLnu,Lnu_eq,Lgamma,Pdir,P0,AB,Aspin,Echem_erg,chain_residual,subsystem_residual,R18,escape,beta_e,beta_mu,depletion,frozen_utilization,Cstar,eta_dot_e,eta_dot_mu,xdot\n";
    std::ofstream endpoints(out/"endpoints.csv");endpoints<<std::setprecision(17)<<"case,tier,x,eta_e_MeV,eta_mu_MeV,accepted,rejected\n";
    for(int k=0;k<4;++k){auto context=cases[k];std::array<State,2> final;
        for(int tier=0;tier<2;++tier){const Tol tol=tier?tiers[2]:tiers[1];
            auto main=main_run(context,0,State{{0,0,0}},Year,1,tol,{0,Year});final[tier]=main.final_state;
            ActualCheckpointContext eval(context,1);RunState r(context->OrdinaryContext());
            auto emit=[&](double time,const State& y){EV::UnpackStateVector(r.state,r.layout,y.data());auto e=context->Evaluate(time,r.state,r.ctx);const auto& d=e.diagnostics;
                double scale=std::max({1.,std::abs(d.A_B_erg_s),std::abs(d.A_spin_erg_s),d.LH_erg_s});
                require(std::abs(d.Echem_chain_residual_erg_s)<1e-11*scale,"trajectory chemical identity");
                require(d.P_dir_actual_erg_s==0&&d.frozen.valid&&d.Tinf_K>=1e6&&d.Tinf_K<=1e9,"trajectory gates");
                require(std::abs(y[1])<4e-7&&std::abs(y[2])<8e-7&&y[0]>=-.03&&y[0]<=.2,"short state outside Stage B reviewed box");
                const double be=d.spin_reference_count_s[0]?d.sigma_count_s[0]/-d.spin_reference_count_s[0]:0;
                const double bm=d.spin_reference_count_s[1]?d.sigma_count_s[1]/-d.spin_reference_count_s[1]:0;
                records<<k<<','<<tier<<','<<time<<','<<y[0]<<','<<y[1]<<','<<y[2]<<','<<d.Tinf_K<<','<<d.Tsurface_inf_K<<','<<d.LH_erg_s<<','<<d.DeltaLnu_erg_s<<','<<d.Lnu_eq_erg_s<<','<<d.Lgamma_erg_s<<','<<d.P_dir_actual_erg_s<<','<<d.P_dir_eq_erg_s<<','<<d.A_B_erg_s<<','<<d.A_spin_erg_s<<','<<RC::MeVToErg*d.Echem_MeV<<','<<d.Echem_chain_residual_erg_s<<','<<d.subsystem_residual_erg_s<<','<<d.R18_residual_erg_s<<','<<d.L_esc_star_inf_erg_s<<','<<be<<','<<bm<<','<<-d.DeltaB_over_B0<<','<<d.frozen.max_utilization<<','<<d.Cstar_erg_K<<','<<e.eta_dot_MeV_s[0]<<','<<e.eta_dot_MeV_s[1]<<','<<e.x_dot_s<<'\n';};
            emit(0,main.initial_state);for(const auto& step:main.accepted_steps)emit(step.t_right_s,step.y_right);
            endpoints<<k<<','<<tier<<','<<main.final_state[0]<<','<<main.final_state[1]<<','<<main.final_state[2]<<','<<main.statistics.accepted_steps<<','<<main.statistics.rejected_steps<<'\n';
            if(tier==0){auto scheduled=main_run(context,0,State{{0,0,0}},Year,1,tol,{0,Year/3,Year});
                require(scheduled.final_state==main.final_state&&scheduled.accepted_steps.size()==main.accepted_steps.size(),"passive schedule regression");
                for(size_t i=0;i<main.accepted_steps.size();++i)require(scheduled.accepted_steps[i].y_right==main.accepted_steps[i].y_right&&scheduled.accepted_steps[i].t_right_s==main.accepted_steps[i].t_right_s,"passive accepted-history regression");}
        }
        require(std::abs(final[1][0]-final[0][0])<1e-12&&std::abs(final[1][1]-final[0][1])<1e-19&&std::abs(final[1][2]-final[0][2])<1e-19,"short tolerance comparison");
        if(k==0)require(final[1][1]==0&&final[1][2]==0,"zero-source exact chemistry");
        if(k==1)require(final[1][1]>0&&final[1][2]>0,"spin-only sign");
        if(k==2)require(final[1][1]<0&&final[1][2]<0,"BNV-only sign");
        if(k==3)require(final[1][1]<0&&final[1][2]>0,"combined mixed-sign response");
    }
    records.close();endpoints.close();
    std::ofstream status(out/"status.txt");status<<"PASS short four-case production runs; corner/partition reductions, guards, signs, chain identity, tolerance and passive schedule checks. Independent cross-check pending.\n";
    std::cout<<"S2 SHORT PRODUCTION PASS\n";return 0;
}catch(const std::exception& e){std::cerr<<"STOP "<<e.what()<<'\n';return 1;}}
