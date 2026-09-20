#define PHASE6A1_BA12R_ULTRA_SUPPORT_ONLY
#include "coupled_trajectory.cpp"

int main(int argc,char** argv)
{
 try
 {
    gsl_set_error_handler_off();
    require(argc==11,"profile certificate thermal work entry frozen-certificate coefficients output pretrajectory-record source-or-control");
    std::cout<<std::setprecision(17)<<std::unitbuf;
    const std::filesystem::path profile=argv[1],certificate=argv[2],thermal_path=argv[3],work=argv[4],entry=argv[5],frozen=argv[6],coefficients=argv[7],output=argv[8],pretrajectory=argv[9];
    const std::string requested=argv[10];
    require(requested=="source"||requested=="control","BA12R mode must be source or control");
    require(!std::filesystem::exists(work),"fresh BA12R assembly directory required");
    require(!std::filesystem::exists(output)&&!std::filesystem::exists(output.string()+".steps"),"fresh BA12R output path required");
    {std::ifstream in(pretrajectory);std::string text((std::istreambuf_iterator<char>(in)),{});require(text.find("PRETRAJECTORY PASS")!=std::string::npos&&text.find("no BNV trajectory had been generated")!=std::string::npos,"committed pretrajectory authorization missing");}

    const auto& card=Campaign::RunCards().back();
    require(card.identity=="CPL-P2-LINEAR-QSS-v1"&&card.partition=="P2"&&
      card.fractional_drive_per_year==-1.0e-12&&card.duration_year==5.0e5&&
      card.checkpoints==8193&&!card.reaction_free,"immutable BA12R P2 card changed");
    const RC::ComponentTolerances ultra{1.0e-11,{1.0e-16,1.0e-22,1.0e-22}};
    ultra.Validate();

    std::filesystem::create_directories(output.parent_path());
    auto fixture=Fixture(profile,certificate.string(),work/"owning-star",80000);
    auto channels=Channels(fixture);
    EOS::CompOSE_Thermo::Options options;options.Tmin_for_derivative_MeV=0;options.clamp_to_domain=false;
    auto thermal=std::make_shared<const RC::FrozenThermalSource>(thermal_path,options,"controlled mathematical fixed-background free-gas entropy; qualified radial80000");
    auto tangent=Campaign::Tangent(fixture,profile,work);
    auto bnv_token=std::make_shared<BNV::BnvDependencyToken>();
    auto monitor=Phase6A1Test::LoadFrozenMonitor(frozen,tangent,bnv_token);
    auto run_token=std::make_shared<RC::RunDependencyToken>();
    auto spin=std::make_shared<const BNV::StaticZeroSpinHistory>(run_token);
    auto ordinary=Context(fixture,channels,thermal,spin,Campaign::Qualification(profile,certificate,entry),run_token);
    const double mu_B=Campaign::Column(coefficients,"mu_B_inf");
    const std::string potential="frozen Structure-1 B0 equilibrium mu_B plus governed moving-reference actual-potential correction";

    std::shared_ptr<const BNV::FrozenControlledBnvRunContext> run_context;
    if(requested=="source")
    {
        const double Bdot=card.fractional_drive_per_year*tangent->B0Count()/Year;
        require(Bdot==-2.4136520263641375e37,"immutable BA12R Bdot changed");
        const std::string fate_id="phase6a1-P2-terminal-fate-v1";
        BNV::ProductFateLedger fate(fate_id,{{"controlled-terminal-product",BNV::TerminalProductFate::BoundInert,1,"abstract-neutron-disappearance"}});
        auto partition=Partition(card,fixture,bnv_token,fate_id);
        auto history=std::make_shared<const Phase6A1Test::UniformProperNeutronSinkHistory>(fixture.central,tangent->B0Count(),Bdot,tangent->DomainIdentity(),tangent->StarIdentity(),tangent->SequenceStateIdentity(),partition->Identity(),fate_id,bnv_token);
        run_context=std::make_shared<const BNV::FrozenControlledBnvRunContext>(ordinary,tangent,history,partition,fate,monitor,channels,mu_B,potential,card.identity);
    }
    else
    {
        auto zero_history=std::make_shared<const Campaign::ExactZeroHistory>(tangent->B0Count(),tangent->DomainIdentity(),tangent->StarIdentity(),tangent->SequenceStateIdentity(),bnv_token);
        auto zero_partition=std::make_shared<const Campaign::ExactZeroPartition>(bnv_token);
        BNV::ProductFateLedger zero_fate("zero-control-fate",{{"none",BNV::TerminalProductFate::BoundInert,1,"zero-control"}});
        run_context=std::make_shared<const BNV::FrozenControlledBnvRunContext>(ordinary,tangent,zero_history,zero_partition,zero_fate,monitor,channels,mu_B,potential,card.identity+"-MATCHED-CONTROL");
    }

    const auto summary=Run(run_context,BNV::ControlledBnvSecularDriver::Mode::Coupled,card,output,ultra);
    std::cout<<"BA12R_ULTRA PASS "<<requested<<" rows "<<summary.rows
             <<" accepted "<<summary.accepted<<" rejected "<<summary.rejected
             <<" rhs "<<summary.rhs<<" minimum_step_s "<<summary.min_step
             <<" maximum_step_s "<<summary.max_step<<'\n';
    return 0;
 }
 catch(const std::exception& e){std::cerr<<"STOP "<<e.what()<<'\n';return 1;}
}
