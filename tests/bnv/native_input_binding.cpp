// Candidate engineering identity test; no BNV RHS or evolution is invoked.
#include "campaign_fixture.hpp"
#include <iomanip>

int main(int argc,char** argv)
{
 try {
    require(argc==6,"profile thermal entry coefficients fresh-work");
    const std::filesystem::path original=argv[1],thermal_path=argv[2],entry=argv[3],coefficients=argv[4],work=argv[5];
    require(!std::filesystem::exists(work),"work directory must be fresh");
    std::filesystem::create_directories(work);
    const auto profile=work/"profile";
    std::filesystem::copy(original,profile,std::filesystem::copy_options::recursive);
    // Only this private fresh copy is mutated. The scientific archive is read-only.
    require(std::filesystem::canonical(profile)!=std::filesystem::canonical(original),"private fixture aliases archive");
    require(std::filesystem::canonical(profile/"model.txt").parent_path()==std::filesystem::canonical(profile),"private model escapes fixture");
    require(!std::filesystem::equivalent(profile/"model.txt",original/"model.txt"),"private model aliases archive inode");
    const auto certificate=profile/"certificate.txt";
    auto fixture=Fixture(profile,certificate.string(),work/"owning-star",80000);
    auto channels=Channels(fixture);
    EOS::CompOSE_Thermo::Options options;options.Tmin_for_derivative_MeV=0;options.clamp_to_domain=false;
    auto thermal=std::make_shared<const RC::FrozenThermalSource>(thermal_path,options,
      "controlled mathematical fixed-background free-gas entropy; qualified radial80000");
    auto token=std::make_shared<RC::RunDependencyToken>();
    auto spin=std::make_shared<const Phase6A1Campaign::BNV::StaticZeroSpinHistory>(token);
    auto q=Phase6A1Campaign::Qualification(profile,certificate,entry);
    MustRefuse([&]{Context(fixture,channels,thermal,spin,q,token);},"default historical identity refuses native bytes");
    q.background_inputs=RC::QualifiedBackgroundInputs::NativeStructure1;
    auto context=Context(fixture,channels,thermal,spin,q,token);context->RequireFullCurrent();
    std::cout<<"NATIVE_CONTEXT_EXACT_OWNING_STAR PASS\n";
    for(int which=0;which<4;++which) {
      auto bad=q;std::string* paths[]={&bad.profile_path,&bad.model_path,&bad.certificate_path,&bad.eos_path};
      *paths[which]=(work/"missing").string();
      MustRefuse([&]{Context(fixture,channels,thermal,spin,bad,token);},"missing identity input");
      *paths[which]=entry.string();
      MustRefuse([&]{Context(fixture,channels,thermal,spin,bad,token);},"mixed or foreign identity input");
    }
    std::filesystem::copy_file(profile/"profile.tsv",work/"identical-profile.tsv");
    auto foreign=q;foreign.profile_path=(work/"identical-profile.tsv").string();
    MustRefuse([&]{Context(fixture,channels,thermal,spin,foreign,token);},"equal bytes outside semantic owner");
    auto invalid=q;invalid.background_inputs=static_cast<RC::QualifiedBackgroundInputs>(123);
    MustRefuse([&]{Context(fixture,channels,thermal,spin,invalid,token);},"unknown identity");
    auto production=q;production.purpose=RC::RunPurpose::ControlledTrajectory;
    MustRefuse([&]{Context(fixture,channels,thermal,spin,production,token);},"native identity does not authorize production trajectory");
    auto src=source(profile/"freegas.tsv");auto numbers=ParticleNumbers::Compute(input(*fixture.central,src));
    double B0=0;for(std::size_t i=0;i<species.size();++i)if(species[i].label=="10"||species[i].label=="11")B0+=numbers.Values()[i];
    require(B0==Phase6A1Campaign::Column(coefficients,"B_solved"),"native coefficient B0 differs from owning star");
    const double rate=-1e-12*B0/Year;
    require(rate==-2.413652026364133e37,"native derived Bdot identity changed");
    require(rate!=-2.4136520263641375e37,"historical/native distinction lost");
    std::cout<<std::setprecision(17)<<"NATIVE_B0 "<<B0<<" NATIVE_BDOT "<<rate<<"\n";
    ++token->generation;MustRefuse([&]{context->RequireCheapCurrent();},"revoked source token");
    --token->generation; // test-only restoration, followed by independent byte mutation
    {std::ofstream f(profile/"model.txt",std::ios::app);f<<"\nmutation test\n";}
    MustRefuse([&]{context->RequireFullCurrent();},"modified retained source bytes");
    MustRefuse([&]{Context(fixture,channels,thermal,spin,q,token);},"modified native identity input");
    std::cout<<"PASS native input binding; zero RHS calls, zero trajectories\n";
    return 0;
 } catch(const std::exception& e){std::cerr<<"STOP "<<e.what()<<'\n';return 1;}
}
