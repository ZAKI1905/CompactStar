// Independent Phase-5D tests; existing synthetic G construction supplies semantic ownership.
#define main phase5c_contract_main
#include "../analysis/chemical_production_contract.cpp"
#undef main
#include <CompactStar/Physics/Rotochemical/RotochemicalReactionResponse.hpp>
using namespace CompactStar::Physics::Rotochemical;
void near(double x,double y,double relative,const char* label)
{ require(std::isfinite(x)&&std::abs(x-y)<=relative*std::max(std::abs(y),1e-290),label); }
UrcaIntegrationRequest urca_request()
{
    auto g=std::make_shared<GlobalChemicalNumberResponse>(GlobalChemicalNumberResponse::Compute(request()));
    UrcaIntegrationRequest r;r.chemical_domain=g;r.domain_identity=g->Lifetime()->revision->domain;
    r.metric_identity="independent analytic metric";r.metric=[](double){return UrcaMetric{0,0};};r.radial_partition_km={0,.5,1};
    r.normalizations.emplace_back(UrcaProcess::Me,[](double){return 1e-40;},std::vector<UrcaSupportInterval>{{0,1}},"positive Me benchmark");
    r.normalizations.emplace_back(UrcaProcess::Mmu,[](double){return 3e-40;},std::vector<UrcaSupportInterval>{{0,.75}},"positive Mmu benchmark");
    return r;
}
int main(int argc,char** argv)
{
    try {
        std::cout<<std::setprecision(17);
        if(argc>1 && std::string(argv[1])=="functions") {
            for(double x:{-1e4,-10.,-5.,-1.,-1e-5,0.,1e-5,.1,.5,1.,2.,4.,5.,8.,10.,1e4})
                std::cout<<x<<' '<<UrcaImbalanceFunctions::FD(x)<<' '<<UrcaImbalanceFunctions::HD(x)<<' '<<UrcaImbalanceFunctions::FM(x)<<' '<<UrcaImbalanceFunctions::HM(x)<<'\n';
            return 0;
        }
        if(argc>1 && std::string(argv[1])=="integrals") {
            for(auto p:{UrcaProcess::De,UrcaProcess::Me}) for(int mutant=0;mutant<9;++mutant) {
                auto r=urca_request();r.selection=UrcaProcessSelection{p};r.normalizations.clear();
                const int q=TemperatureExponent(p);
                // At q=6, changing q to 8 equals the extra-two-lapse mutant (3).
                // It is one algebraic case, not an additional mutation credit.
                if(q==6 && mutant==7) continue;
                // Equivalent transformed-input mutants exercise the ACTUAL integrator.
                // Lapse: remove all, missing one, G_y substitution, extra, sign-flipped.
                r.metric=[mutant,q](double x){ double nu=-.4+.1*x*x,lambda=.2*x*x;
                    if(mutant==1)nu=0;
                    if(mutant==2)nu*=double(1-q)/(2-q);
                    if(mutant==3)nu*=double(-q)/(2-q);
                    if(mutant==4)nu*=double(3-q)/(2-q);
                    if(mutant==5)nu=-nu;
                    if(mutant==6)lambda=0;
                    if(mutant==7)nu*=double(2-(q==6?8:6))/(2-q);
                    return UrcaMetric{nu,lambda};};
                auto support=mutant==8?std::vector<UrcaSupportInterval>{{0,1}}:std::vector<UrcaSupportInterval>{{.2,.45},{.7,.9}};
                r.normalizations.emplace_back(p,[](double x){return 1e-40*(2+x*x);},support,"RE10b positive even outside support");
                auto c=GlobalUrcaChannelCoefficient::Compute(r);
                std::cout<<"INTEGRAL "<<q<<' '<<mutant<<' '<<c.LuminosityCoefficient(p)<<'\n';
                if(mutant==0) {
                    r.normalizations.clear();r.normalizations.emplace_back(p,[](double x){return 1e-40*(2+x*x)*std::exp(-.4*x*x);},support,"inverted proper volume transformed S");
                    std::cout<<"INTEGRAL "<<q<<" 9 "<<GlobalUrcaChannelCoefficient::Compute(r).LuminosityCoefficient(p)<<'\n';
                }
            }
            return 0;
        }
        const double pi=std::acos(-1.0);
        const double kb=double(1.380649e-23L/1.602176634e-13L);
        near(BoltzmannMeVPerK,kb,2e-15,"RE1 kB MeV");near(BoltzmannErgPerK,1.380649e-16,2e-15,"RE1 kB erg");
        near(MeVToErg,1.602176634e-6,2e-15,"RE1 energy conversion");
        ChemicalImbalanceState state(.012,-.025);CompactStar::Physics::State::ChemState storage;storage.Resize(2);state.Store(storage);
        require(storage.Eta(0)==.012&&storage.Eta(1)==-.025,"RE1 state ordering");
        for(auto ch:{BetaChannel::Npe,BetaChannel::NpMu}) for(double nu:{-.4,-.2,0.})
            near(state.LocalMeV(ch,nu)/(kb*1e7*std::exp(-nu)),state.Xi(ch,1e7),2e-15,"RE2 xi redshift");
        for(double x:{-1e4,-10.,-1.,-.01,0.,.01,1.,10.,1e4}) {
            near(UrcaImbalanceFunctions::FD(x),UrcaImbalanceFunctions::FD(-x),0,"RE4 FD parity");
            near(UrcaImbalanceFunctions::FM(x),UrcaImbalanceFunctions::FM(-x),0,"RE4 FM parity");
            require(UrcaImbalanceFunctions::HD(x)==-UrcaImbalanceFunctions::HD(-x)&&UrcaImbalanceFunctions::HM(x)==-UrcaImbalanceFunctions::HM(-x),"RE4 H odd");
            require(x*UrcaImbalanceFunctions::HM(x)>=0&&x*UrcaImbalanceFunctions::HD(x)>=0,"RE7 sign");
        }
        near(UrcaImbalanceFunctions::HM(1e-5)/1e-5,14680/(11513*pi*pi),1e-10,"RE5 HM slope");
        near(UrcaImbalanceFunctions::HD(1e-5)/1e-5,714/(457*pi*pi),1e-10,"RE5 HD slope");
        near(UrcaImbalanceFunctions::ModifiedIncrement(1e-5)/1e-10,22020/(11513*pi*pi),1e-10,"RE5 FM slope");
        near(UrcaImbalanceFunctions::DirectIncrement(1e-5)/1e-10,1071/(457*pi*pi),1e-10,"RE5 FD slope");
        near(UrcaImbalanceFunctions::HM(1e4)/std::pow(1e4,7),24/(11513*std::pow(pi,8)),1e-5,"RE6 HM pi8");
        near(UrcaImbalanceFunctions::FM(1e4)/std::pow(1e4,8),9/(11513*std::pow(pi,8)),1e-5,"RE6 FM asymptote");
        near(UrcaImbalanceFunctions::HD(1e4)/std::pow(1e4,5),42/(457*std::pow(pi,6)),1e-5,"RE6 HD asymptote");
        for(auto v:{std::array<double,3>{4.7870134733368985,0,1},{4.909710028924132,1,1},{5.458531594867599,0,0},{5.633717467648343,1,0}}) {
            double x=v[0], f=v[1]?UrcaImbalanceFunctions::FM(x):UrcaImbalanceFunctions::FD(x),h=v[1]?UrcaImbalanceFunctions::HM(x):UrcaImbalanceFunctions::HD(x);
            require(std::abs(x*h-f+v[2])<1e-12,"RE6 defined full/incremental root");
        }
        auto r=urca_request();auto coefficients=std::make_shared<GlobalUrcaChannelCoefficient>(GlobalUrcaChannelCoefficient::Compute(r));
        near(coefficients->LuminosityCoefficient(UrcaProcess::Me),4*pi/3*1e-25,2e-15,"RE10b flat prefactor");
        RotochemicalReactionResponse response(coefficients);
        auto zero=response.Evaluate(1e7,{0,0});auto power0=RotochemicalThermalPower::From(zero);
        require(zero.Rate(BetaChannel::Npe)==0&&zero.Rate(BetaChannel::NpMu)==0&&zero.IncrementErgPerSecond()==0&&power0.heating_erg_s==0,"RE8 exact zero");
        require(response.Coefficients().get()==coefficients.get(),"RE9 one authority");
        for(double t:{1e7,2e7}) {
            auto result=response.Evaluate(t,state);auto power=RotochemicalThermalPower::From(result);
            double expectedH=0;
            for(auto p:{UrcaProcess::Me,UrcaProcess::Mmu}) {
                auto ch=LeptonChannel(p);double x=state.InfinityMeV(ch)/(kb*t),u=x/pi;
                double H=(14680*u+7560*std::pow(u,3)+840*std::pow(u,5)+24*std::pow(u,7))/(11513*pi);
                double rate=coefficients->LuminosityCoefficient(p)/1.380649e-16*std::pow(t,7)*H;
                near(result.Rate(ch),rate,8e-15,"RE12 rate exponent and kB units");
                require(state.InfinityMeV(ch)*result.Rate(ch)>0,"RE7 eta R positive");
                expectedH+=coefficients->LuminosityCoefficient(p)*std::pow(t,8)*x*H;
            }
            near(power.heating_erg_s,expectedH,8e-15,"thermal conversion once");
            near(power.incremental_beta_erg_s-result.EquilibriumErgPerSecond(),power.full_beta_erg_s,2e-15,"RE9 thermal ledger");
        }
        // Explicit disconnected support: inner closed region must never be sampled.
        auto disconnected=r;disconnected.normalizations.clear();
        disconnected.selection=UrcaProcessSelection{UrcaProcess::Me};
        disconnected.normalizations.emplace_back(UrcaProcess::Me,[](double x){require((x>.2&&x<.45)||(x>.7&&x<.9),"RE17 closed region swept");return 1e-40;},std::vector<UrcaSupportInterval>{{.2,.45},{.7,.9}},"disconnected support");
        auto dc=GlobalUrcaChannelCoefficient::Compute(disconnected);
        near(dc.LuminosityCoefficient(UrcaProcess::Me),4*pi/3*1e-25*(std::pow(.45,3)-std::pow(.2,3)+std::pow(.9,3)-std::pow(.7,3)),3e-15,"RE17 disconnected volume");
        require(DirectUrcaTriangle(1e-8,1e-8,1e-8),"RE17 constructible below-guard triangle");
        r.normalizations.emplace_back(UrcaProcess::De,[](double){throw std::runtime_error("disabled DU evaluated");return 1.;},std::vector<UrcaSupportInterval>{{.8,1}},"triangle-open disabled DU");
        require(GlobalUrcaChannelCoefficient::Compute(r).LuminosityCoefficient(UrcaProcess::De)==0,"RE17 process selection");
        r.domain_identity+="wrong";refuse([&]{GlobalUrcaChannelCoefficient::Compute(r);},"domain mismatch");r=urca_request();
        std::reverse(r.radial_partition_km.begin(),r.radial_partition_km.end());refuse([&]{GlobalUrcaChannelCoefficient::Compute(r);},"innermost-first");
        coefficients->ChemicalDomain()->Lifetime()->revision->domain+="changed";
        refuse([&]{response.Evaluate(1e7,state);},"Stale");
        std::cout<<"PASS focused state/units/functions/rates/response-ledger/support/currentness checks; coupled RE ladder pending\n";return 0;
    } catch(const std::exception& e) {std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}
}
