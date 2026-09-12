#include "coupled.hpp"
namespace CompactStar::Physics::Rotochemical {
struct FrozenContextTestAccess {static void Semantic(const std::vector<BetaChannel>& c,const Analysis::ChemicalMatrix& m){FrozenRotochemicalRunContext::ValidateSemanticZ(c,m);}};
struct FrozenThermalTestAccess {static std::unique_ptr<const FrozenThermalSource> Race(const std::filesystem::path& p,EOS::CompOSE_Thermo::Options o,const std::function<void()>& hook){return std::unique_ptr<const FrozenThermalSource>(new FrozenThermalSource(p,o,"construction race control",hook));}};
}

// Independent closed-form exponential for A=Z diag(d), never the production RKF45.
std::array<double,2> LinearOracle(const M& z,std::array<double,2>d,std::array<double,2>x,double t){long double a=z[0][0]*d[0],b=z[0][1]*d[1],c=z[1][0]*d[0],e=z[1][1]*d[1],tr=(a+e)/2,k=std::sqrt((a-e)*(a-e)/4+b*c),u=std::exp(-tr*t),ch=std::cosh(k*t),sh=std::sinh(k*t)/k;return {double(u*(ch*x[0]-sh*((a-tr)*x[0]+b*x[1]))),double(u*(ch*x[1]-sh*(c*x[0]+(e-tr)*x[1])))};}
double Lyapunov(const M&z,std::array<double,2>x){long double det=(long double)z[0][0]*z[1][1]-(long double)z[0][1]*z[1][0];return (z[1][1]*x[0]*x[0]-(z[0][1]+z[1][0])*x[0]*x[1]+z[0][0]*x[1]*x[1])/det;}
struct LateInvalidationSpin final:RC::PrescribedSpinHistory {
 const std::shared_ptr<RC::RunDependencyToken> token;const std::string id="late source invalidation test";
 explicit LateInvalidationSpin(std::shared_ptr<RC::RunDependencyToken> t):token(std::move(t)){}
 RC::SpinHistorySample Sample(double)const override{++token->generation;return {0,0};}
 void RequireCurrent()const override{}const std::string& Identity()const override{return id;}
};
struct LyapunovObserver final:EV::Observers::IObserver {
 const M z;double previous;size_t samples=0;
 LyapunovObserver(M matrix,std::array<double,2>x):z(std::move(matrix)),previous(Lyapunov(z,x)){}
 void OnSample(const EV::Observers::SampleInfo&,const EV::StateVector&s,const EV::DriverContext&)override{auto eta=RC::ChemicalImbalanceState::Read(s.GetChem());double v=Lyapunov(z,{eta.InfinityMeV(RC::BetaChannel::Npe),eta.InfinityMeV(RC::BetaChannel::NpMu)});require(v<=previous*(1+1e-13),"Lyapunov increase at accepted checkpoint");previous=v;++samples;}
};
// Test-only transformed evaluation boundary. Mutants enter the actual EvolutionSystem.
// This is not a source-edit mutation build. Exact aliases receive one family credit.
struct FaultDriver final:P::IDriver {
 const RC::SecularEvolutionDriver& nominal;const RC::FrozenRotochemicalRunContext& context;const std::string fault;
 FaultDriver(const RC::SecularEvolutionDriver&n,const RC::FrozenRotochemicalRunContext&c,std::string f):nominal(n),context(c),fault(std::move(f)){}
 std::string Name()const override{return "TestTransformedRHS";}
 const std::vector<Tag>& DependsOn()const override{return nominal.DependsOn();}const std::vector<Tag>& Updates()const override{return nominal.Updates();}
 void AccumulateRHS(double t,const EV::StateVector&s,EV::RHSAccumulator&r,const EV::DriverContext&ctx)const override{
  auto v=nominal.Evaluate(t,s,ctx);auto rate=v.reaction.rate_count_s;double T=s.GetThermal().Tinf(),power=v.Pnet_erg_s;const double eq=v.reaction.EquilibriumErgPerSecond();
  double factor=1;if(fault=="wrong_kB_energy")factor=RC::MeVToErg;if(fault=="omit_kB")factor=RC::BoltzmannErgPerK;if(fault=="double_kB")factor=1/RC::BoltzmannErgPerK;if(fault=="doubled_kB_denominator")factor=.5;if(fault=="T8_rate")factor=T;
  if(fault=="reaction_sign"||fault=="eta_sign")factor=-1;for(auto&x:rate)x*=factor;if(fault=="channel_swap")std::swap(rate[0],rate[1]);
  if(fault=="state_slot_swap") {auto eta=RC::ChemicalImbalanceState::Read(s.GetChem());for(size_t i=0;i<2;++i)rate[i]=context.Ltilde(i==0?RC::UrcaProcess::Me:RC::UrcaProcess::Mmu)/RC::BoltzmannErgPerK*std::pow(T,7)*RC::UrcaImbalanceFunctions::HM(eta.InfinityMeV(i==0?RC::BetaChannel::NpMu:RC::BetaChannel::Npe)/(RC::BoltzmannMeVPerK*T));}
  if(fault=="omit_MeV_erg")power+=v.beta.heating_erg_s*(1/RC::MeVToErg-1);if(fault=="double_MeV_erg")power+=v.beta.heating_erg_s*(RC::MeVToErg-1);if(fault=="omit_heating")power-=v.beta.heating_erg_s;
  if(fault=="double_equilibrium"||fault=="full_not_increment")power-=eq;
  // Explicit legacy-source injection is a separate foreign authority, not another Ltilde.
  if(fault=="legacy_placeholder"){TH::NeutrinoCooling::Options opts;opts.include_direct_urca=false;opts.include_modified_urca=true;opts.include_pair_breaking=false;TH::NeutrinoCooling legacy(opts);auto historical=TH::Detail::NeutrinoCooling_Details::ComputeDerived(legacy,s,ctx);require(historical.ok&&historical.L_nu_MU_inf_erg_s>0,"legacy mutation inactive");power-=historical.L_nu_MU_inf_erg_s;}
  auto spin=v.spin;if(fault=="OmegaDot_sign")spin.omega_dot_rad_s2=-spin.omega_dot_rad_s2;
  for(auto row:{RC::BetaChannel::Npe,RC::BetaChannel::NpMu}){size_t i=RC::ChannelIndex(row);double drive=2*context.W(row)*spin.omega_rad_s*spin.omega_dot_rad_s2;if(fault=="W_sign")drive=-drive;if(fault=="omit_factor2")drive/=2;
   double d=drive;for(auto col:{RC::BetaChannel::Npe,RC::BetaChannel::NpMu})if(!((fault=="drop_cross_Z"||(fault=="drop_cross_Z_row0"&&i==0)||(fault=="drop_cross_Z_row1"&&i==1))&&row!=col))d-=context.Z(row,col)*rate[RC::ChannelIndex(col)];
   // M6 transposes the full nonsymmetric linear map A=Z diag(slopes), not Z.
   if(fault=="transpose_linear_map"){auto eta=RC::ChemicalImbalanceState::Read(s.GetChem());double slope=context.Ltilde(i==0?RC::UrcaProcess::Me:RC::UrcaProcess::Mmu)*(14680/(11513*M_PI*M_PI))*std::pow(T,6)/(RC::BoltzmannErgPerK*RC::BoltzmannMeVPerK);d=0;for(auto col:{RC::BetaChannel::Npe,RC::BetaChannel::NpMu})d-=slope*context.Z(col,row)*eta.InfinityMeV(col);}
   r.AddTo(Tag::Chem,i,d);}
  r.AddTo(Tag::Thermal,0,power/(T*v.Cstar_erg_K));
 }
};
void CoupledOracles(const ControlledFixture& f,std::shared_ptr<const RC::GlobalUrcaChannelCoefficient> channels,std::shared_ptr<const RC::FrozenThermalSource> thermal,RC::RunQualification q){
 auto dipole_token=std::make_shared<RC::RunDependencyToken>();RC::PrescribedDipoleHistory dipole(dipole_token);
 for(double t:{0.,1e6*Year,1e10*Year}){double h=t==0?.01*Year:(t+Year)*1e-5;auto sample=dipole.Sample(t);double derivative=t>h?(dipole.Sample(t+h).omega_rad_s-dipole.Sample(t-h).omega_rad_s)/(2*h):(dipole.Sample(t+h).omega_rad_s-sample.omega_rad_s)/h;Near(derivative,sample.omega_dot_rad_s2,t==0?1e-2:1e-7,"prescribed dipole derivative oracle");}
 std::cout<<"PRESCRIBED_DIPOLE_DERIVATIVE PASS finite-difference independent check\n";
 // Independent five-point finite difference of F; no production polynomial derivative.
 for(double xi:{-100.,-10.,-1.,-.1,0.,.1,1.,10.,100.})for(bool modified:{false,true}){auto F=[&](double x){return modified?RC::UrcaImbalanceFunctions::FM(x):RC::UrcaImbalanceFunctions::FD(x);};double h=1e-3*std::max(1.,std::abs(xi));double derivative=(-F(xi+2*h)+8*F(xi+h)-8*F(xi-h)+F(xi-2*h))/(12*h);double expected=3*(modified?RC::UrcaImbalanceFunctions::HM(xi):RC::UrcaImbalanceFunctions::HD(xi));require(std::abs(derivative-expected)<=1e-8*std::max(1.,std::abs(expected)),"F prime equals 3H oracle failed");}
 std::cout<<"F_PRIME_3H PASS DU and MU signed independent five-point derivative\n";
 q.purpose=RC::RunPurpose::AnalyticControl;auto token=std::make_shared<RC::RunDependencyToken>();auto zero=std::make_shared<const TestSpin>(0,0);const auto z=f.z->Values();
 auto build=[&](auto c,auto spin){return Context(f,c,thermal,spin,q,token);};
 auto no=Channels(f,16,false,RC::UrcaProcessSelection{});double spin_error=0,reaction_error=0;
 for(double sign:{-1.,1.}){auto context=build(no,std::make_shared<const TestSpin>(10,sign));RunState r(context);RC::ChemicalImbalanceState(.001,-.002).Store(r.chem);auto y=r.Pack();auto driver=std::make_shared<RC::SecularEvolutionDriver>(context,RC::SecularEvolutionDriver::Mode::ChemicalOnly);EV::EvolutionSystem sys(r.ctx,r.state,r.rhs,r.layout,{driver});RC::ScaledRKF45 solver(sys,r.layout,context);RC::IntegrationStatistics stats;solver.Integrate(0,LinearTimes(2),y.data(),stats,false);
  for(auto c:{RC::BetaChannel::Npe,RC::BetaChannel::NpMu}){size_t i=RC::ChannelIndex(c);double expected=(i==0?.001:-.002)+context->W(c)*(std::pow(10+sign*2,2)-100);spin_error=std::max(spin_error,std::abs(y[i+1]-expected)/std::abs(expected));Near(y[i+1],expected,1e-6,"spin-only oracle failed");}}
 std::cout<<"SPIN_ONLY PASS max_relative_error "<<spin_error<<'\n';
 auto reaction_context=build(channels,zero);
 for(double T:{1e8,2e8})for(double sign:{-1.,1.}){RunState r(reaction_context);r.thermal.SetTinf(T);std::array<double,2> initial{sign*RC::BoltzmannMeVPerK*T*1e-5,sign*RC::BoltzmannMeVPerK*T*2e-5},slopes{};
  for(size_t i=0;i<2;++i)slopes[i]=reaction_context->Ltilde(i==0?RC::UrcaProcess::Me:RC::UrcaProcess::Mmu)*(14680/(11513*M_PI*M_PI))*std::pow(T,6)/(RC::BoltzmannErgPerK*RC::BoltzmannMeVPerK);
  double duration=1/std::max(z[0][0]*slopes[0],z[1][1]*slopes[1]);RC::ChemicalImbalanceState(initial[0],initial[1]).Store(r.chem);auto y=r.Pack();auto driver=std::make_shared<RC::SecularEvolutionDriver>(reaction_context,RC::SecularEvolutionDriver::Mode::ChemicalOnly);EV::EvolutionSystem sys(r.ctx,r.state,r.rhs,r.layout,{driver});auto observer=std::make_shared<LyapunovObserver>(z,initial);sys.AddObserver(observer);RC::ScaledRKF45 solver(sys,r.layout,reaction_context);RC::IntegrationStatistics stats;solver.Integrate(0,LinearTimes(duration),y.data(),stats,false);auto expected=LinearOracle(z,slopes,initial,duration);
  for(size_t i=0;i<2;++i){reaction_error=std::max(reaction_error,std::abs(y[i+1]-expected[i])/std::abs(expected[i]));Near(y[i+1],expected[i],1e-6,"reaction-only oracle failed");require(std::abs(y[i+1])<std::abs(initial[i]),"reaction relaxation sign");}require(observer->samples==16,"missing Lyapunov checkpoints");
  auto fault=std::make_shared<FaultDriver>(*driver,*reaction_context,"transpose_linear_map");EV::EvolutionSystem mutated(r.ctx,r.state,r.rhs,r.layout,{fault});y[0]=std::log(T/1e8);y[1]=initial[0];y[2]=initial[1];std::array<double,3> actual;mutated(0,y.data(),actual.data());bool killed=false;for(size_t i=0;i<2;++i){double correct=-z[i][0]*slopes[0]*initial[0]-z[i][1]*slopes[1]*initial[1];killed|=std::abs(actual[i+1]-correct)>1e-6*std::abs(correct);}require(killed,"M6 full map transpose escaped independent linear derivative oracle");}
 std::cout<<"MUTATION transpose_linear_map DETECTED independent small-xi derivative both signs and temperatures\n";
 std::cout<<"REACTION_ONLY PASS max_relative_error "<<reaction_error<<" ACTIVE_LYAPUNOV PASS\n";
 for(bool active:{false,true}){auto c=active?Channels(f,16,false,RC::UrcaProcessSelection{RC::UrcaProcess::Me}):no;auto context=build(c,zero);RunState r(context);std::array<double,2> initial{1e-8,2e-8};RC::ChemicalImbalanceState(initial[0],initial[1]).Store(r.chem);auto y=r.Pack();auto driver=std::make_shared<RC::SecularEvolutionDriver>(context,RC::SecularEvolutionDriver::Mode::ChemicalOnly);EV::EvolutionSystem sys(r.ctx,r.state,r.rhs,r.layout,{driver});sys.AddObserver(std::make_shared<LyapunovObserver>(z,initial));RC::ScaledRKF45 solver(sys,r.layout,context);RC::IntegrationStatistics stats;solver.Integrate(0,LinearTimes(1e12),y.data(),stats,false);if(active)require(y[2]!=initial[1],"dead channel cross Z omitted");else require(y[1]==initial[0]&&y[2]==initial[1],"all dead eta moved");}
 std::cout<<"DEAD_CHANNEL_LYAPUNOV PASS cross_Z motion retained\n";
 RunState r(reaction_context);auto driver=std::make_shared<RC::SecularEvolutionDriver>(reaction_context);auto eq=driver->Evaluate(0,r.state,r.ctx);require(eq.reaction.Rate(RC::BetaChannel::Npe)==0&&eq.reaction.Rate(RC::BetaChannel::NpMu)==0&&eq.beta.heating_erg_s==0&&eq.beta.neutrino_increment_erg_s==0&&eq.x_dot_s==eq.equilibrium_x_dot_s,"same-Ltilde RE9 failed");driver->AccumulateRHS(0,r.state,r.rhs,r.ctx);require(r.rhs.Peek(Tag::Thermal,0)==eq.equilibrium_x_dot_s,"RE9 actual coupled output failed");std::cout<<"COUPLED_RE9 PASS exact same-Ltilde equilibrium\n";
 RC::ChemicalImbalanceState(.04,-.025).Store(r.chem);auto result=driver->Evaluate(0,r.state,r.ctx);
 constexpr double mev_erg=1.602176634e-6,kb_erg=1.380649e-16,kb_mev=kb_erg/mev_erg;std::array<double,2> rates{};double equilibrium=0,increment=0,heating=0;
 for(size_t i=0;i<2;++i){double eta=i==0?.04:-.025,u=eta/(kb_mev*1e8*M_PI),L=reaction_context->Ltilde(i==0?RC::UrcaProcess::Me:RC::UrcaProcess::Mmu);double H=(14680*u+7560*std::pow(u,3)+840*std::pow(u,5)+24*std::pow(u,7))/(11513*M_PI),F=(22020*u*u+5670*std::pow(u,4)+420*std::pow(u,6)+9*std::pow(u,8))/11513;rates[i]=L/kb_erg*1e56*H;equilibrium+=L*1e64;increment+=L*1e64*F;heating+=eta*rates[i]*mev_erg;Near(result.reaction.rate_count_s[i],rates[i],2e-14,"independent SI reaction oracle");}
 Near(result.Pnet_erg_s,heating-increment-equilibrium-result.Lgamma_erg_s,2e-14,"independent thermal ledger");for(size_t i=0;i<2;++i)Near(result.eta_dot_MeV_s[i],-z[i][0]*rates[0]-z[i][1]*rates[1],2e-14,"independent cross-Z derivative");std::cout<<"SI_LEDGER PASS single conversion and cross Z\n";
 // Nonzero spin and unequal rates distinguish every sign/ordering fault.
 auto spin_context=build(channels,std::make_shared<const TestSpin>(10,-3e-10));RunState mr(spin_context);RC::ChemicalImbalanceState(.04,-.025).Store(mr.chem);auto nominal=std::make_shared<RC::SecularEvolutionDriver>(spin_context);auto y=mr.Pack();std::array<double,3> correct;EV::EvolutionSystem nominal_system(mr.ctx,mr.state,mr.rhs,mr.layout,{nominal});nominal_system(0,y.data(),correct.data());
 for(const char* name:{"OmegaDot_sign","W_sign","omit_factor2","reaction_sign","eta_sign","channel_swap","state_slot_swap","drop_cross_Z","drop_cross_Z_row0","drop_cross_Z_row1","wrong_kB_energy","omit_kB","double_kB","doubled_kB_denominator","T8_rate","omit_MeV_erg","double_MeV_erg","full_not_increment","omit_heating","double_equilibrium","legacy_placeholder"}){
  auto fault=std::make_shared<FaultDriver>(*nominal,*spin_context,name);EV::EvolutionSystem sys(mr.ctx,mr.state,mr.rhs,mr.layout,{fault});std::array<double,3> actual;sys(0,y.data(),actual.data());bool detected=false;for(size_t i=0;i<3;++i)detected|=std::abs(actual[i]-correct[i])>1e-9*std::max(std::abs(correct[i]),1e-300);require(detected,"transformed production RHS mutation escaped");std::cout<<"MUTATION "<<name<<" DETECTED boundary-transformed-production\n";}
 std::cout<<"MUTATION_ALIASES W_sign=OmegaDot_sign; eta_sign=reaction_sign; full_not_increment=double_equilibrium at this output boundary; not independent credits\n";
 auto empty=build(no,zero);RunState er(empty);RC::ChemicalImbalanceState(.04,-.025).Store(er.chem);RC::SecularEvolutionDriver ed(empty);auto e=ed.Evaluate(0,er.state,er.ctx);require(e.reaction.Rate(RC::BetaChannel::Npe)==0&&e.reaction.Rate(RC::BetaChannel::NpMu)==0&&e.reaction.EquilibriumErgPerSecond()==0&&e.beta.neutrino_increment_erg_s==0&&e.beta.heating_erg_s==0&&e.Cstar_erg_K>0&&e.x_dot_s<0,"empty-channel guard");std::cout<<"EMPTY_CHANNEL PASS finite photon cooling\n";
 // Sentinel outputs remain bit-identical for stale tokens and foreign authorities.
 r.rhs.Clear();r.rhs.AddTo(Tag::Thermal,0,123);r.rhs.AddTo(Tag::Chem,0,456);r.rhs.AddTo(Tag::Chem,1,789);
 auto unchanged=[&]{require(r.rhs.Peek(Tag::Thermal,0)==123&&r.rhs.Peek(Tag::Chem,0)==456&&r.rhs.Peek(Tag::Chem,1)==789,"dependency failure altered sentinel");};
 token->generation++;MustRefuse([&]{driver->AccumulateRHS(0,r.state,r.rhs,r.ctx);},"stale source identity before RHS");unchanged();
 EV::EvolutionSystem sys(r.ctx,r.state,r.rhs,r.layout,{driver});RC::ScaledRKF45 solver(sys,r.layout,reaction_context);double sentinel[]{123,456,789};auto packed=r.Pack();MustRefuse([&]{solver.Derivative(0,packed.data(),sentinel);},"stale dependency external derivative");require(sentinel[0]==123&&sentinel[1]==456&&sentinel[2]==789,"changed derivative sentinel");unchanged();token->generation--;
 auto saved=r.ctx.thermo;r.ctx.thermo=nullptr;MustRefuse([&]{driver->AccumulateRHS(0,r.state,r.rhs,r.ctx);},"foreign thermal owner");unchanged();r.ctx.thermo=saved;
 EV::StarContext foreign(f.central->Profile());EV::GeometryCache foreigngeo(foreign);auto savedstar=r.ctx.star;auto savedgeo=r.ctx.geo;r.ctx.star=&foreign;r.ctx.geo=&foreigngeo;MustRefuse([&]{driver->AccumulateRHS(0,r.state,r.rhs,r.ctx);},"foreign matched geometry");unchanged();r.ctx.star=savedstar;r.ctx.geo=savedgeo;
 MustRefuse([&]{reaction_context->RequireOwners(r.ctx,std::make_shared<const TestSpin>(0,0).get());},"foreign equal-valued spin");
 auto revision=f.g->Lifetime()->revision;auto domain=revision->domain;revision->domain+=" changed";MustRefuse([&]{driver->AccumulateRHS(0,r.state,r.rhs,r.ctx);},"stale support domain");unchanged();revision->domain=domain;
 EV::RHSAccumulator malformed;malformed.Configure(Tag::Thermal,1);malformed.Configure(Tag::Chem,1);malformed.AddTo(Tag::Thermal,0,123);MustRefuse([&]{driver->AccumulateRHS(0,r.state,malformed,r.ctx);},"malformed accumulator before accumulation");require(malformed.Peek(Tag::Thermal,0)==123,"partial accumulation survived");
 auto spin_token=std::make_shared<RC::RunDependencyToken>();auto prescribed=std::make_shared<const RC::PrescribedDipoleHistory>(spin_token);auto stale_context=build(channels,prescribed);RunState sr(stale_context);RC::SecularEvolutionDriver sd(stale_context);spin_token->alive=false;MustRefuse([&]{sd.AccumulateRHS(0,sr.state,sr.rhs,sr.ctx);},"same-owner stale spin");spin_token->alive=true;
 auto late_token=std::make_shared<RC::RunDependencyToken>();auto late_spin=std::make_shared<const LateInvalidationSpin>(late_token);auto late=Context(f,channels,thermal,late_spin,q,late_token);RunState lr(late);RC::SecularEvolutionDriver ld(late);lr.rhs.AddTo(Tag::Thermal,0,321);lr.rhs.AddTo(Tag::Chem,0,654);lr.rhs.AddTo(Tag::Chem,1,987);
 MustRefuse([&]{ld.AccumulateRHS(0,lr.state,lr.rhs,lr.ctx);},"late dependency failure after local evaluation before accumulation");require(lr.rhs.Peek(Tag::Thermal,0)==321&&lr.rhs.Peek(Tag::Chem,0)==654&&lr.rhs.Peek(Tag::Chem,1)==987,"late failure partially accumulated");
 // Production accepts only the privately owned timing law, never an arbitrary callback subclass.
 auto prodq=q;prodq.purpose=RC::RunPurpose::ControlledTrajectory;MustRefuse([&]{Context(f,channels,thermal,zero,prodq,token);},"unqualified callback spin refused for production");
 auto badq=q;badq.radial_resolution=40000;MustRefuse([&]{Context(f,channels,thermal,zero,badq,token);},"wrong radial qualification");
 auto badfixed=std::make_shared<FixedBaryonNumberResponse>(*f.fixed);Core::NStar foreignstar;badfixed->metadata.sources.front().star=&foreignstar;MustRefuse([&]{std::make_shared<const RC::FrozenRotochemicalRunContext>(f.z,badfixed,f.vi,f.gn,f.gv,channels,thermal,zero,q,token);},"foreign Z/W fixed owner before dereference");
 // Source mutation affects pre-run qualification, but sealed runtime values are unchanged.
 const auto& source=thermal->Sources().front();const auto before=driver->Evaluate(0,r.state,r.ctx);
 {std::ofstream out(source.path,std::ios::app);out<<'\n';}
 MustRefuse([&]{reaction_context->RequireFullCurrent();},"thermal bytes changed before run");auto after=driver->Evaluate(0,r.state,r.ctx);require(after.x_dot_s==before.x_dot_s&&after.Cstar_erg_K==before.Cstar_erg_K,"runtime reread thermal source");{std::ofstream out(source.path,std::ios::binary);out.write(source.bytes.data(),source.bytes.size());}
 reaction_context->RequireFullCurrent();std::cout<<"THERMAL_SOURCE_BYTES PASS pre-run refusal and sealed runtime\n";
 RC::FrozenContextTestAccess::Semantic({RC::BetaChannel::Npe,RC::BetaChannel::NpMu},z);
 MustRefuse([&]{RC::FrozenContextTestAccess::Semantic({RC::BetaChannel::NpMu,RC::BetaChannel::Npe},z);},"semantic Npe/NpMu swap construction gate");
 MustRefuse([&]{RC::FrozenContextTestAccess::Semantic({RC::BetaChannel::Npe,RC::BetaChannel::NpMu},{{1.}});},"semantic Z shape construction gate");
 EOS::CompOSE_Thermo::Options opts;opts.Tmin_for_derivative_MeV=0;opts.clamp_to_domain=false;
 MustRefuse([&]{RC::FrozenThermalTestAccess::Race(std::filesystem::path(source.path).parent_path(),opts,[&]{std::ofstream out(source.path,std::ios::app);out<<'\n';});},"source race during construction");
 {std::ofstream out(source.path,std::ios::binary);out.write(source.bytes.data(),source.bytes.size());}
 // Both disabled DU sources throw if sampled, despite triangle-open kinematics.
 require(RC::DirectUrcaTriangle(1e-30,1e-30,1e-30),"DU triangle control invalid");
 RC::UrcaIntegrationRequest shell;shell.chemical_domain=f.g;shell.metric=f.metric;shell.radial_partition_km=f.g->Partition();shell.domain_identity=f.g->Lifetime()->revision->domain;shell.metric_identity="qualified disconnected support analytic control";
 double edge=channels->Entry(RC::UrcaProcess::Mmu).support.front().right_km;std::vector<RC::UrcaSupportInterval> intervals{{.2*edge,.3*edge},{.7*edge,.9*edge}};
 for(auto p:{RC::UrcaProcess::Me,RC::UrcaProcess::Mmu})shell.normalizations.emplace_back(p,[intervals,p](double x){require(std::any_of(intervals.begin(),intervals.end(),[&](auto a){return x>=a.left_km&&x<=a.right_km;}),"closed interior activated");return p==RC::UrcaProcess::Me?1e-51:2e-51;},intervals,"disconnected controlled test support");
 for(auto p:{RC::UrcaProcess::De,RC::UrcaProcess::Dmu})shell.normalizations.emplace_back(p,[](double)->double{throw std::runtime_error("disabled DU sampled");},intervals,"triangle-open disabled DU test");
 auto shellc=std::make_shared<const RC::GlobalUrcaChannelCoefficient>(RC::GlobalUrcaChannelCoefficient::Compute(shell));auto shellcontext=build(shellc,zero);RunState shr(shellcontext);RC::ChemicalImbalanceState(.04,-.025).Store(shr.chem);RC::SecularEvolutionDriver shd(shellcontext);auto sh=shd.Evaluate(0,shr.state,shr.ctx);
 require(sh.reaction.equilibrium_erg_s[0]==0&&sh.reaction.equilibrium_erg_s[1]==0,"disabled DU coupled power");for(size_t i=0;i<2;++i){auto p=i==0?RC::UrcaProcess::Me:RC::UrcaProcess::Mmu;Near(sh.reaction.rate_count_s[i],result.reaction.rate_count_s[i]*shellcontext->Ltilde(p)/reaction_context->Ltilde(p),2e-14,"disconnected coupled coefficient response");}
 auto wrong_partition=shell;wrong_partition.radial_partition_km.back()=std::nextafter(wrong_partition.radial_partition_km.back(),0.);MustRefuse([&]{RC::GlobalUrcaChannelCoefficient::Compute(wrong_partition);},"changed partition endpoint");
 for(bool profile_fault:{false,true}){auto bad=q;auto original_path=profile_fault?q.profile_path:q.certificate_path;auto copy=std::filesystem::path(source.path).parent_path()/(profile_fault?"altered-profile.txt":"altered-certificate.txt");{std::ofstream out(copy,std::ios::binary);out<<RC::FrozenSource::Read(original_path)<<" changed";}if(profile_fault)bad.profile_path=copy.string();else bad.certificate_path=copy.string();MustRefuse([&]{Context(f,channels,thermal,zero,bad,token);},profile_fault?"altered qualified profile hash":"altered qualification certificate hash");}
 std::cout<<"COUPLED_SUPPORT PASS explicit disjoint intervals and disabled DU\n";
 std::cout<<"ARCHITECTURE_GUARDS PASS exercised downstream sentinel and ownership controls\n";
}
