// Diagnostic-only local replays. No authoritative trajectory entry point is called.
#define PHASE6A1_CAMPAIGN_SUPPORT_ONLY
#include "adr0017_production_qualification.cpp"
#include <gsl/gsl_version.h>
#include <cstring>

namespace Probe {
using State=BNV::CheckpointState;
// Standard explicit-instantiation access exception: read-only test introspection.
// No production header, layout, cache payload, or lookup implementation is changed.
struct CacheTag {friend auto member(CacheTag);};
template<class Tag,auto M> struct Access {friend auto member(Tag){return M;}};
template struct Access<CacheTag,&EV::StarContext::m_cv_cache>;
const auto& Cache(const EV::StarContext& s){return s.*member(CacheTag{});}
double TM(double x){return 1e8*std::exp(x)*RC::BoltzmannMeVPerK;}
struct Case {int index;double left,obs,right;State yl,yr;int knots;std::size_t lo,hi;};
std::vector<Case> ReadCases(const std::filesystem::path& p){std::ifstream f(p);std::string line;std::getline(f,line);std::vector<Case> v;while(std::getline(f,line)){auto s=Split(line);require(s.size()==13,"case schema");v.push_back({std::stoi(s[0]),std::stod(s[1]),std::stod(s[2]),std::stod(s[3]),{{std::stod(s[4]),std::stod(s[5]),std::stod(s[6])}},{{std::stod(s[7]),std::stod(s[8]),std::stod(s[9])}},std::stoi(s[10]),std::stoull(s[11]),std::stoull(s[12])});}return v;}
struct Attempt {double t,h;State before,after,err;int rc;bool accepted=false;};
struct Trace {
 ActualCheckpointContext context;EV::DriverContext owner;std::ofstream* rhs=nullptr;
 std::vector<Attempt> attempts;const gsl_odeiv2_step_type* method=nullptr;double start=0;
 std::exception_ptr failure;std::size_t calls=0;
 Trace(std::shared_ptr<const BNV::FrozenControlledBnvRunContext> p):context(p,1),owner(p->OrdinaryContext()->DriverContext()){}
 static int Fn(double t,const double*y,double*d,void*vp) noexcept{
  auto& a=*static_cast<Trace*>(vp);try{State in{{y[0],y[1],y[2]}},out{};const auto before=Cache(*a.owner.star).last_i;a.context.Derivative(t,in,out);std::copy(out.begin(),out.end(),d);++a.calls;
   if(a.rhs)*a.rhs<<a.attempts.size()<<'\t'<<t<<'\t'<<y[0]<<'\t'<<y[1]<<'\t'<<y[2]<<'\t'<<out[0]<<'\t'<<out[1]<<'\t'<<out[2]<<'\t'<<TM(y[0])<<'\t'<<before<<'\t'<<Cache(*a.owner.star).last_i<<'\n';return GSL_SUCCESS;
  }catch(...){a.failure=std::current_exception();return GSL_EBADFUNC;}}
};
thread_local Trace* active=nullptr;
int Apply(void*opaque,std::size_t n,double t,double h,double*y,double*err,const double*di,double*doo,const gsl_odeiv2_system*sys){
 auto& a=*active;Attempt row{t,h,{{y[0],y[1],y[2]}},{},{},0};a.attempts.push_back(row);
 int rc=a.method->apply(opaque,n,t,h,y,err,di,doo,sys);auto& q=a.attempts.back();q.after={{y[0],y[1],y[2]}};q.err={{err[0],err[1],err[2]}};q.rc=rc;return rc;
}
struct Result {State y{};double t=0;std::size_t accepted=0,rejected=0,rhs=0;int rc=0;double wall=0;};
struct Lab {
 std::shared_ptr<const BNV::FrozenControlledBnvRunContext> physics;
 std::filesystem::path output;std::ofstream summary;std::size_t solve_count=0;
 Lab(std::shared_ptr<const BNV::FrozenControlledBnvRunContext> p,std::filesystem::path o):physics(p),output(o),summary(o/"solutions.tsv"){
  summary<<std::setprecision(17)<<"name\tmethod\tlevel\tt_start\tt_end\tx\teta_e\teta_mu\taccepted\trejected\trhs\trc\twall_s\n";
 }
 Result Solve(std::string name,const gsl_odeiv2_step_type* method,int level,double start,double end,State initial,bool trace=true,bool full=true){
  ++solve_count;require(end>start,"positive bounded replay interval");Trace a(physics);a.method=method;
  if(full)a.context.RequireFullCurrent();
  std::ofstream rhs;if(trace){rhs.open(output/(name+".rhs.tsv"));rhs<<std::setprecision(17)<<"attempt\tt\tx\teta_e\teta_mu\txdot\tedot\tmdot\tT_MeV\thint_before\tcell_after\n";a.rhs=&rhs;}
  auto type=*method;type.apply=Apply;
  // O1/O2 literals match production configuration exactly (pow is not authority).
  const double rels[]={0,1e-12,1e-13,1e-14,1e-15};
  const State ats[]={{},{1e-17,1e-23,1e-23},{1e-18,1e-24,1e-24},{1e-19,1e-25,1e-25},{1e-20,1e-26,1e-26}};
  auto* step=gsl_odeiv2_step_alloc(&type,3);auto* control=gsl_odeiv2_control_scaled_new(1.,rels[level],1.,0.,ats[level].data(),3);auto* evolve=gsl_odeiv2_evolve_alloc(3);require(step&&control&&evolve,"GSL allocation");
  gsl_odeiv2_system system{Trace::Fn,nullptr,3,&a};double t=start,h=end-start;State y=initial;std::size_t count=0;int rc=0;auto begin=Clock::now();active=&a;
  while(t<end){if(++count>20000){rc=GSL_EMAXITER;break;}const auto before=a.attempts.size();rc=gsl_odeiv2_evolve_apply(evolve,control,step,&system,&t,end,&h,y.data());if(a.failure)std::rethrow_exception(a.failure);if(rc!=0)break;require(a.attempts.size()>before,"no GSL attempt");a.attempts.back().accepted=true;}
  Result result{y,t,0,static_cast<std::size_t>(evolve->failed_steps),a.calls,rc,std::chrono::duration<double>(Clock::now()-begin).count()};for(auto&q:a.attempts)result.accepted+=q.accepted;
  if(trace){std::ofstream f(output/(name+".steps.tsv"));f<<std::setprecision(17)<<"attempt\tt\th\taccepted\trc\tx_before\tx_after\teta_e_before\teta_e_after\teta_mu_before\teta_mu_after\terr_x\terr_eta_e\terr_eta_mu\n";std::size_t i=0;for(auto&q:a.attempts)f<<++i<<'\t'<<q.t<<'\t'<<q.h<<'\t'<<q.accepted<<'\t'<<q.rc<<'\t'<<q.before[0]<<'\t'<<q.after[0]<<'\t'<<q.before[1]<<'\t'<<q.after[1]<<'\t'<<q.before[2]<<'\t'<<q.after[2]<<'\t'<<q.err[0]<<'\t'<<q.err[1]<<'\t'<<q.err[2]<<'\n';}
  summary<<name<<'\t'<<method->name<<'\t'<<level<<'\t'<<start<<'\t'<<result.t<<'\t'<<y[0]<<'\t'<<y[1]<<'\t'<<y[2]<<'\t'<<result.accepted<<'\t'<<result.rejected<<'\t'<<result.rhs<<'\t'<<rc<<'\t'<<result.wall<<'\n';summary.flush();
  gsl_odeiv2_evolve_free(evolve);gsl_odeiv2_control_free(control);gsl_odeiv2_step_free(step);active=nullptr;
  if(full)a.context.RequireFullCurrent();return result;
 }
};
void StateLine(std::ostream&f,const State&y){for(double v:y)f<<'\t'<<v;}
}

#ifndef OBS32_SUPPORT_ONLY
int main(int argc,char**argv){try{
 using namespace Probe;gsl_set_error_handler_off();std::cout<<std::setprecision(17)<<std::unitbuf;
 require(argc==10,"profile certificate thermal work entry frozen coefficients cases output");Inputs in;in.profile=argv[1];in.certificate=argv[2];in.thermal=argv[3];in.work=argv[4];in.entry=argv[5];in.frozen=argv[6];in.coefficients=argv[7];in.output=argv[9];require(!std::filesystem::exists(in.work)&&!std::filesystem::exists(in.output),"fresh local context/output required");std::filesystem::create_directories(in.output);
 auto cases=ReadCases(argv[8]);const auto c=cases.front();require(c.index==32&&c.obs==61635937500.&&c.knots==1,"wrong failure bracket");
 auto begin=Clock::now();auto physics=BuildContext(in);std::cout<<"CONTEXT_READY wall "<<std::chrono::duration<double>(Clock::now()-begin).count()<<'\n';
 // Production reconstruction first. A mismatch prevents every further experiment.
 BNV::AcceptedStepRecord record;record.ordinal=1;record.t_left_s=c.left;record.t_right_s=c.right;record.y_left=c.yl;record.y_right=c.yr;record.cstar_knots_crossed=1;
 BNV::ObservationBracket bracket{32,c.obs,0,1,c.left,c.right,false};
 BNV::Rk8pdCheckpointReconstructor recon;std::size_t id=0;
 auto q=recon.Reconstruct(bracket,{record},[&]{return std::make_unique<ActualCheckpointContext>(physics,++id);});
 const State expected1{{0.12472490938574241,-2.0287605965623057e-7,-5.0871537810150439e-7}},expected2{{0.12472490938614615,-2.0287605965622805e-7,-5.0871537810149634e-7}};
 require(q.O1_state==expected1&&q.O2_state==expected2&&q.status==BNV::CheckpointStatus::NumericallyUnresolved,"archived failure did not reproduce; stop provenance investigation");
 std::ofstream reproduction(in.output/"reproduction.tsv");reproduction<<std::setprecision(17)<<"component\tO1\tO2\td\tD_O1\tF_i\tU\td_over_D\tU_over_point2F\n";for(int i=0;i<3;++i)reproduction<<i<<'\t'<<q.O1_state[i]<<'\t'<<q.O2_state[i]<<'\t'<<q.qualification.d_O[i]<<'\t'<<q.qualification.D_O1[i]<<'\t'<<q.qualification.F_i[i]<<'\t'<<q.qualification.U_O[i]<<'\t'<<q.qualification.d_O[i]/q.qualification.D_O1[i]<<'\t'<<q.qualification.U_O[i]/(.20*q.qualification.F_i[i])<<'\n';reproduction.close();std::cout<<"ARCHIVED_FAILURE_REPRODUCED_EXACTLY\n";
 auto ctx=physics->OrdinaryContext()->DriverContext();const auto& cache=Cache(*ctx.star);require(cache.loaded&&cache.Tinf_MeV.size()==160,"cache topology");
 {std::ifstream f(std::filesystem::path(argv[8]).parent_path()/"archived-prefix.tsv");std::string line;std::getline(f,line);auto header=Split(line);auto col=[&](std::string key){auto it=std::find(header.begin(),header.end(),key);require(it!=header.end(),"missing archived cache field");return std::size_t(it-header.begin());};const auto tc=col("Tinf_K"),cc=col("Cstar_erg_K");std::size_t n=0;while(std::getline(f,line)){auto row=Split(line);double T=std::stod(row.at(tc)),C=std::stod(row.at(cc));require(ctx.star->HeatCapacityStar_Tinf(T*RC::BoltzmannMeVPerK,*ctx.thermo,ctx.geo)==C,"archived Cstar value mismatch");++n;}require(n==32,"archived prefix rows");std::cout<<"ARCHIVED_CSTAR_PREFIX_EXACT 32\n";}

 std::ofstream cf(in.output/"cstar-cache.tsv");cf<<std::setprecision(17)<<"i\tT_MeV\tT_K\tCstar\tT_hex\tC_hex\n";for(std::size_t i=0;i<cache.Tinf_MeV.size();++i)cf<<i<<'\t'<<cache.Tinf_MeV[i]<<'\t'<<cache.Tinf_MeV[i]/RC::BoltzmannMeVPerK<<'\t'<<cache.C_star[i]<<'\t'<<std::hexfloat<<cache.Tinf_MeV[i]<<'\t'<<cache.C_star[i]<<std::defaultfloat<<'\n';cf.close();
 std::vector<std::size_t> knots;for(std::size_t i=1;i<159;++i)if(TM(c.yl[0])<cache.Tinf_MeV[i]&&cache.Tinf_MeV[i]<TM(c.yr[0]))knots.push_back(i);require(knots.size()==1,"actual cache knot count differs");auto ki=knots[0];double knot=cache.Tinf_MeV[ki];std::cout<<"KNOT index "<<ki<<" T_MeV "<<knot<<" T_K "<<knot/RC::BoltzmannMeVPerK<<'\n';
 Lab lab(physics,in.output);std::vector<Result> unsplit;
 for(int l=1;l<=4;++l){auto r=lab.Solve("unsplit-O"+std::to_string(l),gsl_odeiv2_step_rk8pd,l,c.left,c.obs,c.yl);unsplit.push_back(r);if(l<=3)require(r.rc==0,"required unsplit level failed");if(l==1)require(r.y==q.O1_state,"instrumented O1 changed result");if(l==2)require(r.y==q.O2_state,"instrumented O2 changed result");std::cout<<"UNSPLIT O"<<l<<" rc "<<r.rc<<" accepted "<<r.accepted<<" rejected "<<r.rejected<<'\n';}
 // Bounded event root: evaluate the local O3 IVP from the preserved left state.
 // Full currentness encloses the whole root experiment; every trial uses fresh GSL
 // and state objects, and ordinary cheap currentness remains in every RHS call.
 require(TM(c.yl[0])<knot&&TM(unsplit[2].y[0])>knot,"knot not inside left-to-observation interval");physics->RequireFullCurrent();double lo=c.left,hi=c.obs;Result root;std::ofstream roots(in.output/"event-root.tsv");roots<<std::setprecision(17)<<"iteration\tlo\thi\tmid\tT_minus_knot_MeV\n";
 for(int i=0;i<64;++i){double mid=lo+(hi-lo)/2.;if(mid==lo||mid==hi)break;root=lab.Solve("root-"+std::to_string(i),gsl_odeiv2_step_rk8pd,3,c.left,mid,c.yl,false,false);require(root.rc==0,"event solve pathology");double f=TM(root.y[0])-knot;roots<<i<<'\t'<<lo<<'\t'<<hi<<'\t'<<mid<<'\t'<<f<<'\n';if(f<0)lo=mid;else hi=mid;}
 const double tk=lo+(hi-lo)/2.;physics->RequireFullCurrent();std::cout<<"EVENT_BRACKET "<<lo<<' '<<hi<<" chosen "<<tk<<'\n';
 std::ofstream event(in.output/"event.tsv");event<<std::setprecision(17)<<"knot_index\tT_MeV\tT_K\tt_lo\tt_hi\tt_chosen\n"<<ki<<'\t'<<knot<<'\t'<<knot/RC::BoltzmannMeVPerK<<'\t'<<lo<<'\t'<<hi<<'\t'<<tk<<'\n';event.close();
 for(int l=1;l<=3;++l){auto a=lab.Solve("split-left-O"+std::to_string(l),gsl_odeiv2_step_rk8pd,l,c.left,tk,c.yl);require(a.rc==0,"split first leg");auto b=lab.Solve("split-right-O"+std::to_string(l),gsl_odeiv2_step_rk8pd,l,tk,c.obs,a.y);require(b.rc==0,"split second leg");std::cout<<"SPLIT O"<<l<<" event_T_residual_MeV "<<TM(a.y[0])-knot<<'\n';}
 for(int l=2;l<=3;++l){auto a=lab.Solve("independent-unsplit-O"+std::to_string(l),gsl_odeiv2_step_rkf45,l,c.left,c.obs,c.yl);require(a.rc==0,"independent unsplit");auto b=lab.Solve("independent-left-O"+std::to_string(l),gsl_odeiv2_step_rkf45,l,c.left,tk,c.yl);require(b.rc==0,"independent first leg");auto d=lab.Solve("independent-right-O"+std::to_string(l),gsl_odeiv2_step_rkf45,l,tk,c.obs,b.y);require(d.rc==0,"independent second leg");}
 for(std::size_t i=1;i<cases.size();++i){auto a=cases[i];require(a.knots==0,"control contains knot");for(int l=1;l<=3;++l){auto r=lab.Solve("smooth-"+std::to_string(a.index)+"-O"+std::to_string(l),gsl_odeiv2_step_rk8pd,l,a.left,a.obs,a.yl);require(r.rc==0,"smooth control pathology");}}
 // Cell ownership and complete RHS audit. Read actual cache hint, never change it.
 std::ofstream cells(in.output/"cell-ownership.tsv");cells<<std::setprecision(17)<<"warm_side\tpoint\tT_MeV\tCstar\thint_before\tcell_after\n";
 for(int side:{-1,1})for(int p:{-1,0,1}){double warm=cache.Tinf_MeV[ki+side];ctx.star->HeatCapacityStar_Tinf(warm,*ctx.thermo,ctx.geo);double t=p==0?knot:std::nextafter(knot,p<0?0.:INFINITY);auto before=Cache(*ctx.star).last_i;double v=ctx.star->HeatCapacityStar_Tinf(t,*ctx.thermo,ctx.geo);cells<<side<<'\t'<<p<<'\t'<<t<<'\t'<<v<<'\t'<<before<<'\t'<<Cache(*ctx.star).last_i<<'\n';}
 std::ofstream rhsout(in.output/"rhs-continuity.tsv");rhsout<<std::setprecision(17)<<"offset_x\tx\tT_MeV\tcell\txdot\teta_dot_e\teta_dot_mu\tCstar\tPnet\tPdir\tLH\tDeltaLnu\tLnu_eq\tLnu_full\tLgamma\tsigma_e\tsigma_mu\tR_e\tR_mu\n";
 double xk=std::log((knot/RC::BoltzmannMeVPerK)/1e8);auto ev=lab.Solve("event-final-O3",gsl_odeiv2_step_rk8pd,3,c.left,tk,c.yl,false);require(ev.rc==0,"event final");
 for(double offset:{-1e-6,-1e-7,-1e-8,-1e-9,-1e-10,-1e-12,0.,1e-12,1e-10,1e-9,1e-8,1e-7,1e-6}){State y=ev.y;y[0]=xk+offset;ActualCheckpointContext isolated(physics,++id);State dy;isolated.Derivative(tk,y,dy);auto d=isolated.EvaluateDiagnostics(tk,y);rhsout<<offset<<'\t'<<y[0]<<'\t'<<TM(y[0])<<'\t'<<Cache(*ctx.star).last_i;StateLine(rhsout,dy);rhsout<<'\t'<<d.Cstar_erg_K<<'\t'<<d.Pnet_erg_s<<'\t'<<d.P_dir_actual_erg_s<<'\t'<<d.LH_erg_s<<'\t'<<d.DeltaLnu_erg_s<<'\t'<<d.Lnu_eq_erg_s<<'\t'<<d.Lnu_full_erg_s<<'\t'<<d.Lgamma_erg_s<<'\t'<<d.sigma_count_s[0]<<'\t'<<d.sigma_count_s[1]<<'\t'<<d.R_count_s[0]<<'\t'<<d.R_count_s[1]<<'\n';}
 // A separately constructed cold thermal cache must reproduce all 160 payload entries.
 EV::StarContext cold(*ctx.star->Provenance().source);cold.HeatCapacityStar_Tinf(knot,*ctx.thermo,ctx.geo);
 require(Cache(cold).Tinf_MeV==cache.Tinf_MeV&&Cache(cold).C_star==cache.C_star,"cold cache payload changed");std::cout<<"COLD_CACHE_PAYLOAD_EXACT 160\n";
 // Route/history check: identical input in fresh disposable states after opposite warms.
 std::ofstream route(in.output/"context-route.tsv");route<<std::setprecision(17)<<"sample\tcomponent\tmain_route\treplay_cold\treplay_warm\n";
 RunState state(physics->OrdinaryContext());auto driver=std::make_shared<BNV::ControlledBnvSecularDriver>(physics,BNV::ControlledBnvSecularDriver::Mode::Coupled);EV::EvolutionSystem main_route(state.ctx,state.state,state.rhs,state.layout,{driver});
 int sample=0;for(State y:{c.yl,c.yr,q.O1_state,q.O2_state,ev.y}){State a,b,d;require(main_route(c.obs,y.data(),a.data())==0,"main-style RHS");ActualCheckpointContext fresh(physics,++id);fresh.Derivative(c.obs,y,b);ctx.star->HeatCapacityStar_Tinf(cache.Tinf_MeV[ki+1],*ctx.thermo,ctx.geo);fresh.Derivative(c.obs,y,d);require(a==b&&a==d,"context/hint route mismatch");for(int i=0;i<3;++i)route<<sample<<'\t'<<i<<'\t'<<a[i]<<'\t'<<b[i]<<'\t'<<d[i]<<'\n';++sample;}
 physics->RequireFullCurrent();std::ofstream done(in.output/"completion.tsv");done<<std::setprecision(17)<<"key\tvalue\nlocal_solver_invocations\t"<<lab.solve_count+2<<"\nfull_trajectories\t0\nproduction_reproduction\tEXACT_FAIL\nwall_s\t"<<std::chrono::duration<double>(Clock::now()-begin).count()<<"\ngsl\t"<<gsl_version<<'\n';std::cout<<"BOUNDED_FORENSICS_COMPLETE local_solves "<<lab.solve_count+2<<'\n';
 return 0;
 }catch(const std::exception&e){std::cerr<<"STOP "<<e.what()<<'\n';return 1;}}

#endif
