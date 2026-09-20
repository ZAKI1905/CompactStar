#define PHASE6A1_BA12R_ULTRA_SUPPORT_ONLY
#include "coupled_trajectory.cpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

namespace
{
constexpr std::size_t ExpectedPositiveObservations=240;
constexpr std::size_t ExpectedStrictInterior=239;
constexpr std::size_t ExpectedTotalSolves=1434;
constexpr std::size_t MaximumApplyCalls=100000;

struct MatrixRow
{
    std::size_t observation=0,left_index=0,right_index=0;
    double t_obs=0,t_left=0,t_right=0;
    std::array<double,3> left{},right{};
    char category='A';
    int knots=0;
    bool deep=false,exact=false,strict=false;
    std::size_t local_integrations=0;
};

struct SolveStatistics
{
    std::size_t accepted=0,rejected=0,rhs=0;
    double wall_s=0,user_s=0,sys_s=0;
};

struct MethodResult
{
    std::array<double,3> state{};
    std::array<double,3> f_left{{0,0,0}},f_right{{0,0,0}};
    bool has_endpoint_rhs=false;
    SolveStatistics solve;
};

struct Inputs
{
    // profile_root is the EOS/profile DIRECTORY (freegas.tsv, model.txt, profile.tsv); it is never a trajectory TSV.
    std::filesystem::path matrix,profile_root,certificate,thermal,work,entry,frozen,coefficients,output,predeclaration,
      pretrajectory,oracle_flag,authorized_solves,execution_ledger;
};

std::vector<std::string> Split(const std::string& line)
{
    std::vector<std::string> values;std::string item;std::istringstream in(line);
    while(std::getline(in,item,'\t'))values.push_back(item);
    return values;
}

std::vector<MatrixRow> ReadMatrix(const std::filesystem::path& path)
{
    std::ifstream in(path);if(!in)throw std::runtime_error("cannot open solve matrix");
    std::string line;if(!std::getline(in,line))throw std::runtime_error("empty solve matrix");
    const std::vector<std::string> expected{
      "observation_index","t_obs_s","left_endpoint_index","t_left_s","x_left","eta_e_left_MeV","eta_mu_left_MeV",
      "right_endpoint_index","t_right_s","x_right","eta_e_right_MeV","eta_mu_right_MeV","cstar_category",
      "cstar_knots_crossed","deep_interior","exact_endpoint","strict_interior","linear_eval","hermite_eval",
      "replay1","replay2","replay1_repeat","replay2_repeat","oracle1","oracle2","local_integrations"};
    if(Split(line)!=expected)throw std::runtime_error("solve matrix header changed");
    std::vector<MatrixRow> rows;std::size_t total=0,a=0,b=0,c=0,deep=0,exact=0,strict=0;
    while(std::getline(in,line))
    {
        if(line.empty())continue;const auto v=Split(line);if(v.size()!=expected.size())throw std::runtime_error("solve matrix row shape changed");
        MatrixRow r;r.observation=std::stoull(v[0]);r.t_obs=std::stod(v[1]);r.left_index=std::stoull(v[2]);r.t_left=std::stod(v[3]);
        r.left={std::stod(v[4]),std::stod(v[5]),std::stod(v[6])};r.right_index=std::stoull(v[7]);r.t_right=std::stod(v[8]);
        r.right={std::stod(v[9]),std::stod(v[10]),std::stod(v[11])};if(v[12].size()!=1)throw std::runtime_error("invalid Cstar category");
        r.category=v[12][0];r.knots=std::stoi(v[13]);r.deep=std::stoi(v[14])!=0;r.exact=std::stoi(v[15])!=0;r.strict=std::stoi(v[16])!=0;
        for(std::size_t i=17;i<=24;++i)if(std::stoi(v[i])!=(r.strict?1:0) && i>=19)throw std::runtime_error("solve applicability changed");
        if(std::stoi(v[17])!=1||std::stoi(v[18])!=1)throw std::runtime_error("algebraic applicability changed");
        r.local_integrations=std::stoull(v[25]);
        if(r.observation!=rows.size()+1||!(r.t_left<r.t_obs||r.exact)||!(r.t_obs<=r.t_right)||r.right_index!=r.left_index+1)
            throw std::runtime_error("solve matrix bracket changed");
        if(r.strict!=(r.t_left<r.t_obs&&r.t_obs<r.t_right)||r.exact!=(r.t_obs==r.t_right)||r.local_integrations!=(r.strict?6u:0u))
            throw std::runtime_error("solve matrix semantics changed");
        total+=r.local_integrations;deep+=r.deep;exact+=r.exact;strict+=r.strict;a+=r.category=='A';b+=r.category=='B';c+=r.category=='C';rows.push_back(r);
    }
    if(rows.size()!=ExpectedPositiveObservations||a!=237||b!=3||c!=0||deep!=81||exact!=1||strict!=ExpectedStrictInterior||total!=ExpectedTotalSolves)
        throw std::runtime_error("solve matrix authority mismatch");
    return rows;
}

MatrixRow InitialRow()
{
    MatrixRow r;r.category='A';r.exact=true;return r;
}

double CpuSeconds(const timeval& value){return value.tv_sec+1.0e-6*value.tv_usec;}

std::string Sha256(const std::filesystem::path& path)
{
    const std::string command="/usr/bin/shasum -a 256 '"+path.string()+"'";
    std::unique_ptr<FILE,decltype(&pclose)> pipe(popen(command.c_str(),"r"),pclose);
    if(!pipe)throw std::runtime_error("cannot start SHA-256 tool");char buffer[256]{};
    if(!fgets(buffer,sizeof(buffer),pipe.get()))throw std::runtime_error("cannot read SHA-256");
    std::string hash(buffer);hash=hash.substr(0,hash.find_first_of(" \t\r\n"));
    if(hash.size()!=64)throw std::runtime_error("invalid SHA-256 result");return hash;
}

void AtomicWrite(const std::filesystem::path& path,const std::string& contents)
{
    const auto temporary=path.string()+".tmp."+std::to_string(getpid());
    {std::ofstream out(temporary);if(!out)throw std::runtime_error("cannot create temporary output");out<<contents;if(!out)throw std::runtime_error("output write failed");}
    std::filesystem::rename(temporary,path);
}

// Recovery correction: the profile/EOS root is a typed directory argument authenticated by the same
// byte hashes production qualification enforces (FrozenRotochemicalRunContext.hpp:89-92). A trajectory
// TSV, any regular file, a missing file, or a hash mismatch refuses before any context is built.
void ValidateProfileRoot(const std::filesystem::path& root)
{
    if(std::filesystem::is_regular_file(root))throw std::runtime_error("profile root must be the EOS/profile directory, not a trajectory/TSV file: "+root.string());
    if(!std::filesystem::is_directory(root))throw std::runtime_error("profile root is not a directory: "+root.string());
    const std::array<std::pair<const char*,const char*>,3> expected{{
      {"profile.tsv","e9cd03b0b8449806f6c9883d75de1d3dff0cf1481675efc45a56519655d40890"},
      {"model.txt","3ea70de79e15b70c5a6d68f48335d18047ff80e60b55a9acdb78084e9be4d6d4"},
      {"freegas.tsv","7cd44c92e1e7206e0e68e3fed7e3f0ca68e79ab4517d02b96ff78b9be23d3f1a"}}};
    for(const auto& [name,hash]:expected)
    {
        const auto file=root/name;if(!std::filesystem::is_regular_file(file))throw std::runtime_error(std::string("profile root lacks ")+name);
        if(Sha256(file)!=hash)throw std::runtime_error(std::string("profile root ")+name+" differs from the authenticated Structure-1 qualification bytes");
    }
}

struct EvaluationEnvironment
{
    RunState state;
    std::shared_ptr<BNV::ControlledBnvSecularDriver> driver;
    EV::EvolutionSystem system;
    explicit EvaluationEnvironment(const std::shared_ptr<const BNV::FrozenControlledBnvRunContext>& context)
      :state(context->OrdinaryContext()),driver(std::make_shared<BNV::ControlledBnvSecularDriver>(context,BNV::ControlledBnvSecularDriver::Mode::Coupled)),
       system(state.ctx,state.state,state.rhs,state.layout,{driver}){}
    std::array<double,3> Derivative(double t,const std::array<double,3>& y)
    {
        std::array<double,3> out{};if(system(t,y.data(),out.data())!=0)throw std::runtime_error("RHS failure");
        for(double value:out)if(!std::isfinite(value))throw std::runtime_error("nonfinite RHS");return out;
    }
    void Save(double t,const std::array<double,3>& y,const std::filesystem::path& path)
    {
        EV::UnpackStateVector(state.state,state.layout,y.data());
        {Capture capture(path,*driver);capture.Save(t,state.state,state.ctx);if(capture.rows!=1)throw std::runtime_error("diagnostic row count changed");}
    }
};

struct CallbackData
{
    EvaluationEnvironment* environment=nullptr;std::size_t rhs=0;std::exception_ptr failure;
};

int Callback(double t,const double* y,double* out,void* pointer)noexcept
{
    auto& data=*static_cast<CallbackData*>(pointer);
    try
    {
        std::array<double,3> state{{y[0],y[1],y[2]}};const auto value=data.environment->Derivative(t,state);std::copy(value.begin(),value.end(),out);++data.rhs;return GSL_SUCCESS;
    }
    catch(...){data.failure=std::current_exception();return GSL_EBADFUNC;}
}

MethodResult Integrate(EvaluationEnvironment& environment,const MatrixRow& row,
    const gsl_odeiv2_step_type* type,double rtol,const std::array<double,3>& atol)
{
    MethodResult result;result.state=row.left;
    std::unique_ptr<gsl_odeiv2_step,decltype(&gsl_odeiv2_step_free)> step(gsl_odeiv2_step_alloc(type,3),gsl_odeiv2_step_free);
    std::unique_ptr<gsl_odeiv2_control,decltype(&gsl_odeiv2_control_free)> control(
      gsl_odeiv2_control_scaled_new(1.0,rtol,1.0,0.0,atol.data(),3),gsl_odeiv2_control_free);
    std::unique_ptr<gsl_odeiv2_evolve,decltype(&gsl_odeiv2_evolve_free)> evolve(gsl_odeiv2_evolve_alloc(3),gsl_odeiv2_evolve_free);
    if(!step||!control||!evolve)throw std::bad_alloc();CallbackData callback{&environment};gsl_odeiv2_system system{Callback,nullptr,3,&callback};
    double t=row.t_left,h=row.t_obs-row.t_left;if(!(h>0))throw std::runtime_error("invalid local initial step");
    const auto wall_start=std::chrono::steady_clock::now();rusage before{},after{};getrusage(RUSAGE_SELF,&before);
    while(t<row.t_obs)
    {
        if(++result.solve.accepted>MaximumApplyCalls)throw std::runtime_error("local apply-call budget exceeded");
        const double previous=t;const int rc=gsl_odeiv2_evolve_apply(evolve.get(),control.get(),step.get(),&system,&t,row.t_obs,&h,result.state.data());
        result.solve.rejected=evolve->failed_steps;if(callback.failure)std::rethrow_exception(callback.failure);
        if(rc!=GSL_SUCCESS)throw std::runtime_error(std::string("local GSL failure: ")+gsl_strerror(rc));
        if(!(t>previous)||!std::isfinite(t))throw std::runtime_error("local integration failed progress");
        const double T=1.0e8*std::exp(result.state[0]);if(!(T>0)||!std::isfinite(T)||!(T*RC::BoltzmannMeVPerK>1.0e-5&&T*RC::BoltzmannMeVPerK<1.0))
            throw std::runtime_error("local accepted state outside thermal domain");
    }
    getrusage(RUSAGE_SELF,&after);result.solve.wall_s=std::chrono::duration<double>(std::chrono::steady_clock::now()-wall_start).count();
    result.solve.user_s=CpuSeconds(after.ru_utime)-CpuSeconds(before.ru_utime);result.solve.sys_s=CpuSeconds(after.ru_stime)-CpuSeconds(before.ru_stime);result.solve.rhs=callback.rhs;
    if(t!=row.t_obs)throw std::runtime_error("local target not reached exactly");return result;
}

bool StartsWith(const std::string& value,const std::string& prefix){return value.rfind(prefix,0)==0;}
bool IsIntegration(const std::string& method){return StartsWith(method,"oracle")||StartsWith(method,"replay");}

MethodResult EvaluateMethod(EvaluationEnvironment& environment,const MatrixRow& row,const std::string& method)
{
    MethodResult result;
    if(row.observation==0||row.exact){result.state=row.observation==0?std::array<double,3>{{0,0,0}}:row.right;return result;}
    const double s=(row.t_obs-row.t_left)/(row.t_right-row.t_left),H=row.t_right-row.t_left;
    if(StartsWith(method,"linear"))
    {
        for(std::size_t i=0;i<3;++i)result.state[i]=(1.0-s)*row.left[i]+s*row.right[i];return result;
    }
    if(StartsWith(method,"hermite"))
    {
        result.f_left=environment.Derivative(row.t_left,row.left);result.f_right=environment.Derivative(row.t_right,row.right);result.has_endpoint_rhs=true;
        const double h00=2*s*s*s-3*s*s+1,h10=s*s*s-2*s*s+s,h01=-2*s*s*s+3*s*s,h11=s*s*s-s*s;
        for(std::size_t i=0;i<3;++i)result.state[i]=h00*row.left[i]+h10*H*result.f_left[i]+h01*row.right[i]+h11*H*result.f_right[i];
        return result;
    }
    if(StartsWith(method,"replay1"))return Integrate(environment,row,gsl_odeiv2_step_rkf45,1.0e-11,{1.0e-16,1.0e-22,1.0e-22});
    if(StartsWith(method,"replay2"))return Integrate(environment,row,gsl_odeiv2_step_rkf45,1.0e-12,{1.0e-17,1.0e-23,1.0e-23});
    if(method=="oracle1")return Integrate(environment,row,gsl_odeiv2_step_rk8pd,1.0e-12,{1.0e-17,1.0e-23,1.0e-23});
    if(method=="oracle2")return Integrate(environment,row,gsl_odeiv2_step_rk8pd,1.0e-13,{1.0e-18,1.0e-24,1.0e-24});
    throw std::runtime_error("unknown reconstruction method");
}

std::shared_ptr<const BNV::FrozenControlledBnvRunContext> BuildContext(const Inputs& inputs,const std::string& method)
{
    const auto method_work=inputs.work/method;if(std::filesystem::exists(method_work))throw std::runtime_error("fresh method work root required");
    const auto& card=Campaign::RunCards().back();
    require(card.identity=="CPL-P2-LINEAR-QSS-v1"&&card.partition=="P2"&&card.fractional_drive_per_year==-1.0e-12&&
      card.duration_year==5.0e5&&card.checkpoints==8193&&!card.reaction_free,"immutable reconstruction P2 card changed");
    auto fixture=Fixture(inputs.profile_root,inputs.certificate.string(),method_work/"owning-star",80000);auto channels=Channels(fixture);
    EOS::CompOSE_Thermo::Options options;options.Tmin_for_derivative_MeV=0;options.clamp_to_domain=false;
    auto thermal=std::make_shared<const RC::FrozenThermalSource>(inputs.thermal,options,
      "controlled mathematical fixed-background free-gas entropy; qualified radial80000");
    auto tangent=Campaign::Tangent(fixture,inputs.profile_root,method_work);auto bnv_token=std::make_shared<BNV::BnvDependencyToken>();
    auto monitor=Phase6A1Test::LoadFrozenMonitor(inputs.frozen,tangent,bnv_token);auto run_token=std::make_shared<RC::RunDependencyToken>();
    auto spin=std::make_shared<const BNV::StaticZeroSpinHistory>(run_token);
    auto ordinary=Context(fixture,channels,thermal,spin,Campaign::Qualification(inputs.profile_root,inputs.certificate,inputs.entry),run_token);
    const double mu_B=Campaign::Column(inputs.coefficients,"mu_B_inf");
    const std::string potential="frozen Structure-1 B0 equilibrium mu_B plus governed moving-reference actual-potential correction";
    const double Bdot=card.fractional_drive_per_year*tangent->B0Count()/Year;require(Bdot==-2.4136520263641375e37,"immutable reconstruction Bdot changed");
    const std::string fate_id="phase6a1-P2-terminal-fate-v1";BNV::ProductFateLedger fate(fate_id,{{"controlled-terminal-product",BNV::TerminalProductFate::BoundInert,1,"abstract-neutron-disappearance"}});
    auto partition=Partition(card,fixture,bnv_token,fate_id);auto history=std::make_shared<const Phase6A1Test::UniformProperNeutronSinkHistory>(fixture.central,
      tangent->B0Count(),Bdot,tangent->DomainIdentity(),tangent->StarIdentity(),tangent->SequenceStateIdentity(),partition->Identity(),fate_id,bnv_token);
    return std::make_shared<const BNV::FrozenControlledBnvRunContext>(ordinary,tangent,history,partition,fate,monitor,channels,mu_B,potential,card.identity);
}

std::string ObservationName(std::size_t observation)
{
    std::ostringstream out;out<<"obs-"<<std::setw(3)<<std::setfill('0')<<observation;return out.str();
}

void WriteMeta(const std::filesystem::path& path,const std::string& method,const MatrixRow& row,const MethodResult& result,
    double total_wall,double total_user,double total_sys)
{
    std::ostringstream out;out<<std::setprecision(17)
      <<"method\tobservation_index\tt_obs_s\tleft_endpoint_index\tt_left_s\tright_endpoint_index\tt_right_s\tcategory\tdeep\texact\tstrict"
      <<"\taccepted\trejected\trhs\tsolve_wall_s\tsolve_cpu_user_s\tsolve_cpu_sys_s\tjob_wall_s\tjob_cpu_user_s\tjob_cpu_sys_s"
      <<"\tfL_x\tfL_eta_e\tfL_eta_mu\tfR_x\tfR_eta_e\tfR_eta_mu\tpid\n"
      <<method<<'\t'<<row.observation<<'\t'<<row.t_obs<<'\t'<<row.left_index<<'\t'<<row.t_left<<'\t'<<row.right_index<<'\t'<<row.t_right
      <<'\t'<<row.category<<'\t'<<row.deep<<'\t'<<row.exact<<'\t'<<row.strict<<'\t'<<result.solve.accepted<<'\t'<<result.solve.rejected<<'\t'
      <<result.solve.rhs<<'\t'<<result.solve.wall_s<<'\t'<<result.solve.user_s<<'\t'<<result.solve.sys_s<<'\t'<<total_wall<<'\t'<<total_user<<'\t'<<total_sys;
    if(result.has_endpoint_rhs)for(double value:result.f_left)out<<'\t'<<value;else out<<"\tNA\tNA\tNA";
    if(result.has_endpoint_rhs)for(double value:result.f_right)out<<'\t'<<value;else out<<"\tNA\tNA\tNA";
    out<<'\t'<<getpid()<<'\n';AtomicWrite(path,out.str());
}

constexpr std::size_t AuthorizedNewIntegrations=956;

std::vector<std::string> ReadLines(const std::filesystem::path& path)
{
    std::vector<std::string> lines;std::ifstream in(path);std::string line;while(std::getline(in,line))if(!line.empty())lines.push_back(line);return lines;
}

// Atomic execution accounting: a solve runs only if its ID is in the frozen authorization list and absent
// from the global ledger, and only while fewer than 956 new integrations exist. Retries would collide.
void AuthorizeSolve(const Inputs& inputs,const std::string& solve_id)
{
    const auto authorized=ReadLines(inputs.authorized_solves);
    if(authorized.size()!=AuthorizedNewIntegrations)throw std::runtime_error("authorized solve list is not the frozen 956-entry list");
    if(std::find(authorized.begin(),authorized.end(),solve_id)==authorized.end())throw std::runtime_error("solve not authorized: "+solve_id);
    std::size_t executed=0;
    for(const auto& line:ReadLines(inputs.execution_ledger))
    {
        if(StartsWith(line,"solve_id\t"))continue;++executed;
        if(line.substr(0,line.find('\t'))==solve_id)throw std::runtime_error("solve already executed in this recovery: "+solve_id);
    }
    if(executed>=AuthorizedNewIntegrations)throw std::runtime_error("956-integration authorization exhausted");
}

void RecordSolve(const Inputs& inputs,const std::string& line)
{
    const bool fresh=!std::filesystem::exists(inputs.execution_ledger);
    std::ofstream ledger(inputs.execution_ledger,std::ios::app);if(!ledger)throw std::runtime_error("cannot open execution ledger");
    if(fresh)ledger<<"solve_id\tmethod\tobservation_index\tpid\tstart_unix_ns\tend_unix_ns\tresult_sha256\tmeta_sha256\n";
    ledger<<line;ledger.flush();if(!ledger)throw std::runtime_error("execution ledger write failure");
}

void RunMethodBatch(const Inputs& inputs,const std::vector<MatrixRow>& positive,const std::string& method)
{
    const auto method_output=inputs.output/method;if(std::filesystem::exists(method_output))throw std::runtime_error("fresh method output root required");
    std::filesystem::create_directories(method_output);auto context=BuildContext(inputs,method);
    std::ofstream ledger(method_output/"solve_ledger.tsv");if(!ledger)throw std::runtime_error("cannot create method solve ledger");
    ledger<<"solve_id\tobservation_index\tmethod\tleft_endpoint_index\tt_left_s\ttarget_s\tpid\tstart_unix_ns\tend_unix_ns\texit_status\tresult_sha256\tmeta_sha256\n";
    std::vector<MatrixRow> rows;rows.reserve(positive.size()+1);rows.push_back(InitialRow());rows.insert(rows.end(),positive.begin(),positive.end());
    for(const auto& row:rows)
    {
        const auto stem=ObservationName(row.observation);
        const auto result_path=method_output/(stem+".tsv");
        const auto meta_path=method_output/(stem+".meta.tsv");
        if(std::filesystem::exists(result_path)||std::filesystem::exists(meta_path))throw std::runtime_error("hidden retry/output collision");
        const bool integration=row.strict&&IsIntegration(method);const std::string solve_id=method+"-"+stem;
        if(integration)AuthorizeSolve(inputs,solve_id);
        const auto unix_start=std::chrono::system_clock::now().time_since_epoch();const auto wall_start=std::chrono::steady_clock::now();rusage before{},after{};getrusage(RUSAGE_SELF,&before);
        EvaluationEnvironment environment(context);const auto result=EvaluateMethod(environment,row,method);environment.Save(row.t_obs,result.state,result_path.string()+".tmp");
        std::filesystem::rename(result_path.string()+".tmp",result_path);getrusage(RUSAGE_SELF,&after);
        const double wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-wall_start).count();
        WriteMeta(meta_path,method,row,result,wall,CpuSeconds(after.ru_utime)-CpuSeconds(before.ru_utime),CpuSeconds(after.ru_stime)-CpuSeconds(before.ru_stime));
        const auto unix_end=std::chrono::system_clock::now().time_since_epoch();
        if(integration)
        {
            ledger<<solve_id<<'\t'<<row.observation<<'\t'<<method<<'\t'<<row.left_index<<'\t'<<std::setprecision(17)<<row.t_left<<'\t'<<row.t_obs<<'\t'
              <<getpid()<<'\t'<<std::chrono::duration_cast<std::chrono::nanoseconds>(unix_start).count()<<'\t'
              <<std::chrono::duration_cast<std::chrono::nanoseconds>(unix_end).count()<<"\t0\t"<<Sha256(result_path)<<'\t'<<Sha256(meta_path)<<'\n';
            ledger.flush();if(!ledger)throw std::runtime_error("solve-ledger write failure");
            std::ostringstream line;line<<solve_id<<'\t'<<method<<'\t'<<row.observation<<'\t'<<getpid()<<'\t'
              <<std::chrono::duration_cast<std::chrono::nanoseconds>(unix_start).count()<<'\t'
              <<std::chrono::duration_cast<std::chrono::nanoseconds>(unix_end).count()<<'\t'<<Sha256(result_path)<<'\t'<<Sha256(meta_path)<<'\n';
            RecordSolve(inputs,line.str());
        }
    }
    AtomicWrite(method_output/"COMPLETE",method+"\n");
}

const std::vector<std::string>& CandidateStages()
{
    static const std::vector<std::string> stages{"linear","linear-repeat","hermite","hermite-repeat","replay1","replay2","replay1-repeat","replay2-repeat"};
    return stages;
}

// One stage = one candidate method in one fresh child process (concurrency 1). No retry on failure.
void RunStage(const Inputs& inputs,const std::vector<MatrixRow>& rows,const std::string& method)
{
    const auto& stages=CandidateStages();if(std::find(stages.begin(),stages.end(),method)==stages.end())throw std::runtime_error("stage must be a candidate method; oracle rerun is not authorized");
    if(std::filesystem::exists(inputs.output/method)||std::filesystem::exists(inputs.work/method))throw std::runtime_error("fresh stage output/work root required");
    std::filesystem::create_directories(inputs.output);std::filesystem::create_directories(inputs.work);
    const pid_t pid=fork();if(pid<0)throw std::runtime_error("fork failed");
    if(pid==0)
    {
        try{RunMethodBatch(inputs,rows,method);_exit(0);}catch(const std::exception& error){std::ofstream out(inputs.output/(method+".error"));out<<error.what()<<'\n';_exit(2);}
    }
    int status=0;if(waitpid(pid,&status,0)!=pid)throw std::runtime_error("wait failed");
    if(!WIFEXITED(status)||WEXITSTATUS(status)!=0)throw std::runtime_error("stage "+method+" failed; no retry authorized");
    if(!std::filesystem::exists(inputs.output/method/"COMPLETE"))throw std::runtime_error("stage completed without COMPLETE marker");
    std::cout<<"STAGE PASS "<<method<<'\n';
}

// Pre-integration dry run: build the corrected context, evaluate diagnostics/RHS at authenticated accepted
// endpoints for representative observations, and check RHS repeatability. No ODE solve is performed.
void DryRun(const Inputs& inputs,const std::vector<MatrixRow>& positive,const std::vector<std::size_t>& observations)
{
    const std::string method="dryrun";const auto output=inputs.output/method;
    if(std::filesystem::exists(output)||std::filesystem::exists(inputs.work/method))throw std::runtime_error("fresh dry-run root required");
    std::filesystem::create_directories(output);auto context=BuildContext(inputs,method);
    std::ofstream summary(output/"dryrun.tsv");summary<<std::setprecision(17)<<"observation_index\tcategory\tdeep\texact\tstrict\tt_left_s\tt_right_s\tfL_x\tfL_eta_e\tfL_eta_mu\tfL_repeat_identical\tfR_x\tfR_eta_e\tfR_eta_mu\tfR_repeat_identical\n";
    for(const std::size_t index:observations)
    {
        if(index<1||index>positive.size())throw std::runtime_error("dry-run observation outside 1..240");
        const auto& row=positive[index-1];const auto stem=ObservationName(row.observation);
        std::array<double,3> fl{},fr{};bool same_l=false,same_r=false;
        {EvaluationEnvironment environment(context);fl=environment.Derivative(row.t_left,row.left);}
        {EvaluationEnvironment environment(context);same_l=environment.Derivative(row.t_left,row.left)==fl;}
        {EvaluationEnvironment environment(context);fr=environment.Derivative(row.t_right,row.right);}
        {EvaluationEnvironment environment(context);same_r=environment.Derivative(row.t_right,row.right)==fr;}
        {EvaluationEnvironment environment(context);environment.Save(row.t_left,row.left,output/(stem+".left.tsv"));}
        {EvaluationEnvironment environment(context);environment.Save(row.t_right,row.right,output/(stem+".right.tsv"));}
        summary<<row.observation<<'\t'<<row.category<<'\t'<<row.deep<<'\t'<<row.exact<<'\t'<<row.strict<<'\t'<<row.t_left<<'\t'<<row.t_right;
        for(double v:fl)summary<<'\t'<<v;summary<<'\t'<<same_l;for(double v:fr)summary<<'\t'<<v;summary<<'\t'<<same_r<<'\n';
        if(!same_l||!same_r)throw std::runtime_error("endpoint RHS not repeatable in dry run");
    }
    if(!summary)throw std::runtime_error("dry-run summary write failure");AtomicWrite(output/"COMPLETE",method+"\n");
}

Inputs ParseInputs(int argc,char** argv,int first)
{
    std::map<std::string,std::filesystem::path> values;
    for(int i=first;i+1<argc;i+=2)
    {
        const std::string key=argv[i];if(!StartsWith(key,"--"))throw std::runtime_error("expected --name value pairs");
        if(!values.emplace(key.substr(2),argv[i+1]).second)throw std::runtime_error("duplicate argument "+key);
    }
    if((argc-first)%2!=0)throw std::runtime_error("dangling argument");
    const std::vector<std::string> required{"matrix","profile-root","certificate","thermal","work-root","entry-manifest","frozen-certificate",
      "coefficients","output-root","predeclaration","pretrajectory","oracle-qualified-flag","authorized-solves","execution-ledger"};
    for(const auto& key:required)if(!values.count(key))throw std::runtime_error("missing required argument --"+key);
    for(const auto& [key,value]:values)if(std::find(required.begin(),required.end(),key)==required.end())throw std::runtime_error("unknown argument --"+key);
    Inputs inputs{values["matrix"],values["profile-root"],values["certificate"],values["thermal"],values["work-root"],values["entry-manifest"],
      values["frozen-certificate"],values["coefficients"],values["output-root"],values["predeclaration"],values["pretrajectory"],
      values["oracle-qualified-flag"],values["authorized-solves"],values["execution-ledger"]};
    ValidateProfileRoot(inputs.profile_root);
    return inputs;
}
}

int main(int argc,char** argv)
{
 try
 {
    gsl_set_error_handler_off();std::cout<<std::setprecision(17)<<std::unitbuf;
    if(argc==3&&std::string(argv[1])=="matrix-check")
    {
        const auto rows=ReadMatrix(argv[2]);std::cout<<"SOLVE_MATRIX PASS rows "<<rows.size()<<" strict "<<ExpectedStrictInterior<<" total "<<ExpectedTotalSolves<<'\n';return 0;
    }
    const std::string mode=argc>1?argv[1]:"";
    if(mode!="dry-run"&&mode!="stage")throw std::runtime_error("mode must be matrix-check, dry-run or stage <method>");
    const int first=mode=="stage"?3:2;if(argc<first)throw std::runtime_error("stage requires a method name");
    const auto inputs=ParseInputs(argc,argv,first);const auto rows=ReadMatrix(inputs.matrix);
    {std::ifstream in(inputs.predeclaration);std::string text((std::istreambuf_iterator<char>(in)),{});require(
      text.find("32f3277cdf3984318fe2323da825de1d5e337778bb05106725fb0fcfa0529616")!=std::string::npos&&
      text.find("TOTAL AUTHORIZED LOCAL INTEGRATIONS")!=std::string::npos,"committed reconstruction predeclaration missing");}
    {std::ifstream in(inputs.pretrajectory);std::string text((std::istreambuf_iterator<char>(in)),{});require(text.find("PRETRAJECTORY PASS")!=std::string::npos&&text.find("no BNV trajectory had been generated")!=std::string::npos,"committed historical pretrajectory evidence missing");}
    {std::ifstream in(inputs.oracle_flag);std::string text((std::istreambuf_iterator<char>(in)),{});require(text.find("ORACLE QUALIFIED")!=std::string::npos,"independently verified oracle qualification flag missing");}
    if(mode=="dry-run")
    {
        std::vector<std::size_t> observations;
        // Representative set: a no-knot interior observation, each one-knot observation, the deepest-interior observation, the exact endpoint.
        std::size_t deepest=0;double best=0;for(const auto& row:rows){if(!row.strict)continue;const double s=(row.t_obs-row.t_left)/(row.t_right-row.t_left);const double score=std::min(s,1-s)*(row.t_right-row.t_left);if(score>best){best=score;deepest=row.observation;}}
        for(const auto& row:rows)if(row.category=='A'&&row.strict&&!row.deep){observations.push_back(row.observation);break;}
        for(const auto& row:rows)if(row.category=='B')observations.push_back(row.observation);
        observations.push_back(deepest);for(const auto& row:rows)if(row.exact)observations.push_back(row.observation);
        DryRun(inputs,rows,observations);std::cout<<"DRY_RUN PASS observations";for(auto o:observations)std::cout<<' '<<o;std::cout<<'\n';return 0;
    }
    RunStage(inputs,rows,argv[2]);return 0;
 }
 catch(const std::exception& error){std::cerr<<"STOP "<<error.what()<<'\n';return 1;}
}
