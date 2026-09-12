#pragma once
#include <CompactStar/Physics/Rotochemical/FrozenRotochemicalRunContext.hpp>
#include <CompactStar/Physics/Evolution/EvolutionSystem.hpp>
#include <CompactStar/Physics/Evolution/StateLayout.hpp>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <algorithm>
#include <exception>
#include <limits>
namespace CompactStar::Physics::Rotochemical
{
struct ComponentTolerances
{
    double relative=1e-7;
    std::array<double,3> absolute{1e-12,1e-18,1e-18};
    static ComponentTolerances Refined(){return {1e-9,{1e-14,1e-20,1e-20}};}
    void Validate()const{if(!(relative>0)||!std::isfinite(relative))throw std::runtime_error("invalid relative tolerance");for(double a:absolute)if(!(a>0)||!std::isfinite(a))throw std::runtime_error("invalid component tolerance");}
    gsl_odeiv2_control* Allocate()const{Validate();auto* p=gsl_odeiv2_control_scaled_new(1,relative,1,0,absolute.data(),3);if(!p)throw std::bad_alloc();return p;}
};
struct StepOutput {double time_s;size_t accepted,rejected;double last_step_s;};
struct IntegrationStatistics
{
    size_t accepted_steps=0,rejected_steps=0,rhs_evaluations=0;
    double minimum_step_s=std::numeric_limits<double>::infinity(),maximum_step_s=0;
    std::vector<StepOutput> outputs;
};
// Exact D_i = atol_i + rtol |y_i|: eps_abs=1,a_y=1,a_dydt=0.
// This bounded adapter forwards to public EvolutionSystem without modifying it.
class ScaledRKF45 final
{
  public:
    ScaledRKF45(Evolution::EvolutionSystem& system,const Evolution::StateLayout& layout,
      std::shared_ptr<const FrozenRotochemicalRunContext> context,ComponentTolerances tolerances={})
      :system_(system),context_(std::move(context)),tolerances_(tolerances)
    {tolerances_.Validate();if(!context_||layout.TotalSize()!=3||layout.BlockSize(State::StateTag::Thermal)!=1||layout.BlockSize(State::StateTag::Chem)!=2||layout.Offset(State::StateTag::Thermal)!=0||layout.Offset(State::StateTag::Chem)!=1)throw std::runtime_error("scaled RKF45 requires thermal,Npe,NpMu layout");}
    void Derivative(double t,const double* y,double* out)
    {
        context_->RequireCheapCurrent();if(stats_)++stats_->rhs_evaluations;std::array<double,3> local;const int rc=system_(t,y,local.data());if(rc)throw std::runtime_error("EvolutionSystem RHS failure");context_->RequireCheapCurrent();for(double v:local)if(!std::isfinite(v))throw std::runtime_error("nonfinite derivative");std::copy(local.begin(),local.end(),out);
    }
    void Integrate(double start,const std::vector<double>& checkpoints,double* y,IntegrationStatistics& stats,bool enforce_thermal_domain=true)
    {
        context_->RequireFullCurrent();if(!std::isfinite(start))throw std::runtime_error("invalid start time");if(checkpoints.empty()||checkpoints.size()>1000)throw std::runtime_error("invalid checkpoint count");double previous=start;for(double t:checkpoints){if(!std::isfinite(t)||!(t>previous))throw std::runtime_error("unordered checkpoints");previous=t;}
        stats={};stats_=&stats;failure_=nullptr;
        std::unique_ptr<gsl_odeiv2_step,decltype(&gsl_odeiv2_step_free)> step(gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkf45,3),gsl_odeiv2_step_free);
        std::unique_ptr<gsl_odeiv2_control,decltype(&gsl_odeiv2_control_free)> control(tolerances_.Allocate(),gsl_odeiv2_control_free);
        std::unique_ptr<gsl_odeiv2_evolve,decltype(&gsl_odeiv2_evolve_free)> evolve(gsl_odeiv2_evolve_alloc(3),gsl_odeiv2_evolve_free);
        if(!step||!evolve)throw std::bad_alloc();gsl_odeiv2_system sys{Callback,nullptr,3,this};double t=start,h=std::min(1.,(checkpoints.front()-start)*1e-3),last=0;
        try{ValidateAccepted(y,enforce_thermal_domain);system_.NotifyStart(start,checkpoints.back(),y);for(size_t i=0;i<checkpoints.size();++i){size_t count=0;while(t<checkpoints[i]){if(++count>100000)throw std::runtime_error("RKF45 internal step budget exceeded");context_->RequireCheapCurrent();double before=t;int rc=gsl_odeiv2_evolve_apply(evolve.get(),control.get(),step.get(),&sys,&t,checkpoints[i],&h,y);stats.rejected_steps=evolve->failed_steps;if(failure_)std::rethrow_exception(failure_);if(rc!=GSL_SUCCESS)throw std::runtime_error(std::string("RKF45 failure: ")+gsl_strerror(rc));last=t-before;if(!(last>0)||!std::isfinite(last))throw std::runtime_error("RKF45 failed progress");++stats.accepted_steps;stats.minimum_step_s=std::min(stats.minimum_step_s,last);stats.maximum_step_s=std::max(stats.maximum_step_s,last);ValidateAccepted(y,enforce_thermal_domain);std::array<double,3> endpoint;Derivative(t,y,endpoint.data());}
          stats.outputs.push_back({t,stats.accepted_steps,stats.rejected_steps,last});system_.NotifySample(t,y,i);
        }system_.NotifyFinish(t,y,true);}catch(...){auto original=std::current_exception();stats_=nullptr;try{system_.NotifyFinish(t,y,false);}catch(...){}std::rethrow_exception(original);}stats_=nullptr;
    }
  private:
    static int Callback(double t,const double*y,double*out,void*p)noexcept{auto& self=*static_cast<ScaledRKF45*>(p);try{self.Derivative(t,y,out);return GSL_SUCCESS;}catch(...){self.failure_=std::current_exception();return GSL_EBADFUNC;}}
    void ValidateAccepted(const double*y,bool thermal)const{context_->RequireCheapCurrent();for(size_t i=0;i<3;++i)if(!std::isfinite(y[i]))throw std::runtime_error("nonfinite accepted state");double T=1e8*std::exp(y[0]);if(!(T>0)||!std::isfinite(T))throw std::runtime_error("invalid accepted temperature");if(thermal&&!(T*BoltzmannMeVPerK>1e-5&&T*BoltzmannMeVPerK<1))throw std::runtime_error("accepted state outside predeclared thermal cache domain");}
    Evolution::EvolutionSystem& system_;
    const std::shared_ptr<const FrozenRotochemicalRunContext> context_;
    const ComponentTolerances tolerances_;
    IntegrationStatistics* stats_=nullptr;
    std::exception_ptr failure_;
};
}
