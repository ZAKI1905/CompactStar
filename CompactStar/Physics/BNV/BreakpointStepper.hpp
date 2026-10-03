#pragma once

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>
#include <set>

namespace CompactStar::Physics::BNV
{
// Candidate ADR-0017 A1: a small three-state, scalar-interface stepper.
// The RHS callback evaluates the requested cell's smooth extension. No schedule
// enters this object. Event locations belong to this instance and this trajectory.
class BreakpointStepper
{
public:
    using State=std::array<double,3>;
    using RHS=std::function<void(double,const State&,State&,int)>;
    struct Event
    {
        int knot=-1,direction=0;
        double time=0,time_low=0,time_high=0,residual=0,time_ulp=0;
        State state{},rk_error{};
    };
    struct Statistics
    {
        std::size_t rejected=0,event_rejected=0,applies=0,event_trials=0,discarded=0,
                    resets=0,rhs=0,stage_excursions=0,initializations=1;
    };
    struct Step
    {
        double left=0,right=0,h_next=0;
        State before{},after{},rk_error{};
        int cell=-1;
        Event event;
    };
    BreakpointStepper(RHS rhs,std::vector<double> knots,const gsl_odeiv2_step_type* method,
                      double relative,State absolute,double max_step=INFINITY)
      :rhs_(std::move(rhs)),knots_(std::move(knots)),max_step_(max_step),step_(gsl_odeiv2_step_alloc(method,3),gsl_odeiv2_step_free),
       control_(gsl_odeiv2_control_scaled_new(1,relative,1,0,absolute.data(),3),gsl_odeiv2_control_free),
       evolve_(gsl_odeiv2_evolve_alloc(3),gsl_odeiv2_evolve_free)
    {
        Require(step_&&control_&&evolve_,"breakpoint GSL allocation failure");
        Require(relative>0&&max_step>0,"invalid breakpoint configuration");
        for(double a:absolute)Require(a>0&&std::isfinite(a),"invalid absolute tolerance");
        for(std::size_t i=0;i<knots_.size();++i)
            Require(std::isfinite(knots_[i])&&(i==0||knots_[i]>knots_[i-1]),"unordered breakpoint grid");
        Require(knots_.empty()||knots_.size()>=2,"incomplete breakpoint domain");
        system_={Call,nullptr,3,this};
    }
    const Statistics& Stats() const {return stats_;}
    const std::vector<Event>& Events() const {return events_;}
    std::size_t TargetCount() const {return targets_.size();}
    Step Advance(double& t,State& y,double& h,double target)
    {
        Require(target>t&&h>0&&std::isfinite(t)&&std::isfinite(target),"invalid advance interval");
        const double before=t;const State initial=y;const double proposal=std::max(std::min(h,max_step_),std::nextafter(t,INFINITY)-t);
        SelectCell(t,y);
        Step record;record.left=t;record.before=y;record.cell=cell_;
        h=proposal;
        ResetStageAudit();
        Apply(t,y,h,std::min(target,std::max(t+max_step_,std::nextafter(t,INFINITY))));
        const double trial_time=t;
        const int direction=y[0]>initial[0]?1:y[0]<initial[0]?-1:0;
        int knot=-1;
        if(cell_>=0)
        {
            if(y[0]>=knots_[cell_+1]&&y[0]>initial[0])knot=cell_+1;
            else if(y[0]<=knots_[cell_]&&y[0]<initial[0])knot=cell_;
            // A reversing trial that leaves the cell cannot certify event order.
            Require(!(min_derivative_<0&&max_derivative_>0&&stage_outside_),
                    "unresolved reversing breakpoint trial");
        }
        if(knot>=0)
        {
            Require(knot>0&&knot+1<int(knots_.size()),"breakpoint domain exit");
            Require(direction!=0,"unresolved crossing direction");
            ++stats_.discarded;
            // A discarded evolve call has cached the WRONG endpoint derivative.
            t=before;y=initial;h=proposal;Reset();
            const double boundary=knots_[knot];
            Require(direction*(boundary-initial[0])>0,"ambiguous boundary start");
            double low=before,high=trial_time,candidate=high;
            State root=y,error{};double root_t=before;
            // Every trial starts from the same saved physical state and fresh GSL
            // caches, with the same method/tolerances and its OWN event solve.
            for(int iteration=0;iteration<80;++iteration)
            {
                State dy{};rhs_(root_t,root,dy,cell_);
                double next=root_t+(boundary-root[0])/dy[0];
                if(!std::isfinite(next)||!(next>low&&next<high))next=low+(high-low)/2;
                if(!(next>low&&next<high))break;
                candidate=next;double tt=before,hh=std::min(proposal,candidate-before);State yy=initial;
                Reset();ResetStageAudit();
                std::size_t work=0;
                while(tt<candidate)
                {
                    Require(++work<10000,"event IVP work exhausted");
                    const auto rejected_before=stats_.rejected;
                    Apply(tt,yy,hh,candidate);++stats_.event_trials;
                    stats_.event_rejected+=stats_.rejected-rejected_before;
                }
                Require(!(min_derivative_<0&&max_derivative_>0),"non-monotone event approach");
                root=yy;root_t=tt;
                for(int i=0;i<3;++i)error[i]=std::abs(evolve_->yerr[i]);
                const double residual=direction*(root[0]-boundary);
                if(residual<0)low=tt;else high=tt;
                if(residual==0) {low=high=tt;break;}
                const double ulp=std::nextafter(tt,INFINITY)-tt;
                if(high-low<=2*ulp)break;
            }
            Require(root_t>before,"event failed finite-time progress");
            const double residual=root[0]-boundary;
            Require(std::abs(residual)<=64*std::numeric_limits<double>::epsilon()*std::max(1.,std::abs(boundary)),
                    "event localization residual unresolved");
            // Retain all three integrated components. There is NO knot projection.
            t=root_t;y=root;h=proposal;
            record.event={knot,direction,t,low,high,residual,std::nextafter(t,INFINITY)-t,y,error};
            events_.push_back(record.event);
            last_knot_=knot;
            cell_=direction>0?knot:knot-1;
            Reset();
        }
        record.right=t;record.after=y;record.h_next=h;
        if(knot<0)for(int i=0;i<3;++i)record.rk_error[i]=std::abs(evolve_->yerr[i]);
        else record.rk_error=record.event.rk_error;
        Require(t>before,"breakpoint step failed progress");
        return record;
    }
private:
    static void Require(bool ok,const char* msg){if(!ok)throw std::runtime_error(msg);}
    static int Call(double t,const double* y,double* out,void* p) noexcept
    {
        auto& self=*static_cast<BreakpointStepper*>(p);
        try
        {
            State in{{y[0],y[1],y[2]}},d{};self.rhs_(t,in,d,self.cell_);++self.stats_.rhs;
            self.min_derivative_=std::min(self.min_derivative_,d[0]);
            self.max_derivative_=std::max(self.max_derivative_,d[0]);
            if(self.cell_>=0&&(in[0]<self.knots_[self.cell_]||in[0]>self.knots_[self.cell_+1]))
            {self.stage_outside_=true;++self.stats_.stage_excursions;}
            for(int i=0;i<3;++i){Require(std::isfinite(d[i]),"nonfinite breakpoint RHS");out[i]=d[i];}
            return GSL_SUCCESS;
        }
        catch(...){self.failure_=std::current_exception();return GSL_EBADFUNC;}
    }
    void SelectCell(double t,const State& y)
    {
        if(knots_.empty()){cell_=-1;return;}
        Require(y[0]>=knots_.front()&&y[0]<=knots_.back(),"state outside breakpoint domain");
        State d{};rhs_(t,y,d,-1);
        int selected=std::min(int(knots_.size())-2,int(std::upper_bound(knots_.begin(),knots_.end(),y[0])-knots_.begin())-1);
        // A just-landed residual is not snapped. The selected continuation owns
        // that tiny residual strip; record.event accounts for its width.
        if(last_knot_>=0&&std::abs(y[0]-knots_[last_knot_])<=64*std::numeric_limits<double>::epsilon()*std::max(1.,std::abs(y[0])))
        {
            Require(d[0]!=0,"unresolved tangential knot start");
            selected=d[0]>0?last_knot_:last_knot_-1;
        }
        else if(selected>0&&y[0]==knots_[selected])
        {
            Require(d[0]!=0,"unresolved stationary knot start");
            if(d[0]<0)--selected;
        }
        if(cell_!=selected){cell_=selected;Reset();}
    }
    void Reset(){gsl_odeiv2_evolve_reset(evolve_.get());gsl_odeiv2_step_reset(step_.get());++stats_.resets;}
    void ResetStageAudit(){min_derivative_=INFINITY;max_derivative_=-INFINITY;stage_outside_=false;}
    void Apply(double& t,State& y,double& h,double target)
    {
        targets_.insert(target);
        const auto rejected=evolve_->failed_steps;const double old=t;++stats_.applies;
        const int rc=gsl_odeiv2_evolve_apply(evolve_.get(),control_.get(),step_.get(),&system_,&t,target,&h,y.data());
        stats_.rejected+=evolve_->failed_steps-rejected;
        if(failure_)std::rethrow_exception(failure_);
        Require(rc==GSL_SUCCESS&&t>old,"breakpoint GSL failure or no progress");
        for(double v:y)Require(std::isfinite(v),"nonfinite breakpoint state");
    }
    RHS rhs_;std::vector<double> knots_;double max_step_;
    std::unique_ptr<gsl_odeiv2_step,decltype(&gsl_odeiv2_step_free)> step_;
    std::unique_ptr<gsl_odeiv2_control,decltype(&gsl_odeiv2_control_free)> control_;
    std::unique_ptr<gsl_odeiv2_evolve,decltype(&gsl_odeiv2_evolve_free)> evolve_;
    gsl_odeiv2_system system_{};std::exception_ptr failure_;
    int cell_=-1,last_knot_=-1;
    double min_derivative_=INFINITY,max_derivative_=-INFINITY;bool stage_outside_=false;
    Statistics stats_;std::vector<Event> events_;std::set<double> targets_;
};
} // namespace CompactStar::Physics::BNV
