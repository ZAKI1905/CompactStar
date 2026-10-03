#include <CompactStar/Physics/BNV/BreakpointStepper.hpp>
#include <iostream>
#include "checked_evidence.hpp"
#include <cstring>
using B=CompactStar::Physics::BNV::BreakpointStepper;
void need(bool x,const char* m){if(!x)throw std::runtime_error(m);}
B::State solve(int direction,double maxstep,const gsl_odeiv2_step_type* method,double start=-.2)
{
    B solver([=](double,const B::State& y,B::State& d,int cell){
      const bool right=cell<0?y[0]>=0:cell==1;
      d[0]=direction*(1+(right?2.:.5)*y[0]);d[1]=y[0];d[2]=y[1];},
      {-2,0,2},method,1e-12,{{1e-14,1e-14,1e-14}},maxstep);
    double t=0,h=.7;B::State y{{direction>0?start:-start,0,0}};
    while(t<.5)solver.Advance(t,y,h,.5);
    need(solver.Events().size()==(start==0?0:1),"event count");
    if(start!=0){const auto& e=solver.Events()[0];
      double exact=direction>0?-std::log(1+.5*start)/.5:std::log(1-2*start)/2;
      need(std::abs(e.time-exact)<1e-11,"analytic event location");
      need(std::abs(e.residual)<1e-14,"event residual");}
    return y;
}
void rollback()
{
    gsl_odeiv2_system sys{[](double,const double*y,double*d,void*) -> int {d[0]=.1*y[0]*y[0]+y[0];return GSL_SUCCESS;},nullptr,1,nullptr};
    auto* s=gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkf45,1);
    auto* c=gsl_odeiv2_control_standard_new(1e-12,1e-12,1,0);
    auto* e=gsl_odeiv2_evolve_alloc(1);auto* fresh=gsl_odeiv2_evolve_alloc(1);
    double t=0,y=1,h=.1;gsl_odeiv2_evolve_apply(e,c,s,&sys,&t,10,&h,&y);
    const double t0=t,y0=y,h0=h;gsl_odeiv2_evolve_apply(e,c,s,&sys,&t,10,&h,&y);
    double a=t0,b=y0,hb=h0;gsl_odeiv2_evolve_apply(e,c,s,&sys,&a,10,&hb,&b);
    double u=t0,v=y0,hv=h0;gsl_odeiv2_evolve_apply(fresh,c,s,&sys,&u,10,&hv,&v);
    need(b!=v,"stale derivative falsifier failed to discriminate");
    gsl_odeiv2_evolve_reset(e);gsl_odeiv2_step_reset(s);
    a=t0;b=y0;hb=h0;gsl_odeiv2_evolve_apply(e,c,s,&sys,&a,10,&hb,&b);
    need(b==v&&a==u,"rollback reset regression");
    gsl_odeiv2_evolve_reset(e);a=t0;b=y0;hb=h0;
    gsl_odeiv2_evolve_apply(e,c,s,&sys,&a,t0+.3*h0,&hb,&b);need(hb==h0,"clipping lost h proposal");
    gsl_odeiv2_evolve_free(e);gsl_odeiv2_evolve_free(fresh);gsl_odeiv2_control_free(c);gsl_odeiv2_step_free(s);
}
int main(){try{
    gsl_set_error_handler_off();rollback();
    const auto out=std::filesystem::temp_directory_path()/"compactstar-controlled-write-failure.tsv";
    need(!std::filesystem::exists(out),"test output collision");bool failed=false;
    try{PublishTsv(out,"t\tx\n0\t1\n1\t2\n",2,2,0,7);}catch(const std::exception&){failed=true;}
    need(failed&&!std::filesystem::exists(out)&&!std::filesystem::exists(out.string()+".partial"),"write failure published success");
    PublishTsv(out,"t\tx\n0\t1\n1\t2\n",2,2,0);std::filesystem::remove(out);

    for(int d:{-1,1}){auto a=solve(d,INFINITY,gsl_odeiv2_step_rkf45),b=solve(d,.01,gsl_odeiv2_step_rk8pd);
      for(int i=0;i<3;++i)need(std::abs(a[i]-b[i])<1e-11,"coupled independent agreement");
      solve(d,.01,gsl_odeiv2_step_rk8pd,0);}
    // Reversal away from a knot must allow a later crossing in the other direction.
    B turning([](double t,const B::State& y,B::State& d,int){d={{std::cos(t),y[0],y[1]}};},
      {-2,0,2},gsl_odeiv2_step_rk8pd,1e-12,{{1e-14,1e-14,1e-14}},.05);
    double tt=0,hh=.05;B::State yy{{-.2,0,0}};
    while(tt<4)turning.Advance(tt,yy,hh,4);
    need(turning.Events().size()==2&&turning.Events()[0].direction==1&&turning.Events()[1].direction==-1,"direction change crossings");
    need(std::abs(yy[0]-(std::sin(4.)-.2))<1e-12,"turning endpoint");
    bool unresolved=false;
    try {B missed([](double t,const B::State&,B::State& d,int){d={{7*std::cos(7*t),0,0}};},
      {-2,0,2},gsl_odeiv2_step_rkf45,10,{{10,10,10}});
      double t=0,h=1;B::State y{{-.1,0,0}};missed.Advance(t,y,h,1);
    }catch(const std::exception&){unresolved=true;}
    need(unresolved,"missed crossing falsifier not rejected");
    std::cout<<"PASS analytic heating/cooling, knot starts, coupled method comparison, GSL rollback and h carry\n";
    return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
