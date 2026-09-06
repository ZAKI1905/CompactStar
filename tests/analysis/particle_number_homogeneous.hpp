// PB7 TEST-ONLY independent homogeneous sensitivity oracle.
// Different EOS (phase-space quadrature), variable (enthalpy), integrator (RKF45),
// centre expansion and derivative owner. No production number/sequence/Hartle call.
#pragma once
#include "../eos/structure1/tov_oracle.hpp"
#include <array>
namespace phase5b_oracle
{
// PB1 count-only oracle: direct phase-space EOS and independently integrated
// enthalpy TOV plus proper-volume counts. No sensitivity equations are run here.
struct Counts
{
    structure1::EnthalpyOracle background;
    static int rhs(double t,const double y[],double f[],void *ptr)
    {
        auto &o=*static_cast<Counts*>(ptr);
        auto q=o.background.at_H(o.background.Hc-t);
        double e=q.e*o.background.to_geo,p=q.p*o.background.to_geo;
        double u=y[0],v=y[1],metric=1-2*u*v;
        double A=metric/(v+4*M_PI*p);
        f[0]=2*A; f[1]=(4*M_PI*e-3*v)*A/u;
        for(int i=0;i<4;++i) f[2+i]=4*M_PI*std::sqrt(u)/std::sqrt(metric)*q.ns[i]*A;
        for(int i=0;i<6;++i) if(!std::isfinite(f[i])) return GSL_EBADFUNC;
        return GSL_SUCCESS;
    }
    std::array<double,4> solve(double rho,double tolerance,double delta)
    {
        background.Hc=background.H_for_rho(rho); auto q=background.at_H(background.Hc);
        double e=q.e*background.to_geo,p=q.p*background.to_geo;
        double u=3*delta/(2*M_PI*(e+3*p));
        double y[6]={u,4*M_PI*e/3-4*M_PI*(e+p)*q.slope*delta/5};
        for(int i=0;i<4;++i) y[2+i]=4*M_PI/3*q.ns[i]*std::pow(u,1.5);
        gsl_odeiv2_system sys{rhs,nullptr,6,this};
        auto driver=gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_rkf45,delta,tolerance*1e-3,tolerance);
        double t=delta; int status=gsl_odeiv2_driver_apply(driver,&t,background.Hc,y);
        gsl_odeiv2_driver_free(driver);
        if(status!=GSL_SUCCESS) throw std::runtime_error("PB1 independent direct-vacuum count integration failed");
        std::array<double,4> result{}; for(int i=0;i<4;++i) result[i]=y[2+i]*1e54; return result;
    }
};
struct Homogeneous
{
    structure1::EnthalpyOracle background;
    double hc=0;
    std::array<double,4> number_H(const structure1::OraclePoint &q) const
    {
        std::array<double,4> s{},out{};
        for(int i=0;i<4;++i)
        {
            double k=background.eos.b*std::cbrt(3*M_PI*M_PI*q.ns[i]);
            s[i]=std::hypot(background.eos.masses[i],k)*k/(M_PI*M_PI*std::pow(background.eos.b,3));
        }
        if(q.n==0) return out;
        if(q.ns[0]==0) out[1]=out[2]=q.h*s[1]*s[2]/(s[1]+s[2]);
        else
        { double dt_dh=s[1]/(s[1]+s[2]+s[3]); out={q.h*s[0],q.h*(s[2]+s[3])*dt_dh,q.h*s[2]*dt_dh,q.h*s[3]*dt_dh}; }
        return out;
    }
    static int rhs(double t,const double y[],double f[],void *ptr)
    {
        auto &o=*static_cast<Homogeneous*>(ptr);
        auto q=o.background.at_H(o.hc-t); auto nh=o.number_H(q);
        double e=q.e*o.background.to_geo,p=q.p*o.background.to_geo;
        double ph=e+p,eh=q.n>0?ph*q.slope:0;
        double u=y[0],v=y[1],U=y[2],V=y[3],metric=1-2*u*v,den=v+4*M_PI*p;
        if(!(u>0 && metric>0 && den>0)) return GSL_EBADFUNC;
        double A=metric/den,Au=-2*v/den,Av=(-2*u*den-metric)/(den*den),Ah=-metric*4*M_PI*ph/(den*den);
        double C=4*M_PI*e-3*v;
        f[0]=2*A; f[1]=C*A/u;
        f[2]=2*(Au*U+Av*V+Ah);
        f[3]=C*(Au/u-A/(u*u))*U+(-3*A/u+C*Av/u)*V+(4*M_PI*eh*A+C*Ah)/u;
        double weight=4*M_PI*std::sqrt(u)/std::sqrt(metric);
        for(int i=0;i<4;++i)
        {
            double F=weight*q.ns[i]*A;
            f[4+i]=F;
            f[8+i]=F*((.5/u+Au/A+v/metric)*U+(Av/A+u/metric)*V)+weight*(nh[i]*A+q.ns[i]*Ah);
        }
        for(int i=0;i<12;++i) if(!std::isfinite(f[i])) return GSL_EBADFUNC;
        return GSL_SUCCESS;
    }
    struct Result { std::array<double,4> N,B; double R,M,Hc; };
    Result solve(double rho,double tolerance,double delta)
    {
        hc=background.H_for_rho(rho); auto q=background.at_H(hc); auto nh=number_H(q);
        double e=q.e*background.to_geo,p=q.p*background.to_geo,ph=e+p,eh=ph*q.slope;
        double u=3*delta/(2*M_PI*(e+3*p)),U=-u*(eh+3*ph)/(e+3*p);
        double y[12]={u,4*M_PI*e/3-4*M_PI*ph*q.slope*delta/5,U,4*M_PI*eh/3};
        for(int i=0;i<4;++i)
        { y[4+i]=4*M_PI/3*q.ns[i]*std::pow(u,1.5); y[8+i]=4*M_PI/3*(nh[i]*std::pow(u,1.5)+1.5*q.ns[i]*std::sqrt(u)*U); }
        gsl_odeiv2_system system{rhs,nullptr,12,this};
        auto driver=gsl_odeiv2_driver_alloc_y_new(&system,gsl_odeiv2_step_rkf45,delta,tolerance*1e-3,tolerance);
        double t=delta; int status=gsl_odeiv2_driver_apply(driver,&t,hc,y); gsl_odeiv2_driver_free(driver);
        if(status!=GSL_SUCCESS) throw std::runtime_error("PB7 independent homogeneous integration failed");
        Result result{{},{},std::sqrt(y[0]),y[1]*std::pow(y[0],1.5)/background.solar_length,hc};
        // Upper integration endpoint t=Hc moves; F_i at exact vacuum is zero.
        // d/d epsilon_c = (d/d Hc)/(d epsilon_c/d Hc), using this oracle's EOS.
        for(int i=0;i<4;++i) { result.N[i]=y[4+i]*1e54; result.B[i]=y[8+i]/eh*1e54; }
        return result;
    }
};
// PB13 independent direct-vacuum continuation in decreasing enthalpy.
struct TailContinuation
{
    Homogeneous eos_owner;
    double Hs=0;
    static int rhs(double t,const double y[],double f[],void *ptr)
    {
        auto &o=*static_cast<TailContinuation*>(ptr); auto &background=o.eos_owner.background;
        auto q=background.at_H(o.Hs-t); auto nh=o.eos_owner.number_H(q);
        double e=q.e*background.to_geo,p=q.p*background.to_geo;
        double r=y[0],m=y[1],nu=y[2],s=y[3],sp=y[4],mh=y[5],ph=y[6];
        double D=1-2*m/r,g=(m+4*M_PI*r*r*r*p)/(r*(r-2*m)),E=std::exp(-2*nu);
        double lp=4*M_PI*r*(e+p)/D,eh=q.h*(nh[0]+nh[1])*background.to_geo;
        f[0]=1; f[1]=4*M_PI*r*r*e; f[2]=g; f[3]=sp;
        f[4]=-(4/r-lp)*sp+4*lp/r*s;
        f[5]=4*M_PI*r*r*eh*ph+std::pow(r,4)*E*D*sp*sp/12+8*M_PI/3*std::pow(r,4)*(e+p)*E*s*s;
        f[6]=-mh*(1+8*M_PI*r*r*p)/std::pow(r-2*m,2)-4*M_PI*(e+p)*r*r*ph/(r-2*m)
            +r*r*r*E*sp*sp/12+2*r*E*s*(s+r*sp-r*g*s)/3;
        double W=4*M_PI*r*r/std::sqrt(D);
        for(int i=0;i<4;++i)
        {f[7+i]=W*q.ns[i]; f[11+i]=W*(nh[i]*ph+q.ns[i]*(mh/(r-2*m)+r*r*E*s*s/3));}
        for(int i=0;i<15;++i) {f[i]/=g; if(!std::isfinite(f[i])) return GSL_EBADFUNC;}
        return GSL_SUCCESS;
    }
    struct Result {double R; std::array<double,4> N,A;};
    Result solve(const CompactStar::Core::NStar &star,double tolerance)
    {
        const auto &p=star.Profile(); const auto &f=star.RotationResponse(); const auto &h=*star.MonopoleResponse();
        Hs=eos_owner.background.H_for_pressure(CompactStar::RelativityUnits::PressureKmMinus2ToDynCm2((*p.GetPressure())[-1]));
        double y[15]={(*p.GetRadius())[-1],(*p.GetMass())[-1],(*p.GetMetricNu())[-1],f.omega_bar_over_Omega[-1],
            f.domega_bar_over_Omega_dr[-1],h.m0_over_Omega2[-1],h.p0star_over_Omega2[-1]};
        gsl_odeiv2_system sys{rhs,nullptr,15,this};
        auto driver=gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_rkf45,Hs/100,tolerance*1e-12,tolerance);
        gsl_odeiv2_driver_set_hmax(driver,Hs/100);
        double t=0; int status=gsl_odeiv2_driver_apply(driver,&t,Hs,y); gsl_odeiv2_driver_free(driver);
        if(status!=GSL_SUCCESS) throw std::runtime_error("PB13 direct-vacuum continuation failed");
        Result out{y[0],{}, {}}; for(int i=0;i<4;++i) {out.N[i]=y[7+i]*1e54;out.A[i]=y[11+i]*1e54;}
        return out;
    }
};
}
