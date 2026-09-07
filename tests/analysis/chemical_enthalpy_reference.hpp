// TEST ONLY. Independent full-star background + chemical integral reference.
// Reuses the PRE-EXISTING independent phase-space/enthalpy EOS test oracle;
// does not use Core TOV, profile interpolation, Geometry, H or chemical helpers.
#pragma once
#include <tests/eos/structure1/tov_oracle.hpp>
#include <array>
namespace chemical_reference {
struct Reference : structure1::EnthalpyOracle {
    std::array<double,8> response(double H) const {
        const auto q=at_H(H);
        std::array<double,8> out{};
        if (q.n==0) return out;
        std::array<double,4> chi{};
        for (int i=0;i<4;++i) {
            const double pf=eos.b*std::cbrt(3*M_PI*M_PI*q.ns[i]);
            chi[i]=std::hypot(eos.masses[i],pf)*pf/(M_PI*M_PI*std::pow(eos.b,3));
            out[4+i]=chi[i]; // OLD intrinsic route, separate accumulation.
        }
        const double den=chi[1]+chi[2]+chi[3];
        out[0]=chi[0];out[1]=chi[2]*(chi[1]+chi[3])/den;
        out[2]=-chi[2]*chi[3]/den;out[3]=chi[3]*(chi[1]+chi[2])/den;
        return out;
    }
    static int derivative(double x,const double y[],double f[],void* ptr) {
        auto& o=*static_cast<Reference*>(ptr);
        const double H=o.Hc-x;
        auto q=o.at_H(H);
        const double en=q.e*o.to_geo,p=q.p*o.to_geo;
        const double A=(1-2*y[0]*y[1])/(y[1]+4*M_PI*p);
        f[0]=2*A;f[1]=(4*M_PI*en-3*y[1])*A/y[0];
        auto c=o.response(H);
        const double measure=4*M_PI*std::sqrt(y[0])*std::exp(H)*A/std::sqrt(1-2*y[0]*y[1]);
        for (int i=0;i<8;++i) f[2+i]=measure*c[i];
        for (int i=0;i<10;++i) if (!std::isfinite(f[i])) return GSL_EBADFUNC;
        return GSL_SUCCESS;
    }
    std::array<double,10> compute(double tolerance,double start) {
        Hc=H_for_rho(1.10e15);
        const auto q=at_H(Hc);
        const double en=q.e*to_geo,p=q.p*to_geo;
        const double delta=start*Hc;
        double x=delta;
        double y[10]={3*delta/(2*M_PI*(en+3*p)),
            4*M_PI*en/3-4*M_PI*(en+p)*q.slope*delta/5};
        auto c=response(Hc);
        for (int i=0;i<8;++i) y[i+2]=4*M_PI/3*std::pow(y[0],1.5)*std::exp(Hc)*c[i];
        gsl_odeiv2_system sys={derivative,nullptr,10,this};
        auto* driver=gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_rkf45,
                                                  delta,tolerance*1e-4,tolerance);
        const int rc=gsl_odeiv2_driver_apply(driver,&x,Hc,y);
        gsl_odeiv2_driver_free(driver);
        if (rc!=GSL_SUCCESS) throw std::runtime_error("test-only enthalpy chemical reference failed");
        const double r=std::sqrt(y[0]),mass=y[1]*r*r*r;
        const double normalization=1/std::sqrt(1-2*mass/r);
        std::array<double,10> result{};result[0]=r;result[1]=mass;
        for (int i=0;i<8;++i) result[i+2]=y[i+2]*normalization*1e54;
        return result;
    }
};
}
