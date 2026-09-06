// ADR-0012 U1/U2/U5/U11/U13. Independent long-double literal formulas, not
// round trips alone. Negative conventions are test data, never runtime modes.
#include "CompactStar/RelativityUnits.hpp"
#include "CompactStar/Core/NStar.hpp"
#include "CompactStar/Core/TOVSolver.hpp"
#include "CompactStar/Physics/Evolution/StarContext.hpp"
#include "CompactStar/Physics/Driver/Thermal/Boundary/TbDefinition.hpp"
#include <gsl/gsl_version.h>
#include <gsl/gsl_errno.h>
#include <Zaki/Physics/Constants.hpp>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
namespace U = CompactStar::RelativityUnits;
using namespace CompactStar::Core;
constexpr long double G = 6.673e-8L, c = 29979245800.L, Ms = 1.98892e33L;
void require(bool p, const char* what) { if (!p) throw std::runtime_error(what); }
void close(double actual, long double expected, const char* what)
{
    require(std::isfinite(actual), what);
    require(expected == 0 ? actual == 0 : std::abs((static_cast<long double>(actual)-expected)/expected) <= 16*std::numeric_limits<double>::epsilon(), what);
}
int main()
{
 try {
    gsl_set_error_handler_off();
    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
    require(std::string(GSL_VERSION)=="2.7.1" && std::string(gsl_version)=="2.7.1", "U1 dependency vintage");
    require(GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT == 6.673e-8 && GSL_CONST_CGSM_SPEED_OF_LIGHT == 29979245800. && GSL_CONST_CGSM_SOLAR_MASS == 1.98892e33, "U1 TOV constants");
    std::cout << "U1 GSL header=" << GSL_VERSION << " linked=" << gsl_version
      << " G=" << GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT << " c=" << GSL_CONST_CGSM_SPEED_OF_LIGHT
      << " solar_g=" << GSL_CONST_CGSM_SOLAR_MASS << " negative_nonowner_SUN_M_KM=" << Zaki::Physics::SUN_M_KM << '\n';
    for (double x : {0., 1e-200, 1., 1e10, 1e15, 1e33, 1e100, 1e200}) {
      const long double a=x;
      close(U::MassGramsToKm(x), G*a/c/c/1e5L, "U2 mass forward");
      close(U::MassKmToGrams(U::MassGramsToKm(x)), a, "U2 mass round trip");
      close(U::MassKmToGrams(x), a*c*c/G*1e5L, "U2 mass inverse");
      close(U::MassDensityGcm3ToEnergyKmMinus2(x), G*a/c/c*1e10L, "U2 energy forward");
      close(U::EnergyKmMinus2ToMassDensityGcm3(U::MassDensityGcm3ToEnergyKmMinus2(x)), a, "U2 energy round trip");
      close(U::EnergyKmMinus2ToMassDensityGcm3(x), a*c*c/G/1e10L, "U2 energy inverse");
      close(U::PressureDynCm2ToKmMinus2(x), G*a/c/c/c/c*1e10L, "U2 pressure forward");
      close(U::PressureKmMinus2ToDynCm2(U::PressureDynCm2ToKmMinus2(x)), a, "U2 pressure round trip");
      close(U::PressureKmMinus2ToDynCm2(x), a*c*c*c*c/G/1e10L, "U2 pressure inverse");
      close(U::NuPrimeCmInverseToKmInverse(x), a*1e5L, "U2 nuprime forward");
      close(U::NuPrimeKmInverseToCmInverse(x), a/1e5L, "U2 nuprime inverse");
      close(U::NuPrimeKmInverseToCmInverse(U::NuPrimeCmInverseToKmInverse(x)), a, "U2 nuprime round trip");
      close(U::MassGramsToLiteralSolarMass(x), a/Ms, "U2 public mass");
      close(U::LiteralSolarMassToMassKm(x), G*a*Ms/c/c/1e5L, "U2 public to geometry");
      close(U::MassKmToLiteralSolarMass(x), a*c*c*1e5L/G/Ms, "U2 geometry to public");
    }
    std::cout << "U2 PASS (16 epsilon operation-count bound, independent long-double formulas)\n";
    // Exact Schwarzschild interior, independent analytic physical inputs.
    constexpr double M=2., R=13.; const double eps=3*M/(4*M_PI*R*R*R), ys=std::sqrt(1-2*M/R);
    std::vector<TOVPoint> points;
    for (int i=0;i<2001;++i) {
      double r=1e-5+(R-1e-5)*i/2000., y=std::sqrt(1-2*M*r*r/(R*R*R));
      double m=M*r*r*r/(R*R*R), p=eps*(y-ys)/(3*ys-y), np=2*M*r/(R*R*R*y*(3*ys-y));
      points.emplace_back(r, static_cast<double>(m*c*c*1e5L/G/Ms), np/1e5,0.,
         static_cast<double>(p*c*c*c*c/G/1e10L),static_cast<double>(eps*c*c/G/1e10L),.1,std::vector<double>{},0.);
    }
    NStar ns(points); const auto& p=ns.Profile();
    double worst=0.;
    for (size_t i=0;i<points.size();++i) {
      double r=points[i].r, m=M*r*r*r/(R*R*R), y=std::sqrt(1-2*M*r*r/(R*R*R));
      close((*p.GetMass())[i],m,"U5 analytic mass"); close((*p.GetEnergyDensity())[i],eps,"U5 analytic energy");
      close((*p.GetPressure())[i],eps*(y-ys)/(3*ys-y),"U5 analytic pressure");
      double F=((*p.GetMass())[i]+4*M_PI*r*r*r*(*p.GetPressure())[i])/(r*(r-2*(*p.GetMass())[i]));
      worst=std::max(worst,std::abs((*p.GetMetricNuPrime())[i]/F-1));
      require(std::abs((*p.GetMetricNu())[i]-std::log((3*ys-y)/2))<2e-7,"U5 analytic lapse (sampled integration envelope)");
    }
    require(worst<4e-15,"U5 exact nuprime identity");
    require(ns.MassSurface()==points.back().m && ns.GetSequence().ec==points.front().e && ns.GetSequence().pc==points.front().p,"U4 public physical bits");
    std::cout << "U5 analytic TOV PASS worst_nuprime=" << worst << '\n';
    CompactStar::Physics::Evolution::StarContext context(p);
    for(size_t i=0;i<points.size();++i) close((*context.MassDensity_gcm3())[i],points[i].e,"U11 context inverse");
    // Density-selected boundary with values within the old-vs-new inverse mismatch.
    // This synthetic column probes conversion only; it makes no new stellar-physics claim.
    auto profile=ns.Profile();
    { auto edit=profile.Edit(); auto &col=profile.RadialMutable()[profile.GetColumnIndex(StarProfile::Column::EnergyDensity)];
      for(size_t i=0;i<points.size();++i) col[i]=static_cast<double>((1.0001-i*0.0001)*1e10L*G/c/c*1e10L); }
    CompactStar::Physics::Evolution::StarContext layers(profile);
    auto index=CompactStar::Physics::Driver::Thermal::Boundary::FindTbIndex(layers,1e10);
    require(index==1,"U11 density layer");
    const double old_factor=Zaki::Physics::MEV_FM3_2_INV_KM2/Zaki::Physics::MEV_FM3_2_G_CM3;
    require((*layers.MassDensity_gcm3())[0]/((*profile.GetEnergyDensity())[0]/old_factor)-1>1e-4,"U13 wrong inverse detector");
    // Fixed coherent interior node. Each mixed field independently violates F or m'.
    double r=6., m=M*r*r*r/(R*R*R), y=std::sqrt(1-2*m/r), pressure=eps*(y-ys)/(3*ys-y);
    auto F=[r](double mass,double pr){return (mass+4*M_PI*r*r*r*pr)/(r*(r-2*mass));};
    const double beta=Zaki::Physics::SUN_M_KM/static_cast<double>(G*Ms/c/c/1e5L);
    const double alpha=old_factor/static_cast<double>(G/c/c*1e10L);
    require(std::abs(beta-1)>5e-5,"U13 nominal mass negative");
    require(std::abs(alpha-1)>1e-4,"U13 modern G negative");
    require(std::abs(beta/alpha-1)>2e-4,"U13 mixed mass equation");
    require(std::abs(F(m,pressure)/F(beta*m,alpha*pressure)-1)>1e-5,"U13 mixed nuprime");
    require(std::abs(F(m,pressure)/F(beta*m,pressure)-1)>1e-5,"U13 mass-only");
    require(std::abs(F(m,pressure)/F(m,alpha*pressure)-1)>1e-6,"U13 pressure-only");
    require(std::abs(1/alpha-1)>1e-4,"U13 energy-only / both density factors");
    require(std::abs((F(m,pressure)/1e5)/F(m,pressure)-1)>.99,"U13 carried nuprime length");
    std::cout << "U11 thermal inverse PASS; U13 all mixed controls detected; mass_residual=" << beta/alpha-1
      << " nuprime_residual=" << F(m,pressure)/F(beta*m,alpha*pressure)-1 << " wrong_inverse=" << 1/alpha-1 << '\n';
    return 0;
 } catch(const std::exception& e) {std::cerr<<"FAIL: "<<e.what()<<'\n'; return 1;}
}
