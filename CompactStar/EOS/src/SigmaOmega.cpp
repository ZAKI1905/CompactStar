/*
  Sigma-Omega Model Class
*/

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include <Zaki/Physics/Constants.hpp>
#include <Zaki/Math/GSLFuncWrapper.hpp>

#include "CompactStar/EOS/SigmaOmega.hpp"
#include "CompactStar/EOS/Common.hpp"

using namespace Zaki::Physics ;
//==============================================================
//                   Sigma-Omega Model Class
//==============================================================
// Constructor
CompactStar::SigmaOmega::SigmaOmega() 
: Model({0.153, 10*0.153})
{
  SetName("SigmaOmega") ;
}

//--------------------------------------------------------------
/// Destructor
CompactStar::SigmaOmega::~SigmaOmega() { }

//--------------------------------------------------------------
void CompactStar::SigmaOmega::SetPars(const SigmaOmegaPar& in_pars)
{
  params = in_pars ;
}
//--------------------------------------------------------------
double CompactStar::SigmaOmega::I1(const double& in_rho, const double& x) 
{
  double k = pow(1.5*M_PI*M_PI*in_rho, 1./3.) ;
  double a = MN - x ;

  double out = k/2.*sqrt(k*k + a*a) - a*a/2.*atanh(k / sqrt(k*k + a*a));

  return out ;
}

//--------------------------------------------------------------
double CompactStar::SigmaOmega::I2(const double& in_rho, const double& x) 
{
  double k = pow(1.5*M_PI*M_PI*in_rho, 1./3.) ;
  double a = MN - x ;
  
  double out = -pow(a,3) * asinh(k/a) / sqrt(1+pow(k/a,2));
  out       += 2*pow(k,3) ;
  out       += pow(a,2)*k ;
  out       *= sqrt(a*a + k*k)/8. ;

  // double out = k*(2*k*k + a*a)*sqrt(k*k + a*a) ;
  // out += -pow(a,4)/8.*log(k + sqrt(k*k + a*a)) ;
  // out += pow(a,4)/8.*log(abs(a)) ;

  return out ;
}
//--------------------------------------------------------------
double CompactStar::SigmaOmega::I3(const double& in_rho, const double& x) 
{
  double k = pow(1.5*M_PI*M_PI*in_rho, 1./3.) ;
  double a = MN - x ;

  double out = k*(2*k*k - 3*a*a)*sqrt(k*k + a*a) ;
  out += 3.*pow(a,4)* atanh(k/sqrt(k*k + a*a)) ;
  out *= 1./8. ;

  return out ;
}

//--------------------------------------------------------------
double CompactStar::SigmaOmega::gsigEquation(double x)
{
  double out  = 2./M_PI/M_PI*(MN - x)*I1(rho, x) ;
  out  +=  - params.b*MN*pow(x, 2) ;
  out  +=  - params.c*pow(x, 3) ;
  out *= params.gsig_m_sqrd  ;
  out += -x ;
  
  return out ;
}

//--------------------------------------------------------------
double CompactStar::SigmaOmega::SolvegsigEquation()
{
  int status;
  int iter = 0, max_iter = 1000;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r = 0.5 ;
  double x_lo = 0.1, x_hi = 0.98*MN ;

  Zaki::Math::GSLFuncWrapper<SigmaOmega, double (SigmaOmega::*)(double)> 
  func(this, &SigmaOmega::gsigEquation) ;

  gsl_function F = static_cast<gsl_function> (func) ; 
  
  // F.function = &quadratic;
  // F.params = &params;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  // gsl_set_error_handler_off() ;
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);

  // printf ("using %s method\n",
  //         gsl_root_fsolver_name (s));

  // printf ("%5s [%9s, %9s] %9s %10s %9s\n",
  //         "iter", "lower", "upper", "root",
  //         "err", "err(est)");

  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      r = gsl_root_fsolver_root (s);
      x_lo = gsl_root_fsolver_x_lower (s);
      x_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi,
                                       1e-10, 1e-10);

      // if (status == GSL_SUCCESS)
      //   printf ("Converged:\n");

      // printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
      //         iter, x_lo, x_hi,
      //         r, r - r_expected,
      //         x_hi - x_lo);
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free (s);

  // std::cout<<"\n gsig= " << r << "\n" ;
  return r ;
}

//--------------------------------------------------------------
// Energy density (in g/cm^3 )
double CompactStar::SigmaOmega::EDens(const double& in_rho) 
{
  rho = in_rho ;
  double gsig = SolvegsigEquation() ;

  double out = params.b*MN*pow(gsig, 3)/3.;
  out += params.c*pow(gsig, 4)/4. ;
  out += pow(gsig,2)/2./params.gsig_m_sqrd ;
  out += params.gome_m_sqrd*pow(in_rho, 2)/2. ;
  out += 2./M_PI/M_PI*I2(in_rho, gsig);

  // double R = pow(3.*element.A/(4.*M_PI*in_rho), 1./3.);

  // // element.mu() is the atomic mass in GeV
  // double out = (element.mu() - element.Z*ELECTRON_M_GEV)*1e+3 
  //             - 0.9*element.Z*element.Z*Q_E_SQRD_MeVfm / R ;

  // out *= in_rho/element.A ;

  // // Converting the units
  // out *= MEV_FM3_2_G_CM3 ;

  // out += FermiEDens(pow(3.*M_PI*M_PI*in_rho*element.Z/element.A, 1./3.))*
  //         FM4_2_G_CM3 ;
  out *= INV_FM4_2_G_CM3 ;

  return out ;
}

//--------------------------------------------------------------
// Pressure
double CompactStar::SigmaOmega::Press(const double& in_rho) 
{
  rho = in_rho ;
  double gsig = SolvegsigEquation() ;

  double out = -params.b*MN*pow(gsig, 3)/3.;
  out += -params.c*pow(gsig, 4)/4. ;
  out += -pow(gsig,2)/2./params.gsig_m_sqrd ;
  out += params.gome_m_sqrd*pow(in_rho, 2)/2. ;
  out += 2./3./M_PI/M_PI*I3(in_rho, gsig);

  // double out = -0.3*pow(4.*M_PI/3., 1./3.)
  //               *element.Z*element.Z*Q_E_SQRD_MeVfm
  //               *pow(in_rho/element.A, 4./3.) ;

  // // Converting the units
  // out *= MEV_FM3_2_Dyn_CM2 ;

  // out += FermiPress(pow(3.*M_PI*M_PI*in_rho*element.Z/element.A, 1./3.))
  //       * FM4_2_Dyn_CM2 ;

  out *= INV_FM4_2_Dyn_CM2 ;

  return out ;
}

//--------------------------------------------------------------

//==============================================================
