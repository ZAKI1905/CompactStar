/*
  Last edited on Aug 26, 2021
  RotationSolver class
*/
#include <gsl/gsl_math.h>
// #include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_integration.h>  // Aug 6, 2020
#include <gsl/gsl_odeiv2.h>       // Aug 6, 2020
// #include <gsl/gsl_const_cgsm.h>

#include <Zaki/Util/Instrumentor.hpp>
#include <Zaki/File/CSVIterator.hpp>
#include <Zaki/File/VecSaver.hpp>
#include <Zaki/Math/GSLFuncWrapper.hpp>
#include <Zaki/Physics/Constants.hpp>

#include "CompactStar/Core/RotationSolver.hpp"
#include "CompactStar/Core/NStar.hpp"
#include "CompactStar/Core/MixedStar.hpp"

#define R_SOLVER_VERBOSE 0

using namespace CompactStar ;
//==============================================================
//                        RotationSolver class
//==============================================================
// Constructor
RotationSolver::RotationSolver() : Prog("RotationSolver", true) 
{
  // Z_LOG_NOTE("Constructor Called!") ;
  // accel  = gsl_interp_accel_alloc ();
}

//--------------------------------------------------------------
RotationSolver::~RotationSolver()
{
  // Z_LOG_NOTE("Destructor Called!") ;

  // std::stringstream ss;

  // if(m_spline)
  // { 
  //   // ss << "*m_spline = " << static_cast<const void*>(m_spline);
  //   gsl_spline_free (m_spline) ;
  // }
  // if(e_spline)
  // {
  //   // ss << "\n*e_spline = " << static_cast<const void*>(e_spline);
  //   gsl_spline_free (e_spline) ;
  // }
  // if(p_spline)
  // {
  //   // ss << "\n*p_spline = " << static_cast<const void*>(p_spline);
  //   gsl_spline_free (p_spline) ;
  // }

  // if(accel)
  // {
  //   // ss << "\n*accel = " << static_cast<const void*>(accel);
  //   gsl_interp_accel_free (accel) ;
  // }

  // std::cout <<  ss.str() << "\n" ;
}

//--------------------------------------------------------------
// void RotationSolver::ImportTOVSolution(const TOVTable& table) 
// {

//   tov_solution = table ;

//   Z_LOG_INFO("TOV solution imported.") ;

//   m_spline  = gsl_spline_alloc (gsl_interp_cspline, tov_solution.Size());
//   p_spline  = gsl_spline_alloc (gsl_interp_cspline, tov_solution.Size()); 
//   e_spline  = gsl_spline_alloc (gsl_interp_cspline, tov_solution.Size()); 

//   // This function initializes the interpolation object
//   // x has to be strictly increasing
//   gsl_spline_init (m_spline, &tov_solution.r[0], &tov_solution.m[0], 
//                     tov_solution.Size());
//   gsl_spline_init (p_spline, &tov_solution.r[0], &tov_solution.pre[0], 
//                     tov_solution.Size());
//   gsl_spline_init (e_spline, &tov_solution.r[0], &tov_solution.eps[0], 
//                     tov_solution.Size());
  
//   R_Star = tov_solution.r[tov_solution.Size()-1] ;
//   M_Star = tov_solution.m[tov_solution.Size()-1] ;
// }

//--------------------------------------------------------------
// void RotationSolver::ImportTOVSolution(const Zaki::String::Directory& f_name) 
// {
//   std::ifstream     file( (wrk_dir + "/" + f_name).Str());

//   // Error opening the file
//   if (file.fail()) 
//   {
//     Z_LOG_ERROR("File '"+(wrk_dir + "/" + f_name).Str() +"' cannot be opened!") ;
//     Z_LOG_ERROR("Importing TOV solution failed!") ;
//     exit(EXIT_FAILURE) ;
//     return ;
//   }

//   size_t line_num = 0 ;
//   // Reading the input file
//   for(Zaki::File::CSVIterator loop(file, '\t'); loop != Zaki::File::CSVIterator(); ++loop)
//   {
//     if( (*loop).size() < 4)
//     {
//       Z_LOG_ERROR("TOV solution file is not complete!") ;
//       break ;
//     }

//     // First line
//     if (line_num == 0)
//     {     }
//     else // Aug 6, 2020: The issue is that we change the format of the files on December 2020!
//     {
//       // Radius in km
//       tov_solution.r.push_back(std::atof((*loop)[0].c_str())) ;

//       // Mass (in M_sun), which we multiply by M_Sun = 1.47663 km
//       tov_solution.m.push_back(std::atof((*loop)[1].c_str())*Zaki::Physics::SUN_M_KM) ;

//       // Pressure (dyne/cm^2), which we convert to 1/km^2
//       // tov_solution.pre.push_back(std::atof((*loop)[3].c_str())*8.26006*1e-40) ;
//       tov_solution.pre.push_back( std::atof((*loop)[3].c_str()) *
//                                   Zaki::Physics::INV_FM4_2_INV_KM2 /
//                                   Zaki::Physics::INV_FM4_2_Dyn_CM2 ) ;

//       // Energy density (g/cm^3), which we convert to 1/km^2
//       // tov_solution.eps.push_back(std::atof((*loop)[4].c_str())*7.42367*1e-19) ;
//       tov_solution.eps.push_back( std::atof( (*loop)[4].c_str()) * 
//                                   Zaki::Physics::INV_FM4_2_INV_KM2 / 
//                                   Zaki::Physics::INV_FM4_2_G_CM3 ) ;
//     }
//     line_num++ ;
//   }

//   Z_LOG_INFO("TOV solution imported from: "+ (wrk_dir+f_name).Str()+".") ;

//   m_spline  = gsl_spline_alloc (gsl_interp_cspline, tov_solution.Size());
//   p_spline  = gsl_spline_alloc (gsl_interp_cspline, tov_solution.Size()); 
//   e_spline  = gsl_spline_alloc (gsl_interp_cspline, tov_solution.Size()); 

//   // This function initializes the interpolation object
//   // x has to be strictly increasing
//   gsl_spline_init (m_spline, &tov_solution.r[0], &tov_solution.m[0], 
//                     tov_solution.Size());
//   gsl_spline_init (p_spline, &tov_solution.r[0], &tov_solution.pre[0], 
//                     tov_solution.Size());
//   gsl_spline_init (e_spline, &tov_solution.r[0], &tov_solution.eps[0], 
//                     tov_solution.Size());
  
//   R_Star = tov_solution.r[tov_solution.Size()-1] ;
//   M_Star = tov_solution.m[tov_solution.Size()-1] ;
// }
//--------------------------------------------------------------
void RotationSolver::AttachNStar(NStar* ns_ptr)
{
  Z_LOG_INFO("Attaching NStar to RotationSolver class...") ;

  if(!ns_ptr)
  {
    Z_LOG_ERROR("The NStar pointer is a nullptr!") ;
    return ;
  }

  nstar_ptr = ns_ptr ;

}

//--------------------------------------------------------------
void RotationSolver::AttachMixedStar(MixedStar* m_ns_ptr)
{
  Z_LOG_INFO("Attaching MixedStar to RotationSolver class...") ;

  if(!m_ns_ptr)
  {
    Z_LOG_ERROR("The MixedStar pointer is a nullptr!") ;
    return ;
  }

  mixedstar_ptr = m_ns_ptr ;

  // if( !(mixedstar_ptr->nu_r_spline) )
  //   mixedstar_ptr->EvaluateNu() ;
}

//--------------------------------------------------------------
// Input radius (in km), output mass (in km)
// double RotationSolver::GetMass(const double& in_R) 
// {
//   if ( in_R < nstar_ptr->radius[-1] )
//   {
//     // return gsl_spline_eval(m_spline, in_R, accel);
//     return nstar_ptr->GetMass(in_R) ;
//   }
//   // outside the star
//   else
//     return nstar_ptr->mass[-1] ;
// }

//--------------------------------------------------------------
// Input radius (in km), output pressure (in km^-2)
// double RotationSolver::GetPress(double in_R) 
// {
//   // return gsl_spline_eval(p_spline, in_R, accel);
//   return nstar_ptr->GetPress(in_R) ;
// }

//--------------------------------------------------------------
// Input radius (in km), outputs the derivative of pressure (km^-3)
// double RotationSolver::GetPressDer(const double& in_R) 
// {
//   Zaki::Math::GSLFuncWrapper<RotationSolver, double (RotationSolver::*)(double)> 
//     Fp(this, &RotationSolver::GetPress);     

//   gsl_function F = static_cast<gsl_function> (Fp) ; 

//   double result, abserr ;
//   gsl_deriv_forward (&F, in_R, 1e-5, &result, &abserr); // Aug 6, 2020
//   // std::cout << "P' (" << in_R << ") = " << result << " +/- " << abserr << "\n";

//   return result;
// }

//--------------------------------------------------------------
// Input radius, output energy density
double RotationSolver::GetEDens(const double& in_R) 
{
  // return gsl_spline_eval(e_spline, in_R, accel);
  return nstar_ptr->GetEps(in_R) ;
}

//--------------------------------------------------------------
// nu(r) integrand (Aug 6, 2020)
// double RotationSolver::NuIntegrand(double r) 
// {
//   //         km^-3       /   km^-2
//   return -GetPressDer(r) / ( GetPress(r) + GetEDens(r) ); 
// }

//--------------------------------------------------------------
// Input radius, output metric function nu(r) (Aug 6, 2020)
// [ Aug 26, 2021 ]: The output is only valid inside the star
// for outside the star this should simply return:
//  0.5*log( 1 - 2*M_Star/R_Star)
//  !!! Radius is in km here !!!
//  Changed on Mar 23, 2022 !
double RotationSolver::GetNu(double in_R) 
{
  return nstar_ptr->GetNu(in_R) ;

  // // double r_star = tov_solution.r[tov_solution.Size()-1] ;
  // // double m_star = tov_solution.m[tov_solution.Size()-1] ;

  // double result, err, integ_const ;

  // if (in_R < nstar_ptr->radius[-1] )
  // {
  //   gsl_integration_workspace *w = 
  //     gsl_integration_workspace_alloc(2000);

  //   Zaki::Math::GSLFuncWrapper<RotationSolver, double (RotationSolver::*)(double)> 
  //     Fp(this, &RotationSolver::NuIntegrand);     

  //   gsl_function F = static_cast<gsl_function> (Fp) ; 

  //   gsl_integration_qag(&F, tov_solution.r[0], R_Star, 1e-8, 1e-8, 1000, 1, w, &integ_const, &err);

  //   gsl_integration_qag(&F, tov_solution.r[0], in_R, 1e-8, 1e-8, 1000, 1, w, &result, &err);

  //   result = result - integ_const + 0.5*log( 1 - 2*M_Star/R_Star);
  // }
  // else
  // {
  //   result = 0.5*log( 1 - 2*M_Star/in_R) ;
  //   integ_const = -1 ; err = 0 ;
  // }
  // // Printing for debugging purposes:
  // // std::cout << "nu (" << in_R << ") = " << result << " +/- " << err 
  // //           << ", C = " << integ_const <<"\n";

  // return result ;
}

//--------------------------------------------------------------
// Added on Aug 6, 2020
// Returns the value of the j function (Hartle Eq. 6.17, P. 253)
// double RotationSolver::GetHartleJ(double r) 
// {
//   if ( r < nstar_ptr->radius[-1] )
//   {
//     return exp(-GetNu(r)) * sqrt( 1 - 2* GetMass(r) / r) ;
//   }
//   else
//     return 1 ;

// }

//--------------------------------------------------------------
// Added on Aug 6, 2020
// Returns the derivative of the j function (Hartle Eq. 6.18, P. 253)
// double RotationSolver::GetHartleJDer(double r) 
// {
//   if ( r < nstar_ptr->radius[-1] )
//   {
//     double out = exp(-GetNu(r)) / sqrt( 1. - 2.* GetMass(r) / r)  ;
//     out *= -4. * M_PI * r * (GetPress(r) + GetEDens(r)) ;
//     return out ;
//   }
//   else
//     return 0 ;

// }

//--------------------------------------------------------------
// Added on Aug 6, 2020
// int RotationSolver::ODE(double r, const double y[], double f[], void *params)
// {
//   RotationSolver* rot_obj = (RotationSolver *)params ; 

//   // First derivative of bar{omega(r)}
//   f[0] = y[1];

//   // Second derivative of bar{omega(r)}
//   // f[1] = -4. * y[0] * rot_obj->GetHartleJDer(r) / r / rot_obj->GetHartleJ(r) 
//   //       - (4./r + rot_obj->GetHartleJDer(r) / rot_obj->GetHartleJ(r)) * y[1];

//   // Second derivative of bar{omega(r)}
//   // A better way!
//   f[1] = -4.*y[1]/r 
//         + rot_obj->GetHartleOmegaCoeff(r) * y[0] 
//         + y[1] * rot_obj->GetHartleDOmegaCoeff(r) ;

//   // if ( r < rot_obj->GetNStar()->radius[-1])
//   // {
//   //   f[1] = -4.*y[1]/r 
//   //         + 4.*M_PI*(rot_obj->GetPress(r) + rot_obj->GetEDens(r)) 
//   //           * (r * y[1] + 4. * y[0])
//   //         / ( 1. - 2.*rot_obj->GetMass(r) / r );
//   // }
//   // else
//   // {
//   //   f[1] = -4.*y[1]/r ;
//   // }
//   return GSL_SUCCESS;

// }

//--------------------------------------------------------------
// Added on Apr 20, 2022
// Edited on May 3, 2022
int RotationSolver::ODE_Mixed(double r, const double y[], double f[], void *params)
{
  RotationSolver* rot_obj = (RotationSolver *)params ; 

  // First derivative of bar{omega(r)}
  f[0] = y[1];

  // Second derivative of bar{omega(r)}
  f[1] = -4.*y[1]/r 
        + rot_obj->GetHartleOmegaCoeff_Mixed(r) * y[0]
        + y[1] * rot_obj->GetHartleDOmegaCoeff_Mixed(r) ;
  
//  std::cout << "r = " << r << "\n" << std::flush ;
  
  return GSL_SUCCESS;

}

//--------------------------------------------------------------
// Created on June 3, 2022
int RotationSolver::ODE_N_Fast(double r, const double y[], double f[], void *params)
{
  RotationSolver* rot_obj = (RotationSolver *)params ; 

  // First derivative of bar{omega(r)}
  f[0] = y[1];

  // Second derivative of bar{omega(r)}
  f[1] = -4.*y[1]/r 
        + rot_obj->GetHartleOmegaCoeff_N_Fast(r) * y[0]
        + y[1] * rot_obj->GetHartleDOmegaCoeff_N_Fast(r) ;
  
  
  return GSL_SUCCESS;
}

//--------------------------------------------------------------
// Created on May 3, 2022
int RotationSolver::ODE_Mixed_Fast(double r, const double y[], double f[], void *params)
{
  RotationSolver* rot_obj = (RotationSolver *)params ; 

  // First derivative of bar{omega(r)}
  f[0] = y[1];

  // Second derivative of bar{omega(r)}
  f[1] = -4.*y[1]/r 
        + rot_obj->GetHartleOmegaCoeff_Mixed_Fast(r) * y[0]
        + y[1] * rot_obj->GetHartleDOmegaCoeff_Mixed_Fast(r) ;
  
//  std::cout << "r = " << r << "\n" << std::flush ;
  
  return GSL_SUCCESS;

}

//--------------------------------------------------------------
// Created on May 3, 2022
int RotationSolver::ODE_Mixed_Out(double r, const double y[], double f[], void *params)
{
  RotationSolver* rot_obj = (RotationSolver *)params ;
  
    // First derivative of bar{omega(r)}
  f[0] = y[1] ;
  
    // Second derivative of bar{omega(r)}
  f[1] = -4.*y[1]/r ;
  
  
  return GSL_SUCCESS;
  
}


//--------------------------------------------------------------
/// Returns the total baryon number density given pressure
// double RotationSolver::GetRho(const double& in_p) 
// {
//   return gsl_spline_eval(rho_spline, in_p, accel) ;
// }

//--------------------------------------------------------------
// Coefficient for y[0]
// double RotationSolver::GetHartleOmegaCoeff(const double r)
// {
//   if ( r < GetNStar()->radius[-1])
//   {
//     return  16.*M_PI*(GetPress(r) + GetEDens(r)) 
//             / ( 1. - 2.*GetMass(r) / r );
//   }
//   else
//   {
//     return 0 ;
//   }
// }

//--------------------------------------------------------------
// // Coefficient for y[1]
// double RotationSolver::GetHartleDOmegaCoeff(const double r)
// {
//   if ( r < GetNStar()->radius[-1])
//   {
//     return  4.*M_PI*(GetPress(r) + GetEDens(r)) * r
//            / ( 1. - 2.*GetMass(r) / r ) ;
//   }
//   else
//   {
//     return 0 ;
//   }
// }

//--------------------------------------------------------------
// MixedStar: coefficient for y[0]
double RotationSolver::GetHartleOmegaCoeff_Mixed(const double r)
{
  return  16.*M_PI*(
                    mixedstar_ptr->GetPress_Visible(r)
                  + mixedstar_ptr->GetPress_Dark(r)
                  + mixedstar_ptr->GetEps_Visible(r)
                  + mixedstar_ptr->GetEps_Dark(r)
                   )
          / ( 1. - 2.*mixedstar_ptr->GetMass_Total(r) / r );
  
}

//--------------------------------------------------------------
// MixedStar: coefficient for y[1]
double RotationSolver::GetHartleDOmegaCoeff_Mixed(const double r)
{
  return  4.*M_PI*(
                    mixedstar_ptr->GetPress_Visible(r)
                  + mixedstar_ptr->GetPress_Dark(r)
                  + mixedstar_ptr->GetEps_Visible(r)
                  + mixedstar_ptr->GetEps_Dark(r)
                  ) * r
         / ( 1. - 2.*mixedstar_ptr->GetMass_Total(r) / r ) ;
}

//--------------------------------------------------------------
// MixedStar: coefficient for y[0]
double RotationSolver::GetHartleOmegaCoeff_N_Fast(const double r)
{
  return  16.*M_PI*( fast_p + fast_e )
                  / ( 1. - 2.* fast_m / r );
}

//--------------------------------------------------------------
// MixedStar: coefficient for y[1]
double RotationSolver::GetHartleDOmegaCoeff_N_Fast(const double r)
{
  return  4.*M_PI*( fast_p + fast_e ) * r
         / ( 1. - 2.* fast_m / r ) ;
}

//--------------------------------------------------------------
// MixedStar: coefficient for y[0]
double RotationSolver::GetHartleOmegaCoeff_Mixed_Fast(const double r)
{
  return  16.*M_PI*( fast_p_v + fast_p_d + fast_e_v + fast_e_d )
                  / ( 1. - 2.* fast_m_tot / r );
}

//--------------------------------------------------------------
// MixedStar: coefficient for y[1]
double RotationSolver::GetHartleDOmegaCoeff_Mixed_Fast(const double r)
{
  return  4.*M_PI*( fast_p_v + fast_p_d + fast_e_v + fast_e_d ) * r
         / ( 1. - 2.* fast_m_tot / r ) ;
}

//--------------------------------------------------------------
// bar{omega(r)} not omega
double RotationSolver::GetInitOmegaBar() const 
{
  return  init_omega_bar ;
}

//--------------------------------------------------------------
// // Originally added on Aug 6, 2020
// void RotationSolver::Solve(const Zaki::Math::Axis& in_ax,
//                       const Zaki::String::Directory& in_dir)  
// {

//   // Resolution of the solver (division of the radial distance)
//   // unsigned int solver_res = 5000 ;

//   double r_min = nstar_ptr->radius[0] ;
//   double r_max = 1.1 * (nstar_ptr->radius[-1]) ;
//   // double Rstar = tov_solution.r[tov_solution.Size()-1] ;

//   std::cout << "\n\n\t\t ****************************************"
//             <<"*******************************"<<" \n" ;
//   std::cout <<     "\t\t *                   "
//             <<"Rotation Solver Sequence Results"
//             <<"                  * \n" ;
//   std::cout <<     "\t\t ******************************************"
//             <<"*****************************"<<"\n\n" ;

//   for (size_t idx = 0; idx <= in_ax.res ; idx++)
//   {
//     init_omega_bar = in_ax[idx] ;
//     double r = r_min ; 
//     // double r1 = r_max ;

//     double y[2] ; 
//     y[0] = init_omega_bar ;
//     y[1] = 0 ;

//     gsl_odeiv2_system ode_sys = {CompactStar::RotationSolver::ODE, nullptr, 2, this} ;

//     gsl_odeiv2_driver *tmp_driver = gsl_odeiv2_driver_alloc_y_new
//       (&ode_sys, gsl_odeiv2_step_rk8pd,
//         1.e-1, 1.e-10, 1.e-10);


//     double min_log_r  = log10(r_min) ;
//     double max_log_r  = log10(r_max) ;
//     // -------------------------------------------------
//     // Added on March 23, 2022 to fix the bug 
//     // in J, Omega, I calculations.
//     double log_r_star = log10(nstar_ptr->radius[-1]) ;
    
//     // Total angular momentum and velocity:
//     double ang_mom_J, ang_vel_Omega, mom_inertia ;
//     // -------------------------------------------------


//     double step = (max_log_r - min_log_r) / radial_res ;

//     // One extra point for the exact value at r = nstar_ptr->radius[-1].
//     omega_results.reserve(radial_res+1) ;
//     bool surface_reached = false ;

//     for (double log_r_i = min_log_r;  log_r_i <= max_log_r; log_r_i += step)
//     {
//       double ri = pow(10, log_r_i) ;

//       // This function evolves the driver system d from t to t1.
//       // Initially vector y should contain the values of dependent
//       // variables at point t. 
//       int status = gsl_odeiv2_driver_apply (tmp_driver, &r, ri, y);

//       if (status != GSL_SUCCESS)
//       {
//         printf ("error, return value=%d\n", status);
//         break;
//       }

//       // Units : { [ km ], [ M_Sun ], - , [ km^-1 ], [ km^-2 ] }
//       omega_results.emplace_back(r, GetMass(r) / Zaki::Physics::SUN_M_KM, GetNu(r), y[0], y[1]);

//       // March 23, 2022
//       if ( !surface_reached && abs(log_r_i - log_r_star) < step ) 
//       {
//         // std::cout << "R_Star = "<< R_Star << ", r = " 
//         //           << r << ", r_i = " << ri << ", log_r_i - log_r_star = "
//         //           << log_r_i - log_r_star << ", Step = " << step 
//         //           << ", y[0] = " << y[0]<< ", y[1] = " << y[1] << "\n" ;

//         status = gsl_odeiv2_driver_apply (tmp_driver, &r, nstar_ptr->radius[-1], y) ;
//         omega_results.emplace_back(r, GetMass(r) / Zaki::Physics::SUN_M_KM, GetNu(r), y[0], y[1]) ;

//         // std::cout << "R_Star = "<< R_Star << ", r = " 
//         //           << r << ", r_i = " << ri << "\n" ;
//         // Total angular momentum and velocity:
//         ang_mom_J = pow(nstar_ptr->radius[-1],4)*y[1]/6. ;
//         ang_vel_Omega = y[0] + nstar_ptr->radius[-1]*y[1]/3. ;
//         mom_inertia = ang_mom_J / ang_vel_Omega ;

//         // std::cout << "ang_mom_J = "<< ang_mom_J << ", ang_vel_Omega = " 
//         //           << ang_vel_Omega << ", mom_inertia = " << mom_inertia << "\n" ;
//         surface_reached = true ;
//       }
//     }
//     // ............ Radius Loop Ends .........................
    


//     // ............................. omega_results Loop Begins .....................................
//     // We go back now and update the omega_results with the values for omega(r)
//     // We needed ang_vel_Omega for omega(r) since:
//     // omega(r) = ang_vel_Omega - bar{omega(r)}
//     // The unit for omega(r) is changed from km^{-1} to s^{-1}
//     for (size_t i = 0; i < omega_results.size(); i++)
//     {
//       omega_results[i].omega = (ang_vel_Omega - omega_results[i].omega_bar)*Zaki::Physics::LIGHT_C_KM_S ;
//     }
//     // .............................. omega_results Loop Ends ......................................


//     // ...............................................................................................
//     // Printing the sequence point data in terminal
//     printf ("\u2554%-40s \u03A9-Sequence %-3zu out of %-3lu %43s\u2557\n",
//     Zaki::String::Multiply("\u2550", 40).c_str(), idx+1, in_ax.res+1, 
//     Zaki::String::Multiply("\u2550", 43).c_str());
//     printf ("\u2551 %14s %14s %17s %16s %14s %12s %12s %9s\n", "R (km)", "M (km)", "\u03B5_c (km^-2)", 
//     "\u03C9c_bar (s^-1)", "\u03A9 (s^-1) ", "J (km^2)", "I (km^3)", "\u2551");
//     printf ("\u2551 %14le %14le %14le %14le %14le %14le %14le %7s\n",
//     nstar_ptr->radius[-1],
//     nstar_ptr->mass[-1],
//     nstar_ptr->eps[0], 
//     init_omega_bar*Zaki::Physics::LIGHT_C_KM_S,
//     ang_vel_Omega*Zaki::Physics::LIGHT_C_KM_S,
//     ang_mom_J,
//     mom_inertia,
//     "\u2551");
//     printf ("\u255A%110s\u255D\n", Zaki::String::Multiply("\u2550", 110).c_str());

//     // ...............................................................................................

//     omega_seq_pts.emplace_back(init_omega_bar*Zaki::Physics::LIGHT_C_KM_S, 
//                                 nstar_ptr->mass[-1], 
//                                 nstar_ptr->radius[-1], ang_mom_J, ang_vel_Omega*Zaki::Physics::LIGHT_C_KM_S );

//     gsl_odeiv2_driver_free (tmp_driver);

//     // ------------------------------------------------------------
//     //                      Saving to file
//     // ------------------------------------------------------------
//     if(in_dir.Str() != "")
//     {
//       std::vector<std::string> tmp_fname_v = Zaki::String::Pars(in_dir.Str(), "*") ;
//       std::string tmp_fname ;
//       if(tmp_fname_v.size() <2 )
//       {
//         Z_LOG_WARNING("File name patter doesn't match '[]*[]'.") ;
//         tmp_fname = tmp_fname_v[0] + std::to_string(idx) ;
//       }
//       else
//       {
//         tmp_fname = tmp_fname_v[0] + std::to_string(idx) + tmp_fname_v[1] ;
//       }
      
//       Zaki::File::VecSaver vec_saver(wrk_dir + "/" + tmp_fname) ;

//       // .............................................
//       // Header
//       std::stringstream ss ;
//       char res_header[300] ;

//       sprintf (res_header, "\u2554%-40s Sequence %-3zu out of %-3lu %45s\u2557\n",
//       Zaki::String::Multiply("\u2550", 40).c_str(), idx+1, in_ax.res+1, 
//       Zaki::String::Multiply("\u2550", 45).c_str());
//       ss << "# " << res_header ;

//       sprintf (res_header, "\u2551 %14s %12s %17s %16s %14s %12s %12s %11s\n", "R (km)", "M (Sun)", "\u03B5_c (g/cm^3)", 
//       "\u03C9c_bar (s^-1)", "\u03A9 (s^-1) ", "J (km^2)", "I (km^3)", "\u2551");
//       ss << "# " << res_header ;

//       sprintf (res_header, "\u2551 %14le %14le %14le %14le %14le %14le %14le %7s\n",
//       nstar_ptr->radius[-1],
//       nstar_ptr->mass[-1],
//       nstar_ptr->eps[0], 
//       init_omega_bar*Zaki::Physics::LIGHT_C_KM_S,
//       ang_vel_Omega*Zaki::Physics::LIGHT_C_KM_S,
//       ang_mom_J,
//       mom_inertia,
//       "\u2551");
//       ss << "# " << res_header ;

//       sprintf (res_header, "\u255A%110s\u255D\n", Zaki::String::Multiply("\u2550", 110).c_str());
//       ss << "# " << res_header ;

//       sprintf(res_header, "%-19s\t %-19s\t %-19s\t %-19s\t %-19s\t %-19s", 
//               "r [km]", "M [Sun]", "nu", "omega [1/s]", "omega_bar [1/km]",  "omega'_bar [1/km^2]" ) ;

//       ss << res_header ;

//       // std::string tmp_label = res_header ;
//       // tmp_label.c_str()
//       // .............................................

//       vec_saver.SetHeader(ss.str().c_str()) ;
//       vec_saver.Export1D(omega_results) ;
//     }
//     // ------------------------------------------------------------
    
//     // Testing this for freeing memory
//     // omega_results = std::vector<OmegaPoint>() ;
//     // instead of
//     omega_results.clear() ;
//     omega_results.shrink_to_fit() ;
//   } // omega sequence loop ends!
//   std::cout <<"\n\t\t ************************"
//             <<" Rot. Solver Finished *************************"<<"\n\n" ;
//   // ------------------------------------------------------------

// }
//--------------------------------------------------------------
// Added on Mar 23, 2022
// Input an initial omega(0) value
// void RotationSolver::Solve( const double& in_omega_0, 
//                             const Zaki::String::Directory& in_dir) 
// {
//   PROFILE_FUNCTION() ;

//   // Resolution of the solver (division of the radial distance)
//   // unsigned int solver_res = 5000 ;

//   double r_min = nstar_ptr->radius[0] ;
//   double r_max = 1.1*nstar_ptr->radius[-1] ;
//   // double Rstar = tov_solution.r[tov_solution.Size()-1] ;

//   std::cout << "\n\n\t\t ****************************************"
//             <<"*******************************"<<" \n" ;
//   std::cout <<     "\t\t *                   "
//             <<"Rotation Solver Sequence Results"
//             <<"                  * \n" ;
//   std::cout <<     "\t\t ******************************************"
//             <<"*****************************"<<"\n\n" ;


//   init_omega_bar = in_omega_0 ;
//   double r = r_min ; 
//   // double r1 = r_max ;

//   double y[2] ; 
//   y[0] = init_omega_bar ;
//   y[1] = 0 ;

//   gsl_odeiv2_system ode_sys = {CompactStar::RotationSolver::ODE, nullptr, 2, this} ;

//   gsl_odeiv2_driver *tmp_driver = gsl_odeiv2_driver_alloc_y_new
//     (&ode_sys, gsl_odeiv2_step_rk8pd,
//       1.e-1, 1.e-10, 1.e-10);


//   double min_log_r  = log10(r_min) ;
//   double max_log_r  = log10(r_max) ;
//   // -------------------------------------------------
//   // Added on March 23, 2022 to fix the bug 
//   // in J, Omega, I calculations.
//   double log_r_star = log10(nstar_ptr->radius[-1]) ;
  
//   // Total angular momentum and velocity:
//   double ang_mom_J, ang_vel_Omega, mom_inertia ;
//   // -------------------------------------------------


//   double step = (max_log_r - min_log_r) / radial_res ;

//   // One extra point for the exact value at r = R_star.
//   omega_results.reserve(radial_res+1) ;
//   bool surface_reached = false ;

//   for (double log_r_i = min_log_r;  log_r_i <= max_log_r; log_r_i += step)
//   {
//   //   PROFILE_SCOPE("R-Loop") ;

//     double ri = pow(10, log_r_i) ;

//     // This function evolves the driver system d from t to t1.
//     // Initially vector y should contain the values of dependent
//     // variables at point t. 
//     int status = gsl_odeiv2_driver_apply (tmp_driver, &r, ri, y);

//     if (status != GSL_SUCCESS)
//     {
//       printf ("error, return value=%d\n", status);
//       break;
//     }

//     // Units : { [ km ], [ M_Sun ], - , [ km^-1 ], [ km^-2 ] }
//     omega_results.emplace_back(r, GetMass(r) / Zaki::Physics::SUN_M_KM, GetNu(r), y[0], y[1]);

//     // March 23, 2022
//     if ( !surface_reached && abs(log_r_i - log_r_star) < step ) 
//     {
//       // std::cout << "nstar_ptr->radius[-1] = "<< nstar_ptr->radius[-1] << ", r = " 
//       //           << r << ", r_i = " << ri << ", log_r_i - log_r_star = "
//       //           << log_r_i - log_r_star << ", Step = " << step 
//       //           << ", y[0] = " << y[0]<< ", y[1] = " << y[1] << "\n" ;

//       status = gsl_odeiv2_driver_apply (tmp_driver, &r, nstar_ptr->radius[-1], y) ;
//       omega_results.emplace_back(r, GetMass(r) / Zaki::Physics::SUN_M_KM, GetNu(r), y[0], y[1]) ;

//       // std::cout << "nstar_ptr->radius[-1] = "<< nstar_ptr->radius[-1] << ", r = " 
//       //           << r << ", r_i = " << ri << "\n" ;
//       // Total angular momentum and velocity:
//       ang_mom_J = pow(nstar_ptr->radius[-1],4)*y[1]/6. ;
//       ang_vel_Omega = y[0] + nstar_ptr->radius[-1]*y[1]/3. ;
//       mom_inertia = ang_mom_J / ang_vel_Omega ;

//       // std::cout << "ang_mom_J = "<< ang_mom_J << ", ang_vel_Omega = " 
//       //           << ang_vel_Omega << ", mom_inertia = " << mom_inertia << "\n" ;
//       surface_reached = true ;
//     }
//   }
//   // ............ Radius Loop Ends .........................
    


//   // ............................. omega_results Loop Begins .....................................
//   // We go back now and update the omega_results with the values for omega(r)
//   // We needed ang_vel_Omega for omega(r) since:
//   // omega(r) = ang_vel_Omega - bar{omega(r)}
//   // The unit for omega(r) is changed from km^{-1} to s^{-1}
//   for (size_t i = 0; i < omega_results.size(); i++)
//   {
//     omega_results[i].omega = (ang_vel_Omega - omega_results[i].omega_bar)*Zaki::Physics::LIGHT_C_KM_S ;
//   }
//   // .............................. omega_results Loop Ends ......................................


//   // ...............................................................................................
//   // Printing the sequence point data in terminal
//   printf ("\u2554%-40s \u03A9-Sequence %-3zu out of %-3lu %43s\u2557\n",
//   Zaki::String::Multiply("\u2550", 40).c_str(), (size_t)1, (size_t)1, 
//   Zaki::String::Multiply("\u2550", 43).c_str());
//   printf ("\u2551 %14s %14s %17s %16s %14s %12s %12s %9s\n", "R (km)", "M (km)", "\u03B5_c (km^-2)", 
//   "\u03C9c_bar (s^-1)", "\u03A9 (s^-1) ", "J (km^2)", "I (km^3)", "\u2551");
//   printf ("\u2551 %14le %14le %14le %14le %14le %14le %14le %7s\n",
//   nstar_ptr->radius[-1],
//   nstar_ptr->mass[-1],
//   nstar_ptr->eps[0], 
//   init_omega_bar*Zaki::Physics::LIGHT_C_KM_S,
//   ang_vel_Omega*Zaki::Physics::LIGHT_C_KM_S,
//   ang_mom_J,
//   mom_inertia,
//   "\u2551");
//   printf ("\u255A%110s\u255D\n", Zaki::String::Multiply("\u2550", 110).c_str());

//   // ...............................................................................................

//   omega_seq_pts.emplace_back(init_omega_bar*Zaki::Physics::LIGHT_C_KM_S, 
//                               nstar_ptr->mass[-1], 
//                               nstar_ptr->radius[-1], ang_mom_J, ang_vel_Omega*Zaki::Physics::LIGHT_C_KM_S );

//   gsl_odeiv2_driver_free (tmp_driver) ;

//   // ------------------------------------------------------------
//   //                      Saving to file
//   // ------------------------------------------------------------
//   if(false && in_dir.Str() != "")
//   {
//     std::vector<std::string> tmp_fname_v = Zaki::String::Pars(in_dir.Str(), "*") ;
//     std::string tmp_fname ;

//     char tmp_omeg_str[100] ;
//     sprintf(tmp_omeg_str, "%.3e", init_omega_bar*Zaki::Physics::LIGHT_C_KM_S) ;

//     if(tmp_fname_v.size() <2 )
//     {
//       Z_LOG_WARNING("File name patter doesn't match '[]*[]'.") ;
//       tmp_fname = tmp_fname_v[0] + std::string(tmp_omeg_str) ;
//     }
//     else
//     {
//       tmp_fname = tmp_fname_v[0] + std::string(tmp_omeg_str) + tmp_fname_v[1] ;
//     }
    
//     Zaki::File::VecSaver vec_saver(wrk_dir + "/" + tmp_fname) ;

//     // .............................................
//     // Header
//     std::stringstream ss ;
//     char res_header[300] ;

//     sprintf (res_header, "\u2554%-40s Sequence %-3zu out of %-3lu %45s\u2557\n",
//     Zaki::String::Multiply("\u2550", 40).c_str(), (size_t)1, (size_t)1, 
//     Zaki::String::Multiply("\u2550", 45).c_str());
//     ss << "# " << res_header ;

//     sprintf (res_header, "\u2551 %14s %12s %17s %16s %14s %12s %12s %11s\n", "R (km)", "M (Sun)", "\u03B5_c (g/cm^3)", 
//     "\u03C9c_bar (s^-1)", "\u03A9 (s^-1) ", "J (km^2)", "I (km^3)", "\u2551");
//     ss << "# " << res_header ;

//     sprintf (res_header, "\u2551 %14le %14le %14le %14le %14le %14le %14le %7s\n",
//     nstar_ptr->radius[-1],
//     nstar_ptr->mass[-1],
//     nstar_ptr->eps[0], 
//     init_omega_bar*Zaki::Physics::LIGHT_C_KM_S,
//     ang_vel_Omega*Zaki::Physics::LIGHT_C_KM_S,
//     ang_mom_J,
//     mom_inertia,
//     "\u2551");
//     ss << "# " << res_header ;

//     sprintf (res_header, "\u255A%110s\u255D\n", Zaki::String::Multiply("\u2550", 110).c_str());
//     ss << "# " << res_header ;

//     sprintf(res_header, "%-19s\t %-19s\t %-19s\t %-19s\t %-19s\t %-19s", 
//             "r [km]", "M [Sun]", "nu", "omega [1/s]", "omega_bar [1/km]",  "omega'_bar [1/km^2]" ) ;

//     ss << res_header ;

//     // std::string tmp_label = res_header ;
//     // tmp_label.c_str()
//     // .............................................

//     vec_saver.SetHeader(ss.str().c_str()) ;
//     vec_saver.Export1D(omega_results) ;
//   }
//   // ------------------------------------------------------------
  
//   // Testing this for freeing memory
//   // omega_results = std::vector<OmegaPoint>() ;
//   // instead of 
//   omega_results.clear() ;
//   omega_results.shrink_to_fit() ;

//   std::cout <<"\n\t\t ************************"
//             <<" Rot. Solver Finished *************************"<<"\n\n" ;
//   // ------------------------------------------------------------

// }

//--------------------------------------------------------------
// Added on Apr 20, 2022
// Input an initial omega(0) value
// Need to fix this method!!!!
void RotationSolver::Solve_Mixed( const double& in_omega_0,
                            const Zaki::String::Directory& in_dir) 
{
  PROFILE_FUNCTION() ;

  // Resolution of the solver (division of the radial distance)
  // unsigned int solver_res = 5000 ;

  double r_min = mixedstar_ptr->core_region.Min()  ;
  double r_max = 1.1*mixedstar_ptr->mantle_region.Max() ; 

#if R_SOLVER_VERBOSE
  std::cout << "\n\n\t\t ****************************************"
            <<"*******************************"<<" \n" ;
  std::cout <<     "\t\t *                   "
            <<"Rotation Solver Sequence Results"
            <<"                  * \n" ;
  std::cout <<     "\t\t ******************************************"
            <<"*****************************"<<"\n\n" ;
#endif

  init_omega_bar = in_omega_0 ;
  double r = r_min ; 
//  double r1 = r_max ;

  double y[2] ; 
  y[0] = init_omega_bar ;
  y[1] = 0 ;

  // ------------------------------------------------
  // Inside the mixed star:
  gsl_odeiv2_system ode_sys = {CompactStar::RotationSolver::ODE_Mixed, nullptr, 2, this} ;

  gsl_odeiv2_driver *tmp_driver = gsl_odeiv2_driver_alloc_y_new
    (&ode_sys, gsl_odeiv2_step_rk8pd,
      1.e-1, 1.e-10, 1.e-10);
  // ------------------------------------------------
  // Outside the mixed star
  gsl_odeiv2_system ode_sys_out = {CompactStar::RotationSolver::ODE_Mixed_Out, nullptr, 2, this} ;
  
  gsl_odeiv2_driver *tmp_driver_out = gsl_odeiv2_driver_alloc_y_new
  (&ode_sys_out, gsl_odeiv2_step_rk8pd,
   1.e-1, 1.e-10, 1.e-10);
  // ------------------------------------------------
  
//  double min_log_r  = log10(r_min) ;
//  double max_log_r  = log10(r_max) ;
  // -------------------------------------------------
  // Added on Apr 20, 2022
  double r_surface = mixedstar_ptr->mantle_region.Max() ;

  // Added on March 23, 2022 to fix the bug 
  // in J, Omega, I calculations.
//  double log_r_star = log10(r_surface) ;
  
  // Total angular momentum and velocity:
  double ang_mom_J, ang_vel_Omega, mom_inertia ;
  // -------------------------------------------------


//  double step = (max_log_r - min_log_r) / radial_res ;
  double step = (r_max - r_min) / radial_res ;
  
  // One extra point for the exact value at r = R_star.
  omega_results.reserve(radial_res+1) ;
//  bool surface_reached = false ;

  // ---------------------------------------------------------------------
  //                         INNER RADIUS LOOP BEGINS
  // ---------------------------------------------------------------------
  for (double r_i = r_min;  r_i < r_surface; r_i += step)
  {
  //   PROFILE_SCOPE("R-Loop") ;

//    double ri = pow(10, log_r_i) ;

    // This function evolves the driver system d from t to t1.
    // Initially vector y should contain the values of dependent
    // variables at point t. 
    int status = gsl_odeiv2_driver_apply (tmp_driver, &r, r_i, y);

    if (status != GSL_SUCCESS)
    {
#if R_SOLVER_VERBOSE
      printf ("error, return value=%d\n", status);
#endif
      break;
    }

    // Units : { [ km ], [ M_Sun ], - , [ km^-1 ], [ km^-2 ] }
    omega_results_dark.emplace_back(r, 
                  mixedstar_ptr->GetMass_Visible(r) / Zaki::Physics::SUN_M_KM,
                  mixedstar_ptr->GetMass_Dark(r) / Zaki::Physics::SUN_M_KM,
                  mixedstar_ptr->GetNu(r), y[0], y[1]);


    // March 23, 2022
//    if ( !surface_reached && abs(log_r_i - log_r_star) < step )
//    {
//
//      status = gsl_odeiv2_driver_apply (tmp_driver, &r, r_surface, y) ;
//
//      omega_results_dark.emplace_back(r,
//                  mixedstar_ptr->GetMass_Visible(r) / Zaki::Physics::SUN_M_KM,
//                  mixedstar_ptr->GetMass_Dark(r) / Zaki::Physics::SUN_M_KM,
//                  mixedstar_ptr->GetNu(r), y[0], y[1]) ;
//
//      // Total angular momentum and velocity:
//      ang_mom_J = pow(r_surface,4)*y[1]/6. ;
//      ang_vel_Omega = y[0] + r_surface*y[1]/3. ;
//      mom_inertia = ang_mom_J / ang_vel_Omega ;
//
//      surface_reached = true ;
//    }
  }
  // ............ Inner Radius Loop Ends .........................
  gsl_odeiv2_driver_apply (tmp_driver, &r, r_surface, y) ;
    
  omega_results_dark.emplace_back(r,
                                  mixedstar_ptr->GetMass_Visible(r) / Zaki::Physics::SUN_M_KM,
                                  mixedstar_ptr->GetMass_Dark(r) / Zaki::Physics::SUN_M_KM,
                                  mixedstar_ptr->GetNu(r), y[0], y[1]) ;
    
  // Total angular momentum and velocity:
  ang_mom_J = pow(r_surface,4)*y[1]/6. ;
  ang_vel_Omega = y[0] + r_surface*y[1]/3. ;
  mom_inertia = ang_mom_J / ang_vel_Omega ;
    
  // ---------------------------------------------------------------------
  //                         OUTER RADIUS LOOP BEGINS
  // ---------------------------------------------------------------------
  for (double r_i = r_surface+step;  r_i < r_max; r_i += step)
  {
    if ( GSL_SUCCESS != gsl_odeiv2_driver_apply (tmp_driver_out, &r, r_i, y) )
      break ;
    
    // Units : { [ km ], [ M_Sun ], - , [ km^-1 ], [ km^-2 ] }
    omega_results_dark.emplace_back(r,
                                    mixedstar_ptr->GetMass_Visible(r) / Zaki::Physics::SUN_M_KM,
                                    mixedstar_ptr->GetMass_Dark(r) / Zaki::Physics::SUN_M_KM,
                                    mixedstar_ptr->GetNu(r), y[0], y[1]);
  }
  // ............ Outer Radius Loop Ends .........................
  
  // ............................. omega_results_dark Loop Begins .....................................
  // We go back now and update the omega_results_dark with the values for omega(r)
  // We needed ang_vel_Omega for omega(r) since:
  // omega(r) = ang_vel_Omega - bar{omega(r)}
  // The unit for omega(r) is changed from km^{-1} to s^{-1}
  for (size_t i = 0; i < omega_results_dark.size(); i++)
  {
    omega_results_dark[i].omega = (
                                    ang_vel_Omega - omega_results_dark[i].omega_bar
                                  ) * Zaki::Physics::LIGHT_C_KM_S ;
  }
  // .............................. omega_results_dark Loop Ends ......................................

#if R_SOLVER_VERBOSE
  // ...............................................................................................
  // Printing the sequence point data in terminal
  printf ("\u2554%-40s \u03A9-Sequence %-3zu out of %-3lu %43s\u2557\n",
  Zaki::String::Multiply("\u2550", 40).c_str(), (size_t)1, (size_t)1, 
  Zaki::String::Multiply("\u2550", 43).c_str());
  printf ("\u2551 %14s %14s %17s %16s %14s %12s %12s %9s\n", "R (km)", "M (km)", "\u03B5_c (km^-2)", 
  "\u03C9c_bar (s^-1)", "\u03A9 (s^-1) ", "J (km^2)", "I (km^3)", "\u2551");
  printf ("\u2551 %14le %14le %14le %14le %14le %14le %14le %7s\n",
  r_surface,
  mixedstar_ptr->mass[-1],
  mixedstar_ptr->eps[0], 
  init_omega_bar*Zaki::Physics::LIGHT_C_KM_S,
  ang_vel_Omega*Zaki::Physics::LIGHT_C_KM_S,
  ang_mom_J,
  mom_inertia,
  "\u2551");
  printf ("\u255A%110s\u255D\n", Zaki::String::Multiply("\u2550", 110).c_str());
#endif
  // ...............................................................................................

  omega_seq_pts.emplace_back(init_omega_bar*Zaki::Physics::LIGHT_C_KM_S, 
                              mixedstar_ptr->sequence.v.m, 
                              r_surface, ang_mom_J, ang_vel_Omega*Zaki::Physics::LIGHT_C_KM_S );

  gsl_odeiv2_driver_free (tmp_driver) ;

  // ------------------------------------------------------------
  //                      Saving to file
  // ------------------------------------------------------------
  if(false && in_dir.Str() != "")
  {
    std::vector<std::string> tmp_fname_v = Zaki::String::Pars(in_dir.Str(), "*") ;
    std::string tmp_fname ;

    char tmp_omeg_str[100] ;
    snprintf(tmp_omeg_str, sizeof(tmp_omeg_str), 
              "%.3e", init_omega_bar*Zaki::Physics::LIGHT_C_KM_S) ;

    if(tmp_fname_v.size() <2 )
    {
      Z_LOG_WARNING("File name patter doesn't match '[]*[]'.") ;
      tmp_fname = tmp_fname_v[0] + std::string(tmp_omeg_str) ;
    }
    else
    {
      tmp_fname = tmp_fname_v[0] + std::string(tmp_omeg_str) + tmp_fname_v[1] ;
    }
    
    Zaki::File::VecSaver vec_saver(wrk_dir + "/" + tmp_fname) ;

    // .............................................
    // Header
    std::stringstream ss ;
    char res_header[300] ;

    snprintf (res_header, sizeof(res_header), 
              "\u2554%-40s Sequence %-3zu out of %-3lu %45s\u2557\n",
    Zaki::String::Multiply("\u2550", 40).c_str(), (size_t)1, (size_t)1, 
    Zaki::String::Multiply("\u2550", 45).c_str());
    ss << "# " << res_header ;

    snprintf (res_header, sizeof(res_header), 
              "\u2551 %14s %12s %17s %16s %14s %12s %12s %11s\n", "R (km)", "M (Sun)", "\u03B5_c (g/cm^3)", 
    "\u03C9c_bar (s^-1)", "\u03A9 (s^-1) ", "J (km^2)", "I (km^3)", "\u2551");
    ss << "# " << res_header ;

    snprintf (res_header, sizeof(res_header), 
                "\u2551 %14le %14le %14le %14le %14le %14le %14le %7s\n",
    mixedstar_ptr->sequence.v.r,
    mixedstar_ptr->sequence.v.m,
    mixedstar_ptr->sequence.v.ec, 
    init_omega_bar*Zaki::Physics::LIGHT_C_KM_S,
    ang_vel_Omega*Zaki::Physics::LIGHT_C_KM_S,
    ang_mom_J,
    mom_inertia,
    "\u2551");
    ss << "# " << res_header ;

    snprintf (res_header, sizeof(res_header), "\u255A%110s\u255D\n", Zaki::String::Multiply("\u2550", 110).c_str());
    ss << "# " << res_header ;

    snprintf(res_header, sizeof(res_header), "%-19s\t %-19s\t %-19s\t %-19s\t %-19s\t %-19s\t %-19s", 
            "r [km]", "M_V [Sun]", "M_D [Sun]", "nu", "omega [1/s]", 
            "omega_bar [1/km]",  "omega'_bar [1/km^2]" ) ;

    ss << res_header ;

    // std::string tmp_label = res_header ;
    // tmp_label.c_str()
    // .............................................

    vec_saver.SetHeader(ss.str().c_str()) ;
    vec_saver.Export1D(omega_results_dark) ;
  }
  // ------------------------------------------------------------
  
  // Testing this for freeing memory
  // omega_results_dark = std::vector<OmegaPointDark>() ;
  // instead of :
  omega_results_dark.clear() ;
  omega_results_dark.shrink_to_fit() ;

#if R_SOLVER_VERBOSE
  std::cout <<"\n\t\t ************************"
            <<" Rot. Solver Finished *************************"<<"\n\n" ;
#endif
  // ------------------------------------------------------------

}

//--------------------------------------------------------------
// Added on Aug 6, 2020
// Exports the results of solving the rotation equations
void RotationSolver::ExportResults(const Zaki::String::Directory& in_dir) const 
{
  Zaki::File::VecSaver vec_saver_2(wrk_dir + "/" + in_dir);

  char seq_header[200] ;
  snprintf(seq_header, sizeof(seq_header), "%-16s\t %-14s\t %-14s\t %-14s\t %-16s", 
          "omega_bar_c (1/s)", "M",  "R", "J", "Omega (1/s)" ) ;

  vec_saver_2.SetHeader(seq_header) ;
  vec_saver_2.Export1D(omega_seq_pts) ;
}

//--------------------------------------------------------------
/// Returns omega_seq_pts
std::vector<OmegaSeqPoint> RotationSolver::GetOmegaSeq() const 
{
  return omega_seq_pts ;
}

//--------------------------------------------------------------
NStar* RotationSolver::GetNStar() 
{
 return nstar_ptr ; 
}

//--------------------------------------------------------------
MixedStar* RotationSolver::GetMixedStar() 
{
 return mixedstar_ptr ; 
}

//--------------------------------------------------------------
// Resets the containers
void RotationSolver::Reset() 
{

}

//--------------------------------------------------------------
// Added on June 3, 2022
// Evaluates the moment of inertia for the neutron star
void RotationSolver::FindNMomInertia() 
{
  double r_min = nstar_ptr->GetRadius()->operator[](0)  ;
  double r_surface = nstar_ptr->GetRadius()->operator[](-1) ;

  init_omega_bar = 5e-3 ;
  double r = r_min ; 

  double y[2] ; 
  y[0] = init_omega_bar ;
  y[1] = 0 ;

  double ang_mom_J, ang_vel_Omega, mom_inertia ;


gsl_odeiv2_system fast_ode_sys = {CompactStar::RotationSolver::ODE_N_Fast, nullptr, 2, this} ;

  gsl_odeiv2_driver *fast_driver = gsl_odeiv2_driver_alloc_y_new
    (&fast_ode_sys, gsl_odeiv2_step_rk8pd,
      1.e-1, 1.e-10, 1.e-10);

  // Radius loop inside the core
  for (size_t i = 0; i < nstar_ptr->ds[0].Size() ; i++)
  {
    fast_p = nstar_ptr->ds[3][i] ;
    fast_e = nstar_ptr->ds[4][i] ;
    fast_m = nstar_ptr->ds[1][i] ;

    if ( GSL_SUCCESS !=
        gsl_odeiv2_driver_apply (fast_driver, &r, nstar_ptr->ds[0][i], y) )
      break ;
  }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++

  // Total angular momentum and velocity:
  ang_mom_J = pow(r_surface,4)*y[1]/6. ;
  ang_vel_Omega = y[0] + r_surface*y[1]/3. ;
  mom_inertia = ang_mom_J / ang_vel_Omega ;
    
  nstar_ptr->MomI = mom_inertia ;

  gsl_odeiv2_driver_free (fast_driver) ;

  // nstar_ptr->MomI = 0 ;
}

//--------------------------------------------------------------
// Added on Apr 22, 2022
void RotationSolver::FindMixedMomInertia() 
{
  // Resolution of the solver (division of the radial distance)
  // unsigned int solver_res = 5000 ;

  double r_min = mixedstar_ptr->core_region.Min()  ;
//  double r_mantle = mixedstar_ptr->mantle_region.Min() ;
//  double r_max = mixedstar_ptr->mantle_region.Max() ;

  init_omega_bar = 5e-3 ;
  double r = r_min ; 

  double y[2] ; 
  y[0] = init_omega_bar ;
  y[1] = 0 ;

  // gsl_odeiv2_system ode_sys = {CompactStar::RotationSolver::ODE_Mixed, nullptr, 2, this} ;

  // gsl_odeiv2_driver *tmp_driver = gsl_odeiv2_driver_alloc_y_new
  //   (&ode_sys, gsl_odeiv2_step_rk8pd,
  //     1.e-1, 1.e-10, 1.e-10);
  

  // -------------------------------------------------
  // Added on Apr 20, 2022
  double r_surface = mixedstar_ptr->mantle_region.Max() ;

  // Added on March 23, 2022 to fix the bug 
  // in J, Omega, I calculations.
//  double log_r_star = log10(r_surface) ;
  
  // Total angular momentum and velocity:
  double ang_mom_J, ang_vel_Omega, mom_inertia ;
  // -------------------------------------------------

  // double step = (r_surface - r_min) / radial_res ;
  
  // One extra point for the exact value at r = R_star.
  // bool surface_reached = false ;

  // ---------------------------------------------------------------------
  //                         RADIUS LOOP BEGINS
  // ---------------------------------------------------------------------
  // for (double r_i = r_min;  r_i < r_surface; r_i += step)
  // {
  //   if ( GSL_SUCCESS !=
  //       gsl_odeiv2_driver_apply (tmp_driver, &r, r_i, y) )
  //     break ;
  // }
  // ............  Radius Loop Ends .........................

  // gsl_odeiv2_driver_apply (tmp_driver, &r, r_surface, y) ;
  // ++++++++++++++++++++++++++++++++++++++++++++++++++   

  // ++++++++++++++++++++++++++++++++++++++++++++++++++
  // New method:
  gsl_odeiv2_system fast_ode_sys = {CompactStar::RotationSolver::ODE_Mixed_Fast, nullptr, 2, this} ;

  gsl_odeiv2_driver *fast_driver = gsl_odeiv2_driver_alloc_y_new
    (&fast_ode_sys, gsl_odeiv2_step_rk8pd,
      1.e-1, 1.e-10, 1.e-10);

  if(mixedstar_ptr->dark_core)
  {
    // Loop inside the core
    for (size_t i = 0; i < mixedstar_ptr->ds_dar[0].Size() ; i++)
    {
      fast_p_v = mixedstar_ptr->ds_vis[3][i] ;
      fast_p_d = mixedstar_ptr->ds_dar[3][i] ;
      fast_e_v = mixedstar_ptr->ds_vis[4][i] ;
      fast_e_d = mixedstar_ptr->ds_dar[4][i] ;
      fast_m_tot = mixedstar_ptr->ds_vis[1][i] + mixedstar_ptr->ds_dar[1][i] ;

      if ( GSL_SUCCESS !=
          gsl_odeiv2_driver_apply (fast_driver, &r, mixedstar_ptr->ds_dar[0][i], y) )
        break ;
    }
    // Loop inside the mantle
    for (size_t i = mixedstar_ptr->ds_dar[0].Size() ; i < mixedstar_ptr->ds_vis[0].Size() ; i++)
    {
      fast_p_v = mixedstar_ptr->ds_vis[3][i] ;
      fast_p_d = 0 ;
      fast_e_v = mixedstar_ptr->ds_vis[4][i] ;
      fast_e_d = 0 ;
      fast_m_tot = mixedstar_ptr->ds_vis[1][i] ;

      if ( GSL_SUCCESS !=
          gsl_odeiv2_driver_apply (fast_driver, &r, mixedstar_ptr->ds_vis[0][i], y) )
        break ;
    }
  }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++
  else // We have a dark mantle
  {
    // Loop inside the core
    for (size_t i = 0; i < mixedstar_ptr->ds_vis[0].Size() ; i++)
    {
      fast_p_v = mixedstar_ptr->ds_vis[3][i] ;
      fast_p_d = mixedstar_ptr->ds_dar[3][i] ;
      fast_e_v = mixedstar_ptr->ds_vis[4][i] ;
      fast_e_d = mixedstar_ptr->ds_dar[4][i] ;
      fast_m_tot = mixedstar_ptr->ds_vis[1][i] + mixedstar_ptr->ds_dar[1][i] ;

      if ( GSL_SUCCESS !=
          gsl_odeiv2_driver_apply (fast_driver, &r, mixedstar_ptr->ds_vis[0][i], y) )
        break ;
    }
    // Loop inside the mantle
    for (size_t i = mixedstar_ptr->ds_vis[0].Size() ; i < mixedstar_ptr->ds_dar[0].Size() ; i++)
    {
      fast_p_v = 0 ;
      fast_p_d = mixedstar_ptr->ds_dar[3][i] ;
      fast_e_v = 0 ; 
      fast_e_d = mixedstar_ptr->ds_dar[4][i] ;
      fast_m_tot = mixedstar_ptr->ds_dar[1][i] ;

      if ( GSL_SUCCESS !=
          gsl_odeiv2_driver_apply (fast_driver, &r, mixedstar_ptr->ds_dar[0][i], y) )
        break ;
    }
  }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++

  // Total angular momentum and velocity:
  ang_mom_J = pow(r_surface,4)*y[1]/6. ;
  ang_vel_Omega = y[0] + r_surface*y[1]/3. ;
  mom_inertia = ang_mom_J / ang_vel_Omega ;
    
  mixedstar_ptr->MomI = mom_inertia ;

  gsl_odeiv2_driver_free (fast_driver) ;
 
}
// ------------------------------------------------------------

//==============================================================
