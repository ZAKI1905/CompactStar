/*
  NStar class
*/

// #include <gsl/gsl_math.h>

#include <Zaki/Util/Instrumentor.hpp>
#include <Zaki/Math/GSLFuncWrapper.hpp>
#include <Zaki/Physics/Constants.hpp>

#include <gsl/gsl_integration.h>

#include "CompactStar/Core/NStar.hpp"
#include "CompactStar/Core/TOVSolver.hpp"
#include "CompactStar/Core/RotationSolver.hpp"


//==============================================================
//                        NStar class
//==============================================================
// Default Constructor
// Be sure to run SurfaceIsReached if using this constructor!
CompactStar::NStar::NStar() : Prog("NStar")
{
  // radius.label  = "r"   ;
  // mass.label    = "m"   ;
  // rho.label     = "rho" ;
  // eps.label     = "eps" ;
  // press.label   = "p" ;
  // nu_der.label  = "nu'" ;
  // nu.label      = "nu" ;

  rot_solver.AttachNStar(this) ;
}

//--------------------------------------------------------------
// Constructor from TOV Solutions
CompactStar::NStar::NStar(const std::vector<CompactStar::TOVPoint>& in_tov) 
: Prog("NStar", true),
  ds(7+in_tov[0].rho_i.size(), in_tov.size()), 
  B_integrand(2, in_tov.size())
{ 

  B_integrand[0].label = "r(km)" ;       
  B_integrand[1].label = "B_v" ;

  rot_solver.AttachNStar(this) ;

  ds[m_idx].label       = "m" ;
  ds[r_idx].label       = "r" ;
  ds[rho_idx].label     = "rho" ;
  ds[eps_idx].label     = "eps" ;
  ds[pre_idx].label     = "p" ;
  ds[nu_der_idx].label  = "nu'" ;
  ds[nu_idx].label      = "nu" ;

  for (size_t j = 0; j < in_tov[0].rho_i.size(); j++)
  {
    rho_i_idx.emplace_back(6 + j) ;
    ds[ rho_i_idx[j] ].label = "X" ; // !!!
  }

  // ...............................................
  for (auto &&i : in_tov)
  {
    ds[r_idx].vals.emplace_back(i.r) ; // in km
    ds[m_idx].vals.emplace_back(Zaki::Physics::SUN_M_KM*i.m) ; // in km
    ds[rho_idx].vals.emplace_back(i.rho) ; // in fm^{-3}
    ds[eps_idx].vals.emplace_back(i.e * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_G_CM3) ; // in km^{-2}
    ds[pre_idx].vals.emplace_back(i.p * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_Dyn_CM2 ) ; // in km^{-2}
    ds[nu_der_idx].vals.emplace_back(i.nu_der * 1e+5) ; // convert 1/cm to 1/km

    for (size_t j = 0; j < i.rho_i.size(); j++)
    {
      ds[ rho_i_idx[j] ].vals.emplace_back( i.rho_i[j] ) ;
    }
  }
  // ...............................................

  ds.Interpolate(r_idx, {m_idx, nu_der_idx, rho_idx, eps_idx, pre_idx} ) ;

  EvaluateNu() ;

  // ..................................................................  
  //            B Integrand
  // ..................................................................  
  B_integrand[0] = ds[r_idx] ;
  B_integrand[1] = ds[r_idx].pow(2) ;
  B_integrand[1] *= 4*M_PI * ds[rho_idx] ; 

  B_integrand[1] /= (1. - 2.* ds[m_idx] / ds[r_idx]).sqrt() ;

  // converting fm^{-3} to km^{-3}
  B_integrand[1] *= pow(10, 54) ;

  B_integrand.Interpolate(0, 1) ;
  // ..................................................................  


  sequence.ec = ds[eps_idx][0] * Zaki::Physics::INV_FM4_2_G_CM3 /
                                      Zaki::Physics::INV_FM4_2_INV_KM2  ;
  sequence.m  = ds[m_idx][-1] / Zaki::Physics::SUN_M_KM ;
  sequence.r  = ds[r_idx][-1] ;
  sequence.pc = ds[pre_idx][0] * Zaki::Physics::INV_FM4_2_Dyn_CM2 /
                                      Zaki::Physics::INV_FM4_2_INV_KM2 ;
  // sequence.b  = Find_BaryonNum_Visible() ;
  sequence.b  = B_integrand.Integrate(1, {ds[r_idx][0], ds[r_idx][-1]}) ;

  sequence.I  = Find_MomInertia() ;
  // ..................................................................  

  // accel  = gsl_interp_accel_alloc ();

  // mass_r_spline = gsl_spline_alloc (gsl_interp_cspline, in_tov.size()) ;
  // rho_r_spline  = gsl_spline_alloc (gsl_interp_cspline, in_tov.size()) ;
  // eps_r_spline  = gsl_spline_alloc (gsl_interp_cspline, in_tov.size()) ;
  // press_r_spline = gsl_spline_alloc (gsl_interp_cspline, in_tov.size()) ;
  // nu_der_r_spline = gsl_spline_alloc (gsl_interp_cspline, in_tov.size()) ;
  // nu_r_spline = gsl_spline_alloc (gsl_interp_cspline, in_tov.size()) ;

  // radius.label  = "r"   ;
  // mass.label    = "m"   ;
  // rho.label     = "rho" ;
  // eps.label     = "eps" ;
  // press.label   = "p" ;
  // nu_der.label  = "nu'" ;
  // nu.label      = "nu" ;


  // radius.Reserve(in_tov.size()) ;
  // mass.Reserve(in_tov.size()) ;
  // rho.Reserve(in_tov.size()) ;
  // eps.Reserve(in_tov.size()) ;
  // press.Reserve(in_tov.size()) ;
  // nu_der.Reserve(in_tov.size()) ;
  // nu.Reserve(in_tov.size()) ;

  // for (auto &&i : in_tov)
  // {
  //   radius.vals.emplace_back(i.r) ; // in km
  //   mass.vals.emplace_back(Zaki::Physics::SUN_M_KM*i.m) ; // in km
  //   rho.vals.emplace_back(i.rho) ; // in fm^{-3}
  //   eps.vals.emplace_back(i.e * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_G_CM3) ; // in km^{-2}
  //   press.vals.emplace_back(i.p * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_Dyn_CM2 ) ; // in km^{-2}
  //   nu_der.vals.emplace_back(i.nu_der * 1e+5) ; // convert 1/cm to 1/km
  // }

  // gsl_spline_init (mass_r_spline, &radius.vals[0], &mass.vals[0], in_tov.size()) ;
  // gsl_spline_init (rho_r_spline, &radius.vals[0], &rho.vals[0], in_tov.size()) ;
  // gsl_spline_init (eps_r_spline, &radius.vals[0], &eps.vals[0], in_tov.size()) ;
  // gsl_spline_init (press_r_spline, &radius.vals[0], &press.vals[0], in_tov.size()) ;
  // gsl_spline_init (nu_der_r_spline, &radius.vals[0], &nu_der.vals[0], in_tov.size()) ;
}

//--------------------------------------------------------------
// Initializes the dataset
void CompactStar::NStar::Init(const TOVSolver* in_tov_solver) 
{
// Resizes the columns to '7', and
  // reserves 'radial_res' space for each column
  ds.Reserve(7+in_tov_solver->eos_tab.rho_i.size(), 
                    in_tov_solver->radial_res) ;

  ds[m_idx].label       = "m(km)" ;
  ds[r_idx].label       = "r(km)" ;
  ds[rho_idx].label     = "rho(fm^-3)" ;
  ds[eps_idx].label     = "eps(km^-2)" ;
  ds[pre_idx].label     = "p(km^-2)" ;
  ds[nu_der_idx].label  = "nu'(km^-1)" ;
  ds[nu_idx].label      = "nu" ;

  for (size_t i = 0; i < in_tov_solver->eos_tab.rho_i.size() ; i++)
  {
    rho_i_idx.emplace_back( i + 7 ) ;
    ds[ rho_i_idx[i] ].label =
                            in_tov_solver->eos_tab.extra_labels[i] ;
  }

  B_integrand.Reserve(2, in_tov_solver->radial_res) ;

  B_integrand[0].label = "r(km)" ;       
  B_integrand[1].label = "B" ;
  
}

//--------------------------------------------------------------
// This has to be run so the class
// knows when to initialize all the splines
void CompactStar::NStar::SurfaceIsReached() 
{
  PROFILE_FUNCTION() ;

  ds.Interpolate(r_idx, {m_idx, nu_der_idx, rho_idx, eps_idx, pre_idx} ) ;
  
  EvaluateNu() ;

  // ..................................................................  
  //                                B Integrand
  // ..................................................................  
  B_integrand[0] = ds[r_idx] ;
  B_integrand[1] = ds[r_idx].pow(2) ;
  B_integrand[1] *= 4*M_PI * ds[rho_idx] ; 

  B_integrand[1] /= (1. - 2.* ds[m_idx] / ds[r_idx]).sqrt() ;

  // converting fm^{-3} to km^{-3}
  B_integrand[1] *= pow(10, 54) ;

  B_integrand.Interpolate(0, 1) ;
  // ..................................................................  

  // Units must be converted from km
  sequence.ec = ds[eps_idx][0] * Zaki::Physics::INV_FM4_2_G_CM3 /
                                      Zaki::Physics::INV_FM4_2_INV_KM2  ;
  sequence.m  = ds[m_idx][-1] / Zaki::Physics::SUN_M_KM ;
  sequence.r  = ds[r_idx][-1] ;
  sequence.pc = ds[pre_idx][0] * Zaki::Physics::INV_FM4_2_Dyn_CM2 /
                                      Zaki::Physics::INV_FM4_2_INV_KM2 ;
  // sequence.b  = Find_BaryonNum_Visible() ;
  sequence.b  = B_integrand.Integrate(1, {ds[r_idx][0], ds[r_idx][-1]}) ;

  sequence.I  = Find_MomInertia() ;
  // accel  = gsl_interp_accel_alloc ();

  // mass_r_spline = gsl_spline_alloc (gsl_interp_cspline, radius.Size()) ;
  // rho_r_spline  = gsl_spline_alloc (gsl_interp_cspline, radius.Size()) ;
  // eps_r_spline  = gsl_spline_alloc (gsl_interp_cspline, radius.Size()) ;
  // press_r_spline = gsl_spline_alloc (gsl_interp_cspline, radius.Size()) ;
  // nu_der_r_spline = gsl_spline_alloc (gsl_interp_cspline, radius.Size()) ;
  // nu_r_spline = gsl_spline_alloc (gsl_interp_cspline, radius.Size()) ;

 

  // gsl_spline_init (mass_r_spline, &radius.vals[0], &mass.vals[0], radius.Size()) ;
  // gsl_spline_init (rho_r_spline, &radius.vals[0], &rho.vals[0], radius.Size()) ;
  // gsl_spline_init (eps_r_spline, &radius.vals[0], &eps.vals[0], radius.Size()) ;
  // gsl_spline_init (press_r_spline, &radius.vals[0], &press.vals[0], radius.Size()) ;
  // gsl_spline_init (nu_der_r_spline, &radius.vals[0], &nu_der.vals[0], radius.Size()) ;
}
//--------------------------------------------------------------
// Reserving the needed space for all datacolumns
// void CompactStar::NStar::Reserve(const size_t& space_size) 
// {
//   radius.Reserve(space_size) ;
//   mass.Reserve(space_size) ;
//   rho.Reserve(space_size) ;
//   eps.Reserve(space_size) ;
//   press.Reserve(space_size) ;
//   nu_der.Reserve(space_size) ;
//   nu.Reserve(space_size) ;
// }
//--------------------------------------------------------------
// Appends tov points to the NStar
void CompactStar::NStar::Append(const TOVPoint& in_tov) 
{
  // radius.vals.emplace_back(in_tov.r) ; // in km
  // mass.vals.emplace_back(Zaki::Physics::SUN_M_KM*in_tov.m) ; // in km
  // rho.vals.emplace_back(in_tov.rho) ; // in fm^{-3}
  // eps.vals.emplace_back(in_tov.e * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_G_CM3) ; // in km^{-2}
  // press.vals.emplace_back(in_tov.p * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_Dyn_CM2 ) ; // in km^{-2}
  // nu_der.vals.emplace_back(in_tov.nu_der * 1e+5) ; // convert 1/cm to 1/km

  ds[r_idx].vals.emplace_back(in_tov.r) ; // in km
  ds[m_idx].vals.emplace_back(Zaki::Physics::SUN_M_KM*in_tov.m) ; // in km
  ds[rho_idx].vals.emplace_back(in_tov.rho) ; // in fm^{-3}
  ds[eps_idx].vals.emplace_back(in_tov.e * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_G_CM3) ; // in km^{-2}
  ds[pre_idx].vals.emplace_back(in_tov.p * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_Dyn_CM2 ) ; // in km^{-2}
  ds[nu_der_idx].vals.emplace_back(in_tov.nu_der * 1e+5) ; // convert 1/cm to 1/km

  for (size_t i = 0; i < in_tov.rho_i.size(); i++)
  {
    // in fm^{-3}
    ds[ rho_i_idx[i] ].vals.emplace_back( in_tov.rho_i[i] ) ;
  }
}

//--------------------------------------------------------------
// Sets the work directory for the member objects
CompactStar::Prog* CompactStar::NStar::SetMemWrkDir(
  const Zaki::String::Directory& in_dir) 
{
  ds.SetWrkDir( in_dir ) ;

  return this ;
}

//--------------------------------------------------------------

// void CompactStar::NStar::SetSharedPtr(std::shared_ptr<CompactStar::NStar> in_ptr) 
// {
//   ns_shared_ptr = in_ptr ;
// }

// //--------------------------------------------------------------

// std::shared_ptr<CompactStar::NStar> CompactStar::NStar::GetSharedPtr() 
// {
//   return ns_shared_ptr ;
// }

//--------------------------------------------------------------
/// Initializing the class
// void CompactStar::NStar::Init(const std::vector<TOVPoint>& in_tov) 
// {

//   mass_r_spline = gsl_spline_alloc (gsl_interp_cspline, in_tov.size()) ;
//   rho_r_spline  = gsl_spline_alloc (gsl_interp_cspline, in_tov.size()) ;
//   eps_r_spline  = gsl_spline_alloc (gsl_interp_cspline, in_tov.size()) ;
//   press_r_spline = gsl_spline_alloc (gsl_interp_cspline, in_tov.size()) ;
//   nu_der_r_spline = gsl_spline_alloc (gsl_interp_cspline, in_tov.size()) ;
//   nu_r_spline = gsl_spline_alloc (gsl_interp_cspline, in_tov.size()) ;

//   if (radius.Size() != 0 )
//     Reset() ;

//   radius.Reserve(in_tov.size()) ;
//   mass.Reserve(in_tov.size()) ;
//   rho.Reserve(in_tov.size()) ;
//   eps.Reserve(in_tov.size()) ;
//   press.Reserve(in_tov.size()) ;
//   nu_der.Reserve(in_tov.size()) ;
//   nu.Reserve(in_tov.size()) ;

//   for (auto &&i : in_tov)
//   {
//     radius.vals.emplace_back(i.r) ; // in km
//     mass.vals.emplace_back(Zaki::Physics::SUN_M_KM*i.m) ; // in km
//     rho.vals.emplace_back(i.rho) ; // in fm^{-3}
//     eps.vals.emplace_back(i.e * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_G_CM3) ; // in km^{-2}
//     press.vals.emplace_back(i.p * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_Dyn_CM2 ) ; // in km^{-2}
//     nu_der.vals.emplace_back(i.nu_der * 1e+5) ; // convert 1/cm to 1/km
//   }

//   gsl_spline_init (mass_r_spline, &radius.vals[0], &mass.vals[0], in_tov.size()) ;
//   gsl_spline_init (rho_r_spline, &radius.vals[0], &rho.vals[0], in_tov.size()) ;
//   gsl_spline_init (eps_r_spline, &radius.vals[0], &eps.vals[0], in_tov.size()) ;
//   gsl_spline_init (press_r_spline, &radius.vals[0], &press.vals[0], in_tov.size()) ;
//   gsl_spline_init (nu_der_r_spline, &radius.vals[0], &nu_der.vals[0], in_tov.size()) ;
// }

//--------------------------------------------------------------
/// Similar to the destructor
void CompactStar::NStar::Reset() 
{
  ds.ClearRows() ;
  B_integrand.ClearRows() ;
  sequence.Reset() ;
  // if(mass_r_spline)
  //   gsl_spline_free (mass_r_spline);
  
  // if(rho_r_spline)
  //   gsl_spline_free (rho_r_spline);
  
  // if(eps_r_spline)
  //   gsl_spline_free (eps_r_spline);

  // if(press_r_spline)
  //   gsl_spline_free (press_r_spline);
    
  // if(nu_der_r_spline)
  //   gsl_spline_free (nu_der_r_spline);
    
  // if(nu_r_spline)
  //   gsl_spline_free (nu_r_spline);

  // if(accel)
  //   gsl_interp_accel_free (accel);

  // radius.vals.clear() ;
  // mass.vals.clear() ;
  // rho.vals.clear() ;
  // eps.vals.clear() ;
  // press.vals.clear() ;
  // nu_der.vals.clear() ;
  // nu.vals.clear() ;
}

//--------------------------------------------------------------
/// Constructor from file
// CompactStar::NStar::NStar(const Zaki::String::Directory& in_tov_file) 
// : Prog("NStar")
// {
//   accel  = gsl_interp_accel_alloc ();
// }

//--------------------------------------------------------------
/// Destructor
CompactStar::NStar::~NStar() 
{
  // if(mass_r_spline)
  //   gsl_spline_free (mass_r_spline);
  
  // if(rho_r_spline)
  //   gsl_spline_free (rho_r_spline);
  
  // if(eps_r_spline)
  //   gsl_spline_free (eps_r_spline);

  // if(press_r_spline)
  //   gsl_spline_free (press_r_spline);
    
  // if(nu_der_r_spline)
  //   gsl_spline_free (nu_der_r_spline);

  // if(nu_r_spline)
  //   gsl_spline_free (nu_r_spline);

  // if(accel)
  //   gsl_interp_accel_free (accel);
}

//--------------------------------------------------------------
/// Returns the nu_der value given the radius input
// r is in km ! ----> !!!
// double CompactStar::NStar::GetNuDerSpline(const double& in_r) 
// {
//   return gsl_spline_eval(nu_der_r_spline, in_r, accel);
// }

//--------------------------------------------------------------
/// Evaluate the metric function
void CompactStar::NStar::EvaluateNu() 
{
  PROFILE_FUNCTION() ;

  // -----------------------------------
  // Integrate to find nu(r) :
  // -----------------------------------
  ds[nu_idx].Resize( ds[r_idx].Size() ) ;

  // Boundary condition
  double nu_at_R = 0.5*log(
                  1 - 2 * ds[m_idx][-1] / ds[r_idx].Max()
                          ) ;

  // if (dark_core)
  // {
    // ds.AddColumn("nu") ;

  // Calculate the surface term first to find out delta_nu_r 
  ds[nu_idx][-1] = ds.Integrate(nu_der_idx, { ds[r_idx][0], ds[r_idx][-1] } ) ;
  double delta_nu_r = ds[nu_idx][-1] - nu_at_R ;

  ds[nu_idx][-1] = nu_at_R ;

  // We don't need to calculate the surface term anymore,
  //  therefore we have 'i < ds[r_idx].Size()-1'
  for(size_t i = 0 ; i < ds[r_idx].Size()-1 ; ++i)
  {
    ds[nu_idx][i] = ds.Integrate(nu_der_idx, { ds[r_idx][0], ds[r_idx][i] } ) 
                        - delta_nu_r ;
  }

  ds.Interpolate(0, nu_idx) ;

  // // -----------------------------------
  // // Integrate to find nu(r) :
  // // -----------------------------------

  // gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000*radius.Size());
  // double err;

  // Zaki::Math::GSLFuncWrapper<NStar, double (NStar::*)( const double& )> 
  //   Fp(this, &NStar::GetNuDerSpline);     

  // gsl_function F = static_cast<gsl_function> (Fp) ; 

  // double tmp_nu ;

  // for(size_t i = 0 ; i < radius.Size() ; ++i)
  // {
  //   gsl_integration_qag(&F, radius[0], radius[i], 1e-9, 1e-9, 1000, 1, w, &tmp_nu, &err);
  //   nu.vals.push_back(tmp_nu) ;
  //   // tmp_nu = 0 ;
  // }
  
  // // std::cout << "tmp_nu: " << tmp_nu << "\n" ;
  // gsl_integration_workspace_free(w) ;
  
  // // -------------------------------------------
  // // Saving the results for [ r, M(r), nu(r) ] 
  // // -------------------------------------------
  
  // // Matching the boundary condition for nu(r)
  // double nu_at_R = 0.5*log(1 - 2*mass[-1] / radius[-1]) ;

  // double delta_nu_r = nu[-1] - nu_at_R ;

  // // std::cout << "1 - 2*mass[-1] / radius[-1]: " << 1 - 2*mass[-1] / radius[-1] << "\n" ;
  // // std::cout << "log(1 - 2*mass[-1] / radius[-1]): " << log(1 - 2*mass[-1] / radius[-1]) << "\n" ;
  // // std::cout << "nu[-1]: " << nu[-1] << "\n" ;
  // // std::cout << "nu_at_R: " << nu_at_R << "\n" ;

  // // std::cout << "delta_nu_r: " << delta_nu_r << "\n" ;

  // nu = nu - delta_nu_r ;

  // // std::cout << "nu[-1]: " << nu[-1] << "\n" ;

  // gsl_spline_init (nu_r_spline, &radius.vals[0], &nu.vals[0], nu.Size()) ;
}

//--------------------------------------------------------------
/// Metric function as a function of radius (in km)
double CompactStar::NStar::GetNu(const double& in_r) const 
{
  // return gsl_spline_eval(nu_r_spline, in_r, accel) ;
  return ds.Evaluate(nu_idx, in_r) ;
}

//--------------------------------------------------------------
// Metric function as a function of radius (in km)
Zaki::Vector::DataColumn* CompactStar::NStar::GetNu() 
{
  return &ds[nu_idx] ;
}

//--------------------------------------------------------------
/// Mass (in km) as a function of radius
double CompactStar::NStar::GetMass(const double& in_r) const 
{
  // if (in_r < radius[0] || in_r > radius[-1])
  // {
  //   Z_LOG_ERROR("Input radius is out of range [ " 
  //               + std::to_string(radius[0]) + ", "
  //               + std::to_string(radius[-1]) +" ]") ;

  //   return -1 ;
  // }
  
  // return gsl_spline_eval(mass_r_spline, in_r, accel) ;
  if (in_r < 0)
  {
    Z_LOG_ERROR("Input radius is must be positive!") ;

    return -1 ;
  }
  
  if (in_r < ds[r_idx][0])
    return 0 ;

  if( in_r > ds[r_idx][-1])
    return ds[m_idx][-1] ;

  // return gsl_spline_eval(mass_r_spline, in_r, visi_r_accel) ;
  return ds.Evaluate(m_idx, in_r) ;
}

//--------------------------------------------------------------
// Mass (in km) as a function of radius
Zaki::Vector::DataColumn* 
CompactStar::NStar::GetMass()  
{
  return &ds[m_idx] ;
}

//--------------------------------------------------------------
/// Returns the radius dataset
Zaki::Vector::DataColumn* CompactStar::NStar::GetRadius() 
{
  return &ds[r_idx] ;
}

//--------------------------------------------------------------
/// Baryon number density (fm^{-3}) as a function of radius
double CompactStar::NStar::GetRho(const double& in_r) const 
{
  // if (in_r < radius[0] || in_r > radius[-1])
  // {
  //   Z_LOG_ERROR("Input radius is out of range [ " 
  //               + std::to_string(radius[0]) + ", "
  //               + std::to_string(radius[-1]) +" ]") ;

  //   return -1 ;
  // }
  
  // return gsl_spline_eval(rho_r_spline, in_r, accel) ;

  if (in_r < 0)
  {
    Z_LOG_ERROR("Input radius is must be positive!") ;

    return -1 ;
  }

  if (in_r < ds[r_idx][0] || in_r > ds[r_idx][-1])
  {
    return 0 ;
  }
  
  // return gsl_spline_eval(rho_r_spline, in_r, visi_r_accel) ;
  return ds.Evaluate(rho_idx, in_r) ;
}

//--------------------------------------------------------------
// Baryon number density (fm^{-3}) as a function of radius
Zaki::Vector::DataColumn* 
CompactStar::NStar::GetRho()  
{
  return &ds[rho_idx] ;
}

//--------------------------------------------------------------
// Visible baryon number density (fm^{-3}) as a function of radius
// for a specific species labeled as (in_label)
Zaki::Vector::DataColumn* 
CompactStar::NStar::GetRho_i(const std::string& in_label)  
{
  
  for (auto &&c : rho_i_idx)
  {
    if ( ds[ c ].label ==  in_label )
    {
      return  &ds[ c ] ;
    }
    
  }

  Z_LOG_ERROR("Species density with label '"+ 
                in_label +"' was not found. Returning nullptr.") ;
                
  return nullptr ;
}

//--------------------------------------------------------------
/// Energy density (in km^{-2}) as a function of radius
double CompactStar::NStar::GetEps(const double& in_r) const 
{
  // if (in_r < radius[0] || in_r > radius[-1])
  // {
  //   Z_LOG_ERROR("Input radius is out of range [ " 
  //               + std::to_string(radius[0]) + ", "
  //               + std::to_string(radius[-1]) +" ]") ;

  //   return -1 ;
  // }
  
  // return gsl_spline_eval(eps_r_spline, in_r, accel) ;
  if (in_r < 0)
  {
    Z_LOG_ERROR("Input radius is must be positive!") ;

    return -1 ;
  }

  if (in_r < ds[r_idx][0] || in_r > ds[r_idx][-1])
  {
    return 0 ;
  }
  
  return ds.Evaluate(eps_idx, in_r) ;
}

//--------------------------------------------------------------
// Energy density (in km^{-2}) as a function of radius
Zaki::Vector::DataColumn* 
CompactStar::NStar::GetEps()  
{
  return &ds[eps_idx] ;
}

//--------------------------------------------------------------
/// Pressure (in km^{-2}) as a function of radius
double CompactStar::NStar::GetPress(const double& in_r) const 
{
  // if (in_r < radius[0] || in_r > radius[-1])
  // {
  //   Z_LOG_ERROR("Input radius is out of range [ " 
  //               + std::to_string(radius[0]) + ", "
  //               + std::to_string(radius[-1]) +" ]") ;

  //   return -1 ;
  // }
  
  // return gsl_spline_eval(press_r_spline, in_r, accel) ;
  if (in_r < 0)
  {
    Z_LOG_ERROR("Input radius is must be positive!") ;

    return -1 ;
  }

  if (in_r < ds[r_idx][0] || in_r > ds[r_idx][-1])
  {
    return 0 ;
  }

  // return gsl_spline_eval(press_r_spline, in_r, visi_r_accel) ;
  return ds.Evaluate(pre_idx, in_r) ;
}

//--------------------------------------------------------------
// Pressure (in km^{-2}) as a function of radius
Zaki::Vector::DataColumn* 
CompactStar::NStar::GetPress()  
{
  return &ds[pre_idx] ;
}

//--------------------------------------------------------------
// Returns 'sequence'
CompactStar::SeqPoint CompactStar::NStar::GetSequence() const 
{
  return sequence ;
}

//--------------------------------------------------------------
// Baryon number inegrand
double CompactStar::NStar::BaryonNumIntegrand(double in_r) 
{
  double out = in_r * in_r ;
        out *= 4*M_PI*GetRho(in_r) ; 

        // // converting fm^{-3} to km^{-3}
        // out *= pow(10, 54) ;

        out /= sqrt(1 - 2*GetMass(in_r) / in_r) ;

  return out ;
}
//--------------------------------------------------------------
/// Total baryon number inegrand
// double CompactStar::NStar::BaryonNumIntegrand(double in_r) 
// {
//   double out = in_r * in_r ;
//         out *= 4*M_PI*GetRho(in_r) ; 

//         // converting fm^{-3} to km^{-3}
//         out *= pow(10, 54) ;

//         out /= sqrt(1 - 2*GetMass(in_r) / in_r) ;

//   return out ;
// }

//--------------------------------------------------------------
/// Total baryon number as a function of radius
// double CompactStar::NStar::GetBaryonNum(const double& in_r) 
// {
//   PROFILE_FUNCTION() ;
//   double result, err ;

//   gsl_integration_workspace *w = 
//    gsl_integration_workspace_alloc(2500);

//   Zaki::Math::GSLFuncWrapper<NStar, double (NStar::*)(double)> 
//       Fp(this, &NStar::BaryonNumIntegrand);     

//   gsl_function F = static_cast<gsl_function> (Fp) ; 

//   gsl_integration_qag(&F, radius[0], in_r, 1e-9, 1e-9, 2500, 1, w, &result, &err);

//   return result ;
// }

//--------------------------------------------------------------
/// Total baryon number
// double CompactStar::NStar::GetBaryonNum() 
// {
//   return GetBaryonNum(radius[-1]) ;
// }

//--------------------------------------------------------------
/// Moment of inertia (in km^3) inegrand
// This is wrong and missing a factor of 
//  [ 1 - \omega(r) / \Omega ]
// double CompactStar::NStar::MomInertiaIntegrand(double in_r) 
// {
//   double out = pow(in_r, 4) ;
//         out *= 8*M_PI/3. ; 
//         out *= ( GetEps(in_r) + GetPress(in_r) ) ; 

//         // // converting fm^{-3} to km^{-3}
//         // out *= pow(10, 54) ;

//         out *= ::exp(-GetNu(in_r)) ;
//         out /= sqrt(1 - 2*GetMass(in_r) / in_r) ;

//   return out ;
// }

//--------------------------------------------------------------
/// Total moment of inertia (in km^3)
// MomInertiaIntegrand is wrong!
double CompactStar::NStar::Find_MomInertia() 
{
  PROFILE_FUNCTION() ;

  // double result = 0 ;

  // double err ;
  // gsl_integration_workspace *w = 
  //  gsl_integration_workspace_alloc(2000);

  // Zaki::Math::GSLFuncWrapper<NStar, double (NStar::*)(double)> 
  //     Fp(this, &NStar::MomInertiaIntegrand);     

  // gsl_function F = static_cast<gsl_function> (Fp) ; 

  // gsl_integration_qag(&F, radius[0], radius[-1], 1e-9, 1e-9, 2000, 1, w, &result, &err);



  // CompactStar::RotationSolver r_solver ;

  // r_solver.SetWrkDir(wrk_dir) ;

  // r_solver.AttachNStar(ns_shared_ptr) ;

  // r_solver.ImportTOVSolution({radius.vals, mass.vals, press.vals, eps.vals}) ;

  // r_solver.Solve({{5e-3, 5e-2}, 1, "Log"}) ;
  // r_solver.Solve(5e-3) ;

  // OmegaSeqPoint omega_seq_pt = r_solver.GetOmegaSeq()[0] ;
  // result = omega_seq_pt.J / ( omega_seq_pt.Omega / Zaki::Physics::LIGHT_C_KM_S ) ;

  // std::cout 
  //           << "J = " << omega_seq_pt.J 
  //           << ", Omega = " << omega_seq_pt.Omega << ","
  //           << " I = " << result << "\n" ;

  rot_solver.FindNMomInertia() ;

  return MomI ;
}

//--------------------------------------------------------------
/// Precision for printing the profile
/// by default it is set to '9' digits
void CompactStar::NStar::SetProfilePrecision(const int& in_prec) 
{
  profile_precision = in_prec ;
}

//--------------------------------------------------------------
// Exports the star profile
void CompactStar::NStar::Export(const Zaki::String::Directory& in_dir) 
{
  char seq_header[200] ;
  snprintf(seq_header, sizeof(seq_header),
          "    %-14s\t %-14s\t %-14s\t %-14s\t %-14s\t %-14s", 
          "ec(g/cm^3)", "M",  "R(km)", "pc(dyne/cm^2)", "B", "I(km^3)") ;
  
  std::time_t end_time = std::chrono::system_clock::to_time_t(
    std::chrono::system_clock::now()) ;

  ds.AddHead("# --------------------------------------------------------"
                 "--------------------------------------------------------\n") ;
  
  ds.AddHead("# Profile generated on ") ;
  ds.AddHead(std::string( std::ctime(&end_time)) ) ;
  ds.AddHead("# --------------------------------------------------------"
                 "--------------------------------------------------------\n") ;
  ds.AddHead("# Sequence point info: \n");
  ds.AddHead("#         " + std::string(seq_header) + "\n" ) ;
  ds.AddHead("#         " + sequence.Str() + "\n");
  ds.AddHead("# --------------------------------------------------------"
                 "--------------------------------------------------------\n") ;
  ds.AddFoot("# --------------------------------------------------------"
                 "--------------------------------------------------------") ;

  // std::cout << "\n\t profile_precision = " << profile_precision << "\n" ;
  ds.SetPrecision(profile_precision) ;
  ds.Export(in_dir.ThisFileDir() + "/" + in_dir.ThisFile().Str()) ;
  
  ds.ClearHeadFoot() ;
}

//--------------------------------------------------------------

//==============================================================
