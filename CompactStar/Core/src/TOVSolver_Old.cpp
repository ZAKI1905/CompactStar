/*
  Last edited on ?, 2022

  TOVSolver class
*/

// Creating directories
#include <sys/stat.h>

#include <array>

#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_const_cgsm.h>

#include <Zaki/Util/Instrumentor.hpp>
#include <Zaki/File/CSVIterator.hpp>
#include <Zaki/File/VecSaver.hpp>
#include <Zaki/Math/GSLFuncWrapper.hpp>
#include <Zaki/Physics/Constants.hpp>

#include "CompactStar/Core/NStar.hpp"
#include "CompactStar/Core/MixedStar.hpp"
#include "CompactStar/Core/TOVSolver.hpp"

#define TOV_SOLVER_VERBOSE 0

using namespace CompactStar ;

//==============================================================
//                        MixedSequence class
//==============================================================
// Constructor
MixedSequence::MixedSequence() : Prog("MixedSequence") 
{}
//--------------------------------------------------------------

void MixedSequence::Add(const MixedStar& in_star)
{
  seq.emplace_back ( in_star.sequence ) ;
}

//--------------------------------------------------------------
// Exports the mixed star sequence
void MixedSequence::Export(const Zaki::String::Directory& in_dir) 
const 
{
  Zaki::File::VecSaver vec_saver(wrk_dir + in_dir) ;

  char seq_header[400] ;
  snprintf(seq_header, sizeof(seq_header), "%-14s\t %-14s\t %-14s\t %-14s\t %-14s"
          "\t %-14s\t %-14s\t %-14s\t %-14s\t %-14s\t %-14s\t %-14s", 
          "ec(g/cm^3)", "M",  "R(km)", "pc(dyne/cm^2)", "B",
          "I(km^3)", "ec_d(g/cm^3)", "M_d",  "R_d(km)",
          "pc_d(dyne/cm^2)", "B_d", "I_d(km^3)" ) ;

  vec_saver.SetHeader(seq_header) ;

  vec_saver.Export1D(seq) ;
}
//--------------------------------------------------------------
// Combines two sequences
void MixedSequence::Combine(const MixedSequence& other) 
{
  // std::lock_guard<std::mutex> lock(m_mutex) ;

  seq.reserve( seq.size() + other.seq.size()) ;

  for (auto &&s : other.seq)
  {
    seq.emplace_back(s) ;
  }
  std::cout << "[ Thread = " << std::this_thread::get_id()
            << " ] " << name << " Size = " 
            << seq.size() << "\n" ;
}

//--------------------------------------------------------------

//==============================================================
//                        TOVSolver class
//==============================================================
// Constructor
TOVSolver::TOVSolver() : Prog("TOVSolver")
{
  mixed_r_accel = gsl_interp_accel_alloc () ;
  visi_p_accel  = gsl_interp_accel_alloc () ;
  dark_p_accel  = gsl_interp_accel_alloc () ;

  // tov_counter++ ;
}

//--------------------------------------------------------------
TOVSolver::~TOVSolver()
{

  if(visi_eps_p_spline)
    gsl_spline_free (visi_eps_p_spline);
  
  if(dark_eps_p_spline)
    gsl_spline_free (dark_eps_p_spline);
  
  if(visi_rho_p_spline)
    gsl_spline_free (visi_rho_p_spline);
  
  if(dark_rho_p_spline)
    gsl_spline_free (dark_rho_p_spline);
  
  // if(rho_r_spline)
  //   gsl_spline_free (rho_r_spline);
  
  // if(rho_r_spline_dark)
  //   gsl_spline_free (rho_r_spline_dark);
  
  if(nu_der_r_spline)  
    gsl_spline_free (nu_der_r_spline);
  
  if(mixed_r_accel)
    gsl_interp_accel_free (mixed_r_accel);

  // if(dark_accel)
  //   gsl_interp_accel_free (dark_accel);

  if(visi_p_accel)
    gsl_interp_accel_free (visi_p_accel);

  if(dark_p_accel)
    gsl_interp_accel_free (dark_p_accel);

  for(auto sp : visi_rho_i_p_spline)
  {
    if(sp)
      gsl_spline_free(sp) ;
  }
  for(auto sp : dark_rho_i_p_spline)
  {
    if(sp)
      gsl_spline_free(sp) ;
  }

  // tov_counter-- ;

  // {
  //   std::lock_guard<std::mutex> lock(m_mutex) ;

  //   mixed_seq_static.Combine(mixed_sequence) ;
  //   // std::cout << "\n\t tov_counter = " << tov_counter << "\n";
  // }

  // if( tov_counter == 0 )
  // {
  //   mixed_seq_static.Export("Full_Mixed_Sequence.txt") ;
  // }
}

//--------------------------------------------------------------
/// Sets the work directory for the member objects
Prog* TOVSolver::SetMemWrkDir(const Zaki::String::Directory& in_dir) 
{
  mixed_star.SetWrkDir(in_dir) ;
  mixed_sequence.SetWrkDir(in_dir) ;

  // {
  //   std::lock_guard<std::mutex> lock(m_mutex) ;

  //   if ( !mixed_seq_static.IsWrkDirSet() )
  //   {
  //     mixed_seq_static.SetWrkDir(in_dir) ;
  //   }
  // } 

  return this ;
}

//--------------------------------------------------------------
void TOVSolver::ImportEOS(const Zaki::String::Directory& f_name) 
{
  std::ifstream     file( (wrk_dir + "/" + f_name).Str());

  // Error opening the file
  if (file.fail()) 
  {
    Z_LOG_ERROR("File '"+(wrk_dir + "/" + f_name).Str() +"' cannot be opened!") ;
    Z_LOG_ERROR("Importing EOS data failed!") ;
    exit(EXIT_FAILURE) ;
    return ;
  }

  size_t line_num = 0 ;
  // Reading the input file
   ;
  for(Zaki::File::CSVIterator loop(file, '\t'); loop != Zaki::File::CSVIterator(); ++loop)
  {
    if( (*loop).size() < 3)
    {
      std::cout << "\n (*loop)[0]: " << (*loop)[0] << "\n";
      std::cout << "(*loop).size(): " << (*loop).size() << "\n";
      Z_LOG_ERROR("EOS file is not complete!") ;
      break ;
    }

    // First line
    if (line_num == 0)
    { 
      eos_tab.SetLabels((*loop)[0], (*loop)[1], (*loop)[2]) ;
      
      // We want to do this once only
      if ( (*loop).size() > 3)
      {
        for (size_t i = 3; i < (*loop).size(); i++)
        {
          eos_tab.rho_i.push_back({}) ;
          eos_tab.AddExtraLabels((*loop)[i]) ;
        }
      }
    }
    else
    {
      eos_tab.eps.push_back(std::atof((*loop)[0].c_str())) ;
      eos_tab.pre.push_back(std::atof((*loop)[1].c_str())) ;
      eos_tab.rho.push_back(std::atof((*loop)[2].c_str())) ;
      
      for( size_t i=3 ; i < (*loop).size() ; i++)
        eos_tab.rho_i[i-3].push_back(
          std::atof((*loop)[i].c_str())
          ) ;
    }
    line_num++ ;
  }

  Z_LOG_INFO("EOS data imported from: "+ (wrk_dir+f_name).Str()+".") ;
  visi_eps_p_spline  = gsl_spline_alloc (gsl_interp_steffen, eos_tab.Size());
  visi_rho_p_spline = gsl_spline_alloc (gsl_interp_steffen, eos_tab.Size()); 

  for(size_t i=0 ; i < eos_tab.rho_i.size() ; i++)
  {
    visi_rho_i_p_spline.emplace_back(gsl_spline_alloc(gsl_interp_steffen, eos_tab.Size())) ; 
    gsl_spline_init(visi_rho_i_p_spline[i], &eos_tab.pre[0], &eos_tab.rho_i[i][0], eos_tab.Size() ) ;
  }

  // This function initializes the interpolation object
  // x has to be strictly increasing
  gsl_spline_init (visi_eps_p_spline, &eos_tab.pre[0], &eos_tab.eps[0], eos_tab.Size());
  gsl_spline_init (visi_rho_p_spline, &eos_tab.pre[0], &eos_tab.rho[0], eos_tab.Size());

  mixed_star.InitVisible(this) ;
}

//--------------------------------------------------------------
void TOVSolver::ImportEOS_Dark(const Zaki::String::Directory& f_name) 
{
  std::ifstream     file( (wrk_dir + "/" + f_name).Str());

  // Error opening the file
  if (file.fail()) 
  {
    Z_LOG_ERROR("File '"+(wrk_dir + "/" + f_name).Str() +"' cannot be opened!") ;
    Z_LOG_ERROR("Importing EOS data failed!") ;
    exit(EXIT_FAILURE) ;
    return ;
  }

  size_t line_num = 0 ;
  // Reading the input file
   ;
  for(Zaki::File::CSVIterator loop(file, '\t'); loop != Zaki::File::CSVIterator(); ++loop)
  {
    if( (*loop).size() < 3)
    {
      std::cout << "\n (*loop)[0]: " << (*loop)[0] << "\n";
      std::cout << "(*loop).size(): " << (*loop).size() << "\n";
      Z_LOG_ERROR("EOS file is not complete!") ;
      break ;
    }

    // First line
    if (line_num == 0)
    { 
      eos_tab_dark.SetLabels((*loop)[0], (*loop)[1], (*loop)[2]) ;
      
      // We want to do this once only
      if ( (*loop).size() > 3)
      {
        for (size_t i = 3; i < (*loop).size(); i++)
        {
          eos_tab_dark.rho_i.push_back({}) ;
          eos_tab_dark.AddExtraLabels((*loop)[i]) ;
        }
      }
    }
    else
    {
      eos_tab_dark.eps.push_back(std::atof((*loop)[0].c_str())) ;
      eos_tab_dark.pre.push_back(std::atof((*loop)[1].c_str())) ;
      eos_tab_dark.rho.push_back(std::atof((*loop)[2].c_str())) ;
      
      for( size_t i=3 ; i < (*loop).size() ; i++)
        eos_tab_dark.rho_i[i-3].push_back(
          std::atof((*loop)[i].c_str())
          ) ;
    }
    line_num++ ;
  }

  Z_LOG_INFO("Dark EOS data imported from: "+ (wrk_dir+f_name).Str()+".") ;

  dark_eps_p_spline  = gsl_spline_alloc (gsl_interp_steffen, eos_tab_dark.Size());
  dark_rho_p_spline = gsl_spline_alloc (gsl_interp_steffen, eos_tab_dark.Size()); 

  for(size_t i=0 ; i < eos_tab_dark.rho_i.size() ; i++)
  {
    dark_rho_i_p_spline.emplace_back(gsl_spline_alloc(gsl_interp_steffen, eos_tab_dark.Size())) ; 
    gsl_spline_init(dark_rho_i_p_spline[i], &eos_tab_dark.pre[0], &eos_tab_dark.rho_i[i][0], eos_tab_dark.Size() ) ;
  }

  // This function initializes the interpolation object
  // x has to be strictly increasing
  gsl_spline_init (dark_eps_p_spline, &eos_tab_dark.pre[0], &eos_tab_dark.eps[0], eos_tab_dark.Size());
  gsl_spline_init (dark_rho_p_spline, &eos_tab_dark.pre[0], &eos_tab_dark.rho[0], eos_tab_dark.Size());

  mixed_star.InitDark(this) ;
}

//--------------------------------------------------------------
void TOVSolver::PrintEOSTable() const 
{
  eos_tab.Print() ;
}

//--------------------------------------------------------------
//            Added on December 15, 2020
/// Returns the nu_der value given the radius input
// r is in cm!
double TOVSolver::GetNuDerSpline(const double& in_r) 
{
  return gsl_spline_eval(nu_der_r_spline, in_r, mixed_r_accel);
}

//--------------------------------------------------------------
//            Added on December 15, 2020
/// Evaluates & exports the nu(r) function table
/// input must be a TOV solution file (with nu' column)
// void TOVSolver::ExportNu(const Zaki::String::Directory& f_name)
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

//   std::vector<double> tov_radius ;
//   std::vector<double> tov_nu_der ;
//   std::vector<double> tov_mass   ;

//   size_t line_num = 0 ;
//   for(Zaki::File::CSVIterator loop(file, '\t'); loop != Zaki::File::CSVIterator(); ++loop)
//   {
//     // First line
//     if (line_num == 0)
//     {
//       if( Zaki::String::Strip((*loop)[2], ' ') != "nu'")
//       {
//         Z_LOG_ERROR("TOV solution doesn't include nu'(r) !") ;
//         return ;
//       }
//     }
//     // Other lines:
//     else
//     {
//       // Converting the radius into cm --> Why? 
//       // Mar-2022: Because nu' is in (1/cm).
//       tov_radius.push_back(std::atof((*loop)[0].c_str()) * 1e5) ;
//       tov_mass.push_back(std::atof((*loop)[1].c_str())) ;
//       tov_nu_der.push_back(std::atof((*loop)[2].c_str())) ;
//     }
//     line_num++ ;
//   }

//   Z_LOG_INFO("TOV soluton imported from: "+ (wrk_dir+f_name).Str()+".") ;

//   nu_der_r_spline  = gsl_spline_alloc (gsl_interp_steffen, tov_radius.size());

//   // This function initializes the interpolation object
//   // x has to be strictly increasing
//   gsl_spline_init (nu_der_r_spline, &tov_radius[0], &tov_nu_der[0], tov_radius.size());

//   // -----------------------------------
//   // Integrate to find nu(r) :
//   // -----------------------------------

//   gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000*tov_radius.size());
//   double err;

//   Zaki::Math::GSLFuncWrapper<TOVSolver, double (TOVSolver::*)( const double& )> 
//     Fp(this, &TOVSolver::GetNuDerSpline);     

//   gsl_function F = static_cast<gsl_function> (Fp) ; 

//   double tmp_nu = 0;
//   std::vector<double> nu_result ;
//   nu_result.reserve(tov_radius.size()) ;

//   for(size_t i = 0 ; i < tov_radius.size() ; ++i)
//   {
//     gsl_integration_qag(&F, tov_radius[0], tov_radius[i], 1e-9, 1e-9, 1000, 1, w, &tmp_nu, &err);
//     nu_result.push_back(tmp_nu) ;
//     tmp_nu = 0 ;
//   }
  
//   gsl_integration_workspace_free(w);

  
//   // -------------------------------------------
//   // Saving the results for [ r, M(r), nu(r) ] 
//   // -------------------------------------------
  
//   // Matching the boundary condition for nu(r)
//   double nu_at_R = 0.5*log(1 
//                   - 2*Zaki::Physics::SUN_M_KM*tov_mass[tov_mass.size()-1]
//                     / (tov_radius[tov_mass.size()-1]*1e-5)
//                   ) ;

//   double delta_nu_r = nu_result[tov_mass.size()-1] - nu_at_R ;
//   std::vector<TOV_Nu_Point> rmnu_results ;
//   rmnu_results.reserve(tov_radius.size()) ;

//   for (size_t i = 0; i < tov_radius.size(); i++)
//   {
//     rmnu_results.emplace_back(tov_radius[i]*1e-5, 
//                     tov_mass[i],
//                     nu_result[i] - delta_nu_r) ;
//   }
  
//   std::string out_f_name = Zaki::String::Pars((wrk_dir + "/" + f_name).Str(), ".tsv")[0] ;

//   Zaki::File::VecSaver vec_saver(out_f_name + "_nu.tsv");

//   char seq_header[200] ;
//   sprintf(seq_header, "%-14s\t %-14s\t %-14s", 
//           "r(km)", "m",  "nu") ;
//   vec_saver.SetHeader(seq_header) ;
//   vec_saver.Export1D(rmnu_results) ;

// }

//--------------------------------------------------------------
// Input pressure, output energy density
double TOVSolver::GetEDens(const double& in_pres) 
{
  return gsl_spline_eval(visi_eps_p_spline, in_pres, visi_p_accel) ; 
}

//--------------------------------------------------------------
// Input pressure, output energy density (dark sector)
double TOVSolver::GetEDens_Dark(const double& in_pres) 
{
  return gsl_spline_eval(dark_eps_p_spline, in_pres, dark_p_accel);
}

//--------------------------------------------------------------
double TOVSolver::cost_p_of_e(const double in_p)
{
  return GetEDens(in_p) - cost_p_of_e_input ;
}

//--------------------------------------------------------------
double TOVSolver::cost_p_of_e_dark(const double in_p)
{
  return GetEDens_Dark(in_p) - cost_p_of_e_input_dark ;
}

//--------------------------------------------------------------
// Inverse function of "GetEDens"
double TOVSolver::p_of_e(const double& in_e) 
{
  double p_min = eos_tab.pre[0] ;
  double p_max = eos_tab.pre[eos_tab.pre.size()-1] ;

  // std::cout << "p_min = " << p_min 
  //           << ", p_max = " << p_max ;

  // std::cout << ", p_of_e(" <<in_e << ") = " ;

  cost_p_of_e_input = in_e ;

  Zaki::Math::GSLFuncWrapper<TOVSolver, double (TOVSolver::*)(double)> 
  func(this, &TOVSolver::cost_p_of_e) ;

  gsl_function F = static_cast<gsl_function> (func) ; 

  const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
  gsl_root_fsolver *s = gsl_root_fsolver_alloc (T);

//  gsl_set_error_handler_off() ;
  gsl_root_fsolver_set (s, &F, p_min, p_max);

  int iter = 0, max_iter = 250 ;
  int status ; double out_p = 0 ;
  do
  {
    iter++;
    status = gsl_root_fsolver_iterate (s);
    out_p = gsl_root_fsolver_root (s);
    p_min = gsl_root_fsolver_x_lower (s);
    p_max = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (p_min, p_max,
                                      0.0001, 0.0001);
  }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free (s);

  // std::cout << out_p << ", GetEDens(" 
  //           << out_p << ") = " 
  //           << GetEDens(out_p) << ", iter = " << iter 
  //           << ", p_min = " << p_min 
  //           << ", p_max = " << p_max << "\n" ;
  
  return out_p ;
}

//--------------------------------------------------------------
// Inverse function of "GetEDens_Dark"
double TOVSolver::p_of_e_dark(const double& in_e) 
{
  double p_min = eos_tab_dark.pre[0] ;
  double p_max = eos_tab_dark.pre[eos_tab_dark.pre.size()-1] ;

  cost_p_of_e_input_dark = in_e ;

  // std::cout << "p_min = " << p_min 
  //           << ", p_max = " << p_max ;

  // std::cout << ", p_of_e(" <<in_e << ") = " ;

  Zaki::Math::GSLFuncWrapper<TOVSolver, double (TOVSolver::*)(double)> 
  func(this, &TOVSolver::cost_p_of_e_dark) ;

  gsl_function F = static_cast<gsl_function> (func) ; 

  const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
  gsl_root_fsolver *s = gsl_root_fsolver_alloc (T);

  gsl_set_error_handler_off() ;
  gsl_root_fsolver_set (s, &F, p_min, p_max);

  int iter = 0, max_iter = 250 ;
  int status ; double out_p = 0 ;
  do
  {
    iter++;
    status = gsl_root_fsolver_iterate (s);
    out_p = gsl_root_fsolver_root (s);
    p_min = gsl_root_fsolver_x_lower (s);
    p_max = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (p_min, p_max,
                                      0.0001, 0.0001);
  }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free (s);

  // std::cout << out_p << ", GetEDens(" 
  //           << out_p << ") = " 
  //           << GetEDens_Dark(out_p) << ", iter = " << iter 
  //           << ", p_min = " << p_min 
  //           << ", p_max = " << p_max << "\n" ;

  return out_p ;
}

//--------------------------------------------------------------
double TOVSolver::PressureCutoff() const
{
  // return std::max(1.e-15 * GetInitPress(), eos_tab.pre[0]) ;
  return eos_tab.pre[0] ;
}

//--------------------------------------------------------------
double TOVSolver::PressureCutoff_Dark() const
{
  return eos_tab_dark.pre[0] ;
}

//--------------------------------------------------------------
// Modified ODE in the presence of a dark core 
// ......................
// Dictionary :
// y[0] = visible mass(r)
// y[1] = dark mass(r)
// y[2] = visible pressure(r)
// y[3] = dark pressure(r)
// f[0] = visible m'(r)
// f[1] = dark m'(r)
// f[2] = visible p'(r)
// f[3] = dark p'(r)
// ......................
int TOVSolver::ODE_Dark_Core(double r, const double y[], double f[], void *params)
{
  TOVSolver* tov_obj = (TOVSolver *)params ; 
  
  // Visible surface reached
  if (y[2] < tov_obj->PressureCutoff())
  {
#if TOV_SOLVER_VERBOSE
    printf ("\u2554----------------- Visible Core Surface reached ----------------\u2557\n");
    printf ("\u2551 %14s %14s %15s %20s\n", "R (km)", "M (M_sun)", "\u03B5_c (g/cm^3)", "\u2551");
    printf ("\u2551 %14le %14le %14le %20s\n", r / 1.e+5,
      y[0] / GSL_CONST_CGSM_SOLAR_MASS, tov_obj->GetEDens(tov_obj->GetInitPress()), "\u2551");
    printf ("\u255A---------------------------------------------------------------\u255D\n");

    printf ("\n --------------------------------------------------\n");
      std::cout << " Vis. P_0 =" << tov_obj->GetInitPress()<< "\n";
      std::cout << " Vis. P Cut-Off =" << tov_obj->PressureCutoff()<< "\n";
      std::cout << " Current Vis. P = y[2] = " << y[2] << "\n" ;
    printf (" --------------------------------------------------\n\n");
#endif
    tov_obj->dark_core = false ;

    return GSL_EBADFUNC ;	
  }

  // Dark core surface reached
  if (y[3] < tov_obj->PressureCutoff_Dark())
  {
#if TOV_SOLVER_VERBOSE
    printf ("\u2554----------------- Dark Core Surface reached ----------------\u2557\n");
    printf ("\u2551 %14s %14s %15s %17s\n", "R (km)", "M (M_sun)", "\u03B5_c (g/cm^3)", "\u2551");
    printf ("\u2551 %14le %14le %14le %17s\n", r / 1.e+5,
      y[1] / GSL_CONST_CGSM_SOLAR_MASS, tov_obj->GetEDens_Dark(tov_obj->GetInitPress_Dark()), "\u2551");
    printf ("\u255A------------------------------------------------------------\u255D\n");

    printf ("\n --------------------------------------------------\n");
      std::cout << " Dark P_0 =" << tov_obj->GetInitPress_Dark()<< "\n";
      std::cout << " Dark P Cut-Off =" << tov_obj->PressureCutoff_Dark()<< "\n";
      std::cout << " Current Dark P = y[3] = " << y[3] << "\n" ;
    printf (" --------------------------------------------------\n\n");
#endif
    tov_obj->dark_core = true ;

    return GSL_EBADFUNC ;
  }

  // Mass Continuity Equations
  f[0] = 4. * M_PI * r * r * tov_obj->GetEDens(y[2]) ;
  f[1] = 4. * M_PI * r * r * tov_obj->GetEDens_Dark(y[3]) ;

  // TOV Equation
  // visible pressure derivative
  f[2] = -(GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT / pow (r, 2.))
    * (tov_obj->GetEDens(y[2])
       + (y[2] / pow (GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.)))
    * (y[0] + y[1] + 4 * M_PI * pow (r, 3.)
       * (y[2] + y[3] ) / pow (GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.))
    / (1. - (2. * GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT *
              (y[0] + y[1])
	      / (pow (GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.) * r)));

  // dark pressure derivative
  f[3] = -(GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT / pow (r, 2.))
    * (tov_obj->GetEDens_Dark(y[3])
       + (y[3] / pow (GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.)))
    * (y[0] + y[1] + 4 * M_PI * pow (r, 3.)
       * (y[2] + y[3] ) / pow (GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.))
    / (1. - (2. * GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT *
              (y[0] + y[1])
	      / (pow (GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.) * r)));

  return GSL_SUCCESS;
}

//--------------------------------------------------------------
// Modified ODE in the presence of a dark core 
// ......................
// Dictionary :
// y[0] = mantle mass(r)
// y[1] = mantle pressure(r)
// f[0] = mantle m'(r)
// f[1] = mantle p'(r)
// ......................
int TOVSolver::ODE_Dark_Mantle(double r, const double y[], double f[], void *params)
{
  TOVSolver* tov_obj = (TOVSolver *)params ; 
  double m_c = tov_obj->m_core ;

  // Surface reached
  if (tov_obj->dark_core && y[1] < tov_obj->PressureCutoff())
  {
#if TOV_SOLVER_VERBOSE
    printf ("\u2554----------------- Vis. Mantle Surface reached ----------------\u2557\n");
    printf ("\u2551 %14s %14s %15s %19s\n", "R (km)", "M (M_sun)", "\u03B5_c (g/cm^3)", "\u2551");
    printf ("\u2551 %14le %14le %14le %19s\n", r / 1.e+5,
      y[0] / GSL_CONST_CGSM_SOLAR_MASS, tov_obj->GetEDens(tov_obj->GetInitPress()), "\u2551");
    printf ("\u255A--------------------------------------------------------------\u255D\n");

    printf ("\n --------------------------------------------------\n");
      std::cout << " Vis. P_0 =" << tov_obj->GetInitPress()<< "\n";
      std::cout << " Vis. P Cut-Off =" << tov_obj->PressureCutoff()<< "\n";
      std::cout << " Current Vis. P = y[1] = " << y[1] << "\n" ;
    printf (" --------------------------------------------------\n\n");
#endif
    return GSL_EBADFUNC ;
  }
  if (!(tov_obj->dark_core) && y[1] < tov_obj->PressureCutoff_Dark())
  {
#if TOV_SOLVER_VERBOSE
    printf ("\u2554----------------- Dark Mantle Surface reached ----------------\u2557\n");
    printf ("\u2551 %14s %14s %15s %19s\n", "R (km)", "M (M_sun)", "\u03B5_c (g/cm^3)", "\u2551");
    printf ("\u2551 %14le %14le %14le %19s\n", r / 1.e+5,
      y[0] / GSL_CONST_CGSM_SOLAR_MASS, tov_obj->GetEDens_Dark(tov_obj->GetInitPress_Dark()), "\u2551");
    printf ("\u255A--------------------------------------------------------------\u255D\n");
    
    printf ("\n --------------------------------------------------\n");
      std::cout << " Dark P_0 =" << tov_obj->GetInitPress_Dark()<< "\n";
      std::cout << " Dark P Cut-Off =" << tov_obj->PressureCutoff_Dark()<< "\n";
      std::cout << " Current Dark P = y[1] = " << y[1] << "\n" ;
    printf (" --------------------------------------------------\n\n");
#endif
    return GSL_EBADFUNC ;
  }

  double e_den = 0 ;

  if(tov_obj->dark_core)
    e_den = tov_obj->GetEDens(y[1]) ;
  else
    e_den = tov_obj->GetEDens_Dark(y[1]) ;

  // Mass Continuity Equations
  f[0] = 4. * M_PI * r * r * e_den ;

  // TOV Equation
  // visible pressure derivative
  f[1] = -(GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT / pow (r, 2.))
    * (e_den
       + (y[1] / pow (GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.)))
    * (y[0] + m_c + 4 * M_PI * pow (r, 3.)
       * (y[1] + 0 ) / pow (GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.))
    / (1. - (2. * GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT *
              (y[0] + m_c)
	      / (pow (GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.) * r)));


  return GSL_SUCCESS;
}

//--------------------------------------------------------------
// ......................
// Dictionary :
// y[0] = mass(r)
// y[1] = pressure(r)
// f[0] = m'(r)
// f[1] = p'(r)
// ......................
int TOVSolver::ODE(double r, const double y[], double f[], void *params)
{
  TOVSolver* tov_obj = (TOVSolver *)params ; 

  // set a minimum pressure cutoff. if we don't, the ODE solver will wobble all
  // over the surface and crash, or if you make the error tolerance really strict
  // it'll integrate forever
  if (y[1] < tov_obj->PressureCutoff())
  {
#if TOV_SOLVER_VERBOSE
    printf ("\u2554----------------- Surface reached ----------------\u2557\n");
    printf ("\u2551 %14s %14s %15s %7s\n", "R (km)", "M (M_sun)", "\u03B5_c (g/cm^3)", "\u2551");
    printf ("\u2551 %14le %14le %14le %7s\n", r / 1.e+5,
      y[0] / GSL_CONST_CGSM_SOLAR_MASS, tov_obj->GetEDens(tov_obj->GetInitPress()), "\u2551");
    printf ("\u255A--------------------------------------------------\u255D\n");
      std::cout << "\n init_press =" << tov_obj->GetInitPress()<< "\n";
      std::cout << "\n Pressure Cut Off =" << tov_obj->PressureCutoff()<< "\n";
      std::cout << " y[1] = " << y[1] << "\n" ;
#endif
    return GSL_EBADFUNC;	// this flag tells GSL integrator to quit
  }


  // Mass Continuity Equation
  f[0] = 4. * M_PI * r * r * tov_obj->GetEDens(y[1]);

  // TOV Equation
  //      - [ eps(r) + p/c^2 ]
  // f[1]  = - eos_eps(y[1])
  //       + (y[1] / pow (GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.));
  
  // //    [ M(r) + 4.pi.r^3 p(r) ]
  // f[1] *= y[0] + 4*M_PI*r*r*r*y[1] / pow(GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.) ;

  // //    / [ r.(r - 2M(r)) ]
  // f[1]  *= 1./r ;
  // f[1]  *= 1./( r/GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT 
  //               - 2*y[0]/pow (GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.)) ;

  f[1] = -(GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT / pow (r, 2.))
    * (tov_obj->GetEDens(y[1])
       + (y[1] / pow (GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.)))
    * (y[0] + 4 * M_PI * pow (r, 3.)
       * y[1] / pow (GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.))
    / (1. - (2. * GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * y[0]
	      / (pow (GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.) * r)));

  return GSL_SUCCESS;
}

//--------------------------------------------------------------
//               Added on December 15, 2020
/// Returns the derivative of the metric nu(r) function
// r must be in cm!
double TOVSolver::GetNuDer(const double r, const std::vector<double>& y) 
{
  // (dp/dr) has units of  [ g / (cm^2 s^2) ] 
  double dpdr = -(GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT / pow (r, 2.))
    * (GetEDens(y[1])
       + (y[1] / pow (GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.)))
    * (y[0] + 4 * M_PI * pow (r, 3.)
       * y[1] / pow (GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.))
    / (1. - (2. * GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * y[0]
	      / (pow (GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.) * r)));

  return - dpdr / ( 
          
          y[1] +  pow (GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.)*GetEDens(y[1])

                  ) ;
}

//--------------------------------------------------------------
double TOVSolver::GetNuDer_Dark(const double r, const std::vector<double>& y) 
{

  double out = (GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT / pow (r, 2.))
    * (y[0] + 4 * M_PI * pow (r, 3.)
       * y[1] / pow (GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.))
    / (1. - (2. * GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * y[0]
	      / (pow (GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.) * r))) ;

  out /= pow (GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.) ;

  return out ;
}

//--------------------------------------------------------------
/// Returns the total baryon number density given pressure
double TOVSolver::GetRho(const double& in_p) 
{
  return gsl_spline_eval(visi_rho_p_spline, in_p, visi_p_accel) ;
}

//--------------------------------------------------------------
/// Returns the total baryon number density given pressure
double TOVSolver::GetRho_Dark(const double& in_p) 
{
  return gsl_spline_eval(dark_rho_p_spline, in_p, dark_p_accel) ;
}

//--------------------------------------------------------------
/// Returns the specific number density given pressure
std::vector<double> TOVSolver::GetRho_i(const double& in_p) 
{
  std::vector<double> out ;
  out.reserve(visi_rho_i_p_spline.size()) ;

  for (auto sp : visi_rho_i_p_spline)
  {
    out.push_back(gsl_spline_eval(sp, in_p, visi_p_accel)) ;
  }
  
  return out ;
}

//--------------------------------------------------------------
/// Returns the specific number density given pressure
std::vector<double> TOVSolver::GetRho_i_Dark(const double& in_p) 
{
  std::vector<double> out ;
  out.reserve(dark_rho_i_p_spline.size()) ;

  for (auto sp : dark_rho_i_p_spline)
  {
    out.push_back(gsl_spline_eval(sp, in_p, dark_p_accel)) ;
  }
  
  return out ;
}

//--------------------------------------------------------------
double TOVSolver::GetInitPress() const 
{
  return  init_press ;
}

//--------------------------------------------------------------
double TOVSolver::GetInitPress_Dark() const 
{
  return  init_press_dark ;
}

//--------------------------------------------------------------
double TOVSolver::GetInitEDens() const 
{
  return  init_edens ;
}

//--------------------------------------------------------------
// /// Returns the total baryon number density given radius
// double TOVSolver::GetRho_r(const double& in_r) 
// {
//   if (results.size() == 0)
//   {
//     Z_LOG_ERROR("Call 'GetRho_r' only after the"
//                 " TOV is fully solved for a star.") ;
//     return 0 ;
//   }
  
//   std::vector<double> r_list ;
//   std::vector<double> rho_list ;

//   rho_list.reserve(results.size()) ;
//   r_list.reserve(results.size()) ;

//   for (auto &&i : results)
//   {
//     r_list.emplace_back(i.r) ;
//     rho_list.emplace_back(i.rho) ;
//   }

//   rho_r_spline = gsl_spline_alloc (gsl_interp_cspline, results.size()); 

//   gsl_spline_init (rho_r_spline, &r_list[0], &rho_list[0], results.size()) ;

//   return gsl_spline_eval(rho_r_spline, in_r, accel) ;
// }

//--------------------------------------------------------------
/// Adds the condition for printing the mixed star profile
void TOVSolver::AddCondition(bool (*func)(const MixedStar&)) 
{
  exp_cond_f = func ;
}

//--------------------------------------------------------------
void TOVSolver::Solve(const Zaki::Math::Axis& in_ax,
                      const Zaki::String::Directory& in_dir)  
{
  // PROFILE_FUNCTION() ;

  std::cout << "\n\n\t\t ****************************************"
            <<"*******************************"<<" \n" ;
  std::cout <<     "\t\t *                     "
            <<"TOV Solver Sequence Results"
            <<"                     * \n" ;
  std::cout <<     "\t\t ******************************************"
            <<"*****************************"<<"\n\n" ;

  // This star object saves the TOV results
  auto tmp_star = std::make_shared<NStar> () ;
  tmp_star->Reserve(1000) ;

  for (size_t idx = 0; idx <= in_ax.res ; idx++)
  {
    Z_LOG_INFO("Sequence " + std::to_string(idx+1) 
                + " out of "+ std::to_string(in_ax.res+1) +".") ;

    // init_edens = in_ax[idx] ;
    init_press = p_of_e(in_ax[idx])  ;

    double r = r_min ; 
    // double r1 = r_max ;
    double y[2] ; 
    y[1] = init_press;
    y[0] = (4. / 3.) * M_PI * pow (r, 3.) * GetEDens(y[1]) ;

    gsl_odeiv2_system ode_sys = {CompactStar::TOVSolver::ODE, nullptr, 2, this} ;

    gsl_odeiv2_driver *tmp_driver = gsl_odeiv2_driver_alloc_y_new
      (&ode_sys, gsl_odeiv2_step_rk8pd,
      1.e-1, 1.e-10, 1.e-10);


    double min_log_r = log10(r_min) ;
    double max_log_r = log10(r_max) ;

    double step = (max_log_r - min_log_r) / 1000 ;
    double step_scale = 1 ; // Adaptive steps (Aug 6, 2020)

    // results.reserve(1000) ; // (deleted on Apr 22, 2022)

    for (double log_r_i = min_log_r;  log_r_i <= max_log_r; log_r_i += step*step_scale)
    {
      double ri = pow(10, log_r_i) ;

      double tmp_delta_p = y[1] ;     // Adaptive steps (Aug 6, 2020)

      // This function evolves the driver system d from t to t1.
      // Initially vector y should contain the values of dependent
      // variables at point t. 
      int status = gsl_odeiv2_driver_apply (tmp_driver, &r, ri, y);

      if (status != GSL_SUCCESS)
      {
        printf ("error, return value=%d\n", status);
        printf("Pressure = %2.2e.\n", y[1]) ;
        break;
      }

      // ................................................................
      //       Determining the step size adaptively (Aug 6, 2020)
      //       The steps are becoming too large, cap them! (Nov 8, 2021)
      // ................................................................
      tmp_delta_p = tmp_delta_p - y[1] ;
      // Adapting the step size if changes are
      if( tmp_delta_p/y[1] < 1e-8 ) // too small
      {
        // double the step-size
        step_scale *= 1.5 ;
      }
      // going back to normal scaling 
      else if ( tmp_delta_p/y[1] < 1e-6 )
      {
        step_scale = 1 ;
      }
      // too big
      else 
      {
        // scale step-size
        step_scale = pow(10.,-6.-log10(tmp_delta_p/y[1])) ;
      }

      // Checking if the next radius is closer than 1e-8 to this ri
      while(1e-3 > (pow(10, log_r_i + step*step_scale) - ri ) )
      {
        // std::cout << " -> ri = " << ri << "\n" ;
        step_scale *= 1.5 ;
      }

      // Added on Nov 8, 2021
      // Checking if the next radius is further than 5 (m) to this ri
      while(5e+2 < (pow(10, log_r_i + step*step_scale) - ri ) )
      {
        // std::cout << " -> ri = " << ri << "\n" ;
        // std::cout << " -> Delta_r = " << pow(10, log_r_i + step*step_scale) - ri  << "\n" ;
        step_scale /= 1.5 ;
      }
      // ................................................................

      // We set nu(r) = 0, and calculate it after the loop
      // results.emplace_back( r/1.e+5, y[0]/GSL_CONST_CGSM_SOLAR_MASS, 
      //                       GetNuDer(r, {y[0], y[1]}), 0,
      //                       y[1], GetEDens(y[1]), 
      //                       GetRho(y[1]), GetRho_i(y[1]));
      tmp_star->Append({
                        r/1.e+5, y[0]/GSL_CONST_CGSM_SOLAR_MASS, 
                        GetNuDer(r, {y[0], y[1]}), 0,
                        y[1], GetEDens(y[1]), 
                        GetRho(y[1]), GetRho_i(y[1])
                        });
    }
    // ----------------------------------------------------------------
    //                    RADIUS LOOP ENDS HERE
    // ----------------------------------------------------------------
    // tmp_star->SurfaceIsReached() ;
    // tmp_star->EvaluateNu() ;
    // tmp_star->SetSharedPtr(tmp_star) ;

    // Updating results values for nu(r) :
    // for (auto &&r_i : results)
    // {
    //   r_i.nu = tmp_star->GetNu(r_i.r) ;
    // }

    gsl_odeiv2_driver_free (tmp_driver);

    sequence.emplace_back(GetEDens(init_press), y[0]/GSL_CONST_CGSM_SOLAR_MASS,
                          r/1.e+5, init_press,
                          0, 0
                          // tmp_star->GetBaryonNum(),
                          // tmp_star->GetMomInertia()
                          );


    std::vector<std::string> tmp_fname_v = Zaki::String::Pars(in_dir.Str(), "*") ;
    std::string tmp_fname ;
    if(tmp_fname_v.size() <2 )
    {
      Z_LOG_WARNING("File name pattern doesn't match '[]*[]'.") ;
      tmp_fname = tmp_fname_v[0] + std::to_string(idx) ;
    }
    else
    {
      tmp_fname = tmp_fname_v[0] + std::to_string(idx) + tmp_fname_v[1] ;
    }
     
    Zaki::File::VecSaver vec_saver(wrk_dir + "/" + tmp_fname) ;
    char res_header[200] ;
    snprintf(res_header, sizeof(res_header), "%-14s\t %-14s\t %-14s\t %-15s\t %-14s\t %-14s\t %-14s", 
            "r(km)", "m", "nu'(1/cm)", " nu", "p(dyne/cm^2)", "e(g/cm^3)", "rho(1/fm^3)" ) ;

    std::string tmp_label = res_header ;
    for(auto& lab : eos_tab.extra_labels)
    {
      snprintf(res_header, sizeof(res_header), "\t %-14s", Zaki::String::Strip(lab, ' ').c_str() ) ;
      tmp_label += res_header ;
    }
      
    vec_saver.SetHeader(tmp_label.c_str()) ;
    // vec_saver.Export1D(results) ;

    // Testing this for freeing memory
    // results = std::vector<TOVPoint>();
    // results.clear() ;
    // results.shrink_to_fit() ;

    // tmp_star->Reset() ;
  
  } // TOV sequence loop ends!

  std::cout <<"\n\t\t *************************"
            <<" TOV Solver Finished *************************"<<"\n\n" ;
}

//--------------------------------------------------------------
void TOVSolver::Solve_Mixed( const Zaki::Math::Axis& in_v_ax,
                            const Zaki::Math::Axis& in_d_ax,
                      const Zaki::String::Directory& in_dir,
                      const Zaki::String::Directory& in_file)  
{
  PROFILE_FUNCTION() ;

#if TOV_SOLVER_VERBOSE
  std::cout << "\n\n\t\t ****************************************"
            <<"*******************************"<<" \n" ;
  std::cout <<     "\t\t *                     "
            <<"TOV Solver Sequence Results"
            <<"                     * \n" ;
  std::cout <<     "\t\t ******************************************"
            <<"*****************************"<<"\n\n" ;
#endif

  // ----------------------------------------------------------------
  //                  TOV Dark sequence loop begins
  // ----------------------------------------------------------------
  for (size_t d_idx = 0; d_idx <= in_d_ax.res ; d_idx++)
  {
    Z_LOG_INFO("Dark sequence " + std::to_string(d_idx+1) 
                  + " out of "+ std::to_string(in_d_ax.res+1) +".") ;

    init_press_dark = p_of_e_dark(in_d_ax[d_idx])  ;
    // ----------------------------------------------------------------
    //                  TOV Visible sequence loop begins
    // ----------------------------------------------------------------
    for (size_t v_idx = 0; v_idx <= in_v_ax.res ; v_idx++)
    {
      Z_LOG_INFO("Mixed sequence (" + std::to_string(d_idx+1) + ", " 
                  + std::to_string(v_idx+1) 
                  + ") out of ("+ std::to_string(in_d_ax.res+1) 
                  + ", " + std::to_string(in_v_ax.res+1) +").") ;

      if (v_idx % 5 == 0 )
      {
        char tmp_term[150] ;
        snprintf(tmp_term, sizeof(tmp_term), "Mixed sequence (%3.lu, %3.lu) out "
          "of (%3.lu, %3.lu).\r", d_idx+1, v_idx+1, in_d_ax.res+1, in_v_ax.res+1) ;
        std::cout << tmp_term << std::flush ;
      }

      init_press = p_of_e(in_v_ax[v_idx])  ;

      double r = r_min ;
      double y[4], y_mantle[2] ; 
    
      y[2] = init_press ;
      y[3] = init_press_dark ;
      y[0] = (4. / 3.) * M_PI * pow (r, 3.) * GetEDens(y[2]) ;
      y[1] = (4. / 3.) * M_PI * pow (r, 3.) * GetEDens_Dark(y[3]) ;

      // ------------------------------------------
      //                RADIUS LOOP
      // ------------------------------------------
      RadiusLoopDark(r, y, y_mantle) ;
      // ------------------------------------------


      // ----------------------------------------------
      mixed_star.SurfaceIsReached() ;
      // ----------------------------------------------
#if DarkCore_Analysis
      double m_n = Zaki::Physics::NEUTRON_M_FM ;
      double m_chi = 0.8*Zaki::Physics::NEUTRON_M_FM ;

      double critical_rho_val = pow(m_n - m_chi*m_chi/m_n, 3)
                                / ( 24 * M_PI*M_PI ) ;
      int critical_idx = mixed_star.ds_dar[
                              mixed_star.rho_idx 
                              ].GetIdx(critical_rho_val) ;
      std::cout << "Critical Rho = " << critical_rho_val 
                << ", Index = " << critical_idx 
                << ", Size = " 
                << mixed_star.ds_dar.Dim()[0] 
                << ", rho_d [v_idx-1] = " << mixed_star.ds_dar[
                              mixed_star.rho_idx 
                              ][critical_idx-1]
                << ", rho_d [v_idx] = " << mixed_star.ds_dar[
                              mixed_star.rho_idx 
                              ][critical_idx]
                << "\n\t r [v_idx-1] = " << mixed_star.ds_dar[
                              mixed_star.r_idx 
                              ][critical_idx-1] 
                << ", r [v_idx] = " << mixed_star.ds_dar[
                              mixed_star.r_idx 
                              ][critical_idx] << "\n" ;

    // std::map<std::string, std::string> baryons = {
    //   {"10", "Neutron", 6},
    //   {"100", "Lambda", 7}, {"110", "Sigma-", 10}, {"111", "Sigma0"},
    //   {"112", "Sigma+"}, {"120","Xi-"}, {"121", "Xi0"}
    // } ;

      Zaki::Vector::DataColumn b_den = mixed_star.ds_vis[
                                        mixed_star.rho_i_v_idx[7] 
                                        ] * mixed_star.ds_vis[
                                        mixed_star.rho_idx 
                                        ];
      Zaki::Vector::DataColumn r_set = mixed_star.ds_vis[
                                        mixed_star.r_idx 
                                        ];
      Zaki::Vector::DataColumn M_r = mixed_star.ds_vis[
                                        mixed_star.m_idx 
                                        ];
      Zaki::Vector::DataColumn nu_r = mixed_star.ds_vis[
                                        mixed_star.nu_idx 
                                        ];

      Zaki::Vector::DataColumn integ_fr = 4*M_PI*b_den ;
      integ_fr *= r_set.pow(2) * (1. - 2*M_r / r_set ).pow(-0.5) ;

      Zaki::Vector::DataColumn integ_b_dot = integ_fr ;
      integ_b_dot *= exp( nu_r ) ;

      Zaki::Vector::DataSet integrand({r_set, integ_fr, integ_b_dot}) ;

      std::cout << "\n integ_fr.Size() = " << integ_fr.Size() ;
      std::cout << "\n integ_b_dot.Size() = " << integ_b_dot.Size() ;

      integrand.Interpolate(0, {1, 2}) ;
      double fr_result = integrand.Integrate(1, {r_set[0], r_set[-1]}) ;
      double b_dot_result = integrand.Integrate(2, {r_set[0], r_set[-1]}) ;
      
      double fr_result_choked = integrand.Integrate(1, {r_set[critical_idx], r_set[-1]}) ;
      double b_dot_result_choked = integrand.Integrate(2, {r_set[critical_idx], r_set[-1]}) ;

      std::cout << "\nFraction = " << fr_result* 1e54 / mixed_star.sequence.v.b 
                << ", [B_dot/B] = " 
                << b_dot_result* 1e54 / mixed_star.sequence.v.b  
                <<"\n Fraction choked = "
                <<  fr_result_choked* 1e54 / mixed_star.sequence.v.b  
                << ", [B_dot/B] choked = "
                << b_dot_result_choked* 1e54 / mixed_star.sequence.v.b  << "\n\n" ; 
#endif
      // ----------------------------------------------
      // Saving the results
      // ----------------------------------------------
      // mixed_star.SetWrkDir( wrk_dir ) ;
      if(exp_cond_f)
      {
        if(exp_cond_f(mixed_star))
          ExportMixedStarProfile(v_idx, d_idx, in_dir + in_file) ;
          // mixed_star.Export(in_dir + "/Mixed_" + 
          //     std::to_string(d_idx) + "_" +
          //     std::to_string(v_idx) + ".tsv") ;
      }
      mixed_sequence.Add(mixed_star) ;
      // ----------------------------------------------
      
#if 0
      // ----------------------
      sequence.emplace_back(
                      GetEDens(init_press), 
                      tmp_m_star->mass[-1]/Zaki::Physics::SUN_M_KM,
                      tmp_m_star->radius[-1], init_press, 
                      tmp_m_star->GetBaryonNum(), 
                      tmp_m_star->GetMomInertia()
                            );
      // ----------------------------------------------

      std::vector<std::string> tmp_fname_v = Zaki::String::Pars(in_file.Str(), "*") ;
      std::string tmp_fname ;
      if(tmp_fname_v.size() <2 )
      {
        Z_LOG_WARNING("File name pattern doesn't match '[]*[]'.") ;
        tmp_fname = tmp_fname_v[0] + std::to_string(v_idx) ;
      }
      else
      {
        tmp_fname = tmp_fname_v[0] + std::to_string(v_idx) + tmp_fname_v[1] ;
      }
      
    // ............ Creating a directory ............
    if (mkdir((wrk_dir + in_dir + "/"
                      + std::to_string(d_idx)
                      ).Str().c_str(), ACCESSPERMS) == -1) 
    {
      Z_LOG_INFO("Directory '"+ (wrk_dir + in_dir + "/" 
                                      + std::to_string(d_idx)
                                      ).Str()
                                      +"' wasn't created, because: "
                                      +strerror(errno)+".") ;    
    }
    else
      Z_LOG_INFO("Directory '"+(wrk_dir + in_dir + "/" 
                                      + std::to_string(d_idx) 
                                      ).Str()
                                      +"' created.") ; 
    // .................................................

      Zaki::File::VecSaver vec_saver(wrk_dir + in_dir + "/" 
                                      + std::to_string(d_idx)
                                      + "/" + tmp_fname) ;

      // -------------------- Visible Header Begins --------------------
      char res_header[200] ;
      snprintf(res_header, sizeof(res_header), "%-14s\t %-14s\t %-14s\t %-15s\t %-14s\t %-14s\t %-14s", 
              "r(km)", "m", "nu'(1/cm)", " nu", "p(dyne/cm^2)", "e(g/cm^3)", "rho(1/fm^3)" ) ;

      std::string tmp_label = res_header ;
      for(auto& lab : eos_tab.extra_labels)
      {
        snprintf(res_header, sizeof(res_header), "\t %-14s", Zaki::String::Strip(lab, ' ').c_str() ) ;
        tmp_label += res_header ;
      }
      // -------------------- Visible Header Ends--------------------

      vec_saver.SetHeader(tmp_label.c_str()) ;
      vec_saver.Export1D(results) ;
      
      // -------------------- Dark Header Begins --------------------
      snprintf(res_header, sizeof(res_header), "%-14s\t %-14s\t %-14s\t %-15s\t %-14s\t %-14s\t %-14s", 
              "r(km)", "m", "nu'(1/cm)", " nu", "p(dyne/cm^2)", "e(g/cm^3)", "rho(1/fm^3)" ) ;

      tmp_label = res_header ;
      for(auto& lab : eos_tab_dark.extra_labels)
      {
        snprintf(res_header, sizeof(res_header), "\t %-14s", Zaki::String::Strip(lab, ' ').c_str() ) ;
        tmp_label += res_header ;
      }
      // -------------------- Dark Header Ends--------------------
      vec_saver.SetHeader(tmp_label.c_str()) ;

      tmp_fname_v = Zaki::String::Pars(tmp_fname, ".") ;
      if(tmp_fname_v.size() == 2)
        vec_saver.SetFileName(wrk_dir + in_dir +  "/" 
                                      + std::to_string(d_idx)
                                      + "/" + tmp_fname_v[0] 
                                      + "_dark." +  tmp_fname_v[1] ) ;
      else
      {
        vec_saver.SetFileName(wrk_dir + in_dir + "/" 
                                      + std::to_string(d_idx)
                                      + "/" + tmp_fname + "_dark" ) ;
      }
      
      vec_saver.Export1D(results_dark) ;


#endif

      mixed_star.Reset() ;

#if TOV_SOLVER_VERBOSE
      printf("\n========================== %s (%lu) "
            "==========================\n\n", "End of Seq.", v_idx+1) ;
#endif
    }
    // ----------------------------------------------------------------
    //                  TOV Visible sequence loop ends!
    // ----------------------------------------------------------------
#if TOV_SOLVER_VERBOSE
    std::cout <<"\n\t\t *************************"
              <<" TOV Solver Finished *************************"<<"\n\n" ;
#endif

  }
  // ----------------------------------------------------------------
  //                  TOV Dark sequence loop ends!
  // ----------------------------------------------------------------
  ExportMixedSequence(in_dir + in_file + "_Sequence.tsv") ;
  // mixed_sequence.Export(in_dir + "/Mixed_Sequence.tsv") ;
}

//--------------------------------------------------------------
void TOVSolver::ExportMixedStarProfile
  ( const size_t& v_idx, const size_t& d_idx, 
    const Zaki::String::Directory& in_dir)
{
  mixed_star.Export(in_dir + 
              std::to_string(d_idx) + "_" +
              std::to_string(v_idx) + ".tsv") ;
}

//--------------------------------------------------------------
void TOVSolver::ExportMixedSequence
  (const Zaki::String::Directory& in_dir)
{
  mixed_sequence.Export(in_dir) ;
}

//--------------------------------------------------------------
// The radius iteration in the mixed star scenario
void TOVSolver::RadiusLoopDark(double& in_r, double* in_y,
                                  double* in_y_mantle)
{
  PROFILE_FUNCTION() ;

  //----------------------------------------
  //          GSL ODE SYSTEM SETUP
  //----------------------------------------
  gsl_odeiv2_system ode_sys_core = 
    {CompactStar::TOVSolver::ODE_Dark_Core, nullptr, 4, this} ;

  gsl_odeiv2_driver *tmp_driver_core = gsl_odeiv2_driver_alloc_y_new
    (&ode_sys_core, gsl_odeiv2_step_rk8pd,
    1.e-1, 1.e-10, 1.e-10);

  gsl_odeiv2_system ode_sys_mantle = 
    {CompactStar::TOVSolver::ODE_Dark_Mantle, nullptr, 2, this} ;

  gsl_odeiv2_driver *tmp_driver_mantle = gsl_odeiv2_driver_alloc_y_new
    (&ode_sys_mantle, gsl_odeiv2_step_rk8pd,
    1.e-1, 1.e-10, 1.e-10);
  //----------------------------------------

  double min_log_r = log10(r_min) ;
  double max_log_r = log10(r_max) ;

  // double min_log_r = log(r_min) ;
  // double max_log_r = log(r_max) ;

  double step = (max_log_r - min_log_r) / radial_res ;
  double step_scale = 1 ; // Adaptive steps (Aug 6, 2020)
  // mixed_star.Reserve(radial_res) ;
  bool CORE_REGION = true ;

  for (double log_r_i = min_log_r;  log_r_i <= max_log_r; log_r_i += step*step_scale)
  {
    double ri = pow(10, log_r_i) ;
    // double ri = exp(log_r_i) ;
//    double ri = log_r_i ;

    // double tmp_delta_p = 0 ;     // Adaptive steps (Aug 6, 2020)

    int status = -1 ;

    if(CORE_REGION)
    {
      // tmp_delta_p = in_y[2] ;
      status = gsl_odeiv2_driver_apply (tmp_driver_core, &in_r, ri, in_y);
      // tmp_delta_p = (tmp_delta_p - in_y[2]) / in_y[2] ;
    }
    else
    {
      // tmp_delta_p = in_y_mantle[1] ;

      status = gsl_odeiv2_driver_apply (tmp_driver_mantle, &in_r, ri, in_y_mantle);

      // tmp_delta_p = (tmp_delta_p - in_y_mantle[1]) / in_y_mantle[1] ;
    }

    if (status != GSL_SUCCESS)
    {
#if TOV_SOLVER_VERBOSE
      printf("\t-------------------%s-------------------\n", "GSL") ;
      printf ("\t GSL error, return value=%d\n", status);
#endif
      if(CORE_REGION)
      {
#if TOV_SOLVER_VERBOSE
        printf("\t Visible Pressure = %2.2e.\n", in_y[2]) ;
        printf("\t Dark Pressure    = %2.2e.\n", in_y[3]) ;
#endif
        if(dark_core)
        { 
          m_core = in_y[1] ;
          
          // Initial condition for ODE_Mantle
          in_y_mantle[0] = in_y[0] ;
          in_y_mantle[1] = in_y[2] ;
        }
        else
        {
          m_core = in_y[0] ;
          
          // Initial condition for ODE_Mantle
          in_y_mantle[0] = in_y[1] ;
          in_y_mantle[1] = in_y[3] ;
#if TOV_SOLVER_VERBOSE
          std::cout << "\n\t m_core = "<<in_y[0]/GSL_CONST_CGSM_SOLAR_MASS<<"\n" ;
          std::cout << "\t in_y_mantle[0] = "<< in_y_mantle[0] << "\n";
          std::cout << "\t in_y_mantle[1] = "<< in_y_mantle[1] << "\n";
#endif
        }

        CORE_REGION = false ;
        gsl_odeiv2_driver_free (tmp_driver_core) ;
      }
      else // Mantle's surface reached!
      {
#if TOV_SOLVER_VERBOSE
        printf("\t  Surface Pressure = %2.2e.\n", in_y_mantle[1]) ;
        printf("\t-------------------%s-------------------\n", "GSL") ;
#endif
        break ;
      }
#if TOV_SOLVER_VERBOSE
      printf("\t-------------------%s-------------------\n", "GSL") ;
#endif

      // Experimental : NOT SURE !!!!!!!!!!!!!!!!!!!
      continue ; // Jump over the boundary to avoid duplicate values
    }

    if (CORE_REGION)
    {
      mixed_star.Append_Core(
        { in_r/1.e+5, in_y[0]/GSL_CONST_CGSM_SOLAR_MASS, 
          GetNuDer_Dark(in_r, {in_y[0] + in_y[1], in_y[2] + in_y[3]}), 0,
          in_y[2], GetEDens(in_y[2]), 
          GetRho(in_y[2]), GetRho_i(in_y[2])
        },
        { in_r/1.e+5, in_y[1]/GSL_CONST_CGSM_SOLAR_MASS, 
          GetNuDer_Dark(in_r, {in_y[0] + in_y[1], in_y[2] + in_y[3]}), 0,
          in_y[3], GetEDens_Dark(in_y[3]), 
          GetRho_Dark(in_y[3]), GetRho_i_Dark(in_y[3])
        }
                              ) ;

    }
    else if (dark_core) // dark core with a visible mantle
    {
      mixed_star.Append_Visible_Mantle(
        { in_r/1.e+5, in_y_mantle[0]/GSL_CONST_CGSM_SOLAR_MASS, 
          GetNuDer_Dark(in_r, {in_y_mantle[0] + m_core, in_y_mantle[1]}), 0,
          in_y_mantle[1], GetEDens(in_y_mantle[1]), 
          GetRho(in_y_mantle[1]), GetRho_i(in_y_mantle[1])
        }
                                        ) ;
    }
    else // visible core, with a dark mantle
    {
      mixed_star.Append_Dark_Mantle(
        { in_r/1.e+5, in_y_mantle[0]/GSL_CONST_CGSM_SOLAR_MASS, 
          GetNuDer_Dark(in_r, {in_y_mantle[0] + m_core, in_y_mantle[1]}), 0,
          in_y_mantle[1], GetEDens_Dark(in_y_mantle[1]), 
          GetRho_Dark(in_y_mantle[1]), GetRho_i_Dark(in_y_mantle[1])
        }
                                    ) ;
    }
    
  }

  gsl_odeiv2_driver_free (tmp_driver_mantle) ;

}

//--------------------------------------------------------------
// Exports the star sequence
void TOVSolver::ExportSequence(const Zaki::String::Directory& in_dir) const 
{
  Zaki::File::VecSaver vec_saver_2(wrk_dir + "/" + in_dir);

  char seq_header[200] ;
  snprintf(seq_header, sizeof(seq_header), "%-14s\t %-14s\t %-14s\t %-14s\t %-14s\t %-14s", 
          "ec(g/cm^3)", "M",  "R(km)", "pc(dyne/cm^2)", "B", "I(km^3)" ) ;
  vec_saver_2.SetHeader(seq_header) ;
  vec_saver_2.Export1D(sequence) ;
}
//--------------------------------------------------------------
// Exports the mixed star sequence
// void TOVSolver::ExportMixedSequence(const Zaki::String::Directory& in_dir) const 
// {
//   Zaki::File::VecSaver vec_saver_2(wrk_dir + "/" + in_dir);

//   char seq_header[400] ;
//   sprintf(seq_header, "%-14s\t %-14s\t %-14s\t %-14s\t %-14s\t %-14s"
//                   "\t %-14s\t %-14s\t %-14s\t %-14s\t %-14s\t %-14s", 
//           "ec(g/cm^3)", "M",  "R(km)", "pc(dyne/cm^2)", "B", "I(km^3)",
//           "ec_d(g/cm^3)", "M_d",  "R_d(km)", "pc_d(dyne/cm^2)", "B_d", "I_d(km^3)" ) ;

//   vec_saver_2.SetHeader(seq_header) ;

//   vec_saver_2.Export1D(mixed_sequence) ;
// }
//--------------------------------------------------------------

//==============================================================
