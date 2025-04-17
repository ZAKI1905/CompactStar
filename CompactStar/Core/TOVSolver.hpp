// -*- lsst-c++ -*-
/*
* CompactStar
* See License file at the top of the source tree.
*
* Copyright (c) 2023 Mohammadreza Zakeri
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*/

/**
 * @file TOVSolver.hpp
 *
 * @brief Solves TOV equations.
 *
 * @ingroup Core
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/
// Last edit on Dec 15 2020
#ifndef CompactStar_TOVSolver_H
#define CompactStar_TOVSolver_H

// #include <thread>
#include <vector>
#include <gsl/gsl_spline.h>

#include <Zaki/Math/Math_Core.hpp>

#include "CompactStar/Core/Prog.hpp"
#include "CompactStar/Core/MixedStar.hpp"


//==============================================================
namespace CompactStar
{

class Analysis ;
//==============================================================
class eps_pair : public Zaki::Math::Coord2D
{
  // double e_v, e_d ;
  public:
  eps_pair(const double& in_e_v, const double& in_e_d)
        : Coord2D(in_e_v, in_e_d) {}
  
  eps_pair(const Coord2D& in_c)
        : Coord2D(in_c.x, in_c.y) {}

  double e_v() const
  {
    return x ;
  }

  double e_d() const
  {
    return y ;
  }
};
//==============================================================


//==============================================================

class Contour 
// : public Zaki::Math::Curve2D<eps_pair>
{
  public:
    Zaki::Math::Curve2D curve ;
    double val = 0 ;
    double precision = 1e-8 ;
    size_t max_steps = 35 ;
  
  Contour() {} 
  Contour(const std::string& in_label) : curve(in_label) {} 
  // Contour(const Zaki::Math::Curve2D<Zaki::Math::Coord2D>& in_curve) 
  // {
  //   for (size_t i = 0; i < in_curve.Size() ; i++)
  //   {
  //     pts.emplace_back(in_curve.pts[i].x, in_curve.pts[i].y) ;
  //   }
    
  //   SetLabel(in_curve.label) ;
  //   val = std::atof(in_curve.label.c_str()) ;
  // } 

  // Contour(const Zaki::Math::Curve2D<eps_pair>& in_curve) 
  // {
  //   for (size_t i = 0; i < in_curve.Size() ; i++)
  //   {
  //     pts.emplace_back(in_curve.pts[i].x, in_curve.pts[i].y) ;
  //   }
    
  //   SetLabel(in_curve.label) ;
  //   val = std::atof(in_curve.label.c_str()) ;
  // } 
  // std::vector<eps_pair> guide ;

  // Contour(const Zaki::String::Directory& in_file) ;

  size_t Size() const 
  {
    return curve.Size() ;
  }

  void Import(const Zaki::String::Directory& in_file) 
  {
    curve.Import(in_file) ;
    val = std::atof(curve.GetLabel().c_str()) ;
  }

};

//==============================================================
struct EOSTable
{
  private:
    std::string eps_label ;
    std::string pre_label ;
    std::string rho_label ;

  public:
    std::vector<double> eps ;
    std::vector<double> pre ;
    std::vector<double> rho ;

    std::vector<std::vector<double> > rho_i ;
    std::vector<std::string> extra_labels   ;

    size_t Size()
    {
      return eps.size() ;
    }

    void SetLabels( const std::string& in_eps_label, 
                    const std::string& in_pre_label,
                    const std::string& in_rho_label)
    {
      eps_label = Zaki::String::Strip(in_eps_label, ' ') ;
      pre_label = Zaki::String::Strip(in_pre_label, ' ') ;
      rho_label = Zaki::String::Strip(in_rho_label, ' ') ;
    }

    void AddExtraLabels(const std::string& in_label)
    {
      extra_labels.emplace_back(in_label) ;
    }

    void Print() const
    {
      std::cout << " *-------------------------------------------* " << "\n" ;
      std::cout <<  " | "<<eps_label<<"   | "<<pre_label
                <<")   | "<<rho_label<<"  |\n" ;
      std::cout << " *-------------------------------------------* " << "\n" ;
      std::string tmp_str;
      for (size_t i = 0; i < eps.size() ; i++)
      {
        std::cout << " | " << eps[i] << "\t      | " 
                  << pre[i] << "\t      | " << rho[i] <<"\n" ;
      }
      std::cout << " *--------------------------------------------* " << "\n" ;
    }
};
//==============================================================
//            nu_der added on December 15, 2020
// Struct representing TOV solution points
struct TOVPoint
{
  double r, m, nu_der, nu, p, e, rho ;
  std::vector<double> rho_i ;

  // TOVPoint(const double& in_r, const double& in_m,
  //           const double& in_p, const double& in_e)
  //           : r(in_r), m(in_m), p(in_p), e(in_e), rho(-1)
  //         {}

  TOVPoint(const double& in_r, const double& in_m, 
            const double& in_nu_der, const double& in_nu,
            const double& in_p, const double& in_e,
            const double& in_rho, const std::vector<double>& in_rho_i)
            : r(in_r), m(in_m), nu_der(in_nu_der), nu(in_nu),
              p(in_p), e(in_e), rho(in_rho),
              rho_i(in_rho_i)
          {}

  std::string Str() const
  {
    std::stringstream ss;
    char tmp[200] ;
    snprintf(tmp, sizeof(tmp), "%.8e\t %.8e\t %.8e\t %.8e\t %.8e\t %.8e\t %.8e", 
                     r,     m,  nu_der,  nu,     p,     e,   rho) ;
    ss << tmp ;
    for (double tmp_rho : rho_i)
    {
      snprintf(tmp, sizeof(tmp), "\t %.8e", tmp_rho) ;
      ss << tmp ;
    }
    
    return ss.str() ;
  }
};
//==============================================================
//            nu_der added on December 15, 2020
// Struct representing TOV solution points
struct TOV_Nu_Point
{
  double r, m, nu ;

  TOV_Nu_Point(const double& in_r, const double& in_m, const double& in_nu)
            : r(in_r), m(in_m), nu(in_nu)
          {}

  std::string Str() const
  {
    std::stringstream ss;
    char tmp[150] ;
    snprintf(tmp, sizeof(tmp), "%.8e\t %.8e\t %.8e", r, m, nu) ;
    ss << tmp ;
    
    return ss.str() ;
  }
};

//==============================================================
//                        Sequence class
//==============================================================
class Sequence : public Prog
{
private:
  std::vector<SeqPoint> seq ;
  // std::mutex m_mutex ;

public:
  Sequence() ;
  // Sequence(const std::string& in_name ) 
  // : Prog(in_name) {}

  ~Sequence() 
  {
    // std::cout << "[ Thread = " << std::this_thread::get_id()
    //           << "] Destructor called for: " 
    //           << name << "\n" ;
  }

  void Add(const NStar& in_star) ;

  // Exports the star sequence
  void Export(const Zaki::String::Directory& in_dir="") const ;

  // Combines two sequences
  void Combine(const Sequence& other) ;

  // Clears the sequence
  void Clear() ;
};

//==============================================================
//                        MixedSequence class
//==============================================================
class MixedSequence : public Prog
{
private:
  std::vector<MixedSeqPoint> seq ;
  // std::mutex m_mutex ;

public:
  MixedSequence() ;
  // MixedSequence(const std::string& in_name ) 
  // : Prog(in_name) {}

  ~MixedSequence() 
  {
    // std::cout << "[ Thread = " << std::this_thread::get_id()
    //           << "] Destructor called for: " 
    //           << name << "\n" ;
  }

  void Add(const MixedStar& in_star) ;

  // Exports the mixed star sequence
  void Export(const Zaki::String::Directory& in_dir="") const ;

  // Combines two sequences
  void Combine(const MixedSequence& other) ;

  // Clears the sequence
  void Clear() ;
};

//==============================================================
class TOVSolver : public Prog
{
  friend class MixedStar ;
  friend class NStar ;
  //--------------------------------------------------------------
  protected:
    
    /// The number of TOV object instances
    // static inline std::atomic<size_t> tov_counter = 0;

    // std::mutex m_mutex ;

    /// The radial resolution for the solver
    size_t radial_res = 10000 ;

    EOSTable eos_tab ;
    EOSTable eos_tab_dark ;
    
    MixedStar mixed_star ;
    NStar n_star ;

    // This workspace stores state variables for interpolation lookups.
    // It caches the previous value of an index lookup. 
    // When the subsequent interpolation point falls in the same 
    // interval its index value can be returned immediately.
    // There has to be separate accelerators for variables
    // with different domains.
    //
    // EOS splines ( Domain = pressure )
    gsl_interp_accel *visi_p_accel  = nullptr ;
    gsl_interp_accel *dark_p_accel  = nullptr ;
    gsl_interp_accel *mixed_r_accel = nullptr ;

    const gsl_interp_type* TOV_gsl_interp_type = gsl_interp_steffen ;
    // const gsl_interp_type* TOV_gsl_interp_type = gsl_interp_linear ;

    // gsl_interp_accel *accel = nullptr ;
    // gsl_interp_accel *dark_accel = nullptr ;


    // We use cubic spline
    // Cubic spline with natural boundary conditions. 
    // The resulting curve is piecewise cubic on each interval, 
    // with matching first and second derivatives at the supplied data-points.
    //  The second derivative is chosen to be zero at the first point and last point.
    gsl_spline *visi_eps_p_spline = nullptr ;
    gsl_spline *dark_eps_p_spline = nullptr ;
    
    /// Total baryon number density (function of p)
    gsl_spline *visi_rho_p_spline = nullptr ;
    gsl_spline *dark_rho_p_spline = nullptr ;

    /// Total baryon number density (function of r)
    // gsl_spline *visi_rho_r_spline = nullptr ;
    // gsl_spline *dark_rho_r_spline = nullptr ;

    /// This is for the extra rho's (individual species)
    std::vector<gsl_spline*> visi_rho_i_p_spline ;
    std::vector<gsl_spline*> dark_rho_i_p_spline ;

    // Added on December 15, 2020
    /// nu(r) derivative spline
    gsl_spline *nu_der_r_spline = nullptr ;


    /// Initial pressure (at r = r(i=0))
    double init_press = -1 ;
    double init_press_dark= -1 ;

    /// Initial energy density (at r = r(i=0))
    double init_edens = -1 ;
    double init_edens_dark = -1 ;

    /// Solution to TOV are saved here
    // std::vector<TOVPoint> results ;
    // std::vector<TOVPoint> results_dark ;

    // std::vector<SeqPoint> sequence ;
    Sequence sequence ;
    // std::vector<MixedSeqPoint> mixed_sequence ;
    MixedSequence mixed_sequence ;

    // static inline MixedSequence mixed_seq_static ;

    /// Returns the pressure corresponding to in_e
    /// It's the inverse function of "GetEDens"
    double p_of_e(const double& in_e) ;
    double p_of_e_dark(const double& in_e) ;

    /// Cost function for finding the pressure
    /// for "p_of_e" method.
    double cost_p_of_e(const double in_e) ;
    double cost_p_of_e_dark(const double in_e) ;
    double cost_p_of_e_input = 0 ;
    double cost_p_of_e_input_dark = 0 ;

    /// The value of pressure cut-off is the pressure 
    /// at the surface of the star
    /// Theoretically it's zero, but we choose 
    /// 1e-5 times smaller than the central pressure
    /// or the lowest possible pressure given by the Eq. 
    /// of state.
    double PressureCutoff() const ;
    double PressureCutoff_Dark() const ;

    /// Minimum radius in the solver (cm)
    double r_min = 1. ;

    /// Maximum radius in the solver (cm)
    // double r_max = 20e+5 ; // 20 km
    // double r_max = 50e+5 ; // 50 km
    double r_max = 70e+5 ; // 70 km

    // dark core or visible core :
    bool dark_core = true ;

    double m_core = -1 ;
    enum class TOVSolverStatus
    {
      Vis_Surf_Reached = 100, Dark_Surf_Reached = 101, Mantle_Surf_Reached = 102
    };

    /// Sets the work directory for the member objects
    Prog* SetMemWrkDir(const Zaki::String::Directory&) override ;

    /// Pointer to the condition for exporting
    ///  the mixed star profile
    bool (*mix_exp_cond_f)(const MixedStar&) = nullptr ;

    /// Pointer to the condition for exporting 
    /// the star profile
    bool (*n_exp_cond_f)(const NStar&) = nullptr ;

    /// Pointer to the analysis object
    Analysis* analysis = nullptr ;

    /// Exclusion region in (ec_v, ec_d) space
    Zaki::Math::Cond_Polygon c_poly ;

    /// Number of points that are ignored 
    ///  because they are in the exclusion region.
    size_t ignored_counter = 0 ;

    void Hidden_ImportEOS_Vis(const Zaki::String::Directory&) ;
    void Hidden_ImportEOS_Dar(const Zaki::String::Directory&) ;

    // The precision in printing the profiles
    int profile_precision = 9 ;

    /// The precision in evaluation pressure as a function of density 
    double p_of_e_prec = 1e-4 ;

  //--------------------------------------------------------------
  public:
    
    TOVSolver();
    ~TOVSolver() ;
    
    /// Copy Constructor 
    TOVSolver(const TOVSolver &) = delete ;

    /// Assignment operator
    TOVSolver& operator= (const TOVSolver&) = delete ;

    void ImportEOS(const Zaki::String::Directory&) ;
    void ImportEOS(const Zaki::String::Directory& vis_eos, 
                  const Zaki::String::Directory& dar_eos) ;

    // void ImportEOS_Visible(const Zaki::String::Directory&) ;
    // void ImportEOS_Dark(const Zaki::String::Directory&) ;
    void PrintEOSTable() const ;

    /// Returns the energy density given pressure
    double GetEDens(const double&) ;
    double GetEDens_Dark(const double&) ;

    /// Returns the total baryon number density given pressure
    double GetRho(const double&) ;
    double GetRho_Dark(const double&) ;


    //........... NOV 3, 2021 Begins ...........
    /// Returns the total baryon number density given radius
    double GetRho_r(const double&) ;
    //........... NOV 3, 2021 Ends  ............

    /// Returns the specific number density given pressure
    std::vector<double> GetRho_i(const double& in_p) ;
    std::vector<double> GetRho_i_Dark(const double& in_p) ;

    double GetInitPress() const ;
    double GetInitEDens() const ;
    double GetInitPress_Dark() const ;

    //               Added on December 15, 2020
    /// Returns the derivative of the metric nu(r) function
    double GetNuDer(const double r, const std::vector<double>& y) ;
    double GetNuDer_Dark(const double r, const std::vector<double>& y) ;

    //            Added on December 15, 2020
    // Returns the nu_der value given the radius input
    double GetNuDerSpline(const double& in_r) ;

    /// Adds the condition for printing the mixed star profile
    void AddMixCondition(bool (*func)(const MixedStar&)) ;

    /// Adds the condition for printing the star profile
    void AddNCondition(bool (*func)(const NStar&)) ;

    /// Attaches a pointer to analysis 
    void AddAnalysis(Analysis*) ;

    /// Input a range of initial pressure
    void Solve( const Zaki::Math::Axis&, 
                const Zaki::String::Directory&, 
                const Zaki::String::Directory& file_name) ;

    void Solve_Mixed(  const Zaki::Math::Axis& vis_ax, 
                      const Zaki::Math::Axis& dark_ax, 
                      const Zaki::String::Directory& dir,
                      const Zaki::String::Directory& file_name) ;

    void Solve_Mixed( const Contour& eps_cont, 
                      const Zaki::String::Directory& dir,
                      const Zaki::String::Directory& file_name) ;

    // The radius iteration in the neutron star scenario
    void RadiusLoop(double& r, double* y) ;

    // The radius iteration in the mixed star scenario
    void RadiusLoopMixed(double& r, double* y_core,
                        double* y_mantle) ;

    /// Exports the star sequence
    void ExportSequence(const Zaki::String::Directory&) const ;
    
    /// Exports the mixed sequence
    virtual void ExportMixedSequence(const Zaki::String::Directory&) ;

    /// Exports the mixed star profile
    virtual void ExportMixedStarProfile(const size_t& v_idx, 
      const size_t& d_idx, const Zaki::String::Directory&) ;

    /// Exports the neutron star profile
    virtual void ExportNStarProfile(const size_t& idx, 
                            const Zaki::String::Directory&) ;

    virtual void SurfaceIsReached(const size_t& v_idx,
                                  const size_t& d_idx) ;

    void SurfaceIsReached() ;


    // For mixed stars
    virtual void PrintStatus(const size_t& v_idx,
                             const size_t& d_idx,
                             const size_t& v_res, 
                             const size_t& d_res) ;
    // For neutron stars
    void PrintStatus(const size_t& in_idx, const size_t& in_res) ;


    //            Added on December 15, 2020
    /// Evaluates & exports the nu(r) function table
    // void ExportNu(const Zaki::String::Directory& in_dir) ;


    static int ODE(double r, const double y[], double f[], void *params) ;
    static int ODE_Dark_Core(double r, const double y[], double f[], void *params) ;
    static int ODE_Dark_Mantle(double r, const double y[], double f[], void *params) ;

    void SetExclusionRegion(const Zaki::Math::Cond_Polygon&) ;
    void SetRadialRes(const size_t&) ;

    /// Sets the printing precision for the NStar profiles
    void SetProfilePrecision(const int& prec) ;

    /// Sets maximum value for radius in the solver
    /// if a maximum is known, say 15 km for NS,
    /// it would increase the radial resolution
    /// defined by:
    ///  delta R = scale(R) * (R_max - R_min) / radial_res
    void SetMaxRadius(const double&) ; 

    /// Empties the sequence
    void ClearSequence() ;

    /// @brief It generates a sequence of NS by varying radial resolution
    /// @param e_c central energy density
    /// @param dir directory for saving the results
    /// @param file result file name
    void GenTestSequence(const double& e_c, 
                  const Zaki::String::Directory& dir,
                  const Zaki::String::Directory& file) ;
};

//==============================================================
} // CompactStar namespace
//==============================================================
// std::ostream& operator << ( std::ostream &, const Coord2D&);
// std::ostream& operator << ( std::ostream &, const Segment&);
//==============================================================
#endif /*CompactStar_TOVSolver_H*/
