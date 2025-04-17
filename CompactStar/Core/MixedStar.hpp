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
 * @file MixedStar.hpp
 *
 * @brief Dark+Visible admixed star.
 *
 * @ingroup Core
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/

#ifndef CompactStar_MixedStar_H
#define CompactStar_MixedStar_H

// #include <gsl/gsl_spline.h>
// #include <gsl/gsl_integration.h>

// #include <Zaki/Math/GSLFuncWrapper.hpp>
#include <Zaki/String/Directory.hpp>
#include <Zaki/Vector/DataSet.hpp>

#include "CompactStar/Core/Prog.hpp"
#include "CompactStar/Core/NStar.hpp"
// #include "CompactStar/Core/RotationSolver.hpp"

//==============================================================
// struct gsl_integration_workspace ;
//==============================================================
namespace CompactStar
{

// forward declaration 
struct TOVPoint ;     
class TOVSolver ;

//==============================================================
//             Mixed Star Sequence Points Class
//==============================================================
class MixedSeqPoint
{
  public:
    
    SeqPoint  v, d ;
    size_t v_idx = 0 ;
    size_t d_idx = 0 ;
    // double ec_d, ec, m, m_d, r_d, r, pc_d, pc, b, I ;

    MixedSeqPoint()
      { } ;

    MixedSeqPoint(const size_t& in_v_idx, 
                  const SeqPoint& in_v,
                  const size_t& in_d_idx, 
                  const SeqPoint& in_d)
            : v(in_v), d(in_d), v_idx(in_v_idx), 
              d_idx(in_d_idx)
      { }

    MixedSeqPoint(const size_t& in_v_idx, 
                  const size_t& in_d_idx, 
                  const std::vector<double>& in_seq_row)
                  : v_idx(in_v_idx), d_idx(in_d_idx)
      {
        if ( in_seq_row.size() != 12 )
          Z_LOG_ERROR("The mixed sequence data point is incomplete.") ;
        else
          {
            v  = { in_seq_row[0], in_seq_row[1], in_seq_row[2],
                  in_seq_row[3], in_seq_row[4], in_seq_row[5] } ;
            d  = { in_seq_row[6], in_seq_row[7], in_seq_row[8],
                  in_seq_row[9], in_seq_row[10], in_seq_row[11] } ;
          }
      }

    std::string Str() const
    {
      // May 10, 2022: For some reason 
      // %lu and %zu specifiers don't 
      // print zeros of v_idx and d_idx!
      // Fixed: Removed "." in "%5.lu".
      std::stringstream ss ;
      char tmp[75] ;
      snprintf(tmp, sizeof(tmp), "%5lu\t %5lu\t ", v_idx, d_idx ) ;
      // std::cout << "\n\t\t" << tmp  << "\n" ;
      ss << tmp << v.Str() << "\t " << d.Str() ;
      return ss.str() ;
    }

    void Reset() 
    {
      v_idx = 0 ;
      d_idx = 0 ;
      v.Reset() ;
      d.Reset() ;
    }

    //..................................................
    // Operator overloading
    // Addition (+)
    MixedSeqPoint operator+(const MixedSeqPoint& seq) const
    {
      MixedSeqPoint out_seq(0, v+seq.v, 0, d+seq.d) ;
      return out_seq ;
    }
    // Multiplication (*)
    MixedSeqPoint operator*(const double& num) const
    {
      MixedSeqPoint out_seq( 0, v*num, 0, d*num ) ;
      return out_seq ;
    }
    //..................................................
};
//==============================================================
//                      MixedStar Class
//==============================================================
class MixedStar : public Prog
{
  friend class RotationSolver ;
  friend class TOVSolver ;
  friend class MixedSequence ;
  //--------------------------------------------------------------
  private:
    
    struct Region
    {
      Zaki::Math::Range<double> r ;
      std::string name ;

      double Max() const 
      {
        return r.max ;
      }

      double Min() const
      {
        return r.min ;
      }

      Region(const std::string& in_name) 
      : name(in_name) 
      {}

      bool SetRange(const double& a, const double& b)
      {
        if ( a < b )
        {
          r.min = a ;
          r.max = b ; 
          return true ;
        }
        else
        {
          r.max = a ; 
          r.min = b ;
          return false ;
        }
      }
    };
    
    /// The MixedSeqPoint that corresponds to this mixed star
    /// This will be initialized after SurfaceIsReached 
    /// is called, or if the whole TOV is imported from a file
    MixedSeqPoint sequence ;

    /// Rotation solver
    RotationSolver rot_solver ;

    /// Moment of inertia
    double MomI = 0 ;

    /// Do we have a dark core or a dark halo? 
    bool dark_core = true ;

    Region core_region ;
    Region mantle_region ;

    Zaki::Vector::DataSet ds_vis ;
    Zaki::Vector::DataSet ds_dar ;
    Zaki::Vector::DataSet B_vis_integrand ;
    Zaki::Vector::DataSet B_dar_integrand ;
    Zaki::Vector::DataColumn mass_tot_dc ;

    int r_idx = 0 ;
    int m_idx = 1 ;
    int nu_der_idx = 2 ;
    int pre_idx = 3 ;
    int eps_idx = 4 ;
    int rho_idx = 5 ;
    int nu_idx = 6 ;

    std::vector<int> rho_i_v_idx ;
    std::vector<int> rho_i_d_idx ;

    /// Sets the work directory for the member objects
    Prog* SetMemWrkDir(const Zaki::String::Directory&) override ;

    // /// GSL Workspace for performing integrations
    // gsl_integration_workspace *integ_wrk_space = nullptr ;
    
    // /// GSL function for dark baryon integration
    // gsl_function F_B_Dark ;

    // /// GSL function for visible baryon integration
    // gsl_function F_B_Visi ;

    // Zaki::Math::GSLFuncWrapper<MixedStar, double (MixedStar::*)(double)> 
    //   FW_B_Dark;     

    // Zaki::Math::GSLFuncWrapper<MixedStar, double (MixedStar::*)(double)> 
    //   FW_B_Visi;  
  //--------------------------------------------------------------
  public:
    
    /// Default Constructor 
    MixedStar() ;

    /// Constructor from TOV Solutions
    MixedStar(const size_t& in_v_idx, const size_t& in_d_idx,
              const std::vector<TOVPoint>& in_tov_results, 
              const std::vector<TOVPoint>& in_dark_tov_results) ;

    /// Initializes the visible dataset
    void InitVisible(const TOVSolver* in_tov_solver) ;

    /// Initializes the dark dataset
    void InitDark(const TOVSolver* in_tov_solver) ;

    /// Appends tov points to the MixedStar in the core region
    void Append_Core(const TOVPoint& in_tov_pts, 
                const TOVPoint& in_dark_tov_pts) ;

    /// Appends tov points to the MixedStar in the mantle
    /// from the dark solutions
    void Append_Dark_Mantle(const TOVPoint& in_dark_tov_pts) ;

    /// Appends tov points to the MixedStar in the mantle
    /// from the visible solutions
    void Append_Visible_Mantle(const TOVPoint& in_tov_pts) ;

    /// This has to be run so the class
    /// knows when to initialize all the splines
    void SurfaceIsReached(const size_t& in_v_idx, 
                          const size_t& in_d_idx) ;

    /// Reserving the needed space for all datacolumns
    // [ Removed on May 1, 2022 because of the changes 
    //    to the constructors! ]
    // void Reserve(const size_t& space_size) ;

    /// Similar to the destructor
    void Reset() ;

    /// Constructor from file
    // MixedStar(const Zaki::String::Directory& in_tov_file) ;

    /// Destructor
    ~MixedStar() ;

    /// Copy Constructor 
    MixedStar(const MixedStar &) = delete ;

    /// Assignment operator
    MixedStar& operator= (const MixedStar&) = delete;

    /// Visible mass as a function of radius
    double GetMass_Visible(const double& in_r) const ;
    Zaki::Vector::DataColumn* GetMass_Visible() ;
    
    /// Dark mass as a function of radius
    double GetMass_Dark(const double& in_r) const ;
    Zaki::Vector::DataColumn* GetMass_Dark() ;

    /// Total mass as a function of radius
    double GetMass_Total(const double& in_r) const ;

    /// Returns a data column holding total mass values
    /// extending from the origin to the surface of the star 
    Zaki::Vector::DataColumn* GetMass_Total() ;

    /// Visible baryon number density as a function of radius
    double GetRho_Visible(const double& in_r) const ;
    Zaki::Vector::DataColumn* GetRho_Visible() ;

    // Visible baryon number density (fm^{-3}) as a function 
    // of radius for a specific species labeled as (in_label)
    Zaki::Vector::DataColumn* 
      GetRho_i_Visible(const std::string& in_label) ;

    /// Dark baryon number density as a function of radius
    double GetRho_Dark(const double& in_r) const ;
    Zaki::Vector::DataColumn* GetRho_Dark() ;
    
    /// Energy density as a function of radius
    double GetEps_Visible(const double& in_r) const ;
    Zaki::Vector::DataColumn* GetEps_Visible() ;
    
    /// Dark energy density as a function of radius
    double GetEps_Dark(const double& in_r) const ;
    Zaki::Vector::DataColumn* GetEps_Dark() ;

    /// Pressure as a function of radius
    double GetPress_Visible(const double& in_r) const ;
    Zaki::Vector::DataColumn* GetPress_Visible() ;

    /// Dark pressure as a function of radius
    double GetPress_Dark(const double& in_r) const ;
    Zaki::Vector::DataColumn* GetPress_Dark() ;

    /// Returns 'sequence' 
    MixedSeqPoint GetSequence() const ;

    /// Evaluate the metric function
    void EvaluateNu() ;

    /// Metric function as a function of radius
    double GetNu(const double& in_r) const ;
    Zaki::Vector::DataColumn* GetNu() ;

    /// Returns the visible radius dataset
    Zaki::Vector::DataColumn* GetRadius_Visible() ;

    /// Returns the dark radius dataset
    Zaki::Vector::DataColumn* GetRadius_Dark() ;

    /// Returns the larger radius dataset
    Zaki::Vector::DataColumn* GetRadius() ;

    /// Derivative of the metric function as a function of radius
    // double GetNuDerSpline(const double& in_r) ;

    /// Visible baryon number inegrand
    double BaryonNumIntegrand(double in_r) ;

    /// Visible baryon number as a function of radius
    double Find_BaryonNum_Visible(const double& in_r) ;

    /// Visible baryon number
    double Find_BaryonNum_Visible() ;

    /// Dark baryon number inegrand
    double BaryonNumIntegrand_Dark(double in_r) ;

    /// Dark baryon number as a function of radius
    double Find_BaryonNum_Dark(const double& in_r) ;

    /// Dark baryon number
    double Find_BaryonNum_Dark() ;

    // /// Moment of inertia inegrand
    // double MomInertiaIntegrand(double in_r) ;

    /// Total moment of inertia
    double Find_MomInertia() ;

    /// Exports the mixed star profile
    void Export(const Zaki::String::Directory& in_dir) ;
};

//==============================================================
} // CompactStar namespace
//==============================================================
#endif /*CompactStar_MixedStar_H*/
