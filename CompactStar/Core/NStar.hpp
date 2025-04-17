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
 * @file NStar.hpp
 *
 * @brief Typical neutron star.
 *
 * @ingroup Core
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/
// Last edit Nov 3, 2021
#ifndef CompactStar_NStar_H
#define CompactStar_NStar_H

// #include <gsl/gsl_spline.h>

#include <Zaki/String/Directory.hpp>
#include <Zaki/Vector/DataSet.hpp>

#include "CompactStar/Core/RotationSolver.hpp"
#include "CompactStar/Core/Prog.hpp"

//==============================================================
namespace CompactStar
{

// forward declaration 
struct TOVPoint ;     
class RotationSolver ;
class TOVSolver ;


//==============================================================
//             Usual Star Sequence Points Class
//==============================================================
class SeqPoint
{
  public:
    double ec, m, r, pc, b, I ;

    SeqPoint() : ec(0), m(0), r(0), pc(0),
              b(0), I(0) 
      { } ;

    SeqPoint(const double& in_ec, const double& in_m,
            const double& in_r, const double& in_pc,
            const double& in_b, const double& in_I)
            : ec(in_ec), m(in_m), r(in_r), pc(in_pc),
              b(in_b), I(in_I)
      { }

    SeqPoint(const std::vector<double>& in_seq_row)
      {
        if ( in_seq_row.size() != 6 )
          Z_LOG_ERROR("The sequence data point is incomplete.") ;
        else
          {
            ec  = in_seq_row[0] ;
            m   = in_seq_row[1] ;
            r   = in_seq_row[2] ;
            pc  = in_seq_row[3] ;
            b   = in_seq_row[4] ;
            I   = in_seq_row[5] ;
          }
      }

    std::string Str() const
    {
      std::stringstream ss;
      char tmp[150] ;
      snprintf(tmp, sizeof(tmp), "%.8e\t %.8e\t %.8e\t %.8e\t %.8e\t %.8e",
              ec, m, r, pc, b, I) ;
      ss << tmp ;

      return ss.str() ;
    }

    void Reset()
    {
      ec  = 0;
      m   = 0;
      r   = 0;
      pc  = 0;
      b   = 0;
      I   = 0;
    }
    //..................................................
    // Operator overloading
    // Addition (+)
    SeqPoint operator+(const SeqPoint& seq) const
    {
      SeqPoint out_seq( ec+seq.ec, m+seq.m, r+seq.r, 
                        pc+seq.pc, b+seq.b, I+seq.I) ;
      return out_seq ;
    }
    // Multiplication (*)
    SeqPoint operator*(const double& num) const
    {
      SeqPoint out_seq( ec*num, m*num, r*num, 
                        pc*num, b*num, I*num) ;
      return out_seq ;
    }
    //..................................................
};

//==============================================================
class NStar : public Prog
{
  friend class RotationSolver ;
  friend class TOVSolver ;
  friend class Sequence ;
  //--------------------------------------------------------------
  private:

    Zaki::Vector::DataSet ds ;
    Zaki::Vector::DataSet B_integrand ;

    int r_idx = 0 ;
    int m_idx = 1 ;
    int nu_der_idx = 2 ;
    int pre_idx = 3 ;
    int eps_idx = 4 ;
    int rho_idx = 5 ;
    int nu_idx = 6 ;

    std::vector<int> rho_i_idx ;

    // Sets the work directory for the member objects
    CompactStar::Prog* SetMemWrkDir(const 
                  Zaki::String::Directory& in_dir)  ;

    // Zaki::Vector::DataColumn radius ;
    // Zaki::Vector::DataColumn mass ;
    // Zaki::Vector::DataColumn rho ;
    // Zaki::Vector::DataColumn eps ;
    // Zaki::Vector::DataColumn press ;
    // Zaki::Vector::DataColumn nu_der ;
    // Zaki::Vector::DataColumn nu ;

    // This workspace stores state variables for interpolation lookups.
    // It caches the previous value of an index lookup. 
    // When the subsequent interpolation point falls in the same 
    // interval its index value can be returned immediately.
    // gsl_interp_accel *accel = nullptr ;

    // We use cubic spline
    // Cubic spline with natural boundary conditions. 
    // The resulting curve is piecewise cubic on each interval, 
    // with matching first and second derivatives at the supplied data-points.
    //  The second derivative is chosen to be zero at the first point and last point.
    // gsl_spline *mass_r_spline = nullptr ;
    // gsl_spline *rho_r_spline = nullptr ;
    // gsl_spline *eps_r_spline = nullptr ;
    // gsl_spline *press_r_spline = nullptr ;
    // gsl_spline *nu_der_r_spline = nullptr ;
    // gsl_spline *nu_r_spline = nullptr ;

    // std::shared_ptr<NStar> GetSharedPtr() ;
    // void SetSharedPtr(std::shared_ptr<NStar> ) ; 

    /// The SeqPoint that corresponds to this star
    /// This will be initialized after SurfaceIsReached 
    /// is called, or if the whole TOV is imported from a file
    SeqPoint sequence ;

    /// Rotation solver
    RotationSolver rot_solver ;

    /// Moment of inertia
    double MomI = 0 ;

    /// Precision for printing the profile
    int profile_precision = 9 ;

  //--------------------------------------------------------------
  public:
    
    /// Default Constructor 
    NStar();

    /// Constructor from TOV Solutions
    NStar(const std::vector<TOVPoint>& in_tov_results) ;

    /// Initializes the dataset
    void Init(const TOVSolver* in_tov_solver) ;

    /// Appends tov points to the NStar
    void Append(const TOVPoint& in_tov_pts) ;

    /// This has to be run so the class
    /// knows when to initialize all the splines
    void SurfaceIsReached() ;

    // /// Reserving the needed space for all datacolumns
    // void Reserve(const size_t& space_size) ;

    /// Similar to the destructor
    void Reset() ;

    /// Constructor from file
    // NStar(const Zaki::String::Directory& in_tov_file) ;

    /// Destructor
    ~NStar() ;

    /// Copy Constructor 
    NStar(const NStar &) = delete ;

    /// Assignment operator
    NStar& operator= (const NStar&) = delete;

    // /// Mass as a function of radius
    // double GetMass(const double& in_r) const ;

    // /// Baryon number density as a function of radius
    // double GetRho(const double& in_r) const ;

    // /// Energy density as a function of radius
    // double GetEps(const double& in_r) const ;
    
    // /// Pressure as a function of radius
    // double GetPress(const double& in_r) const ;

    // /// Metric function as a function of radius
    // double GetNu(const double& in_r) const ;

    // /// Derivative of the metric function as a function of radius
    // double GetNuDerSpline(const double& in_r) ;

    // /// Total baryon number inegrand
    // double BaryonNumIntegrand(double in_r) ;

    // /// Total baryon number as a function of radius
    // double GetBaryonNum(const double& in_r) ;

    // /// Total baryon number
    // double GetBaryonNum() ;

    // // /// Moment of inertia inegrand
    // // double MomInertiaIntegrand(double in_r) ;

    // /// Total moment of inertia
    // double GetMomInertia() ;

    // ------------------------------
    ///  Mass as a function of radius
    double GetMass(const double& in_r) const ;
    Zaki::Vector::DataColumn* GetMass() ;

    ///  Baryon number density (fm^{-3}) as a function 
    /// of radius for a specific species labeled as (in_label)
    Zaki::Vector::DataColumn* GetRho_i(const 
                                    std::string& in_label) ;

    /// Baryon number density as a function of radius
    double GetRho(const double& in_r) const ;
    Zaki::Vector::DataColumn* GetRho() ;
    
    /// Energy density as a function of radius
    double GetEps(const double& in_r) const ;
    Zaki::Vector::DataColumn* GetEps() ;

    /// Pressure as a function of radius
    double GetPress(const double& in_r) const ;
    Zaki::Vector::DataColumn* GetPress() ;

    /// Returns 'sequence' 
    SeqPoint GetSequence() const ;

    /// Evaluate the metric function
    void EvaluateNu() ;

    /// Metric function as a function of radius
    double GetNu(const double& in_r) const ;
    Zaki::Vector::DataColumn* GetNu() ;

    /// Returns the radius dataset
    Zaki::Vector::DataColumn* GetRadius() ;

    /// Baryon number inegrand
    double BaryonNumIntegrand(double in_r) ;

    ///  Baryon number as a function of radius
    double Find_BaryonNum(const double& in_r) ;

    /// Baryon number
    double Find_BaryonNum() ;

    // /// Moment of inertia inegrand
    // double MomInertiaIntegrand(double in_r) ;

    /// Total moment of inertia
    double Find_MomInertia() ;

    /// Precision for printing the profile
    /// by default it is set to '9' digits
    void SetProfilePrecision(const int& profile_precision) ;

    /// Exports the star profile
    void Export(const Zaki::String::Directory& in_dir) ;
    //------------------------------------------
};

//==============================================================
} // CompactStar namespace
//==============================================================
#endif /*CompactStar_NStar_H*/
