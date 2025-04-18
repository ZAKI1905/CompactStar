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
 * @brief Dark+Visible admixed star.
 *
 * This class represents a mixed compact object consisting of a visible (baryonic)
 * and dark matter component, modeled as concentric spheres in hydrostatic equilibrium.
 * It supports post-processing TOV solutions, metric evaluation, and property export.
 *
 * @ingroup Core
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@eku.edu
 */

#ifndef CompactStar_MixedStar_H
#define CompactStar_MixedStar_H

#include <Zaki/String/Directory.hpp>
#include <Zaki/Vector/DataSet.hpp>

#include "CompactStar/Core/Prog.hpp"
#include "CompactStar/Core/NStar.hpp"

//==============================================================
// struct gsl_integration_workspace ;
//==============================================================
namespace CompactStar
{

// Forward declarations
struct TOVPoint ;     
class TOVSolver ;

//==============================================================
//             Mixed Star Sequence Points Class
//==============================================================
/// @class MixedSeqPoint
/// @brief Represents a data point in the sequence of mixed stars (visible + dark).
class MixedSeqPoint
{
  public:
    
    SeqPoint  v, d ; ///< Visible and dark sequence points
    size_t v_idx = 0, d_idx = 0 ; ///< Indices of visible and dark points in the sequence

    /// Default constructor
    MixedSeqPoint()
      { } ;

    /// Constructor from visible and dark sequence points
    /// @param in_v_idx Index of the visible sequence point
    /// @param in_v Visible sequence point
    /// @param in_d_idx Index of the dark sequence point
    /// @param in_d Dark sequence point
    /// @note The indices are used to identify the sequence points
    ///       in the mixed star sequence.
    MixedSeqPoint(const size_t& in_v_idx, 
                  const SeqPoint& in_v,
                  const size_t& in_d_idx, 
                  const SeqPoint& in_d)
            : v(in_v), d(in_d), v_idx(in_v_idx), 
              d_idx(in_d_idx)
      { }

    /// Constructor from visible and dark sequence points
    /// @param in_v_idx Index of the visible sequence point
    /// @param in_d_idx Index of the dark sequence point
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

    /// @brief Prints the string representation of the sequence point.
    /// @return A string representation of the sequence point.
    /// @note The string format is:
    ///       "v_idx d_idx v_str d_str"
    ///       where v_str and d_str are the string representations of
    ///       the visible and dark sequence points, respectively.
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

    /// @brief Resets the sequence point to default.
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
    /// @brief Adds two mixed sequence points.
    /// @param seq The mixed sequence point to add.
    /// @return A new mixed sequence point that is the sum of this and the given point.
    MixedSeqPoint operator+(const MixedSeqPoint& seq) const
    {
      MixedSeqPoint out_seq(0, v+seq.v, 0, d+seq.d) ;
      return out_seq ;
    }

    // Multiplication (*)
    /// @brief Multiplies a sequence point by a scalar.
    /// @param num The scalar to multiply with.
    /// @return A new mixed sequence point that is the product of this and the given scalar.    
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
/// @class MixedStar
/// @brief Represents a star composed of both visible and dark matter.
///
/// Encapsulates the core and mantle regions of an admixed star, stores TOV solutions,
/// supports data analysis (e.g., moment of inertia), and exports computed profiles.
///
/// The class is designed to work with the TOVSolver and RotationSolver classes.
/// It provides methods to initialize, append TOV points, and compute various properties
/// such as mass, baryon number, and moment of inertia.
/// It also allows for the evaluation of the metric function and its derivative.
/// The class is capable of exporting the mixed star profile to a file.
/// @ingroup Core
/// @author Mohammadreza Zakeri
class MixedStar : public Prog
{
  friend class RotationSolver ;
  friend class TOVSolver ;
  friend class MixedSequence ;
  //--------------------------------------------------------------
  private:

    /// @struct Region
    /// @brief Defines a radial region within the star (e.g., core or mantle).
    /// @details Each region has a range defined by a minimum and maximum radius.
    /// @note The region is defined by a name (e.g., "core", "mantle") and a range.
    /// @note The region is used to determine the boundaries for TOV point appending.
    struct Region
    {
      /// @brief Range of the region
      /// @details The range is defined by a minimum and maximum radius.
      Zaki::Math::Range<double> r ;

      /// @brief Name of the region
      /// @details The name is used to identify the region (e.g., "core", "mantle").
      std::string name ;
      
      /// @brief Returns the maximum radius of the region.
      /// @return The maximum radius of the region.
      double Max() const 
      {
        return r.max ;
      }

      /// @brief Returns the minimum radius of the region.
      /// @return The minimum radius of the region.
      double Min() const
      {
        return r.min ;
      }

      /// @brief Constructor for the region.
      /// @param in_name The name of the region (e.g., "core", "mantle").
      Region(const std::string& in_name) 
      : name(in_name) 
      {}

      /// @brief Sets the range of the region.
      /// @param a The minimum radius of the region.
      /// @param b The maximum radius of the region.
      /// @return True if the range was set successfully, false otherwise.
      /// @note The method ensures that the minimum radius is less than the maximum radius.
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
    //--------------------------------------------------------------

    /// @brief The sequence point for the mixed star
    /// @details This point is used to store the sequence point
    /// for the mixed star, which includes both visible and dark components.
    /// @note The MixedSeqPoint is initialized after the surface is reached
    ///       or if the whole TOV is imported from a file.
    MixedSeqPoint sequence ;


    /// Rotation solver
    /// @brief The rotation solver computes the first order corrections due to rotation using Hartle's relativistic equations.
    /// @details The rotation solver is used to compute the moment of inertia of the mixed star.
    RotationSolver rot_solver ;

    /// @brief The moment of inertia of the mixed star.
    double MomI = 0 ;

    /// @brief Indicates whether the core is dark or visible.
    /// @details If true, the core is dark. If false, the core is visible.
    bool dark_core = true ;

    /// @brief The core region of the mixed star.
    Region core_region ;
    
    /// @brief The mantle region of the mixed star.
    Region mantle_region ;


    /// @brief Dataset for visible matter sector.
    Zaki::Vector::DataSet ds_vis ;
    
    /// @brief Dataset for dark matter sector.
    Zaki::Vector::DataSet ds_dar ;

    /// @brief Dataset storing the visible baryon number density integrand.
    Zaki::Vector::DataSet B_vis_integrand ;

    /// @brief Dataset storing the dark baryon number density integrand.
    Zaki::Vector::DataSet B_dar_integrand ;

    /// @brief Datacolumn storing the total mass of the mixed star.
    Zaki::Vector::DataColumn mass_tot_dc ;

    // Index bookkeeping
    /// @brief Index of radius.
    int r_idx = 0 ;

    /// @brief Index of mass.
    int m_idx = 1 ;

    /// @brief Index metric function derivative d[nu(r)]/dr.
    int nu_der_idx = 2 ;

    /// @brief Index of pressure.
    int pre_idx = 3 ;

    /// @brief Index of energy density.
    int eps_idx = 4 ;

    /// @brief Index of baryon number density.
    int rho_idx = 5 ;

    /// @brief Index of the metric function nu(r).
    int nu_idx = 6 ;

    std::vector<int> rho_i_v_idx ;
    std::vector<int> rho_i_d_idx ;

    /// @brief Sets the working directory for the mixed star.
    /// @details This method sets the working directory for the mixed star and its member objects.
    /// @param in_dir The directory to set.
    Prog* SetMemWrkDir(const Zaki::String::Directory&) override ;
  //--------------------------------------------------------------
  public:
    
    /// @brief Default constructor for the MixedStar class.
    MixedStar() ;

    /// Constructor from TOV Solutions
    /// @brief Initializes the MixedStar with TOV solutions.
    /// @param in_v_idx Index of the visible TOV solution.
    /// @param in_d_idx Index of the dark TOV solution.
    /// @param in_tov_results Vector of TOV points for the visible solution.
    /// @param in_dark_tov_results Vector of TOV points for the dark solution.
    MixedStar(const size_t& in_v_idx, const size_t& in_d_idx,
              const std::vector<TOVPoint>& in_tov_results, 
              const std::vector<TOVPoint>& in_dark_tov_results) ;

    /// @brief Initializes the visible dataset with TOV solutions.
    void InitVisible(const TOVSolver* in_tov_solver) ;

    /// @brief Initializes the dark dataset with TOV solutions.
    void InitDark(const TOVSolver* in_tov_solver) ;

    /// @brief Appends TOV points to the core region of the mixed star.
    /// @param in_tov_pts Vector of TOV points for the visible solution.
    /// @param in_dark_tov_pts Vector of TOV points for the dark solution.
    void Append_Core(const TOVPoint& in_tov_pts, 
                const TOVPoint& in_dark_tov_pts) ;

    /// @brief Appends TOV points to the mantle region of the mixed star from the dark solutions
    /// @param in_dark_tov_pts Vector of TOV points for the dark solution.
    void Append_Dark_Mantle(const TOVPoint& in_dark_tov_pts) ;

    /// @brief Appends TOV points to the mantle region of the mixed star from the visible solutions.
    /// @param in_tov_pts Vector of TOV points for the visible solution.
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
    /// @brief Resets the mixed star to its initial state.
    void Reset() ;

    /// Constructor from file
    // MixedStar(const Zaki::String::Directory& in_tov_file) ;

    /// @brief Destructor for the MixedStar class.
    ~MixedStar() ;

    /// @brief Copy constructor for the MixedStar class.
    /// @details This constructor is deleted to prevent copying of the MixedStar object.
    MixedStar(const MixedStar &) = delete ;

    /// @brief Assignment operator for the MixedStar class.
    /// @details This operator is deleted to prevent assignment of the MixedStar object.
    MixedStar& operator= (const MixedStar&) = delete;

    /// @brief Returns the mass of the visible component at a given radius.
    double GetMass_Visible(const double& in_r) const ;

    /// @brief Returns the mass of the visible component as a function of radius.
    Zaki::Vector::DataColumn* GetMass_Visible() ;
    
    /// @brief Returns the mass of the dark component at a given radius.
    double GetMass_Dark(const double& in_r) const ;

    /// @brief Returns the mass of the dark component as a function of radius.
    Zaki::Vector::DataColumn* GetMass_Dark() ;

    /// @brief Returns the total mass of the mixed star at a given radius.
    /// @details The total mass is the sum of the visible and dark masses.
    /// @param in_r The radius at which to evaluate the total mass.
    double GetMass_Total(const double& in_r) const ;

    /// Returns a data column holding total mass values
    /// extending from the origin to the surface of the star 
    /// @brief Returns the total mass of the mixed star as a function of radius.
    Zaki::Vector::DataColumn* GetMass_Total() ;

    /// Visible baryon number density as a function of radius
    /// @brief Returns the visible baryon number density at a given radius.
    /// @param in_r The radius at which to evaluate the baryon number density.
    double GetRho_Visible(const double& in_r) const ;

    /// @brief Returns the visible baryon number density as a function of radius.
    /// @details The baryon number density is the number of baryons per unit volume.
    Zaki::Vector::DataColumn* GetRho_Visible() ;

    /// @brief Returns the visible baryon number density fm^{-3} for a specific species at a given radius.
    /// @param in_r The radius at which to evaluate the baryon number density.
    /// @param in_label The label of the baryon species (e.g., "n", "p", "e").
    Zaki::Vector::DataColumn* 
      GetRho_i_Visible(const std::string& in_label) ;

    /// @brief Returns the dark baryon number density at a given radius.
    /// @param in_r The radius at which to evaluate the baryon number density.
    double GetRho_Dark(const double& in_r) const ;

    /// @brief Returns the dark baryon number density as a function of radius.
    Zaki::Vector::DataColumn* GetRho_Dark() ;
    
    /// @brief Returns the energy density at a given radius.
    /// @param in_r The radius at which to evaluate the energy density.
    double GetEps_Visible(const double& in_r) const ;

    /// @brief Returns the visible energy density as a function of radius.
    Zaki::Vector::DataColumn* GetEps_Visible() ;
    
    /// @brief Returns the dark energy density at a given radius.
    /// @param in_r The radius at which to evaluate the energy density.
    double GetEps_Dark(const double& in_r) const ;

    /// @brief Returns the dark energy density as a function of radius.
    Zaki::Vector::DataColumn* GetEps_Dark() ;

    /// @brief Returns the visible pressure at a given radius.
    /// @param in_r The radius at which to evaluate the pressure.
    double GetPress_Visible(const double& in_r) const ;

    /// @brief Returns the visible pressure as a function of radius.
    Zaki::Vector::DataColumn* GetPress_Visible() ;

    /// @brief Returns the dark pressure at a given radius.
    /// @param in_r The radius at which to evaluate the pressure.
    double GetPress_Dark(const double& in_r) const ;

    /// @brief Returns the dark pressure as a function of radius.
    Zaki::Vector::DataColumn* GetPress_Dark() ;

    /// Returns 'sequence' 
    /// @brief Returns the sequence point of the mixed star.
    MixedSeqPoint GetSequence() const ;

    /// @brief Evaluates the metric function by integrating nu'(r) over the radius.
    void EvaluateNu() ;

    /// @brief Returns the metric function at a given radius.
    /// @param in_r The radius at which to evaluate the metric function.
    double GetNu(const double& in_r) const ;

    /// @brief Returns the metric function as a function of radius.
    Zaki::Vector::DataColumn* GetNu() ;

    /// @brief Returns a pointer to the visible radius dataset.
    Zaki::Vector::DataColumn* GetRadius_Visible() ;

    /// @brief Returns a pointer to the dark radius dataset.
    Zaki::Vector::DataColumn* GetRadius_Dark() ;

    /// @brief Returns a pointer to the radius dataset of the mixed star.
    Zaki::Vector::DataColumn* GetRadius() ;

    /// @brief Returns the value for visible baryon number integrand at a given radius.
    /// @param in_r The radius at which to evaluate the baryon number integrand.
    double BaryonNumIntegrand(double in_r) ;

    /// @brief Returns the visible baryon number at a given radius.
    /// @param in_r The radius at which to evaluate the baryon number.
    double Find_BaryonNum_Visible(const double& in_r) ;

    /// Visible baryon number
    /// @brief Returns the total visible baryon number.
    double Find_BaryonNum_Visible() ;

    /// @brief Returns the value for dark baryon number integrand at a given radius.
    /// @param in_r The radius at which to evaluate the baryon number integrand.
    double BaryonNumIntegrand_Dark(double in_r) ;

    /// @brief Returns the dark baryon number at a given radius.
    /// @param in_r The radius at which to evaluate the baryon number.
    double Find_BaryonNum_Dark(const double& in_r) ;

    /// @brief Returns the total dark baryon number.
    double Find_BaryonNum_Dark() ;


    /// @brief Returns the total moment of inertia of the mixed star.
    double Find_MomInertia() ;

    /// @brief Exports the mixed star profile to a file.
    /// @param in_dir The directory where the profile should be exported.
    void Export(const Zaki::String::Directory& in_dir) ;
};

//==============================================================
} // CompactStar namespace
//==============================================================
#endif /*CompactStar_MixedStar_H*/
