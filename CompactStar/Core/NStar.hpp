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
 * @brief Definition of neutron star profile and sequence handling.
 *
 * This file declares the SeqPoint struct for storing star sequence points
 * and the NStar class for managing TOV solutions, rotation, and profile export.
 *
 * @ingroup Core
 * @author Mohammadreza Zakeri
 * @contact M.Zakeri@eku.edu
 */
// Last edit Nov 3, 2021
#ifndef CompactStar_NStar_H
#define CompactStar_NStar_H


#include <Zaki/String/Directory.hpp>
#include <Zaki/Vector/DataSet.hpp>

#include "CompactStar/Core/RotationSolver.hpp"
#include "CompactStar/Core/Prog.hpp"

//==============================================================
namespace CompactStar
{

/**
 * @struct TOVPoint
 * @brief Forward declaration for Tolman-Oppenheimer-Volkoff solution point.
 */
struct TOVPoint ;     
class RotationSolver ;
class TOVSolver ;


//==============================================================
//             Usual Star Sequence Points Class
//==============================================================
/**
 * @class SeqPoint
 * @brief Holds a single data point in a star sequence.
 *
 * Stores energy density, mass, radius, central pressure, baryon number,
 * and moment of inertia for a neutron star model.
 */
class SeqPoint
{
  public:
    /** @brief Energy density (EC). */
    double ec;

    /** @brief Gravitational mass (M). */
    double m;

    /** @brief Radius (R). */
    double r;

    /** @brief Central pressure (PC). */
    double pc;

    /** @brief Baryon number integral (B). */
    double b;

    /** @brief Moment of inertia (I). */
    double I;

    /**
     * @brief Default constructor, initializes all values to zero.
     */
    SeqPoint() : ec(0), m(0), r(0), pc(0),
              b(0), I(0) 
      { } ;

    /**
     * @brief Parameterized constructor.
     *
     * @param in_ec  Energy density.
     * @param in_m   Mass.
     * @param in_r   Radius.
     * @param in_pc  Central pressure.
     * @param in_b   Baryon number.
     * @param in_I   Moment of inertia.
     */
    SeqPoint(const double& in_ec, const double& in_m,
            const double& in_r, const double& in_pc,
            const double& in_b, const double& in_I)
            : ec(in_ec), m(in_m), r(in_r), pc(in_pc),
              b(in_b), I(in_I)
      { }

    /**
     * @brief Construct from a sequence row vector.
     *
     * Expects a vector of length 6: [ec, m, r, pc, b, I].
     * Logs an error if size mismatch.
     *
     * @param in_seq_row  Input row vector.
     */
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

    /**
     * @brief Format the sequence point as a string.
     *
     * @return Tab-delimited string of values in scientific notation.
     */
    std::string Str() const
    {
      std::stringstream ss;
      char tmp[150] ;
      snprintf(tmp, sizeof(tmp), "%.8e\t %.8e\t %.8e\t %.8e\t %.8e\t %.8e",
              ec, m, r, pc, b, I) ;
      ss << tmp ;

      return ss.str() ;
    }

    /**
     * @brief Reset all sequence values to zero.
     */
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

    /**
     * @brief Addition operator.
     *
     * Adds corresponding fields of two SeqPoint objects.
     *
     * @param seq  Sequence point to add.
     * @return     New SeqPoint representing the sum.
     */
    SeqPoint operator+(const SeqPoint& seq) const
    {
      SeqPoint out_seq( ec+seq.ec, m+seq.m, r+seq.r, 
                        pc+seq.pc, b+seq.b, I+seq.I) ;
      return out_seq ;
    }
    
    /**
     * @brief Scalar multiplication operator.
     *
     * Multiplies all fields by a scalar.
     *
     * @param num  Scalar multiplier.
     * @return     New SeqPoint scaled by num.
     */
    SeqPoint operator*(const double& num) const
    {
      SeqPoint out_seq( ec*num, m*num, r*num, 
                        pc*num, b*num, I*num) ;
      return out_seq ;
    }
    //..................................................
};

//==============================================================
//                        NStar Class
//==============================================================

/**
 * @class NStar
 * @brief Main neutron star class for handling TOV solutions and rotation.
 *
 * Inherits from Prog to provide utility functions, manages interpolation
 * datasets, computes baryon number and moment of inertia, and exports profiles.
 */
class NStar : public Prog
{
  friend class RotationSolver ;
  friend class TOVSolver ;
  friend class Sequence ;
  //--------------------------------------------------------------
  private:

    /** @brief Dataset holding radius, mass, pressure, etc. */
    Zaki::Vector::DataSet ds ;

    /** @brief Integrand dataset for baryon number calculation. */
    Zaki::Vector::DataSet B_integrand ;

    int r_idx      = 0; /**< Radius index in dataset. */
    int m_idx      = 1; /**< Mass index.          */
    int nu_der_idx = 2; /**< Derivative of metric function index. */
    int pre_idx    = 3; /**< Pressure index.      */
    int eps_idx    = 4; /**< Energy density index.*/
    int rho_idx    = 5; /**< Baryon density index.*/
    int nu_idx     = 6; /**< Metric function index.*/

    std::vector<int> rho_i_idx ; /**< Indices for individual species densities. */

    /**
     * @brief Set working directory for internal objects.
     *
     * @param in_dir  Input directory path.
     * @return        Pointer to Prog base for chaining calls.
     */
    CompactStar::Prog* SetMemWrkDir(const 
                  Zaki::String::Directory& in_dir)  ;

    /// The SeqPoint that corresponds to this star
    /// This will be initialized after SurfaceIsReached 
    /// is called, or if the whole TOV is imported from a file
    SeqPoint sequence ; /**< Cached sequence point after surface is reached. */

    RotationSolver rot_solver ; /**< Rotation solver instance. */

    double MomI = 0 ; /**< Cached moment of inertia. */

    int profile_precision = 9 ; /**< Digits of precision when printing profiles. */

  //--------------------------------------------------------------
  public:
    
    /**
     * @brief Default constructor.
     */ 
    NStar();

    /**
     * @brief Construct from TOV solution points.
     *
     * Initializes dataset from provided TOV results.
     *
     * @param in_tov_results  Vector of TOVPoint structures.
     */
    NStar(const std::vector<TOVPoint>& in_tov_results) ;

    /**
     * @brief Initialize internal datasets from a TOV solver.
     *
     * Must be called before appending points.
     *
     * @param in_tov_solver  Pointer to initialized TOVSolver.
     */
    void Init(const TOVSolver* in_tov_solver) ;

    /**
     * @brief Append a single TOV solution point.
     *
     * @param in_tov_pts  TOVPoint to append.
     */
    void Append(const TOVPoint& in_tov_pts) ;

    /// This has to be run so the class
    /// knows when to initialize all the splines
    /**
     * @brief Signal that surface has been reached, enabling spline init.
     */
    void SurfaceIsReached() ;

    /**
     * @brief Reset the NStar object to initial state.
     */
    void Reset() ;

    /**
     * @brief Destructor cleans up any allocated resources.
     */
    ~NStar() ;

    // Deleted copy operations to enforce unique ownership
    NStar(const NStar &) = delete ;

    /// Assignment operator
    NStar& operator= (const NStar&) = delete;

    // ------------------------------
    //  Getters
    // ------------------------------

    /**
     * @brief Get interpolated mass at a given radius.
     *
     * @param in_r  Radius value.
     * @return      Mass at radius.
     */
    double GetMass(const double& in_r) const ;

    /**
     * @brief Get the DataColumn of mass values.
     *
     * @return Pointer to mass DataColumn.
     */
    Zaki::Vector::DataColumn* GetMass() ;

    ///  Baryon number density (fm^{-3}) as a function 
    /// of radius for a specific species labeled as (in_label)

    /**
     * @brief Get baryon density (fm^{-3}) DataColumn
     *   as a function of radius for a labeled species.
     *
     * @param in_label  Species label string.
     * @return          Pointer to density DataColumn.
     */
    Zaki::Vector::DataColumn* GetRho_i(const 
                                    std::string& in_label) ;

    /**
     * @brief Get interpolated total baryon density at radius.
     *
     * @param in_r  Radius value.
     * @return      Baryon density.
     */
    double GetRho(const double& in_r) const ;

    /**
     * @brief Get DataColumn of total baryon density.
     *
     * @return Pointer to density DataColumn.
     */
    Zaki::Vector::DataColumn* GetRho() ;
    
    /**
     * @brief Get interpolated energy density at radius.
     *
     * @param in_r  Radius value.
     * @return      Energy density.
     */
    double GetEps(const double& in_r) const ;

    /**
     * @brief Get DataColumn of energy density.
     *
     * @return Pointer to energy density DataColumn.
     */
    Zaki::Vector::DataColumn* GetEps() ;

    /**
     * @brief Get interpolated pressure at radius.
     *
     * @param in_r  Radius value.
     * @return      Pressure.
     */
    double GetPress(const double& in_r) const ;

    /**
     * @brief Get DataColumn of pressure values.
     *
     * @return Pointer to pressure DataColumn.
     */
    Zaki::Vector::DataColumn* GetPress() ;

    /**
     * @brief Retrieve the computed sequence point.
     *
     * @return SeqPoint instance.
     */
    SeqPoint GetSequence() const ;

    /**
     * @brief Compute metric function spline values.
     */
    void EvaluateNu() ;

    /**
     * @brief Get interpolated metric function at radius.
     *
     * @param in_r  Radius value.
     * @return      Metric function value.
     */      
    double GetNu(const double& in_r) const ;

    /**
     * @brief Get DataColumn of metric function values.
     *
     * @return Pointer to metric DataColumn.
     */
    Zaki::Vector::DataColumn* GetNu() ;

    /**
     * @brief Get DataColumn of radius values.
     *
     * @return Pointer to radius DataColumn.
     */
    Zaki::Vector::DataColumn* GetRadius() ;

    /**
     * @brief Compute integrand for baryon number at radius.
     *
     * @param in_r  Radius value.
     * @return      Integrand value.
     */
    double BaryonNumIntegrand(double in_r) ;

    /**
     * @brief Compute baryon number as a function of radius.
     *
     * @param in_r  Radius up to which to integrate.
     * @return      Baryon number.
     */
    double Find_BaryonNum(const double& in_r) ;

    /**
     * @brief Compute total baryon number.
     *
     * @return Baryon number.
     */
    double Find_BaryonNum() ;

    /**
     * @brief Compute total moment of inertia.
     *
     * @return Moment of inertia.
     */
    double Find_MomInertia() ;

    /**
     * @brief Set precision for exporting profile values.
     *
     * @param precision  Number of digits.
     * @note Default is 9.
     */
    void SetProfilePrecision(const int& profile_precision) ;

    /**
     * @brief Export the star profile to directory.
     *
     * @param in_dir  Output directory path.
     */
    void Export(const Zaki::String::Directory& in_dir) ;
    //------------------------------------------
};

//==============================================================
} // CompactStar namespace
//==============================================================
#endif /*CompactStar_NStar_H*/
